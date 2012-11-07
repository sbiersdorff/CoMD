/*

   Copyright (c) 2011, Los Alamos National Security, LLC All rights
   reserved.  Copyright 2011. Los Alamos National Security, LLC. This
   software was produced under U.S. Government contract DE-AC52-06NA25396
   for Los Alamos National Laboratory (LANL), which is operated by Los
   Alamos National Security, LLC for the U.S. Department of Energy. The
   U.S. Government has rights to use, reproduce, and distribute this
   software.

   NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY
   WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF
   THIS SOFTWARE.

   If software is modified to produce derivative works, such modified
   software should be clearly marked, so as not to confuse it with the
   version available from LANL.

   Additionally, redistribution and use in source and binary forms, with
   or without modification, are permitted provided that the following
   conditions are met:

   · Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

   · Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   · Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
   CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
   BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
   FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS
   ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY
   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
   GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
   IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
   OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
   ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#include "pmd.h"
#include "cheby.h"
#include "ic_fcc.h"
#include "docompute.h"
#include "read.h"

void printArray(real_t* array, int n, char *name)
{
    int i;
    printf("%s;\n", name);
    for (i=0;i<n;i++){
        printf("%d, %17.9e\n", i, array[i]);
    }
}


simflat_t *initSimFromCmdLine(int argc, char **argv) {
    command_t cmd;
    struct pmd_base_potential_t *pot;
    simflat_t *sim;
    SimThreadData* threadData;

    printf("    Double Precision=%s\n", (sizeof(real_t)==4?"false":"true") );

    /* get command line params */
    parseCommandLine(&cmd, argc, argv);

    printCmd(&cmd);

    /* decide whether to get LJ or EAM potentials */
    if(cmd.doeam) 
    {
        printf("Setting eam potential\n");
        pot = setEamPot(cmd.potdir, cmd.potname);
    } else {
        pot = (pmd_base_potential_t *) getLJPot();
    }

    if ( ! pot ) simulationAbort(-2,(char *) "Unable to initialize potential");

    printf("Simulation data\n");
    if (cmd.doeam) {
        eampotential_t *new_pot;
        new_pot = (eampotential_t*) pot;
        printf("EAM potential values:\n");
        printf("cutoff = "EMT1"\n", new_pot->cutoff);
        printf("mass = "EMT1"\n", new_pot->mass);
        printf("phi potential:\n");
        printf("n = %d\n", new_pot->phi->n);
        cmd.lat = new_pot->lat;
        printf("cmd.lat = "EMT1"\n", cmd.lat);
        printf("%e, %e\n", new_pot->phi->x0, new_pot->phi->xn);
#if (DIAG_LEVEL > 1)
        printf("x0 = "EMT1"\n", new_pot->phi->x0);
        printf("xn = "EMT1"\n", new_pot->phi->xn);
        printArray(new_pot->phi->values, new_pot->phi->n, "phi ref");
        printArray(new_pot->rho->values, new_pot->rho->n, "rho ref");
        printArray(new_pot->f->values, new_pot->f->n, "f ref");
#endif
    } else {
        ljpotential_t *new_pot;
        new_pot = (ljpotential_t*) pot;
        printf("LJ potential values:\n");
        printf("cutoff = "EMT1"\n", new_pot->cutoff);
        printf("mass = "EMT1"\n", new_pot->mass);
        cmd.lat = 1.122462048*new_pot->cutoff ;// * 1.53;      // This is for Lennard-Jones
    }

    if(strcmp(cmd.filename, "")) {
        /* Read in file */
        printf("\n\n    File is __%s__\n\n", cmd.filename);
        sim = fromFileASCII(cmd, pot);
        if ( ! sim ) simulationAbort(-3,(char *) "Input file does not exist");
    } else {
        printf("cmd.lat = "EMT1"\n", cmd.lat);
        sim = create_fcc_lattice(cmd, pot);
    }

    if (cmd.doeam) {
        eampotential_t *new_pot = (eampotential_t*) pot;
        printf("%e, %e\n", new_pot->phi->x0, new_pot->phi->xn);
        sim->ch_pot = setChebPot(new_pot, 32);
#if (DIAG_LEVEL > 0)
        printf("Chebychev coefficients:\n");
        fflush(stdout);
        printf("%d, %d, %d\n", 
                sim->ch_pot->phi->n,
                sim->ch_pot->rho->n,
                sim->ch_pot->f->n);
        fflush(stdout);

        printArray(sim->ch_pot->phi->values, sim->ch_pot->phi->n, "phi");
        printArray(sim->ch_pot->rho->values, sim->ch_pot->rho->n, "rho");
        printArray(sim->ch_pot->f->values,   sim->ch_pot->f->n,   "f");
        printArray(sim->ch_pot->dphi->values, sim->ch_pot->dphi->n, "dphi");
        printArray(sim->ch_pot->drho->values, sim->ch_pot->drho->n, "drho");
        printArray(sim->ch_pot->df->values,   sim->ch_pot->df->n,   "df");
#endif
    }

    printf("    total atoms is: %d\n", sim->ntot);
    printf("box factor is ("EMT1", "EMT1", "EMT1")\n", 
            sim->boxsize[0]/pot->cutoff,
            sim->boxsize[1]/pot->cutoff,
            sim->boxsize[2]/pot->cutoff);

    /* initial output for consistency check */
    reBoxAll(sim);

    printf("Initialization finished\n");
    return sim;
}


void *do_compute_work(SimThreadData *data)
{
    int iter;
    double eone;
    double ts, te;
    int niter = 20;
    int nsteps = 10;
    double dt;

#if (DIAG_LEVEL > 1)
    printf("Chebychev coefficients:\n");
    fflush(stdout);
    printf("%d, %d, %d\n", 
            data->sim->ch_pot->phi->n,
            data->sim->ch_pot->rho->n,
            data->sim->ch_pot->f->n);
    fflush(stdout);
#endif

    printf("Starting simulation\n");
    fflush(stdout);
    ts = timeNow();
    computeForce(data->sim); 
    te = timeNow();

    /* convert the timestep */
    dt = 1.0e-15;

    printf("Starting timesteps\n");
    fflush(stdout);
    for(iter=0; iter<niter; iter++) {
        double ns;
        ns = (iter==0?1.0:(double)(1+nsteps));
        /* note the 'nsteps+1' below.
         * this is because we have an extra
         * force call at the end of timesteps()
         */
        printf(" %3d %30.20g computed in %.3fs (%8.4f us/atom for %d atoms)\n",
                iter*nsteps, data->sim->e,(te-ts),1.0e6*(te-ts)/(double)ns/(double)data->sim->ntot,data->sim->ntot);
        printf("Virial stress = %g\n", data->sim->stress);
        ts = timeNow();
        eone = nTimeSteps(nsteps,data->sim,dt);
        te = timeNow();

#ifdef USE_IN_SITU_VIZ
        data->viz->updateData(data->sim);
#endif
    }

    printf(" %3d %30.20f computed in %.3fs (%8.4f us/atom for %d atoms)\n",
            iter*nsteps, data->sim->e,(te-ts),1.0e6*(te-ts)/(double)(nsteps+1)/(double)data->sim->ntot,data->sim->ntot);
    return NULL;
}
