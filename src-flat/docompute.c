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

simflat_t *initSimFromCmdLine(int argc, char **argv) {
    command_t cmd;
    struct pmd_base_potential_t *pot;
    simflat_t *sim;
    SimThreadData* threadData;

    printf("    Double Precision=%s\n", (sizeof(real_t)==4?"false":"true") );

    /* get command line params */
    parseCommandLine(&cmd, argc, argv);

    /* decide whether to get LJ or EAM potentials */
    if(cmd.doeam) 
    {
        printf("Setting eam potential\n");
        pot = setEamPot(cmd.potdir, cmd.potname);
    } else {
        pot = (pmd_base_potential_t *) getLJPot();
    }

    if ( ! pot ) simulationAbort(-2,(char *) "Unable to initialize potential");

    /* Read in file */
    printf("\n\n    File is __%s__\n\n", cmd.filename);
    sim = fromFileASCII(cmd.filename, pot);
    if ( ! sim ) simulationAbort(-3,(char *) "Input file does not exist");

    printf("    total atoms is: %d\n", sim->ntot);

    printf("Simulation data\n");
    if (cmd.doeam) {
        eampotential_t *new_pot;
        new_pot = (eampotential_t*) pot;
        printf("EAM potential values:\n");
        printf("cutoff = %e\n", new_pot->cutoff);
        printf("mass = %e\n", new_pot->mass);
        printf("phi potential:\n");
        printf("n = %d\n", new_pot->phi->n);
    } else {
        ljpotential_t *new_pot;
        new_pot = (ljpotential_t*) pot;
        printf("LJ potential values:\n");
        printf("cutoff = %e\n", new_pot->cutoff);
        printf("mass = %e\n", new_pot->mass);
    }

    /* initial output for consistency check */
    reBoxAll(sim);

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

    ts = timeNow();
    computeForce(data->sim); 
    te = timeNow();

    /* convert the timestep */
    dt = 1.0e-15;

    for(iter=0; iter<niter; iter++) {
        double ns;
        ns = (iter==0?1.0:(double)(1+nsteps));
        /* note the 'nsteps+1' below.
         * this is because we have an extra
         * force call at the end of timesteps()
         */
        printf(" %3d %30.20g computed in %.3fs (%8.4f us/atom for %d atoms)\n",
                iter*nsteps, data->sim->e,(te-ts),1.0e6*(te-ts)/(double)ns/(double)data->sim->ntot,data->sim->ntot);
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
