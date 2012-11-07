#include "mytype.h"
#include "pmdOCL.h"
#include "cheby.h"
#include "ic_fcc.h"
#include "read.h"

simflat_t *initSimFromCmdLine(int argc, char **argv, int *eam_flag, int *gpu_flag) {
    command_t cmd;
    struct pmd_base_potential_t *pot;
    simflat_t *sim;
    SimThreadData* threadData;

    printf("    Double Precision=%s\n", (sizeof(real_t)==4?"false":"true") );

    /* get command line params */
    parseCommandLine(&cmd, argc, argv);

    printCmd(&cmd);

    *eam_flag = cmd.doeam;
    *gpu_flag = cmd.usegpu;

    /* decide whether to get LJ or EAM potentials */
    if(cmd.doeam) 
    {
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
        printf("cutoff = %e\n", new_pot->cutoff);
        printf("mass = %e\n", new_pot->mass);
        printf("phi potential:\n");
        printf("n = %d\n", new_pot->phi->n);
        cmd.lat = new_pot->lat;
        printf("cmd.lat = %e\n", cmd.lat);
    } else {
        ljpotential_t *new_pot;
        new_pot = (ljpotential_t*) pot;
        printf("LJ potential values:\n");
        printf("cutoff = %e\n", new_pot->cutoff);
        printf("mass = %e\n", new_pot->mass);
        cmd.lat = 1.122462048*new_pot->cutoff ;// * 1.53;      // This is for Lennard-Jones
        printf("cmd.lat = %e\n", cmd.lat);
    }

    if(strcmp(cmd.filename, "")) {
        /* Read in file */
        printf("\n\n    File is __%s__\n\n", cmd.filename);
        sim = fromFileASCII(cmd, pot);
        if ( ! sim ) simulationAbort(-3,(char *) "Input file does not exist");
    } else {
        sim = create_fcc_lattice(cmd, pot);
    }
    printf("Initial condition set\n");

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
    printf("box factor is (%e, %e, %e)\n", 
            sim->boxsize[0]/pot->cutoff,
            sim->boxsize[1]/pot->cutoff,
            sim->boxsize[2]/pot->cutoff);

    /* initial output for consistency check */
    reBoxAll(sim);

    return sim;
}


void *do_compute_work(simflat_t* sim)
{
    int iter;
    double eone;
    double ts, te;
    int niter = 10;
    int nsteps = 10;
    double dt;

    ts = timeNow();
    computeForce(sim); 
    te = timeNow();

    printIt(sim, stdout);

    /* convert the timestep */
    dt = 1.0e-15;

    for(iter=0; iter<niter; iter++) {
        double ns;
        ns = (iter==0?1.0:(double)(1+nsteps));
        /* note the 'nsteps+1' below.
         * this is because we have an extra
         * force call at the end of timesteps()
         */
        printf(" %3d %30.20f computed in %.3es (%8.4f us/atom for %d atoms)\n",
                iter*nsteps, sim->e,(te-ts),1.0e6*(te-ts)/(double)ns/(double)sim->ntot,sim->ntot);
        printf("Virial stress = %g\n", sim->stress);
        ts = timeNow();
        eone = nTimeSteps(nsteps,sim,dt);
        te = timeNow();
        //printIt(sim, stdout);
    }

    printf(" %3d %30.20f computed in %.3es (%8.4f us/atom for %d atoms)\n",
            iter*nsteps, sim->e,(te-ts),1.0e6*(te-ts)/(double)(nsteps+1)/(double)sim->ntot,sim->ntot);
    return NULL;
}
