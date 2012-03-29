#include "pmdOCL.h"

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

    /* Read in file */
    printf("\n\n    File is __%s__\n\n", cmd.filename);
    sim = fromFileASCII(cmd.filename, pot);
    if ( ! sim ) simulationAbort(-3,(char *) "Input file does not exist");

    printf("    total atoms is: %d\n", sim->ntot);

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
        ts = timeNow();
        eone = nTimeSteps(nsteps,sim,dt);
        te = timeNow();
    }

    printf(" %3d %30.20f computed in %.3fs (%8.4f us/atom for %d atoms)\n",
            iter*nsteps, sim->e,(te-ts),1.0e6*(te-ts)/(double)(nsteps+1)/(double)sim->ntot,sim->ntot);
    return NULL;
}
