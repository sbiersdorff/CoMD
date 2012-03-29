/**
 * a simple md simulator
 **/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

#include "pmd.h"
#include "eamTypes.h"

int main(int argc, char **argv) {
  command_t cmd;
  parseCommandLine(&cmd, argc, argv);
  return cmdPPU(&cmd);
}

void printIt(simulation_t *sim,FILE *fp) {
  int i;
  for(i=0; i<sim->boxes.nboxes; i++) {
    int j;
    real_t *x, *y, *z;
    real_t *fx, *fy, *fz;
    real_t *px, *py, *pz;
    domain_t *idm;
    int *id;
    idm  = sim->boxes.domains[i];
    x = idm->x;
    y = idm->y;
    z = idm->z;
    fx = idm->fx;
    fy = idm->fy;
    fz = idm->fz;
    px = idm->px;
    py = idm->py;
    pz = idm->pz;
    id = idm->id;
    for(j=0; j<idm->nAtoms[0]; j++,x++,y++,z++,px++,py++,pz++,fx++,fy++,fz++) {
      if ( id[j] < 10) {
	fprintf(fp,
		"%02d %+020.12e %+020.12e %+020.12e 1 %+020.12e %+020.12e %+020.12e F= %+020.12e %+020.12e %+020.12e\n",
		id[j]+1,*x,*y,*z,*px,*py,*pz,*fx,*fy,*fz);
      }
    }
  }
  return;
}
int cmdPPU(command_t *pcmd) {
  struct pmd_base_potential_t *pot;
  // variables;
  int rc;  
  simulation_t *sim;
  double ts, te;
  extern ljpotential_t *getLJPot();
  extern double timeNow();

  
  if(pcmd->doeam) {
    /* EAM potential */
    struct pmd_base_potential_t *p;
    pot = setEamPot(pcmd->potdir, pcmd->potname);
  }
  else {
    /* LJ potential */
    pot = (pmd_base_potential_t *) getLJPot();
  }
    
  if ( ! pot ) simulationAbort(-2,"Unable to read eam potential");
  printf("\n\n    File is __%s__\n\n", pcmd[0].filename);
  
  // get simulation info.  Simulation contains points to domain
  // point info and general metadata;
  sim = fromFileASCII(pcmd[0].filename, pot);

  if ( ! sim ) {
    printf("   \n\n"
	   "   ERROR:\n"
	   "     File %s does not exist\n"
	   "   \n\n", pcmd[0].filename);
    pot->destroy(&pot);
    exit(1);
  }

  printf("    Double Precision=%s\n",
	 (sizeof(real_t)==4?"false":"true")
	 );
  printf("    Domainheadersize=%d\n", (int) DOMAINHEADERSIZE);
  printf("    total atoms is: %d\n", sim->nAtoms);
  printf("    Max atoms in any box is: %d\n", getJmax());
  if (!sim) {
    fprintf(stderr,"Unable to get simulation");
    destroySimulation(&sim);
    pot->destroy((void *) &pot);
    exit(-1); 
  }	

  printf("    Computing the ppuforce\n");
  ts = timeNow();
  sim->force(sim);
  te = timeNow();
  //x  printIt(sim,stdout);
  printf("    %.3f s to compute the ppuforce (e = %30.20g)\n", te-ts,sim->e);
  if(pcmd[0].debug) {
    extern void debugOutput(int iv, char *file, simulation_t *s);
    debugOutput(pcmd[0].debug, "d_ppu.txt",sim);
  }

  {
    double eone;
    int itsa;
    extern double nTimeSteps(int n, simulation_t *s, real_t dt);

    for(itsa=0; itsa<200; itsa++) {
      ts = timeNow();
      sim->force(sim);
      printf(" %3d %30.20g computed in %.3fs\n", itsa, sim->e,(te-ts));
      eone = nTimeSteps(1,sim,1.0e-15);
      te = timeNow();
    }
    sim->force(sim);
    printf(" %3d %30.20g computed in %.3fs\n", itsa, sim->e,(te-ts));
    //    printIt(sim,stdout);
  }
  destroySimulation(&sim); 	

  return 0;
}


void debugOutput(int iv, char *filename, simulation_t *sim) {
  FILE *fpn;
  int i, j;
  domain_t *d;
  real_t maxdiff[3],a,b,c;
  printf("    Writing forces to file .. ");fflush(stdout);
  fpn = fopen(filename,"w");
  fprintf(fpn,"id ea fx fy fz\n");
  j = 0;
  maxdiff[0] = maxdiff[1] = maxdiff[2] = 0.0;

  for(i=0; i<sim->boxes.nboxes; i++) {
    d = sim->boxes.domains[i];
    if ( !d ) continue;
    if ( ! d->nAtoms) continue;
    if(iv > 1 ) fprintf(fpn,"Box %10d with %d atoms\n",i,*(d->nAtoms));
    for(j=0; j<*(d->nAtoms); j++) {
      fprintf(fpn," %10d % 19.12e (% 19.12e,% 19.12e,% 19.12e)\n",
	      d->id[j],d->ea[j],d->fx[j],d->fy[j],d->fz[j]);
    }
    if ( i%(sim->boxes.nboxes/10) == 9){
      printf(" .");
      fflush(stdout);
    }
  }
  printf(" Done\n");
  fclose(fpn);
    
}

