/**
 * The file for timestep subroutines.
 * Initially we will only have verlet
 *
 **/
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include "pmd.h"

static void advanceVelocity(simulation_t *s, const int nStart, const int nTodo, real_t dt) {
  int iBox;
  domain_t *iDomain;
  const real_t cfac = (real_t)(
			       (double)0.529/
			       (double)2.418e-17
			       );
  dt = dt*cfac;

  for(iBox=nStart; iBox<nStart+nTodo; iBox++) {
    /* loop over my boxes in system */
    real_t *px, *py, *pz;
    real_t *fx, *fy, *fz;
    int ii,nIBox;
    iDomain = s->boxes.domains[iBox];
    nIBox = *(iDomain->nAtoms);
    if ( ! nIBox) continue;
    px = iDomain->px;
    py = iDomain->py;
    pz = iDomain->pz;
    fx = iDomain->fx;
    fy = iDomain->fy;
    fz = iDomain->fz;
    for(ii=0; ii<nIBox; ii++) {
      /* loop over atoms in iBox */
      px[ii] -= dt*fx[ii];
      py[ii] -= dt*fy[ii];
      pz[ii] -= dt*fz[ii];
    }
  }
}

static void advancePositions(simulation_t *s, const int nStart, const int nTodo, real_t dt) {
  int iBox;
  domain_t *iDomain;
  const real_t cfac = (real_t)(
			       (double)0.529/
			       (double)2.418e-17
			       );
  real_t rmass;  

  //  printf("==================mass=%f\n",s->pot->mass);
  rmass = (real_t)(((double)1822.83*(double)s->pot->mass));

  dt = dt*cfac/rmass;

  for(iBox=nStart; iBox<nStart+nTodo; iBox++) {
    /* loop over my boxes in system */
    real_t *px, *py, *pz;
    real_t *x, *y, *z;
    int ii,nIBox;
    real_t dpx, dpy, dpz;
    dpx = s->periodicBounds[3];
    dpy = s->periodicBounds[4];
    dpz = s->periodicBounds[5];
    iDomain = s->boxes.domains[iBox];
    nIBox = *(iDomain->nAtoms);
    if ( ! nIBox) continue;
    px = iDomain->px;
    py = iDomain->py;
    pz = iDomain->pz;
    x = iDomain->x;
    y = iDomain->y;
    z = iDomain->z;
    for(ii=0; ii<nIBox; ii++) {
      /* loop over atoms in iBox */
      x[ii] += dt*px[ii];
      y[ii] += dt*py[ii];
      z[ii] += dt*pz[ii];
      if ( x[ii] < 0.0 ) x[ii] += dpx;
      else if ( x[ii] >= dpx ) x[ii] -= dpx;
      if ( y[ii] < 0.0 ) y[ii] += dpy;
      else if ( y[ii] >= dpy ) y[ii] -= dpy;
      if ( z[ii] < 0.0 ) z[ii] += dpz;
      else if ( z[ii] >= dpz ) z[ii] -= dpz;
    }
  }
}


static void reBox(simulation_t *s) {
  extern int getBoxID(simulation_t *s, real_t x, real_t y, real_t z);
  extern void moveAtom(simulation_t *s, int iId, int iBox, int jBox);
  int iBox;
  for(iBox=s->boxes.nboxes-1; iBox>=0; iBox--) {
    int id;
    int n;
    int jBox;
    domain_t *idm;
    real_t *x, *y, *z;
    idm = s->boxes.domains[iBox];
    n = idm->nAtoms[0];
    x = idm->x; y = idm->y; z = idm->z;
    id = 0;
    while ( id < idm->nAtoms[0]) {
      jBox = getBoxID(s,x[id],y[id],z[id]);
      if(jBox == iBox) id++;
      else moveAtom(s,id,iBox,jBox);
    }
  }
  return;
}

double nTimeSteps(int n, simulation_t *s, real_t dt) {
  extern void printIt(simulation_t *sim,FILE *fp);
  int i;
  /* for now no repositioning */
  /* half step forward for positions */
  //  advancePositions(s,0,s->boxes.nboxes,(dt/2.0));

  for(i=0; i<n; i++) {
    advancePositions(s,0,s->boxes.nboxes,(dt/2.0));
    reBox(s); /* reposition particles if needed */
    s->force(s); /* compute force */
    advanceVelocity(s,0,s->boxes.nboxes,dt); /* advance velocity */
    advancePositions(s,0,s->boxes.nboxes,(dt/2.0));
    //    advancePositions(s,0,s->boxes.nboxes,dt); /* advance positions */
  }
  reBox(s); /* reposition particles if needed */
  /* half step backward for positions */
  return s->e;
}

      
