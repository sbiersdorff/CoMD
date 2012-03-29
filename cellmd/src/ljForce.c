/**
 * The file for force subroutines.
 * Initially we will only have lennard jones potential here.
 *
 **/
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include "pmd.h"

extern int LJ(void *s);

void ljDestroy(void **inppot) {
  ljpotential_t *pot;
  if ( ! inppot ) return;
  pot = *(ljpotential_t **)inppot;
  if ( ! pot ) return;
  suAlignedFree(pot);
  *inppot = NULL;
  return;
}


ljpotential_t *getLJPot() {
  ljpotential_t *pot;
  pot = suAlignedMalloc(sizeof(ljpotential_t));
  pot->force = LJ;
  pot->destroy = ljDestroy;
  pot->sigma = 3.41;
  pot->epsilon = 0.012606;
  pot->cutoff = 5.0*pot->sigma;
  pot->mass = 1.0;
  return pot;
}


int LJ(void *inS) {
  /**
   * calculates forces for the 12-6 lennard jones potential **/
  /**
   * Notes on LJ:
   *
   * http://en.wikipedia.org/wiki/Lennard_Jones_potential
   *
   * LJ is a simple potential of the form:
   *
   * e_lj(r) = 4*epsilon*((sigma/r)^12 - (sigma/r)^6)
   * F(r) = 4*epsilon*(12*sigma^12/r^13 - 6*sigma^6/r^7)
   *
   * epsilon and sigma are the adjustable parameters in the potential.
   *    epsilon = well depth
   *    sigma   = hard sphere diameter
   *
   * You can also adjust the 12 & 6, but few people do.
   *
   * Some typical values for epsilon and sigma:
   *
   *   material            epsilon            sigma   
   *     N                 36.20 K             3.299 A
   *     O                 44.06 K             2.956 A
   **/
  simulation_t *s = (simulation_t *) inS;
  int ii, ij, i, j, jTmp;
  real_t r2, r, r6, r7;
  real_t s6;
  real_t f;
  int nIBox;
  int nJBox;
  int *nbrBoxes;
  domain_t *iDomain;
  domain_t *jDomain;
  int iBox, jBox;
  real_t dx, dy, dz;
  ljpotential_t *pot;
  real_t sigma = -1;
  real_t epsilon = -1;
  real_t rcut=-1;
  real_t r2cut = -1;
  double etot;
  
  pot = (ljpotential_t *) s->pot;
  sigma = pot->sigma;
  epsilon = pot->epsilon;
  rcut = pot->cutoff;
  r2cut = rcut*rcut;

  /**
   * point to energy and force and zero them out **/
  etot = 0.0;
  s->e = (real_t) 0.0;

  for(i=0; i<s->boxes.nboxes; i++) {
    zeroForcesInBox(s,i);
  }
  s6 = sigma*sigma*sigma*sigma*sigma*sigma;


  /**
   * The more I read this, the more it seems that I should transfer
   * atom position / force data to the box structure so that the
   * atoms in a box will be contiguous in memory.
   *
   * Right now the atoms are part of the overall simulation structure.
   * This means that potentially two atoms within a box may be widely
   * separated in main memory.  
   **/
  
  for(iBox=0; iBox<s->boxes.nboxes; iBox++) /* loop over all boxes in system */
    {
      iDomain = s->boxes.domains[iBox];
      nIBox = *(iDomain->nAtoms);
      if ( ! nIBox) continue;

      nbrBoxes = getNeighborBoxes(s,iBox);

      for(jTmp=0; jTmp<NUMNEIGHBORS; jTmp++) /* loop over neighbor boxes */
	{
	  jBox = nbrBoxes[jTmp];
	  if(jBox<0) break; /* done with neighbor boxes */

	  jDomain = s->boxes.domains[jBox];
	  nJBox = *(jDomain->nAtoms);
	  if ( ! nJBox) continue;
	  
	  for(ii=0; ii<nIBox; ii++) /* loop over atoms in iBox */
	    {
	      i = iDomain->id[ii];  /* the ij-th atom in iBox */
	      for(ij=0; ij<nJBox; ij++)  /* loop over atoms in iBox */
		{
		  j = jDomain->id[ij];  /* the ij-th atom in iBox */
		  if ( j <= i ) continue;
	  
		  /** dist2(r2,iDomain,jDomain,ii,ij,dx,dy,dz);**/

		  dx = iDomain->x[ii]-jDomain->x[ij];
		  dy = iDomain->y[ii]-jDomain->y[ij];
		  dz = iDomain->z[ii]-jDomain->z[ij];
		  r2 = dx*dx + dy*dy + dz*dz;
		  if ( r2 > r2cut) continue;
	  
		  /**
		   * Important note:
		   *
		   * from this point on r actually refers to 1.0/r
		   *
		   **/
		  r2 =(real_t)1.0/r2;
		  r = sqrt(r2);
		  r6 = (r2*r2*r2);
		  r7 = r6*r;
		  etot += r6*(s6*r6 - 1.0); 
		  f = 4.0*epsilon*s6*r7*(12.0*s6*r6 - 6.0);
		  iDomain->fx[ii] += dx*f*r;
		  iDomain->fy[ii] += dy*f*r;
		  iDomain->fz[ii] += dz*f*r;

		  jDomain->fx[ij] -= dx*f*r;
		  jDomain->fy[ij] -= dy*f*r;
		  jDomain->fz[ij] -= dz*f*r;

		} /* loop over atoms in iBox */
	    } /* loop over atoms in iBox */
	} /* loop over neighbor boxes */

    } /* loop over all boxes in system */

  etot = etot*4.0*epsilon*s6;
  s->e = (real_t) etot;
  return 0;
}
      
