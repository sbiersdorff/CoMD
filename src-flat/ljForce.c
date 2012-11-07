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
extern int PERIODIC;
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
  /* return a new potential with LJ for copper */
  ljpotential_t *pot;
  pot = (ljpotential_t*)suAlignedMalloc(sizeof(ljpotential_t));
  pot->force = LJ;
  pot->destroy = ljDestroy;
  pot->sigma = 1.53;
  pot->epsilon = 0.0085;
  pot->cutoff = 3.0*pot->sigma;
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
  simflat_t *s = (simflat_t *) inS;
  int ii, ij, i, j, jTmp;
  real_t r2, r, r6, r7;
  real_t s6;
  real_t f, fr;
  int nIBox;
  int nJBox;
  int *nbrBoxes;
  int ibox;
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

  /* zero forces / energy */
  memset(s->f,0,s->nboxes*MAXATOMS*sizeof(real4));

  s6 = sigma*sigma*sigma*sigma*sigma*sigma;

  //virial stress added here
  s->stress = 0.0;


  for(ibox=0; ibox<s->nboxes; ibox++) { /* loop over all boxes in system */
    nIBox = s->natoms[ibox];
    if ( ! nIBox) continue;
    nbrBoxes = getNeighborBoxes(s,ibox);
    for(jTmp=0; jTmp<nbrBoxes[-1]; jTmp++) { /* loop over neighbor boxes */
      real3 drbox;
      int ioff;
      int k;
      int jbox;
      jbox = nbrBoxes[jTmp];
      if(jbox<0) break; /* done with neighbor boxes */
      for(j=0; j<3; j++) {
	drbox[j] = s->dcenter[ibox][j]-s->dcenter[jbox][j];
	if(PERIODIC) {
	  if(drbox[j]<-0.5*s->bounds[j]) drbox[j] += s->bounds[j];
	  else if (drbox[j] > 0.5*s->bounds[j] ) drbox[j] -= s->bounds[j];
	}
      }

      nJBox = s->natoms[jbox];
      if ( ! nJBox) continue;
	  
      for(ioff=ibox*MAXATOMS,ii=0; ii<nIBox; ii++,ioff++) { /* loop over atoms in ibox */
	int joff;
        s->stress -= s->p[ioff][0]*s->p[ioff][0]/s->mass[ioff];
	int i = s->id[ioff];  /* the ij-th atom in ibox */
	for(joff=MAXATOMS*jbox,ij=0; ij<nJBox; ij++,joff++) { /* loop over atoms in ibox */
	  int m;
	  real_t dr[3];
	  int j = s->id[joff];  /* the ij-th atom in ibox */
	  if ( j <= i ) continue;
	  r2 = 0.0;
	  for(m=0; m<3; m++) {
	    dr[m] = drbox[m]+s->r[ioff][m]-s->r[joff][m];
	    r2+=dr[m]*dr[m];
	  }
	    
	  if ( r2 > r2cut) continue;
	  
	  /**
	   * Important note:
	   *
	   * from this point on r actually refers to 1.0/r
	   *
	   **/
	  r2 =(real_t)1.0/r2;
	  //r = sqrt(r2);
	  r6 = (r2*r2*r2);
	  //r7 = r6*r;
	  s->f[ioff][3] += 0.5*r6*(s6*r6 - 1.0); 
	  s->f[joff][3] += 0.5*r6*(s6*r6 - 1.0); 
	  etot += r6*(s6*r6 - 1.0); 
	  //f = 4.0*epsilon*s6*r7*(12.0*s6*r6 - 6.0);
          // different formulation to avoid sqrt computation
	  fr = 4.0*epsilon*s6*r6*r2*(12.0*s6*r6 - 6.0);
	  for(m=0; m<3; m++) {
	    s->f[ioff][m] += dr[m]*fr;
	    s->f[joff][m] -= dr[m]*fr;
	  }
          s->stress += 2.0*fr*dr[0]*dr[0];
	} /* loop over atoms in jbox */
      } /* loop over atoms in ibox */
    } /* loop over neighbor boxes */

  } /* loop over all boxes in system */

  etot = etot*4.0*epsilon*s6;
  s->e = (real_t) etot;

  // renormalize stress
  s->stress = s->stress/(s->bounds[0]*s->bounds[1]*s->bounds[2]);

  return 0;
}
      
