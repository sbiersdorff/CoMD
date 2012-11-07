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
 * An interface for reading afv style potential files
 *
 * Written by Sriram Swaminarayan 9/11/2006
 **/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>

#include "memUtils.h"
#include "pmd.h"
#include "cheby.h"
#include "eam.h"
#include "utility.h"

/**
 * An endian swapping utility
 * Currently undefined
 **/
#define endianSwapIfNeeded(p, nElements, elementSize, swapFlag) {	\
    int  i, j;								\
    char *ptr;								\
    char *c;								\
    char *d;								\
    char e;								\
    if ( swapFlag ) {							\
      ptr = (char *) (p);						\
      for(i=0; i<(nElements); i++,ptr+=(elementSize)) {			\
	c = ptr;							\
	d = c + (elementSize)-1;					\
	for(j=0; j<(elementSize)/2; j++,c++,d--) {			\
	  e = *c;							\
	  *c = *d;							\
	  *d = e;							\
	}								\
      }									\
    }									\
}


extern int PERIODIC;

struct eampotential_t *myPot = NULL;

struct pmd_base_potential_t *setEamPotFromPotential(struct pmd_base_potential_t *inPot) {
  myPot = (struct eampotential_t *) inPot;
  return (struct pmd_base_potential_t *) myPot;
}
  
struct pmd_base_potential_t *setEamPot(char *dir, char *file) {
  
  if(myPot) eamDestroy((void **) &myPot);
  myPot = eamReadASCII(dir,file);
  myPot->destroy=eamDestroy;
  if ( ! myPot) simulationAbort(-2,(char *) "Unable to read potential file\n");
  return (struct pmd_base_potential_t *) myPot;
}

struct eampotential_t *getEamPot() {
  if ( ! myPot) setEamPot((char *) "pots",(char *) "ag");
  myPot->destroy=eamDestroy;
  return myPot;
}

static void destroyPotentialArray(struct potentialarray_t **a, int doubleFlag) {
  if ( ! a ) return;
  if ( ! *a ) return;
  if ( (*a)->values) {
    (*a)->values--;
    (*a)->values-=doubleFlag;
    suAlignedFree((*a)->values);
  }
  suAlignedFree(*a);
  *a = NULL;
  return;
}

static struct potentialarray_t *allocPotentialArray(int n, real_t x0, real_t xn, real_t invDx) {
  struct potentialarray_t *a;
  int is;
  is = (sizeof(struct potentialarray_t)+15 ) & ~15;
  a = (struct potentialarray_t*)suAlignedCalloc(is);
  
  if ( ! a ) return NULL;

  // Always assume double precision arrays!
  is = ((n+3)*sizeof(real_t)+15 ) & ~15;
  a->values = (real_t*)suAlignedCalloc(is);

  if ( ! a->values) {
    suAlignedFree(a);
    return NULL;
  }
  a->values++; 
  a->n = n;
  a->invDx = invDx;
  a->xn = xn;
  a->x0 = x0 + (xn-x0)/(double)n;
  return a;

}


static struct potentialarray_t *getPotentialArrayFromBinaryFile(char *file) {
  struct potentialarray_t *retArray;
  FILE *fp;
  int n;
  int recSize;
  real_t *vals;
  double x0, xn, invDx;
  double *inData; 
  char swapFlag = 0;
  int eightByteHeader = 0;
  int itmp;
  int iflo, ifhi;
  
  fp = fopen(file,"rb");
  if ( ! fp ) {
    return NULL;
  }

  /* read record header and decide swap or not */
  fread(&recSize,sizeof(int),1,fp);
  swapFlag = (recSize > 4096);
  endianSwapIfNeeded(&recSize,1,sizeof(int),swapFlag);

  fread(&n,sizeof(int),1,fp);
  endianSwapIfNeeded(&n,1,sizeof(int),swapFlag);

  
  fread(&x0,sizeof(double),1,fp);
  endianSwapIfNeeded(&x0,1,sizeof(double),swapFlag);
  
  fread(&xn,sizeof(double),1,fp);
  endianSwapIfNeeded(&xn,1,sizeof(double),swapFlag);
  
  fread(&invDx,sizeof(double),1,fp);
  endianSwapIfNeeded(&invDx,1,sizeof(double),swapFlag);

  /* discard two integers */
  fread(&iflo,sizeof(int),1,fp);
  endianSwapIfNeeded(&iflo,1,sizeof(int),swapFlag);

  fread(&ifhi,sizeof(int),1,fp);
  endianSwapIfNeeded(&ifhi,1,sizeof(int),swapFlag);

  retArray = allocPotentialArray(n,x0,xn,invDx);
  if ( ! retArray ) {
        fclose(fp);
	return NULL;
  }

  /* read record trailer */
  fread(&recSize,sizeof(int),1,fp);
  endianSwapIfNeeded(&recSize,1,sizeof(int),swapFlag);
  
  /* read next record header and confirm size */
  fread(&recSize,sizeof(int),1,fp);
  endianSwapIfNeeded(&recSize,1,sizeof(int),swapFlag);

  if ( recSize != n*sizeof(double) ) {
    fclose(fp);
    destroyPotentialArray(&retArray,0);
    printf("sizes are: %d,%d\n",recSize, (int)(n*sizeof(double)));
    simulationAbort(-500,(char *) "Size mismatch error reading binary potential file.");
  }

  /* allocate space and read in potential data */
  inData = (double*)suAlignedCalloc(n*sizeof(double));
  if ( ! inData) {
    fclose(fp);
    destroyPotentialArray(&retArray,0);
    return NULL;
  }
  fread(inData,sizeof(double),n,fp);
  endianSwapIfNeeded(inData,n,sizeof(double),swapFlag);
  vals = retArray->values;
  for(n=0; n<retArray->n; n++, vals++) *vals = (real_t) inData[n];
  suAlignedFree(inData);
  {
    /* array values for robust interpolation */
    real_t lo, hi;
    if(iflo==0) lo=0.0;
    else lo = retArray->values[0];

    if(ifhi==0) hi=0.0;
    else hi = retArray->values[n-1];

    n = retArray->n;
    retArray->values[-1] = lo;
    retArray->values[n] = hi;
    retArray->values[n+1] = hi;

  }
  fclose(fp);
  return retArray;
}


static struct potentialarray_t *getPotentialArrayFromFile(char *file) {
  struct potentialarray_t *retArray;
  char tmp[4096];
  FILE *fp;
  int n;
  real_t *vals;
  real_t x0, xn, invDx;

  /* check on binary file */
  if(file[strlen(file)-1] != 't') return getPotentialArrayFromBinaryFile(file);

  fp = fopen(file,"r");
  if ( ! fp ) {
    return NULL;
  }

  /**
   * read first line **/
  fgets(tmp,sizeof(tmp),fp);
  sscanf(tmp,"%10d " FMT1 " " FMT1 " " FMT1 ,
	 &n, &x0, &xn, &invDx);

  retArray = allocPotentialArray(n,x0,xn,invDx);
  if ( ! retArray ) return NULL;

  vals = retArray->values;
  for(n=0; n<retArray->n; n++, vals++) fscanf(fp,FMT1,vals);
  {
    /* array values for robust interpolation */
    n = retArray->n;
    retArray->values[-1] = retArray->values[0];
    retArray->values[n] = retArray->values[n-1];
    retArray->values[n+1] = retArray->values[n-1];
  }
  fclose(fp);
  return retArray;
}

static double getMassFromFile(char *file) {
  double mass;
  char tmp[4096];
  FILE *fp;
  int n;
  fp = fopen(file,"r");
  if ( ! fp ) {
    return -1.0;
  }

  /**
   * read first line **/
  fgets(tmp,sizeof(tmp),fp);
  /**
   * get mass from second line **/
  fgets(tmp,sizeof(tmp),fp);
  sscanf(tmp,"%lf ",&mass);
  return mass;
}

static double getLatFromFile(char *file) {
  double lat;
  char tmp[4096];
  FILE *fp;
  int n;
  fp = fopen(file,"r");
  if ( ! fp ) {
    return -1.0;
  }

  /**
   * read first two lines **/
  fgets(tmp,sizeof(tmp),fp);
  fgets(tmp,sizeof(tmp),fp);
  /**
   * get lat from third line **/
  fgets(tmp,sizeof(tmp),fp);
  sscanf(tmp,"%lf ",&lat);
  return lat;
}

void eamDestroy(void **inppot) {
  pmd_base_potential_t **pPot = (pmd_base_potential_t **) inppot;
  eampotential_t *pot;
  if ( ! pPot ) return;
  pot = *(eampotential_t **)pPot;
  if ( pot == myPot) myPot = NULL;
  if ( ! pot ) return;
  if(pot->phi) destroyPotentialArray(&(pot->phi),0);
  if(pot->rho) destroyPotentialArray(&(pot->rho),0);
  if(pot->f) destroyPotentialArray(&(pot->f),0);
  suAlignedFree(pot);
  *pPot = NULL;
  myPot = NULL;
  return;
}

eampotential_t *eamReadASCII(char *dir, char *potname) {
  /**
   * reads potential potname from directory dir.
   * returns a poitner to an eampotential_t struct.
   **/
  eampotential_t *retPot;
  char tmp[4096];
  int is;
  is = (sizeof(struct eampotential_t)+15 ) & ~15;
  retPot = (eampotential_t*)suAlignedCalloc(is);
  if ( ! retPot ) return NULL;

  /**
   * read the phi component **/
  sprintf(tmp,"%s/%s.phi",dir,potname);
  retPot->phi = getPotentialArrayFromFile(tmp);

  /**
   * read the rho component **/
  sprintf(tmp,"%s/%s.rho",dir,potname);
  retPot->rho = getPotentialArrayFromFile(tmp);

  /**
   * read the F component **/
  sprintf(tmp,"%s/%s.f",dir,potname);
  retPot->f = getPotentialArrayFromFile(tmp);

  sprintf(tmp,"%s/%s.doc",dir,potname);
  retPot->mass = (real_t) getMassFromFile(tmp);
  retPot->lat = (real_t) getLatFromFile(tmp);
  printf("lattice constant: %e\n", retPot->lat);

  if ( (retPot->mass < 0.0 ) || (! (retPot->phi && retPot->rho && retPot->f )) ) {
     printf("\n\n"
	    "    ****  Unable to open potential file %s.  **** \n\n"
	    "    Did you untar pots.tgz (tar zxvf pots.tgz)?"
	    "\n\n"
	    ,potname);
     eamDestroy((void **) &retPot);
     
    return NULL;
  }

  /**
   * set the cutoff from the phi component **/
  retPot->cutoff = retPot->phi->xn;

  retPot->force = eamForce;

  return retPot;
  
}


static inline void eamInterpolateDeriv(struct potentialarray_t *a, real_t r, int iType, int jType, real_t *value1, real_t *f1) {
  /**
   *
   * This routine will not crash if r is out of range.
   *
   * if ( r < a->x0) r = a->x0;
   * if ( r > a->xn)   r = a->xn;
   **/
   
  int i1;
  real_t gi, gi1;

  if ( r<a->x0) r = a->x0;
  else if (r>a->xn) r = a->xn;
  
  r = (r-a->x0)*(a->invDx) ;
  i1 = (int)floor(r);

  /* reset r to fractional distance */
  r = r - floor(r);

  gi  = a->values[i1+1] - a->values[i1-1];
  gi1 = a->values[i1+2] - a->values[i1];


  *value1 = a->values[i1] + 0.5*r*(
				   r*(a->values[i1+1]+ a->values[i1-1] -2.0*a->values[i1]) +
				   gi
				   );
  if(i1<=0) *f1 = 0.0;
  else *f1 = 0.5*(gi + r*(gi1-gi))*a->invDx;

  return;
}



int eamForce(void *inS) {
  /**
   * calculates forces for the EAM potential **/
  simflat_t *s = (simflat_t *) inS;
  eampotential_t *pot = NULL;
  eam_cheby_t *ch_pot = s->ch_pot;
  int ii, ij, i, j, jTmp;
  real_t r2,r;
  int nIBox;
  int nJBox;
  int *nbrBoxes;
  int iBox, jBox;
  real_t rcut2;
  real_t rhoTmp;
  real_t phiTmp;
  real_t dPhi, dRho;
  real_t bndX, bndY, bndZ;
  int ncc = 0;
  double etot;

  pot = (eampotential_t *) s->pot;
  rcut2 = pot->cutoff*pot->cutoff;


  /* zero forces / energy / rho /rhoprime */
  etot = 0.0;
  memset(s->f,0,s->nboxes*MAXATOMS*sizeof(real4));
  memset(s->fi,0,s->nboxes*MAXATOMS*sizeof(real_t));
  memset(s->rho,0,s->nboxes*MAXATOMS*sizeof(real_t));
  
  // virial stress computation added here
  s-> stress = 0.0;

/*
#if (USE_CHEBY) 
    //Test to see if the data is there
    printf("Chebychev coefficients:\n");
    fflush(stdout);
    printf("%d, %d, %d\n", 
    ch_pot->phi->n,
    ch_pot->rho->n,
    ch_pot->f->n);
    fflush(stdout);
#endif
    */

  for(iBox=0; iBox<s->nboxes; iBox++) { /* loop over all boxes in system */ 

    nIBox = s->natoms[iBox];
    nbrBoxes = getNeighborBoxes(s,iBox);
      
    for(jTmp=0; jTmp<nbrBoxes[-1]; jTmp++) { /* loop over neighbor boxes */
      int ioff;
      real3 drbox;
      jBox = nbrBoxes[jTmp];
      if(jBox<0) break;
      if (jBox < iBox ) continue;

      for(j=0; j<3; j++) {
	drbox[j] = s->dcenter[iBox][j]-s->dcenter[jBox][j];
	if(PERIODIC) {
	  if(drbox[j]<-0.5*s->bounds[j]) drbox[j] += s->bounds[j];
	  else if (drbox[j] > 0.5*s->bounds[j] ) drbox[j] -= s->bounds[j];
	}
      }

      nJBox = s->natoms[jBox];
      for(ioff=MAXATOMS*iBox,ii=0; ii<nIBox; ii++,ioff++) { /* loop over atoms in iBox */
	int joff;
	for(joff=MAXATOMS*jBox,ij=0; ij<nJBox; ij++,joff++) {  /* loop over atoms in iBox */
	  int k;
	  real3 dr;

	  if ( (iBox==jBox) &&(ij <= ii) ) continue;

	  /** dist2(r2,atomsInIBox,jDomain,ii,ij,dx,dy,dz);**/
	  r2 = 0.0;
	  for(k=0; k<3; k++) {
	    dr[k]=drbox[k]+s->r[ioff][k]-s->r[joff][k];
	    r2+=dr[k]*dr[k];
	  }
	  if(r2>rcut2) continue;

	  r = sqrt(r2);

#if (USE_CHEBY)
	phiTmp = eamCheb(s->ch_pot->phi, r);
	dPhi = eamCheb(s->ch_pot->dphi, r);

	rhoTmp = eamCheb(s->ch_pot->rho, r);
	dRho = eamCheb(s->ch_pot->drho, r);
#else 
	  eamInterpolateDeriv(pot->phi,r,0,0,&phiTmp,&dPhi);
	  eamInterpolateDeriv(pot->rho,r,0,0,&rhoTmp,&dRho);
#endif

	  for(k=0; k<3; k++) {
	    s->f[ioff][k] += dPhi*dr[k]/r;
	    s->f[joff][k] -= dPhi*dr[k]/r;
	  }
	  /* update energy terms */
	  etot += (double) phiTmp;
	  s->f[ioff][3] += phiTmp/2.0;
	  s->f[joff][3] += phiTmp/2.0;

	  /* update rho terms */
	  s->rho[ioff] += rhoTmp;
	  s->rho[joff] += rhoTmp;
	} /* loop over atoms in jBox */
      } /* loop over atoms in iBox */
    } /* loop over neighbor boxes */
  } /* loop over all boxes in system */

/*
  for(i=0; i<s->nboxes; i++) {
    int j;
    int ioff;
    int *id;

    for(ioff=i*MAXATOMS,j=0; j<s->natoms[i]; j++,ioff++) {
      if ( s->id[ioff] < 10) {
	fprintf(stdout,
		"%02d %02d %+020.12e %+020.12e %+020.12e 1 %+020.12e %+020.12e %+020.12e F= %+020.12e %+020.12e %+020.12e\n",
                i,
		s->id[ioff]+1,
		s->r[ioff][0],s->r[ioff][1],s->r[ioff][2],
		s->rho[ioff],s->fi[ioff],s->p[ioff][2],
		s->f[ioff][0],s->f[ioff][1],s->f[ioff][2]
		);
      }
    }
  }
  */

  /**
   * Now loop again and add in the F(rho(ij))
   **/

  for(iBox=0; iBox<s->nboxes; iBox++) { /* loop over all boxes in system */
    int ioff;
    nIBox =  s->natoms[iBox];
      
    for(ioff=MAXATOMS*iBox,ii=0; ii<nIBox; ii++,ioff++) { /* loop over atoms in iBox */
      real_t fi, fiprime;
#if (USE_CHEBY)
	fi = eamCheb(s->ch_pot->f, s->rho[ioff]);
	fiprime = eamCheb(s->ch_pot->df, s->rho[ioff]);
#else
      eamInterpolateDeriv(pot->f,s->rho[ioff],0,0,&fi,&fiprime);
#endif
      s->fi[ioff] = fiprime; /* update rhoprime */
      /* update energy terms */
      etot += (double) fi; 
      s->f[ioff][3] += fi;
      //s->stress -= s->p[ioff][0]*s->p[ioff][0]/s->mass[ioff];
    }
  }
  for(iBox=0; iBox<s->nboxes; iBox++) { /* loop over all boxes in system */
    int ioff;
    nIBox =  s->natoms[iBox];

    nbrBoxes = getNeighborBoxes(s,iBox);
      
    for(ioff=MAXATOMS*iBox,ii=0; ii<nIBox; ii++,ioff++) { /* loop over atoms in iBox */
      for(jTmp=0; jTmp<NUMNEIGHBORS; jTmp++) { /* loop over neighbor boxes */
	real3 drbox;
	int joff;
	
	jBox = nbrBoxes[jTmp];
	if(jBox<0) break;

	if(jBox < iBox) continue;

	for(j=0; j<3; j++) {
	  drbox[j] = s->dcenter[iBox][j]-s->dcenter[jBox][j];
	  if(PERIODIC) {
	    if(drbox[j]<-0.5*s->bounds[j]) drbox[j] += s->bounds[j];
	    else if (drbox[j] > 0.5*s->bounds[j] ) drbox[j] -= s->bounds[j];
	  }
	}

	nJBox = s->natoms[jBox];

	for(joff=MAXATOMS*jBox,ij=0; ij<nJBox; ij++,joff++)  { /* loop over atoms in iBox */
	  real3 dr;
	  real_t rhoijprime;
	  int k;
	  if ((iBox==jBox) && (ij <= ii))  continue;
		  
	  /** dist2(r2,atomsInIBox,jDomain,ii,ij,dx,dy,dz);**/
	  r2 = 0.0;
	  for(k=0; k<3; k++) {
	    dr[k]=drbox[k]+s->r[ioff][k]-s->r[joff][k];
	    r2+=dr[k]*dr[k];
	  }
	  if(r2>=rcut2) continue;

	  r = sqrt(r2);

#if (USE_CHEBY)
	dRho = eamCheb(s->ch_pot->drho, r);
#else
	  eamInterpolateDeriv(pot->rho,r,0,0,&rhoTmp,&dRho);
#endif
	  rhoijprime = dRho;

	  for(k=0; k<3; k++) {
	    s->f[ioff][k] += (s->fi[ioff]+s->fi[joff])*rhoijprime*dr[k]/r;
	    s->f[joff][k] -= (s->fi[ioff]+s->fi[joff])*rhoijprime*dr[k]/r;
	  }

          s->stress += 2.0*(s->fi[ioff]+s->fi[joff])*rhoijprime*dr[0]/r*dr[0];
	} /* loop over atoms in jBox */
      } /* loop over atoms in iBox */
    } /* loop over neighbor boxes */
  } /* loop over all boxes in system */

  s->e = (real_t) etot;

  s->stress = s->stress/(s->bounds[0]*s->bounds[1]*s->bounds[2]);

  return 0;
}


/**
 * utility comparison routine **/
static void adiffpot(char *name,potentialarray_t *a, potentialarray_t *b) {
  int i;
  printf("---------------------------------------\n");
  printf("  comparison of %s\n", name);
  printf("    n = %4d   /  %4d\n", a->n, b->n);	 
  printf("   x0 = %10.2g  /  %10.2g\n", a->x0, b->x0);	 
  printf("   xn = %10.2g  /  %10.2g\n", a->xn, b->xn);	 
  printf("   dx = %10.2g  /  %10.2g\n", a->invDx, b->invDx);
  for(i=-1; i<a->n+2;i++) {
    if ( a->values[i] != b->values[i]) {
      printf("   v[%d] = %10.2g  /  %10.2g\n", i, a->values[i],b->values[i]);
    }
  }
  printf("---------------------------------------\n");
  return;
}

void eamComparePots(eampotential_t *a, eampotential_t *b) {
  adiffpot((char *) "phi", a->phi, b->phi);
  adiffpot((char *) "rho", a->rho, b->rho);
  adiffpot((char *) "f", a->f, b->f);
  return;
}
