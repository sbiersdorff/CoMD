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

#include "eamTypes.h"

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

extern void simulationAbort(int ineCode, char *inmsg);


int PERIODIC = 1;
#ifdef INCLUDE_R_IN_DERIVATIVE

#if ( ( EAMLINEAR != 1 ) && defined(INCLUDE_R_IN_DERIVATIVE) )
#error You cannot have cubic interpolation and r included in derivative
#endif

#endif

#ifndef CELLCODE
#define CELLCODE 1
#endif

int icalls;
int icount;
struct eampotential_t *myPot = NULL;
extern int eamForce(void *s);


struct pmd_base_potential_t *setEamPotFromPotential(struct pmd_base_potential_t *inPot) {
  myPot = (struct eampotential_t *) inPot;
  return (struct pmd_base_potential_t *) myPot;
}
  
struct pmd_base_potential_t *setEamPot(char *dir, char *file) {
  if(myPot) eamDestroy((void **) &myPot);
  myPot = eamReadASCII(dir,file);
  if ( ! myPot) simulationAbort(-2,"Unable to read potential file\n");
  return (struct pmd_base_potential_t *) myPot;
}
struct eampotential_t *getEamPot() {
  if ( ! myPot) setEamPot("pots","ag");
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
  a = suAlignedCalloc(is);
  
  if ( ! a ) return NULL;

  // Always assume double precision arrays!
  is = ((n+3)*sizeof(real_t)+15 ) & ~15;
  a->values = suAlignedCalloc(is);

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

static struct potentialarray_t *eamCombinePhiAndRho(struct potentialarray_t *phi, potentialarray_t *rho) {
  /**
   * combines the phi and rho into an interleaved
   * single potential array
   **/
  struct potentialarray_t *retPot;
  int i;
  int is;
  double eps = 1.0e-14;

  if ( (abs(phi->n - rho->n) > 0)   ||
       (fabs(phi->xn - rho->xn) > eps)   ||
       (fabs(phi->x0 - rho->x0) > eps)   ||
       (fabs(phi->invDx - rho->invDx) > eps)
       ) {
	
    printf("%g %g %g %d\n", phi->x0-rho->x0, phi->xn-rho->xn, phi->invDx-rho->invDx, phi->n-rho->n);
    simulationAbort(-600,"phi and rho must have identical representations in potential");
  }

  /* allocate space for the structure */
  is = (sizeof(struct potentialarray_t)+15 ) & ~15;
  retPot = suAlignedCalloc(is);
  if ( ! retPot ) return NULL;
       
  
  /** shallow copy **/
  memcpy(retPot,phi,sizeof(struct potentialarray_t));

  /** resetting the values **/
  is = ((6+2*phi->n)*sizeof(real_t)+15 ) & ~15;
  retPot->values = suAlignedCalloc(is);
  if ( ! retPot->values) {
    suAlignedFree(retPot);
    return NULL;
  }
  retPot->values += 2;
  for(i=-1; i<phi->n+2; i++) {
    retPot->values[2*i] = phi->values[i];
    retPot->values[2*i+1] = rho->values[i];
  }
  return retPot;
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
    simulationAbort(-500,"Size mismatch error reading binary potential file.");
  }

  /* allocate space and read in potential data */
  inData = suAlignedCalloc(n*sizeof(double));
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
  printf("==================mass=%f\n",mass);
  return mass;
}


#if 1
struct pmd_base_potential_t *createEAMPotentialFromData
                             (int ftnIndexing,
			      int n_phi,double x0_phi, double xN_phi, double invDx_phi, double *values_phi,
			      int n_rho,double x0_rho, double xN_rho, double invDx_rho, double *values_rho,
			      int n_f,  double x0_f,   double xN_f,   double invDx_f,   double *values_f ) {
  /**
   * creates an eam potential in cellMD form from standard arrays **/
  struct eampotential_t *pot;
  extern int        suAbort(int status, char *message);
  int is;

  /**
     * Since I diddle my values, I need to reset them
     * You should not need to do this for your potentials
     **/

  if ( ftnIndexing == 2) {
    n_phi++;
    n_rho++;
    ftnIndexing = 0;
  }

  if(! ftnIndexing) {
    x0_phi = (x0_phi*n_phi-xN_phi)/(double)(n_phi-1.0);
    x0_rho = (x0_rho*n_rho-xN_rho)/(double)(n_rho-1.0);
    x0_f   = (x0_f*n_f-xN_f)/(double)(n_f-1.0);
  }

  /** allocate main struct **/
  is = (sizeof(struct eampotential_t)+15 ) & ~15;
  pot = suAlignedCalloc(is);
  if ( ! pot ) suAbort(-1,"out of memory 1 in createEAMPotential");


  /** do phi **/
  pot->phi = allocPotentialArray(n_phi,x0_phi,xN_phi,invDx_phi);
  if ( ! pot->phi ) {
    suAbort(-1,"out of memory 2 in createEAMPotential");
    suAlignedFree(pot);
  }
  memcpy(pot->phi->values,values_phi,n_phi*sizeof(double));
  pot->phi->values[-1] = values_phi[0];
  pot->phi->values[n_phi] = pot->phi->values[n_phi+1] = 0.0;

  /** do rho **/
  pot->rho = allocPotentialArray(n_rho,x0_rho,xN_rho,invDx_rho);
  if ( ! pot->rho ) {
    destroyPotentialArray(&(pot->phi),1);
    suAbort(-1,"out of memory 2 in createEAMPotential");
    suAlignedFree(pot);
  }
  memcpy(pot->rho->values,values_rho,n_rho*sizeof(double));
  pot->rho->values[-1] = values_rho[0];
  pot->rho->values[n_rho] = pot->rho->values[n_rho+1] = 0.0;

  /** do f **/
  pot->f = allocPotentialArray(n_f,x0_f,xN_f,invDx_f);
  if ( ! pot->f ) {
    suAbort(-1,"out of memory 2 in createEAMPotential");
    destroyPotentialArray(&(pot->phi),1);
    destroyPotentialArray(&(pot->rho),1);
    suAlignedFree(pot);
  }
  memcpy(pot->f->values,values_f,n_f*sizeof(double));
  pot->f->values[-1] = 0.0;
  pot->f->values[n_f] = pot->f->values[n_f+1] = values_f[n_f-1];


  pot->cutoff = xN_phi;
  pot->phirho = eamCombinePhiAndRho(pot->phi, pot->rho);
  pot->force = eamForce;
  
  return (pmd_base_potential_t *) pot;
}
  
#endif
void eamDestroy(void **inppot) {
  pmd_base_potential_t **pPot = (pmd_base_potential_t **) inppot;
  eampotential_t *pot;
  if ( ! pPot ) return;
  pot = *(eampotential_t **)pPot;
  if ( pot == myPot) myPot = NULL;
  if ( ! pot ) return;
  if(pot->phirho) destroyPotentialArray(&(pot->phirho),1);
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
  retPot = suAlignedCalloc(is);
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

  if ( (retPot->mass < 0.0 ) || (! (retPot->phi && retPot->rho && retPot->f )) ) {
     printf("\n\n"
	    "    ****  Unable to open potential file %s.  **** \n\n"
	    "    Did you untar pots.tgz (tar zxvf pots.tgz)?"
	    "\n\n"
	    ,potname);
     eamDestroy((void **) &retPot);
     
    return NULL;
  }

  retPot->phirho = eamCombinePhiAndRho(retPot->phi, retPot->rho);

  /**
   * set the cutoff from the phi component **/
  retPot->cutoff = retPot->phi->xn;

  retPot->force = eamForce;
  return retPot;
  
}

real_t eamInterpolate(struct potentialarray_t *a, real_t r, int itype) {
  /**
   *
   * This routine will not crash.
   *
   * if ( r < a->x0) r = a->x0;
   * if ( r > a->xn)   r = a->xn;
   **/
   
  int i1, i2;
  real_t ret;
  icalls++;
  ret = 0.0;
  r = (r<a->x0?a->x0:r); /* r = min(r,a->x0)  */
  r = (r>a->xn?a->xn:r);     /* r = max(r,a->xn)    */

  
  r = (r-a->x0)*(a->invDx);

  i1 = (int)(floor(r));
  r = r-floor(r);
  i2 = i1+1; /* no checks on i2 because we are guaranteed existence */


  ret = a->values[i1] + 0.5*r*(
			       (a->values[i1+1]+ a->values[i1-1] -2.0*a->values[i1])*r +
			       a->values[i1+1] - a->values[i1-1]
			       );
  return ret;
}
void eamInterpolate2(struct potentialarray_t *a, real_t r, int itype, real_t *p1, real_t *p2) {
  /**
   *
   * This routine will not crash.
   *
   * if ( r < a->x0) r = a->x0;
   * if ( r > a->xn)   r = a->xn;
   **/
   
  int i1, i2;
  icalls+=2;
  
  r = (r<a->x0?a->x0:r); /* r = min(r,a->x0)  */
  r = (r>a->xn?a->xn:r);     /* r = max(r,a->xn)    */

  
  r = (r-a->x0)*(a->invDx);

  i1 = 2*(int)(floor(r));
  i2 = i1+2; /* no checks on i2 because we are guaranteed existence */
  r = r-floor(r);
  *p1 = a->values[i1] + r*(a->values[i2] - a->values[i1]);
  *p2 = a->values[i1+1] + r*(a->values[i2+1] - a->values[i1+1]);

  *p1 = a->values[i1] + 0.5*r*(
			       (a->values[i1+2]+ a->values[i1-2] -2.0*a->values[i1])*r +
			       a->values[i1+2] - a->values[i1-2]
			       );
  i1++;
  *p1 = a->values[i1] + 0.5*r*(
			       (a->values[i1+2]+ a->values[i1-2] -2.0*a->values[i1])*r +
			       a->values[i1+2] - a->values[i1-2]
			       );



  return;
}

void eamInterpolateDeriv2(struct potentialarray_t *a, real_t r,
			  int iType, int jType,
			  real_t *value1, real_t *value2,
			  real_t *f1, real_t *f2) {
  /**
   *
   * This routine will not crash if r is out of range.
   *
   * if ( r < a->x0) r = a->x0;
   * if ( r > a->xn)   r = a->xn;
   **/
   
  int i1, i1_minus_1, i1_plus_1,i1_plus_2;
  real_t gi, gi1;
  
  icalls+=2;
  r = (r<a->x0?a->x0:r); /* r = min(r,a->x0)  */
  r = (r>a->xn?a->xn:r);     /* r = max(r,a->xn)    */
  
  r = (r-a->x0)*(a->invDx) ;
  i1         = 2* (int)floor(r);
  i1_plus_1  = i1+2;
  i1_plus_2  = i1+4;
  i1_minus_1 = i1-2;

  /* reset r to fractional distance */
  r = r - floor(r);

#if 0  
  *value1 = a->values[i1] + r*(a->values[i1_plus_1] - a->values[i1]);
#endif
  *value1 = a->values[i1] + 0.5*r*(
			       (a->values[i1+2]+ a->values[i1-2] -2.0*a->values[i1])*r +
			       a->values[i1+2] - a->values[i1-2]
			       );
  gi  = a->values[i1_plus_1] - a->values[i1_minus_1];
  gi1 = a->values[i1_plus_2] - a->values[i1];
  *f1 = 0.5*(gi + r*(gi1-gi))*a->invDx;
  if(i1<=0) *f1 = 0.0;

  /**
   * reset indices to point to phi **/
  i1++; 
  i1_minus_1++;
  i1_plus_1++;
  i1_plus_2++;

#if 0  
  *value2 = a->values[i1] + r*(a->values[i1_plus_1] - a->values[i1]);
#endif
  *value2 = a->values[i1] + 0.5*r*(
			       (a->values[i1+2]+ a->values[i1-2] -2.0*a->values[i1])*r +
			       a->values[i1+2] - a->values[i1-2]
			       );
  gi  = a->values[i1_plus_1] - a->values[i1_minus_1];
  gi1 = a->values[i1_plus_2] - a->values[i1];
  *f2 = 0.5*(gi + r*(gi1-gi))*a->invDx;
  if(i1<=1) *f2 = 0.0;

  
  return;
}

void eamInterpolateDeriv(struct potentialarray_t *a, real_t r, int iType, int jType, real_t *value1, real_t *f1) {
  /**
   *
   * This routine will not crash if r is out of range.
   *
   * if ( r < a->x0) r = a->x0;
   * if ( r > a->xn)   r = a->xn;
   **/
   
  int i1, i1_minus_1, i1_plus_1,i1_plus_2;
  real_t gi, gi1;
  
  icalls++;
  r = (r<a->x0?a->x0:r); /* r = min(r,a->x0)  */
  r = (r>a->xn?a->xn:r);     /* r = max(r,a->xn)    */
  
  r = (r-a->x0)*(a->invDx) ;
  i1         = (int)floor(r);
  i1_plus_1  = i1+1;
  i1_plus_2  = i1+2;
  i1_minus_1 = i1-1;

  /* reset r to fractional distance */
  r = r - floor(r);

#if 0  
  *value1 = a->values[i1] + r*(a->values[i1_plus_1] - a->values[i1]);
#endif
  *value1 = a->values[i1] + 0.5*r*(
			       (a->values[i1+1]+ a->values[i1-1] -2.0*a->values[i1])*r +
			       a->values[i1+1] - a->values[i1-1]
			       );
  gi  = a->values[i1_plus_1] - a->values[i1_minus_1];
  gi1 = a->values[i1_plus_2] - a->values[i1];
  *f1 = 0.5*(gi + r*(gi1-gi))*a->invDx;
  if(i1<=0) *f1 = 0.0;

  return;
}


void eamInterpolateDeriv2Linear(struct potentialarray_t *a, real_t r,
			  int iType, int jType,
			  real_t *value1, real_t *value2,
			  real_t *f1, real_t *f2) {
  /**
   * simple linear interpolation routine
   *
   * This routine will not crash if r is out of range.
   *
   * if ( r < a->x0) r = a->x0;
   * if ( r > a->xn)   r = a->xn;
   **/
   
  int i0, i1, i2, i3, i4, i5, mi1, mi2;
  real_t g1, g2, g1a, g2a, g1b, g2b;
  double r1a;
  double r1b;
#if 0
    eamInterpolateDeriv2(a,r,iType,jType,value1,value2,f1,f2);
    return;
#endif
  icalls+=2;
  r = (r<a->x0?a->x0:r); /* r = min(r,a->x0)  */
  r = (r>a->xn?a->xn:r);     /* r = max(r,a->xn)    */
  
  r = (r-a->x0)*(a->invDx) ;

#ifdef INCLUDE_R_IN_DERIVATIVE		  
  r1a = (floor(r))/a->invDx + a->x0;
  r1b = (r1a + 1.0/a->invDx);
#endif

  i0 = 2 * ((int)floor(r));

  mi2 = i0-2;
  mi1 = i0-1;
  i1 = i0+1;
  i2 = i0+2;
  i3 = i0+3;
  i4 = i0+4;
  i5 = i0+5;

  /* reset r to fractional distance */
  r = r - floor(r);


  g1 = (a->values[i2] - a->values[i0]);
  *value1 = a->values[i0] + r*g1;

  g2 = (a->values[i3] - a->values[i1]);
  *value2 = a->values[i1] + r*g2;

#ifdef INCLUDE_R_IN_DERIVATIVE		  

  g1b = (a->values[i4] - a->values[i0])*a->invDx/2.0/r1b;
  if(i0<=0) g1a = g1b;
  else g1a = (a->values[i2] - a->values[mi2])*a->invDx/2.0/r1a;
  *f1 = (g1a + r * (g1b-g1a));

  g2b = (a->values[i5] - a->values[i1])*a->invDx/2.0/r1b;
  if(i1<=1) g2a = g2b;
  else g2a = (a->values[i3] - a->values[mi1])*a->invDx/2.0/r1a;
  *f2 = (g2a + r * (g2b-g2a));


#else

  g1b = (a->values[i4] - a->values[i0])*a->invDx/2.0;
  if(i0<=0) g1a = g1b;
  else g1a = (a->values[i2] - a->values[mi2])*a->invDx/2.0;
  *f1 = (g1a + r * (g1b-g1a));

  g2b = (a->values[i5] - a->values[i1])*a->invDx/2.0;
  if(i1<=1) g2a = g2b;
  else g2a = (a->values[i3] - a->values[mi1])*a->invDx/2.0;
  *f2 = (g2a + r * (g2b-g2a));

#endif
  return;
}

void eamInterpolateDerivLinear(struct potentialarray_t *a, real_t r, int iType, int jType, real_t *value1, real_t *f1) {
  /**
   * smiple linear interpolation routine
   *
   * This routine will not crash if r is out of range.
   *
   * if ( r < a->x0) r = a->x0;
   * if ( r > a->xn)   r = a->xn;
   **/
   
  int i0, i1, mi1, i2;
  real_t g1,g1a, g1b;

#if 0
  eamInterpolateDeriv(a,r,iType,jType,value1,f1);
  return;
#endif
  icalls++;
  r = (r<a->x0?a->x0:r); /* r = min(r,a->x0)  */
  r = (r>a->xn?a->xn:r);     /* r = max(r,a->xn)    */

  r = (r-a->x0)*(a->invDx) ;

  i0 = (int)floor(r);
  mi1 = i0-1;
  i1 = i0+1;
  i2 = i0+2;

  /* reset r to fractional distance */
  r = r - floor(r);

  g1 = (a->values[i1] - a->values[i0]);
  *value1 = a->values[i0] + r*g1;

  g1b = (a->values[i2] - a->values[i0])*a->invDx/2.0;
  if(i0<=0) g1a = g1b;
  else g1a = (a->values[i1] - a->values[mi1])*a->invDx/2.0;
    
  *f1 = (g1a + r * (g1b-g1a));

  return;
}


#if ( EAMLINEAR > 0 )
#define mymacro_eamInterpolateDeriv  eamInterpolateDerivLinear
#define mymacro_eamInterpolateDeriv2 eamInterpolateDeriv2Linear
#else
#define mymacro_eamInterpolateDeriv  eamInterpolateDeriv
#define mymacro_eamInterpolateDeriv2 eamInterpolateDeriv2
#endif

int eamForce(void *inS) {
  /**
   * calculates forces for the EAM potential **/
  simulation_t *s = (simulation_t *) inS;
  eampotential_t *pot = NULL;
  int ii, ij, i, j, jTmp;
  real_t r2,r;
  int nIBox;
  int nJBox;
  int *nbrBoxes;
  domain_t *iDomain;
  domain_t *jDomain;
  int iBox, jBox;
  real_t dx, dy, dz;
  real_t rcut2;
  real_t *rhosum;
  real_t *fprime;
  /*
  real_t *dRhoI2Sum;
  */
  real_t rhoTmp;
  real_t phiTmp;
  real_t dPhi, dRho;
  real_t bndX, bndY, bndZ;
  real_t *eatom;
  int ncc = 0;
  int *ncount;
  double etot;

  eatom = suAlignedCalloc(s->nAtoms*sizeof(double));
  if ( ! eatom) simulationAbort(-2,"Unable to allocate eatom.");
    
  ncount = suAlignedCalloc(s->nAtoms*sizeof(int));
  if ( ! ncount) simulationAbort(-2,"Unable to allocate ncount.");
    
    

  icount = 0;
  pot = getEamPot();
  if ( ! pot ) {
    simulationAbort(-2,"Unable to read in potential file for Ag (did you run 'tar zxvf pots.tgz' ?");
  }
  
  rcut2 = pot->cutoff*pot->cutoff;

  /**
   * allocate space **/
  rhosum = suAlignedCalloc(s->nAtoms*sizeof(real_t));
  fprime = suAlignedCalloc(s->nAtoms*sizeof(real_t));

  /**
   * point to energy and force and zero them out **/
  etot = 0.0;

  for(i=0; i<s->boxes.nboxes; i++) {
    zeroForcesInBox(s,i);
  }

  /**
   * The more I read this, the more it seems that I should transfer
   * atom position / force data to the box structure so that the
   * atoms in a box will be contiguous in memory.
   *
   * Right now the atoms are part of the overall simulation structure.
   * This means that potentially two atoms within a box may be widely
   * separated in main memory.  
   **/


  bndX = (s->periodicBounds[3]-s->periodicBounds[0])/2.;
  bndY = (s->periodicBounds[4]-s->periodicBounds[1])/2.;
  bndZ = (s->periodicBounds[5]-s->periodicBounds[2])/2.;

  memset(eatom,0,sizeof(real_t)*s->nAtoms);
  memset(ncount,0,sizeof(int)*s->nAtoms);
  
  for(iBox=0; iBox<s->boxes.nboxes; iBox++) /* loop over all boxes in system */
    {
      int inn;
      inn=0;

      iDomain = s->boxes.domains[iBox];
      nIBox = *(iDomain->nAtoms);
      nbrBoxes = getNeighborBoxes(s,iBox);
      
      for(jTmp=0; jTmp<NUMNEIGHBORS; jTmp++) /* loop over neighbor boxes */
	{
	  jBox = nbrBoxes[jTmp];
	  if((jBox<0) || (jBox < iBox ) ) continue;
	  jDomain = s->boxes.domains[jBox];
	  nJBox = *(jDomain->nAtoms);

	  for(ii=0; ii<nIBox; ii++) /* loop over atoms in iBox */
	    {
	      i = iDomain->id[ii];  /* the ij-th atom in iBox */

	      for(ij=0; ij<nJBox; ij++)  /* loop over atoms in iBox */
		{
		  j = jDomain->id[ij];  /* the ij-th atom in iBox */

		  if ( (iBox==jBox) &&(ij <= ii) ) continue;

		  /** dist2(r2,atomsInIBox,jDomain,ii,ij,dx,dy,dz);**/

		  dx = (iDomain->x[ii]-jDomain->x[ij]);
		  dy = (iDomain->y[ii]-jDomain->y[ij]);
		  dz = (iDomain->z[ii]-jDomain->z[ij]);

		  if(PERIODIC) {
		    /**
		     * find the periodic image **/
		    dx = (dx>bndX?dx-2.*bndX:dx);
		    dy = (dy>bndY?dy-2.*bndY:dy);
		    dz = (dz>bndZ?dz-2.*bndZ:dz);

		    dx = (-dx>bndX?dx+2.*bndX:dx);
		    dy = (-dy>bndY?dy+2.*bndY:dy);
		    dz = (-dz>bndZ?dz+2.*bndZ:dz);

		  }
		  r2 = dx*dx + dy*dy + dz*dz;

		  if(r2>rcut2) continue;
		  icount++;
		  ncount[i]++;
		  ncount[j]++;

		  r = sqrt(r2);

		  mymacro_eamInterpolateDeriv2(pot->phirho,r,0,0,&phiTmp,&rhoTmp,&dPhi,&dRho);
#ifdef INCLUDE_R_IN_DERIVATIVE		  
		  iDomain->fx[ii] += dPhi*dx;
		  iDomain->fy[ii] += dPhi*dy;
		  iDomain->fz[ii] += dPhi*dz;

		  jDomain->fx[ij] -= dPhi*dx;
		  jDomain->fy[ij] -= dPhi*dy;
		  jDomain->fz[ij] -= dPhi*dz;
#else
		  iDomain->fx[ii] += dPhi*dx/r;
		  iDomain->fy[ii] += dPhi*dy/r;
		  iDomain->fz[ii] += dPhi*dz/r;

		  jDomain->fx[ij] -= dPhi*dx/r;
		  jDomain->fy[ij] -= dPhi*dy/r;
		  jDomain->fz[ij] -= dPhi*dz/r;
#endif
		  etot += (double) phiTmp;
		  eatom[i] += phiTmp/2.0;
		  eatom[j] += phiTmp/2.0;
		  rhosum[i] += rhoTmp;
		  rhosum[j] += rhoTmp;
		} /* loop over atoms in jBox */

	    } /* loop over atoms in iBox */

	} /* loop over neighbor boxes */

    } /* loop over all boxes in system */

  if (0) {
    int iAtom;
    FILE *fpnc =NULL;
    fpnc = fopen("counterPPU.txt","w");
    for(iBox=0; iBox<s->boxes.nboxes; iBox++) {
      iDomain = s->boxes.domains[iBox];
      nIBox =  *(iDomain->nAtoms);
      for(iAtom=0; iAtom<nIBox; iAtom++) {
	i = iDomain->id[iAtom];  /* the ij-th atom in iBox */
	fprintf(fpnc,"%10d %10d\n",i+1,ncount[i]);
      }
    }
    fclose(fpnc);
  }

  /**
   * Now loop again and add in the F(rho(ij))
   **/

  for(iBox=0; iBox<s->boxes.nboxes; iBox++) /* loop over all boxes in system */
    {
      int inn;
      real_t fi, fiprime;
      real_t rhoij, rhoijprime;
      inn=0;

      iDomain = s->boxes.domains[iBox];
      nIBox =  *(iDomain->nAtoms);

      nbrBoxes = getNeighborBoxes(s,iBox);
      
      for(ii=0; ii<nIBox; ii++) /* loop over atoms in iBox */
	{
	  i = iDomain->id[ii];  /* the ij-th atom in iBox */
	  mymacro_eamInterpolateDeriv(pot->f,rhosum[i],0,0,&fi,&fiprime);
	  fprime[i] = fiprime;
	  etot += (double) fi;
	  eatom[i] += fi;
	}
    }
  for(iBox=0; iBox<s->boxes.nboxes; iBox++) /* loop over all boxes in system */
    {
      int inn;
      real_t fi, fiprime;
      real_t rhoij, rhoijprime;
      inn=0;

      iDomain = s->boxes.domains[iBox];
      nIBox =  *(iDomain->nAtoms);

      nbrBoxes = getNeighborBoxes(s,iBox);
      
      for(ii=0; ii<nIBox; ii++) /* loop over atoms in iBox */
	{
	  i = iDomain->id[ii];  /* the ij-th atom in iBox */
			
	  for(jTmp=0; jTmp<NUMNEIGHBORS; jTmp++) /* loop over neighbor boxes */
	    {
	      jBox = nbrBoxes[jTmp];
	      if((jBox < iBox) || (jBox < 0) ) continue;
	      
	      jDomain = s->boxes.domains[jBox];
	      nJBox = *(jDomain->nAtoms);

	      for(ij=0; ij<nJBox; ij++)  /* loop over atoms in iBox */
		{
		  j = jDomain->id[ij];  /* the ij-th atom in iBox */

		  if ((iBox==jBox) && (ij <= ii))  continue;
		  
		  dx = (iDomain->x[ii]-jDomain->x[ij]);
		  dy = (iDomain->y[ii]-jDomain->y[ij]);
		  dz = (iDomain->z[ii]-jDomain->z[ij]);
		  
		  if(PERIODIC) {
		    /**
		     * find the periodic image **/
		    dx = (dx>bndX?dx-2.0*bndX:dx);
		    dy = (dy>bndY?dy-2.0*bndY:dy);
		    dz = (dz>bndZ?dz-2.0*bndZ:dz);

		    dx = (-dx>bndX?dx+2.0*bndX:dx);
		    dy = (-dy>bndY?dy+2.0*bndY:dy);
		    dz = (-dz>bndZ?dz+2.0*bndZ:dz);
		  }
		  r2 = dx*dx + dy*dy + dz*dz;

		  if(r2>=rcut2) continue;
		  icount++;
		  r = sqrt(r2);

		  mymacro_eamInterpolateDeriv2(pot->phirho,r,0,0,&phiTmp,&rhoTmp,&dPhi,&dRho);
		  rhoij = rhoTmp;
		  rhoijprime = dRho;

#ifdef INCLUDE_R_IN_DERIVATIVE		  
		  iDomain->fx[ii] += (fprime[i]+fprime[j])*rhoijprime*dx;
		  iDomain->fy[ii] += (fprime[i]+fprime[j])*rhoijprime*dy;
		  iDomain->fz[ii] += (fprime[i]+fprime[j])*rhoijprime*dz;

		  jDomain->fx[ij] -= (fprime[i]+fprime[j])*rhoijprime*dx;
		  jDomain->fy[ij] -= (fprime[i]+fprime[j])*rhoijprime*dy;
		  jDomain->fz[ij] -= (fprime[i]+fprime[j])*rhoijprime*dz;
#else
		  iDomain->fx[ii] += (fprime[i]+fprime[j])*rhoijprime*dx/r;
		  iDomain->fy[ii] += (fprime[i]+fprime[j])*rhoijprime*dy/r;
		  iDomain->fz[ii] += (fprime[i]+fprime[j])*rhoijprime*dz/r;

		  jDomain->fx[ij] -= (fprime[i]+fprime[j])*rhoijprime*dx/r;
		  jDomain->fy[ij] -= (fprime[i]+fprime[j])*rhoijprime*dy/r;
		  jDomain->fz[ij] -= (fprime[i]+fprime[j])*rhoijprime*dz/r;
#endif

		} /* loop over atoms in jBox */
	    } /* loop over atoms in iBox */
	} /* loop over neighbor boxes */
    } /* loop over all boxes in system */


  suAlignedFree(rhosum);
  suAlignedFree(fprime);
  suAlignedFree(eatom);
  suAlignedFree(ncount);
  rhosum = fprime = NULL;

  s->e = (real_t) etot;


  return 0;
}


/**
 * utility comparison routined **/
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
  adiffpot("phi", a->phi, b->phi);
  adiffpot("rho", a->rho, b->rho);
  adiffpot("f", a->f, b->f);
  adiffpot("phirho", a->phirho, b->phirho);
  return;
}
