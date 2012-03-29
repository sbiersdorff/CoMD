#ifndef __EAM_TYPES_H_
#define __EAM_TYPES_H_

#include <pmd.h>

typedef struct potentialarray_t {
  int n;
  real_t x0;
  real_t xn;
  real_t invDx;
  real_t *values;
} potentialarray_t;

typedef struct eampotential_t {
  real_t cutoff;       /* the potential cutoff distance */
  real_t mass;           /* mass of atoms */
  int (*force)(void *s);
  void (*destroy)(void **pot);
  struct potentialarray_t *phi; /* the phi array */
  struct potentialarray_t *rho; /* the rho array */
  struct potentialarray_t *f;   /* the F array   */
  struct potentialarray_t *phirho; /* the phi array */
} eampotential_t;


extern void eamDestroy(void **pPot);
extern eampotential_t *eamReadASCII(char *dir, char *potname);
extern real_t eamInterpolate(potentialarray_t *a, real_t r, int iType);
extern void eamInterpolateDeriv2(struct potentialarray_t *a, real_t r,
				 int iType, int jType,
				 real_t *value1, real_t *value2,
				 real_t *f1, real_t *f2);
extern void eamInterpolateDeriv(struct potentialarray_t *a, real_t r,int iType, int jType, real_t *v1, real_t *f1);

extern struct eampotential_t *getEamPot();
extern struct pmd_base_potential_t *setEamPot(char *dir, char *file);
extern struct pmd_base_potential_t *setEamPotFromPotential(struct pmd_base_potential_t *inPot);
extern struct pmd_base_potential_t *getEamPotentialFromClient();

/**
 *  createEAMPotentialFromData(..) creates a cellMD type potential from
 *  standard eam tables you may have lying around.  
 **/
struct pmd_base_potential_t *createEAMPotentialFromData(int ftnIndexing,
           int n_phi,double x0_phi, double xN_phi, double invDx_phi, double *values_phi,
	   int n_rho,double x0_rho, double xN_rho, double invDx_rho, double *values_rho,
           int n_f,  double x0_f,   double xN_f,   double invDx_f,   double *values_f );
  

#define eamGetPhiRho(pot,r,iType,phiRet,rhoRet) {eamInterpolate2(pot->phirho,r,iType,phiRet,rhoRet);}
#define eamGetPhi(pot,r,iType) (eamInterpolate(pot->phi,r,iType))
#define eamGetRho(pot,r,iType) (eamInterpolate(pot->rho,r,iType))
#define eamGetF(pot,rho,iType) (eamInterpolate(pot->f,rho,iType))

#endif
