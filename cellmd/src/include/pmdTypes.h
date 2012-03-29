#ifndef __PMDTYPES_H_
#define __PMDTYPES_H_

#include "mytype.h"
#include "memUtils.h"


// the base form off of which all potentials will be set
typedef struct pmd_base_potential_t {
  real_t cutoff;         /* potential cutoff distance */
  real_t mass;           /* mass of atoms */
  int (*force)(void *s); /**< the actual parameter is "struct simulation_t *s" **/
  void (*destroy)(void *pot); /* destruction of the potential */
} pmd_base_potential_t;
  

typedef struct ljpotential_t {
  real_t cutoff;       
  real_t mass;           /* mass of atoms */
  int (*force)(void *s);
  void (*destroy)(void **pot); 
  real_t sigma;
  real_t epsilon;
} ljpotential_t;


typedef struct domain_t {
  void *data;        /** the actual data
		      * isize, natoms, npad, nextid, 
		      * 27 neighbors
		      * 27 neighbor EA
		      * x1--xn, xpad
		      * ditto for y,z,dx,dy,dz,fx,fy,fz,px,py,pz
		      * itype1..itypen,itypepad
		      * id1...idn, idpad,
		      * moveflag1, moveflagn, movepad
		      */

  /**
   * all elements below point into data above **/

  int *iSize;   /* current size of atom arrays  **/
  int *nAtoms;  /* total number of atoms in box **/
  int *nPad;   /* Is total number of atoms an odd number **/
  int *nextID;  /* nextid for the spu           **/
  
  /** the three positions **/
  real_t *x;
  real_t *y;
  real_t *z;

  /** the value of rhobar **/
  real_t *fi;

  /** the three forces **/
  real_t *fx;
  real_t *fy;
  real_t *fz;

  /* energy per atom **/
  real_t *ea;

  /** the three momenta **/
  real_t *px;
  real_t *py;
  real_t *pz;

  int *Move;  /* Can the atom move?           **/

  int *iType;  /* the type of atoms, not initially used. **/  
  int *id;     /* The original ID of the atom  **/

} domain_t;


typedef struct boxes_t {
  /**
   * a struct to hold information on space partitioning */
  unsigned int size_b;
  int nbx;      /* number of boxes in x-direction       */
  int nby;      /* number of boxes in y-direction       */
  int nbz;      /* number of boxes in z-direction       */

  int nboxes;   /* total number of boxes                        */


  real_t rcut;  /* cutoff for box size, ignored                 */
  int    idmax; /* the highest id in this box                   */
  domain_t **domains;

} boxes_t;


typedef struct simulation_t {
  /**
   * a (short) description of the simulation **/

  /** the total energy of the simulation **/

  unsigned size_s;
  unsigned size_dh;
  unsigned nBytes;
  unsigned nBytesSwapped;
  unsigned contiguous;
  unsigned offset_Domain_data;
  unsigned offset_Domain_array;
  unsigned offset_Domain_pointers;

  int nAtoms;  /* The total number of atoms in the simulation  **/
  int nMove;   /* The number of moving atoms in the simulation **/

  real_t e;
  

  real_t periodicBounds[6]; /* The periodic boundaries*/

  /**
   * a superfluous structure, but I believe it will
   * make remote communication easier.
   *
   * I might make this an pointer in later versions.
   **/
  
  boxes_t boxes; 
  char *comment;

  /* the potential */
  pmd_base_potential_t *pot;
  
  /**
   * the force function **/
  int (*force)(void *s);
} simulation_t;

#endif
