/**
 * The include file for boxit.c
 * contains function declarations
 **/
#ifndef __BOXIT_H_
#define __BOXIT_H_

#include "pmdTypes.h"

#ifndef advancePtr
#define advancePtr(type,p,is) ((type *)(((char *)(p))+(is)))
#endif

#define NUMNEIGHBORS 27
#define DMA_PAD(iSize) (((iSize)+15) & ~15)
#define DOMAINHEADERSIZE DMA_PAD(		\
				 4*sizeof(int)	\
				  )

#define DOMAINSIZE(isize) DOMAINHEADERSIZE +	\
  (						\
   11*DMA_PAD(isize*sizeof(real_t))		\
   + 3*DMA_PAD(isize*sizeof(int))		\
    )

extern void zeroForcesInBox(simulation_t *s, int iBox);

extern void allocDomains(simulation_t *s);
extern void destroyBoxes(simulation_t *s);
extern void redoBoxes(simulation_t *s);

extern int *getNeighborBoxes(simulation_t *s, int iBox);
extern void reBoxAtom(simulation_t *s, int iatom);
extern void putAtomInBox(simulation_t *s,
			 int id, char flagMove, int iType,
			 real_t x,real_t y,real_t z,
			 real_t px,real_t py,real_t pz,
			 real_t fx,real_t fy,real_t fz);
extern int      findAtomByID(simulation_t *s, int id, int *boxID, int *iInBox);
extern void     eliminateBlankBoxes(simulation_t *s);
extern int      getIBoxFromMerge(simulation_t *s, int iBoxM);
extern void updateSPUDomainInfo(simulation_t *s);
extern void updateSPUServerDomainInfo(simulation_t *s);
extern void updateDomainPointers(domain_t *d);
extern void *updateDomainPointersWithSizeAndData(domain_t *d, int isz, void *data);

/**
 * Some useful defines **/
#define  getAtomsInBox(s, iBox) (s->boxes.atoms[iBox])
#define  getNAtomsInBox(s, iBox) (s->boxes.atoms[iBox]->nAtoms)

#endif
