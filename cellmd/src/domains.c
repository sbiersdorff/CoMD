/**
 * this file will split space up into domains so that
 * we will have fast neighbor listing **/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include "pmd.h"
#include "memUtils.h"

// if you define ONEBOXONLY, this should force an n^2 computation
// note that this will fail for _large_ number of atoms

#ifdef ONEBOXONLY
#define NOBOXES -1000
#else
#define NOBOXES 0
#endif
int nstack = 0;
void **stacker[100000] = {0};
static int jmax=0;
static int jmax2=0;


int getJmax() {return jmax;}

int checkStack(void *ptr) {
  int i;
  return 0;
  if ( ! ptr ) return 0;
  for(i=0; i<nstack;i++) if(ptr == stacker[i]) return 1;
  stacker[nstack] = ptr;
  nstack++;
  return 0;
}


#define check_ATOM_FREE(atom,element)					\
  {								\
    if(atom->element) suAlignedFree(atom->element);			\
    atom->element = NULL;					\
    if(checkStack(atom->element)) {				\
      simulationAbort(129,"double free of atom->" #element);	\
    }								\
    else printf("Successfully freed atom->" #element "\n");	\
    fflush(stdout);						\
  }

#define COPYARRAY(dest, src, element, size) { \
    if(src->element) memcpy(dest->element,src->element,size);\
    else memset(dest->element,-1,size); \
  }

void updateDomainPointers(domain_t *d) {

  int isz;
  int ioff;
  isz = *((int *) d->data);

  /**
   * repoints pointers within the domain structure **/
  d->iSize = (int *) d->data;
  d->nAtoms = d->iSize + 1;
  d->nPad = d->iSize + 2;
  d->nextID = d->iSize + 3;

  d->x = (real_t *)( d->data + DOMAINHEADERSIZE);
  d->y = d->x + isz;
  d->z = d->y + isz;

  d->fi = d->z + isz;

  d->fx = d->fi + isz;
  d->fy = d->fx + isz;
  d->fz = d->fy + isz;

  d->ea = d->fz + isz;

  d->px = d->ea + isz;
  d->py = d->px + isz;
  d->pz = d->py + isz;
  d->Move  = (int *)(d->pz + isz);
  d->iType = d->Move + isz;
  d->id    = d->iType + isz;
  return;
}

void *updateDomainPointersWithSizeAndData(domain_t *d, int isz, void *data) {

  int ioff;

  //  printf("update size is %d\n", isz);

  /**
   * repoints pointers within the domain structure **/
  d->data = data;
  d->iSize = (int *) d->data;
  d->nAtoms = d->iSize + 1;
  d->nPad = d->iSize + 2;
  d->nextID = d->iSize + 3;

  d->x = (real_t *)( d->data + DOMAINHEADERSIZE);
  d->y = d->x + isz;
  d->z = d->y + isz;

  d->fi = d->z + isz;

  d->fx = d->fi + isz;
  d->fy = d->fx + isz;
  d->fz = d->fy + isz;

  d->ea = d->fz + isz;

  d->px = d->ea + isz;
  d->py = d->px + isz;
  d->pz = d->py + isz;

  d->Move  = d->id + isz;
  d->iType = (int *)(d->pz+isz);
  d->id    = d->iType + isz;

  data = (advancePtr(void,data,DOMAINSIZE(isz)));
  return data;
}





static void *allocDomainData(int isize) {
  void *ptr;
  int is;
  int *iptr;

  is = DOMAINSIZE(isize);
			   
  ptr = suAlignedMalloc(is);
  if ( ! ptr ) simulationAbort(111,"Unable to allocate space for domain data");
  memset(ptr,0,is);
  iptr = (int *) ptr;
  iptr[0] = isize;
  return ptr;
}

domain_t *newDomain(int n) {
  domain_t *d;
  if ( n<0) n = 16;
  else n = n+n%2;
  d = suAlignedMalloc(sizeof(domain_t));
  if( ! d ) return NULL;
  d->data = allocDomainData(n);
  updateDomainPointers(d);
  return d;
}

#define aCOPYARRAY(dest, src, element, size) {			\
    if(src->element) memcpy(dest->element,src->element,size);\
    else {printf("no element src->%s\n",#element);memset(dest->element,-1,size);} \
  }


static void copyDomainArrays(domain_t *dest, domain_t *src, int isz) {
  *(dest->nAtoms) = *(src->nAtoms);

  COPYARRAY(dest, src, x, (isz*sizeof(real_t)));
  COPYARRAY(dest, src, y, (isz*sizeof(real_t)));
  COPYARRAY(dest, src, z, (isz*sizeof(real_t)));
  
  COPYARRAY(dest, src, fi, (isz*sizeof(real_t)));

  COPYARRAY(dest, src, fx, (isz*sizeof(real_t)));
  COPYARRAY(dest, src, fy, (isz*sizeof(real_t)));
  COPYARRAY(dest, src, fz, (isz*sizeof(real_t)));

  COPYARRAY(dest, src, ea, (isz*sizeof(real_t)));
  COPYARRAY(dest, src, px, (isz*sizeof(real_t)));
  COPYARRAY(dest, src, py, (isz*sizeof(real_t)));
  COPYARRAY(dest, src, pz, (isz*sizeof(real_t)));
  COPYARRAY(dest, src, Move,  (isz*sizeof(int)));

  COPYARRAY(dest, src, iType, (isz*sizeof(int)));
  COPYARRAY(dest, src, id,    (isz*sizeof(int)));

  return;
}

void expandDomain(domain_t *d) {

  domain_t oldDomain;
  int sizeold, sizenew;
  void *oldData;
  void *newData;

  if ( ! d ) return;

  // copy old pointers
  memcpy(&oldDomain, d, sizeof(domain_t));

  // allocate new pointers
  sizeold = *(d->iSize);
  sizenew = sizeold + 4;
  newData = allocDomainData(sizenew);

  memcpy(newData, oldDomain.data, DOMAINHEADERSIZE);
  updateDomainPointers(&oldDomain);
  
  // update new pointers
  d->data = newData;
  *((int *)(d->data)) = sizenew;

  updateDomainPointers(d);
#if 0
  printf("Domain expanded from %d to %d\n", *(d->nAtoms), *(oldDomain.nAtoms));
  printf("  sizeDomain expanded to %d from %d\n", *(d->iSize), *(oldDomain.iSize));
#endif
  copyDomainArrays(d,&oldDomain,sizeold);

  suAlignedFree(oldDomain.data);
  return;
  
}  

void destroyDomain(domain_t *d) {

  if ( ! d  ) return;
  if ( d->data) suAlignedFree(d->data);
  suAlignedFree(d);
  return;
}


/**
 * static box data **/
static real_t pdx=-1;
static real_t pdy=-1;
static real_t pdz=-1;

/**
 * static box functions **/
static int getIBoxFromxyz(simulation_t *s, int ix, int iy, int iz) {
  int ibox;
  boxes_t *b;

  b = &(s->boxes);
  ix = (ix>=0?ix:ix+s->boxes.nbx);
  ix = (ix>=s->boxes.nbx?ix-s->boxes.nbx:ix);

  iy = (iy>=0?iy:iy+s->boxes.nby);
  iy = (iy>=s->boxes.nby?iy-s->boxes.nby:iy);

  iz = (iz>=0?iz:iz+s->boxes.nbz);
  iz = (iz>=s->boxes.nbz?iz-s->boxes.nbz:iz);

  ibox = ix + iy*s->boxes.nbx + iz*s->boxes.nbx*s->boxes.nby;

  return ibox;

}

static void getIxyz(simulation_t *s, int ibox, int *ix, int *iy, int *iz) {
  *ix = ibox%(s->boxes.nbx);
  *iy = (ibox/s->boxes.nbx)%(s->boxes.nby);
  *iz = (ibox/s->boxes.nbx/s->boxes.nby)%(s->boxes.nbz);
  return;
}

static int inegs= 0;

int getBoxID2(simulation_t *s, real_t x, real_t y, real_t z) {
  int ibox,ix, iy, iz;
  ix = (int)((x-s->periodicBounds[0])/pdx);
  iy = (int)((y-s->periodicBounds[1])/pdy);
  iz = (int)((z-s->periodicBounds[2])/pdz);

  ix = (ix+s->boxes.nbx)%(s->boxes.nbx);
  iy = (iy+s->boxes.nby)%(s->boxes.nby);
  iz = (iz+s->boxes.nbz)%(s->boxes.nbz);

  inegs += (((ix<0) ||(iz<0) ||(iy<0) )>0?1:0);
  if(ix<0) ix += s->boxes.nbx;
  if(iy<0) iy += s->boxes.nby;
  if(iz<0) iz += s->boxes.nbz;
  
  ibox = ix + iy*s->boxes.nbx + iz*s->boxes.nbx*s->boxes.nby;

  return ibox;
}
int getBoxID(simulation_t *s, real_t x, real_t y, real_t z) {
  int ibox,ix, iy, iz;
  ix = (int)((x-s->periodicBounds[0])/pdx);
  iy = (int)((y-s->periodicBounds[1])/pdy);
  iz = (int)((z-s->periodicBounds[2])/pdz);


  inegs += (((ix<0) ||(iz<0) ||(iy<0) )>0?1:0);
  if(ix<0) ix += s->boxes.nbx;
  if(iy<0) iy += s->boxes.nby;
  if(iz<0) iz += s->boxes.nbz;

  ix = (ix+3*s->boxes.nbx)%(s->boxes.nbx);
  iy = (iy+3*s->boxes.nby)%(s->boxes.nby);
  iz = (iz+3*s->boxes.nbz)%(s->boxes.nbz);
  
  ibox = ix + iy*s->boxes.nbx + iz*s->boxes.nbx*s->boxes.nby;

  return ibox;
}

void copyAtom(int i, domain_t *idm, int j, domain_t *jdm) {
  /* copy atom iId in domain idm to atom jId in domain jdm */
  jdm->id[j] = idm->id[i];
  jdm->Move[j] = idm->Move[i];
  jdm->iType[j] = idm->iType[i];
  jdm->x[j] = idm->x[i];
  jdm->y[j] = idm->y[i];
  jdm->z[j] = idm->z[i];
  jdm->fx[j] = idm->fx[i];
  jdm->fy[j] = idm->fy[i];
  jdm->fz[j] = idm->fz[i];
  jdm->px[j] = idm->px[i];
  jdm->py[j] = idm->py[i];
  jdm->pz[j] = idm->pz[i];
  jdm->ea[j] = idm->ea[i];
  jdm->fi[j] = idm->fi[i];
}

void moveAtom(simulation_t *s, int iId, int iBox, int jBox) {
  int nj,ni;
  domain_t *idm, *jdm;
  idm = s->boxes.domains[iBox];
  jdm = s->boxes.domains[jBox];

  nj = jdm->nAtoms[0];
  copyAtom(iId, idm, nj, jdm);
  jdm->nAtoms[0]++;

  ni = idm->nAtoms[0] = idm->nAtoms[0] - 1;
  copyAtom(ni,idm,iId, idm);
  return;
}
  
void putAtomInBox(simulation_t *s,
		  int id, char flagMove, int iType,
		  real_t x,real_t y,real_t z,
		  real_t px,real_t py,real_t pz,
		  real_t fx,real_t fy,real_t fz) {

  /**
   * finds an appropriate box for an atom based on the
   * spatial cooridnates and puts it in there.
   *
   * reallocates memory if needed.
   **/

  domain_t *d;
  int ibox;
  int i;

  if(pdx<0.0 || pdy<0.0 || pdz<0.0) simulationAbort(-200,"Boxes not inited");
  /**
   * Find correct box.
   * for now, only one box **/
  ibox = getBoxID(s,x,y,z);
  //printf("Server %7d %7d\n",id,ibox);

  /**
   * point to correct box **/
  d = s->boxes.domains[ibox];

  /**
   * Reallocate if needed **/
  if(d->nAtoms[0] == d->iSize[0]) expandDomain(d);

  
  /**
   * push atom into primary period **/
  x = (x>=s->periodicBounds[3]?x-s->periodicBounds[3]:x);
  x = (x<0.0?x+s->periodicBounds[3]:x);
  
  y = (y>=s->periodicBounds[4]?y-s->periodicBounds[4]:y);
  y = (y<0.0?y+s->periodicBounds[4]:y);
  
  z = (z>=s->periodicBounds[5]?z-s->periodicBounds[5]:z);
  z = (z<0.0?z+s->periodicBounds[5]:z);
  /**
   * assign values to array elements **/

  i = d->nAtoms[0];
  d->nAtoms[0] = i+1;
  jmax = (jmax>i+1?jmax:i+1);
  //printf("id %d in box %d with isz=%d, natoms=%d\n", id, ibox, *(d->nAtoms), *(d->iSize));
  d->id[i] = id;
  d->Move[i] = flagMove;
  d->iType[i] = iType;
  d->x[i] = x;
  d->y[i] = y;
  d->z[i] = z;
  d->px[i] = px;
  d->py[i] = py;
  d->pz[i] = pz;
  d->fx[i] = fx;
  d->fy[i] = fy;
  d->fz[i] = fz;

  if ( id > s->boxes.idmax)    s->boxes.idmax = id;

  return;
}


static int *p_dcounts=NULL;
static int myndx, myndy,myndz,myndoms;

void setDomainSizes(int *domainCounts, int ndx, int ndy, int ndz, int nDomains, real_t *perb) {
//int i;
  p_dcounts = domainCounts;
  myndx = ndx;
  myndy = ndy;
  myndz = ndz;
  myndoms = nDomains;
  pdx = perb[3]-perb[0];
  pdy = perb[4]-perb[1];
  pdz = perb[5]-perb[2];
  pdx = pdx/(real_t) ndx;
  pdy = pdy/(real_t) ndy;
  pdz = pdz/(real_t) ndz;
  return;
}

void allocDomains(simulation_t *s) {
  /**
   * for now only one box with all atoms **/

  boxes_t *b;
  int i;

  real_t cutoff;


  if ( ! s ) return;
  cutoff = s->boxes.rcut;
  
  /**
   * grab a pointer to the boxes structure **/
  b = &(s->boxes);
  

  if ( p_dcounts) {
    b->nbx = myndx;
    b->nby = myndy;
    b->nbz = myndz;
    b->nboxes = myndoms;
    (b->domains) = suAlignedMalloc((b->nboxes*sizeof(domain_t *)));
    if ( ! b->domains) {
      destroyBoxes(s);
      simulationAbort(102,"Unable to allocate atoms in allocBoxes()");
    }
    for(i=0; i<b->nboxes; i++) b->domains[i] = newDomain(p_dcounts[i]);
  }
  else {
      

    /**
     * The number of boxes.
     * This is initially set to unity, but will be replaced
     * by boxes of size 1.5*cutoff when we are done.
     **/
    if(NOBOXES+cutoff <= 0.0 ) {
      b->nbx = b->nby = b->nbz = 1;
      pdx = s->periodicBounds[3];
      pdy = s->periodicBounds[4];
      pdz = s->periodicBounds[5];
    }
    else {
      int isz = 1024;
      float factor = 1.03;
      pdx = s->periodicBounds[3]-s->periodicBounds[0];
      pdy = s->periodicBounds[4]-s->periodicBounds[1];
      pdz = s->periodicBounds[5]-s->periodicBounds[2];
      b->nbx = (int)floor(pdx/factor/cutoff);
      b->nby = (int)floor(pdy/factor/cutoff);
      b->nbz = (int)floor(pdz/factor/cutoff);

      b->nbx = (b->nbx>0?b->nbx:1);
      b->nby = (b->nby>0?b->nby:1);
      b->nbz = (b->nbz>0?b->nbz:1);
      b->nboxes = b->nbx*b->nby*b->nbz;
      if(b->nbx<3 || b->nby<3 || b->nbz<3 ) b->nbx = b->nby = b->nbz = 1;
      b->nboxes = b->nbx*b->nby*b->nbz;
      pdx = pdx/(real_t)b->nbx;
      pdy = pdy/(real_t)b->nby;
      pdz = pdz/(real_t)b->nbz;
    
    }
    b->nboxes = b->nbx*b->nby*b->nbz;

    (b->domains) = suAlignedMalloc((b->nboxes*sizeof(domain_t *)));
    if ( ! b->domains) {
      destroyBoxes(s);
      simulationAbort(102,"Unable to allocate atoms in allocBoxes()");
    }
    for(i=0; i<b->nboxes; i++) b->domains[i] = newDomain(-1);
  }


  return;
}

void destroyDomains(simulation_t *s) {
  /**
   * release all memory associated with domains
   **/
  int i;
  domain_t **ppd;
  if ( ! s ) return;
  if ( ! s->boxes.domains) return;
  ppd = s->boxes.domains;
  for(i=0; i<s->boxes.nboxes; i++) destroyDomain(ppd[i]);
  suAlignedFree(ppd);
  s->boxes.domains = NULL;
  return;
}


void redoBoxes(simulation_t *s) {
  simulationAbort(-3,"Placeholder routine redoBoxes() called");
  return;
}

void reBoxAtom(simulation_t *s, int iatom) {
  simulationAbort(-3,"Placeholder routine reBoxAtom() called");
  return;
}

void destroyBoxes(simulation_t *s) {
  boxes_t *b;
  b = &(s->boxes);
  destroyDomains(s);
  memset(b,0,sizeof(boxes_t));
  return;
}

#define ADDIF 1

int *getNeighborBoxes(simulation_t *s, int iBox) {
  static int actualNbrs[1+NUMNEIGHBORS];
  int *nbrs;
  int ix, iy, iz;
  int in;
  boxes_t *b;
  b = &(s->boxes);

  memset(actualNbrs,-1,1+NUMNEIGHBORS*sizeof(int));
  nbrs = actualNbrs+1;

  if(b->nboxes == 1) {
    nbrs[-1] = 1;
    nbrs[0] = 0;
    return nbrs;
  }

  getIxyz(s, iBox, &ix, &iy, &iz);

  in = 0;
  
  ix--;iy--;iz--;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;

  iz++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;

  iz++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;

  
  iz-=2; iy++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;

  iz++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;

  iz++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;

  iz-=2; iy++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;

  iz++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;

  iz++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;

  iz-=2; iy-=2;ix++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;
  iz++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;
  iz++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);

  if(ADDIF) in++;
  iz-=2; iy++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;
  iz++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;
  iz++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);

  if(ADDIF) in++;
  iz-=2; iy++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;
  iz++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;
  iz++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);

  if(ADDIF) in++;
  iz-=2; iy-=2;ix++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;
  iz++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;
  iz++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);

  if(ADDIF) in++;
  iz-=2; iy++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;
  iz++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;
  iz++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);

  if(ADDIF) in++;
  iz-=2; iy++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;
  iz++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(ADDIF) in++;
  iz++;
  nbrs[in] = getIBoxFromxyz(s,ix,iy,iz);
  if(! ADDIF) nbrs[in] = -1;

  nbrs[-1] = in;
  return nbrs;
}

void zeroForcesInBox(simulation_t *s, int i) {
  int n;
  n = *(s->boxes.domains[i]->nAtoms);
  if ( n ) {
    memset(s->boxes.domains[i]->fx,0,n*sizeof(real_t));
    memset(s->boxes.domains[i]->fy,0,n*sizeof(real_t));
    memset(s->boxes.domains[i]->fz,0,n*sizeof(real_t));
  }

  return;
}






