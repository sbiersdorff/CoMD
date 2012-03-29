#include "pmd.h"
#include <time.h>
#define DOREADTIMERS 4
#ifndef GZIPSUPPORT
simulation_t *fromFileASCII(char *filename, struct pmd_base_potential_t *pot) {
  /**
   * Reads in a simulation from a file written by clsman.
   *
   * Initially only mode 0 files are read
   **/
  simulation_t *s = NULL;
  FILE *fp;
  int i,itype;
  real_t x,y,z,px,py,pz;
  int natoms, nmove;
#ifdef DOREADTIMERS  
  clock_t start,old;
  clock_t count=0;
#endif

  PMDDEBUGPRINTF(0,"tmp file in read is: %s\n", filename);
  fp = fopen(filename,"r");
  if( ! fp ) return NULL;

#ifdef DOREADTIMERS  
    start = clock();
#endif
  s = blankSimulation(pot);
  if ( ! s ) {
    fclose(fp);
    return s;
  }


  /* read in natom and nmove */
  fscanf(fp,"%d%d\n",&natoms,&nmove);
  PMDDEBUGPRINTF(0,"\n   natoms = %d\n   nmove = %d\n\n",natoms,nmove);
  s->nAtoms = natoms;
  s->nMove = nmove;

  /* read in comment */
  s->comment = suAlignedCalloc(1024*sizeof(char));
  if(0) {
    char fmta[128];
    sprintf(fmta,"%%%d\n",s->nMove);
    fscanf(fp,fmta,s->comment);
  }
  else if ( nmove) {
    fgets(s->comment,1024,fp);
  }
  PMDDEBUGPRINTF(0,"\n   comment = %s\n",s->comment);

  /* read in periodic  boundaries */
  s->periodicBounds[0] = s->periodicBounds[1] = s->periodicBounds[2] = (real_t) 0.0;

  fscanf(fp,FMT1 " " FMT1 " " FMT1 "\n",s->periodicBounds+3,s->periodicBounds+4,s->periodicBounds+5);
  PMDDEBUGPRINTF(0,"Periodic bounds: " FMT1 " " FMT1 " " FMT1 "\n",
		 s->periodicBounds[3],s->periodicBounds[4],s->periodicBounds[5]);

#if DOREADTIMERS  >2
    old = clock();
#endif
  allocDomains(s);
#if DOREADTIMERS  >2
    count = clock()-old+count;
#endif

  for(i=0; i<s->nAtoms; i++,x++,y++,z++,px++,py++,pz++,itype++) {
    PMDDEBUGPRINTF(0,"Reading %d\r", i);
    fscanf(fp,FMT1 " " FMT1 " " FMT1 " %d " FMT1 " " FMT1 " " FMT1 "\n",
	   &x,&y,&z,&itype,&px,&py,&pz);
    PMDDEBUGPRINTF(0,"%4d: " FMT1 " " FMT1 " " FMT1 " %d " FMT1 " " FMT1 " " FMT1 "\n",
		   i,x,y,z,itype,px,py,pz);
#if DOREADTIMERS  >2
    old = clock();
#endif
    putAtomInBox(s,i,(char)(i<nmove),itype,
		 x,y,z,
		 px,py,pz,
		 (real_t) 0.0,(real_t) 0.0,(real_t) 0.0);
#if DOREADTIMERS  >2
    count = clock()-old+count;
#endif
  }
#ifdef DOREADTIMERS  
    old = clock();
#endif
  if(0) printf("\n    ---- Read file took %.2fs of which boxingtime = %.2fs\n\n",	 
	 (float)(old-start)/(float)(CLOCKS_PER_SEC),(float)(count)/(float)(CLOCKS_PER_SEC));
  fclose(fp);
    

  return s;
}
simulation_t *fromFileGzip(char *filename, struct pmd_base_potential_t *pot) {
  printf("    fromFileGzip(): Trying fromFileASCII()n");
  return (fromFileASCII(filename,pot));
}
#else
#include "zlib.h"
simulation_t *fromFileGzip(char *filename, struct pmd_base_potential_t *pot) {
  /**
   * Reads in a simulation from a file written by clsman.
   *
   * Initially only mode 0 files are read
   **/
#ifdef DOREADTIMERS  
  clock_t start,old;
  int count;
#endif
  simulation_t *s = NULL;
  gzFile *fp;
  int i,itype;
  real_t x,y,z,px,py,pz;
  int natoms, nmove;
  char zstr[16348];

#ifdef DOREADTIMERS  
  start = clock();
  count = 0;
#endif
  PMDDEBUGPRINTF(0,"tmp file in read is: %s\n", filename);
  fp = gzopen(filename,"rb");
  if( ! fp ) return NULL;

  s = blankSimulation(pot);
  if ( ! s ) {
    gzclose(fp);
    return s;
  }


  /* read in natom and nmove */
  gzgets(fp,zstr,sizeof(zstr));
  sscanf(zstr,"%d%d\n",&natoms,&nmove);
  PMDDEBUGPRINTF(0,"\n   -0natoms = %d\n   nmove = %d\n\n",natoms,nmove);
  s->nAtoms = natoms;
  s->nMove = nmove;

  /* read in comment */
  s->comment = suAlignedCalloc(1024*sizeof(char));
  if(0) {
    char fmta[128];
    sprintf(fmta,"%%%d\n",s->nMove);
    gzgets(fp,zstr,sizeof(zstr));
    sscanf(zstr,fmta,s->comment);
  }
  else {
      gzgets(fp,zstr,sizeof(zstr));
      if(1)strncpy(s->comment,zstr,1024);
      else s->comment = "no comment";
  }

  PMDDEBUGPRINTF(0,"\n   comment = %s\n",s->comment);

  /* read in periodic  boundaries */
  s->periodicBounds[0] = s->periodicBounds[1] = s->periodicBounds[2] = (real_t) 0.0;

  gzgets(fp,zstr,sizeof(zstr));
  sscanf(zstr,FMT1 " " FMT1 " " FMT1 "\n",s->periodicBounds+3,s->periodicBounds+4,s->periodicBounds+5);
  PMDDEBUGPRINTF(0,"Periodic bounds: " FMT1 " " FMT1 " " FMT1 "\n",
		 s->periodicBounds[3],s->periodicBounds[4],s->periodicBounds[5]);

  allocDomains(s);

  for(i=0; i<s->nAtoms; i++) {
    PMDDEBUGPRINTF(0,"Reading %d\r", i);
    gzgets(fp,zstr,sizeof(zstr));
    sscanf(zstr,FMT1 " " FMT1 " " FMT1 " %d " FMT1 " " FMT1 " " FMT1 "\n",
	   &x,&y,&z,&itype,&px,&py,&pz);
    PMDDEBUGPRINTF(0,"%4d: " FMT1 " " FMT1 " " FMT1 " %d " FMT1 " " FMT1 " " FMT1 "\n",
		   i,x,y,z,itype,px,py,pz);
#ifdef DOREADTIMERS  
    old = clock();
#endif
    putAtomInBox(s,i,(char)(i<nmove),itype,
		 x,y,z,
		 px,py,pz,
		 (real_t) 0.0,(real_t) 0.0,(real_t) 0.0);
#ifdef DOREADTIMERS  
    count = clock()-old+count;
#endif
  }
  gzclose(fp);
#ifdef DOREADTIMERS  
  old=clock();
  if(0) printf("\n    ---- Read file took %.2fs of which boxingtime = %.2fs\n\n",	 
		 (float)(old-start)/(float)(CLOCKS_PER_SEC),(float)(count)/(float)(CLOCKS_PER_SEC));
#endif
  return s;
}
simulation_t *fromFileASCII(char *filename, struct pmd_base_potential_t *pot) {
  return(fromFileGzip(filename,pot));
}
#endif



