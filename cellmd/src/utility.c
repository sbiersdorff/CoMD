/**
 * The data for our simulation will be held in a
 * structure of arrays.  To make it go fast on the
 * CBE we have decided to separate out the three
 * positions and three momenta into separate
 * arrays.  This will allow for more overlap
 * between memory access and computation.
 **/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include "pmdTypes.h"
#include "pmd.h"
#include "memUtils.h"
#define offsetPtr(ptr,offset) (((char *)ptr) + offset)

/* some utility routines */
double timeNow() {
  struct timeval tnow;
  double t;
  gettimeofday(&tnow,NULL);
  t = (double) tnow.tv_sec + ((double)tnow.tv_usec)*1.0e-6;
  return t;
}

void breakinme() {/** a blank break**/return;}
int suAbort(int status,char *message) {
  fprintf(stderr,"\n\n\n   ABORT!: \n   %s\n\n\n",message);
  breakinme();				\
  exit(status);
}
char *dupString(char *str) {
  char *s;
  s = suAlignedCalloc((strlen(str)+1)*sizeof(char));
  strcpy(s,str);
  return s;
}

void simulationAbort(int ineCode, char *inmsg) {	
  int eCode;					
  char *msg;					
  eCode = ineCode;				
  msg = inmsg;				
  if(! msg) {				
    msg = "Unknown error";		
  }					
  fprintf(stderr,"\n\n"			
	  "    Error code: %d\n"	
	  "    %s\n\n",			
	  eCode, msg);			
  exit(eCode);				
}

struct simulation_t *blankSimulation(struct pmd_base_potential_t *pot) {
  simulation_t *s;

  s = suAlignedMalloc((1*sizeof(simulation_t)));

  if(pot) {
    s->pot = pot;
    s->force = pot->force;
    s->boxes.rcut = pot->cutoff;
  }
  else {
    extern ljpotential_t *getLJPot();
    s->pot = (pmd_base_potential_t *) getLJPot();
    if(s->pot) s->force = s->pot->force;
    s->boxes.rcut = s->pot->cutoff;
  }
  
  if ( s->boxes.rcut < 0.0) simulationAbort(-1,"Blanksimulation() got a negative cutofff");
  
  
  return s;

}


void destroySimulation(simulation_t **ps) {
  /**
   * frees all data associated with *s and frees *s **/
  simulation_t *s;
  pmd_base_potential_t *pot;

  if ( ! ps ) return;

  s = *ps;

  if ( ! s ) return;
  pot = s->pot;
  if ( pot) pot->destroy(&pot);
  if ( ! s->contiguous) {
    if(s->comment) suAlignedFree(s->comment);
    destroyBoxes(s);
  }
  suAlignedFree(s);
  *ps = NULL;
  return;
}

