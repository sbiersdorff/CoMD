#ifndef __PMD_H_
#define __PMD_H_

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "mytype.h"
#include "mycommand.h"
#include "pmdTypes.h"
#include "eamTypes.h"
#include "ljTypes.h"
#include "domains.h"
#include "constants.h"

#define DEBUGLEVEL 0
#define PMDDEBUGPRINTF(xxx,...) {if(xxx>DEBUGLEVEL) printf(__VA_ARGS__);}
#define fPMDDEBUGPRINTF(xxx,...) {if(xxx>DEBUGLEVEL) fprintf(__VA_ARGS__);}

extern simflat_t *blankSimulation(struct pmd_base_potential_t *pot);
extern void destroySimulation(simflat_t **s);

extern simflat_t *fromFileASCII(char *filename, struct pmd_base_potential_t *pot);
extern simflat_t *fromFileGzip(char *filename, struct pmd_base_potential_t *pot);

extern void writeClsman(simflat_t *s, char *fname);
extern void printIt(simflat_t *s,FILE *fp);

extern double nTimeSteps(int n, simflat_t *s, real_t dt);

/* utility routines */
extern void breakinme();
extern void simulationAbort(int ineCode, char *inmsg);
extern double timeNow();
extern int computeForce(simflat_t *s);
extern void *do_compute_work(simflat_t *sim);
extern simflat_t *initSimFromCmdLine(int argc, char **argv, int *eam_flag, int *gpu_flag);

#endif
