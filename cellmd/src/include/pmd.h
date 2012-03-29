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
#include "domains.h"

#define DEBUGLEVEL 0
#define PMDDEBUGPRINTF(xxx,...) {if(xxx>DEBUGLEVEL) printf(__VA_ARGS__);}
#define fPMDDEBUGPRINTF(xxx,...) {if(xxx>DEBUGLEVEL) fprintf(__VA_ARGS__);}

extern int getJmax();
extern simulation_t *blankSimulation(struct pmd_base_potential_t *pot);
extern void destroySimulation(simulation_t **s);
extern simulation_t *fromFileASCII(char *filename, struct pmd_base_potential_t *pot);
extern simulation_t *fromFileGzip(char *filename, struct pmd_base_potential_t *pot);
extern real_t getLJCutoff();
extern int cmdPPU(command_t *pcmd);
extern void doPPU(simulation_t *sim);
extern void breakinme();
extern int suAbort(int status,char *message);
extern void simulationAbort(int ineCode, char *inmsg);
extern double nowTime();
#endif
