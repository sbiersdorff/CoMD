
#ifndef __UTILITY_H
#define __UTILITY_H

#include "pmdTypes.h"

double timeNow();
void breakinme();
char *dupString(char *str);
void simulationAbort(int ineCode, char *inmsg);
struct simflat_t *blankSimulation(struct pmd_base_potential_t *pot);
void destroySimulation(simflat_t **ps);

#endif
