
#ifndef __DOCOMPUTE_H
#define __DOCOMPUTE_H

#include "pmdTypes.h"
#include "mytype.h"

void printArray(real_t* array, int n, char *name);
struct simflat_t *initSimFromCmdLine(int argc, char **argv);
void *do_compute_work(SimThreadData *data);

#endif
