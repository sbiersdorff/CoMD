/*

Copyright (c) 2011, Los Alamos National Security, LLC All rights
reserved.  Copyright 2011. Los Alamos National Security, LLC. This
software was produced under U.S. Government contract DE-AC52-06NA25396
for Los Alamos National Laboratory (LANL), which is operated by Los
Alamos National Security, LLC for the U.S. Department of Energy. The
U.S. Government has rights to use, reproduce, and distribute this
software.

NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY
WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF
THIS SOFTWARE.

If software is modified to produce derivative works, such modified
software should be clearly marked, so as not to confuse it with the
version available from LANL.

Additionally, redistribution and use in source and binary forms, with
or without modification, are permitted provided that the following
conditions are met:

· Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

· Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

· Neither the name of Los Alamos National Security, LLC, Los Alamos
  National Laboratory, LANL, the U.S. Government, nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS
ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef __PMD_H_
#define __PMD_H_

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "mytype.h"
#include "constants.h"
#include "mycommand.h"
#include "pmdTypes.h"
#include "eamTypes.h"
#include "ljTypes.h"
#include "domains.h"

#define DEBUGLEVEL 0
#define PMDDEBUGPRINTF(xxx,...) {if(xxx>DEBUGLEVEL) printf(__VA_ARGS__);}
#define fPMDDEBUGPRINTF(xxx,...) {if(xxx>DEBUGLEVEL) fprintf(__VA_ARGS__);}

extern simflat_t *blankSimulation(struct pmd_base_potential_t *pot);
extern void destroySimulation(simflat_t **s);

extern simflat_t *fromFileASCII(char *filename, struct pmd_base_potential_t *pot);
extern simflat_t *fromFileGzip(char *filename, struct pmd_base_potential_t *pot);
extern simflat_t *fromFileTim(char *filename, struct pmd_base_potential_t *pot);

extern void writeClsman(simflat_t *s, char *fname);
extern void printIt(simflat_t *s,FILE *fp);

extern double nTimeSteps(int n, simflat_t *s, real_t dt);

/* utility routines */
extern void breakinme();
extern void simulationAbort(int ineCode, char *inmsg);
extern double timeNow();
extern int computeForce(simflat_t *s);
extern void *do_compute_work(SimThreadData *data);
extern simflat_t *initSimFromCmdLine(int argc, char **argv);

#endif
