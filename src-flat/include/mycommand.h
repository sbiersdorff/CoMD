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

/**
 * \brief a set of utility routines for parsing input
 *
 * We read input and fillout a command_t structure
 **/
#ifndef MYCOMMAND_H
#define MYCOMMAND_H

#include "mytype.h"

typedef struct command_t {
  char filename[1024];/**< the input file name **/
  char potdir[1024];  /**< the directory where eam potentials reside **/
  char potname[1024]; /**< the name of the potential **/
  int debug;          /**< a flag to determine whether to debugt or not **/
  int doeam;          /**< a flag to determine whether we're running EAM potentials **/
  int periodic;       /**< a flag that controls whether periodic boundary conditions are used **/
  int usegpu;         /**< a flag that controls whether OpenCL tries to target a gpu **/
  int nx;             /**< number of unit cells in x **/
  int ny;             /**< number of unit cells in y **/
  int nz;             /**< number of unit cells in z **/
  double bf;          /**< ratio of box size to cutoff radius **/
  double lat;         /**< lattice constant **/
  double defgrad;     /**< deformation gradient **/
  double temp;        /**< simulation temperature **/
} command_t;

extern void usage(char *exename);
extern void parseCommandLine(command_t *cmd, int argc, char **argv);
extern void printCmd(command_t *cmd);

#endif
