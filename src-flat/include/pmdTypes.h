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

#ifndef __PMDTYPES_H_
#define __PMDTYPES_H_

#include "mytype.h"
#include "memUtils.h"

#define MAXATOMS 64
#define SIM_NOSTATE 0
#define SIM_ALLOCED 1


/**
 * \brief the base form off of which all potentials will be set
 *
 * - All potentials are expected to conform to the
 * following units:
 *   -# atomic distances are in Angstroms
 *   -# atomic masses are in AMUs (atomic mass units)
 *   -# forces are returned in hartrees/angstrom
 *   -# energies are computed in hartrees
 * Note that phi, rho, and f are of type potentialarray_t.
 *
 **/
typedef struct pmd_base_potential_t {
  real_t cutoff;         /**< potential cutoff distance in Angstroms **/
  real_t mass;           /**< mass of atoms in atomic mass units **/
  int (*force)(void *s); /**< the actual parameter is "struct simulation_t *s" **/
  void (*destroy)(void *pot); /**< destruction of the potential **/
} pmd_base_potential_t;



#ifdef USE_IN_SITU_VIZ
typedef struct {   /**< Neighbor info used for centrosymmetry parameter **/
  real3    r;
  real_t   rsq;
  int      index;
} Neighbor;
typedef Neighbor Neighbor12[12];
#endif

typedef struct fileAtom {
  float  x, y, z;       
  float  bond_order;
  float  centrosymmetry;
} FileAtom;
  

/** \def array offset by n boxes **/
#define arrayAtBox(array,boxid) ((array)+MAXATOMS*boxid)

/**
 * The basic flat simulation data structure with MAXATOMS in every box **/
typedef struct simflat_t {
  int stateflag; /**< unused for now **/
  int nbx[3];    /**< number of boxes in each dimension **/
  int nboxes;    /**< total number of boxes **/
  int ntot; /**< total number of atoms**/
  real3 bounds; /**< periodic bounds**/
  real3 boxsize; /**< size of domains**/

  real3 *dcenter; /**< an array that contains the center of each box **/
  int *natoms;    /**< the total number of atoms in the simulation **/
  int *id;     /**< The original ID of the atom  **/
  int *iType;  /**< the type of atoms**/
  real_t *mass; /**< mass of the atoms**/
  real3 *r; /**< positions**/
  real3 *p; /**< momenta of atoms**/
  real4 *f; /**< fx, fy, fz, energy**/
  real_t *rho;   /**< rhosum for EAM potential**/
  real_t *fi;    /**< rhobar for EAM potential**/

#ifdef USE_IN_SITU_VIZ
  real_t *centro;
#endif

  /** the total potential energy of the simulation **/
  real_t e;    /**< the total energy of the system **/
  int nAtoms;  /**< The total number of atoms in the simulation  **/


  pmd_base_potential_t *pot; /**< the potential**/

  char *comment; /**< free form string that describes the simulation **/
} simflat_t;

/*
  add-onns for visualization
*/

#ifdef USE_IN_SITU_VIZ
#include "AtomVisualize.h"

#else
typedef void AtomVisualize;
#endif

typedef struct simThreadData
{
	AtomVisualize* viz;
	simflat_t* sim;
	int niter;
	int nsteps;
        int eam_flag;
        int gpu_flag;
} SimThreadData;

#endif
