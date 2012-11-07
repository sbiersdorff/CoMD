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


#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <float.h>
#include <time.h>

#include "pmd.h"
#include "utility.h"
#include "ic_fcc.h"

#define DOICTIMERS 4

simflat_t *create_fcc_lattice(command_t cmd, struct pmd_base_potential_t *pot) {
    /**
     * Creates an fcc lattice with nx * ny * nz unit cells and lattice constant lat
     *
     **/
    int nx = cmd.nx;
    int ny = cmd.ny;
    int nz = cmd.nz;
    real_t lat = cmd.lat;
    real_t defgrad = cmd.defgrad;
    real_t boxfactor = cmd.bf;

    simflat_t *s = NULL;
    int     i, j, k, n, itype, natoms;
    real_t  x, y, z, halflat;
#ifdef DOICTIMERS  
    clock_t start,old,count=0;
    real2_t fx,fy,fz, px,py,pz;
    fx=fy=fz=px=py=pz=0.0;
#endif

#ifdef DOICTIMERS  
    start = clock();
#endif
    s = blankSimulation(pot);
    if ( ! s ) simulationAbort(-50,(char *) "Unable to create Simulation data structure");

    /* Optional simulation comment (blank for now) */
    s->comment = (char*)suAlignedCalloc(1024*sizeof(char));

    natoms = 4*nx*ny*nz;
    halflat = lat / 2.0;

    /* periodic  boundaries  */
    memset(s->bounds,0,sizeof(real3));
    memset(s->boxsize,0,sizeof(real3));
    // stretch in x direction
    s->defgrad = defgrad;
    s->bounds[0] = nx * lat * defgrad;
    s->bounds[1] = ny * lat;
    s->bounds[2] = nz * lat;

    s->bf = boxfactor;

#if DOICTIMERS  >2
    old = clock();
#endif
    allocDomains(s);

#if DOICTIMERS  >2
    count = clock()-old+count;
#endif
    i = j = k = n = 0;
    z = lat / 4.0;
    while (z < s->bounds[2]) {
	y = lat / 4.0;
	while (y < s->bounds[1]) {
	    x = lat * defgrad / 4.0;
	    while (x < s->bounds[0]) {
		if ((i+j+k) % 2)
		    putAtomInBox(s,n++,1,1,s->pot->mass,x,y,z,px,py,pz,fx,fy,fz);
		x += halflat * defgrad;
		i++;
	    }
	    y += halflat;
	    j++;
	}
	z += halflat;
	k++;
    }

#ifdef DOICTIMERS  
    old = clock();
#endif
    if(DOICTIMERS) printf("\n    ---- FCC initial condition took %.2gs for %d atoms\n\n",
	    (float)(old-start)/(float)(CLOCKS_PER_SEC), n);

    return s;
}



