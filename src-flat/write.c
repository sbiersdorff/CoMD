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

#include <stdio.h>
#include "pmd.h"

typedef struct myrl_t {
  int box;
  int id;
} myrl_t;

static myrl_t *getRL(simflat_t *s) {
  int *idone = NULL;
  myrl_t *rl = NULL;
  int ib, i;
  rl = (myrl_t*)calloc(s->ntot,sizeof(myrl_t));
  idone = (int*)calloc(s->ntot,sizeof(int));
  for(ib=0; ib<s->nboxes; ib++) {
    int ioff;
    for(ioff=MAXATOMS*ib,i=0; i<s->natoms[ib]; i++,ioff++) {
      int id = s->id[ioff];
      rl[id].id = i;
      rl[id].box = ib;
      if(idone[id]) {
	fprintf(stderr,"_____regot id=%d box %d,atom %d\n",id,ib,i);
      }
      else {
	idone[id]++;
      }
    }
  }
  for(i=0; i<s->ntot; i++) {
    if ( ! idone[i]) {
      fprintf(stderr,"_____didnt get id=%d \n",i);
    }
  }
  free(idone);
  fflush(stderr);
      
  return rl;
}

void printIt(simflat_t *s,FILE *fp) {
  int i;
  for(i=0; i<s->nboxes; i++) {
    int j;
    int ioff;
    int *id;

    for(ioff=i*MAXATOMS,j=0; j<s->natoms[i]; j++,ioff++) {
      if ( s->id[ioff] < 10) {
	fprintf(fp,
		"%02d, %02d, X=(%+020.12e %+020.12e %+020.12e) 1 P=(%+020.12e %+020.12e %+020.12e) F=(%+020.12e %+020.12e %+020.12e)\n",
                i,
		s->id[ioff]+1,
		s->r[ioff][0],s->r[ioff][1],s->r[ioff][2],
		s->p[ioff][0],s->p[ioff][1],s->p[ioff][2],
		s->f[ioff][0],s->f[ioff][1],s->f[ioff][2]
		);
      }
    }
  }
  return;
}

void writeClsman(simflat_t *s, char *fname) {
  myrl_t *rl;
  FILE *fp;
  int i;
  int ti[16];
  double td[16];
  

  fp = fopen(fname,"wb");

  ti[0] = 3*sizeof(int);/*header*/
  ti[1]=s->ntot;/*natoms*/
  ti[2]=0;/*nmove*/
  ti[3]=0;/**/
  ti[4] = 3*sizeof(int);/*header*/
  ti[5]=0;/*header*/
  ti[6]=0;/*header*/
  fwrite(ti,7,sizeof(int),fp);

  rl = getRL(s);

  ti[0] = 3*sizeof(double);/*header*/
  fwrite(ti,1,sizeof(int),fp);

  /* bounds */
  td[0] = (double)s->bounds[0];
  td[1] = (double)s->bounds[1];
  td[2] = (double)s->bounds[2];
  fwrite(td,3,sizeof(double),fp); 

  ti[0] = 3*sizeof(double);/*header*/
  fwrite(ti,1,sizeof(int),fp);
  
  ti[0] = s->ntot*(6*sizeof(double) + 1*sizeof(int));
  fwrite(ti,1,sizeof(int),fp);
  
  for(i=0; i<s->ntot; i++) {
    const int jtype=1;
    int id = MAXATOMS*rl[i].box+rl[i].id;
    int j;
    double r[3],p[3];
    for(j=0; j<3; j++) {
      td[j] = (double)s->dcenter[rl[i].box][j] + (double)s->r[id][j];
      td[j+3] = (double)s->p[id][j];
    }
    fwrite(td,3,sizeof(double),fp);
    fwrite(&jtype,1,sizeof(int),fp);
    fwrite(td+3,3,sizeof(double),fp);
  }
  ti[0] = s->ntot*(6*8 + 1*4);
  fwrite(ti,1,sizeof(int),fp);
  fflush(fp);
  fclose(fp);
  free(rl);
  return;
}
  
