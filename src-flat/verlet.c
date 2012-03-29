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
 * verlet.c: containes routines for verlet time stepping
 *
 * Author: Sriram Swaminarayan
 * Date: 11/4/11
 *
 **/

/**
 * implements the timestepping for verlet and thermostated schemes **/

#include "pmd.h"

static real_t dt = (real_t) 1.0;
static real_t dthalf = (real_t) 0.5;
static real_t dtminushalf = -0.5;

static inline void verlet_velocity(simulation_t *s, const int nStart, const int nTodo) {
  int i;
  real_t dvt;
  real_t *vx,*vy,*vz;
  real_t *fx,*fy,*fz;
  real_t *invMass;

  const real_t cfac = (real_t)(double)0.529/(double)2.418e-17;

  fx = s->f.x.d2;
  fy = s->f.y.d2;
  fz = s->f.z.d2;

  vx = s->p.x.d2;
  vy = s->p.y.d2;
  vz = s->p.z.d2;
  invMass = s->p.w.d2;

  dvt   = dt*cfac;

  for(i=0; i<nTodo; i++) {
    vx[i] = cikNMSub64fp(dvt,fx[i],vx[i]);
    vy[i] = cikNMSub64fp(dvt,fy[i],vy[i]);
    vz[i] = cikNMSub64fp(dvt,fz[i],vz[i]);
  }

  return;
}
static inline void verlet_velocity_half_forward(simulation_t *s, int nStart, int nTodo) {
  real_t dtbackup;
  dtbackup = dt;
  dt = dthalf;
  verlet_velocity(s,nStart,nTodo);
  dt = dtbackup;
  return;
}
static inline void verlet_velocity_half_backward(simulation_t *s, int nStart, int nTodo) {
  real_t dtbackup;
  real_t dtnew;
  dtbackup = dt;
  dt = dtminushalf;
  verlet_velocity(s,nStart,nTodo);
  dt = dtbackup;
  return;
}
static inline void verlet_positions(simulation_t *s, int nStart, int nTodo) {
  int i;
  real_t *x,*y,*z;
  real_t *vx,*vy,*vz;
  real_t *tdx,*tdy,*tdz;
  real_t dxMax = zerosfp64;
  real_t halfnpad;
  const real_t *im = s->p.w.d2;
  const int fjtag = 5;
  const double cfac_2 = (double)0.529/(double)2.418e-17/(double)1822.83;
  const real_t cfac_v= ((real_t){cfac_2,cfac_2});;
  const real_t dxt = cikMul64fp(dt,cfac_v);
  real_t cfac;


  x = s->r.x.d2 + nStart;
  y = s->r.y.d2 + nStart;
  z = s->r.z.d2 + nStart;
  
  vx = s->p.x.d2;
  vy = s->p.y.d2;
  vz = s->p.z.d2;

  halfnpad = cikMul64fp(padNL,halffp64);
  //  pn5("positions 1");
  for(i=0; i<nTodo; i++) {
    cfac = cikMul64fp(dxt,im[i]);
    
    x[i] = cikMAdd64fp(cfac,vx[i],x[i]);
    y[i] = cikMAdd64fp(cfac,vy[i],y[i]);
    z[i] = cikMAdd64fp(cfac,vz[i],z[i]);

  }
  //  pn5("positions 2");
  return;
    
}
static inline void verlet_positions_half_forward(simulation_t *s, int nStart, int nTodo) {
  real_t dtbackup;
  dtbackup = dt;
  dt = dthalf;
  verlet_positions(s,nStart,nTodo);
  dt = dtbackup;
  return;
}
static inline void verlet_positions_half_backward(simulation_t *s, int nStart, int nTodo) {
  real_t dtbackup;
  real_t dtnew;
  dtbackup = dt;
  dt = dtminushalf;
  verlet_positions(s,nStart,nTodo);
  dt = dtbackup;
  return;
}

#define force interact_1
extern double interact_1(int mode,simulation_t *s, int nStart, int nTodo);

void verlet_params(double dtin, double Tin, double frictionin) {
  g_ts_dt = dtin;
  g_ts_ttherm = Tin;
  g_ts_fcoef = frictionin;
  dt = (real_t) {dtin, dtin};
  dthalf = cikMul64fp(halffp64,dt);
  dtminushalf = cikMul64fp(minushalffp64,dt);
  //  PRINTCPUTYPE; printf("dt=%g / %g %g ignored\n",dtin,Tin, frictionin);fflush(stdout);
}

extern void speSimGetXYZ(int tag);
extern void speSimPutXYZ(int tag);
extern void speSimPutMomenta(int tag);

static double e[4] __attribute__ ((aligned(128)));
static double *verlet_one_timestep(simulation_t *s, int nStart, int nTodo) {
  int tag = 30;
  real_t r2_tmp_v, r2_v;

  r2_tmp_v = cikAdd64fp(g_pot->phirho->xn_v,padNL);
  r2_v = cikMul64fp(r2_tmp_v, r2_tmp_v);

  // advance positions
  verlet_positions_half_forward(s,nStart,nTodo);

  // compute new neighbor list
  nlUpdate(g_s_spu,cikExtract32(g_s_spu->isize,0),cikExtract32(g_s_spu->isize,1),r2_v);

  // compute new force
  e[0] = force(-1,s,nStart,nTodo);

  // step velocities by half
  (void) verlet_velocity(s,nStart,nTodo);

  // move positions one more half forward
  verlet_positions_half_forward(s,nStart,nTodo);

  // put positions to main memory
  speSimReloadXYZ(tag);

  e[0] = force(-1,s,nStart,nTodo);

  return e;
}
double *verlet_n_timesteps(int n, simulation_t *s, int nStart, int nTodo) {

  int i;
  int tag = 7;
  real_t r2_tmp_v, r2_v;

  r2_tmp_v = cikAdd64fp(g_pot->phirho->xn_v,padNL);
  r2_v = cikMul64fp(r2_tmp_v, r2_tmp_v);

  // half step forward bty positions
  (void) verlet_positions_half_forward(s,nStart,nTodo);

  for(i=0; i<n; i++) {

    speSimReloadXYZ(tag);
    nlUpdate(g_s_spu,cikExtract32(g_s_spu->isize,0),cikExtract32(g_s_spu->isize,1),r2_v);
    
    // compute new force
    e[0] = force(-1,s,nStart,nTodo);

    // step velocities by a full timestep
    (void) verlet_velocity(s,nStart,nTodo);

    // advance positions
    (void) verlet_positions(s,nStart,nTodo);

  }

  // step back velocities by half
  (void) verlet_positions_half_backward(s,nStart,nTodo);
  speSimReloadXYZ(tag);
  nlUpdate(g_s_spu,cikExtract32(g_s_spu->isize,0),cikExtract32(g_s_spu->isize,1),r2_v);
  // compute new force
  e[0] = force(-1,s,nStart,nTodo);

  speBarrier();

  return e;
}
void verlet_do() {
 ts_n = verlet_n_timesteps;
 ts_one = verlet_one_timestep;
 ts_params = verlet_params;
 return;
}
