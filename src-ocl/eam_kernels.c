/** Kernels for computing the EAM potential
  Since we can only block on kernel completion, the 3 sweeps done in the original code 
  need to be implemented as 3 separate kernels. 
  Note also that the potential arrays are large enough to require accessing them from 
  global memory.
  Since OpenCL doesn't pick up #include properly, we need to manually switch real_t from 
  float to double in each kernel file individually.

  Note: More careful analysis shows we can consolidate kernels 1 and 2 into a single pass;
  there is a flag PASS_2 in this file which should be set to match the flag PASS_2 in 
  helpers.c. Switching this flag allows you to test that a) the results match and 2) evaluate
  the overhead of adding the extra kernel which only loops over all particles.
 **/

//Initial implementation of the MD code
#define N_MAX_ATOMS 64
#define N_MAX_NEIGHBORS 27
#define PERIODIC 1

#define KERN_DIAG 0
#define USE_SPLINE 0
#define USE_CHEBY 0
#define PASS_2 0 // this should match the setting in helpers.c

#if defined(cl_khr_fp64)  // Khronos extension available?
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)  // AMD extension available?
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

/* CL_REAL_T is set to single or double depending on compile time flags */
typedef CL_REAL_T real_t; 
typedef CL_REAL4_T real3_t;

#if (USE_CHEBY)
// given a list of Chebyshev coefficients c, compute the value at x
// x must be in the range [a, b]
// call signature might change to put it in line with the regular EAM call
real_t chebev(
        real_t a, 
        real_t b, 
        __global const real_t *c,
        int m, 
        real_t x) 
{
   real_t d, dd, sv, y, y2;
   real_t ch;
   int j;
   if ((x-a)*(x-b) > 0.0) {
      printf("x not in range in chebev, %f\n", x);
   }
   d=0.0;
   dd=0.0;
   y=(2.0*x-a-b)/(b-a);
   y2=2.0*y;
   for(j=m-1;j>0;j--) {
      sv=d;
      d=y2*d-dd+c[j];
      dd=sv;
   }
   ch=y*d-dd+0.5*c[0];
   return ch;
}
#endif

#if (USE_SPLINE)
void PhiSpline(real_t r, real_t *f, real_t *df)
{
    // values for copper
    /* 
    // original values from 87 paper
    real_t a_k[6] = { 29.059214 , -140.05681 , 130.07331 , -17.48135 , 31.82546 , 71.58749};
    real_t r_k[6] = { 1.2247449 , 1.1547054 , 1.1180065 , 1.0000000 , 0.8660254 , 0.7071068};
    */
    // new values for smoother potential
    real_t a_k[6] = {61.73525861, -108.18467800, 57.00053948,-12.88796578, 39.16381901, 0.0};
    real_t r_k[6] = {1.225, 1.202, 1.154, 1.050, 0.866, 0.707};
    real_t az = 1.0*3.615;
    real_t az3 = az*az*az;
    real_t az2 = az*az;
    int k;
    // set output values to zero
    *f=0.0;
    *df=0.0;
    //r = 1.0; //3.615;
    r = r/3.615;
    //printf("r = %e\n", r);
    // sum over all coefficients
    for (k=0;k<6;k++) {
	r_k[k] = r_k[k]*3.615;
	if (r < r_k[k]) {
	    *f += (r_k[k]-r)*(r_k[k]-r)*(r_k[k]-r)*a_k[k]/az3;
	    *df -= 3.0*(r_k[k]-r)*(r_k[k]-r)*a_k[k]/az2;
	}
    }
}

void RhoSpline(real_t r, real_t *f, real_t *df)
{
    // values for copper
    /*
    // original values from 87 paper
    real_t R_k[2] = { 1.2247449 , 1.0000000 };
    real_t A_k[2] = { 9.806694 , 16.774638 };
    */
    // new values for smoother potential
    real_t R_k[2] = { 1.225, 0.990 };
    real_t A_k[2] = { 10.03718305, 17.06363299 };
    real_t az = 1.0*3.615;
    real_t az3 = az*az*az;
    real_t az2 = az*az;
    int k;
    // set output values to zero
    *f=0.0;
    *df=0.0;
    //r = 1.0; //3.615;
    r = r/3.615;
    // sum over all coefficients
    for (k=0;k<2;k++) {
	R_k[k] = R_k[k]*3.615;
	if (r < R_k[k]) {
	    *f += (R_k[k]-r)*(R_k[k]-r)*(R_k[k]-r)*A_k[k]/az3;
	    *df -= 3.0*(R_k[k]-r)*(R_k[k]-r)*A_k[k]/az2;
	}
    }
}

void FSpline(real_t rho, real_t *f, real_t *df)
{
    *f = -1.0*sqrt(rho);
    *df = -0.5/(*f);
}

#endif

void eamInterpolateDeriv(real_t r,
	__global const real_t* values,
	int n_values,
	real_t *value1, 
	real_t *f1)
{
    int i1;
    int i;
    real_t gi, gi1;

    // extract values from potential 'struct'
    real_t x0 = values[n_values+3];
    real_t xn = values[n_values+4];
    real_t invDx = values[n_values+5];

    // identical to Sriram's loop in eam.c
    if ( r<x0) r = x0;
    else if (r>xn) r = xn;

    r = (r-x0)*(invDx) ;
    i1 = (int)floor(r);

    /* reset r to fractional distance */
    r = r - floor(r);

    // all indices shifted up by one compared to the original code
    gi  = values[i1+2] - values[i1];
    gi1 = values[i1+3] - values[i1+1];

    // Note the shift removes [i1-1] as a possibility
    // values[i1-1] is guaranteed(?) inbounds because 
    // a->x0 = x0 + (xn-x0)/(double)n; 
    // appears in allocPotentialArray
    *value1 = values[i1+1] + 0.5*r*(
	    r*(values[i1+2]+ values[i1] -2.0*values[i1+1]) +
	    gi
	    );
    if(i1<=1) 
	*f1 = 0.0;
    else 
	*f1 = 0.5*(gi + r*(gi1-gi))*invDx;

    return;

}

// Simple version without local blocking to check for correctness
__kernel void EAM_Force_1(
	__global real_t* x_pos,
	__global real_t* y_pos,
	__global real_t* z_pos,

	__global real_t* fx,
	__global real_t* fy,
	__global real_t* fz,

	__global real_t* energy,
	__global real_t* rho,
	__global real_t* rhobar,

	__global real_t* dcx,
	__global real_t* dcy,
	__global real_t* dcz,

	__global real_t* bounds,
	__global const int* neighbor_list,
	__global const int* n_neighbors,
	__global const int* natoms,

	__global const int* n_values,

	__global const real_t* phi_pot, // the potentials are passed in as real arrays: x0, xn, invDx, values[?]
	__global const real_t* rho_pot,
	__global const real_t* F_pot,

	const real_t cutoff) 
{
#if USE_SPLINE
    // values for copper
    real_t a_k[6] = { 29.059214 , -140.05681 , 130.07331 , -17.48135 , 31.82546 , 71.58749};
    real_t r_k[6] = { 1.2247449 , 1.1547054 , 1.1180065 , 1.0000000 , 0.8660254 , 0.7071068};

    real_t R_k[2] = { 1.2247449 , 1.0000000 };
    real_t A_k[2] = { 9.806694 , 16.774638 };
#endif

    int iatom = get_global_id(0);
    int ibox = get_global_id(1);

    real_t dx, dy, dz;
    real_t r, r2, r6;
    real_t fr, e_i;
    real_t rho_i;

    real_t dxbox, dybox, dzbox;

    // accumulate local force value
    real_t fx_i, fy_i, fz_i;

    real_t rcut = cutoff;
    real_t r2cut = rcut*rcut;
    real_t rhoTmp;
    real_t phiTmp;
    real_t dPhi, dRho;
    real_t fi, fiprime;

    int i;
    int j, j_local;

    int i_offset;
    int i_particle;

    // zero out forces on particles
    i_offset = ibox*N_MAX_ATOMS;
    i_particle = i_offset + iatom;

    fx_i = 0.0;
    fy_i = 0.0;
    fz_i = 0.0;

    rho_i = 0.0;

    e_i = 0.0;


    if (iatom < natoms[ibox]) { // each thread executes on a single atom in the box

#if(KERN_DIAG > 0) 
	printf("i = %d, %f, %f, %f\n", i_particle, x_pos[i_particle], y_pos[i_particle], z_pos[i_particle]);

	printf("ibox = %d, n_neighbors = %d\n", ibox, n_neighbors[ibox]);
#endif

	for (j = 0; j<n_neighbors[ibox]; j++) { // loop over neighbor cells
	    int jbox = neighbor_list[ibox*N_MAX_NEIGHBORS + j];
	    int j_offset = jbox*N_MAX_ATOMS;

	    // compute box center offsets
	    dxbox = dcx[ibox] - dcx[jbox];
	    dybox = dcy[ibox] - dcy[jbox];
	    dzbox = dcz[ibox] - dcz[jbox];

	    // correct for periodic 
	    if(PERIODIC) {
		if (dxbox<-0.5*bounds[0]) dxbox += bounds[0];
		else if (dxbox > 0.5*bounds[0] ) dxbox -= bounds[0];
		if (dybox<-0.5*bounds[1]) dybox += bounds[1];
		else if (dybox > 0.5*bounds[1] ) dybox -= bounds[1];
		if (dzbox<-0.5*bounds[2]) dzbox += bounds[2];
		else if (dzbox > 0.5*bounds[2] ) dzbox -= bounds[2];
	    }

	    //printf("dxbox, dybox, dzbox = %f, %f, %f\n", dxbox, dybox, dzbox);

	    int jatom;
	    for (jatom = 0; jatom<natoms[jbox]; jatom++) { // loop over all groups in neighbor cell 

		int j_particle = j_offset + jatom; // global offset of particle

		dx = dxbox + x_pos[i_particle] - x_pos[j_particle];
		dy = dybox + y_pos[i_particle] - y_pos[j_particle];
		dz = dzbox + z_pos[i_particle] - z_pos[j_particle];

#if(KERN_DIAG > 0) 
		printf("dx, dy, dz = %f, %f, %f\n", dx, dy, dz);
		printf("i = %d, j = %d, %f, %f, %f\n", i_particle, j_particle, x_pos[j_particle], y_pos[j_particle], z_pos[j_particle]);
#endif

		r2 = dx*dx + dy*dy + dz*dz;

		if ( r2 <= r2cut && r2 > 0.0) { // no divide by zero

#if(KERN_DIAG > 0) 
		    printf("%d, %d, %f\n", i_particle, j_particle, r2);
		    printf("r2, r2cut = %f, %f\n", r2, r2cut);
#endif

		    r = sqrt(r2);
		    //r = 1.0;

#if(USE_SPLINE)
		    //r = r*r_k[0]/cutoff;
		    //r = r/3.615;
		    PhiSpline(r, &phiTmp, &dPhi);
		    RhoSpline(r, &rhoTmp, &dRho);
#else
		    eamInterpolateDeriv(r, phi_pot, n_values[0], &phiTmp, &dPhi);
		    eamInterpolateDeriv(r, rho_pot, n_values[1], &rhoTmp, &dRho);
#endif

#if(KERN_DIAG > 0) 
		    printf("%d, %d, %f\n", i_particle, j_particle, r2);
		    printf("i_particle = %d, j_particle = %d, phiTmp = %f, dPhi = %f\n", i_particle, j_particle, phiTmp, dPhi);
		    printf("i_particle = %d, j_particle = %d, rhoTmp = %f, dRho = %f\n", i_particle, j_particle, rhoTmp, dRho);
#endif

#if(USE_SPLINE)
		    fx_i += (dRho+dPhi)*dx/r;
		    fy_i += (dRho+dPhi)*dy/r;
		    fz_i += (dRho+dPhi)*dz/r;
#else
		    fx_i += dPhi*dx/r;
		    fy_i += dPhi*dy/r;
		    fz_i += dPhi*dz/r;
#endif

		    e_i += phiTmp;

		    rho_i += rhoTmp;

		} else {
		}


	    } // loop over all atoms
	} // loop over neighbor cells

	fx[i_particle] = fx_i;
	fy[i_particle] = fy_i;
	fz[i_particle] = fz_i;

	// since we loop over all particles, each particle contributes 1/2 the pair energy to the total
	//energy[i_particle] = e_i*0.5;

	rho[i_particle] = rho_i;

	 // we can actually include the Force_2 kernel here, to save some time!
	 // skip the next 4 lines if PASS_2 = 1 in helpers.c
#if(PASS_2 == 0)
	eamInterpolateDeriv(rho_i,F_pot,n_values[2],&fi,&fiprime);
	rhobar[i_particle] = fiprime; // update rhoprime 
	//update energy terms 
	energy[i_particle] = e_i*0.5 + fi;
#else
	energy[i_particle] = e_i*0.5;
#endif
    } else { // zero out the energy of the other particles for safety
	energy[i_particle] = 0.0;
    }
}

__kernel void EAM_Force_2(
	__global real_t* rhobar,
	__global real_t* energy,
	__global real_t* rho,
	__global const int* natoms,
	__global const real_t* F_pot,
	__global const int* n_values
	)
{

    int iatom = get_global_id(0);
    int ibox = get_global_id(1);

    real_t fi, fiprime;


    int i;

    int i_offset;
    int i_particle;

    i_offset = ibox*N_MAX_ATOMS;
    i_particle = i_offset + iatom;

    /*
    // local copy of F potential
    real_t F_local[n_values+3];

    // load values into local potentials
    for (i=0;i<n_values+3;i++) {
    F_local[i] = F_pot[i];
    }
     */

    if (iatom < natoms[ibox]) { // each thread executes on a single atom in the box


	//printf("ibox = %d, iatom = %d\n", ibox, iatom);

#if(USE_SPLINE)
	FSpline(rho[i_particle], &fi, &fiprime);
#else
	eamInterpolateDeriv(rho[i_particle],F_pot,n_values[2],&fi,&fiprime);
#endif
	rhobar[i_particle] = fiprime; // update rhoprime 
	//update energy terms 
	energy[i_particle] += fi;

    }

}

__kernel void EAM_Force_3(
	__global real_t* x_pos,
	__global real_t* y_pos,
	__global real_t* z_pos,

	__global real_t* fx,
	__global real_t* fy,
	__global real_t* fz,

	__global real_t* fi,

	__global real_t* dcx,
	__global real_t* dcy,
	__global real_t* dcz,

	__global const real_t* bounds,
	__global const int* neighbor_list,
	__global const int* n_neighbors,
	__global const int* natoms,
	__global const int* n_values,

	__global const real_t* rho_pot,
	const real_t cutoff) 
{

    int iatom = get_global_id(0);
    int ibox = get_global_id(1);

    real_t dx, dy, dz;
    real_t r, r2;

    real_t dxbox, dybox, dzbox;

    // accumulate local force value
    real_t fx_i, fy_i, fz_i;

    real_t rcut = cutoff;
    real_t r2cut = rcut*rcut;
    real_t rhoTmp, dRho;

    int i;
    int j, j_local;
    int jbox, jatom;

    real_t rTmp, rhoijprime;

    // global offset of local thread
    int i_offset = ibox*N_MAX_ATOMS;
    int i_particle = i_offset + iatom;

    if (iatom < natoms[ibox]) { // each thread executes on a single atom in the box

	// zero out forces on particles
	fx_i = fx[i_particle];
	fy_i = fy[i_particle];
	fz_i = fz[i_particle];

#if(KERN_DIAG > 0) 

	printf("i = %d, %f, %f, %f\n", i_particle, x_pos[i_particle], y_pos[i_particle], z_pos[i_particle]);

	printf("ibox = %d, n_neighbors = %d\n", ibox, n_neighbors[ibox]);
#endif

	for (j = 0; j<n_neighbors[ibox]; j++) { // loop over neighbor cells
	    int jbox = neighbor_list[ibox*N_MAX_NEIGHBORS + j];
	    int j_offset = jbox*N_MAX_ATOMS;

	    // compute box center offsets
	    dxbox = dcx[ibox] - dcx[jbox];
	    dybox = dcy[ibox] - dcy[jbox];
	    dzbox = dcz[ibox] - dcz[jbox];

	    // correct for periodic 
	    if(PERIODIC) {
		if (dxbox<-0.5*bounds[0]) dxbox += bounds[0];
		else if (dxbox > 0.5*bounds[0] ) dxbox -= bounds[0];
		if (dybox<-0.5*bounds[1]) dybox += bounds[1];
		else if (dybox > 0.5*bounds[1] ) dybox -= bounds[1];
		if (dzbox<-0.5*bounds[2]) dzbox += bounds[2];
		else if (dzbox > 0.5*bounds[2] ) dzbox -= bounds[2];
	    }

	    // printf("dxbox, dybox, dzbox = %f, %f, %f\n", dxbox, dybox, dzbox);

	    for (jatom = 0; jatom<natoms[jbox]; jatom++) { // loop over all groups in neighbor cell 

		int j_particle = j_offset + jatom; // global offset of particle

		/*
		   dx = x_pos[i_particle] - x_pos[j_particle] + dxbox;
		   dy = y_pos[i_particle] - y_pos[j_particle] + dybox;
		   dz = z_pos[i_particle] - z_pos[j_particle] + dzbox;
		 */

		dx = dxbox + x_pos[i_particle] - x_pos[j_particle];
		dy = dybox + y_pos[i_particle] - y_pos[j_particle];
		dz = dzbox + z_pos[i_particle] - z_pos[j_particle];

#if(KERN_DIAG > 0) 
		printf("dx, dy, dz = %f, %f, %f\n", dx, dy, dz);
		printf("i = %d, j = %d, %f, %f, %f\n", i_particle, j_particle, x_pos[j_particle], y_pos[j_particle], z_pos[j_particle]);
#endif

		r2 = dx*dx + dy*dy + dz*dz;

		if ( r2 < r2cut && r2 > 0.0) { // no divide by zero

#if(KERN_DIAG > 0) 
		    printf("r2, r2cut = %f, %f\n", r2, r2cut);
		    printf("%d, %d, %f\n", i_particle, j_particle, r2);
#endif

		    r = sqrt(r2);

#if (USE_SPLINE)
		    RhoSpline(r, &rhoTmp, &dRho);
#else
		    eamInterpolateDeriv(r, rho_pot, n_values[1], &rhoTmp, &dRho);
#endif
		    rhoijprime = dRho;

#if(KERN_DIAG > 0) 
		    printf("%d, %d, %f\n", i_particle, j_particle, r2);
#endif

		    rTmp = (fi[i_particle]+fi[j_particle])*rhoijprime/r;

		    fx_i += (fi[i_particle]+fi[j_particle])*rhoijprime*dx/r;
		    fy_i += (fi[i_particle]+fi[j_particle])*rhoijprime*dy/r;
		    fz_i += (fi[i_particle]+fi[j_particle])*rhoijprime*dz/r;
		    /*
		       fx_i += rTmp*dx;
		       fy_i += rTmp*dy;
		       fz_i += rTmp*dz;
		     */

		} else {
		}


	    } // loop over all atoms in jbox
	} // loop over neighbor cells

#if (USE_SPLINE)
#else
	fx[i_particle] = fx_i;
	fy[i_particle] = fy_i;
	fz[i_particle] = fz_i;
#endif

    } // loop over all atoms in ibox
}

// AoS Versions
// Simple version without local blocking to check for correctness
__kernel void EAM_Force_1_AoS(
	__global real3_t* pos,

	__global real3_t* f,

	__global real_t* energy,
	__global real_t* rho,
	__global real_t* rhobar,

	__global real3_t* dc,

	__global real3_t* bounds,
	__global const int* neighbor_list,
	__global const int* n_neighbors,
	__global const int* natoms,

	__global const int* n_values,

	__global const real_t* phi_pot, // the potentials are passed in as real arrays: x0, xn, invDx, values[?]
	__global const real_t* rho_pot,
	__global const real_t* F_pot,

	const real_t cutoff)
{

    int iatom = get_global_id(0);
    int ibox = get_global_id(1);

    real_t r, r2, r6;
    real_t fr, e_i;
    real_t rho_i;

    real3_t dr;
    real3_t drbox;
    real3_t f_i;

    real_t rcut = cutoff;
    real_t r2cut = rcut*rcut;
    real_t rhoTmp;
    real_t phiTmp;
    real_t dPhi, dRho;
    real_t fi, fiprime;

    int i;
    int j, j_local;

    int i_offset;
    int i_particle;

    i_offset = ibox*N_MAX_ATOMS;
    i_particle = i_offset + iatom;

    // zero out forces on particles
    f_i.x = 0.0;
    f_i.y = 0.0;
    f_i.z = 0.0;

    rho_i = 0.0;

    e_i = 0.0;


    if (iatom < natoms[ibox]) { // each thread executes on a single atom in the box

#if(KERN_DIAG > 0) 
	printf("i = %d, %f, %f, %f\n", i_particle, pos[i_particle].x, pos[i_particle].y, pos[i_particle].z);

	printf("ibox = %d, n_neighbors = %d\n", ibox, n_neighbors[ibox]);
#endif

	for (j = 0; j<n_neighbors[ibox]; j++) { // loop over neighbor cells
	    int jbox = neighbor_list[ibox*N_MAX_NEIGHBORS + j];
	    int j_offset = jbox*N_MAX_ATOMS;

	    // compute box center offsets
	    drbox = dc[ibox] - dc[jbox];

	    // correct for periodic 
	    if(PERIODIC) {
		if (drbox.x<-0.5*bounds[0].x) drbox.x += bounds[0].x;
		else if (drbox.x > 0.5*bounds[0].x ) drbox.x -= bounds[0].x;
		if (drbox.y<-0.5*bounds[0].y) drbox.y += bounds[0].y;
		else if (drbox.y > 0.5*bounds[0].y ) drbox.y -= bounds[0].y;
		if (drbox.z<-0.5*bounds[0].z) drbox.z += bounds[0].z;
		else if (drbox.z > 0.5*bounds[0].z ) drbox.z -= bounds[0].z;
	    }

#if(KERN_DIAG > 0) 
	    printf("dxbox, dybox, dzbox = %f, %f, %f\n", dxbox, dybox, dzbox);
#endif

	    int jatom;
	    for (jatom = 0; jatom<natoms[jbox]; jatom++) { // loop over all groups in neighbor cell 

		int j_particle = j_offset + jatom; // global offset of particle

		dr = pos[i_particle] - pos[j_particle] + drbox;

#if(KERN_DIAG > 0) 
		printf("dx, dy, dz = %f, %f, %f\n", dr.x, dr.y, dr.z);
		printf("i = %d, j = %d, %f, %f, %f\n", 
		i_particle, j_particle, pos[j_particle].x, pos[j_particle].y, pos[j_particle].z);
#endif

		r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;

		if ( r2 <= r2cut && r2 > 0.0) { // no divide by zero

#if(KERN_DIAG > 0) 
		    printf("%d, %d, %f\n", i_particle, j_particle, r2);
		    printf("r2, r2cut = %f, %f\n", r2, r2cut);
#endif

		    r = sqrt(r2);

		    eamInterpolateDeriv(r, phi_pot, n_values[0], &phiTmp, &dPhi);
		    eamInterpolateDeriv(r, rho_pot, n_values[1], &rhoTmp, &dRho);

#if(KERN_DIAG > 0) 
		    printf("%d, %d, %f\n", i_particle, j_particle, r2);
		    printf("i_particle = %d, j_particle = %d, phiTmp = %f, dPhi = %f\n", i_particle, j_particle, phiTmp, dPhi);
		    printf("i_particle = %d, j_particle = %d, rhoTmp = %f, dRho = %f\n", i_particle, j_particle, rhoTmp, dRho);
#endif

		    f_i.x += dPhi*dr.x/r;
		    f_i.y += dPhi*dr.y/r;
		    f_i.z += dPhi*dr.z/r;

		    e_i += phiTmp;

		    rho_i += rhoTmp;

		} else {
		}


	    } // loop over all atoms
	} // loop over neighbor cells

	f[i_particle] = f_i;

	// since we loop over all particles, each particle contributes 1/2 the pair energy to the total
	//energy[i_particle] = e_i*0.5;

	rho[i_particle] = rho_i;

#if(PASS_2 == 0)
	eamInterpolateDeriv(rho_i,F_pot,n_values[2],&fi,&fiprime);
	rhobar[i_particle] = fiprime; // update rhoprime 
	//update energy terms 
	energy[i_particle] = e_i*0.5 + fi;
#else
	energy[i_particle] = e_i*0.5;
#endif
    }
}

__kernel void EAM_Force_2_AoS(
	__global real_t* rhobar,
	__global real_t* energy,
	__global real_t* rho,
	__global const int* natoms,
	__global const real_t* F_pot,
	__global const int* n_values
	)
{

    int iatom = get_global_id(0);
    int ibox = get_global_id(1);

    real_t fi, fiprime;


    int i;

    int i_offset;
    int i_particle;

    i_offset = ibox*N_MAX_ATOMS;
    i_particle = i_offset + iatom;

    if (iatom < natoms[ibox]) { // each thread executes on a single atom in the box


#if(KERN_DIAG > 0) 
	printf("ibox = %d, iatom = %d\n", ibox, iatom);
#endif

	eamInterpolateDeriv(rho[i_particle],F_pot,n_values[2],&fi,&fiprime);
	rhobar[i_particle] = fiprime; // update rhoprime 
	//update energy terms 
	energy[i_particle] += fi;

    }

}
__kernel void EAM_Force_3_AoS(
	__global real3_t* pos,

	__global real3_t* f,

	__global real_t* fi,

	__global real3_t* dc,

	__global const real3_t* bounds,
	__global const int* neighbor_list,
	__global const int* n_neighbors,
	__global const int* natoms,

	__global const int* n_values,
	__global const real_t* rho_pot,
	const real_t cutoff) 
{

    int iatom = get_global_id(0);
    int ibox = get_global_id(1);

    real3_t dr;
    real_t r, r2;

    real3_t drbox;

    // accumulate local force value
    real3_t f_i;

    real_t rcut = cutoff;
    real_t r2cut = rcut*rcut;
    real_t rhoTmp, dRho;

    int i;
    int j, j_local;
    int jbox, jatom;

    real_t rTmp, rhoijprime;

    // global offset of local thread
    int i_offset = ibox*N_MAX_ATOMS;
    int i_particle = i_offset + iatom;

    if (iatom < natoms[ibox]) { // each thread executes on a single atom in the box

	// zero out forces on particles
	f_i.x = f[i_particle].x;
	f_i.y = f[i_particle].y;
	f_i.z = f[i_particle].z;

#if(KERN_DIAG > 0) 
	printf("i = %d, %f, %f, %f\n", i_particle, pos[i_particle].x, pos[i_particle].y, pos[i_particle].z);
	printf("ibox = %d, n_neighbors = %d\n", ibox, n_neighbors[ibox]);
#endif

	for (j = 0; j<n_neighbors[ibox]; j++) { // loop over neighbor cells
	    int jbox = neighbor_list[ibox*N_MAX_NEIGHBORS + j];
	    int j_offset = jbox*N_MAX_ATOMS;

	    // compute box center offsets
	    drbox = dc[ibox] - dc[jbox];

	    // correct for periodic 
	    if(PERIODIC) {
		if (drbox.x<-0.5*bounds[0].x) drbox.x += bounds[0].x;
		else if (drbox.x > 0.5*bounds[0].x ) drbox.x -= bounds[0].x;
		if (drbox.y<-0.5*bounds[0].y) drbox.y += bounds[0].y;
		else if (drbox.y > 0.5*bounds[0].y ) drbox.y -= bounds[0].y;
		if (drbox.z<-0.5*bounds[0].z) drbox.z += bounds[0].z;
		else if (drbox.z > 0.5*bounds[0].z ) drbox.z -= bounds[0].z;
	    }

	    // printf("dxbox, dybox, dzbox = %f, %f, %f\n", drbox.x, drbox.y, drbox.z);

	    for (jatom = 0; jatom<natoms[jbox]; jatom++) { // loop over all groups in neighbor cell 

		int j_particle = j_offset + jatom; // global offset of particle

		dr = pos[i_particle] - pos[j_particle] + drbox;

#if(KERN_DIAG > 0) 
		printf("dx, dy, dz = %f, %f, %f\n", dr.x, dr.y, dr.z);
		printf("i = %d, j = %d, %f, %f, %f\n", i_particle, j_particle, pos[j_particle].x, pos[j_particle].y, pos[j_particle].z);
#endif

		r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;

		if ( r2 < r2cut && r2 > 0.0) { // no divide by zero

#if(KERN_DIAG > 0) 
		    printf("r2, rcut = %f, %f\n", r2, rcut);
		    printf("%d, %d, %f\n", i_particle, j_particle, r2);
#endif

		    r = sqrt(r2);

		    eamInterpolateDeriv(r, rho_pot, n_values[1], &rhoTmp, &dRho);
		    rhoijprime = dRho;

#if(KERN_DIAG > 0) 
		    printf("%d, %d, %f\n", i_particle, j_particle, r2);
#endif

		    rTmp = (fi[i_particle]+fi[j_particle])*rhoijprime/r;

		    /*
		    f_i.x += (fi[i_particle]+fi[j_particle])*rhoijprime*dr.x/r;
		    f_i.y += (fi[i_particle]+fi[j_particle])*rhoijprime*dr.y/r;
		    f_i.z += (fi[i_particle]+fi[j_particle])*rhoijprime*dr.z/r;
		     */
		       f_i.x += rTmp*dr.x;
		       f_i.y += rTmp*dr.y;
		       f_i.z += rTmp*dr.z;

		} else {
		}


	    } // loop over all atoms in jbox
	} // loop over neighbor cells

	f[i_particle] = f_i;

    } // loop over all atoms in ibox
}

