/** Kernels for computing the EAM potential
  Since we can only block on kernel completion, the 3 sweeps done in the original code 
  need to be implemented as 3 separate kernels. 
  Note also that the potential arrays are large enough to require accessing them from 
  global memory.
  Since OpenCL doesn't pick up #include properly, we need to manually switch real_t from 
  float to double in each kernel file individually.
  **/ 

#define N_MAX_ATOMS 64
#define N_MAX_NEIGHBORS 27
#define PERIODIC 1

#define ALLOW_PRINTF 0

#if defined(cl_khr_fp64)  // Khronos extension available?
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)  // AMD extension available?
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

//Initial implementation of the MD code
/* CL_REAL_T is set to single or double depending on compile time flags */
typedef CL_REAL_T real_t;

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
    real_t x0 = values[n_values+0];
    real_t xn = values[n_values+1];
    real_t invDx = values[n_values+2];

    // identical to Sriram's loop in eam.c
    if ( r<x0) r = x0;
    else if (r>xn) r = xn;

    r = (r-x0)*(invDx) ;
    i1 = (int)floor(r);

    /* reset r to fractional distance */
    r = r - floor(r);

    gi  = values[i1+1] - values[i1-1];
    gi1 = values[i1+2] - values[i1];


    // values[i1-1] is guaranteed(?) inbounds because 
    // a->x0 = x0 + (xn-x0)/(double)n; 
    // appears in allocPotentialArray
    *value1 = values[i1] + 0.5*r*(
	    r*(values[i1+1]+ values[i1-1] -2.0*values[i1]) +
	    gi
	    );
    if(i1<=0) *f1 = 0.0;
    else *f1 = 0.5*(gi + r*(gi1-gi))*invDx;

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

	const real_t cutoff) 
{

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

    int i;
    int j, j_local;
    int jbox, jatom;

    int i_offset, j_offset;
    int i_particle, j_particle;

    // zero out forces on particles
    i_offset = ibox*N_MAX_ATOMS;
    i_particle = i_offset + iatom;

    fx_i = 0.0;
    fy_i = 0.0;
    fz_i = 0.0;

    rho_i = 0.0;

    e_i = 0.0;


    if (iatom < natoms[ibox]) { // each thread executes on a single atom in the box

#if(ALLOW_PRINTF) 
	printf("i = %d, %f, %f, %f\n", i_particle, x_pos[i_particle], y_pos[i_particle], z_pos[i_particle]);

	printf("ibox = %d, n_neighbors = %d\n", ibox, n_neighbors[ibox]);
#endif

	for (j = 0; j<n_neighbors[ibox]; j++) { // loop over neighbor cells
	    int jbox = neighbor_list[ibox*N_MAX_NEIGHBORS + j];
	    j_offset = jbox*N_MAX_ATOMS;

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

	    for (jatom = 0; jatom<natoms[jbox]; jatom++) { // loop over all groups in neighbor cell 

		j_particle = j_offset + jatom; // global offset of particle

		dx = x_pos[i_particle] - x_pos[j_particle] + dxbox;;
		dy = y_pos[i_particle] - y_pos[j_particle] + dybox;;
		dz = z_pos[i_particle] - z_pos[j_particle] + dzbox;;

#if(ALLOW_PRINTF) 
		//printf("dx, dy, dz = %f, %f, %f\n", dx, dy, dz);
		//printf("i = %d, j = %d, %f, %f, %f\n", i_particle, j_particle, x_pos[j_particle], y_pos[j_particle], z_pos[j_particle]);
#endif

		r2 = dx*dx + dy*dy + dz*dz;

		if ( r2 <= r2cut && r2 > 0.0) { // no divide by zero

#if(ALLOW_PRINTF) 
		    printf("%d, %d, %f\n", i_particle, j_particle, r2);
		    //printf("r2, rcut = %f, %f\n", r2, rcut);
#endif
		    //printf("%d, %d, %f\n", i_particle, j_particle, r2);

		    r = sqrt(r2);

		    eamInterpolateDeriv(r, phi_pot, n_values[0], &phiTmp, &dPhi);
		    eamInterpolateDeriv(r, rho_pot, n_values[1], &rhoTmp, &dRho);

#if(ALLOW_PRINTF) 
		    //printf("%d, %d, %f\n", i_particle, j_particle, r2);
		    //printf("i_particle = %d, j_particle = %d, i_b = %d, r6 = %f\n", i_particle, j_particle, i_b, r6);
#endif

		    fx_i += dPhi*dx/r;
		    fy_i += dPhi*dy/r;
		    fz_i += dPhi*dz/r;

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
	energy[i_particle] = e_i*0.5;

	rho[i_particle] = rho_i;

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

    real_t dx, dy, dz;
    real_t r, r2, r6;
    real_t fr, e;

    real_t dxbox, dybox, dzbox;

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

	eamInterpolateDeriv(rho[i_particle],F_pot,n_values[2],&fi,&fiprime);
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

	__global real_t* energy,
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

    int i;
    int j, j_local;
    int jbox, jatom;

    int i_offset, j_offset;
    int i_particle, j_particle;

    real_t rTmp, rhoijprime;

    // zero out forces on particles
    fx_i = 0.0;
    fy_i = 0.0;
    fz_i = 0.0;

    // global offset of local thread
    i_offset = ibox*N_MAX_ATOMS;
    i_particle = i_offset + iatom;

    if (iatom < natoms[ibox]) { // each thread executes on a single atom in the box

#if(ALLOW_PRINTF) 
	printf("i = %d, %f, %f, %f\n", i_particle, x_pos[i_particle], y_pos[i_particle], z_pos[i_particle]);

	printf("ibox = %d, n_neighbors = %d\n", ibox, n_neighbors[ibox]);
#endif

	for (j = 0; j<n_neighbors[ibox]; j++) { // loop over neighbor cells
	    int jbox = neighbor_list[ibox*N_MAX_NEIGHBORS + j];
	    j_offset = jbox*N_MAX_ATOMS;

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

		j_particle = j_offset + jatom; // global offset of particle

		dx = x_pos[i_particle] - x_pos[j_particle] + dxbox;;
		dy = y_pos[i_particle] - y_pos[j_particle] + dybox;;
		dz = z_pos[i_particle] - z_pos[j_particle] + dzbox;;

#if(ALLOW_PRINTF) 
		//printf("dx, dy, dz = %f, %f, %f\n", dx, dy, dz);
		//printf("i = %d, j = %d, %f, %f, %f\n", i_particle, j_particle, x_pos[j_particle], y_pos[j_particle], z_pos[j_particle]);
#endif

		r2 = dx*dx + dy*dy + dz*dz;

		if ( r2 < r2cut && r2 > 0.0) { // no divide by zero

#if(ALLOW_PRINTF) 
		    //printf("r2, rcut = %f, %f\n", r2, rcut);
#endif
		    //printf("%d, %d, %f\n", i_particle, j_particle, r2);

		    r = sqrt(r2);

		    eamInterpolateDeriv(r, rho_pot, n_values[0], &rhoTmp, &dRho);
		    rhoijprime = dRho;

#if(ALLOW_PRINTF) 
		    //printf("%d, %d, %f\n", i_particle, j_particle, r2);
		    //printf("i_particle = %d, j_particle = %d, i_b = %d, r6 = %f\n", i_particle, j_particle, i_b, r6);
#endif

		    rTmp = (fi[i_particle]+fi[j_particle])*rhoijprime/r;

		    fx_i += rTmp*dx;
		    fy_i += rTmp*dy;
		    fz_i += rTmp*dz;

		} else {
		}


	    } // loop over all atoms in jbox
	} // loop over neighbor cells

	   fx[i_particle] += fx_i;
	   fy[i_particle] += fy_i;
	   fz[i_particle] += fz_i;

    }
}

