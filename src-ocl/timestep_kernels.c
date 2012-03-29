// timestep subroutines

/**
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

/* CL_REAL_T is set to single or double depending on compile time flags */
typedef CL_REAL_T real_t; 
//#pragma OPENCL EXTENSION cl_amd_fp64 : enable

__kernel void AdvanceVelocity (
        __global real_t* px,
        __global real_t* py,
        __global real_t* pz,
        const __global real_t* fx,
        const __global real_t* fy,
        const __global real_t* fz,
        const __global int* n_atoms,
        const real_t dt
        )

{
    int i_atom = get_global_id(0);
    int ibox = get_global_id(1);

    int offset = get_global_size(0);

    real_t dt_local = dt;

#if(ALLOW_PRINTF) 
    if (n_atoms[ibox] > 0) 
        printf("%d, %d, %d\n", ibox, i_atom, n_atoms[ibox]);
#endif

    if (i_atom < n_atoms[ibox]) {
        px[i_atom + offset*ibox] -= dt_local*fx[i_atom + offset*ibox];
        py[i_atom + offset*ibox] -= dt_local*fy[i_atom + offset*ibox];
        pz[i_atom + offset*ibox] -= dt_local*fz[i_atom + offset*ibox];
    }

}

__kernel void AdvancePositions (
        const __global real_t* px,
        const __global real_t* py,
        const __global real_t* pz,
        __global real_t* x_pos,
        __global real_t* y_pos,
        __global real_t* z_pos,
        const __global real_t* mass,
        const __global int* n_atoms,
        const real_t rmass,
        const real_t dt
        )

{
    int i_atom = get_global_id(0);
    int ibox = get_global_id(1);

    int offset = N_MAX_ATOMS; //get_global_size(0);

#if(ALLOW_PRINTF) 
    if (n_atoms[ibox] > 0) 
        printf("%d, %d, %d\n", ibox, i_atom, n_atoms[ibox]);
#endif

    real_t dt_local = dt/rmass; 

#if(ALLOW_PRINTF) 
    printf("%f, %f, %f\n", dt_local, dt, rmass);
#endif

    if (i_atom < n_atoms[ibox]) {
        x_pos[i_atom + offset*ibox] += dt_local*px[i_atom + offset*ibox];
        y_pos[i_atom + offset*ibox] += dt_local*py[i_atom + offset*ibox];
        z_pos[i_atom + offset*ibox] += dt_local*pz[i_atom + offset*ibox];

#if(ALLOW_PRINTF) 
        printf("%d, %d, %f, %f, %f\n", ibox, i_atom, x_pos[i_atom + offset*ibox], y_pos[i_atom + offset*ibox], z_pos[i_atom + offset*ibox]);
#endif
    }

}


