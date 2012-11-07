// timestep subroutines

/**
  Since OpenCL doesn't pick up #include properly, we need to manually switch real_t from 
  float to double in each kernel file individually.
  **/

#define N_MAX_NEIGHBORS 27
#define PERIODIC 1

#define KERN_DIAG 0

#if defined(cl_khr_fp64)  // Khronos extension available?
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)  // AMD extension available?
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

/* CL_REAL_T is set to single or double depending on compile time flags */
typedef CL_REAL_T real_t; 
typedef CL_REAL4_T real3_t; 


__kernel void AdvanceVelocity (
        __global real_t* px,
        __global real_t* py,
        __global real_t* pz,
        __global const real_t* fx,
        __global const real_t* fy,
        __global const real_t* fz,
        __global const int* n_atoms,
        const real_t dt
        )

{
    int iatom = get_global_id(0);
    int ibox = get_global_id(1);

    int offset = get_global_size(0);
    int ioff = iatom + offset*ibox;

    real_t dt_local = dt;

#if(KERN_DIAG > 0) 
    if (iatom == 0 && ibox == 0) printf(" AV dt = %e\n", dt);

    if (n_atoms[ibox] > 0) printf("%d, %d, %d\n", ibox, iatom, n_atoms[ibox]);
#endif

    if (iatom < n_atoms[ibox]) {
        px[ioff] -= dt_local*fx[ioff];
        py[ioff] -= dt_local*fy[ioff];
        pz[ioff] -= dt_local*fz[ioff];
    }

}

__kernel void AdvancePositions (
        __global const real_t* px,
        __global const real_t* py,
        __global const real_t* pz,
        __global real_t* x_pos,
        __global real_t* y_pos,
        __global real_t* z_pos,
        //__global const real_t* mass,
        __global const int* n_atoms,
        const real_t rmass,
        const real_t dt
        )

{
    int iatom = get_global_id(0);
    int ibox = get_global_id(1);

    int offset = get_global_size(0);
    int ioff = iatom + offset*ibox;

#if(KERN_DIAG > 0) 

    if (n_atoms[ibox] > 0) printf("%d, %d, %d\n", ibox, iatom, n_atoms[ibox]);
#endif

    real_t dt_local = dt/rmass; 

#if(KERN_DIAG > 0) 
    if (iatom == 0 && ibox == 0) printf(" AP dt = %e\n", dt_local);
    if (iatom == 0 && ibox == 0) printf("%f, %f, %f\n", dt_local, dt, rmass);
#endif

    if (iatom < n_atoms[ibox]) {
        x_pos[ioff] += dt_local*px[ioff];
        y_pos[ioff] += dt_local*py[ioff];
        z_pos[ioff] += dt_local*pz[ioff];

#if(KERN_DIAG > 0) 
    if (iatom == 0 && ibox == 0) printf("%d, %d, %f, %f, %f\n", ibox, iatom, x_pos[ioff], y_pos[ioff], z_pos[ioff]);
#endif
    }

}

__kernel void AdvanceVelocityAoS (
        __global real3_t* p,
        __global const real3_t* f,
        __global const int* n_atoms,
        const real_t dt
        )

{
    int iatom = get_global_id(0);
    int ibox = get_global_id(1);

    int offset = get_global_size(0);
    int ioff = iatom + offset*ibox;

    real3_t dt_local = {dt, dt, dt, 0.0};
    //real3_t dt_local = {0.0, 0.0, 0.0, 0.0};

#if(KERN_DIAG > 0) 
    if (iatom == 0 && ibox == 0) printf("AoS AV dt = %e\n", dt);

    if (n_atoms[ibox] > 0) printf("%d, %d, %d\n", ibox, iatom, n_atoms[ibox]);
#endif

    if (iatom < n_atoms[ibox]) {
        p[ioff] -= dt_local*f[ioff];
    }

}

__kernel void AdvancePositionsAoS (
        __global const real3_t* p,
        __global real3_t* pos,
        //__global const real_t* mass,
        __global const int* n_atoms,
        const real_t rmass,
        const real_t dt
        )

{
    int iatom = get_global_id(0);
    int ibox = get_global_id(1);

    int offset = get_global_size(0);
    int ioff = iatom + offset*ibox;

#if(KERN_DIAG > 0) 
    if (iatom == 0 && ibox == 0) printf("AoS AP dt = %e\n", dt);

    if (n_atoms[ibox] > 0) printf("%d, %d, %d\n", ibox, iatom, n_atoms[ibox]);
#endif

    real_t rdt = dt/rmass;

    real3_t dt_local = {rdt, rdt, rdt, 0.0};

#if(KERN_DIAG > 0) 
    if (iatom == 0 && ibox == 0) printf("%e, %e, %e\n", dt_local.x, dt, rmass);
#endif

    if (iatom < n_atoms[ibox]) {
        pos[ioff] += dt_local*p[ioff];

#if(KERN_DIAG > 0) 
    if (iatom == 0 && ibox == 0) printf("%d, %d, %f, %f, %f\n", ibox, iatom, pos[ioff].x, pos[ioff].y, pos[ioff].z);
#endif
    }

}


