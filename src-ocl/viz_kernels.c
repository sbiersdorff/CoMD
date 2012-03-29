
#define N_MAX_ATOMS 64
#define N_MAX_NEIGHBORS 27
#define PERIODIC 1

#define ALLOW_PRINTF 0

#if defined(cl_khr_fp64)  // Khronos extension available?
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)  // AMD extension available?
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

typedef CL_REAL_T real_t; 

__kernel void Viz(
        __global real_t* x_pos,
        __global real_t* y_pos,
        __global real_t* z_pos,
        const __global int* n_atoms,
        __global real_t* x_cen,
        __global real_t* y_cen,
        __global real_t* z_cen, 
        __global float* vertices) 
{ 
  int i_atom = get_global_id(0);
  int ibox = get_global_id(1);

  int offset = N_MAX_ATOMS; 

  if (i_atom < n_atoms[ibox]) 
  {
    vertices[i_atom*3 + offset*ibox*3 + 0] = x_cen[ibox] + x_pos[i_atom + offset*ibox];
    vertices[i_atom*3 + offset*ibox*3 + 1] = y_cen[ibox] + y_pos[i_atom + offset*ibox];
    vertices[i_atom*3 + offset*ibox*3 + 2] = z_cen[ibox] + z_pos[i_atom + offset*ibox];
  }
  else
  {
    vertices[i_atom*3 + offset*ibox*3 + 0] = 0.0f;
    vertices[i_atom*3 + offset*ibox*3 + 1] = 0.0f;
    vertices[i_atom*3 + offset*ibox*3 + 2] = 0.0f;
  }
}
