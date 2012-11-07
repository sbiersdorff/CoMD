//Initial implementation of the MD code

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

// Simple version without local blocking to check for correctness
__kernel void LJ_Force(
        __global real_t* x_pos,
        __global real_t* y_pos,
        __global real_t* z_pos,
        __global real_t* fx,
        __global real_t* fy,
        __global real_t* fz,
        __global real_t* energy,
        __global real_t* dcx,
        __global real_t* dcy,
        __global real_t* dcz,
        __global real_t* bounds,
        __global int* neighbor_list,
        __global int* n_neighbors,
        __global int* natoms,
        const real_t sigma,
        const real_t epsilon,
        const real_t cutoff) 
{

    int iatom = get_global_id(0);
    int ibox = get_global_id(1);
    int maxatoms = get_global_size(0);

    real_t dx, dy, dz;
    real_t r, r2, r6;
    real_t fr, e;

    real_t dxbox, dybox, dzbox;

    // accumulate local force value
    real_t fx_i, fy_i, fz_i;

    real_t rcut = cutoff;
    real_t r2cut = rcut*rcut;
    real_t s6 = sigma*sigma*sigma*sigma*sigma*sigma;

    int j, j_local;
    int jbox, jatom;

    int i_offset, j_offset;
    int i_particle, j_particle;

    // zero out forces on particles
    fx_i = 0.0;
    fy_i = 0.0;
    fz_i = 0.0;

    e = 0.0;

    i_offset = ibox*maxatoms; //N_MAX_ATOMS;
    i_particle = i_offset + iatom;

    if (iatom < natoms[ibox]) { // each thread executes on a single atom in the box

#if(KERN_DIAG) 
        //if (ibox < 2) printf("i = %d, %f, %f, %f\n", i_particle, x_pos[i_particle], y_pos[i_particle], z_pos[i_particle]);

        //printf("ibox = %d, n_neighbors = %d\n", ibox, n_neighbors[ibox]);
#endif

        for (j = 0; j<n_neighbors[ibox]; j++) { // loop over neighbor cells
            int jbox = neighbor_list[ibox*N_MAX_NEIGHBORS + j];
            j_offset = jbox*maxatoms; //N_MAX_ATOMS;

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

#if(KERN_DIAG) 
                //printf("dx, dy, dz = %f, %f, %f\n", dx, dy, dz);
                //printf("i = %d, j = %d, %f, %f, %f\n", i_particle, j_particle, x_pos[j_particle], y_pos[j_particle], z_pos[j_particle]);
#endif

                r2 = dx*dx + dy*dy + dz*dz;

                if ( r2 <= r2cut && r2 > 0.0) { // no divide by zero

#if(KERN_DIAG) 
                    printf("%d, %d, %f\n", i_particle, j_particle, r2);
                    //printf("r2, rcut = %f, %f\n", r2, rcut);
#endif

                    // reciprocal of r2 now
                    r2 = (real_t)1.0/r2;

                    r6 = r2*r2*r2;

                    e += r6*(s6*r6 - 1.0);

#if(KERN_DIAG) 
                    //printf("%d, %d, %f\n", i_particle, j_particle, r2);
                    //printf("i_particle = %d, j_particle = %d, i_b = %d, r6 = %f\n", i_particle, j_particle, i_b, r6);
#endif

                    fr = 4.0*epsilon*s6*r2*r6*(12.0*r6*s6 - 6.0);

                    fx_i += dx*fr;
                    fy_i += dy*fr;
                    fz_i += dz*fr;

                } else {
                }


            } // loop over all atoms
        } // loop over neighbor cells

        fx[i_particle] = fx_i;
        fy[i_particle] = fy_i;
        fz[i_particle] = fz_i;

        // since we loop over all particles, each particle contributes 1/2 the pair energy to the total
        energy[i_particle] = e*2.0*epsilon*s6;

    }
}


__kernel void LJ_Force_AoS(
        __global real3_t* pos,
        __global real3_t* f,
        __global real_t* energy,
        __global real3_t* dc,
        __global real3_t* bounds,
        __global int* neighbor_list,
        __global int* n_neighbors,
        __global int* natoms,
        const real_t sigma,
        const real_t epsilon,
        const real_t cutoff) 
{

    int iatom = get_global_id(0);
    int ibox = get_global_id(1);
    int maxatoms = get_global_size(0);

    real3_t dr;
    real_t r, r2, r6;
    real_t fr, e;

    real3_t drbox;

    // accumulate local force value
    real3_t f_i;

    real_t rcut = cutoff;
    real_t r2cut = rcut*rcut;
    real_t s6 = sigma*sigma*sigma*sigma*sigma*sigma;

    int j, j_local;
    int jbox, jatom;

    int i_offset, j_offset;
    int i_particle, j_particle;

    // zero out forces on particles
    f_i.x = 0.0;
    f_i.y = 0.0;
    f_i.z = 0.0;

    e = 0.0;

    i_offset = ibox*maxatoms; //N_MAX_ATOMS;
    i_particle = i_offset + iatom;

    if (iatom < natoms[ibox]) { // each thread executes on a single atom in the box

#if(KERN_DIAG) 
        //if (ibox < 2) printf("i = %d, %e, %e, %e\n", i_particle, pos[i_particle].x, pos[i_particle].y, pos[i_particle].z);

        //printf("ibox = %d, n_neighbors = %d\n", ibox, n_neighbors[ibox]);
        //printf("ibox = %d, natoms = %d\n", ibox, natoms[ibox]);
#endif

        for (j = 0; j<n_neighbors[ibox]; j++) { // loop over neighbor cells
            int jbox = neighbor_list[ibox*N_MAX_NEIGHBORS + j];
            j_offset = jbox*maxatoms; //N_MAX_ATOMS;

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

                j_particle = j_offset + jatom; // global offset of particle

                dr = pos[i_particle] - pos[j_particle] + drbox;

#if(KERN_DIAG) 
                //printf("dx, dy, dz = %f, %f, %f\n", dr.x, dr.y, dr.z);
                //printf("i = %d, j = %d, %f, %f, %f\n", i_particle, j_particle, pos.x[j_particle], pos.y[j_particle], pos.z[j_particle]);
#endif

                r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;

                if ( r2 <= r2cut && r2 > 0.0) { // no divide by zero

#if(KERN_DIAG) 
                    //printf("%d, %d, %f\n", i_particle, j_particle, r2);
                    //printf("r2, rcut = %f, %f\n", r2, rcut);
#endif

                    // reciprocal of r2 now
                    r2 = (real_t)1.0/r2;

                    r6 = r2*r2*r2;

                    e += r6*(s6*r6 - 1.0);

#if(KERN_DIAG) 
                    //printf("%d, %d, %f\n", i_particle, j_particle, r2);
                    //printf("i_particle = %d, j_particle = %d, i_b = %d, r6 = %f\n", i_particle, j_particle, i_b, r6);
#endif

                    fr = 4.0*epsilon*s6*r2*r6*(12.0*r6*s6 - 6.0);

                    f_i.x += dr.x*fr;
                    f_i.y += dr.y*fr;
                    f_i.z += dr.z*fr;

                } else {
                }


            } // loop over all atoms
        } // loop over neighbor cells

        f[i_particle] = f_i;

        // since we loop over all particles, each particle contributes 1/2 the pair energy to the total
        energy[i_particle] = e*2.0*epsilon*s6;

    }
}

/*
   __kernel void LJ_Force(
   __global real_t* x_pos,
   __global real_t* y_pos,
   __global real_t* z_pos,
   __global real_t* fx,
   __global real_t* fy,
   __global real_t* fz,
   __global real_t* energy,
   __global int* neighbor_list,
   __global int* n_neighbors,
   __local real_t* x_ii,
   __local real_t* y_ii,
   __local real_t* z_ii,
   __local real_t* x_ij,
   __local real_t* y_ij,
   __local real_t* z_ij,
   const real_t sigma,
   const real_t epsilon,
   const int n_cells) 
   {

   int i_atom = get_global_id(0);
   int ibox = get_global_id(1);
   int i_local = get_local_id(0);
   int n_groups = get_num_groups(0);
   int n_items = get_local_size(0);



#if(KERN_DIAG) 
printf("Number of work groups: %d\n", n_groups);
printf("Number of work items: %d\n", n_items);
#endif

real_t dx, dy, dz;
real_t r, r2, r6;
real_t fr, e;
real_t fx_ii, fy_ii, fz_ii;

real_t rcut = 5.0*sigma;
real_t r2cut = rcut*rcut;
real_t s6 = sigma*sigma*sigma*sigma*sigma*sigma;

int j, j_local;
int i_b, i_p;
int j_b;

int cell_offset, j_offset;
int group_offset;
int i_particle, j_particle;

// zero out forces on particles
fx_ii = 0.0;
fy_ii = 0.0;
fz_ii = 0.0;

e = 0.0;

i_b = get_group_id(0);

cell_offset = ibox*N_MAX_ATOMS;
group_offset = n_items*i_b;
i_particle = group_offset + i_local;

#if(KERN_DIAG) 
//printf("i_particle = %d\n", i_particle);
//printf("i_global = %d, i_local = %d, i_b = %d, n_items = %d\n", i_global, i_local, i_b, n_items);
#endif

// load particle data into local arrays
x_ii[i_local] = x_pos[i_particle + cell_offset];
y_ii[i_local] = y_pos[i_particle + cell_offset];
z_ii[i_local] = z_pos[i_particle + cell_offset];

barrier(CLK_LOCAL_MEM_FENCE);

#if(KERN_DIAG) 
//printf("x_ii, y_ii, z_ii = %f, %f, %f\n", x_ii[i_local], y_ii[i_local], z_ii[i_local]);
printf("%d, %f, %f, %f\n", i_particle, x_ii[i_local], y_ii[i_local], z_ii[i_local]);
#endif

for (j = 0; j<n_neighbors[ibox]; j++) { // loop over neighbor cells
    j_offset = neighbor_list[ibox*N_MAX_NEIGHBORS + j]*N_MAX_ATOMS;
    for (j_b = 0; j_b<n_groups; j_b++) { // loop over all groups in neighbor cell 

        // use i_local to load data in blocks of size n_items
        x_ij[i_local] = x_pos[i_local + j_b*n_items + j_offset];
        y_ij[i_local] = y_pos[i_local + j_b*n_items + j_offset];
        z_ij[i_local] = z_pos[i_local + j_b*n_items + j_offset];

        barrier(CLK_LOCAL_MEM_FENCE);

        for (j_local=0;j_local < n_items; j_local ++) { // loop over all atoms in group

            j_particle = j_local+ j_b*n_items; // global offset of particle

            dx = x_ii[i_local] - x_ij[j_local];
            dy = y_ii[i_local] - y_ij[j_local];
            dz = z_ii[i_local] - z_ij[j_local];

#if(KERN_DIAG) 
            //printf("dx, dy, dz = %f, %f, %f\n", dx, dy, dz);
            printf("%d, %f, %f, %f\n", j_particle, x_ij[j_local], y_ij[j_local], z_ij[j_local]);
            printf("%d, %d, %f, %f, %f\n", i_particle, j_particle, dx, dy, dz);
#endif

            r2 = dx*dx + dy*dy + dz*dz;

#if(KERN_DIAG) 
            printf("%d, %d, %f\n", i_particle, j_particle, r2);
            //printf("r2, rcut = %f, %f\n", r2, rcut);
#endif

            if ( r2 <= r2cut && r2 > 0.0) { // no divide by zero

                // reciprocal of r2 now
                r2 = (real_t)1.0/r2;

                r6 = r2*r2*r2;

                e += r6*(s6*r6 - 1.0);

#if(KERN_DIAG) 
                //printf("%d, %d, %f\n", i_particle, j_particle, r2);
                //printf("i_particle = %d, j_particle = %d, i_b = %d, r6 = %f\n", i_particle, j_particle, i_b, r6);
#endif

                fr = 4.0*epsilon*s6*r2*r6*(12.0*r6*s6 - 6.0);

                fx_ii += dx*fr;
                fy_ii += dy*fr;
                fz_ii += dz*fr;

            } else {
            }

        } // loop over all atoms in group

    } // loop over all groups in neighbor cell
} // loop over neighbor cells

fx[i_particle + cell_offset] = fx_ii;
fy[i_particle + cell_offset] = fy_ii;
fz[i_particle + cell_offset] = fz_ii;

barrier(CLK_LOCAL_MEM_FENCE);

// since we loop over all particles, each particle contributes 1/2 the pair energy to the total
energy[i_particle + cell_offset] = e*2.0*epsilon*s6;

}
*/
