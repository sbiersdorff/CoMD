
#include "helpers.h"

#ifdef INTEROP_VIZ

#if defined (__APPLE__) || defined(MACOSX)
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#include <OpenGL/OpenGL.h>
#include "OpenCL/cl.h"
#include "OpenCL/cl_gl.h"
#include "OpenCL/cl_gl_ext.h"
#include "OpenCL/cl_ext.h"
#include "OpenGL/CGLDevice.h"
#define GL_SHARING_EXTENSION "cl_APPLE_gl_sharing"
#else
#include <GL/gl.h>
#include <GL/glx.h>
#include "CL/cl_gl.h"
#endif

GLuint vboBuffers[3];
cl_mem vboResources[3];
int g_ncells;
Quaternion q;
float cameraFOV;
float centerX, centerY, centerZ;

#endif

// Make sure PASS_1, PASS_2, PASS_3 flags are correctly set!
// Make sure PASS_2 matches that is eam_kernels.c

#define PASS_1 1
#define PASS_2 0
#define PASS_3 1

/** This chunk of code exists only because older OpenCL implementations
  on Apple didn't treat vector types correctly. If you see errors during build time
  about subscripted values, then you have a newer version of OpenCL, and these 'if'
  clauses are no longer needed. They appear a couple more times in this file.
  **/
void dummy_test()
{
    cl_float4 dummy;

#if defined (__APPLE__) || defined(MACOSX)
    //dummy[0] = 1.0;
    //dummy[1] = 1.0;
#else
    dummy.x = 1.0;
    dummy.y = 1.0;
#endif
}

void printArray(real_t* array, int n, char *name)
{
    int i;
    printf("%s;\n", name);
    for (i=0;i<n;i++){
        printf("%d, %17.9e\n", i, array[i]);
    }
}

void computePrintEnergy(
        dev_sim_t sim_D,
        host_sim_t sim_H)
{

    /**
      Copy the array of particle energies from device to host. 
      The total energy is summed and returned in the sim_H.energy variable
     **/

    int ibox, iatom;
    double local_e;

    oclCopyToHost(sim_D.e, sim_H.e, sim_H.array_size);

    local_e = 0.0;
    for (ibox=0;ibox<sim_H.n_cells;ibox++) {
        for (iatom=0;iatom<sim_H.grid.n_atoms[ibox];iatom++) {
            local_e += sim_H.e[ibox*MAXATOMS + iatom];
        }
    }
    sim_H.energy = (cl_real)local_e;

    //printf("System energy = %30.20f\n", sim_H.energy);
    printf(" %30.20f", sim_H.energy);

}

void computePrintEnergyAoS(
        dev_sim_AoS_t sim_D,
        host_sim_AoS_t sim_H)
{

    /**
      Copy the array of particle energies from device to host. 
      The total energy is summed and returned in the sim_H.energy variable
     **/

    int ibox, iatom;
    double local_e;

    oclCopyToHost(sim_D.e, sim_H.e, sim_H.array_size);

    local_e = 0.0;
    for (ibox=0;ibox<sim_H.n_cells;ibox++) {
        for (iatom=0;iatom<sim_H.grid.n_atoms[ibox];iatom++) {
            local_e += sim_H.e[ibox*MAXATOMS + iatom];
        }
    }
    sim_H.energy = (cl_real)local_e;

    //printf("System energy = %30.20f\n", sim_H.energy);
    printf(" %30.20f", sim_H.energy);

}

void printState(
        host_sim_t sim_H,
        int n_atoms)
{

    /**
      Print the box index, atom index, position, momentum and force for the 
      first n_cells boxes of the simulation
     **/

    int i, ibox, iatom;
    int atom_count = 0;
    ibox = 0;
    printf("System state:\n");
    while (atom_count < n_atoms) {
        for (iatom=0;iatom<sim_H.grid.n_atoms[ibox];iatom++) {

            i = ibox*MAXATOMS + iatom;

            printf("%02d, %02d, "
	    "X=(%+020.12e %+020.12e %+020.12e) 1 "
	    "P=(%+020.12e %+020.12e %+020.12e) "
	    "F=(%+020.12e %+020.12e %+020.12e)\n",
                    ibox, iatom, 
                    sim_H.r.x[i],sim_H.r.y[i],sim_H.r.z[i],
                    sim_H.p.x[i],sim_H.p.y[i],sim_H.p.z[i],
                    sim_H.f.x[i],sim_H.f.y[i],sim_H.f.z[i]);
	    atom_count ++;

        }
	ibox ++;
    }

}

void printStateAoS(
        host_sim_AoS_t sim_H,
        int n_atoms)
{

    /**
      Print the box index, atom index, position, momentum and force for the 
      first n_cells boxes of the simulation
     **/

    int i, ibox, iatom;
    int atom_count = 0;
    ibox = 0;
    printf("System state:\n");
    while (atom_count < n_atoms) {
        for (iatom=0;iatom<sim_H.grid.n_atoms[ibox];iatom++) {

            i = ibox*MAXATOMS + iatom;
#if defined (__APPLE__) || defined(MACOSX)
            printf("%02d, %02d, "
	    "X=(%+020.12e %+020.12e %+020.12e) 1 "
	    "P=(%+020.12e %+020.12e %+020.12e) "
	    "F=(%+020.12e %+020.12e %+020.12e)\n",
                    ibox, iatom, 
                    sim_H.r[i][0],sim_H.r[i][1],sim_H.r[i][2],
                    sim_H.p[i][0],sim_H.p[i][1],sim_H.p[i][2],
                    sim_H.f[i][0],sim_H.f[i][1],sim_H.f[i][2]);
#else
            printf("%02d, %02d, "
	    "X=(%+020.12e %+020.12e %+020.12e) 1 "
	    "P=(%+020.12e %+020.12e %+020.12e) "
	    "F=(%+020.12e %+020.12e %+020.12e)\n",
                    ibox, iatom, 
                    sim_H.r[i].x,sim_H.r[i].y,sim_H.r[i].z,
                    sim_H.p[i].x,sim_H.p[i].y,sim_H.p[i].z,
                    sim_H.f[i].x,sim_H.f[i].y,sim_H.f[i].z);
#endif
	    atom_count ++;
        }
	ibox ++;
    }

}

void createDevGrid(dev_grid_t *grid_D, int n_r_size, int n_nl_size, int n_n_size) 
{
    /** Create the device buffers to hold the grid data:
      box centers, neighbor lists info and bounds
     **/

    oclCreateReadWriteBuffer(&grid_D->r_box.x, n_r_size);
    oclCreateReadWriteBuffer(&grid_D->r_box.y, n_r_size);
    oclCreateReadWriteBuffer(&grid_D->r_box.z, n_r_size);

    oclCreateReadWriteBuffer(&grid_D->neighbor_list, n_nl_size);
    oclCreateReadWriteBuffer(&grid_D->n_neighbors, n_n_size);
    oclCreateReadWriteBuffer(&grid_D->n_atoms, n_n_size);

    oclCreateReadWriteBuffer(&grid_D->bounds, sizeof(cl_real)*3);
}

void createDevGridAoS(dev_grid_AoS_t *grid_D, int n_r_size, int n_nl_size, int n_n_size) 
{
    /** Create the device buffers to hold the grid data:
      box centers, neighbor lists info and bounds
     **/

    oclCreateReadWriteBuffer(&grid_D->r_box, n_r_size*r3);

    oclCreateReadWriteBuffer(&grid_D->neighbor_list, n_nl_size);
    oclCreateReadWriteBuffer(&grid_D->n_neighbors, n_n_size);
    oclCreateReadWriteBuffer(&grid_D->n_atoms, n_n_size);

    oclCreateReadWriteBuffer(&grid_D->bounds, sizeof(real3_t));
}

void createDevVec(dev_vec_t *a_D, int array_size) 
{
    oclCreateReadWriteBuffer(&a_D->x, array_size);
    oclCreateReadWriteBuffer(&a_D->y, array_size);
    oclCreateReadWriteBuffer(&a_D->z, array_size);
}

void getVecAoS(cl_mem a_D,
        real3_t* a_H,
        int array_size)
{
    oclCopyToHost(a_D, a_H, r3*array_size);
}

void getVec(dev_vec_t a_D,
        host_vec_t a_H,
        int array_size)
{
    oclCopyToHost(a_D.x, a_H.x, array_size);
    oclCopyToHost(a_D.y, a_H.y, array_size);
    oclCopyToHost(a_D.z, a_H.z, array_size);
}

void getVector(cl_mem ax_D, cl_mem ay_D, cl_mem az_D,
        cl_real* ax_H, cl_real* ay_H, cl_real* az_H,
        int array_size)
{
    oclCopyToHost(ax_D, ax_H, array_size);
    oclCopyToHost(ay_D, ay_H, array_size);
    oclCopyToHost(az_D, az_H, array_size);
}

void putVector(
        cl_real* ax_H, cl_real* ay_H, cl_real* az_H,
        cl_mem ax_D, cl_mem ay_D, cl_mem az_D,
        int array_size)
{
    oclCopyToDevice(ax_H, ax_D, array_size);
    oclCopyToDevice(ay_H, ay_D, array_size);
    oclCopyToDevice(az_H, az_D, array_size);
}

void putVec(
        host_vec_t a_H,
        dev_vec_t a_D,
        int array_size)
{
    oclCopyToDevice(a_H.x, a_D.x, array_size);
    oclCopyToDevice(a_H.y, a_D.y, array_size);
    oclCopyToDevice(a_H.z, a_D.z, array_size);
}

void putGrid(host_grid_t grid_H, dev_grid_t grid_D)
{
    putVec(grid_H.r_box, grid_D.r_box, grid_H.n_r_size);

    oclCopyToDevice(grid_H.n_neighbors, grid_D.n_neighbors, grid_H.n_n_size);
    oclCopyToDevice(grid_H.n_atoms, grid_D.n_atoms, grid_H.n_n_size);
    oclCopyToDevice(grid_H.neighbor_list, grid_D.neighbor_list, grid_H.n_nl_size);

    oclCopyToDevice(grid_H.bounds, grid_D.bounds, sizeof(cl_real)*3);
}

void putGridAoS(host_grid_AoS_t grid_H, dev_grid_AoS_t grid_D)
{
    oclCopyToDevice(grid_H.r_box, grid_D.r_box, grid_H.n_r_size*r3);

    oclCopyToDevice(grid_H.n_neighbors, grid_D.n_neighbors, grid_H.n_n_size);
    oclCopyToDevice(grid_H.n_atoms, grid_D.n_atoms, grid_H.n_n_size);
    oclCopyToDevice(grid_H.neighbor_list, grid_D.neighbor_list, grid_H.n_nl_size);

    oclCopyToDevice(grid_H.bounds, grid_D.bounds, sizeof(real3_t));
}

void putEamPot(host_eam_pot_t eam_pot_H, dev_eam_pot_t eam_pot_D)
{
    oclCopyToDevice(eam_pot_H.rho, eam_pot_D.rho, eam_pot_H.n_p_rho_size);
    oclCopyToDevice(eam_pot_H.phi, eam_pot_D.phi, eam_pot_H.n_p_phi_size);
    oclCopyToDevice(eam_pot_H.F, eam_pot_D.F, eam_pot_H.n_p_F_size);

    oclCopyToDevice(eam_pot_H.n_values, eam_pot_D.n_values, sizeof(cl_int)*3);
}

void putEamCh(host_eam_ch_t eam_ch_H, dev_eam_ch_t eam_ch_D)
{
    oclCopyToDevice(eam_ch_H.rho, eam_ch_D.rho, eam_ch_H.n_ch_rho_size);
    oclCopyToDevice(eam_ch_H.phi, eam_ch_D.phi, eam_ch_H.n_ch_phi_size);
    oclCopyToDevice(eam_ch_H.F, eam_ch_D.F, eam_ch_H.n_ch_F_size);

    oclCopyToDevice(eam_ch_H.n_values, eam_ch_D.n_values, sizeof(cl_int)*3);
}


void oclRunKernel(cl_kernel kernel, cl_event *event, size_t* n_global, size_t* n_local)
{
    int err = clEnqueueNDRangeKernel(commandq, kernel, 2, NULL, n_global, n_local, 0, NULL, event);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to enqueue kernel! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    }
    clWaitForEvents(1, event);
}

void printSim(simflat_t *s,FILE *fp) 
{
    /** Print the base simulation data
      Note this is in a slightly different order than the OpenCL code returns
     **/

    int i;
    for(i=0; i<s->nboxes; i++) {
        int j;
        int ioff;
        int *id;

        for(ioff=i*MAXATOMS,j=0; j<s->natoms[i]; j++,ioff++) {
            if ( s->id[ioff] < 10) {
                fprintf(fp,
                        "%02d %02d "
			"X=(%+020.12e %+020.12e %+020.12e) 1 "
			"P=(%+020.12e %+020.12e %+020.12e) "
			"F=(%+020.12e %+020.12e %+020.12e)\n",
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

void setEAMArgs(
        cl_kernel *force_kernels,
        dev_sim_t sim_D)
{ 
    /** Set the kernel arguments for the three EAM force computation kernels
     **/

    printf("Setting EAM kernel arguments\n");
    printf("Kernel 1...");
    fflush(stdout);
    // set kernel arguments for EAM_Force_1
    int err = 0;
    int n_arg = 0;
    err  = clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.r.x);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.r.y);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.r.z);

    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.f.x);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.f.y);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.f.z);

    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.e);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.rho);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.fi);

    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.grid.r_box.x);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.grid.r_box.y);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.grid.r_box.z);

    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.grid.bounds);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.grid.neighbor_list);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.grid.n_neighbors);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.grid.n_atoms);

    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.eam_pot.n_values);

    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.eam_pot.phi);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.eam_pot.rho);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.eam_pot.F);

    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_real), &sim_D.eam_pot.cutoff);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set EAM_Force_1 arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
    fflush(stdout);
        exit(1);
    } else {
        printf("done\n");
    fflush(stdout);
    }

    // set kernel arguments for EAM_Force_2
    printf("Kernel 2...");
    fflush(stdout);
    err = 0;
    n_arg = 0;
    err  = clSetKernelArg(force_kernels[1], n_arg++, sizeof(cl_mem), &sim_D.fi);
    err |= clSetKernelArg(force_kernels[1], n_arg++, sizeof(cl_mem), &sim_D.e);
    err |= clSetKernelArg(force_kernels[1], n_arg++, sizeof(cl_mem), &sim_D.rho);

    err |= clSetKernelArg(force_kernels[1], n_arg++, sizeof(cl_mem), &sim_D.grid.n_atoms);

    err |= clSetKernelArg(force_kernels[1], n_arg++, sizeof(cl_mem), &sim_D.eam_pot.F);
    err |= clSetKernelArg(force_kernels[1], n_arg++, sizeof(cl_mem), &sim_D.eam_pot.n_values);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set EAM_Force_2 arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
    fflush(stdout);
        exit(1);
    } else {
        printf("done\n");
    fflush(stdout);
    }

    // set kernel arguments for EAM_Force_3
    printf("Kernel 3...");
    err = 0;
    n_arg = 0;
    err  = clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.r.x);
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.r.y);
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.r.z);

    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.f.x);
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.f.y);
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.f.z);

    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.fi);

    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.grid.r_box.x);
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.grid.r_box.y);
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.grid.r_box.z);

    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.grid.bounds);
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.grid.neighbor_list);
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.grid.n_neighbors);
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.grid.n_atoms);

    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.eam_pot.n_values);
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.eam_pot.rho);
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_real), &sim_D.eam_pot.cutoff);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set EAM_Force_3 arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    } else {
        printf("done\n");
    }

}
void setEAMArgsAoS(
        cl_kernel *force_kernels,
        dev_sim_AoS_t sim_D)
{ 
    /** Set the kernel arguments for the three EAM force computation kernels
     **/

    int err, n_arg;

    printf("Setting EAM kernel arguments\n");
    printf("Kernel 1...");
    fflush(stdout);
    // set kernel arguments for EAM_Force_1
    err = 0;
    n_arg = 0;
    err  = clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.r);

    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.f);

    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.e);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.rho);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.fi);

    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.grid.r_box);

    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.grid.bounds);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.grid.neighbor_list);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.grid.n_neighbors);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.grid.n_atoms);

    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.eam_pot.n_values);

    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.eam_pot.phi);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.eam_pot.rho);
    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_mem), &sim_D.eam_pot.F);

    err |= clSetKernelArg(force_kernels[0], n_arg++, sizeof(cl_real), &sim_D.eam_pot.cutoff);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set EAM_Force_1 arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    fflush(stdout);
    } else {
        printf("done\n");
    fflush(stdout);
    }

    // set kernel arguments for EAM_Force_2
    printf("Kernel 2...");
    err = 0;
    n_arg = 0;
    err  = clSetKernelArg(force_kernels[1], n_arg++, sizeof(cl_mem), &sim_D.fi);
    err |= clSetKernelArg(force_kernels[1], n_arg++, sizeof(cl_mem), &sim_D.e);
    err |= clSetKernelArg(force_kernels[1], n_arg++, sizeof(cl_mem), &sim_D.rho);

    err |= clSetKernelArg(force_kernels[1], n_arg++, sizeof(cl_mem), &sim_D.grid.n_atoms);

    err |= clSetKernelArg(force_kernels[1], n_arg++, sizeof(cl_mem), &sim_D.eam_pot.F);
    err |= clSetKernelArg(force_kernels[1], n_arg++, sizeof(cl_mem), &sim_D.eam_pot.n_values);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set EAM_Force_2 arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    fflush(stdout);
    } else {
        printf("done\n");
    fflush(stdout);
    }

    // set kernel arguments for EAM_Force_3
    printf("Kernel 3...");
    err = 0;
    n_arg = 0;
    // field arrays
    err  = clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.r);
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.f);
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.fi);
    // grid data
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.grid.r_box);
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.grid.bounds);
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.grid.neighbor_list);
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.grid.n_neighbors);
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.grid.n_atoms);
    //potential data
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.eam_pot.n_values);
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_mem), &sim_D.eam_pot.rho);
    err |= clSetKernelArg(force_kernels[2], n_arg++, sizeof(cl_real), &sim_D.eam_pot.cutoff);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set EAM_Force_3 arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
    fflush(stdout);
        exit(1);
    } else {
        printf("done\n");
	fflush(stdout);
    }

}

void setLJArgs(cl_kernel LJ_Force,
        dev_sim_t sim_D)
{
    /** Set the kernel arguments for the LJ force kernel **/

    printf("Setting LJ kernel arguments\n");
    // set kernel arguments for LJ_Force
    int err = 0;
    int n_arg = 0;
    // field arrays
    err  = clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.r.x);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.r.y);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.r.z);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.f.x);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.f.y);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.f.z);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.e);
    // grid data
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.grid.r_box.x);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.grid.r_box.y);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.grid.r_box.z);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.grid.bounds);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.grid.neighbor_list);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.grid.n_neighbors);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.grid.n_atoms);
    // potential data
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_real), &sim_D.lj_pot.sigma);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_real), &sim_D.lj_pot.epsilon);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_real), &sim_D.lj_pot.cutoff);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set LJ_Force arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    } else {
        printf("LJ_Force arguments set\n");
    }
}

void setLJArgsAoS(cl_kernel LJ_Force,
        dev_sim_AoS_t sim_D)
{
    /** Set the kernel arguments for the LJ force kernel **/

    printf("Setting LJ kernel arguments\n");
    // set kernel arguments for LJ_Force
    int err = 0;
    int n_arg = 0;
    // field arrays
    err  = clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.r);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.f);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.e);
    // grid data
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.grid.r_box);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.grid.bounds);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.grid.neighbor_list);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.grid.n_neighbors);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_mem), &sim_D.grid.n_atoms);
    // potential data
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_real), &sim_D.lj_pot.sigma);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_real), &sim_D.lj_pot.epsilon);
    err |= clSetKernelArg(LJ_Force, n_arg++, sizeof(cl_real), &sim_D.lj_pot.cutoff);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set LJ_ForceAoS arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    } else {
        printf("LJ_ForceAoS arguments set\n");
    }
}

void setAVArgs(cl_kernel AdvanceVelocity,
        dev_sim_t sim_D,
        cl_real dt)
{
    /** Set the arguments for the AdvanceVelocity kernel.
      Because of the Verlet timestepping scheme we keep the timestep as a separate argument
     **/

    int err = 0;
    int n_arg = 0;
    // field arrays
    err  = clSetKernelArg(AdvanceVelocity, n_arg++, sizeof(cl_mem), &sim_D.p.x);
    err |= clSetKernelArg(AdvanceVelocity, n_arg++, sizeof(cl_mem), &sim_D.p.y);
    err |= clSetKernelArg(AdvanceVelocity, n_arg++, sizeof(cl_mem), &sim_D.p.z);
    err |= clSetKernelArg(AdvanceVelocity, n_arg++, sizeof(cl_mem), &sim_D.f.x);
    err |= clSetKernelArg(AdvanceVelocity, n_arg++, sizeof(cl_mem), &sim_D.f.y);
    err |= clSetKernelArg(AdvanceVelocity, n_arg++, sizeof(cl_mem), &sim_D.f.z);
    // grid data
    err |= clSetKernelArg(AdvanceVelocity, n_arg++, sizeof(cl_mem), &sim_D.grid.n_atoms);
    // timestep
    err |= clSetKernelArg(AdvanceVelocity, n_arg++, sizeof(cl_real), &dt);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set AdvanceVelocity arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    } else {
        printf("AdvanceVelocity arguments set\n");
        printf("dt = %e\n", dt);
    }
}

void setAPArgs(cl_kernel AdvancePosition,
        dev_sim_t sim_D,
        cl_real dt)
{
    /** Set the arguments for the AdvancePosition kernel.
      Because of the Verlet timestepping scheme we keep the timestep as a separate argument
     **/

    int err = 0;
    int n_arg = 0;
    // field arrays
    err  = clSetKernelArg(AdvancePosition, n_arg++, sizeof(cl_mem), &sim_D.p.x);
    err |= clSetKernelArg(AdvancePosition, n_arg++, sizeof(cl_mem), &sim_D.p.y);
    err |= clSetKernelArg(AdvancePosition, n_arg++, sizeof(cl_mem), &sim_D.p.z);
    err |= clSetKernelArg(AdvancePosition, n_arg++, sizeof(cl_mem), &sim_D.r.x);
    err |= clSetKernelArg(AdvancePosition, n_arg++, sizeof(cl_mem), &sim_D.r.y);
    err |= clSetKernelArg(AdvancePosition, n_arg++, sizeof(cl_mem), &sim_D.r.z);
    //err |= clSetKernelArg(AdvancePosition, n_arg++, sizeof(cl_mem), &sim_D.m);
    // grid data
    err |= clSetKernelArg(AdvancePosition, n_arg++, sizeof(cl_mem), &sim_D.grid.n_atoms);
    // timestep
    err |= clSetKernelArg(AdvancePosition, n_arg++, sizeof(cl_real), &sim_D.rmass);
    err |= clSetKernelArg(AdvancePosition, n_arg++, sizeof(cl_real), &dt);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set AdvancePosition arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    } else {
        printf("AdvancePosition arguments set\n");
        printf("dt = %e\n", dt);
    }
}

void setAVArgsAoS(cl_kernel AdvanceVelocity,
        dev_sim_AoS_t sim_D,
        cl_real dt)
{
    /** Set the arguments for the AdvanceVelocity kernel.
      Because of the Verlet timestepping scheme we keep the timestep as a separate argument
     **/

    int err = 0;
    int n_arg = 0;
    // field arrays
    err  = clSetKernelArg(AdvanceVelocity, n_arg++, sizeof(cl_mem), &sim_D.p);
    err |= clSetKernelArg(AdvanceVelocity, n_arg++, sizeof(cl_mem), &sim_D.f);
    err |= clSetKernelArg(AdvanceVelocity, n_arg++, sizeof(cl_mem), &sim_D.grid.n_atoms);
    err |= clSetKernelArg(AdvanceVelocity, n_arg++, sizeof(cl_real), &dt);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set AdvanceVelocityAoS arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    } else {
        printf("AdvanceVelocityAoS arguments set\n");
        printf("dt = %e\n", dt);
    }
}

void setAPArgsAoS(cl_kernel AdvancePosition,
        dev_sim_AoS_t sim_D,
        cl_real dt)
{
    /** Set the arguments for the AdvancePosition kernel.
      Because of the Verlet timestepping scheme we keep the timestep as a separate argument
     **/

    int err = 0;
    int n_arg = 0;
    // field arrays
    err  = clSetKernelArg(AdvancePosition, n_arg++, sizeof(cl_mem), &sim_D.p);
    err |= clSetKernelArg(AdvancePosition, n_arg++, sizeof(cl_mem), &sim_D.r);
    //err |= clSetKernelArg(AdvancePosition, n_arg++, sizeof(cl_mem), &sim_D.m);
    err |= clSetKernelArg(AdvancePosition, n_arg++, sizeof(cl_mem), &sim_D.grid.n_atoms);
    err |= clSetKernelArg(AdvancePosition, n_arg++, sizeof(cl_real), &sim_D.rmass);
    err |= clSetKernelArg(AdvancePosition, n_arg++, sizeof(cl_real), &dt);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to set AdvancePositionAoS arguments! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    } else {
        printf("AdvancePositionAoS arguments set\n");
        printf("dt = %e\n", dt);
    }
}

void getElapsedTime(cl_event event, cl_real* elapsed_time, cl_real* enqueued_time)
{
    /** Helper routine to return the start-to-finish time (elapsed_time) 
      and the time from enqueueing to finish (enqueued_time)
     **/

    cl_ulong t_start, t_end, t_enqueue;
    int err;
    size_t param_size;


    err = clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &t_end, &param_size);
    if (err != CL_SUCCESS)
    {
        printf("Error: %s\n", print_cl_errstring(err));
        printf("t_end = %llu\n", t_end);
        //exit(1);
    }

    err = clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &t_start, &param_size);
    if (err != CL_SUCCESS)
    {
        printf("Error: %s\n", print_cl_errstring(err));
        printf("t_start = %llu\n", t_start);
        //exit(1);
    }

    err = clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &t_enqueue, &param_size);

    *elapsed_time = (t_end - t_start)*1.0e-9 ;
    *enqueued_time = (t_end - t_enqueue)*1.0e-9 ;
}

void MakeHostArrays(
        host_vec_t *r_H,
        host_vec_t *v_H,
        host_vec_t *f_H,
        host_grid_t *grid_H,
        int array_size,
        int n_r_size,
        int n_nl_size,
        int n_n_size
        )
{
    /** Allocate all the host arrays needed for the base simulation data **/

    // host memory
    // location
    r_H->x = malloc(array_size);
    r_H->y = malloc(array_size);
    r_H->z = malloc(array_size);

    // momenta
    v_H->x = malloc(array_size);
    v_H->y = malloc(array_size);
    v_H->z = malloc(array_size);

    // forces
    f_H->x = malloc(array_size);
    f_H->y = malloc(array_size);
    f_H->z = malloc(array_size);

    // box locations
    grid_H->r_box.x = malloc(n_r_size);
    grid_H->r_box.y = malloc(n_r_size);
    grid_H->r_box.z = malloc(n_r_size);

    grid_H->neighbor_list = malloc(n_nl_size);
    grid_H->n_neighbors = malloc(n_n_size);
    grid_H->n_atoms = malloc(n_n_size);
    grid_H->bounds = malloc(3*sizeof(cl_real));

    grid_H->n_n_size = n_n_size;
    grid_H->n_nl_size = n_nl_size;
    grid_H->n_r_size = n_r_size;

}

void getPrintState(dev_sim_t sim_D, host_sim_t sim_H) 
{
    getVec(sim_D.r, sim_H.r, sim_H.array_size);
    getVec(sim_D.p, sim_H.p, sim_H.array_size);
    getVec(sim_D.f, sim_H.f, sim_H.array_size);
    printState(sim_H, 2);
}

void getPrintStateAoS(dev_sim_AoS_t sim_D, host_sim_AoS_t sim_H) 
{
    getVecAoS(sim_D.r, sim_H.r, sim_H.array_size);
    getVecAoS(sim_D.p, sim_H.p, sim_H.array_size);
    getVecAoS(sim_D.f, sim_H.f, sim_H.array_size);
    printStateAoS(sim_H, 2);
}

void initHostEAM(host_eam_pot_t *eam_pot_H, simflat_t *sim) 
{
    /** Allocate and initialize all the EAM potential data needed **/

    int i;
    int n_v_rho;
    int n_v_phi;
    int n_v_F;

    // assign eam potential values
    printf("Using eam potential\n");
    eampotential_t *new_pot;
    new_pot = (eampotential_t*) sim->pot;
    eam_pot_H->cutoff  = new_pot->cutoff;
    printf("Cutoff = %e\n", eam_pot_H->cutoff);

    n_v_rho = new_pot->rho->n;
    eam_pot_H->n_p_rho_size = (6 + new_pot->rho->n)*sizeof(cl_real);
    printf("rho potential size = %d\n", eam_pot_H->n_p_rho_size);

    n_v_phi = new_pot->phi->n;
    eam_pot_H->n_p_phi_size = (6 + new_pot->phi->n)*sizeof(cl_real);
    printf("phi potential size = %d\n", eam_pot_H->n_p_phi_size);

    n_v_F = new_pot->f->n;
    eam_pot_H->n_p_F_size = (6 + new_pot->f->n)*sizeof(cl_real);
    printf("F potential size = %d\n", eam_pot_H->n_p_F_size);

    eam_pot_H->rho = malloc(eam_pot_H->n_p_rho_size);
    eam_pot_H->phi = malloc(eam_pot_H->n_p_phi_size);
    eam_pot_H->F = malloc(eam_pot_H->n_p_F_size);
    eam_pot_H->n_values = malloc(3*sizeof(cl_int));

    // Note the EAM array has 3 extra values to account for over/under flow
    // We also add another 3 values to store x0, xn, invDx
    eam_pot_H->rho[n_v_rho+3] = new_pot->rho->x0;
    eam_pot_H->rho[n_v_rho+4] = new_pot->rho->xn;
    eam_pot_H->rho[n_v_rho+5] = new_pot->rho->invDx;

    for (i=0;i<n_v_rho+3;i++)
    {
        eam_pot_H->rho[i] = new_pot->rho->values[i-1];
    }

    eam_pot_H->phi[n_v_phi+3] = new_pot->phi->x0;
    eam_pot_H->phi[n_v_phi+4] = new_pot->phi->xn;
    eam_pot_H->phi[n_v_phi+5] = new_pot->phi->invDx;

    for (i=0;i<n_v_phi+3;i++)
    {
        eam_pot_H->phi[i] = new_pot->phi->values[i-1];
    }

    eam_pot_H->F[n_v_F+3] = new_pot->f->x0;
    eam_pot_H->F[n_v_F+4] = new_pot->f->xn;
    eam_pot_H->F[n_v_F+5] = new_pot->f->invDx;

    for (i=0;i<n_v_F+3;i++)
    {
        eam_pot_H->F[i] = new_pot->f->values[i-1];
    }

    eam_pot_H->n_values[0] = n_v_phi;
    eam_pot_H->n_values[1] = n_v_rho;
    eam_pot_H->n_values[2] = n_v_F;
}

void initHostEAMCh(host_eam_ch_t *eam_ch_H, simflat_t *sim) 
{
    /** Allocate and initialize all the EAM potential data needed **/

    int i;
    int n_v_rho;
    int n_v_phi;
    int n_v_F;

    // assign eam potential values
    printf("Using eam potential\n");
    eampotential_t *new_pot;
    new_pot = (eampotential_t*) sim->pot;
    eam_ch_H->cutoff  = new_pot->cutoff;
    printf("Cutoff = %e\n", eam_ch_H->cutoff);

    n_v_rho = new_pot->rho->n;
    eam_ch_H->n_ch_rho_size = (3 + new_pot->rho->n)*sizeof(cl_real);
    printf("rho potential size = %d\n", eam_ch_H->n_ch_rho_size);

    n_v_phi = new_pot->phi->n;
    eam_ch_H->n_ch_phi_size = (3 + new_pot->phi->n)*sizeof(cl_real);
    printf("phi potential size = %d\n", eam_ch_H->n_ch_phi_size);

    n_v_F = new_pot->f->n;
    eam_ch_H->n_ch_F_size = (3 + new_pot->f->n)*sizeof(cl_real);
    printf("F potential size = %d\n", eam_ch_H->n_ch_F_size);

    eam_ch_H->rho = malloc(eam_ch_H->n_ch_rho_size);
    eam_ch_H->phi = malloc(eam_ch_H->n_ch_phi_size);
    eam_ch_H->F = malloc(eam_ch_H->n_ch_F_size);
    eam_ch_H->n_values = malloc(3*sizeof(cl_int));

    eam_ch_H->rho[n_v_rho+0] = new_pot->rho->x0;
    eam_ch_H->rho[n_v_rho+1] = new_pot->rho->xn;
    eam_ch_H->rho[n_v_rho+2] = new_pot->rho->invDx;

    for (i=0;i<n_v_rho;i++)
    {
        eam_ch_H->rho[i] = new_pot->rho->values[i];
    }

    eam_ch_H->phi[n_v_phi+0] = new_pot->phi->x0;
    eam_ch_H->phi[n_v_phi+1] = new_pot->phi->xn;
    eam_ch_H->phi[n_v_phi+2] = new_pot->phi->invDx;

    for (i=0;i<n_v_phi;i++)
    {
        eam_ch_H->phi[i] = new_pot->phi->values[i];
    }

    eam_ch_H->F[n_v_F+0] = new_pot->f->x0;
    eam_ch_H->F[n_v_F+1] = new_pot->f->xn;
    eam_ch_H->F[n_v_F+2] = new_pot->f->invDx;

    for (i=0;i<n_v_F;i++)
    {
        eam_ch_H->F[i] = new_pot->f->values[i];
    }

    eam_ch_H->n_values[0] = n_v_phi;
    eam_ch_H->n_values[1] = n_v_rho;
    eam_ch_H->n_values[2] = n_v_F;
}

void initHostLJ(host_lj_pot_t *lj_pot_H, simflat_t *sim)
{
    /** Allocate and initialize all the LJ potential data needed **/

    ljpotential_t *new_pot;
    new_pot = (ljpotential_t*) sim->pot;
    lj_pot_H->sigma  = new_pot->sigma;
    lj_pot_H->epsilon  = new_pot->epsilon;
    lj_pot_H->cutoff  = new_pot->cutoff;

    printf("Using lj potential\n");
    printf("Sigma = %e\n", lj_pot_H->sigma);
    printf("Epsilon = %e\n", lj_pot_H->epsilon);
    printf("Cutoff = %e\n", lj_pot_H->cutoff);

}

void initHostSim (host_sim_t *sim_H, simflat_t *sim)
{
    /** Allocate and initialize all the host-side simulation data needed, 
      including the appropriate potential data 
     **/

    int ibox, iatom, ioff;

    sim_H->n_cells = sim->nboxes;
    sim_H->cfac = bohr_per_atu_to_A_per_s;
    sim_H->rmass = (cl_real)((double)1822.83);
    sim_H->dt = 1.0e-15*sim_H->cfac;
    printf("dt = %e\n", sim_H->dt);


    printf("Max atom count per box:\n");
    int n_max_in_box = 0;
    int box_with_max = 0;
    for (ibox=0;ibox<sim_H->n_cells;ibox++) {
        if (sim->natoms[ibox] > n_max_in_box) {
            n_max_in_box = sim->natoms[ibox];
            box_with_max = ibox;
        }
    }
    printf("%d, %d\n", box_with_max, n_max_in_box);

    sim_H->array_size = MAXATOMS*sim->nboxes*sizeof(cl_real);
    sim_H->n_nl_size = sim->nboxes*NUMNEIGHBORS*sizeof(cl_int);
    sim_H->n_n_size = sim->nboxes*sizeof(cl_int);
    sim_H->n_r_size = sim->nboxes*sizeof(cl_real);
    sim_H->ntot = sim->ntot;

    // location
    sim_H->r.x = malloc(sim_H->array_size);
    sim_H->r.y = malloc(sim_H->array_size);
    sim_H->r.z = malloc(sim_H->array_size);

    // momenta
    sim_H->p.x = malloc(sim_H->array_size);
    sim_H->p.y = malloc(sim_H->array_size);
    sim_H->p.z = malloc(sim_H->array_size);

    // forces
    sim_H->f.x = malloc(sim_H->array_size);
    sim_H->f.y = malloc(sim_H->array_size);
    sim_H->f.z = malloc(sim_H->array_size);

    // box locations
    sim_H->grid.r_box.x = malloc(sim_H->n_r_size);
    sim_H->grid.r_box.y = malloc(sim_H->n_r_size);
    sim_H->grid.r_box.z = malloc(sim_H->n_r_size);

    sim_H->grid.neighbor_list = malloc(sim_H->n_nl_size);
    sim_H->grid.n_neighbors = malloc(sim_H->n_n_size);
    sim_H->grid.n_atoms = malloc(sim_H->n_n_size);
    sim_H->grid.bounds = malloc(3*sizeof(cl_real));

    sim_H->grid.n_n_size = sim_H->n_n_size;
    sim_H->grid.n_nl_size = sim_H->n_nl_size;
    sim_H->grid.n_r_size = sim_H->n_r_size;

    // mass, energy
    sim_H->m = malloc(sim_H->array_size);
    sim_H->e = malloc(sim_H->array_size);

    if(sim_H->eam_flag) {
        sim_H->fi = malloc(sim_H->array_size);
        sim_H->rho = malloc(sim_H->array_size);

        initHostEAM(&sim_H->eam_pot, sim);
        initHostEAMCh(&sim_H->eam_ch, sim);
    } else {
        initHostLJ(&sim_H->lj_pot, sim);
    }

    for(ibox=0;ibox<sim->nboxes;ibox++) {

        int* nbrBoxes;
        nbrBoxes = getNeighborBoxes(sim,ibox);

        sim_H->grid.n_atoms[ibox] = sim->natoms[ibox];

        sim_H->grid.r_box.x[ibox] = sim->dcenter[ibox][0];
        sim_H->grid.r_box.y[ibox] = sim->dcenter[ibox][1];
        sim_H->grid.r_box.z[ibox] = sim->dcenter[ibox][2];

        int j;
        sim_H->grid.n_neighbors[ibox] = nbrBoxes[-1];
        for (j=0;j<sim_H->grid.n_neighbors[ibox];j++) {
            sim_H->grid.neighbor_list[NUMNEIGHBORS*ibox + j] = nbrBoxes[j];
        }

        for(iatom=0;iatom<sim->natoms[ibox];iatom++) {

            ioff = ibox*MAXATOMS + iatom;

            sim_H->r.x[ioff] = sim->r[ioff][0];
            sim_H->r.y[ioff] = sim->r[ioff][1];
            sim_H->r.z[ioff] = sim->r[ioff][2];

            sim_H->p.x[ioff] = sim->p[ioff][0];
            sim_H->p.y[ioff] = sim->p[ioff][1];
            sim_H->p.z[ioff] = sim->p[ioff][2];

            sim_H->f.x[ioff] = sim->f[ioff][0];
            sim_H->f.y[ioff] = sim->f[ioff][1];
            sim_H->f.z[ioff] = sim->f[ioff][2];

        }

        for (j=0;j<3;j++) {
            sim_H->grid.bounds[j] = sim->bounds[j];
        }
    }

#ifdef INTEROP_VIZ
    centerX = sim->bounds[0]/2.0f;
    centerY = sim->bounds[1]/2.0f;
    centerZ = sim->bounds[2]/2.0f;
    printf("Center: %f %f %f\n", centerX, centerY, centerZ);
#endif

}

void initHostSimAoS (host_sim_AoS_t *sim_H, simflat_t *sim)
{
    /** Allocate and initialize all the host-side simulation data needed, 
      including the appropriate potential data 
     **/

    int ibox, iatom, ioff;

    sim_H->n_cells = sim->nboxes;
    sim_H->cfac = bohr_per_atu_to_A_per_s;
    sim_H->rmass = (cl_real)((double)1822.83);
    sim_H->dt = 1.0e-15*sim_H->cfac;
    printf("dt = %e\n", sim_H->dt);


    printf("Max atom count per box:\n");
    int n_max_in_box = 0;
    int box_with_max = 0;
    for (ibox=0;ibox<sim_H->n_cells;ibox++) {
        if (sim->natoms[ibox] > n_max_in_box) {
            n_max_in_box = sim->natoms[ibox];
            box_with_max = ibox;
        }
    }
    printf("%d, %d\n", box_with_max, n_max_in_box);

    sim_H->array_size = MAXATOMS*sim->nboxes*sizeof(cl_real);
    sim_H->n_nl_size = sim->nboxes*NUMNEIGHBORS*sizeof(cl_int);
    sim_H->n_n_size = sim->nboxes*sizeof(cl_int);
    sim_H->n_r_size = sim->nboxes*sizeof(cl_real);

    // location
    sim_H->r = malloc(sim_H->array_size*r3);

    // momenta
    sim_H->p = malloc(sim_H->array_size*r3);

    // forces
    sim_H->f = malloc(sim_H->array_size*r3);

    // box locations
    sim_H->grid.r_box       = malloc(sim_H->n_r_size*r3);

    sim_H->grid.neighbor_list = malloc(sim_H->n_nl_size);
    sim_H->grid.n_neighbors = malloc(sim_H->n_n_size);
    sim_H->grid.n_atoms     = malloc(sim_H->n_n_size);
    sim_H->grid.bounds      = malloc(sizeof(real3_t));

    sim_H->grid.n_n_size    = sim_H->n_n_size;
    sim_H->grid.n_nl_size   = sim_H->n_nl_size;
    sim_H->grid.n_r_size    = sim_H->n_r_size;
    sim_H->ntot             = sim->ntot;

    // mass, energy
    sim_H->m = malloc(sim_H->array_size);
    sim_H->e = malloc(sim_H->array_size);

    if(sim_H->eam_flag) {
        sim_H->fi = malloc(sim_H->array_size);
        sim_H->rho = malloc(sim_H->array_size);

        initHostEAM(&sim_H->eam_pot, sim);
        initHostEAMCh(&sim_H->eam_ch, sim);
    } else {
        initHostLJ(&sim_H->lj_pot, sim);
    }

    for(ibox=0;ibox<sim->nboxes;ibox++) {

        int* nbrBoxes;
        nbrBoxes = getNeighborBoxes(sim,ibox);

        sim_H->grid.n_atoms[ibox] = sim->natoms[ibox];

#if defined (__APPLE__) || defined(MACOSX)
        sim_H->grid.r_box[ibox][0] = sim->dcenter[ibox][0];
        sim_H->grid.r_box[ibox][1] = sim->dcenter[ibox][1];
        sim_H->grid.r_box[ibox][2] = sim->dcenter[ibox][2];
#else
        sim_H->grid.r_box[ibox].x = sim->dcenter[ibox][0];
        sim_H->grid.r_box[ibox].y = sim->dcenter[ibox][1];
        sim_H->grid.r_box[ibox].z = sim->dcenter[ibox][2];
#endif

        int j;
        sim_H->grid.n_neighbors[ibox] = nbrBoxes[-1];
        for (j=0;j<sim_H->grid.n_neighbors[ibox];j++) {
            sim_H->grid.neighbor_list[NUMNEIGHBORS*ibox + j] = nbrBoxes[j];
        }

        for(iatom=0;iatom<sim->natoms[ibox];iatom++) {

            ioff = ibox*MAXATOMS + iatom;
#if defined (__APPLE__) || defined(MACOSX)
            sim_H->r[ioff][0] = sim->r[ioff][0];
            sim_H->r[ioff][1] = sim->r[ioff][1];
            sim_H->r[ioff][2] = sim->r[ioff][2];

            sim_H->p[ioff][0] = sim->p[ioff][0];
            sim_H->p[ioff][1] = sim->p[ioff][1];
            sim_H->p[ioff][2] = sim->p[ioff][2];

            sim_H->f[ioff][0] = sim->f[ioff][0];
            sim_H->f[ioff][1] = sim->f[ioff][1];
            sim_H->f[ioff][2] = sim->f[ioff][2];
#else
            sim_H->r[ioff].x = sim->r[ioff][0];
            sim_H->r[ioff].y = sim->r[ioff][1];
            sim_H->r[ioff].z = sim->r[ioff][2];

            sim_H->p[ioff].x = sim->p[ioff][0];
            sim_H->p[ioff].y = sim->p[ioff][1];
            sim_H->p[ioff].z = sim->p[ioff][2];

            sim_H->f[ioff].x = sim->f[ioff][0];
            sim_H->f[ioff].y = sim->f[ioff][1];
            sim_H->f[ioff].z = sim->f[ioff][2];
#endif
        }

        for (j=0;j<3;j++) {
            sim_H->grid.bounds[j] = sim->bounds[j];
        }
    }

}

void initDevSim(dev_sim_t *sim_D, host_sim_t *sim_H)
{
    /** Allocate all the device-side arrays needed for the simulation **/

    // allocate memory buffer on device
    printf("Allocating device memory...");

    // positions
    createDevVec(&sim_D->r, sim_H->array_size);
    createDevVec(&sim_D->p, sim_H->array_size);
    createDevVec(&sim_D->f, sim_H->array_size);

    // particle mass
    oclCreateReadWriteBuffer(&sim_D->m, sim_H->array_size);

    // particle energy
    oclCreateReadWriteBuffer(&sim_D->e, sim_H->array_size);

    createDevGrid(&sim_D->grid, sim_H->n_r_size, sim_H->n_nl_size, sim_H->n_n_size);

    if (sim_H->eam_flag){
        oclCreateReadWriteBuffer(&sim_D->fi, sim_H->array_size);
        oclCreateReadWriteBuffer(&sim_D->rho, sim_H->array_size);

	//************************************************************************
	// EAM table data
        oclCreateReadWriteBuffer(&sim_D->eam_pot.rho, sim_H->eam_pot.n_p_rho_size);
        oclCreateReadWriteBuffer(&sim_D->eam_pot.phi, sim_H->eam_pot.n_p_phi_size);
        oclCreateReadWriteBuffer(&sim_D->eam_pot.F, sim_H->eam_pot.n_p_F_size);

        oclCreateReadWriteBuffer(&sim_D->eam_pot.n_values, sizeof(cl_int)*3);

        // add this here to make passing arguments to kernels easier
        sim_D->eam_pot.cutoff = sim_H->eam_pot.cutoff;

	//************************************************************************
	// EAM Chebychev coefficient data
        oclCreateReadWriteBuffer(&sim_D->eam_ch.rho, sim_H->eam_ch.n_ch_rho_size);
        oclCreateReadWriteBuffer(&sim_D->eam_ch.phi, sim_H->eam_ch.n_ch_phi_size);
        oclCreateReadWriteBuffer(&sim_D->eam_ch.F, sim_H->eam_ch.n_ch_F_size);

        oclCreateReadWriteBuffer(&sim_D->eam_ch.n_values, sizeof(cl_int)*3);

        // add this here to make passing arguments to kernels easier
        sim_D->eam_ch.cutoff = sim_H->eam_ch.cutoff;
    } else {
        sim_D->lj_pot.cutoff = sim_H->lj_pot.cutoff;
        sim_D->lj_pot.sigma = sim_H->lj_pot.sigma;
        sim_D->lj_pot.epsilon = sim_H->lj_pot.epsilon;
    }
    printf("device memory allocated\n");

}

void initDevSimAoS(dev_sim_AoS_t *sim_D, host_sim_AoS_t *sim_H)
{
    /** Allocate all the device-side arrays needed for the simulation **/

    // allocate memory buffer on device
    printf("Allocating device memory (AoS)...");

    // positions, momenta, force
    oclCreateReadWriteBuffer(&sim_D->r, sim_H->array_size*r3);
    oclCreateReadWriteBuffer(&sim_D->p, sim_H->array_size*r3);
    oclCreateReadWriteBuffer(&sim_D->f, sim_H->array_size*r3);

    // particle mass
    oclCreateReadWriteBuffer(&sim_D->m, sim_H->array_size);

    // particle energy
    oclCreateReadWriteBuffer(&sim_D->e, sim_H->array_size);

    createDevGridAoS(&sim_D->grid, sim_H->n_r_size, sim_H->n_nl_size, sim_H->n_n_size);

    if (sim_H->eam_flag){
        oclCreateReadWriteBuffer(&sim_D->fi, sim_H->array_size);
        oclCreateReadWriteBuffer(&sim_D->rho, sim_H->array_size);

	//************************************************************************
	// EAM table data
        oclCreateReadWriteBuffer(&sim_D->eam_pot.rho, sim_H->eam_pot.n_p_rho_size);
        oclCreateReadWriteBuffer(&sim_D->eam_pot.phi, sim_H->eam_pot.n_p_phi_size);
        oclCreateReadWriteBuffer(&sim_D->eam_pot.F, sim_H->eam_pot.n_p_F_size);

        oclCreateReadWriteBuffer(&sim_D->eam_pot.n_values, sizeof(cl_int)*3);

        // add this here to make passing arguments to kernels easier
        sim_D->eam_pot.cutoff = sim_H->eam_pot.cutoff;

	//************************************************************************
	// EAM Chebychev coefficient data
        oclCreateReadWriteBuffer(&sim_D->eam_ch.rho, sim_H->eam_ch.n_ch_rho_size);
        oclCreateReadWriteBuffer(&sim_D->eam_ch.phi, sim_H->eam_ch.n_ch_phi_size);
        oclCreateReadWriteBuffer(&sim_D->eam_ch.F, sim_H->eam_ch.n_ch_F_size);

        oclCreateReadWriteBuffer(&sim_D->eam_ch.n_values, sizeof(cl_int)*3);

        // add this here to make passing arguments to kernels easier
        sim_D->eam_ch.cutoff = sim_H->eam_ch.cutoff;
    } else {
        sim_D->lj_pot.cutoff = sim_H->lj_pot.cutoff;
        sim_D->lj_pot.sigma = sim_H->lj_pot.sigma;
        sim_D->lj_pot.epsilon = sim_H->lj_pot.epsilon;
    }
    printf("device memory allocated\n");

}

void putSim(host_sim_t sim_H, dev_sim_t sim_D)
{
    /** Copy all the host-side simulation data to the corresponding device arrays **/

    printf("Copying data to device...");

    // copy the input arrays to the device
    // positions

    putVec(sim_H.r, sim_D.r, sim_H.array_size);
    putVec(sim_H.p, sim_D.p, sim_H.array_size);
    putVec(sim_H.f, sim_D.f, sim_H.array_size);

    // mass
    oclCopyToDevice(sim_H.m, sim_D.m, sim_H.array_size);

    // simulation data
    putGrid(sim_H.grid, sim_D.grid);

    if(sim_H.eam_flag) {
        putEamPot(sim_H.eam_pot, sim_D.eam_pot);
    }
    printf("data copied\n");

}

void putSimAoS(host_sim_AoS_t sim_H, dev_sim_AoS_t sim_D)
{
    /** Copy all the host-side simulation data to the corresponding device arrays **/

    printf("Copying data to device (AoS)...");
    fflush(stdout);

    // copy the input arrays to the device
    // positions

    oclCopyToDevice(sim_H.r, sim_D.r, sim_H.array_size*r3);
    oclCopyToDevice(sim_H.p, sim_D.p, sim_H.array_size*r3);
    oclCopyToDevice(sim_H.f, sim_D.f, sim_H.array_size*r3);

    printf("positions...");
    fflush(stdout);

    // mass
    oclCopyToDevice(sim_H.m, sim_D.m, sim_H.array_size);

    printf("mass...");
    fflush(stdout);

    // simulation data
    putGridAoS(sim_H.grid, sim_D.grid);

    printf("grid...");
    fflush(stdout);

    if(sim_H.eam_flag) {
        putEamPot(sim_H.eam_pot, sim_D.eam_pot);
    }
    printf("data copied\n");
    fflush(stdout);

}

void FreeSims(host_sim_t sim_H, dev_sim_t sim_D)
{
    /** clean up all the host and device memory objects **/

    free(sim_H.r.x);
    free(sim_H.r.y);
    free(sim_H.r.z);

    free(sim_H.p.x);
    free(sim_H.p.y);
    free(sim_H.p.z);

    free(sim_H.f.x);
    free(sim_H.f.y);
    free(sim_H.f.z);

    free(sim_H.m);
    free(sim_H.e);

    clReleaseMemObject(sim_D.r.x);
    clReleaseMemObject(sim_D.r.y);
    clReleaseMemObject(sim_D.r.z);

    clReleaseMemObject(sim_D.p.x);
    clReleaseMemObject(sim_D.p.y);
    clReleaseMemObject(sim_D.p.z);

    clReleaseMemObject(sim_D.f.x);
    clReleaseMemObject(sim_D.f.y);
    clReleaseMemObject(sim_D.f.z);

    clReleaseMemObject(sim_D.m);
    clReleaseMemObject(sim_D.e);

    if(sim_H.eam_flag) {
        free(sim_H.rho);
        free(sim_H.fi);

        clReleaseMemObject(sim_D.rho);
        clReleaseMemObject(sim_D.fi);

    }
}

void FreeSimsAoS(host_sim_AoS_t sim_H, dev_sim_AoS_t sim_D)
{
    /** clean up all the host and device memory objects **/

    free(sim_H.r);
    free(sim_H.p);
    free(sim_H.f);
    free(sim_H.m);
    free(sim_H.e);

    clReleaseMemObject(sim_D.r);
    clReleaseMemObject(sim_D.p);
    clReleaseMemObject(sim_D.f);
    clReleaseMemObject(sim_D.m);
    clReleaseMemObject(sim_D.e);

    if(sim_H.eam_flag) {
        free(sim_H.rho);
        free(sim_H.fi);

        clReleaseMemObject(sim_D.rho);
        clReleaseMemObject(sim_D.fi);

    }
}

void buildModules(cl_kernel *force_kernels, 
        cl_kernel *AdvancePosition, 
        cl_kernel *AdvanceVelocity,
        cl_kernel *Viz, 
        host_sim_t sim_H, 
        size_t *n_local, 
        size_t *n_global)
{
    /** Build the kernels to compute force, and advance position and velocity.
      Return the appropriate global and local sizes
     **/

    cl_program timestep_module;
    cl_program lj_module;
    cl_program eam_module;
    cl_program viz_module;

    int err;

    // build the program from the kernel source file

    BuildProgramFromFile(&timestep_module, "./src-ocl/timestep_kernels.c", context, device_id);
    // only build the modules needed for the potential chosen
    if(sim_H.eam_flag) {
        BuildProgramFromFile(&eam_module, "./src-ocl/eam_kernels.c", context, device_id);
    }else{
        BuildProgramFromFile(&lj_module, "./src-ocl/lj_kernels.c", context, device_id);
    }
#ifdef INTEROP_VIZ
    BuildProgramFromFile(&viz_module, "./src-ocl/viz_kernels.c", context, device_id);
#endif

    printf("Program built\n");

    if(sim_H.eam_flag) {
        // create the EAM_Force_x kernels from the program
        force_kernels[0] = clCreateKernel(eam_module, "EAM_Force_1", &err);
        if (!force_kernels[0] || err != CL_SUCCESS)
        {
            printf("Error: Failed to create compute kernel EAM_Force_1!\n");
            exit(1);
        } else {
            printf("Kernel EAM_Force_1 built\n");
        }
        force_kernels[1] = clCreateKernel(eam_module, "EAM_Force_2", &err);
        if (!force_kernels[1] || err != CL_SUCCESS)
        {
            printf("Error: Failed to create compute kernel EAM_Force_2!\n");
            exit(1);
        } else {
            printf("Kernel EAM_Force_2 built\n");
        }
        force_kernels[2] = clCreateKernel(eam_module, "EAM_Force_3", &err);
        if (!force_kernels[2] || err != CL_SUCCESS)
        {
            printf("Error: Failed to create compute kernel EAM_Force_3!\n");
            exit(1);
        } else {
            printf("Kernel EAM_Force_3 built\n");
        }

    } else {
        // create the LJ_Force kernel from the program
        force_kernels[0] = clCreateKernel(lj_module, "LJ_Force", &err);
        if (!force_kernels[0] || err != CL_SUCCESS)
        {
            printf("Error: Failed to create compute kernel LJ_Force!\n");
            exit(1);
        } else {
            printf("Kernel LJ_Force built\n");
        }
    }
    // create the AdvanceVelocity kernel from the program
    *AdvanceVelocity = clCreateKernel(timestep_module, "AdvanceVelocity", &err);
    if (!*AdvanceVelocity || err != CL_SUCCESS)
    {
        printf("Error: Failed to create compute kernel AdvanceVelocity!\n");
        exit(1);
    } else {
        printf("Kernel AdvanceVelocity built\n");
    }
    // create the AdvancePosition kernel from the program
    *AdvancePosition = clCreateKernel(timestep_module, "AdvancePositions", &err);
    if (!*AdvancePosition || err != CL_SUCCESS)
    {
        printf("Error: Failed to create compute kernel AdvancePosition!\n");
        exit(1);
    } else {
        printf("Kernel AdvancePosition built\n");
    }

#ifdef INTEROP_VIZ
    // create the Viz kernel from the program
    *Viz = clCreateKernel(viz_module, "Viz", &err);
    if (!*Viz || err != CL_SUCCESS)
    {
        printf("Error: Failed to create compute kernel Viz!\n");
        exit(1);
    } else {
        printf("Kernel Viz built\n");
    }
#endif

    // determine allowable local work sizes for the device we chose
    if(sim_H.eam_flag) {
    } else {
        err = clGetKernelWorkGroupInfo(force_kernels[0], device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(n_local), n_local, NULL);
        if (err != CL_SUCCESS)
        {
            printf("Error: Failed to retrieve LJ_Force work group info! %d\n", err);
            printf("Error: %s\n", print_cl_errstring(err));
            exit(1);
        }
        printf("Maximum local size for LJ_Force is (%lu, %lu)\n", n_local[0], n_local[1]);
    }

    err = clGetKernelWorkGroupInfo(*AdvanceVelocity, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(n_local), n_local, NULL);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to retrieve AdvanceVelocity work group info! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    }
    printf("Maximum local size for AdvanceVelocity is (%lu, %lu)\n", n_local[0], n_local[1]);

    err = clGetKernelWorkGroupInfo(*AdvancePosition, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(n_local), n_local, NULL);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to retrieve AdvancePosition work group info! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    }
    printf("Maximum local size for AdvancePosition is (%lu, %lu)\n", n_local[0], n_local[1]);

    n_local[1] += 1;

    // set the global size equal to the numer of atoms
    n_global[0] = MAXATOMS;
    n_global[1] = sim_H.n_cells;

    // if the local size is greater than the total size, set local size to global size
    int i;
    for (i=0;i<2;i++) {
        if (n_global[i] < n_local[i]) {
            n_local[i] = n_global[i];
        }
    }
    printf("Global and local sizes are (%lu, %lu), (%lu, %lu)\n", n_global[0], n_global[1], n_local[0], n_local[1]);

}
void buildModulesAoS(cl_kernel *force_kernels, 
        cl_kernel *AdvancePosition, 
        cl_kernel *AdvanceVelocity,
        cl_kernel *Viz, 
        host_sim_AoS_t sim_H, 
        size_t *n_local, 
        size_t *n_global)
{
    /** Build the kernels to compute force, and advance position and velocity.
      Return the appropriate global and local sizes
     **/

    cl_program timestep_module;
    cl_program lj_module;
    cl_program eam_module;
    cl_program viz_module;

    int err;

    // build the program from the kernel source file

    BuildProgramFromFile(&timestep_module, "./src-ocl/timestep_kernels.c", context, device_id);
    // only build the modules needed for the potential chosen
    if(sim_H.eam_flag) {
        BuildProgramFromFile(&eam_module, "./src-ocl/eam_kernels.c", context, device_id);
    }else{
        BuildProgramFromFile(&lj_module, "./src-ocl/lj_kernels_aos.c", context, device_id);
    }
#ifdef INTEROP_VIZ
    BuildProgramFromFile(&viz_module, "./src-ocl/viz_kernels.c", context, device_id);
#endif

    printf("Program built\n");

    if(sim_H.eam_flag) {
        // create the EAM_Force_x kernels from the program
        force_kernels[0] = clCreateKernel(eam_module, "EAM_Force_1_AoS", &err);
        if (!force_kernels[0] || err != CL_SUCCESS)
        {
            printf("Error: Failed to create compute kernel EAM_Force_1!\n");
            exit(1);
        } else {
            printf("Kernel EAM_Force_1 built\n");
        }
        force_kernels[1] = clCreateKernel(eam_module, "EAM_Force_2_AoS", &err);
        if (!force_kernels[1] || err != CL_SUCCESS)
        {
            printf("Error: Failed to create compute kernel EAM_Force_2!\n");
            exit(1);
        } else {
            printf("Kernel EAM_Force_2 built\n");
        }
        force_kernels[2] = clCreateKernel(eam_module, "EAM_Force_3_AoS", &err);
        if (!force_kernels[2] || err != CL_SUCCESS)
        {
            printf("Error: Failed to create compute kernel EAM_Force_3!\n");
            exit(1);
        } else {
            printf("Kernel EAM_Force_3 built\n");
        }

    } else {
        // create the LJ_Force kernel from the program
        force_kernels[0] = clCreateKernel(lj_module, "LJ_Force_AoS", &err);
        if (!force_kernels[0] || err != CL_SUCCESS)
        {
            printf("Error: Failed to create compute kernel LJ_Force_AoS!\n");
            exit(1);
        } else {
            printf("Kernel LJ_Force_AoS built\n");
        }
    }
    // create the AdvanceVelocity kernel from the program
    *AdvanceVelocity = clCreateKernel(timestep_module, "AdvanceVelocityAoS", &err);
    if (!*AdvanceVelocity || err != CL_SUCCESS)
    {
        printf("Error: Failed to create compute kernel AdvanceVelocityAoS!\n");
        exit(1);
    } else {
        printf("Kernel AdvanceVelocity built\n");
    }
    // create the AdvancePosition kernel from the program
    *AdvancePosition = clCreateKernel(timestep_module, "AdvancePositionsAoS", &err);
    if (!*AdvancePosition || err != CL_SUCCESS)
    {
        printf("Error: Failed to create compute kernel AdvancePositionAoS!\n");
        exit(1);
    } else {
        printf("Kernel AdvancePosition built\n");
    }

#ifdef INTEROP_VIZ
    // create the Viz kernel from the program
    *Viz = clCreateKernel(viz_module, "Viz", &err);
    if (!*Viz || err != CL_SUCCESS)
    {
        printf("Error: Failed to create compute kernel Viz!\n");
        exit(1);
    } else {
        printf("Kernel Viz built\n");
    }
#endif

    // determine allowable local work sizes for the device we chose
    if(sim_H.eam_flag) {
    } else {
        err = clGetKernelWorkGroupInfo(force_kernels[0], device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(n_local), n_local, NULL);
        if (err != CL_SUCCESS)
        {
            printf("Error: Failed to retrieve LJ_Force work group info! %d\n", err);
            printf("Error: %s\n", print_cl_errstring(err));
            exit(1);
        }
        printf("Maximum local size for LJ_Force is (%lu, %lu)\n", n_local[0], n_local[1]);
    }

    err = clGetKernelWorkGroupInfo(*AdvanceVelocity, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(n_local), n_local, NULL);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to retrieve AdvanceVelocity work group info! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    }
    printf("Maximum local size for AdvanceVelocity is (%lu, %lu)\n", n_local[0], n_local[1]);

    err = clGetKernelWorkGroupInfo(*AdvancePosition, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(n_local), n_local, NULL);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to retrieve AdvancePosition work group info! %d\n", err);
        printf("Error: %s\n", print_cl_errstring(err));
        exit(1);
    }
    printf("Maximum local size for AdvancePosition is (%lu, %lu)\n", n_local[0], n_local[1]);

    //n_local[1] += 1;

    // set the global size equal to the numer of atoms
    n_global[0] = MAXATOMS;
    n_global[1] = sim_H.n_cells;

    // if the local size is greater than the total size, set local size to global size
    int i;
    for (i=0;i<2;i++) {
        if (n_global[i] < n_local[i]) {
            n_local[i] = n_global[i];
        }
    }
    printf("Global and local sizes are (%lu, %lu), (%lu, %lu)\n", n_global[0], n_global[1], n_local[0], n_local[1]);

}

void computeForceOCL(cl_kernel* force_kernels, 
cl_event *Force_event, 
size_t* n_global, 
size_t* n_local, 
int eam_flag, 
int ntot,
cl_real* t_kern)
{
    /** Execute the appropriate force kernels **/

    int err;
    cl_real t_simple, t_overall;
    cl_real t_total = 0.0;
    if (eam_flag) {
#if (PASS_1)
        if (DIAG_LEVEL > 1) {
            printf("Running EAM kernel 1..");
            fflush(stdout);
        }
        oclRunKernel(force_kernels[0], Force_event, n_global, n_local);
        err = clWaitForEvents(1, Force_event);
        if (DIAG_LEVEL > 1) {
            printf("done\n");
            fflush(stdout);
        }
        getElapsedTime(*Force_event, &t_simple, &t_overall);
        t_total += t_simple;
#endif
#if (PASS_2)
        if (DIAG_LEVEL > 1) {
            printf("Running EAM kernel 2..");
            fflush(stdout);
        }
        oclRunKernel(force_kernels[1], Force_event, n_global, n_local);
        err = clWaitForEvents(1, Force_event);
        if (DIAG_LEVEL > 1) {
            printf("done\n");
            fflush(stdout);
        }
        getElapsedTime(*Force_event, &t_simple, &t_overall);
        t_total += t_simple;
#endif
#if (PASS_3)
        if (DIAG_LEVEL > 1) {
            printf("Running EAM kernel 3..");
            fflush(stdout);
        }
        oclRunKernel(force_kernels[2], Force_event, n_global, n_local);
        err = clWaitForEvents(1, Force_event);
        if (DIAG_LEVEL > 1) {
            printf("done\n");
            fflush(stdout);
        }
        getElapsedTime(*Force_event, &t_simple, &t_overall);
        t_total += t_simple;
#endif
	*t_kern = t_total;
        if (DIAG_LEVEL > 0)
            printf("Kernel EAM_Force executed in %.3e secs. (%e us/atom for %d atoms)\n", 
                    t_total, 1.0e6*t_total/ntot, ntot);

    } else {
        if (DIAG_LEVEL > 1) {
            printf("Running LJ kernel..");
        }
        oclRunKernel(force_kernels[0], Force_event, n_global, n_local);
        err = clWaitForEvents(1, Force_event);
        if (DIAG_LEVEL > 1) {
            printf("done\n");
        }
        getElapsedTime(*Force_event, &t_simple, &t_overall);
        if (DIAG_LEVEL > 0)
            printf("Kernel LJ_Force executed in %.3e secs. (%e us/atom for %d atoms)\n", 
                    t_simple, 1.0e6*t_simple/ntot, ntot);
	*t_kern = t_simple;
    }

}

#ifdef INTEROP_VIZ
// write atom positions to VBO
void oclGraphics(cl_kernel vizKernel, dev_sim_t sim_D, size_t* n_global, size_t* n_local)
{
    // Set the arguments of our kernel
    cl_int err;  unsigned int i;
    const n = g_ncells*MAXATOMS*3;    

    cl_mem cl_vertices = vboResources[0];   

    int n_arg = 0;
    err  = clSetKernelArg(vizKernel, n_arg++, sizeof(cl_mem), &sim_D.r.x);
    err |= clSetKernelArg(vizKernel, n_arg++, sizeof(cl_mem), &sim_D.r.y);
    err |= clSetKernelArg(vizKernel, n_arg++, sizeof(cl_mem), &sim_D.r.z);
    err |= clSetKernelArg(vizKernel, n_arg++, sizeof(cl_mem), &sim_D.grid.n_atoms);
    err |= clSetKernelArg(vizKernel, n_arg++, sizeof(cl_mem), &sim_D.grid.r_box.x);
    err |= clSetKernelArg(vizKernel, n_arg++, sizeof(cl_mem), &sim_D.grid.r_box.y);
    err |= clSetKernelArg(vizKernel, n_arg++, sizeof(cl_mem), &sim_D.grid.r_box.z);
    err |= clSetKernelArg(vizKernel, n_arg++, sizeof(cl_mem), (void *) &cl_vertices);  
    clFinish(commandq);

    // Execute the kernel 
    err = clEnqueueNDRangeKernel(commandq, vizKernel, 2, NULL, n_global, NULL, 0, NULL, 0);  
    clFinish(commandq);

    // Read back and output results
    /*float* data = (float*)malloc(sizeof(float)*n); 
      err = clEnqueueReadBuffer(commandq, cl_vertices, CL_TRUE, 0, n*sizeof(float), data, 0, NULL, 0);
      clFinish(commandq);
      printf("Ran graphics kernel\n");
      for (i=0; i<n; i++) printf("%f ", data[i]);
      printf("\n");*/  
}

// render VBOs
void oclRender()
{
    cl_int err = clEnqueueReleaseGLObjects(commandq, 1, &vboResources[0], 0, 0, 0);
    clFinish(commandq);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(cameraFOV, 1.0f, 1.0f, 2000.0f);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0.0f, 0.0f, 150.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);

    glPushMatrix();
    //float halfGrid = 75.0f/2.0f;
    glTranslatef(-centerX, -centerY, -centerZ); //-halfGrid, -halfGrid, -halfGrid/2.0f);

    float rotationMatrix[16];
    QuaternionGetRotMat(rotationMatrix, q);
    glMultMatrixf(rotationMatrix);

    //float centerX = halfGrid;  float centerY = halfGrid;  float centerZ = halfGrid/2.0f;
    GLfloat matrix[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
    float offsetX = matrix[0]*centerX + matrix[1]*centerY + matrix[2]*centerZ; 
    float offsetY = matrix[4]*centerX + matrix[5]*centerY + matrix[6]*centerZ;
    float offsetZ = matrix[8]*centerX + matrix[9]*centerY + matrix[10]*centerZ;
    offsetX = centerX - offsetX; offsetY = centerY - offsetY; offsetZ = centerZ - offsetZ;
    glTranslatef(-offsetX, -offsetY, -offsetZ);

    //glPointSize(5.0f);
    glColor4f(1.0f, 0.0f, 0.0f, 1.0f);

    glEnableClientState(GL_VERTEX_ARRAY);
    //glEnableClientState(GL_COLOR_ARRAY);
    //glEnableClientState(GL_NORMAL_ARRAY);

    glBindBuffer(GL_ARRAY_BUFFER, vboBuffers[0]);
    glVertexPointer(3, GL_FLOAT, 0, 0);

    //glBindBuffer(GL_ARRAY_BUFFER, vboBuffers[1]);
    //glColorPointer(4, GL_FLOAT, 0, 0);

    //glBindBuffer(GL_ARRAY_BUFFER, vboBuffers[2]);
    //glNormalPointer(GL_FLOAT, 4*sizeof(float), 0);

    glDrawArrays(GL_POINTS, 0, g_ncells*MAXATOMS); 

    glPopMatrix();

    clEnqueueAcquireGLObjects(commandq, 1, &(vboResources[0]), 0, 0, NULL); 
    clFinish(commandq);
}


// initialize vertex buffer objects
void oclInitInterop(int ncells)
{
    g_ncells = ncells;
    int vboSize = ncells*MAXATOMS*3*sizeof(float);
    int numBuffers = 1; 
    int i;
    glGenBuffers(numBuffers, vboBuffers);
    for (i=0; i<numBuffers; i++)
    {
        unsigned int buffer_size = vboSize*sizeof(float);
        glBindBuffer(GL_ARRAY_BUFFER, vboBuffers[i]);
        glBufferData(GL_ARRAY_BUFFER, buffer_size, 0, GL_DYNAMIC_DRAW);
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    for (i=0; i<numBuffers; i++) 
    {
        cl_int status = CL_SUCCESS;
        vboResources[i] = clCreateFromGLBuffer(context, CL_MEM_WRITE_ONLY, vboBuffers[i], &status);
    }

    clEnqueueAcquireGLObjects(commandq, 1, &(vboResources[0]), 0, 0, NULL); 
    clFinish(commandq);
}

#endif

