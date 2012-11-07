#ifndef HELPER_H
#define HELPER_H

#include "cl_utils.h"
#include "pmdOCL.h"
#include "pmdTypes.h"
#include "domains.h"
#include "constants.h"
#include "quaternion.h"

#define DIAG_LEVEL 0

typedef struct host_vec_t {
    cl_real* x;
    cl_real* y;
    cl_real* z;
} host_vec_t;

typedef struct dev_vec_t {
    cl_mem x;
    cl_mem y;
    cl_mem z;
} dev_vec_t;

typedef struct host_grid_t {
    host_vec_t r_box;
    cl_int* neighbor_list;
    cl_int* n_neighbors;
    cl_int* n_atoms;
    cl_real* bounds;
    int n_n_size;
    int n_nl_size;
    int n_r_size;
} host_grid_t;

typedef struct host_grid_AoS_t {
    real3_t* r_box;
    cl_int* neighbor_list;
    cl_int* n_neighbors;
    cl_int* n_atoms;
    cl_real* bounds;
    int n_n_size;
    int n_nl_size;
    int n_r_size;
} host_grid_AoS_t;

typedef struct dev_grid_t {
    dev_vec_t r_box;
    cl_mem neighbor_list;
    cl_mem n_neighbors;
    cl_mem n_atoms;
    cl_mem bounds;
} dev_grid_t;

typedef struct dev_grid_AoS_t {
    cl_mem r_box;
    cl_mem neighbor_list;
    cl_mem n_neighbors;
    cl_mem n_atoms;
    cl_mem bounds;
} dev_grid_AoS_t;

typedef struct host_eam_pot_t {
    cl_real* rho;
    cl_real* phi;
    cl_real* F;
    cl_int* n_values;
    cl_real cutoff;
    int n_p_rho_size;
    int n_p_phi_size;
    int n_p_F_size;
} host_eam_pot_t;

typedef struct dev_eam_pot_t {
    cl_mem rho;
    cl_mem phi;
    cl_mem F;
    cl_mem n_values;
    cl_real cutoff;
} dev_eam_pot_t;

typedef struct host_eam_ch_t {
    cl_real* rho;
    cl_real* phi;
    cl_real* F;
    cl_int* n_values;
    cl_real cutoff;
    int n_ch_rho_size;
    int n_ch_phi_size;
    int n_ch_F_size;
} host_eam_ch_t;

typedef struct dev_eam_ch_t {
    cl_mem rho;
    cl_mem phi;
    cl_mem F;
    cl_mem n_values;
    cl_real cutoff;
} dev_eam_ch_t;

typedef struct host_lj_pot_t {
    cl_real cutoff;
    cl_real sigma;
    cl_real epsilon;
} host_lj_pot_t;

typedef struct dev_lj_pot_t {
    cl_real cutoff;
    cl_real sigma;
    cl_real epsilon;
} dev_lj_pot_t;

typedef struct host_sim_t {
    // grid array values
    host_vec_t r;
    host_vec_t p;
    host_vec_t f;
    cl_real* e;
    cl_real* m;
    cl_real* fi;
    cl_real* rho;
    // grid data 
    host_grid_t grid;
    // real scalars
    cl_real dt;
    cl_real rmass;
    cl_real cfac;
    cl_real energy;
    // integer values
    int array_size;
    int n_n_size;
    int n_r_size;
    int n_nl_size;
    int n_cells;
    int ntot;
    // eam flag
    int eam_flag;
    host_eam_pot_t eam_pot;
    host_eam_ch_t eam_ch;
    host_lj_pot_t lj_pot;
} host_sim_t;

typedef struct host_sim_AoS_t {
    // grid array values
    real3_t* r;
    real3_t* p;
    real3_t* f;
    cl_real* e;
    cl_real* m;
    cl_real* fi;
    cl_real* rho;
    // grid data 
    host_grid_AoS_t grid;
    // real scalars
    cl_real dt;
    cl_real rmass;
    cl_real cfac;
    cl_real energy;
    // integer values
    int array_size;
    int n_n_size;
    int n_r_size;
    int n_nl_size;
    int n_cells;
    int ntot;
    // eam flag
    int eam_flag;
    host_eam_pot_t eam_pot;
    host_eam_ch_t eam_ch;
    host_lj_pot_t lj_pot;
} host_sim_AoS_t;

typedef struct dev_sim_t {
    // grid array values
    dev_vec_t r;
    dev_vec_t p;
    dev_vec_t f;
    cl_mem e;
    cl_mem m;
    cl_mem fi;
    cl_mem rho;
    // grid data 
    dev_grid_t grid;
    // note all scalars can be passed directly to the kernels
    // real scalars
    cl_real dt;
    cl_real rmass;
    cl_real cfac;
    cl_real energy;
    // integer values
    int array_size;
    int n_n_size;
    int n_r_size;
    int n_nl_size;
    int n_cells;
    dev_eam_pot_t eam_pot;
    dev_eam_ch_t eam_ch;
    dev_lj_pot_t lj_pot;
} dev_sim_t;

typedef struct dev_sim_AoS_t {
    // grid array values
    cl_mem r;
    cl_mem p;
    cl_mem f;
    cl_mem e;
    cl_mem m;
    cl_mem fi;
    cl_mem rho;
    // grid data 
    dev_grid_AoS_t grid;
    // note all scalars can be passed directly to the kernels
    // real scalars
    cl_real dt;
    cl_real rmass;
    cl_real cfac;
    cl_real energy;
    // integer values
    int array_size;
    int n_n_size;
    int n_r_size;
    int n_nl_size;
    int n_cells;
    dev_eam_pot_t eam_pot;
    dev_eam_ch_t eam_ch;
    dev_lj_pot_t lj_pot;
} dev_sim_AoS_t;
    

/* General helper utils */

void printArray(real_t* array, int n, char *name);

void printSim(simflat_t *s,FILE *fp);

void createDevVec(dev_vec_t *a_D, int array_size);

void getVector(cl_mem ax_D, cl_mem ay_D, cl_mem az_D,
	cl_real* ax_H, cl_real* ay_H, cl_real* az_H,
	int array_size);

void putVec(host_vec_t a_H, dev_vec_t a_D, int array_size);

void putVector(cl_real* ax_H, cl_real* ay_H, cl_real* az_H, cl_mem ax_D, cl_mem ay_D, cl_mem az_D, int array_size);

void putEamPot(host_eam_pot_t eam_pot_H, dev_eam_pot_t eam_pot_D);

void oclRunKernel(cl_kernel kernel, cl_event *event, size_t* n_global, size_t* n_local);

void getElapsedTime(cl_event event, cl_real* elapsed_time, cl_real* enqueued_time);

void computeForceOCL(cl_kernel* force_kernels, cl_event* Force_event, size_t* n_global, size_t* n_local, int eam_flag, int ntot, cl_real* t_kern);

/* SoA (default) variants) */
void computePrintEnergy(dev_sim_t sim_D, host_sim_t sim_H);

void printState(host_sim_t sim_H, int n_cells);

void createDevGrid(dev_grid_t *grid_D, int n_r_size, int n_nl_size, int n_n_size);

void getVec(dev_vec_t a_D, host_vec_t a_H, int array_size);

void buildModules(cl_kernel *force_kernels, cl_kernel *AdvancePosition, cl_kernel *AdvanceVelocity, cl_kernel *Viz, 
        host_sim_t sim_H, size_t *n_local, size_t *n_global);

void initHostSim (host_sim_t *sim_H, simflat_t *sim);

void initDevSim(dev_sim_t *sim_D, host_sim_t *sim_H);

void putSim(host_sim_t sim_H, dev_sim_t sim_D);

void putGrid(host_grid_t grid_H, dev_grid_t grid_D);

void setLJArgs(cl_kernel LJ_Force, dev_sim_t sim_D);

void setEAMArgs(cl_kernel *force_kernels, dev_sim_t sim_D);

void setAVArgs(cl_kernel AdvanceVelocity, dev_sim_t sim_D, cl_real dt);

void setAPArgs(cl_kernel AdvancePosition, dev_sim_t sim_D, cl_real dt);

/* AoS variants */
void computePrintEnergyAoS(dev_sim_AoS_t sim_D, host_sim_AoS_t sim_H);

void printStateAoS(host_sim_AoS_t sim_H, int n_cells);

void createDevGridAoS(dev_grid_AoS_t *grid_D, int n_r_size, int n_nl_size, int n_n_size);

void getVecAoS(cl_mem a_D, real3_t* a_H, int array_size);

void buildModulesAoS(cl_kernel *force_kernels, cl_kernel *AdvancePosition, cl_kernel *AdvanceVelocity, cl_kernel *Viz, 
        host_sim_AoS_t sim_H, size_t *n_local, size_t *n_global);

void initHostSimAoS (host_sim_AoS_t *sim_H, simflat_t *sim);

void initDevSimAoS(dev_sim_AoS_t *sim_D, host_sim_AoS_t *sim_H);

void putSimAoS(host_sim_AoS_t sim_H, dev_sim_AoS_t sim_D);

void putGridAoS(host_grid_AoS_t grid_H, dev_grid_AoS_t grid_D);

void setLJArgsAoS(cl_kernel LJ_Force, dev_sim_AoS_t sim_D);

void setEAMArgsAoS(cl_kernel *force_kernels, dev_sim_AoS_t sim_D);

void setAVArgsAoS(cl_kernel AdvanceVelocity, dev_sim_AoS_t sim_D, cl_real dt);

void setAPArgsAoS(cl_kernel AdvancePosition, dev_sim_AoS_t sim_D, cl_real dt);

/* Graphics kernels */

void oclGraphics(cl_kernel vizKernel, dev_sim_t sim_D, size_t* n_global, size_t* n_local);

void oclRender();

void oclInitInterop(int ncells);

#endif
