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
    real_t* x;
    real_t* y;
    real_t* z;
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
    real_t* bounds;
    int n_n_size;
    int n_nl_size;
    int n_r_size;
} host_grid_t;

typedef struct dev_grid_t {
    dev_vec_t r_box;
    cl_mem neighbor_list;
    cl_mem n_neighbors;
    cl_mem n_atoms;
    cl_mem bounds;
} dev_grid_t;

typedef struct host_eam_pot_t {
    real_t* rho;
    real_t* phi;
    real_t* F;
    cl_int* n_values;
    real_t cutoff;
    int n_p_rho_size;
    int n_p_phi_size;
    int n_p_F_size;
} host_eam_pot_t;

typedef struct dev_eam_pot_t {
    cl_mem rho;
    cl_mem phi;
    cl_mem F;
    cl_mem n_values;
    real_t cutoff;
} dev_eam_pot_t;

typedef struct host_lj_pot_t {
    real_t cutoff;
    real_t sigma;
    real_t epsilon;
} host_lj_pot_t;

typedef struct dev_lj_pot_t {
    real_t cutoff;
    real_t sigma;
    real_t epsilon;
} dev_lj_pot_t;

typedef struct host_sim_t {
    // grid array values
    host_vec_t r;
    host_vec_t p;
    host_vec_t f;
    real_t* e;
    real_t* m;
    real_t* fi;
    real_t* rho;
    // grid data 
    host_grid_t grid;
    // real scalars
    real_t dt;
    real_t rmass;
    real_t cfac;
    real_t energy;
    // integer values
    int array_size;
    int n_n_size;
    int n_r_size;
    int n_nl_size;
    int n_cells;
    // eam flag
    int eam_flag;
    host_eam_pot_t eam_pot;
    host_lj_pot_t lj_pot;
} host_sim_t;

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
    real_t dt;
    real_t rmass;
    real_t cfac;
    real_t energy;
    // integer values
    int array_size;
    int n_n_size;
    int n_r_size;
    int n_nl_size;
    int n_cells;
    dev_eam_pot_t eam_pot;
    dev_lj_pot_t lj_pot;
} dev_sim_t;
    

void ComputePrintEnergy(dev_sim_t sim_D, host_sim_t sim_H);

void PrintState( host_sim_t sim_H, int n_cells);

void printSim(simflat_t *s,FILE *fp);

void CreateDevGrid(dev_grid_t *grid_D, int n_r_size, int n_nl_size, int n_n_size);

void CreateDevVec(dev_vec_t *a_D, int array_size);

void GetVec(dev_vec_t a_D, host_vec_t a_H, int array_size);

void GetVector(cl_mem ax_D, cl_mem ay_D, cl_mem az_D,
	real_t* ax_H, real_t* ay_H, real_t* az_H,
	int array_size);

void PutVec( host_vec_t a_H, dev_vec_t a_D, int array_size);

void PutVector( real_t* ax_H, real_t* ay_H, real_t* az_H,
        cl_mem ax_D, cl_mem ay_D, cl_mem az_D,
        int array_size);

void PutGrid(host_grid_t grid_H, dev_grid_t grid_D);

void PutEamPot(host_eam_pot_t eam_pot_H, dev_eam_pot_t eam_pot_D);

void oclRunKernel(cl_kernel kernel, cl_event *event, size_t* n_global, size_t* n_local);

void SetLJArgs(cl_kernel LJ_Force, dev_sim_t sim_D);

void SetAVArgs(cl_kernel AdvanceVelocity, dev_sim_t sim_D, real_t dt);

void SetAPArgs(cl_kernel AdvancePosition, dev_sim_t sim_D, real_t dt);

void GetElapsedTime(cl_event event, real_t* elapsed_time, real_t* enqueued_time);

void SetEAMArgs(cl_kernel *force_kernels, dev_sim_t sim_D);

#endif
