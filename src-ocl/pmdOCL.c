/**
 * a simple md simulator
 **/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <pthread.h>

#include "mytype.h"
#include "pmdOCL.h"
#include "helpers.h"

#ifdef INTEROP_VIZ

#if defined (__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#else
#include <glut.h>
#endif

int mouse_old_x, mouse_old_y, mouse_buttons;
int g_threadStatus, g_clFinished;
extern Quaternion q;
extern float cameraFOV;
#endif

cl_kernel AdvanceVelocity;
cl_kernel AdvancePosition;
cl_kernel Viz;
cl_kernel *Force_kernels;

cl_kernel AdvanceVelocityAoS;
cl_kernel AdvancePositionAoS;
cl_kernel VizAoS;
cl_kernel *Force_kernelsAoS;

int err;
size_t n_local[2];
size_t n_global[2];

cl_event Force_event;
cl_event AP_event;
cl_event AV_event;

real_t t_exec;
real_t t_overall;
real_t t_local;
real_t t_dummy;
cl_real dt;

host_sim_t sim_H;
dev_sim_t sim_D;

host_sim_AoS_t sim_AoS_H;
dev_sim_AoS_t sim_AoS_D;

double ts, te;

int eam_flag;
int gpu_flag;
simflat_t *sim;
int iter, ns;
int niter, nsteps;


void computeIterationSoA()
{
    real_t t_kern, t_enq;
    real_t t_acc = 0.0;
    cl_event loc_event;
    for (ns = 0;ns < nsteps; ns++) {
	// advance particle positions dt/2
	oclRunKernel(AdvancePosition, &AP_event, n_global, n_local);
	clWaitForEvents(1, &AP_event);
	getElapsedTime(AP_event, &t_kern, &t_enq);
	t_acc += t_kern;

	if (DIAG_LEVEL > 1)
	    getPrintState(sim_D, sim_H);

	// compute force
	computeForceOCL(Force_kernels, &Force_event, n_global, n_local, eam_flag, sim_H.ntot, &t_kern);
	t_acc += t_kern;

	if (DIAG_LEVEL > 1)
	    getPrintState(sim_D, sim_H);

	// advance velocity a full timestep
	oclRunKernel(AdvanceVelocity, &AV_event, n_global, n_local);
	clWaitForEvents(1, &AV_event);
	getElapsedTime(AV_event, &t_kern, &t_enq);
	t_acc += t_kern;

	if (DIAG_LEVEL > 1)
	    getPrintState(sim_D, sim_H);

	// advance particle positions dt/2
	oclRunKernel(AdvancePosition, &AP_event, n_global, n_local);
	clWaitForEvents(1, &AP_event);
	getElapsedTime(AP_event, &t_kern, &t_enq);
	t_acc += t_kern;

	if (DIAG_LEVEL > 1)
	    getPrintState(sim_D, sim_H);
    }

#ifdef INTEROP_VIZ 
    oclGraphics(Viz, sim_D, n_global, n_local);
#endif

    // compute force
    computeForceOCL(Force_kernels, &Force_event, n_global, n_local, eam_flag, sim_H.ntot, &t_kern);

    if (DIAG_LEVEL > 0)
	getPrintState(sim_D, sim_H);

    printf(" %4d ", nsteps*(iter + 1));
    computePrintEnergy(sim_D, sim_H);
    printf("    computed in  = %e (%e us/atom for %d atoms)\n", 
	    t_acc,1.0e6*t_acc/(double)(sim->ntot*nsteps), sim->ntot);
}

void computeIterationAoS()
{
    real_t t_kern, t_enq;
    real_t t_acc = 0.0;
    for (ns = 0;ns < nsteps; ns++) {
	// advance particle positions dt/2
	oclRunKernel(AdvancePositionAoS, &AP_event, n_global, n_local);
	clWaitForEvents(1, &AP_event);
	getElapsedTime(AP_event, &t_kern, &t_enq);
	t_acc += t_kern;

	if (DIAG_LEVEL > 1)
	    getPrintStateAoS(sim_AoS_D, sim_AoS_H);

	// compute force
	computeForceOCL(Force_kernelsAoS, &Force_event, n_global, n_local, eam_flag, sim_H.ntot, &t_kern);
	t_acc += t_kern;

	if (DIAG_LEVEL > 1)
	    getPrintStateAoS(sim_AoS_D, sim_AoS_H);

	// advance velocity a full timestep
	oclRunKernel(AdvanceVelocityAoS, &AV_event, n_global, n_local);
	clWaitForEvents(1, &AV_event);
	getElapsedTime(AV_event, &t_kern, &t_enq);
	t_acc += t_kern;

	if (DIAG_LEVEL > 1)
	    getPrintStateAoS(sim_AoS_D, sim_AoS_H);

	// advance particle positions dt/2
	oclRunKernel(AdvancePositionAoS, &AP_event, n_global, n_local);
	clWaitForEvents(1, &AP_event);
	getElapsedTime(AP_event, &t_kern, &t_enq);
	t_acc += t_kern;

	if (DIAG_LEVEL > 1)
	    getPrintStateAoS(sim_AoS_D, sim_AoS_H);
    }

#ifdef INTEROP_VIZ 
    //oclGraphics(Viz, sim_D, n_global, n_local);
#endif

    // compute force
    computeForceOCL(Force_kernelsAoS, &Force_event, n_global, n_local, eam_flag, sim_H.ntot, &t_kern);

    if (DIAG_LEVEL > 0)
	getPrintStateAoS(sim_AoS_D, sim_AoS_H);

    //printf(" after %4d steps: ", nsteps*(iter + 1));
    printf(" %4d ", nsteps*(iter + 1));
    computePrintEnergyAoS(sim_AoS_D, sim_AoS_H);
    printf("    computed in  = %e (%e us/atom for %d atoms)\n", 
	    t_acc,1.0e6*t_acc/(double)(sim->ntot*nsteps), sim->ntot);
}

void finishOclSoA()
{
    getElapsedTime(AP_event, &t_exec, &t_overall);
    printf("Kernel AdvancePosition executed in %e secs.\n", t_exec);
    getElapsedTime(AV_event, &t_exec, &t_overall);
    printf("Kernel AdvanceVelocity executed in %e secs.\n", t_exec);

    real_t t_ref = (real_t)(te - ts);

    printf("Reference wallclock elapsed time = %e (%e us/atom for %d atoms)\n", 
	    t_ref,1.0e6*t_ref/(double)(sim->ntot*niter*nsteps), sim->ntot);

    // copy result back from device
    // note this utility forces synchronous copy so 
    // it does not return until the copy is complete
    getVec(sim_D.f, sim_H.f, sim_H.array_size);
    oclCopyToHost(sim_D.e, sim_H.e, sim_H.array_size);

    computePrintEnergy(sim_D, sim_H);
    printf("\n");

    // clean up the OpenCL memory 
    FreeSims(sim_H, sim_D);

}


void finishOclAoS()
{
    getElapsedTime(AP_event, &t_exec, &t_overall);
    printf("Kernel AdvancePositionAoS executed in %e secs.\n", t_exec);
    getElapsedTime(AV_event, &t_exec, &t_overall);
    printf("Kernel AdvanceVelocityAoS executed in %e secs.\n", t_exec);

    real_t t_ref = (real_t)(te - ts);

    printf("Reference wallclock elapsed time = %e (%e us/atom for %d atoms)\n", 
	    t_ref,1.0e6*t_ref/(double)(sim->ntot*niter*nsteps), sim->ntot);

    // copy result back from device
    // note this utility forces synchronous copy so 
    // it does not return until the copy is complete
    getVecAoS(sim_AoS_D.f, sim_AoS_H.f, sim_AoS_H.array_size);
    oclCopyToHost(sim_AoS_D.e, sim_AoS_H.e, sim_AoS_H.array_size);

    computePrintEnergyAoS(sim_AoS_D, sim_AoS_H);
    printf("\n");

    // clean up the OpenCL memory 
    FreeSimsAoS(sim_AoS_H, sim_AoS_D);

    clReleaseEvent(AP_event);
    clReleaseEvent(AV_event);

    oclCleanup();
}

void runRef()
{
    printf("************************************************************************\n");
    printf("Running reference simulation\n");
    // write initial state 
    writeClsman(sim,(char *) "init.bin");

    // do the computation 
    (void)do_compute_work(sim); 

    // write final configuration 
    writeClsman(sim,(char *) "final.bin");

    // free memory 
    destroySimulation(&sim);
}

void computeInitSoA()
{
    printf("************************************************************************\n");
    printf("Initializinging SoA OpenCL simulation\n");

    sim_D.rmass = (cl_real)(amu_to_m_e)*(double)(sim->pot->mass);
    sim_H.eam_flag = eam_flag;

    if(eam_flag) {
	Force_kernels = malloc(sizeof(cl_kernel)*3);
    } else {
	Force_kernels = malloc(sizeof(cl_kernel)*1);
    }

    // set up arrays for OpenCL

    initHostSim(&sim_H, sim);

    initDevSim(&sim_D, &sim_H);

    putSim(sim_H, sim_D);

    // build the program from the kernel source file

    buildModules(Force_kernels, &AdvancePosition, &AdvanceVelocity, &Viz, sim_H, n_local, n_global);

    cl_real dthalf = 0.5*sim_H.dt;
    cl_real dtminushalf = -0.5*sim_H.dt;

    if (eam_flag) {
	// set the arguments for all 3 EAM_Force kernels
	setEAMArgs(Force_kernels, sim_D);

    } else {
	// set kernel arguments for LJ_Force
	setLJArgs(Force_kernels[0], sim_D);
    }

    // set kernel arguments for AdvanceVelocity
    setAVArgs(AdvanceVelocity, sim_D, sim_H.dt);

    // set kernel arguments for AdvancePosition
    setAPArgs(AdvancePosition, sim_D, dthalf);

    // Start the simulation here;
    printf("************************************************************************\n");
    printf("Starting SoA OpenCL simulation\n");

    //ts = timeNow();

    cl_real t_kern;
    computeForceOCL(Force_kernels, &Force_event, n_global, n_local, eam_flag, sim_H.ntot, &t_kern);

    getVec(sim_D.f, sim_H.f, sim_H.array_size);

    printf(" Initial:\n");
    computePrintEnergy(sim_D, sim_H);

    getPrintState(sim_D, sim_H);
}

void computeInitAoS()
{
    printf("************************************************************************\n");
    printf("Initializinging AoS OpenCL simulation\n");

    sim_AoS_D.rmass = (cl_real)(amu_to_m_e)*(double)(sim->pot->mass);
    sim_AoS_H.eam_flag = eam_flag;

    if(eam_flag) {
	Force_kernelsAoS = malloc(sizeof(cl_kernel)*3);
    } else {
	Force_kernelsAoS = malloc(sizeof(cl_kernel)*1);
    }

    // set up arrays for OpenCL

    initHostSimAoS(&sim_AoS_H, sim);

    initDevSimAoS(&sim_AoS_D, &sim_AoS_H);

    putSimAoS(sim_AoS_H, sim_AoS_D);

    // build the program from the kernel source file

    buildModulesAoS(Force_kernelsAoS, &AdvancePositionAoS, &AdvanceVelocityAoS, &VizAoS, sim_AoS_H, n_local, n_global);

    cl_real dthalf = 0.5*sim_AoS_H.dt;
    cl_real dtminushalf = -0.5*sim_AoS_H.dt;

    // set kernel arguments for AdvanceVelocity
    setAVArgsAoS(AdvanceVelocityAoS, sim_AoS_D, sim_AoS_H.dt);

    // set kernel arguments for AdvancePosition
    setAPArgsAoS(AdvancePositionAoS, sim_AoS_D, dthalf);

    if (eam_flag) {
	// set the arguments for all 3 EAM_Force kernels
	setEAMArgsAoS(Force_kernelsAoS, sim_AoS_D);

    } else {
	// set kernel arguments for LJ_Force
	setLJArgsAoS(Force_kernelsAoS[0], sim_AoS_D);
    }

    // Start the simulation here;
    printf("************************************************************************\n");
    printf("Starting AoS OpenCL simulation\n");

    //ts = timeNow();

    cl_real t_kern;
    computeForceOCL(Force_kernelsAoS, &Force_event, n_global, n_local, eam_flag, sim_AoS_H.ntot, &t_kern);

    getVecAoS(sim_AoS_D.f, sim_AoS_H.f, sim_AoS_H.array_size);

    printf(" Initial:\n");
    computePrintEnergyAoS(sim_AoS_D, sim_AoS_H);

    getPrintStateAoS(sim_AoS_D, sim_AoS_H);
}


#ifdef INTEROP_VIZ
void keyboard(unsigned char key, int x, int y) { glutPostRedisplay(); }
void idle() { glutPostRedisplay(); }
void mouse(int button, int state, int x, int y) 
{
    if (state == GLUT_DOWN) mouse_buttons |= 1<<button;
    else if (state == GLUT_UP) mouse_buttons = 0;

    mouse_old_x = x;
    mouse_old_y = y;
    glutPostRedisplay();
}

void motion(int x, int y) 
{
    float dx = x - mouse_old_x;
    float dy = y - mouse_old_y;

    if (mouse_buttons == 1)
    {
	Quaternion newRotX;
	QuaternionSetEulerAngles(&newRotX, -0.2*dx*3.14159/180.0, 0.0, 0.0);
	QuaternionMul(&q, q, newRotX);

	Quaternion newRotY;
	QuaternionSetEulerAngles(&newRotY, 0.0, 0.0, -0.2*dy*3.14159/180.0);
	QuaternionMul(&q, q, newRotY);
    }
    else if (mouse_buttons == 4)
    {
	cameraFOV += dy/25.0f;
    }

    mouse_old_x = x;
    mouse_old_y = y;
    glutPostRedisplay();
}

void renderData() 
{
    if (iter == 0) ts = timeNow();
    if (iter++ < niter) computeIterationSoA();
    if (iter == niter) { te = timeNow(); finishOclSoA(); }

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    oclRender();
    glutSwapBuffers();
}
#endif


int main(int argc, char **argv) {

#ifdef INTEROP_VIZ
    q.x = q.y = q.z = 0.0f;  q.w = 1.0f;  
    mouse_old_x = 0; mouse_old_y = 0; mouse_buttons = 0; cameraFOV = 60.0f;
    g_threadStatus = 0;  g_clFinished = 0;
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1024, 1024);
    glutCreateWindow("CoMD"); 
    glutDisplayFunc(renderData);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutIdleFunc(idle);

    glewInit();
#endif

    /* get sim_flat from cmdline */
    sim = initSimFromCmdLine(argc, argv, &eam_flag, &gpu_flag);

    oclInit(gpu_flag);

#ifdef INTEROP_VIZ
    oclInitInterop(sim->nboxes);
    computeInitSoA();
    iter = 0;  niter = 10;  nsteps = 10;
    glutMainLoop();
#else
    niter = 10;  nsteps = 10;
    // SoA test loop
    computeInitSoA();
    ts = timeNow();
    for (iter=0; iter<niter; iter++) computeIterationSoA();
    te = timeNow();
    finishOclSoA();
    // AoS test loop
    computeInitAoS();
    ts = timeNow();
    for (iter=0; iter<niter; iter++) computeIterationAoS();
    te = timeNow();
    finishOclAoS();
    runRef();
#endif
}

