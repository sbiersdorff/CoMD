/**
 * a simple md simulator
 **/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <pthread.h>

#include "pmdOCL.h"
#include "helpers.h"


#ifdef INTEROP_VIZ

#include "glew.h"

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
cl_kernel *force_kernels;

int err;
size_t n_local[2];
size_t n_global[2];

cl_event Force_event;
cl_event AP_event;
cl_event AV_event;

real_t t_simple;
real_t t_overall;
real_t t_local;
real_t t_dummy;
real_t dt;

host_sim_t sim_H;
dev_sim_t sim_D;

double ts, te;

int eam_flag;
int gpu_flag;
simflat_t *sim;
int iter, ns;
int niter, nsteps;


void compute_iteration()
{
    for (ns = 0;ns < nsteps; ns++) {
        // advance particle positions dt/2
        oclRunKernel(AdvancePosition, &AP_event, n_global, n_local);

        if (DIAG_LEVEL > 0)
            GetPrintState(sim_D, sim_H);

        // compute force
        computeForceOCL(force_kernels, Force_event, n_global, n_local, eam_flag);

        if (DIAG_LEVEL > 0)
            GetPrintState(sim_D, sim_H);

        // advance velocity a full timestep
        oclRunKernel(AdvanceVelocity, &AV_event, n_global, n_local);

        if (DIAG_LEVEL > 0)
            GetPrintState(sim_D, sim_H);

        // advance particle positions dt/2
        oclRunKernel(AdvancePosition, &AP_event, n_global, n_local);

        if (DIAG_LEVEL > 0)
            GetPrintState(sim_D, sim_H);
    }

#ifdef INTEROP_VIZ 
    oclGraphics(Viz, sim_D, n_global, n_local);
#endif

    // compute force
    computeForceOCL(force_kernels, Force_event, n_global, n_local, eam_flag);

    printf(" after dt (-%d-): ", nsteps*iter + 1);
    ComputePrintEnergy(sim_D, sim_H);
}


void cleanup_ocl()
{
    GetElapsedTime(AP_event, &t_simple, &t_overall);
    printf("Kernel AdvancePosition executed in %e secs.\n", t_simple);
    GetElapsedTime(AV_event, &t_simple, &t_overall);
    printf("Kernel AdvanceVelocity executed in %e secs.\n", t_simple);

    real_t t_ref = (real_t)(te - ts);

    printf("Reference elapsed time = %e (%e us/atom for %d atoms)\n", t_ref,1.0e6*t_ref/(double)(sim->ntot*niter*nsteps), sim->ntot);

    // copy result back from device
    // note this utility forces synchronous copy so 
    // it does not return until the copy is complete
    GetVec(sim_D.f, sim_H.f, sim_H.array_size);
    oclCopyToHost(sim_D.e, sim_H.e, sim_H.array_size);

    ComputePrintEnergy(sim_D, sim_H);

    // clean up the OpenCL memory 
    FreeSims(sim_H, sim_D);

    clReleaseEvent(AP_event);
    clReleaseEvent(AV_event);

    oclCleanup();

    // write initial state 
    writeClsman(sim,(char *) "init.bin");

    // do the computation 
    (void)do_compute_work(sim); 

    // write final configuration 
    writeClsman(sim,(char *) "final.bin");

    // free memory 
    destroySimulation(&sim);
}



void compute_init()
{
    sim_D.rmass = (real_t)(amu_to_m_e)*(double)(sim->pot->mass);
    sim_H.eam_flag = eam_flag;

    if(eam_flag) {
        force_kernels = malloc(sizeof(cl_kernel)*3);
    } else {
        force_kernels = malloc(sizeof(cl_kernel)*1);
    }

    // set up arrays for OpenCL

    HostSimInit(&sim_H, sim);

    DevSimInit(&sim_D, sim_H);

    PutSim(sim_H, sim_D);

    // build the program from the kernel source file

    BuildModules(force_kernels, &AdvancePosition, &AdvanceVelocity, &Viz, sim_H, n_local, n_global);

    real_t dthalf = sim_H.dt/2.0;
    real_t dtminushalf = -0.5*sim_H.dt;

    if (eam_flag) {
        // set the arguments for all 3 EAM_Force kernels
        SetEAMArgs(force_kernels, sim_D);

    } else {
        // set kernel arguments for LJ_Force
        SetLJArgs(force_kernels[0], sim_D);
    }

    // set kernel arguments for AdvanceVelocity
    SetAVArgs(AdvanceVelocity, sim_D, sim_H.dt);

    // set kernel arguments for AdvancePosition
    SetAPArgs(AdvancePosition, sim_D, dthalf);

    // Start the simulation here;
    printf("Starting simulation\n");

    //ts = timeNow();

    computeForceOCL(force_kernels, Force_event, n_global, n_local, eam_flag);

    GetVec(sim_D.f, sim_H.f, sim_H.array_size);

    printf(" Initial:\n");
    ComputePrintEnergy(sim_D, sim_H);

    GetPrintState(sim_D, sim_H);
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
  if (iter++ < niter) compute_iteration();
  if (iter == niter) { te = timeNow(); cleanup_ocl(); }

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
    glutCreateWindow("cruft"); 
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
    compute_init();
    iter = 0;  niter = 10;  nsteps = 10;
    glutMainLoop();
#else
    compute_init();
    niter = 10;  nsteps = 10;
    ts = timeNow();
    for (iter=0; iter<niter; iter++) compute_iteration();
    te = timeNow();
    cleanup_ocl();
#endif
}

