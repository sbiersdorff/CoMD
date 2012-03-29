#define RUNTIME_API 1

#ifndef CL_UTILS_H
#define CL_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if defined (__APPLE__) || defined(MACOSX)
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "mytype.h"

#ifdef SINGLE
#define cl_real cl_float
#else
#define cl_real cl_double
#endif


// OpenCL stuff to make syntax simpler for 'Runtime API' constructs
//
#if (RUNTIME_API)

extern cl_uint          num_platforms;
extern cl_platform_id*  platform_id; 
extern cl_device_id     device_id;
extern cl_context       context; 
extern cl_command_queue commandq;
extern cl_program       program;

extern int oclInit(int gpu_request);

extern void oclInitInterop();

extern int oclCleanup();

extern int oclBuildProgramFromFile(char* filename);

extern void oclCreateReadBuffer(cl_mem* device_buffer, int size);

extern void oclCreateWriteBuffer(cl_mem* device_buffer, int size);

extern void oclCopyToDevice(void* host_buffer, cl_mem device_buffer, int size);

extern void oclCopyToHost(cl_mem device_buffer, void* host_buffer, int size);

#else 

extern int InitCL(cl_uint* num_platforms, 
        cl_platform_id* platform_id, 
        cl_device_id* device_id,
        cl_context* context, 
        cl_command_queue* commandq);

extern int BuildProgramFromFile(cl_program* program, 
        char* filename, 
        cl_context context, 
        cl_device_id device_id);

#endif

int ParseGPURequest(int argc, char* argv[]);

extern char* print_cl_errstring(int err);

#endif
