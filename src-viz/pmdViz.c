/**
 * a simple md simulator
 **/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

#include "pmd.h"

#include <pthread.h>


void *computation_thread(void *ptr)
/* do the computation work */
{
  SimThreadData *data = (SimThreadData*)ptr;
  return do_compute_work(data);
}


int main(int argc, char **argv) {
  command_t cmd;
  simflat_t *sim;
  SimThreadData* threadData;

  /* get sim_flat from cmdline */
  sim = initSimFromCmdLine(argc, argv);

  writeClsman(sim,(char *) "init.bin");

  /* start computation thread */
  if(1) {
    pthread_t computationThread;
    pthread_t visualizationThread;
    AtomVisualize* viz = new AtomVisualize();

    viz->initialize(sim);
    viz->updateData(sim);
    viz->renderData();
      
    threadData = new SimThreadData;
    threadData->sim = sim;
    threadData->viz = viz;
		  
    pthread_create(&computationThread, NULL, computation_thread, (void*) threadData);
    viz->interact();
  }

  writeClsman(sim,(char *) "final.bin");

  /* release sim_flat */
  destroySimulation(&sim); 	

  return 0;
}
