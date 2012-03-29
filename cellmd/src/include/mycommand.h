#ifndef MYCOMMAND_H
#define MYCOMMAND_H
#include "mytype.h"
typedef struct command_t {
  char filename[1024];
  char potdir[1024];
  char potname[1024];
  int nSPUs;
  int debug;
  int doppu;
  int doeam;
  int series;
  int periodic;
} command_t;

extern void usage(char *exename);
extern void parseCommandLine(command_t *cmd, int argc, char **argv);
extern void printCmd(command_t *cmd);

#endif
