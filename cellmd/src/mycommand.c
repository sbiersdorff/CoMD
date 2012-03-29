
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include "pmd.h"
#include <time.h>
#include <cmdLineParser.h>

void usage(char *exename) {
  printf("\n\n"
	 "    Usage: %s \n"
	 "         [-e --eam ] \n"
	 "         [-v --verbose <0|1>] \n"
	 "         [-c --doppu ] \n"
	 "         [-d --potdir  potsdir] \n"
	 "         [-p --potname potname]\n"
	 "         [-n --nspus   number_spus]\n"
	 "         [-s --series  ] // runs a series of spus\n"
	 "         [-f --file    file_name]\n"
	 "         [-z --noperiodic no periodic boundary conditions]\n"
	 "\n",
	 exename);
  exit(1);
}

void printCmd(command_t *cmd) {
  printf("\n"
	 "   File = %s\n"
	 "   doeam = %d\n"
	 "   potdir = %s\n"
	 "   potname = %s\n"
	 "   nSPUS = %d\n"
	 "   debug = %d\n"
	 "   doppu = %d\n"
	 "   series = %d\n"
	 "   Periodic = %d\n",
	 cmd->filename,
	 cmd->doeam,
	 cmd->potdir,
	 cmd->potname,
	 cmd->nSPUs,
	 cmd->debug,
	 cmd->doppu,
	 cmd->series,
	 cmd->periodic);
  return;
}
void parseCommandLine(command_t *cmd, int argc, char **argv) {
  extern int PERIODIC;
  char spus[1024];
  char debug[1024];
  int digit_optind = 0;
  int c;
  int noperiodic=0;
  int help=0;

  // fill up cmd with defaults
  strcpy(cmd->filename,"data/ljtest.inp");
  strcpy(cmd->potdir,"pots");
  strcpy(cmd->potname,"ag");
  cmd->doeam = 0;
  cmd->debug = 0;
  cmd->series = 0;
  cmd->periodic = 1;


  // add arguments for processing
  addArg("help",    'h',  0,  'i',  &(help), 0,                     "print this message");
  addArg("infile",  'f',  1,  's',  cmd->filename,  sizeof(cmd->filename), "input file name");
  addArg("potdir",  'd',  1,  's',  cmd->potdir,    sizeof(cmd->potdir),   "potential directory");
  addArg("potname", 'p',  1,  's',  cmd->potname,   sizeof(cmd->potname),  "potential name");
  addArg("nspus",   'n',  1,  'i',  &(cmd->nSPUs),  0,                     "number of SPUs");
  addArg("debug",   'v',  0,  'i',  &(cmd->debug),  0,                     "verbosity level");
  addArg("doppu",   'c',  0,  'i',  &(cmd->doppu),  0,                     "compute force on PPU");
  addArg("doeam",   'e',  0,  'i',  &(cmd->doeam),  0,                     "compute force on PPU");
  addArg("series",  's',  0,  'i',  &(cmd->series), 0,                     "run a series of spus");
  addArg("noperiodic", 'z',  0,  'i',  &noperiodic, 0,                     "do not use periodic bc");

  processArgs(argc,argv);
  if(help) {
    printArgs();
    freeArgs();
    exit(2);
  }
  freeArgs();
  cmd->periodic = !(noperiodic);
  PERIODIC=cmd->periodic;
  return;
}

