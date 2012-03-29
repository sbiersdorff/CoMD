/******************************************
 * cmdLineParser.c
 *
 * A general purpose command line parser
 * written by Sriram Swaminarayan
 *            July 24, 2007
 *****************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <memory.h>

#define nextOption(o) ((myOption_t *) o->next)

typedef struct myOption_t {
  char *help;
  char *longArg;
  unsigned char  shortArg[2];
  int   argFlag;
  char  type;
  int   sz;
  void *ptr;
  void *next;
} myOption_t;

static int longest = 1;
static struct  myOption_t *myargs=NULL;

static char *dupString(char *s) {
  char *d;
  if ( ! s ) s = "";
  d = calloc(strlen(s)+1,sizeof(char));
  strcpy(d,s);
  return d;
}
  
static myOption_t *myOptionAlloc(char *longOption,char shortOption, int has_arg, char type, void *dataPtr, int dataSize, char *help) {
  myOption_t *o;
  static int iBase=129;
  o = calloc(sizeof(myOption_t),1);
  o->help = dupString(help);
  o->longArg = dupString(longOption);
  if(shortOption) o->shortArg[0] = (unsigned char)shortOption;
  else {o->shortArg[0] = iBase; iBase++;}
  o->argFlag = has_arg;
  o->type = type;
  o->ptr = dataPtr;
  o->sz = dataSize;
  if(longOption) longest = (longest>strlen(longOption)?longest:strlen(longOption));
  return o;
}

static myOption_t *myOptionFree(myOption_t *o) {
  myOption_t *r;
  if(!o) return NULL;
  r = nextOption(o);
  if(o->longArg)free(o->longArg);
  if(o->help)free(o->help);
  free(o);
  return r;
}

static myOption_t *lastOption(myOption_t *o) {
  if ( ! o) return o;
  while(nextOption(o)) o = nextOption(o);
  return o;
}

static myOption_t *findOption(myOption_t *o, unsigned char shortArg) {
  while(o) {
    if (o->shortArg[0] == shortArg) return o;
    o = nextOption(o);
  }
  return o;
}
  

int addArg(char *longOption,char shortOption, int has_arg, char type, void *dataPtr, int dataSize, char *help) {
  myOption_t *o,*p;
  o = myOptionAlloc(longOption,shortOption,has_arg,type,dataPtr,dataSize, help);
  if ( ! o ) return 1;
  if ( ! myargs) myargs = o;
  else {
    p = lastOption(myargs);
    p->next = (void *)o;
  }
  return 0;
}
void freeArgs() {while(myargs) {myargs = myOptionFree(myargs);} return;}
void printArgs() {
  myOption_t *o = myargs;
  char s[4096];
  unsigned char *shortArg;
  printf("\n"
	 "  Arguments are: \n");
  sprintf(s,"   --%%-%ds",longest);
  while(o) {
    if(o->shortArg[0]<0xFF) shortArg = o->shortArg;
    else shortArg = (unsigned char *) "---";
    printf(s,o->longArg);
    printf(" -%c  arg=%1d type=%c  %s\n",shortArg[0],o->argFlag,o->type,o->help);
    o = nextOption(o);
    
  }
  printf("\n\n");
  return;
}

void processArgs(int argc, char **argv) {
  myOption_t *o;
  int n=0;
  int i;
  struct option *opts;
  char *sArgs;
  int c;
  if ( ! myargs) return;
  o = myargs;
  while(o) {n++,o=nextOption(o);}

  o = myargs;
  sArgs= calloc(2*(n+2),sizeof(char));
  opts = calloc(n,sizeof(struct option));
  for(i=0; i<n; i++) {
    opts[i].name = o->longArg;
    opts[i].has_arg = o->argFlag;
    opts[i].flag    = 0;
    opts[i].val     = o->shortArg[0];

    strcat(sArgs,(char *) o->shortArg);
    if(o->argFlag) strcat(sArgs,":");
    o = nextOption(o);

  }

  while(1) {

    int option_index = 0;

    c = getopt_long (argc, argv, sArgs, opts, &option_index);
    if ( c == -1) break;
    o = findOption(myargs,c);
    if ( ! o ) {
      printf("\n\n"
	     "    invalid switch : -%c in getopt()\n"
	     "\n\n",
	     c);
      break;
    }      
    if(! o->argFlag) {
      int *i = (int *)o->ptr;
      *i = 1;
    }
    else {
      switch(o->type) {
      case 'i':
	sscanf(optarg,"%d",(int *)o->ptr);
	break;
      case 'f':
	sscanf(optarg,"%f",(float *)o->ptr);
	break;
      case 'd':
	sscanf(optarg,"%lf",(double *)o->ptr);
	break;
      case 's':
	strncpy(o->ptr,optarg,o->sz);
	((char *)o->ptr)[o->sz-1] = '\0';
	break;
      case 'c':
	sscanf(optarg,"%c",(char *)o->ptr);
	break;
      default:
	printf("\n\n"
	       "    invalid type : %c in getopt()\n"
	       "    valid values are 'e', 'z'. 'i','d','f','s', and 'c'\n"
	       "\n\n",
	       c);      
      }
    }
  }

  free(opts);
  free(sArgs);
  //  freeArgs();
  return;
}







#if 0
/* a tester */
int main(int argc, char **argv) {
  int n=0;
  char infile[9]="infile value";
  char outfile[256]="outfile value";
  float fl=-1.0;
  double dbl = -1.0;
  int    flag = 0;

  //Add arguments
  addArg("infile",  'i',  1,  's',  infile,  sizeof(infile), "input file name");
  addArg("outfile", 'o',  1,  's',  outfile, sizeof(infile), "output file name");
  addArg("nSPUs",   'n',  1,  'i',  &n,      0,              "number of SPUs");
  addArg("floater",  0,   1,  'f',  &fl,     0,              "floating number");
  addArg("doubler", 'k',  1,  'd',  &dbl,    0,              "double  number");
  addArg("justFlag",'F',  0,    0,  &flag,   0,              "just a flag");

  // print the argument help
  printArgs();

  // process (and free) arguments
  processArgs(argc,argv);

  // print the variables in the code
  printf("Got n = %d\n",n);
  printf("Got fl = %f\n",fl);
  printf("Got dbl = %lf\n",dbl);
  printf("Got flag = %d\n",flag);
  printf("Got infile = %s\n",infile);
  printf("Got outfile = %s\n",outfile);

  printf("\n\n");
  
  return 0;
}
#endif
