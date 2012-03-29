#ifndef __MYTYPE_H_
#define __MYTYPE_H_

#ifdef SINGLE
typedef float real_t;
  #define FMT1 "%g"
#else
typedef double real_t;
  #define FMT1 "%lg"
#endif

#endif
