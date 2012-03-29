#ifndef _MEMUTILS_H_

#define _MEMUTILS_H_

#include <unistd.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>

#ifndef DMA_PAD 
#define DMA_PAD(iSize) (((iSize)+15) & ~15)
#endif

extern void breakinme();
extern int suAbort(int status, char *msg);
  
    
/**
 * use posix_memalign to alloc **/
    
#define ALIGNED_MALLOC(ptr,size) {					\
    if(posix_memalign((void **) &ptr,128,DMA_PAD(size))){		\
      char t[1024];sprintf(t,"unable to dingo malloc %d\n",size);	\
      breakinme();							\
      suAbort(195,t);							\
    }									\
    memset(ptr,0,size);							\
  }
    
#define ALIGNED_FREE(ptr) {if(ptr) free(ptr);}
  
  
#define ALIGNED_CALLOC(ptr,size) {ALIGNED_MALLOC(ptr,size);memset(ptr,0,size);}
  
static void *suAlignedCalloc(int isize) {
  void *b;
  ALIGNED_CALLOC(b,isize);
  return b;
}
static void *suAlignedMalloc(int isize) {
  void *b;
  ALIGNED_MALLOC(b,isize);
  return b;
}
static void suAlignedFree(void *d) {
  ALIGNED_FREE(d);
}

#endif
