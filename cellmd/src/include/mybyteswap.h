#ifndef _MYBYTESWAP_H_
#define _MYBYTESWAP_H_

#ifndef NO_ENDIAN_SWAP
#include <byteswap.h>
static inline double bswap_64d(double b) {
  typedef union bswap_64d_u {
    double d;
    unsigned long long u;
  
  } bswap_64d_u;
  bswap_64d_u a;
  a.u = bswap_64(((bswap_64d_u)b).u);
  return a.d;
}
#else
#define bswap_32(x)  x
#define bswap_64(x)  x
#define bswap_64d(x) x
#endif

#endif
