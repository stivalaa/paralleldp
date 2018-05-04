#ifndef TBBHASHMAP_H
#define TBBHASHMAP_H
/*****************************************************************************
 * 
 * File:    tbbhashmap.h
 * Author:  Alex Stivala
 * Created: May 2009
 *
 * C callable interface to Intel Thread Building Blocks conccurrent hashmap.
 * 
 * $Id: tbbhashmap.h 2399 2009-05-16 04:31:56Z astivala $
 *
 * Developed on:
 * Linux charikar.csse.unimelb.edu.au 2.6.9-78.0.5.ELsmp #1 SMP Wed Sep 24 05:40:24 EDT 2008 x86_64 x86_64 x86_64 GNU/Linux
 * TBB21_INSTALL_DIR=/home/charikar/pgrad/astivala/tbb21_20080605oss
 * . ~/tbb21_20080605oss/em64t/cc3.4.3_libc2.3.4_kernel2.6.9/bin/tbbvars.sh
 * TBB_ARCH_PLATFORM=em64t/cc3.4.3_libc2.3.4_kernel2.6.9
 *
 *****************************************************************************/

typedef struct _twothings {
  long long high, low;
} _SET;

#ifdef __cplusplus
#define EXTERN_C extern "C"
#else
#define EXTERN_C
#endif

EXTERN_C void tbbhashmap_insert(_SET key, int value);

EXTERN_C int tbbhashmap_lookup(_SET key);

EXTERN_C int tbbhashmap_haskey(_SET key);

#endif /* TBBHASHMAP_H */

