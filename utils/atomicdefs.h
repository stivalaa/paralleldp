#ifndef ATOMICDEFS_H
#define ATOMICDEFS_H
/*****************************************************************************
 * 
 * File:    atomicdefs.h
 * Author:  Alex Stivala
 * Created: April 2009
 *
 * Deefinitions of macros for atomic operations on different platforms
 * 
 * gcc version 4.1.0 or greater is required, in order to use the
 * __sync_val_compare_and_swap()
 * and __sync_val_fetch_and_add() builtins
 * (this module was developed on Linux 2.6.22 (x86) with gcc 4.1.3),
 * except on Solaris, where we use SUNWSpro compiler and
 * cas32()/caslong()/cas64()/casptr() in 
 * /usr/include/sys/atomic.h)
 *
 *
 * $Id: atomicdefs.h 2506 2009-06-11 08:12:43Z astivala $
 *
 *****************************************************************************/

#ifdef SOLARIS
#include <atomic.h>
#define CASPTR(ptr,oldval,newval) atomic_cas_ptr(ptr, oldval, newval)
#define CAS64(ptr,oldval,newval) atomic_cas_64(ptr, oldval, newval)
#define CAS32(ptr,oldval,newval) atomic_cas_32(ptr, oldval, newval)
#define ATOMIC_OR_64(ptr, x) atomic_or_64(ptr, x)
#else
#define CASPTR(ptr,oldval,newval) __sync_val_compare_and_swap(ptr, oldval, newval)
#define CAS64(ptr,oldval,newval) __sync_val_compare_and_swap(ptr, oldval, newval)
#define CAS32(ptr,oldval,newval) __sync_val_compare_and_swap(ptr, oldval, newval)
#define ATOMIC_OR_64(ptr ,x) __sync_fetch_and_or(ptr, x)
#endif

#endif /* ATOMICDEFS_H */

