#ifndef BPAUTILS_H
#define BPAUTILS_H
/*****************************************************************************
 * 
 * File:    bpautils.h
 * Author:  Alex Stivala
 * Created: June 2007
 *
 * Declarations of various utility functions for message logging, etc.
 * 
 *
 * $Id: bpautils.h 2632 2009-07-11 03:59:11Z astivala $
 *
 *****************************************************************************/

#include <stdlib.h>
#include <sys/time.h>

#ifdef __cplusplus
extern "C" {
#endif

#define bool int

#define FALSE 0
#define TRUE  1

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))


#define EPSILON  1E-13  /* epsilon for double comparisons */


/* maximum number of threads that can be used */
#ifndef MAX_NUM_THREADS
#define MAX_NUM_THREADS 256
#endif


/*
 * We will use this to mark unset cells in the matrix and for 'lowest' bound
 */
#define INF 999999     /* dodgy */
#define NEGINF (-INF)  /* dodgy */


/* definition of a type for a 4-tuple */
/* 
 * The key for the hashtable is a 4-tuple, i.e. we look up by an (i,j,k,l)
 * tuple, where i,j,k,l would have been the indices of our 4-dimensional matrix
 */
typedef struct tuple4_s 
{
  unsigned short i;
  unsigned short j;
  unsigned short k;
  unsigned short l;
} tuple4_t;



/* write error message with varargs, and exit with failure code */
void bpa_fatal_error(const char *function_name, const char *format, ...);

/* write error message with varargs */
void bpa_error_msg(const char *function_name, const char *format, ...);

/* write message to stderr, if in verbose mode */
void bpa_log_msg(const char *function_name, const char *format, ...);

/* malloc(), exiting on failure */
void *bpa_malloc(size_t size);

/* calloc(), exiting on failure */
void *bpa_calloc(size_t nelem, size_t size);

/* realloc(), exiting on failure */
void *bpa_realloc(void *ptr, size_t size);

/* set/clear the verbose flag for bpa_log_msg */
void bpa_set_verbose(bool v);


/* timeval subtraction: result <- x - y */     
int
timeval_subtract (struct timeval *result, struct timeval *x, 
                  struct timeval *y) ;

/* return number of processors available on machine */
int get_num_cores(void);

/* arrange the n elements of array in random order  (seed is for rand_r()) */
void shuffle(int *array, size_t n, unsigned int *seed);

/* return a random permutation of 0 ... n-1 in array (allocated by caller) */
void random_permutation(int *array, size_t n, unsigned int *seed);

#ifdef __cplusplus
}
#endif

#endif /* BPAUTILS_H */
