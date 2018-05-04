/*****************************************************************************
 * 
 * File:    bpautils.c
 * Author:  Alex Stivala
 * Created: June 2007
 *
 * Various utility functions for message logging, memory allocation, etc.
 * 
 *
 * $Id: bpautils.c 2632 2009-07-11 03:59:11Z astivala $
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <errno.h>

#include "bpautils.h"

static int verbose = FALSE;


/*
 * write error message with varargs 
 */
void bpa_error_msg(const char *function_name, const char *format, ...)
{
  va_list ap;
  va_start(ap, format);

  fprintf(stderr, "%s: ", function_name);
  vfprintf(stderr, format, ap);
  va_end(ap);
}

/*
 * write error message with variable args, and exit with failure code 
 */
void bpa_fatal_error(const char *function_name, const char *format, ...)
{
  va_list ap;
  va_start(ap, format);

  fprintf(stderr, "%s: ", function_name);
  vfprintf(stderr, format, ap);
  va_end(ap);
  
  exit(EXIT_FAILURE);
}


/*
 * write message with variable args to stderr, if in verbose mode
 */
void bpa_log_msg(const char *function_name, const char *format, ...)
{
  va_list ap;
  va_start(ap, format);

  if (verbose)
  {
    fprintf(stderr, "%s: ", function_name);
    vfprintf(stderr, format, ap);
  }
  va_end(ap);
}

/*
 * malloc(), exiting on failure 
 */
void *bpa_malloc(size_t size)
{
  void *p = NULL;

  if (!(p = malloc(size)))
    bpa_fatal_error("bpa_malloc", "malloc failed\n");

  return p;
}

   

/*
 * calloc(), exiting on failure 
 */
void *bpa_calloc(size_t nelem, size_t size)
{
  void *p = NULL;

  if (!(p = calloc(nelem, size)))
    bpa_fatal_error("bpa_calloc", "calloc failed\n");

  return p;
}

   
/*
 * realloc(), exiting on failure 
 */
void *bpa_realloc(void *ptr, size_t size)
{
  void *p = NULL;

  if (!(p = realloc(ptr, size)))
    bpa_fatal_error("bpa_realloc", "realloc failed\n");
  
  return p;
}

/*
 * set_verbose() - set/clear the vebose flag
 *
 */
void bpa_set_verbose(bool v)
{
  verbose = v;
}


/* Subtract the `struct timeval' values X and Y,
   storing the result in RESULT.
   Return 1 if the difference is negative, otherwise 0.  
(from GNU libc manual) */
     
int
timeval_subtract (struct timeval *result, struct timeval *x, 
                  struct timeval *y)
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }
     
  /* Compute the time remaining to wait.
     tv_usec is certainly positive. */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;
     
  /* Return 1 if result is negative. */
  return x->tv_sec < y->tv_sec;
}

/*
 * Return the number of processors online
 */
int get_num_cores(void)
{
  long ncores;
  if ((ncores = sysconf(_SC_NPROCESSORS_ONLN)) < 0)
    bpa_fatal_error("sysconf() failed (%d)\n", errno);
  return (int)ncores;
}

/* Arrange the N elements of ARRAY in random order.
   Only effective if N is much smaller than RAND_MAX;
   if this may not be the case, use a better random
   number generator. */
/* http://benpfaff.org/writings/clc/shuffle.html */
void shuffle(int *array, size_t n, unsigned int *seed)
{
  if (n > 1) {
    size_t i;
    for (i = 0; i < n - 1; i++) {
      size_t j = i + rand_r(seed) / (RAND_MAX / (n - i) + 1);
      int t = array[j];
      array[j] = array[i];
      array[i] = t;
    }
  }
}



/* return a random permutation of 0 ... n-1 in array (allocated by caller) 
 * Parameters:
 *    array - (out) array of n ints allocated by caller
 *    n     - number of integers to permute (i.e. permute 0..n-1), 
 *              number of integers in array
 *    seed - (in/out) seed for rand_r()
 *  Return value: None
 */
void random_permutation(int *array, size_t n, unsigned int *seed)
{
  int i;
  for (i = 0; i < n; i++)
    array[i] = i;
  shuffle(array, n, seed);
}
