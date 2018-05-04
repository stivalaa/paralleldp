#ifndef BPACOMMON_H
#define BPACOMMON_H
/*****************************************************************************
 * 
 * File:    bpacommon.h
 * Author:  Alex Stivala
 * Created: April 2009
 *
 * Definitions of types, constants etc. used commonly.
 * 
 *
 * $Id: bpacommon.h 2605 2009-07-04 06:55:41Z astivala $
 *
 *****************************************************************************/

/*#include <values.h>*/

#ifdef SOLARIS
#include <inttypes.h>
typedef int64_t myint64_t;
#else
typedef long long  myint64_t;
#endif

#define MAX_INT64 9223372036854775807LL

/* Index into 4d array stored in contiguous memory */
#define INDEX4D(i,k,j,l,m,n) ( (((i)*(n) + (j))*(m) + (k))*(n) + (l) )


/*
   #define INDEX4D(i,j,k,l,zwidth,wwidth,xwidth,ywidth) ((i) + ((xwidth) * (j)) + ((xwidth) * (ywidth) * (k)) + ((xwidth) * (ywidth) * (zwidth) * (l)))
*/


/* integer absolute value */
#define INTEGER_ABS(x) ((x < 0) ? -(x) : x)

/* the master thead has id 0 (don't change this; we assume first in array) */
#define MASTER_THREAD_ID 0

/* 
 * BPA_SIGMA
 *
 * return the sigma (unpaired base substitution) score for the two bases b1,b2 
 *
 * Parameters:
 *    b1  -  one base as a character (A,C,G,U)
 *    b2  -  second base as above
 *
 * Return value:
 *   The unpaired base substitution score for b1, b2.
 *   This is global parameter value sigma_match if they are same,
 *   and sigma_mismatch if they differ.
 *
 */
#define SIGMA_MATCH    1
#define SIGMA_MISMATCH 0
#define BPA_SIGMA(b1,b2) ((b1) == (b2) ? SIGMA_MATCH : SIGMA_MISMATCH)



typedef struct thread_data_s
{
    int thread_id; /* id of this thread to index into per-thread arrays */

    /* (i,j,k,l) for this thread to start at 
     *   i     - left co-ord in first sequence
     *   j     - right co-ord in first sequence
     *   k     - left co-ord in second sequence
     *   l     - right co-ord in second sequence
     *   0 <= i < j <= n1 - 1
     *   0 <= k < l <= n2 - 1
     */
    short i;
    short j;
    short k;
    short l;
    myint64_t *S;  /* The 4d d.p. matrix S allocated linearly 
                   (used for array version  only)*/
    myint64_t score;  /* OUTPUT score computed by this thread */
} thread_data_t;



#endif /* BPACOMMON_H */
