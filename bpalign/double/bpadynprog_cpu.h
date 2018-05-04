#ifndef BPADYNPROG_CPU_H
#define BPADYNPROG_CPU_H
/*****************************************************************************
 * 
 * File:    bpadynprog_cpu.h
 * Author:  Alex Stivala
 * Created: April 2009, 
 *          from phd/bpalign/src/bpadynprog.c June 2007
 *
 * Declarations for simple CPU dp implemnentation
 * 
 *
 * $Id: bpadynprog_cpu.h 2584 2009-06-25 07:54:41Z astivala $
 *
 *****************************************************************************/

#include "bpaipsilist.h"

/*
 *  The dynamic programming (bottom up) array computation for base pair
 *  probability matrix alignment. See Hofacker et al 2004 (p. 2223).
 *  This is the basic CPU implementation that simply computes
 *  the entire dp matrix serially.
 *
 *  
 *  Parameters:   n1    - length of first sequence
 *                n2    - length of second sequence
 *                seqA    - first sequence
 *                seqB    - second sequence
 *                seripsiA  - serialized ipsilist for first seq
 *                ld_seripsiA - leading dimension of seripsiA
 *                seripsiB  - serialized ipsilist for second seq
 *                ld_seripsiB - leading dimension of seripsiB
 *                gamma - gap penalty (<= 0)
 *                M     - minimum size of hairpin loop 
 *                S     - (workarea) the dp matrix S
 *                score - (output) final score computed
 *
 */
void dynprog_cpu(int n1, int n2, 
                                const char *seqA, const char *seqB, 
                                const ipsi_element_t *seripsiA, int ld_seripsiA,
                                const ipsi_element_t *seripsiB, int ld_seripsiB,
                                double gamma, int M, double * S,
                                double *score);



/* dynamic programming (memoization) with no bounding to compute S(i,j,k,l) */
double bpa_dynprogm(int i, int j, int k, int l);

/* dynamic programming (memoization) with no bounding to compute S(i,j,k,l) 
  with array rather than hashtable to store S */
double bpa_dynprogm_array(int i, int j, int k, int l, double *S);


/* dynamic programming (memoization) with no bounding to compute S(i,j,k,l) 
 with threads */
double bpa_dynprogm_thread_master(int i, int j, int k, int l);
void *bpa_dynprogm_thread(void *threadarg);


/* dynamic programming (memoization) with no bounding to compute S(i,j,k,l) 
 with threads, using array not hashtable */
double bpa_dynprogm_thread_array_master(int i, int j, int k, int l, double *S);
void *bpa_dynprogm_thread_array(void *threadarg);

#endif /* BPADYNPROG_CPU_H */
