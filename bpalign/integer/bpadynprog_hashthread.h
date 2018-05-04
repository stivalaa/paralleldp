#ifndef BPADYNPROG_HASHTHREAD_H
#define BPADYNPROG_HASHTHREAD_H
/*****************************************************************************
 * 
 * File:    bpadynprog_hashthread.h
 * Author:  Alex Stivala
 * Created: April 2009, 
 *
 * Declarations for implementations of threaded d.p. using hash table
 * (used for multiple different hashtable implementations).
 * 
 *
 * $Id: bpadynprog_hashthread.h 2584 2009-06-25 07:54:41Z astivala $
 *
 *****************************************************************************/


/* dynamic programming (memoization) with no bounding to compute S(i,j,k,l) 
 with threads */
int64_t bpa_dynprogm_thread_master(int i, int j, int k, int l);
void *bpa_dynprogm_thread(void *threadarg);

#endif /* BPADYNPROG_HASHTHREAD_H */
