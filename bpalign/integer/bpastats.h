#ifndef BPASTATS_H
#define BPASTATS_H
/*****************************************************************************
 * 
 * File:    bpastats.h
 * Author:  Alex Stivala
 * Created: June 2007
 *
 * Declaration of global data structure for keeping instrumentation of
 * calls, cell counts etc. for dynamic programming in bpalign program.
 * 
 *
 * $Id: bpastats.h 2584 2009-06-25 07:54:41Z astivala $
 *
 *****************************************************************************/

typedef unsigned long counter_t;

typedef struct bpastats_s 
{
    counter_t count_dynprogm_entry;  /* calls of dynprogm */
    counter_t count_dynprogm_entry_notmemoed; /* where not memo value return */
    counter_t count_S;               /* matrix S cells computed */
    counter_t count_U;               /* matrix U cells computed */
    counter_t count_recomputations; /* global bounding recomputations */
} bpastats_t;



#endif /* BPASTATS_H */
