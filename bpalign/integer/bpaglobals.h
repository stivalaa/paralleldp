#ifndef BPAGLOBALS_H
#define BPAGLOBALS_H
/*****************************************************************************
 * 
 * File:    bpaglobals.h
 * Author:  Alex Stivala
 * Created: June 2007
 *
 * Definition of globals data structure for the bpalign program.
 * 
 *
 * $Id: bpaglobals.h 2647 2009-07-13 04:17:19Z astivala $
 *
 *****************************************************************************/

#include <stdio.h>

#include "bpautils.h"
#include "bpaipsilist.h"
#include "bpaparse.h"

#define PMIN 1e-04 /* minimum base pairing probability considered significant */
#define MINLOOP 5  /* minimum size of hairpin loop */

typedef struct bpaglobals_s 
{
    /* set from command line options */

    bool   verbose;        /* if true, write trace information etc. to stderr */
    bool   printstats;     /* if true, print instrumentation counts to stderr */
    bool   useglobalbounding; /* if true, use global bounding (else none) */
    bool   exactseqscore;  /* if true, use exact sequence score as ubound */
    bool   useordering;    /* if true, try to order evaluations in inner loop */
    bool   use_bottomup;   /* if true use bottomup not topdown implemetation */
    bool   use_threading;  /* if true use threaded implementation */
    int    num_threads;    /* number of threads to use if use_threading */
    bool   use_array;      /* use array not hashtable for top-down */
    bool   use_random;     /* randomize choices in multithread version */
    FILE  *ubounddata_fp;  /* file to write ubound data for gnuplot to */

    /* constants which should probably be settable from command line (TODO) */

    myint64_t gamma;          /* gap penalty (<=0) */
    myint64_t sigma_match;    /* score for matching base */
    myint64_t sigma_mismatch; /* score for mismatched base */

    /* global data structures */

    char        *seqA;      /* first sequence */
    char        *seqB;      /* second sequence */
    int          seqlenA;   /* length of seqA */
    int          seqlenB;   /* length of seqB */
    ipsi_list_t *ipsilistA; /* (j,psi) lists indexed by i for 1st sequence */
    ipsi_list_t *ipsilistB; /* (j,psi) lists indexed by i for 2nd sequence */
    basepair_t  *pairlistA; /* list of (i,j,p) for 1st seq */
    basepair_t  *pairlistB; /* list of (i,j,p) for 2ns seq */
    int          paircountA;/* length of pairlistA */
    int          paircountB;/* length of pairlistB */
} bpaglobals_t;

extern bpaglobals_t bpaglobals;

#endif /* BPAGLOBALS_H */
