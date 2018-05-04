/*****************************************************************************
 * 
 * File:    bpaglobals.c
 * Author:  Alex Stivala
 * Created: June 2007
 *
 * Globals data definition for the bpalign program.
 * 
 *
 * $Id: bpaglobals.c 2647 2009-07-13 04:17:19Z astivala $
 *
 *****************************************************************************/

#include "bpaglobals.h"

bpaglobals_t bpaglobals = 
{
  FALSE  /* verbose */
  ,FALSE /* prinstats */
  ,FALSE /* useglobalbounding */
  ,FALSE /* exactseqscore */
  ,FALSE /* useordering */
  ,FALSE /* use_bottomup */
  ,FALSE /* use_threading */
  ,0     /* num_threads */
  ,FALSE /* use_array */
  ,TRUE  /* use_random */
  ,NULL  /* ubounddata_fp */

  ,-60*SIGMA_MATCH     /* gamma */  
  ,SIGMA_MATCH   /* sigma_match */ /* FIXME unused */
  ,0      /* sigma_mismatch */
  
  ,NULL   /* first sequence */
  ,NULL   /* 2nd sequence */
  ,0      /* seqlenA */
  ,0      /* seqlenB */
  ,NULL   /* ipsilist A */
  ,NULL   /* ipsilist B */
  ,NULL   /* pairlistA */
  ,NULL   /* pairlistB */
  ,0      /* paircountA */
  ,0      /* paircountB */
};
