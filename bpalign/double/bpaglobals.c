/*****************************************************************************
 * 
 * File:    bpaglobals.c
 * Author:  Alex Stivala
 * Created: June 2007
 *
 * Globals data definition for the bpalign program.
 * 
 *
 * $Id: bpaglobals.c 2584 2009-06-25 07:54:41Z astivala $
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
  ,NULL  /* ubounddata_fp */

  ,-3     /* gamma */  
  ,0.05   /* sigma_match */
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
