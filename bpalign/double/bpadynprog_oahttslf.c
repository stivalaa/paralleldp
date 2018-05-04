/*****************************************************************************
 *
 * File:    bpadynprog_oahttslf.c
 * Author:  Alex Stivala
 * Created: April 2009, 
 *          from phd/bpalign/src/bpadynprog.c June 2007
 *
 * $Id: bpadynprog_oahttslf.c 2693 2009-07-24 01:48:03Z astivala $
 *
 * pthreads implementation using the oahttslf lockfree hashtable of
 * RNA base pair probability matrix alignment by dynamic programming.
 * 
 * Algorithm from Hofacker et al 2004
 * "Alignment of RNA base pairing probability matrices"
 * Bioinformatics 20(14):2222-2227
 * and Hofacker's prototype implementation of PMcomp and PMmulti
 * pmcomp.pl pmmulti.pl from http://www.tbi.univie.ac.at/~ivo/RNA/PMcomp/
 *
 * For a more recent and very efficient implementation (in C++) see LocARNA
 * http://www.bioinf.uni-freiburg.de/Software/LocARNA/
 * and the paper
 * Will et al 2007 "Inferring Noncoding RNA Families and Classes by
 * Means of Genome-Scale Structure-Based Clustering"
 * PLoS Computational Biology 3(4):e65
 *
 *****************************************************************************/

#include <assert.h>
#include <math.h>
#include <pthread.h>

#include "oahttslf.h"
#include "bpacommon.h"
#include "bpaglobals.h"
#include "bpastats.h"
#include "bpautils.h"
#include "bpaipsilist.h"
#include "bpadynprog_hashthread.h"


/* if the SYNCH preprocessor symbol is defined then we use threads
   but just join them immediately so it was just as if called synchronously,
   for testing */


/*****************************************************************************
 *
 * static data
 *
 *****************************************************************************/

#ifdef USE_INSTRUMENT
/* instrumentation is per-thread, each thread only writes to its own element */
static volatile bpastats_t bpastats[MAX_NUM_THREADS];
#endif


/* structures passed to each thread  as parameter */
static thread_data_t thread_data[MAX_NUM_THREADS];

/* pthreads thread handles */
static pthread_t threads[MAX_NUM_THREADS];

/* the number of active threads. This is only used by the master thread */
static int num_active_threads = 1; /* we count the master thread */ 



/*****************************************************************************
 *
 * static functions
 *
 *****************************************************************************/


/* insert by (i,j,k,l) into table */
static void oahttslf_insert_indices(uint16_t i, uint16_t j, 
                                    uint16_t k, uint16_t l,
                                    double value, int thread_id);

/* lookup by (i,j,k,l) */
static double oahttslf_lookup_indices(uint16_t i, uint16_t j,
                                    uint16_t k, uint16_t l);

/* dodgy: 0 is the OAHTTSLF_EMPTY_KEY and OAHTTSLF_EMPTY_VALUE value,
   so if key or value is 0 set it to MAGIC_ZERO instead */
#define MAGIC_ZERO 0xffffffffffffffff




/*
 * oahttslf_insert_indices()
 *
 * Insert value for (i,j,k,l) into the hashtable
 *
 * Parameters:
 *    i,j,k,l - indices to build insertion key
 *    value - value to insert for the key
 *    thread_id - id (0,...n, not pthread id) of this thread
 *
 * Return value:
 *    None.
 */
static void oahttslf_insert_indices(uint16_t i, uint16_t j, 
                                    uint16_t k, uint16_t l,
                                    double value,
                                    int thread_id)

{
  uint64_t key;
  double val;

  key = (i == 0 && j == 0 && k == 0 && l == 0 ? MAGIC_ZERO : 
         ((uint64_t)i << 47) | ((uint64_t)(j & 0xffff) << 31) |
         ((uint64_t)(k & 0xffff) << 15) | (uint64_t)(l & 0xffff));
//  printf("zz %d,%d,%d,%d\t\t%016lX\n",i,j,k,l,key);
  val = (value == 0 ? NEGINF : value);
  oahttslf_insert_double(key, val, thread_id);
}



/*
 * oahttslf_lookup_indices()
 *
 * Get the value for (i,j,k,l) from the hashtable
 *
 * Parameters:
 *     i,j,k,l - indices to build key for lookup
 * 
 * Return value:
 *     value if key found, else NEGINF
 */
static double oahttslf_lookup_indices(uint16_t i, uint16_t j,
                                    uint16_t k, uint16_t l)
{
  uint64_t key;
  double val;
  bool found;
  key = (i == 0 && j == 0 && k == 0 && l == 0 ? MAGIC_ZERO : 
         ((uint64_t)i << 47) | ((uint64_t)(j & 0xffff) << 31) |
         ((uint64_t)(k & 0xffff) << 15) | (uint64_t)(l & 0xffff));

  found =  oahttslf_lookup_double(key, &val);
  if (found)
    return (val <= NEGINF ? 0 : val);
  else
    return NEGINF;
}



/*
 * bpa_dynprogm_thread_call()
 *
 * Utility function to start a new thread running
 * bpa_dynprogm_thread() if we we are master thread and have not
 * reached the thread limit, otherwise call the function in this
 * thread
 *
 *      Parameters:   thread_id - our thread id for calling thread
 *                    num_active (read/write) - 
 *                                         current number of threads active
 *                    new_thread_id (r/w) - thread id for new thread (or  own)
 *                    i     - left co-ord in first sequence
 *                    j     - right co-ord in first sequence
 *                    k     - left co-ord in second sequence
 *                    l     - right co-ord in second sequence
 *                            0 <= i < j <= n1 - 1
 *                            0 <= k < l <= n2 - 1
 *      Uses global data:
 *                  read/write:
 *                     thread_data - structures passed to each thread
 *                     threads     - thread handles
 *                  readonly:
 *                    num_threads - number of threads we can use
 *
 *
 * Return value: None.
 * 
 *    The caller has to get the results from the hashtable at (i,j,k,l)
 *    when the thread (if there is one) has finished.
 */
static void bpa_dynprogm_thread_call(int thread_id, 
                                     int *num_active,
                                     int *new_thread_id,
                                     int i, int j, int k, int l)
{
  static const char *funcname = "bpa_dynprogm_thread_call";

  thread_data_t dummy_thread_data; /* for passing parameter in same thread */
  int rc;
  bool is_master = (thread_id == MASTER_THREAD_ID);

  if (is_master && *num_active < bpaglobals.num_threads)
  {
    fprintf(stderr,"starting thread id %d\n", *num_active);
    *new_thread_id = *num_active;
    thread_data[*num_active].thread_id = *num_active;
    thread_data[*num_active].i = i;
    thread_data[*num_active].j = j;
    thread_data[*num_active].k = k;
    thread_data[*num_active].l = l;
    if ((rc = pthread_create(&threads[*num_active], NULL,
                             bpa_dynprogm_thread, 
                             (void *)&thread_data[*num_active])))
      bpa_fatal_error(funcname, "pthread_create() failed (%d)\n", rc);
    (*num_active)++;
  }
  else
  {
    *new_thread_id = thread_id;
    dummy_thread_data.thread_id = thread_id;
    dummy_thread_data.i = i;
    dummy_thread_data.j = j;
    dummy_thread_data.k = k;
    dummy_thread_data.l = l;
    bpa_dynprogm_thread((void *)&dummy_thread_data);
  }
}


/*****************************************************************************
 *
 * external functions
 *
 *****************************************************************************/


/*
 * bpa_dynprogm_thread()
 *
 *
 *      This version is multi-threaded, sharing hashtable used to
 *      store computed S values between the threads.
 *      If our thread_id is MASTER_THREAD_ID then we are the master
 *      thread, which starts at the final value to compute.
 *      Only the master thread creates other threads, and it waits
 *      for them to finish.
 *
 *      The dynamic programming (top down memoization)
 *      computation for base pair probability matrix alignment.
 *      See Hofacker et al 2004 (p. 2223).
 *
 *      This version is the memory function version; instead of computing
 *      the whole S array bottom-up, it is computed recursively top-down
 *      and values stored in it, and reused (memoization) if already
 *      computed. 
 *
 *      This version uses no bounding.
 *
 *
 *      Is is also required that the ipsilist[i] lists are sorted
 *      by right co-ord ascending. This is assured by input processing.
 *
 *
 *      Parameters:   threadarg - thread data for this thread
 *
 *
 *      Uses global data:
 *                  read/write:
 *                    thread_data - structures passed to each thread
 *                    threads     - thread handles
 *                    hashtable in ht.c
 *                    count_dynprogm_entry[thread_id] - count of calls here
 *                    count_dynprogm_notmemoed[thread_id] -
 *                                              where not return memoed value
 *
 *                  readonly:
 *                    seqA    - first sequence
 *                    seqB    - second sequence
 *                    seqlenA - length of first sequence
 *                    seqlenB - length of second sequence
 *                    ipsilistA - (j,psi) lists indexed by i for 1st seq
 *                    ipsilistB - (j,psi) lists indexed by i for 2nd seq
 *                    gamma   - gap penalty (<= 0)
 *                    num_threads - number of threads we can use
 *
 *      Return value: NULL
 *                     (Declared void * for pthreads)
 *
 */
void *bpa_dynprogm_thread(void *threadarg)
{
  static const char *funcname = "bpa_dynprogm_thread";
  double score = NEGINF;
  double gapA, gapB, unpaired, gapmax, pairedscore;
  double sigma_ik;
  int h,q; /* h and q are the pairing co-ords used in the recurrence */
  int x,y; /* just loop indices, no meaning */
  double psiA_ih, psiB_kq; /* psi value at seqA[i,h] and seqB[k,q] */
  double sm, shq, max_shq;

  int i,j,k,l;
  thread_data_t *mydata = (thread_data_t *)threadarg;
  bool is_master = (mydata->thread_id == MASTER_THREAD_ID);
  bool comp_gapA, comp_gapB, comp_unpaired;




  int rc;
  int t;
  int active_threadids[MAX_NUM_THREADS];
  int actindex;

  i = mydata->i;
  j = mydata->j;
  k = mydata->k;
  l = mydata->l;


  assert(i >= 0);
  assert(i < bpaglobals.seqlenA);
  assert(j >= 0);
  assert(j < bpaglobals.seqlenA);
  assert(i <= j);
  assert(k >= 0);
  assert(k < bpaglobals.seqlenB);
  assert(l >= 0);
  assert(l < bpaglobals.seqlenB);
  assert(k <= l);

  bpa_log_msg(funcname, "%d\t\t%d\t%d\t%d\t%d\n",mydata->thread_id,i,j,k,l);

#ifdef USE_INSTRUMENT
  bpastats[mydata->thread_id].count_dynprogm_entry++;
#endif

  /* memoization: if value here already computed then do nothing */
  if (oahttslf_lookup_indices(i, j, k, l) > NEGINF)
    return NULL;


#ifdef USE_INSTRUMENT
  bpastats[mydata->thread_id].count_dynprogm_entry_notmemoed++;
#endif

  /*
   *  Initialization cases for the d.p. matrix S:
   *    S(i,j,k,l) = |(j-i) - (l-k)| * gamma, for j - i <= MinLoop + 1 or
   *                                              l - k <= MinLoop + 1
   */
  if ((j - i) <= MINLOOP + 1 || (l - k) <= MINLOOP + 1)
  {
    score = fabs((j - i) - (l - k)) * bpaglobals.gamma;
    bpa_log_msg(funcname, "%d\tI\t%d\t%d\t%d\t%d\t%lg\n",mydata->thread_id,i,j,k,l,score);
    oahttslf_insert_indices(i, j, k, l, score);
/*    assert(oahttslf_lookup_indices(i,j,k,l) == score); */
#ifdef USE_INSTRUMENT
    bpastats[mydata->thread_id].count_S++;
#endif
    return NULL;
  }

  /*
   *
   *     Recursive cases for computing elements of the dp matrix S
   *
   *     Compute each of the four cases over which we
   *     choose the max:
   *      1. gap in second sequence (B): S(i+1,j,k,l) + gamma
   *      2. gap in first sequence (A):  S(i,j+1,k,l) + gamma
   *      3. extension of both subsequences with unpaired position:
   *         S(i+1,j,k+1,l) + sigma(Ai, Bk)
   *      4. (more complex case with another max):
   *           match of pair (i, h) in A with (k,q) in B
   */

  actindex = 0;

  if (i + 1 < bpaglobals.seqlenA && i + 1 < j)
  {
    bpa_dynprogm_thread_call(mydata->thread_id,
                             &num_active_threads,
                             &active_threadids[actindex],
                             i + 1, j, k, l);
    if (active_threadids[actindex] != mydata->thread_id)
    {
      fprintf(stderr, "B %d: %d\n", actindex, active_threadids[actindex]);
      actindex++;
    }
    
    comp_gapB = TRUE;
    
#ifdef SYNCH
    if (comp_gapB_threadid != mydata->thread_id)
      if ((rc = pthread_join(threads[comp_gapB_threadid], NULL)))
        bpa_fatal_error(funcname, "pthread_join failed (%d)\n", rc);
#endif

  }
  else
  {
    gapB = NEGINF;
    comp_gapB = FALSE;
  }

  if (k + 1 < bpaglobals.seqlenB && k + 1 < l)
  {
    bpa_dynprogm_thread_call(mydata->thread_id,
                             &num_active_threads,
                             &active_threadids[actindex],
                             i, j, k + 1, l);
    if (active_threadids[actindex] != mydata->thread_id)
    {
      fprintf(stderr, "A %d: %d\n", actindex, active_threadids[actindex]);
      actindex++;
    }

    comp_gapA = TRUE;

#ifdef SYNCH
    if (comp_gapA_threadid != mydata->thread_id)
      if ((rc = pthread_join(threads[comp_gapA_threadid], NULL)))
        bpa_fatal_error(funcname, "pthread_join failed (%d)\n", rc);
#endif

  }
  else
  {
    comp_gapA = FALSE;
    gapA = NEGINF;
  }

  if (i+1 < bpaglobals.seqlenA && i+1 < j && k+1 < bpaglobals.seqlenB && k+1 < l)
  {
     bpa_dynprogm_thread_call(mydata->thread_id,
                              &num_active_threads,
                              &active_threadids[actindex],
                              i+1, j, k+1, l);
    if (active_threadids[actindex] != mydata->thread_id)
    {
      fprintf(stderr, "C %d: %d\n", actindex, active_threadids[actindex]);
      actindex++;
    }

     comp_unpaired = TRUE;

#ifdef SYNCH
    if (comp_unpaired_threadid != mydata->thread_id)
      if ((rc = pthread_join(threads[comp_unpaired_threadid], NULL)))
        bpa_fatal_error(funcname, "pthread_join failed (%d)\n", rc);
#endif

  }
  else
  {
    comp_unpaired = FALSE;
    unpaired = NEGINF;
  }
  

  /*
   * max_shq = max{h<=j,q<=l}( S^M[i,h,k,q] + S[h+1,j,q+1,l] )
   *           where S^M[i,j,k,l] = S[i+1,j+1,k+1,l+1] +
   *                                psiA[i,j] +psiB[k,l] +
   *                                tau[Ai,Aj,Bk,Bl]
   */
  max_shq = NEGINF;
  /* In this loop, we just start threads (or call synchronously in this
     thread) to compute the values, they are used in another loop
     after this one has done and threads have completed */
  for (x = 0; x < bpaglobals.ipsilistA[i].num_elements; x++)
  {
    h = bpaglobals.ipsilistA[i].ipsi[x].right;
    if (h >= j)
      break;  /* sprted ascending so done in (i,j) interval when past j */
    for (y = 0; y < bpaglobals.ipsilistB[k].num_elements; y++)
    {
      q = bpaglobals.ipsilistB[k].ipsi[y].right;
      if (q >= l)
        break; /* similarly have finished in (k,l) interval */

      /* note there appears to be an error in the paper
       * in this equation; it should be h-1 and q-1 as
       * here not +1 and +1 (j+1 and l+1 in the paper's
       * formulation).
       */
/*
      sm = bpa_dynprogm_thread(i+1, h-1, k+1, q-1) + pairedscore;
      shq = sm + bpa_dynprogm_thread(h+1, j, q+1, l);
*/

/*      fprintf(stderr, "%d\t%d\t%d\t%d\n", i, h, k, q);*/

      bpa_dynprogm_thread_call(mydata->thread_id, &num_active_threads,
                               &active_threadids[actindex],
                               i+1, h-1, k+1, q-1);
      if (active_threadids[actindex] != mydata->thread_id)
      {
        fprintf(stderr, "z %d: %d\n", actindex, active_threadids[actindex]);
        actindex++;
      }
      bpa_dynprogm_thread_call(mydata->thread_id, &num_active_threads,
                               &active_threadids[actindex],
                               h+1, j, q+1, l);
      if (active_threadids[actindex] != mydata->thread_id)
      {
        fprintf(stderr, "z %d: %d\n", actindex, active_threadids[actindex]);
        actindex++;
      }

    }
  }

#ifndef SYNCH
  /* FIXME TODO for now we will wait for all threads here.
     We should change this so that instead of crude starting of threads
     we maintain a pool of threads that signal when done, then we can
     use those values and start them with new work. PRobably not even
     start/end threads, just start them all, and use mutex/condition variables
     to co-ordinate */
  if (is_master)
  {
    /* in the master thread, wait for all other threads to finish */
    for (t = 0; t < actindex; t++)
    {
      fprintf(stderr, "%d / %d\n", active_threadids[t],num_active_threads);
      if ((rc = pthread_join(threads[active_threadids[t]], NULL)))
        bpa_fatal_error(funcname, "pthread_join failed (%d)\n", rc);
      num_active_threads--;
    }
  }
#endif

  /* get values from hashtable. In threads other than master,
     these must be here as calls made synchronously in thread */
  if (comp_gapB)
    gapB = oahttslf_lookup_indices(i + 1, j, k, l) + bpaglobals.gamma;
  if (comp_gapA)
    gapA = oahttslf_lookup_indices(i, j, k + 1, l) + bpaglobals.gamma;
  if (comp_unpaired)
  {
    sigma_ik = BPA_SIGMA(bpaglobals.seqA[i], bpaglobals.seqB[k]);
    unpaired = oahttslf_lookup_indices(i+1, j, k+1, l) + sigma_ik;
  }

/*   assert(gapB > NEGINF); */
/*   assert(gapA > NEGINF); */
/*   assert(unpaired > NEGINF); */

  gapmax = MAX(gapA, gapB);
  score = MAX(gapmax, unpaired); /* max of first 3 cases */

  /* Now get all the sm and shq values from hashtable and find max */
  /* All the values must be present as we have either joined the thread
     that set them, or it was called synchronously in this thread */
  for (x = 0; x < bpaglobals.ipsilistA[i].num_elements; x++)
  {
    h = bpaglobals.ipsilistA[i].ipsi[x].right;
    if (h >= j)
      break;  /* sprted ascending so done in (i,j) interval when past j */
    psiA_ih = bpaglobals.ipsilistA[i].ipsi[x].psi;
    for (y = 0; y < bpaglobals.ipsilistB[k].num_elements; y++)
    {
      q = bpaglobals.ipsilistB[k].ipsi[y].right;
      if (q >= l)
        break; /* similarly have finished in (k,l) interval */
      psiB_kq = bpaglobals.ipsilistB[k].ipsi[y].psi;

/*      fprintf(stderr, "retr %d\t%d\t%d\t%d\n", i, h, k, q); */
      
      pairedscore = psiA_ih + psiB_kq; /* TODO: add sigma_tau() score too */
      assert(pairedscore >= 0);
      sm = oahttslf_lookup_indices( i+1, h-1, k+1, q-1) + pairedscore;
      shq = sm + oahttslf_lookup_indices( h+1, j, q+1, l);
/*      assert(sm > NEGINF); */
      if (shq > max_shq)
        max_shq = shq;
    }
  }

  score = MAX(score, max_shq);

  bpa_log_msg(funcname, "%d\tS\t%d\t%d\t%d\t%d\t%lg\n",mydata->thread_id,i,j,k,l,score);
  oahttslf_insert_indices(i, j, k, l, score, thread_id);
/*  assert(oahttslf_lookup_indices(i,j,k,l) == score); */
#ifdef USE_INSTRUMENT
  bpastats[mydata->thread_id].count_S++;
#endif
  return NULL;
}




/*
 * bpa_dynprogm_thread_master()
 *
 *   Caller interface to the multithreaded version:
 *   just calls the actual implementation after setting up thread
 *   parameter block so that the instance called in this thread
 *   is the "master" thread (thread_id 0).
 *
 *      The dynamic programming (top down memoization)
 *      computation for base pair probability matrix alignment.
 *      See Hofacker et al 2004 (p. 2223).
 *
 *      This version is the memory function version; instead of computing
 *      the whole S array bottom-up, it is computed recursively top-down
 *      and values stored in it, and reused (memoization) if already
 *      computed. 
 *
 *      This version is multi-threaded, sharing hashtable used to
 *      store computed S values between the threads.
 *
 *      This version uses no bounding.
 *
 *
 *      Is is also required that the ipsilist[i] lists are sorted
 *      by right co-ord ascending. This is assured by input processing.
 *
 *
 *      Parameters:   i     - left co-ord in first sequence
 *                    j     - right co-ord in first sequence
 *                    k     - left co-ord in second sequence
 *                    l     - right co-ord in second sequence
 *                            0 <= i < j <= n1 - 1
 *                            0 <= k < l <= n2 - 1
 *
 *      Uses global data:
 *                  read/write:
 *                    hashtable in ht.c
 *                    count_dynprogm_entry - count of calls here
 *                    count_dynprogm_notmemoed - where not return memoed value
 *
 *                  readonly:
 *                    seqA    - first sequence
 *                    seqB    - second sequence
 *                    seqlenA - length of first sequence
 *                    seqlenB - length of second sequence
 *                    ipsilistA - (j,psi) lists indexed by i for 1st seq
 *                    ipsilistB - (j,psi) lists indexed by i for 2nd seq
 *                    gamma   - gap penalty (<= 0)
 *                    printstats -flag to print instrumentation data
 *                    num_threads -max number of threads allowed
 *
 *      Return value: value of d.p. at (i,j,k,l).
 *
 */
double bpa_dynprogm_thread_master(int i, int j, int k, int l)
{
  thread_data_t master_thread_data;

#ifdef USE_INSTRUMENT
  int t;
  counter_t total_count_S = 0, total_count_dynprogm_entry = 0, 
    total_count_dynprogm_entry_notmemoed = 0;
#endif

  master_thread_data.thread_id = MASTER_THREAD_ID;
  master_thread_data.i = i;
  master_thread_data.j = j;
  master_thread_data.k = k;
  master_thread_data.l = l;
  /* run master in this thread */
  (void)bpa_dynprogm_thread(&master_thread_data);

  if (bpaglobals.printstats)
  {
#ifdef USE_INSTRUMENT
    for (t = 0; t < bpaglobals.num_threads; t++)
    {
      printf("stats for thread %d:\n", t);
      printf("  S cells computed = %lu\n", bpastats[t].count_S);
      printf("  calls to dynprogm = %lu\n", bpastats[t].count_dynprogm_entry);
      printf("  calls to dynprogm where not memoed = %lu\n",
             bpastats[t].count_dynprogm_entry_notmemoed);
      total_count_S += bpastats[t].count_S;
      total_count_dynprogm_entry += bpastats[t].count_dynprogm_entry;
      total_count_dynprogm_entry_notmemoed += bpastats[t].count_dynprogm_entry_notmemoed;
    }
    printf("totals:\n");
    printf("  S cells computed = %lu\n", total_count_S);
    printf("  calls to dynprogm = %lu\n", total_count_dynprogm_entry);
    printf("  calls to dynprogm where not memoed = %lu\n",
           total_count_dynprogm_entry_notmemoed);
#else
    printf("COMPILED WITHOUT -DUSE_INSTRUMENT\n");
#endif
  }
  return oahttslf_lookup_indices(i, j, k, l);
}



