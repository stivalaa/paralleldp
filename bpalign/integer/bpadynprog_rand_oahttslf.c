/*****************************************************************************
 *
 * File:    bpadynprog_rand_httslf.c
 * Author:  Alex Stivala
 * Created: April 2009, 
 *          from phd/bpalign/src/bpadynprog.c June 2007
 *
 * $Id: bpadynprog_rand_oahttslf.c 3037 2009-12-14 04:05:16Z alexs $
 *
 * pthreads implementation, with dp recursive call ordering random choice,
 * using the httslf lockfree hashtable of RNA base pair
 * probability matrix alignment by dynamic programming.
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

#define MAX_IPSILIST_LEN 600 /* TODO: make this dynamic */

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

/* mutex and condition variable to signal termination so we can wait on
   any thread to terminate in the master, not have to join specific thread */
static pthread_mutex_t term_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t term_cond = PTHREAD_COND_INITIALIZER;
static int term_thread_id = -1;  /* thread_id of finished thread */

/* the number of active threads. This is only used by the master thread */
static unsigned int num_active_threads = 0; /* we do not count the master thread */ 

/*****************************************************************************
 *
 * static functions
 *
 *****************************************************************************/


/* insert by (i,j,k,l) into table */
static void oahttslf_insert_indices(uint16_t i, uint16_t j, 
                                    uint16_t k, uint16_t l,
                                    myint64_t value, int thread_id);

/* lookup by (i,j,k,l) */
static myint64_t oahttslf_lookup_indices(uint16_t i, uint16_t j,
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
                                    myint64_t value,
                                    int thread_id)

{
  uint64_t key;
  myint64_t val;

  key = (i == 0 && j == 0 && k == 0 && l == 0 ? MAGIC_ZERO : 
         ((uint64_t)i << 47) | ((uint64_t)(j & 0xffff) << 31) |
         ((uint64_t)(k & 0xffff) << 15) | (uint64_t)(l & 0xffff));
//  printf("zz %d,%d,%d,%d\t\t%016lX\n",i,j,k,l,key);
  val = (value == 0 ? NEGINF : value);
  oahttslf_insert(key, val, thread_id);
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
static myint64_t oahttslf_lookup_indices(uint16_t i, uint16_t j,
                                    uint16_t k, uint16_t l)
{
  uint64_t key;
  myint64_t val;
  bool found;
  key = (i == 0 && j == 0 && k == 0 && l == 0 ? MAGIC_ZERO : 
         ((uint64_t)i << 47) | ((uint64_t)(j & 0xffff) << 31) |
         ((uint64_t)(k & 0xffff) << 15) | (uint64_t)(l & 0xffff));

  found =  oahttslf_lookup(key, &val);
  if (found)
    return (val <= NEGINF ? 0 : val);
  else
    return NEGINF;
}





/*****************************************************************************
 *
 * external functions
 *
 *****************************************************************************/


/*------------------------------------------------------------------------*/
/*                        httslf version                                  */


void *bpa_dynprogm_thread_wrapper(void *threadarg);

static myint64_t bpa_dynprogm(int i, int j, int k, int l,
                    int thread_id, unsigned int *seed);


/*
 * bpa_dynprogm_thread_wrapper() - thread interface to bpa_dynprogm()
 *
 *      This version is multi-threaded, sharing hashtable used to
 *      store computed values between the threads. This function
 *      is just the pthreads interface to the recursive bpa_dynprogm()
 *      function itself. We compute the value and then signal
 *      a condition variable to indicate the thread has finished,
 *      so the master thread can wait for ANY thread to terminate,
 *      not having to explicitly join a  particular thread.
 *
 *      Parameters:
 *         threadarg - thread data for this thread
 *
 *      Return value:
 *         pointer to myint64_t containing score computed
 *         (declared void* for pthreads).
 *
 */
void *bpa_dynprogm_thread_wrapper(void *threadarg)
{
  thread_data_t *mydata = (thread_data_t *)threadarg;

  unsigned int seed = (unsigned int)pthread_self() * time(NULL);
  mydata->score = bpa_dynprogm(mydata->i, mydata->j, mydata->k, mydata->l,
                               mydata->thread_id, &seed);

  /* signal thread termination so master can detect a thread finished */
  pthread_mutex_lock(&term_mutex);
  if (term_thread_id == -1) /* master only cares about first thread to exit */
    term_thread_id = mydata->thread_id;
  pthread_cond_signal(&term_cond);
  pthread_mutex_unlock(&term_mutex);
  return &mydata->score;
}


/*
 * bpa_dynprogm()
 *
 *      The dynamic programming (top down memoization)
 *      computation for base pair probability matrix alignment.
 *      See Hofacker et al 2004 (p. 2223).
 *
 *      This version is the memory function version; instead of computing
 *      the whole S array bottom-up, it is computed recursively top-down
 *      and values stored in it, and reused (memoization) if already
 *      computed. A hash table is used to store the computed values,
 *      so only storage for actually computed values is allocated.
 *
 *
 *      Choice of subproblem order is randomized so that all threads
 *      run this same code, but diverge in thei rprocessing randomly.
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
 *                  thread_id - id (0,...n, not pthread id) of this thread
 *                  seed   - seed for rand_r()
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
 *
 *      Return value: The value of the dp at i,j,k,l
 *
 */
static myint64_t bpa_dynprogm(int i, int j, int k, int l,
                    int thread_id, unsigned int *seed)
{
  static const char *funcname = "bpa_dynprogm";
  myint64_t score = NEGINF;
  myint64_t gapA, gapB, unpaired, gapmax, pairedscore;
  myint64_t sigma_ik;
  int h,q; /* h and q are the pairing co-ords used in the recurrence */
  int x,y,xprime,yprime,z; /* just loop indices, no meaning */
  myint64_t psiA_ih, psiB_kq; /* psi value at seqA[i,h] and seqB[k,q] */
  myint64_t sm, shq, max_shq;
  myint64_t value;
  int ipsilistA_permutation[MAX_IPSILIST_LEN];
  int ipsilistB_permutation[MAX_IPSILIST_LEN];
  bool done_gapA = FALSE, done_gapB = FALSE, done_unpaired = FALSE;

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

  bpa_log_msg(funcname, "\t%d\t%d\t%d\t%d\n",i,j,k,l);

#ifdef USE_INSTRUMENT
  bpastats[thread_id].count_dynprogm_entry++;
#endif

  /* memoization: if value here already computed then just return it */
  if ((value =  oahttslf_lookup_indices(i, j, k, l)) > NEGINF)
    return value;

#ifdef USE_INSTRUMENT
  bpastats[thread_id].count_dynprogm_entry_notmemoed++;
#endif

  /*
   *  Initialization cases for the d.p. matrix S:
   *    S(i,j,k,l) = |(j-i) - (l-k)| * gamma, for j - i <= MinLoop + 1 or
   *                                              l - k <= MinLoop + 1
   */
  if ((j - i) <= MINLOOP + 1 || (l - k) <= MINLOOP + 1)
  {
    score = fabs((j - i) - (l - k)) * bpaglobals.gamma;
    bpa_log_msg(funcname, "I\t%d\t%d\t%d\t%d\t%lld\n",i,j,k,l,score);
    oahttslf_insert_indices(i, j, k, l, score, thread_id);
#ifdef USE_INSTRUMENT
    bpastats[thread_id].count_S++;
#endif
    return score;
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
  
  /*
   * max_shq = max{h<=j,q<=l}( S^M[i,h,k,q] + S[h+1,j,q+1,l] )
   *           where S^M[i,j,k,l] = S[i+1,j+1,k+1,l+1] +
   *                                psiA[i,j] +psiB[k,l] +
   *                                tau[Ai,Aj,Bk,Bl]
   */
  max_shq = NEGINF;
  /* iterate over  ipsilistA  elements in random order */
  /* The additional 3 indices stand for gapB, gapA and unpaired
     respectively, so we are randomly ordering not just the ipsilist
     (pair match) subproblems, but also the gap and unpaired subproblems */
  if (bpaglobals.ipsilistA[i].num_elements + 3 > MAX_IPSILIST_LEN)
    bpa_fatal_error(funcname, "increase MAX_IPSILIST_LEN to %d\n",
                    bpaglobals.ipsilistA[i].num_elements + 3);
  
  if (bpaglobals.use_random)
    random_permutation(ipsilistA_permutation,
                       bpaglobals.ipsilistA[i].num_elements + 3, seed);
  else
  {
    ipsilistA_permutation[0] = bpaglobals.ipsilistA[i].num_elements;    
    ipsilistA_permutation[1] = ipsilistA_permutation[0] + 1;
    ipsilistA_permutation[2] = ipsilistA_permutation[1] + 1;
    for (z = 3; z < bpaglobals.ipsilistA[i].num_elements + 3; z++)
      ipsilistA_permutation[z] = z - 3;
  }
  for (x = 0; x < bpaglobals.ipsilistA[i].num_elements + 3; x++)
  {
    xprime = ipsilistA_permutation[x];
    if (xprime >= bpaglobals.ipsilistA[i].num_elements)
    {
      /* one of the two gap cases or the unpaired cases, not an ipsilist case*/
      switch (xprime - bpaglobals.ipsilistA[i].num_elements)
      {
        case 0:
          if (i + 1 < bpaglobals.seqlenA && i + 1 < j)
            gapB = bpa_dynprogm(i + 1, j, k, l, thread_id, seed) + bpaglobals.gamma;
          else
            gapB = NEGINF;
          done_gapB = TRUE;
          break;

        case 1:
          if (k + 1 < bpaglobals.seqlenB && k + 1 < l)
            gapA = bpa_dynprogm(i, j, k + 1, l, thread_id, seed) + bpaglobals.gamma;
          else
            gapA = NEGINF;
          done_gapA = TRUE;
          break;

        case 2:
          if (i+1 < bpaglobals.seqlenA && i+1 < j && k+1 < bpaglobals.seqlenB && k+1 < l)
          {
            sigma_ik = BPA_SIGMA(bpaglobals.seqA[i], bpaglobals.seqB[k]);
            unpaired = bpa_dynprogm(i+1, j, k+1, l, thread_id, seed) + sigma_ik;
          }
          else
            unpaired = NEGINF;
          done_unpaired = TRUE;
          break;

        default:
          bpa_fatal_error(funcname, "impossible case %d\n", 
                          xprime - bpaglobals.ipsilistA[i].num_elements);
          break;
      }
      continue; /* done with this case */
    }

    /* processing  one of the ipsilistA elements */
    h = bpaglobals.ipsilistA[i].ipsi[xprime].right;
    if (h >= j) continue;
//    if (h >= j)
//      break;  /* sprted ascending so done in (i,j) interval when past j */
    psiA_ih = bpaglobals.ipsilistA[i].ipsi[xprime].psi;
    /* iterate over ipsilistB elements in random order */
    if (bpaglobals.ipsilistB[k].num_elements > MAX_IPSILIST_LEN)
      bpa_fatal_error(funcname, "increase MAX_IPSILIST_LEN to %d\n",
                      bpaglobals.ipsilistA[k].num_elements);
    if (bpaglobals.use_random)
      random_permutation(ipsilistB_permutation,
                         bpaglobals.ipsilistB[k].num_elements, seed);
    else
      for (z = 0; z < bpaglobals.ipsilistB[k].num_elements; z++)
        ipsilistB_permutation[z] = z;
    for (y = 0; y < bpaglobals.ipsilistB[k].num_elements; y++)
    {
      yprime = ipsilistB_permutation[y];
      q = bpaglobals.ipsilistB[k].ipsi[yprime].right;
//      if (q >= l)
//        break; /* similarly have finished in (k,l) interval */
      if (q >= l) continue;
      psiB_kq = bpaglobals.ipsilistB[k].ipsi[yprime].psi;
      
      pairedscore = psiA_ih + psiB_kq; /* TODO: add sigma_tau() score too */
      
      assert(pairedscore >= 0);
      
      /* note there appears to be an error in the paper
       * in this equation; it should be h-1 and q-1 as
       * here not +1 and +1 (j+1 and l+1 in the paper's
       * formulation).
       */
      sm = bpa_dynprogm(i+1, h-1, k+1, q-1, thread_id, seed) + pairedscore;
      shq = sm + bpa_dynprogm(h+1, j, q+1, l, thread_id, seed);
      if (shq > max_shq)
        max_shq = shq;
    }
  }
  
  assert(done_gapA && done_gapB && done_unpaired);

  gapmax = MAX(gapA, gapB);
  score = MAX(gapmax, unpaired); /* max of first 3 cases */
  score = MAX(score, max_shq);

  bpa_log_msg(funcname, "S\t%d\t%d\t%d\t%d\t%lld\n",i,j,k,l,score);
  oahttslf_insert_indices(i, j, k, l, score, thread_id);
#ifdef USE_INSTRUMENT
  bpastats[thread_id].count_S++;
#endif
  return score;
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
 *                    printstats -flag to compute instrumentation data
 *                    verbose -flag to print instrumentation data
 *                    num_threads -max number of threads allowed
 *
 *      Return value: value of d.p. at (i,j,k,l).
 *
 */
int64_t bpa_dynprogm_thread_master(int i, int j, int k, int l)
{
  static const char *funcname = "dp_dynprogm_thread_master";

  int new_thread_id;
  int finished_thread_id;
  int rc;
  myint64_t *score;
#ifdef USE_INSTRUMENT
  int t;
  extern  counter_t total_count_S , total_count_dynprogm_entry , 
    total_count_dynprogm_entry_notmemoed, num_keys ;
#endif
  int q;
  
  term_thread_id = -1;
  num_active_threads = 0; /* do not count master thread */

  while (num_active_threads < (unsigned)bpaglobals.num_threads)
  {
#ifdef DEBUG
    fprintf(stderr, "starting thread id %d\n", num_active_threads);
#endif
    new_thread_id = num_active_threads;
    thread_data[num_active_threads].thread_id = num_active_threads;
    thread_data[num_active_threads].i = i;
    thread_data[num_active_threads].j = j;
    thread_data[num_active_threads].k = k;
    thread_data[num_active_threads].l = l;
    if ((rc = pthread_create(&threads[num_active_threads], NULL,
                             bpa_dynprogm_thread_wrapper,
                             (void *)&thread_data[num_active_threads])))
      bpa_fatal_error(funcname, "pthread_create() failed (%d)\n", rc);
    num_active_threads++;
  }
  
  
  /* in the master thread, wait for first thread to finish (don't care
     about rest) */
  pthread_mutex_lock(&term_mutex);
  /* only wait if no thread has finsished: in POSIX threads, singalling
     a condition before waiting (which happens here if worker thread finishes
     extremely quickly) is a l5Bogic error (results in master waiting
     forever on my system) */
  if (term_thread_id == -1)
    pthread_cond_wait(&term_cond, &term_mutex); /*unlocks mutex while waiting*/
  finished_thread_id = term_thread_id;
  pthread_mutex_unlock(&term_mutex);
#ifdef DEBUG
  fprintf(stderr, "thread id %d finished\n", finished_thread_id);
#endif
  if ((rc = pthread_join(threads[finished_thread_id],&score)))
      bpa_fatal_error(funcname, "pthread_join failed (%d)\n", rc);
  num_active_threads--;

  /* cleanup threads (some will still be running) */
  for (q = 0 ; q < bpaglobals.num_threads; q++)
  {
    if (q == finished_thread_id)
      continue;
    if ((rc = pthread_cancel(threads[q])))
      /* bpa_error_msg(funcname, "pthread_cancel failed (%d)\n", rc); */ /* safe to ignore this error --- thread may have finished by now */ ;
  }
  for (q = 0 ; q < bpaglobals.num_threads; q++)
  {
    if (q == finished_thread_id)
      continue;
    if ((rc = pthread_join(threads[q], NULL)))
      bpa_fatal_error(funcname, "pthread_join [2] failed (%d)\n", rc);
    num_active_threads--;
  }
  
  if (num_active_threads != 0)
    bpa_fatal_error(funcname, "%d active threads at termination\n", num_active_threads);

  if (bpaglobals.printstats)
  {
#ifdef USE_INSTRUMENT
    for (t = 0; t < bpaglobals.num_threads; t++)
    {
      if (bpaglobals.verbose)
      {
        printf("stats for thread %d:\n", t);
        printf("  S cells computed = %lu\n", bpastats[t].count_S);
        printf("  calls to dynprogm = %lu\n", bpastats[t].count_dynprogm_entry);
        printf("  calls to dynprogm where not memoed = %lu\n",
               bpastats[t].count_dynprogm_entry_notmemoed);
      }
      total_count_S += bpastats[t].count_S;
      total_count_dynprogm_entry += bpastats[t].count_dynprogm_entry;
      total_count_dynprogm_entry_notmemoed += bpastats[t].count_dynprogm_entry_notmemoed;
    }
    num_keys = (counter_t)oahttslf_total_key_count();
    if (bpaglobals.verbose) 
    {
      printf("totals:\n");
      printf("  S cells computed = %lu\n", total_count_S);
      printf("  calls to dynprogm = %lu\n", total_count_dynprogm_entry);
      printf("  calls to dynprogm where not memoed = %lu\n",
             total_count_dynprogm_entry_notmemoed);
    }
#else
    if (bpaglobals.verbose)
      printf("COMPILED WITHOUT -DUSE_INSTRUMENT\n");
#endif
  }

  return *score;
}



/*------------------------------------------------------------------------*/
/*                         array version                                  */


void *bpa_dynprogm_thread_array_wrapper(void *threadarg);

static myint64_t bpa_dynprogm_thread_array(int i, int j, int k, int l,
                                        myint64_t *S,
                                        int thread_id, unsigned int *seed);


/*
 * bpa_dynprogm_thread_array_wrapper() - thread interface to bpa_dynprogm_array()
 *
 *      This version is multi-threaded, sharing array used to
 *      store computed values between the threads. This function
 *      is just the pthreads interface to the recursive
 *      bpa_dynprogm_thread_array()
 *      function itself. We compute the value and then signal
 *      a condition variable to indicate the thread has finished,
 *      so the master thread can wait for ANY thread to terminate,
 *      not having to explicitly join a  particular thread.
 *
 *      Parameters:
 *         threadarg - thread data for this thread
 *
 *      Return value:
 *         pointer to myint64_t containing score computed
 *         (declared void* for pthreads).
 *
 */
void *bpa_dynprogm_thread_array_wrapper(void *threadarg)
{
  thread_data_t *mydata = (thread_data_t *)threadarg;

  unsigned int seed = (unsigned int)pthread_self() * time(NULL);
  mydata->score = bpa_dynprogm_thread_array(
    mydata->i, mydata->j, mydata->k, mydata->l, mydata->S,
    mydata->thread_id, &seed);

  /* signal thread termination so master can detect a thread finished */
  pthread_mutex_lock(&term_mutex);
  if (term_thread_id == -1) /* master only cares about first thread to exit */
    term_thread_id = mydata->thread_id;
  pthread_cond_signal(&term_cond);
  pthread_mutex_unlock(&term_mutex);
  return &mydata->score;
}


/*
 * bpa_dynprogm_thread_array()
 *
 *      The dynamic programming (top down memoization)
 *      computation for base pair probability matrix alignment.
 *      See Hofacker et al 2004 (p. 2223).
 *
 *      This version is the memory function version; instead of computing
 *      the whole S array bottom-up, it is computed recursively top-down
 *      and values stored in it, and reused (memoization) if already
 *      computed. An array is used to store the computed values.
 *
 *
 *      Choice of subproblem order is randomized so that all threads
 *      run this same code, but diverge in thei rprocessing randomly.
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
 *                    S     - the 4d d.p. array
 *                  thread_id - id (0,...n, not pthread id) of this thread
 *                  seed   - seed for rand_r()
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
 *
 *      Return value: The value of the dp at i,j,k,l
 *
 */
static myint64_t bpa_dynprogm_thread_array(int i, int j, int k, int l,
                                        myint64_t *S,
                                        int thread_id, unsigned int *seed)
{
  static const char *funcname = "bpa_dynprogm_thread_array";
  myint64_t score = NEGINF;
  myint64_t gapA, gapB, unpaired, gapmax, pairedscore;
  myint64_t sigma_ik;
  int h,q; /* h and q are the pairing co-ords used in the recurrence */
  int x,y,xprime,yprime,z; /* just loop indices, no meaning */
  myint64_t psiA_ih, psiB_kq; /* psi value at seqA[i,h] and seqB[k,q] */
  myint64_t sm, shq, max_shq;
  int n1 = bpaglobals.seqlenA;
  int n2 = bpaglobals.seqlenB;
  int ipsilistA_permutation[MAX_IPSILIST_LEN];
  int ipsilistB_permutation[MAX_IPSILIST_LEN];
  bool done_gapA = FALSE, done_gapB = FALSE, done_unpaired = FALSE;


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

  bpa_log_msg(funcname, "\t%d\t%d\t%d\t%d\n",i,j,k,l);

#ifdef USE_INSTRUMENT
  bpastats[thread_id].count_dynprogm_entry++;
#endif

  /* memoization: if value here already computed then just return it */
  if ((score = S[INDEX4D(i,j,k,l,n1,n2)]) != NEGINF)
    return score;

#ifdef USE_INSTRUMENT
  bpastats[thread_id].count_dynprogm_entry_notmemoed++;
#endif

  /*
   *  Initialization cases for the d.p. matrix S:
   *    S(i,j,k,l) = |(j-i) - (l-k)| * gamma, for j - i <= MinLoop + 1 or
   *                                              l - k <= MinLoop + 1
   */
  if ((j - i) <= MINLOOP + 1 || (l - k) <= MINLOOP + 1)
  {
    score = fabs((j - i) - (l - k)) * bpaglobals.gamma;
    bpa_log_msg(funcname, "I\t%d\t%d\t%d\t%d\t%lld\n",i,j,k,l,score);
    S[INDEX4D(i,j,k,l,n1,n2)] = score;
#ifdef USE_INSTRUMENT
    bpastats[thread_id].count_S++;
#endif
    return score;
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
  
  /*
   * max_shq = max{h<=j,q<=l}( S^M[i,h,k,q] + S[h+1,j,q+1,l] )
   *           where S^M[i,j,k,l] = S[i+1,j+1,k+1,l+1] +
   *                                psiA[i,j] +psiB[k,l] +
   *                                tau[Ai,Aj,Bk,Bl]
   */
  max_shq = NEGINF;
  /* iterate over  ipsilistA  elements in random order */
  /* The additional 3 indices stand for gapB, gapA and unpaired
     respectively, so we are randomly ordering not just the ipsilist
     (pair match) subproblems, but also the gap and unpaired subproblems */
  if (bpaglobals.ipsilistA[i].num_elements + 3 > MAX_IPSILIST_LEN)
    bpa_fatal_error(funcname, "increase MAX_IPSILIST_LEN to %d\n",
                    bpaglobals.ipsilistA[i].num_elements + 3);
  if (bpaglobals.use_random)
    random_permutation(ipsilistA_permutation,
                       bpaglobals.ipsilistA[i].num_elements + 3, seed);
  else
  {
    ipsilistA_permutation[0] = bpaglobals.ipsilistA[i].num_elements;    
    ipsilistA_permutation[1] = ipsilistA_permutation[0] + 1;
    ipsilistA_permutation[2] = ipsilistA_permutation[1] + 1;
    for (z = 3; z < bpaglobals.ipsilistA[i].num_elements + 3; z++)
      ipsilistA_permutation[z] = z - 3;
  }
  for (x = 0; x < bpaglobals.ipsilistA[i].num_elements + 3; x++)
  {
    xprime = ipsilistA_permutation[x];
    if (xprime >= bpaglobals.ipsilistA[i].num_elements)
    {
      /* one of the two gap cases or the unpaired cases, not an ipsilist case*/
      switch (xprime - bpaglobals.ipsilistA[i].num_elements)
      {
        case 0:
          if (i + 1 < bpaglobals.seqlenA && i + 1 < j)
            gapB = bpa_dynprogm_thread_array(i + 1, j, k, l, S, thread_id, seed) + bpaglobals.gamma;
          else
            gapB = NEGINF;
          done_gapB = TRUE;
          break;

        case 1:
          if (k + 1 < bpaglobals.seqlenB && k + 1 < l)
            gapA = bpa_dynprogm_thread_array(i, j, k + 1, l, S, thread_id, seed) + bpaglobals.gamma;
          else
            gapA = NEGINF;
          done_gapA = TRUE;
          break;

        case 2:
          if (i+1 < bpaglobals.seqlenA && i+1 < j && k+1 < bpaglobals.seqlenB && k+1 < l)
          {
            sigma_ik = BPA_SIGMA(bpaglobals.seqA[i], bpaglobals.seqB[k]);
            unpaired = bpa_dynprogm_thread_array(i+1, j, k+1, l, S, thread_id, seed) + sigma_ik;
          }
          else
            unpaired = NEGINF;
          done_unpaired = TRUE;
          break;

        default:
          bpa_fatal_error(funcname, "impossible case %d\n", 
                          xprime - bpaglobals.ipsilistA[i].num_elements);
          break;
      }
      continue; /* done with this case */
    }

    /* processing  one of the ipsilistA elements */
    h = bpaglobals.ipsilistA[i].ipsi[xprime].right;
    if (h >= j) continue;
//    if (h >= j)
//      break;  /* sprted ascending so done in (i,j) interval when past j */
    psiA_ih = bpaglobals.ipsilistA[i].ipsi[xprime].psi;
    /* iterate over ipsilistB elements in random order */
    if (bpaglobals.ipsilistB[k].num_elements > MAX_IPSILIST_LEN)
      bpa_fatal_error(funcname, "increase MAX_IPSILIST_LEN to %d\n",
                      bpaglobals.ipsilistB[k].num_elements);
    if (bpaglobals.use_random)
      random_permutation(ipsilistB_permutation,
                         bpaglobals.ipsilistB[k].num_elements, seed);
    else
      for (z = 0; z < bpaglobals.ipsilistB[k].num_elements; z++)
        ipsilistB_permutation[z] = z;
    for (y = 0; y < bpaglobals.ipsilistB[k].num_elements; y++)
    {
      yprime = ipsilistB_permutation[y];
      q = bpaglobals.ipsilistB[k].ipsi[yprime].right;
//      if (q >= l)
//        break; /* similarly have finished in (k,l) interval */
      if (q >= l) continue;
      psiB_kq = bpaglobals.ipsilistB[k].ipsi[yprime].psi;
      
      pairedscore = psiA_ih + psiB_kq; /* TODO: add sigma_tau() score too */
      
      assert(pairedscore >= 0);
      
      /* note there appears to be an error in the paper
       * in this equation; it should be h-1 and q-1 as
       * here not +1 and +1 (j+1 and l+1 in the paper's
       * formulation).
       */
      sm = bpa_dynprogm_thread_array(i+1, h-1, k+1, q-1, S, thread_id, seed) + pairedscore;
      shq = sm + bpa_dynprogm_thread_array(h+1, j, q+1, l, S, thread_id, seed);
      if (shq > max_shq)
        max_shq = shq;
    }
  }
  
  assert(done_gapA && done_gapB && done_unpaired);

  gapmax = MAX(gapA, gapB);
  score = MAX(gapmax, unpaired); /* max of first 3 cases */
  score = MAX(score, max_shq);

  bpa_log_msg(funcname, "S\t%d\t%d\t%d\t%d\t%lld\n",i,j,k,l,score);
  S[INDEX4D(i,j,k,l,n1,n2)] = score;
#ifdef USE_INSTRUMENT
  bpastats[thread_id].count_S++;
#endif
  return score;
}



/*
 * bpa_dynprogm_thread_array_master()
 *
 *   Caller interface to the multithreaded array version:
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
 *      This version is multi-threaded, sharing array used to
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
 *                    S     - the 4d d.p. array
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
 *                    printstats -flag to compute instrumentation data
 *                    verbose -flag to print instrumentation data
 *                    num_threads -max number of threads allowed
 *
 *      Return value: value of d.p. at (i,j,k,l).
 *
 */
int64_t bpa_dynprogm_thread_array_master(int i, int j, int k, int l, myint64_t *S)
{
  static const char *funcname = "dp_dynprogm_thread_array_master";

  int new_thread_id;
  int finished_thread_id;
  int rc;
  myint64_t *score;
  int t;
  extern counter_t total_count_S , total_count_dynprogm_entry , 
    total_count_dynprogm_entry_notmemoed;
  int q;
  
  term_thread_id = -1;
  num_active_threads = 0; /* do not count master thread */

  while (num_active_threads < (unsigned)bpaglobals.num_threads)
  {
#ifdef DEBUG
    fprintf(stderr, "starting thread id %d\n", num_active_threads);
#endif
    new_thread_id = num_active_threads;
    thread_data[num_active_threads].thread_id = num_active_threads;
    thread_data[num_active_threads].i = i;
    thread_data[num_active_threads].j = j;
    thread_data[num_active_threads].k = k;
    thread_data[num_active_threads].l = l;
    thread_data[num_active_threads].S = S;
    if ((rc = pthread_create(&threads[num_active_threads], NULL,
                             bpa_dynprogm_thread_array_wrapper,
                             (void *)&thread_data[num_active_threads])))
      bpa_fatal_error(funcname, "pthread_create() failed (%d)\n", rc);
    num_active_threads++;
  }
  
  
  /* in the master thread, wait for first thread to finish (don't care
     about rest) */
  pthread_mutex_lock(&term_mutex);
  /* only wait if no thread has finsished: in POSIAX threads, singalling
     a condition before waiting (which happens here if worker thread finishes
     extremely quickly) is a l5Bogic error (results in master waiting
     forever on my system) */
  if (term_thread_id == -1)
    pthread_cond_wait(&term_cond, &term_mutex); /*unlocks mutex while waiting*/
  finished_thread_id = term_thread_id;
  pthread_mutex_unlock(&term_mutex);
#ifdef DEBUG
  fprintf(stderr, "thread id %d finished\n", finished_thread_id);
#endif
  if ((rc = pthread_join(threads[finished_thread_id], &score)))
      bpa_fatal_error(funcname, "pthread_join failed (%d)\n", rc);
  num_active_threads--;

  /* cleanup threads (some will still be running) */
  for (q = 0 ; q < bpaglobals.num_threads; q++)
  {
    if (q == finished_thread_id)
      continue;
    if ((rc = pthread_cancel(threads[q])))
      /* bpa_error_msg(funcname, "pthread_cancel failed (%d)\n", rc); */ /* safe to ignore this error --- thread may have finished by now */ ;
  }
  for (q = 0 ; q < bpaglobals.num_threads; q++)
  {
    if (q == finished_thread_id)
      continue;
    if ((rc = pthread_join(threads[q], NULL)))
      bpa_fatal_error(funcname, "pthread_join [2] failed (%d)\n", rc);
    num_active_threads--;
  }
  
  if (num_active_threads != 0)
    bpa_fatal_error(funcname, "%d active threads at termination\n", num_active_threads);

  if (bpaglobals.printstats)
  {
#ifdef USE_INSTRUMENT
    for (t = 0; t < bpaglobals.num_threads; t++)
    {
      if (bpaglobals.verbose)
      {
        printf("stats for thread %d:\n", t);
        printf("  S cells computed = %lu\n", bpastats[t].count_S);
        printf("  calls to dynprogm = %lu\n", bpastats[t].count_dynprogm_entry);
        printf("  calls to dynprogm where not memoed = %lu\n",
               bpastats[t].count_dynprogm_entry_notmemoed);
      }
      total_count_S += bpastats[t].count_S;
      total_count_dynprogm_entry += bpastats[t].count_dynprogm_entry;
      total_count_dynprogm_entry_notmemoed += bpastats[t].count_dynprogm_entry_notmemoed;
    }
    if (bpaglobals.verbose)
    {
      printf("totals:\n");
      printf("  S cells computed = %lu\n", total_count_S);
      printf("  calls to dynprogm = %lu\n", total_count_dynprogm_entry);
      printf("  calls to dynprogm where not memoed = %lu\n",
             total_count_dynprogm_entry_notmemoed);
    }
#else
  if (bpaglobals.verbose)
    printf("COMPILED WITHOUT -DUSE_INSTRUMENT\n");
#endif
  }
  return *score;
}



