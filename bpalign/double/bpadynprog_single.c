/*****************************************************************************
 *
 * File:    bpadynprog_single.c
 * Author:  Alex Stivala
 * Created: April 2009, 
 *          from phd/bpalign/src/bpadynprog.c June 2007
 *
 * $Id: bpadynprog_single.c 2584 2009-06-25 07:54:41Z astivala $
 *
 * CPU (single threaded)  implementations of RNA base pair
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

#include "ht.h"
#include "bpacommon.h"
#include "bpaglobals.h"
#include "bpastats.h"
#include "bpautils.h"
#include "bpaipsilist.h"
#include "bpadynprog_cpu.h"


/*****************************************************************************
 *
 * static data
 *
 *****************************************************************************/

#ifdef USE_INSTRUMENT
/* instrumentation is per-thread, each thread only writes to its own element */
static volatile bpastats_t bpastats[MAX_NUM_THREADS];
#endif

/* insert by (i,j,k,l) into table */
static void ht_insert_indices(unsigned short i, unsigned short j, 
                        unsigned short k, unsigned short l, double value);

/* lookup by (i,j,k,l) */
static ht_entry_t *ht_lookup_indices(unsigned short i, unsigned short j,
                              unsigned short k, unsigned short l);



/*
 * ht_insert_indices()
 *
 * Insert value for (i,j,k,l) into the hashtable
 *
 * Parameters:
 *    i,j,k,l - indices to build insertion key
 *    value - value to insert for the key
 *
 * Return value:
 *    None.
 */
static void ht_insert_indices(unsigned short i, unsigned short j, 
                       unsigned short k, unsigned short l, double value)
{
  tuple4_t key;

  key.i = i;
  key.j = j;
  key.k = k;
  key.l = l;
  ht_insert(&key, &value);
}



/*
 * ht_lookup_indices()
 *
 * Get the value for (i,j,k,l) from the hashtable
 *
 * Parameters:
 *     i,j,k,l - indices to build key for lookup
 * 
 * Return value:
 *     pointer to entry with the key, NULL if not found.
 */
static ht_entry_t *ht_lookup_indices(unsigned short i, unsigned short j,
                              unsigned short k, unsigned short l)
{
  tuple4_t key;

  key.i = i;
  key.j = j;
  key.k = k;
  key.l = l;
  return ht_lookup(&key);
}


/*****************************************************************************
 *
 * external functions
 *
 *****************************************************************************/


/*
 *  The dynamic programming (bottom up) array computation for base
 *  pair probability matrix alignment. See Hofacker et al 2004
 *  (p. 2223).  This is the trivial single-threaded implementation
 *  that simply computes the entire dp matrix (i.e. "bottom-up",
 *  compute all cells, rather than "top-down" or call-based
 *  (recursive)) serially just as a conventional host cpu
 *  implementation, in a single thread.
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
 *                gappenalty - gap penalty (<= 0)
 *                M     - minimum size of hairpin loop 
 *                S     - (workarea) the dp matrix S (n1*n1*n2*n2)
 *                score - (output) final score computed
 *
 */
void dynprog_cpu(int n1, int n2, 
                                const char *seqA, const char *seqB, 
                                const ipsi_element_t *seripsiA, int ld_seripsiA,
                                const ipsi_element_t *seripsiB, int ld_seripsiB,
                                double gappenalty, int M, double *S, 
                                double *score)
{
  int i,j,k,l,x,y;
  double gapA, gapB, unpaired;
  int h,q; /* h and q are the pairing co-ords used in the recurrence */
  double sm, shq, max_shq;
  double psiA_ih, psiB_kq; /* psi value at seqA[i,h] and seqB[k,q] */
  double gapmax,matchmax;

/*  bpa_dump_seripsilist(seripsiA, n1, ld_seripsiA);   */

  /*
   *  Initialization cases for the d.p. matrix S:
   *    S(i,j,k,l) = |(j-i) - (l-k)| * gappenalty, for j - i <= MinLoop + 1 or
   *                                              l - k <= MinLoop + 1
   */
  for (i = 0; i < n1; i++)
    for (j = 0; j < n1; j++)
      for (k = 0; k < n2; k++)
        for (l = 0; l < n2; l++)
          if ((j - i) <= M + 1 || (l - k) <= M + 1)
            S[INDEX4D(i,j,k,l,n1,n2)] =  INTEGER_ABS((j - i) - (l - k)) * gappenalty;
          else
            S[INDEX4D(i,j,k,l,n1,n2)] = 0;
  
  for (j = 0; j < n1; j++)
      for (i = j-1; i >= 0; i--)
          for (l = 0; l < n2; l++)
              for (k = l-1; k >=0; k--)
              {
                /*
                 *     Compute each of the four cases over which we
                 *     choose the max:
                 *      1. gap in second sequence (B): S(i+1,j,k,l) + gappenalty
                 *      2. gap in first sequence (A):  S(i,j+1,k,l) + gappenalty
                 *      3. extension of both subsequences with unpaired position
                 *         S(i+1,j,k+1,l) + sigma(Ai, Bk)
                 *      4. (more complex case with another max):
                 *           match of pair (i, h) in A with (k,q) in B
                 */
                gapA = S[INDEX4D(i+1,j,k,l,n1,n2)] + gappenalty;
                gapB = S[INDEX4D(i,j,k+1,l,n1,n2)] + gappenalty;
                unpaired = S[INDEX4D(i+1,j,k+1,l,n1,n2)] + 
                             (double) BPA_SIGMA(seqA[i], seqB[k]);

                /*
                 * max_shq = max{h<=j,q<=l}( S^M[i,h,k,q] + S[h+1,j,q+1,l] )
                 *           where S^M[i,j,k,l] = S[i+1,j+1,k+1,l+1] +
                 *                                psiA[i,j] +psiB[k,l] +
                 *                                tau[Ai,Aj,Bk,Bl]
                 */
                max_shq = -99999; /* dodgy FIXME */
                for (x = 0; x < ld_seripsiA; x++)
                {
                  h = seripsiA[i * ld_seripsiA + x].right;
                  if (h <= 0)
                    break; /* 0 marks empty entry */
                  if (h >= j)
                    break;  /* sorted so done in (i,j) interval when past j */

                  psiA_ih = seripsiA[i * ld_seripsiA + x].psi;
                  for (y = 0; y < ld_seripsiB; y++)
                  {
                    q = seripsiB[k * ld_seripsiB + y].right;
                    if (q <= 0)
                        break; /* 0 marks empty entry */
                    if (q >= l)
                      break; /* similarly have finished in (k,l) interval */

                    psiB_kq = seripsiB[k * ld_seripsiB + y].psi;

                    sm = S[INDEX4D(i+1,h-1,k+1,q-1,n1,n2)] +
                            psiA_ih + psiB_kq;
                                                          /* TODO add tau */
                    shq = sm + S[INDEX4D(h+1,j,q+1,l,n1,n2)];

                    max_shq = MAX(shq, max_shq);
                  }
                }
                gapmax = MAX(gapA, gapB);
                matchmax = MAX(unpaired, max_shq);
                S[INDEX4D(i,j,k,l,n1,n2)] = MAX(gapmax, matchmax);
              }

  
  *score = S[INDEX4D(0,n1-1,0,n2-1,n1,n2)];
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
 *
 *      Return value: The value of the dp at i,j,k,l
 *
 */
double bpa_dynprogm(int i, int j, int k, int l)
{
  static const char *funcname = "bpa_dynprogm";
  double score = NEGINF;
  double gapA, gapB, unpaired, gapmax, pairedscore;
  double sigma_ik;
  int h,q; /* h and q are the pairing co-ords used in the recurrence */
  int x,y; /* just loop indices, no meaning */
  double psiA_ih, psiB_kq; /* psi value at seqA[i,h] and seqB[k,q] */
  double sm, shq, max_shq;
  double *value;


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
  bpastats[0].count_dynprogm_entry++;
#endif

  /* memoization: if value here already computed then just return it */
  if ((value =  (double *)ht_lookup_indices(i, j, k, l)) != NULL)
    return *value;

#ifdef USE_INSTRUMENT
  bpastats[0].count_dynprogm_entry_notmemoed++;
#endif

  /*
   *  Initialization cases for the d.p. matrix S:
   *    S(i,j,k,l) = |(j-i) - (l-k)| * gamma, for j - i <= MinLoop + 1 or
   *                                              l - k <= MinLoop + 1
   */
  if ((j - i) <= MINLOOP + 1 || (l - k) <= MINLOOP + 1)
  {
    score = fabs((j - i) - (l - k)) * bpaglobals.gamma;
    bpa_log_msg(funcname, "I\t%d\t%d\t%d\t%d\t%lg\n",i,j,k,l,score);
    ht_insert_indices(i, j, k, l, score);
#ifdef USE_INSTRUMENT
    bpastats[0].count_S++;
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

  if (i + 1 < bpaglobals.seqlenA && i + 1 < j)
    gapB = bpa_dynprogm(i + 1, j, k, l) + bpaglobals.gamma;
  else
    gapB = NEGINF;

  if (k + 1 < bpaglobals.seqlenB && k + 1 < l)
    gapA = bpa_dynprogm(i, j, k + 1, l) + bpaglobals.gamma;
  else
    gapA = NEGINF;

  if (i+1 < bpaglobals.seqlenA && i+1 < j && k+1 < bpaglobals.seqlenB && k+1 < l)
  {
    sigma_ik = BPA_SIGMA(bpaglobals.seqA[i], bpaglobals.seqB[k]);
    unpaired = bpa_dynprogm(i+1, j, k+1, l) + sigma_ik;
  }
  else
    unpaired = NEGINF;

  gapmax = MAX(gapA, gapB);
  score = MAX(gapmax, unpaired); /* max of first 3 cases */


  /*
   * max_shq = max{h<=j,q<=l}( S^M[i,h,k,q] + S[h+1,j,q+1,l] )
   *           where S^M[i,j,k,l] = S[i+1,j+1,k+1,l+1] +
   *                                psiA[i,j] +psiB[k,l] +
   *                                tau[Ai,Aj,Bk,Bl]
   */
  max_shq = NEGINF;
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

      pairedscore = psiA_ih + psiB_kq; /* TODO: add sigma_tau() score too */

      assert(pairedscore >= 0);

      /* note there appears to be an error in the paper
       * in this equation; it should be h-1 and q-1 as
       * here not +1 and +1 (j+1 and l+1 in the paper's
       * formulation).
       */
      sm = bpa_dynprogm(i+1, h-1, k+1, q-1) + pairedscore;
      shq = sm + bpa_dynprogm(h+1, j, q+1, l);
      if (shq > max_shq)
        max_shq = shq;
    }
  }
  score = MAX(score, max_shq);

  bpa_log_msg(funcname, "S\t%d\t%d\t%d\t%d\t%lg\n",i,j,k,l,score);
  ht_insert_indices(i, j, k, l, score);
#ifdef USE_INSTRUMENT
  bpastats[0].count_S++;
#endif
  return score;
}



/*
 * bpa_dynprogm_array()
 *
 *      The dynamic programming (top down memoization)
 *      computation for base pair probability matrix alignment.
 *      See Hofacker et al 2004 (p. 2223).
 *
 *      This version is the memory function version; instead of computing
 *      the whole S array bottom-up, it is computed recursively top-down
 *      and values stored in it, and reused (memoization) if already
 *      computed. This version uses an array (4d array but allocated
 *      as linear block of memory so no pointers, just array equation
 *      used to access data cells) so memory allocatd for every possible
 *      cell (even though many won't be used) but no overhead of hashing
 *      etc. as for hash table implementation.
 *
 *      This version uses no bounding.
 *
 *
 *      Is is also required that the ipsilist[i] lists are sorted
 *      by right5A co-ord ascending. This is assured by input processing.
 *
 *
 *      Parameters:   i     - left co-ord in first sequence
 *                    j     - right co-ord in first sequence
 *                    k     - left co-ord in second sequence
 *                    l     - right co-ord in second sequence
 *                            0 <= i < j <= n1 - 1
 *                            0 <= k < l <= n2 - 1
 *                S     - (workarea) the dp matrix S (n1*n1*n2*n2)
 *
 *      Uses global data:
 *                  read/write:
 *                    hasht5Aable in ht.c
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
double bpa_dynprogm_array(int i, int j, int k, int l, double *S)
{
  static const char *funcname = "bpa_dynprogm_array";
  double score = NEGINF;
  double gapA, gapB, unpaired, gapmax, pairedscore;
  double sigma_ik;
  int h,q; /* h and q are the pairing co-ords used in the recurrence */
  int x,y; /* just loop indices, no meaning */
  double psiA_ih, psiB_kq; /* psi value at seqA[i,h] and seqB[k,q] */
  double sm, shq, max_shq;
  int n1 = bpaglobals.seqlenA;
  int n2 = bpaglobals.seqlenB;
  


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
  bpastats[0].count_dynprogm_entry++;
#endif

  /* memoization: if value here already computed then just return it */
  score = S[INDEX4D(i,j,k,l,n1,n2)];
  if (score != NEGINF)
    return score;

#ifdef USE_INSTRUMENT
  bpastats[0].count_dynprogm_entry_notmemoed++;
#endif

  /*
   *  Initialization cases for the d.p. matrix S:
   *    S(i,j,k,l) = |(j-i) - (l-k)| * gamma, for j - i <= MinLoop + 1 or
   *                                              l - k <= MinLoop + 1
   */
  if ((j - i) <= MINLOOP + 1 || (l - k) <= MINLOOP + 1)
  {
    score = fabs((j - i) - (l - k)) * bpaglobals.gamma;
    bpa_log_msg(funcname, "I\t%d\t%d\t%d\t%d\t%lg\n",i,j,k,l,score);
    S[INDEX4D(i,j,k,l,n1,n2)] = score;
#ifdef USE_INSTRUMENT
    bpastats[0].count_S++;
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

  if (i + 1 < bpaglobals.seqlenA && i + 1 < j)
    gapB = bpa_dynprogm_array(i + 1, j, k, l, S) + bpaglobals.gamma;
  else
    gapB = NEGINF;

  if (k + 1 < bpaglobals.seqlenB && k + 1 < l)
    gapA = bpa_dynprogm_array(i, j, k + 1, l, S) + bpaglobals.gamma;
  else
    gapA = NEGINF;

  if (i+1 < bpaglobals.seqlenA && i+1 < j && k+1 < bpaglobals.seqlenB && k+1 < l)
  {
    sigma_ik = BPA_SIGMA(bpaglobals.seqA[i], bpaglobals.seqB[k]);
    unpaired = bpa_dynprogm_array(i+1, j, k+1, l, S) + sigma_ik;
  }
  else
    unpaired = NEGINF;

  gapmax = MAX(gapA, gapB);
  score = MAX(gapmax, unpaired); /* max of first 3 cases */


  /*
   * max_shq = max{h<=j,q<=l}( S^M[i,h,k,q] + S[h+1,j,q+1,l] )
   *           where S^M[i,j,k,l] = S[i+1,j+1,k+1,l+1] +
   *                                psiA[i,j] +psiB[k,l] +
   *                                tau[Ai,Aj,Bk,Bl]
   */
  max_shq = NEGINF;
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

      pairedscore = psiA_ih + psiB_kq; /* TODO: add sigma_tau() score too */

      assert(pairedscore >= 0);

      /* note there appears to be an error in the paper
       * in this equation; it should be h-1 and q-1 as
       * here not +1 and +1 (j+1 and l+1 in the paper's
       * formulation).
       */
      sm = bpa_dynprogm_array(i+1, h-1, k+1, q-1, S) + pairedscore;
      shq = sm + bpa_dynprogm_array(h+1, j, q+1, l, S);
      if (shq > max_shq)
        max_shq = shq;
    }
  }
  score = MAX(score, max_shq);

  bpa_log_msg(funcname, "S\t%d\t%d\t%d\t%d\t%lg\n",i,j,k,l,score);
  S[INDEX4D(i,j,k,l,n1,n2)] = score;
#ifdef USE_INSTRUMENT
  bpastats[0].count_S++;
#endif
  return score;
}

