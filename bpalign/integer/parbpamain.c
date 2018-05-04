/*****************************************************************************
 *
 * File:    parbpamain.c
 * Author:  Alex Stivala
 * Created: April 2009
 *
 * $Id: parbpamain.c 2699 2009-07-25 03:23:06Z astivala $
 *
 * This is the main() for a threaded dp implementation of
 * base pair probability alignment. Pairwise RNA structural alignment by
 * finding maximum-weight common secondary structure between two base pairing
 * probability matrices (as per McCaskill 1990).
 *
 * Note we build the version that uses the nbds library into a separate
 * executable than the version using the httslf.c module
 * in order not to overflow the .bss due with static data (hash table);
 * they both use this module as main().
 *
 * Usage: parbpalign [-avsz] [ -t num_threads | -b ] file1.bplist file2.bplist
 *
 *   Input files are sequence and base pair probability list output from
 *   the rnafold2list.py script (which extracts it from the _dp.ps output
 *   of the Vienna RNA package RNAfold -p program).
 *
 *   Output is to stdout (and stderr for errors/debug).
 *
 *  -s prints statistics; length of sequences, number of arcs, cells computed
 *  -v             : write verbose debugging output to stderr.
 *  -t num_threads : use threaded implementation with num_threads threads
 *  -b             : use bottom-up implementeation rather than top-down
 *  -a             : use top-down implementation but with array not hashtable
 *  -z             : do NOT randomize choices in multithread (-t) version
 *
 *
 * Platform and dependencies:
 *
 *   Use output of Vienna RNA package RNAfold (dot plot matrices)
 *   for base pair probability matrices. Developed with Vienna RNA package
 *   version 1.6.3.
 *
 *   Developed on Linux 2.6.22 (x86) with gcc 4.1.3.
 * 
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "bpacommon.h"
#include "bpautils.h"
#include "bpaglobals.h"
#include "bpaparse.h"
#include "bpaipsilist.h"
#include "bpadynprog_cpu.h"
#include "bpadynprog_hashthread.h"
#include "ht.h"
#include "bpastats.h"



/*****************************************************************************
 *
 * static functions
 *
 *****************************************************************************/





/*****************************************************************************
 *
 * external functions
 *
 *****************************************************************************/

counter_t total_count_S = 0, total_count_dynprogm_entry = 0, 
  total_count_dynprogm_entry_notmemoed = 0, num_keys = 0;

/*
 * bpalign()
 *
 * The main program logic of the bpalign program.
 * Read the sequences and base pair probability lists from the input files
 * and compute their optimal alignment.
 *
 * Parameters:
 *    filename1 - filename of first sequence to align (rnafold2list.py output)
 *    filename2 - filename of 2nd sequence to alignt (rnafold2list.py output)
 *
 *      Uses global data:
 *                  readonly:
 *                    use_threading - use threaded implementation
 *                    use_bottomup  - use bottom-up implementation
 *                    num_threads   - number of threads to use
 *                    use_array     - use array not hashtable on top-down
 *                    printstats    - print stats about data
 *                  read/write:
 *                    seqA    - first sequence
 *                    seqB    - second sequence
 *                    seqlenA - length of first sequence
 *                    seqlenB - length of second sequence
 *                    pairlistA - list of (i,j,p) for 1st seq
 *                    pairlistB - list of (i,j,p) for 2nd seq
 *                    paircountA - length of pairlistA
 *                    paircountB - length of pairlistB
 *                    ipsilistA - (j,psi) lists indexed by i for 1st seq
 *                    ipsilistB - (j,psi) lists indexed by i for 2nd seq
 *                    
 *
 * Return value:
 *    0 if successful,
 *    nonzero on error.
 */
static int bpalign(const char *filename1, const char *filename2)
{
  static const char *funcname = "bpalign";
  basepair_t *bplistA, *bplistB;
  int bplenA, bplenB; /* length of bplists */
  ipsi_element_t *seripsiA, *seripsiB;
  int ld_seripsiA, ld_seripsiB; /* leading dimension of seripsi arrays */
  myint64_t score, *dev_score;
  int *dev_seqA, *dev_seqB;
  ipsi_element_t *dev_seripsiA, *dev_seripsiB;
  myint64_t *dev_S;
  volatile myint64_t *matrixS;

  int otime, ttime, etime;
  struct rusage starttime,totaltime,runtime,endtime,opttime;
  struct timeval start_timeval,end_timeval,elapsed_timeval;

  extern bpastats_t bpastats[]; /* for bpadynprog_single only */

  /*
   * read sequences and base pair probabilities and build data structures
   */

  if (!(bplistA = bpa_read_basepairs(filename1, PMIN, &bpaglobals.seqA,
                                     &bplenA)))
  {
    bpa_error_msg(funcname, "could not read basepairs from %s\n",filename1);
    return -1;
  }

  if (!(bplistB = bpa_read_basepairs(filename2, PMIN,
                                     &bpaglobals.seqB, &bplenB)))
  {
    bpa_error_msg(funcname, "could not read basepairs from %s\n",filename2);
    return -1;
  }

  bpaglobals.seqlenA = strlen(bpaglobals.seqA);
  bpaglobals.seqlenB = strlen(bpaglobals.seqB);
  bpaglobals.pairlistA = bplistA;
  bpaglobals.pairlistB = bplistB;
  bpaglobals.paircountA = bplenA;
  bpaglobals.paircountB = bplenB;

/*  bpa_dump_bp_list(bplenA, bplistA, seqA);  */
/*  bpa_dump_bp_list(bplenB, bplistB, seqB); */

  bpaglobals.ipsilistA = bpa_pairlist_to_ipsilist(bplistA, bplenA,
                                                  bpaglobals.seqlenA);
  bpaglobals.ipsilistB = bpa_pairlist_to_ipsilist(bplistB, bplenB,
                                                  bpaglobals.seqlenB);

/*  bpa_dump_ipsilist(bpaglobals.ipsilistA, bpaglobals.seqlenA);   */
/*  bpa_dump_ipsilist(bpaglobals.ipsilistB, bpaglobals.seqlenB);   */


  if (bpaglobals.verbose)
  {
    fprintf(stderr, "seq A len = %d\n", bpaglobals.seqlenA);
    fprintf(stderr, "seq B len = %d\n", bpaglobals.seqlenB);
    fprintf(stderr, "seq A arcs = %d\n", bplenA);
    fprintf(stderr, "seq B arcs = %d\n", bplenB);
  }

  seripsiA = bpa_serialize_ipsilist(bpaglobals.ipsilistA, bpaglobals.seqlenA,
                                    &ld_seripsiA);
  seripsiB = bpa_serialize_ipsilist(bpaglobals.ipsilistB, bpaglobals.seqlenB,
                                    &ld_seripsiB);

  if (bpaglobals.verbose)
  {
    bpa_dump_seripsilist(seripsiA, bpaglobals.seqlenA, ld_seripsiA);  
    bpa_dump_seripsilist(seripsiB, bpaglobals.seqlenB, ld_seripsiB);  
  }

  if (bpaglobals.use_bottomup || bpaglobals.use_array)
  {
    /* allocate workarea for dp matrix */
    matrixS =  (myint64_t *)bpa_malloc(
      bpaglobals.seqlenA * bpaglobals.seqlenA * 
      bpaglobals.seqlenB * bpaglobals.seqlenB *
      sizeof(myint64_t));

    /* need to set array elements all to NEGINF for top-down */
    if (!bpaglobals.use_bottomup)
    {
      int i,j,k,l;
      for (i = 0; i < bpaglobals.seqlenA; i++)
        for (j = 0; j < bpaglobals.seqlenA; j++)
          for (k = 0; k < bpaglobals.seqlenB; k++)
            for (l = 0; l < bpaglobals.seqlenB; l++)
              matrixS[INDEX4D(i,j,k,l,bpaglobals.seqlenA,bpaglobals.seqlenB)] = NEGINF;

    }
  }

  gettimeofday(&start_timeval, NULL);
  getrusage(RUSAGE_SELF, &starttime);

  if (bpaglobals.use_bottomup)
  {
    /* do the d.p. with basic cpu implementation */
    dynprog_cpu(bpaglobals.seqlenA, bpaglobals.seqlenB,
            bpaglobals.seqA, bpaglobals.seqB,
            seripsiA, ld_seripsiA,
            seripsiB, ld_seripsiB,
            bpaglobals.gamma, MINLOOP,
            matrixS,
            &score);
  }
  else if (bpaglobals.use_threading)
  {
    /* threading, top-down  implementation */
    if (bpaglobals.use_array)
      score = bpa_dynprogm_thread_array_master(0, bpaglobals.seqlenA-1, 0, 
                                               bpaglobals.seqlenB-1,
                                               matrixS);
    else
    {
      score = bpa_dynprogm_thread_master(0, bpaglobals.seqlenA-1, 0, 
                                         bpaglobals.seqlenB-1);
    }
  }
  else 
  {
    /* no threading, top-down implementation */
    if (bpaglobals.use_array)
      score = bpa_dynprogm_array(0, bpaglobals.seqlenA-1, 0, bpaglobals.seqlenB-1, matrixS);
    else
    {
      score = bpa_dynprogm(0, bpaglobals.seqlenA-1, 0, bpaglobals.seqlenB-1);
    }
    total_count_S = bpastats[0].count_S;
    total_count_dynprogm_entry = bpastats[0].count_dynprogm_entry;
    total_count_dynprogm_entry_notmemoed = bpastats[0].count_dynprogm_entry_notmemoed;
#ifdef USE_INSTRUMENT
    num_keys = (counter_t)oahttslf_total_key_count();
#endif
  }

  getrusage(RUSAGE_SELF, &endtime);
  gettimeofday(&end_timeval, NULL);
  timeval_subtract(&elapsed_timeval, &end_timeval, &start_timeval);
  /* timeval_subtract(&endtime,&starttime,&runtime); */
  runtime = endtime;
  ttime = 1000 * runtime.ru_utime.tv_sec + runtime.ru_utime.tv_usec/1000 
          + 1000 * runtime.ru_stime.tv_sec + runtime.ru_stime.tv_usec/1000;
  etime = 1000 * elapsed_timeval.tv_sec + elapsed_timeval.tv_usec/1000;

  if (bpaglobals.printstats)
  {
    /* print in one line space-sepearted format for later parsing */
    /* score user+system_cpu_time elased_time lenA lenB arcsA arcsB re hc hn */
    printf("%lld %d %d %d %d %d %d %lu %lu %lu\n", score, ttime, etime,
            bpaglobals.seqlenA, bpaglobals.seqlenB, 
            bplenA, bplenB,
            total_count_dynprogm_entry - total_count_dynprogm_entry_notmemoed,
           total_count_dynprogm_entry_notmemoed, num_keys);
  }
  else
  {
    printf("score = %lld\n", score);
  }


  /* free memory */
  free(bpaglobals.seqA);
  free(bpaglobals.seqB);
  free(bpaglobals.pairlistA);
  free(bpaglobals.pairlistB);
  bpa_free_ipsilist(bpaglobals.ipsilistA, bpaglobals.seqlenA);
  bpa_free_ipsilist(bpaglobals.ipsilistB, bpaglobals.seqlenB);
  free(seripsiA);
  free(seripsiB);
  /* TODO free hash table entries */

  return 0;
}


/*
 * print usage message and exit
 *
 * Parameters:
 *     program - name of program (usu argv[0])
 */
static void usage(const char *program)
{
  fprintf(stderr,
          "usage: %s  [-svaz] [-t num_threads | -b] file1.bplist file2_bplist\n"
          "   -s  :  write instrumentation data to stdout\n"
          "   -v  :  write verbose debug information to stderr\n"
          "   -t num_threads  : use threaded implementation\n"
          "   -a  :  usee array not hashtable for top-down implementations\n"
          "   -b  :  use bottom-up not top-down dynamic programming\n",
          "   -z  :  do NOT randomize choices in multithreaded version\n",
          program);
  exit(EXIT_FAILURE);
}


/*
 * main
 */
int main(int argc, char *argv[])
{
  int c;
  char *filename1, *filename2;
  int exit_status;
  
  bpa_set_verbose(FALSE);

  /* process command line options */

  while ((c = getopt(argc, argv, "ast:bvzh?")) != -1)
  {
    switch (c)
    {
      case 'a':
        bpaglobals.use_array = TRUE; /* use array not hashtable for top-down */
        break;

      case 'b':
        bpaglobals.use_bottomup = TRUE; /*use bottom-up not top-down */
        break;
        
      case 't':
        bpaglobals.use_threading = TRUE; /* used threads */
        if (atoi(optarg) < 1)
        {
          fprintf(stderr, "number of threads must be >= 1\n");
          usage(argv[0]);
        }
        else if (atoi(optarg) > MAX_NUM_THREADS)
        {
          fprintf(stderr, "maximum number of threads is %d\n", MAX_NUM_THREADS);
          usage(argv[0]);
        }
        bpaglobals.num_threads = atoi(optarg);
        break;

      case 'v':
        bpaglobals.verbose = TRUE;
        bpa_set_verbose(TRUE);
        break;

      case 's':
        bpaglobals.printstats = TRUE;
        break;

      case 'z':
        bpaglobals.use_random = FALSE;
        break;

      case 'h':
      case '?':
        usage(argv[0]);
        break;

      default:
        /* should not happen */
        usage(argv[0]);
        break;
    }
  }

  /* we should have exactly two command line parameters (two input files) */
  if (optind == argc - 2)
  {
    filename1 = argv[optind];
    filename2 = argv[optind+1];
  }
  else
  {
    usage(argv[0]);
  }
  
  if (bpaglobals.use_bottomup && bpaglobals.use_threading)
  {
    fprintf(stderr, "cannot use threading (-t) with bottom-up (-b)\n");
    usage(argv[0]);
  }

  if (bpaglobals.use_array && bpaglobals.use_bottomup)
    fprintf(stderr,
            "WARNING: -a (use array) ignored with -b: bottom-up always uses array\n");

  /* main toplevel program logic is in bpalign() */

  if (bpalign(filename1, filename2) == 0)
    exit_status = EXIT_SUCCESS;
  else
    exit_status = EXIT_FAILURE;

  return exit_status;
}
