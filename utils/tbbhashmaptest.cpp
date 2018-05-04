/*****************************************************************************
 * 
 * File:    tbbhashmaptest.cp
 * Author:  Alex Stivala, SET code taken from Stuckey & Garcia de la Banda,
 *          Intel TBB code based on example in manual.
 * Created: May 2009
 *
 * Test harness for Intel TBB concurrent hash map
 * 
 * Usage:
 *    tbbhashmaptest [numthreads]
 *
 * $Id: tbbhashmaptest.cpp 2398 2009-05-16 04:12:48Z astivala $
 *
 * Tested on:
 * Linux charikar.csse.unimelb.edu.au 2.6.9-78.0.5.ELsmp #1 SMP Wed Sep 24 05:40:24 EDT 2008 x86_64 x86_64 x86_64 GNU/Linux
 * TBB21_INSTALL_DIR=/home/charikar/pgrad/astivala/tbb21_20080605oss
 * . ~/tbb21_20080605oss/em64t/cc3.4.3_libc2.3.4_kernel2.6.9/bin/tbbvars.sh
 * TBB_ARCH_PLATFORM=em64t/cc3.4.3_libc2.3.4_kernel2.6.9
 *
 *****************************************************************************/

#include <string>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include "tbb/concurrent_hash_map.h"
#include "tbb/scalable_allocator.h"
#include <assert.h>
#include <pthread.h>
#include <sys/time.h>
#include "bpautils.h"
#include "tbbhashmap.h"

using namespace tbb;
using namespace std;


#define MAX_NUM_THREADS 256

#define NUM_INSERTIONS  10000000




/***************************************************************************
 *
 * Definition of SET type and operations on it implemented as bit operations
 *
 ***************************************************************************/
#define MAXBIT 64

static int unioncount[MAX_NUM_THREADS];
static int bitscount[MAX_NUM_THREADS];
  
typedef struct twothings {
  long long high, low;
} SET;

#define in(e,S) ( (e) < 64 ? bit[e] & S.low : bit[(e)-64] & S.high ) 

#define union(S1,S2,S) { S.high = S1.high | S2.high ; S.low = S1.low | S2.low; } 

#define nonempty(S) ( S.low | S.high )

#define empty(S) ( !(S.low | S.high) )

#define emptyset(S)  { S.high = 0; S.low = 0; }

#define neg(S,N) { N.high = ~S.high; N.low = ~S.low; }

#define intersect(S1,S2,S) {S.high = S1.high & S2.high; S.low = S1.low & S2.low; }

#define subset(S1,S2) ( (S1.low == (S1.low & S2.low)) && (S1.high == (S1.high & S2.high)))

#define removebit(e,S0,S) { S = S0; if (e < 64) S.low = S.low & ~bit[e]; else S.high = S.high & ~bit[e-64]; } 

#define addbit(e,S0,S) { S = S0; if (e < 64) S.low = S.low | bit[e]; else S.high = S.high | bit[e-64]; } 

#define notequal(S1,S2) (S1.low != S2.low || S1.high != S2.high) 
//#define equal(S1,S2) (S1.low == S2.low && S1.high == S2.high) 

#define intersects(S1,S2) ( (S1.low & S2.low) || (S1.high & S2.high)) 
int totalcost;

#define universe(S) { S.high = -1; S.low = -1; }

long long bit[64];

static void printarr(int SIZE, int a[])
{
  int i;
  for (i = 0; i < SIZE; i++)
    printf("%5d",a[i]);
  printf("\n");
}

static void printbit(int nbits, SET n)
{
  int  i;
  for (i = nbits-1; i >= 0; i--)
    if (in(i,n)) printf("1");
    else printf("0");
}




static int nbitsp(SET as, int thread_id)
{
  int c;
  long long word;

  word = as.high;
  for (c = 0; word; c++)
    {
      word &= word - 1; /* clear the least significant bit set */
    }
  word = as.low;
  for (; word; c++)
    {
      word &= word - 1; /* clear the least significant bit set */
    }
  bitscount[thread_id]++;
  return c;
}



/***************************************************************************
 *
 * thread data
 *
 ***************************************************************************/


typedef struct thread_data_s
{
    int thread_id;
    int num_insertions;
} thread_data_t;


static thread_data_t thread_data[MAX_NUM_THREADS];



static void *insert_random(void *threadarg)
{
  _SET s;
  thread_data_t *mydata = (thread_data_t *)threadarg;
  unsigned int seed = mydata->thread_id;
  int value;
 
  int q;

  for (q = 0; q < mydata->num_insertions; q++)
  {
    s.low = rand_r(&seed) << 16 | rand_r(&seed);
    s.high = 0;


    if (!tbbhashmap_haskey(s))
    {
//      printf("inserting %llu\n", s.low);
      value = (int)s.low;
      tbbhashmap_insert(s, value);
//      assert(*(int *)httslf_lookup(&s) == (int)s.low);
      assert(tbbhashmap_lookup(s) == value);
    }
    else
    {
//      printf("found %llu\n", s.low);
      assert(tbbhashmap_lookup(s) == (int)s.low);
    }
  }
  return NULL;
}


int main(int argc, char *argv[])
{
  int t;
  int num_threads;
  int rc;
  pthread_t threads[MAX_NUM_THREADS];
  struct timeval start_timeval,end_timeval,elapsed_timeval;  
  int etime;
  int c,i;


  c = 1;
  for (i = 0; i < MAXBIT; i++) {
    bit[i] = c;
    c *= 2;
  }
  
  if (argc == 1)
    num_threads = 2;
  else if (argc == 2)
    num_threads = atoi(argv[1]);
  else
  {
    fprintf(stderr, "usage: %s [numthreads]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  gettimeofday(&start_timeval, NULL);

  for (t = 0; t < num_threads; t++)
  {
    thread_data[t].thread_id = t;
    thread_data[t].num_insertions = NUM_INSERTIONS / num_threads;

    if ((rc = pthread_create(&threads[t], NULL, insert_random, 
                             (void *)&thread_data[t])))
    {
      fprintf(stderr, "pthread_create() failed (%d)\n", rc);
      exit(EXIT_FAILURE);
    }
  }

  /* wait for threads */
  for (t = 0; t < num_threads; t++)
  {
    if ((rc = pthread_join(threads[t], NULL)))
    {
      fprintf(stderr, "pthread_join() failed (%d)\n", rc);
      exit(EXIT_FAILURE);
    }
  }

  gettimeofday(&end_timeval, NULL);
  timeval_subtract(&elapsed_timeval, &end_timeval, &start_timeval);
  etime = 1000 * elapsed_timeval.tv_sec + elapsed_timeval.tv_usec/1000;
  printf("elapsed time %d ms\n", etime);

  exit(EXIT_SUCCESS);
}

