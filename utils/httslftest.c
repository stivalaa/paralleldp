/*****************************************************************************
 * 
 * File:    httslftest.c
 * Author:  Alex Stivala, SET code taken from Stuckey & Garcia de la Banda
 * Created: April 2009
 *
 * Test harness for thread-safe lock-free separate chaining hash table.
 * 
 * Usage:
 *    httslftest [numthreads]
 *
 * $Id: httslftest.c 3165 2009-12-29 03:23:30Z alexs $
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <pthread.h>
#include <time.h>
#include <sys/time.h>

#include "httslf.h"

#define MAX_NUM_THREADS 256

#define NUM_INSERTIONS  10000000

#define USE_GOOD_HASH 1


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
#define equal(S1,S2) (S1.low == S2.low && S1.high == S2.high) 

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


#ifdef USE_GOOD_HASH
/*
  hash a 64 bit value into 32 bits. From:
  (Thomas Wang, Jan 1997, Last update Mar 2007, Version 3.1)
  http://www.concentric.net/~Ttwang/tech/inthash.htm
  (found by reference in NIST Dictionary of Algorithms and Data Structures)
*/
static unsigned long hash6432shift(unsigned long long key)
{
  key = (~key) + (key << 18); /* key = (key << 18) - key - 1; */
  key = key ^ (key >> 31);
  key = key * 21; /* key = (key + (key << 2)) + (key << 4); */
  key = key ^ (key >> 11);
  key = key + (key << 6);
  key = key ^ (key >> 22);
  return (unsigned long) key;
}

#endif

/* hash function for httslf - takes SET* and returns hash value for the set */
static unsigned int hash_function(const void *vp) {
  unsigned int i;
  unsigned long highhash,lowhash,q;
  SET p = *(SET *)vp;

#ifdef USE_GOOD_HASH
  highhash = hash6432shift(p.high);
  lowhash = hash6432shift(p.low);
  q = lowhash ^ highhash; /* FIXME: should have a hash12832shift() instead */
#else
  q = p.low;
#endif

  i = q & (HTTSLF_SIZE - 1); /* depends on HTTSLF_SIZE being 2^n */

  return i;
}


/* key comparison function for httslf - return nonzero iff
   both sets are the same */
static int setmatch(const void *vpset1, const void *vpset2)
{
#ifdef SPARC
  /* MUCH faster to memcmp() than proper cast & compare on Solaris SPARC */
  /* TODO: why? It uses a huge amount of system (not user) CPU time
     to do it by casting to SET* and comparing with set1->low==set2.low etc.,
     while this mecmp() does not. can't work out why */
  return !memcmp(vpset1, vpset2, sizeof(SET));
#else
  /* but on x86-64 (Linux) it is faster to do it this way,
     to the extent that doing it with memcmp() like Solaris SPARC results
     in slowdown even for 4 (and probably more, but that's all I have) cores.
     No idea why. */
  const SET *pset1 = (const SET *)vpset1;
  const SET *pset2 = (const SET *)vpset2;
  return pset1->low == pset2->low && pset1->high == pset2->high;
#endif
}


static void *insert_random(void *threadarg)
{
  SET s,snew;
  thread_data_t *mydata = (thread_data_t *)threadarg;
  unsigned int seed = mydata->thread_id * time(NULL);
  int *pvalue,*pnewvalue, value,newvalue;
  int q;

  for (q = 0; q < mydata->num_insertions; q++)
  {
    s.low = rand_r(&seed) << 16 | rand_r(&seed);
    s.high = 0;
    if (!(pvalue = httslf_lookup(&s)))
    {
      snew.low = rand_r(&seed) << 16 | rand_r(&seed);
      snew.high = 0;
      if ((pnewvalue = httslf_lookup(&snew)))
      {
        assert(*pnewvalue == (int)snew.low);
      }

      value = (int)s.low;
      httslf_insert(&s, &value);
      /*assert(*(int *)httslf_lookup(&s) == (int)s.low);*/

    }
    else
    {
      assert(*pvalue == (int)s.low);
    }
  }
/*  pthread_exit(NULL); */
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

  httslf_initialize(sizeof(SET), sizeof(int), hash_function,
                    NULL, setmatch, NULL);


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

  if (!httslf_validate())
  {
      fprintf(stderr, "hash table validation failed\n");
      exit(EXIT_FAILURE);
  }

  httslf_printstats();

  pthread_exit(NULL);
}

