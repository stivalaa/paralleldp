/*****************************************************************************
 * 
 * File:    simpletest.c
 * Author:  Alex Stivala, SET code taken from Stuckey & Garcia de la Banda
 * Created: April 2009
 *
 * Test program that does same as oahttslftest but actrually does NOT
 * do any insertes or lookups (i.e. no memory access) to see how it
 * affects artificats on UltraSPARC T1 scalability
 * 
 * Usage:
 *    simpletest [numthreads]
 *
 *
 * $Id: simpletest.c 3184 2010-01-01 04:10:13Z alexs $
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <sys/time.h>

#include "bpautils.h"

#define MAX_NUM_THREADS 256

#define NUM_INSERTIONS  10000000


/*
 *TODO FIXME 
 *  this actually only uses the low 64 bits, doesn't even store high 
 */


typedef unsigned long long int uint64_t;


/***************************************************************************
 *
 * Definition of SET type and operations on it implemented as bit operations
 *
 ***************************************************************************/
#define MAXBIT 64


static int unioncount[MAX_NUM_THREADS];
static int bitscount[MAX_NUM_THREADS];
  
typedef struct twothings {
  uint64_t high, low;
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



/***************************************************************************
 *
 * static functions
 *
 ***************************************************************************/

static void *integer_ops(void *threadarg)
{
  SET s,snew;
  thread_data_t *mydata = (thread_data_t *)threadarg;
  unsigned int seed = mydata->thread_id * time(NULL);
  uint64_t  value,newvalue;
  uint64_t ivalue;
  int q;
  int i;

  for (q = 0; q < mydata->num_insertions; q++)
  {
    /*s.low = ((uint32_t)rand_r(&seed) << 31 | (uint32_t)rand_r(&seed)) + 1;  */
    for  (i = 0; i < 100; i++)
      s.low += 0xdeadbeef * 87234  + 32;
    if (s.low == 0)
        s.low = 1;
    s.high = 0;
    if (1)
    {

      /*snew.low = ((uint32_t)rand_r(&seed) << 31 | (uint32_t)rand_r(&seed)) + 1;  */
      snew.low = 0xfeebdaed * 23423  - 993;

      if (snew.low == 0)
          snew.low = 1;
      snew.high = 0;
      if (TRUE)
      {
          snew.high= value*q;
      }

      value = (uint64_t)s.low;
    }
    else
    {
        value=9999;
    }
  }
  return NULL;
}



/***************************************************************************
 *
 * main
 *
 ***************************************************************************/

int main(int argc, char *argv[])
{
  int t;
  int num_threads;
  int rc;
  pthread_t threads[MAX_NUM_THREADS];
  struct timeval start_timeval,end_timeval,elapsed_timeval;  
  int etime;
  int c,i;
  bool readstdin = FALSE;
  int num_actions;


  c = 1;
  for (i = 0; i < MAXBIT; i++) {
    bit[i] = c;
    c *= 2;
  }
  
  if (argc == 1)
    num_threads = 1;
  else if (argc == 2)
  {
    num_threads = atoi(argv[1]);
  }
  else
  {
    fprintf(stderr, "usage: %s [numthreads]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  if (num_threads >  MAX_NUM_THREADS)
  {
    fprintf(stderr, "max number of threads is %d\n", MAX_NUM_THREADS);
    exit(1);
  }

  gettimeofday(&start_timeval, NULL);


    for (t = 0; t < num_threads; t++)
    {
      thread_data[t].thread_id = t;
      thread_data[t].num_insertions = NUM_INSERTIONS / num_threads;
  
      if ((rc = pthread_create(&threads[t], NULL, integer_ops, 
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

  pthread_exit(NULL);
}

