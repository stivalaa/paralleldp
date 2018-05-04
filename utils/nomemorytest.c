/*****************************************************************************
 * 
 * File:    nomemorytest.c
 * Author:  Alex Stivala, SET code taken from Stuckey & Garcia de la Banda
 * Created: April 2009
 *
 * Test program that does same as oahttslftest but actrually does NOT
 * do any insertes or lookups (i.e. no memory access) to see how it
 * affects artificats on UltraSPARC T1 scalability
 * 
 * Usage:
 *    nomemorytest [-r] [numthreads]
 *
 *     -r : read actions from stdin
 *
 * $Id: nomemorytest.c 3184 2010-01-01 04:10:13Z alexs $
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <sys/time.h>

#include "oahttslf.h"

#define MAX_NUM_THREADS 256

#define NUM_INSERTIONS  10000000


/*
 *TODO FIXME 
 *  this actually only uses the low 64 bits, doesn't even store high 
 */




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

typedef enum
{
  ACTION_LOOKUP = 0,
  ACTION_INSERT = 1
} action_e;

typedef struct action_s 
{
    action_e action;
    unsigned long long key;
    unsigned long long value;
} action_t;

#define MAX_ACTIONS 25000

static action_t actions[MAX_ACTIONS];

/* read space-separated insert|lookup key (hex) value (dec) pairs from fp one per line
   e.g.

lookup 1002861048
insert 1002061048 46
...
insert 2F3D2EF7A7B 45
*/
static int read_actions(FILE *fp, action_t kv[])
{
  int i = 0;
  char actbuf[8];

  while (i < MAX_ACTIONS && !feof(fp))
  {
    fgets(actbuf, 7, fp);
    actbuf[7] = '\0';
    if (!strcmp(actbuf, "lookup"))
    {
      kv[i].action = ACTION_LOOKUP;
      if (fscanf(fp, " %16llX\n", &kv[i].key) != 1)
      {
          fprintf(stderr, "bad lookup line %d\n", i);
          exit(EXIT_FAILURE);
      }
    }
    else if (!strcmp(actbuf, "insert"))
    {
      kv[i].action = ACTION_INSERT;
      if (fscanf(fp, " %16llX %llu\n", &kv[i].key, &kv[i].value) != 2)
      {
          fprintf(stderr, "bad insert line %d\n", i);
          exit(EXIT_FAILURE);
      }
    }
    i++;
  }
  return i;
}


static void *actions_fromlist(void *threadarg)
{
  SET s;
  thread_data_t *mydata = (thread_data_t *)threadarg;
  unsigned int seed = mydata->thread_id;
  uint64_t  value;
  uint64_t ivalue,dummy;
  int q;
  bool found;

  for (q = 0; q < mydata->num_insertions; q++)
  {
#ifdef DEBUG
      printf("thread %d action %d: ", mydata->thread_id, q);
#endif
    switch (actions[q].action)
    {
      case ACTION_LOOKUP:
#ifdef DEBUG
          printf("lookup(%16llX)", actions[q].key);
#endif
          found = oahttslf_lookup(actions[q].key, &ivalue);
#ifdef DEBUG
          if (found)
              printf(" = %llu\n", ivalue);
          else
              printf(" not found\n");
#endif
          break;

      case ACTION_INSERT:
#ifdef DEBUG
          printf("insert(%16llX, %llu)\n", actions[q].key, actions[q].value);
#endif
          oahttslf_insert(actions[q].key, actions[q].value, mydata->thread_id);
          break;

      default:
          fprintf(stderr, "Bad action %d\n", actions[q].action);
          exit(EXIT_FAILURE);
    }
  }
  return NULL;
}

static void *insert_random(void *threadarg)
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
    if (!strcmp(argv[1], "-r"))
      readstdin = TRUE;
    else
      num_threads = atoi(argv[1]);
  }
  else if (argc == 3)
  {
    if (strcmp(argv[1], "-r"))
    {
      fprintf(stderr, "usage: %s [-r] [numthreads]\n", argv[0]);
      exit(EXIT_FAILURE);
    }
    readstdin = TRUE;
    num_threads = atoi(argv[2]);
  }
  else
  {
    fprintf(stderr, "usage: %s [-r] [numthreads]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  if (num_threads >  MAX_NUM_THREADS)
  {
    fprintf(stderr, "max number of threads is %d\n", MAX_NUM_THREADS);
    exit(1);
  }

  if (readstdin)
    num_actions = read_actions(stdin, actions);


  gettimeofday(&start_timeval, NULL);


  if (readstdin)
  {
    for (t = 0; t < num_threads; t++)
    {
      thread_data[t].thread_id = t;
      thread_data[t].num_insertions = num_actions;
  
      if ((rc = pthread_create(&threads[t], NULL, actions_fromlist, 
                               (void *)&thread_data[t])))
      {
        fprintf(stderr, "pthread_create() failed (%d)\n", rc);
        exit(EXIT_FAILURE);
      }
    }
  }
  else
  {
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

#ifdef DEBUG 
  if (!oahttslf_validate())
  {
      fprintf(stderr, "hash table validation failed\n");
      exit(EXIT_FAILURE);
  }

  oahttslf_printstats();
#endif

#ifdef USE_CONTENTION_INSTRUMENT
  printf("total retry count = %ld\n", oahttslf_total_retry_count());
#endif


  pthread_exit(NULL);
}

