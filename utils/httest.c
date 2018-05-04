/*****************************************************************************
 * 
 * File:    httest.c
 * Author:  Alex Stivala, SET code taken from Stuckey & Garcia de la Banda
 * Created: April 2009
 *
 * Test harness for separate chaining hash table.
 * 
 *
 * $Id: httest.c 3094 2009-12-21 00:51:52Z alexs $
 *
 *****************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>

#include "ht.h"


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

/* hash function for ht - takes SET* and returns hash value for the set */
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

  i = q & (HT_SIZE - 1); /* depends on HT_SIZE being 2^n */

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

static void *insert_random(void)
{
  SET s;
  unsigned int seed = time(NULL);
  int *pvalue, value;
  int q;

  srand(seed);

  for (q = 0; q < NUM_INSERTIONS; q++)
  {
    s.low = rand() << 16 | rand();
    s.high = 0;
    if (!(pvalue = ht_lookup(&s)))
    {
      value = (int)s.low;
      ht_insert(&s, &value);
      /*assert(*(int *)ht_lookup(&s) == (int)s.low);*/

      s.low = rand() << 16 | rand();
      s.high = 0;
      if ((pvalue = ht_lookup(&s)))
      {
        assert(*(int *)ht_lookup(&s) == (int)s.low);
      }
    }
    else
    {
      assert(*pvalue == (int)s.low);
    }
  }
}


int main(int argc, char *argv[])
{
  int t;
  int rc;
  struct timeval start_timeval,end_timeval,elapsed_timeval;  
  int etime;
  int c,i;


  c = 1;
  for (i = 0; i < MAXBIT; i++) {
    bit[i] = c;
    c *= 2;
  }
  
  if (argc != 1)
  {
    fprintf(stderr, "usage: %s\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  gettimeofday(&start_timeval, NULL);

  ht_initialize(sizeof(SET), sizeof(int), hash_function,
                    NULL, setmatch, NULL);

  insert_random();

  gettimeofday(&end_timeval, NULL);
  timeval_subtract(&elapsed_timeval, &end_timeval, &start_timeval);
  etime = 1000 * elapsed_timeval.tv_sec + elapsed_timeval.tv_usec/1000;
  printf("elapsed time %d ms\n", etime);

  if (!ht_validate())
  {
      fprintf(stderr, "hash table validation failed\n");
      exit(EXIT_FAILURE);
  }

  ht_printstats();
  
  exit(EXIT_SUCCESS);
}

