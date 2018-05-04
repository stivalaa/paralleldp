/*****************************************************************************
 * 
 * File:    oahttslf.c
 * Author:  Alex Stivala
 * Created: April 2009
 *
 * Open addressing (closed hashing) thread-safe lock-free hash table.
 * Uses linear probing.
 * 
 * gcc version 4.1.0 or greater is required, in order to use the
 * __sync_val_compare_and_swap()
 * and __sync_val_fetch_and_add() builtins
 * (this module was developed on Linux 2.6.22 (x86) with gcc 4.1.3),
 * except on Solaris, where we use SUNWSpro compiler and
 * cas32()/caslong()/cas64()/casptr() in 
 * /usr/include/sys/atomic.h)
 *
 * $Id: oahttslf.c 3148 2009-12-27 04:15:31Z alexs $
 *
 *
 * Preprocessor symbols:
 *
 *
 * USE_GOOD_HASH  - use mixing hash function rather than trivial one
 * DEBUG          - include extra assertion checks etc.
 * ALLOW_UPDATE  - allow insert to update value of existing key
 * USE_INSTRUMENT - compile in (per-thread) instrumentation counts.
 * USE_CONTENTION_INSTRUMENT - per-thread contention counts (only).
 *
 *****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "bpautils.h"
#include "oahttslf.h"
#include "cellpool.h"
#include "atomicdefs.h"


#define USE_GOOD_HASH
#define ALLOW_UPDATE

/* linear probing step size */
#define OAHTTSLF_PROBE_STEP 1


/*****************************************************************************
 *
 * types
 *
 *****************************************************************************/


typedef struct oahttslf_entry_s
{
    uint64_t key;   /* TODO need to be able to have a SET (128 bit) key */
    uint64_t value;
} oahttslf_entry_t;


/*****************************************************************************
 *
 * static data
 *
 *****************************************************************************/


/* The hash table */
/* Each entry is a key-value pair. Writing is synchornized with CAS logic */
/* Note we depend on the empty key/value being 0 since this is declared */
/* static and therefore initilized to zero */
/* TODO: change to allocate dynamically so can have multiple */
static oahttslf_entry_t hashtable[OAHTTSLF_SIZE];


/***************************************************************************
 *
 * Instrumentation, per thread (each thread only writes to its own array 
 * element), and the total_ values which are computed only in the master
 * thread, by summing over the per-thread arrays.
 * Note that using this can really slow thigns down (e.g. on PPC and Intel)
 * due to 'false sharing' so only compiled in if USE_INSTRUMENT is defined.
 * pthreads "thread local storage" 
 * seems way to complicated to even experiment with
 * (for all I know has the same problem anyway).
 * We have this because we have no global counter of keys inserted in order
 * to avoid sychronizing with CAS on a single counter which will cost us
 * in terms of therad scalability. And oahttslf_num_entries() which
 * iterates through to count things (etherefore costing nothing during
 * processing, but just at then end) worked out way too slow (on open
 * stacks experiments on mundara (SPARC) total time was 8 hours rather
 * than less than 1 hour for example). Hence had to go to all this
 * trouble to add thread_id parameter and this instrumentation.
 *
 ***************************************************************************/

#ifdef USE_INSTRUMENT
static unsigned int key_count[MAX_NUM_THREADS];
#endif
#ifdef USE_CONTENTION_INSTRUMENT
static unsigned int retry_count[MAX_NUM_THREADS];
#endif

/*****************************************************************************
 *
 * static functions
 *
 *****************************************************************************/


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


static unsigned int hash_function(uint64_t key) {
  unsigned int i;
  unsigned long q;

#ifdef USE_GOOD_HASH
  q = hash6432shift(key);
#else
  q = key;
#endif

  i = q & (OAHTTSLF_SIZE - 1); /* depends on OAHTTSLF_SIZE being 2^n */

  return i;
}



/*****************************************************************************
 *
 * external functions
 *
 *****************************************************************************/

volatile oahttslf_entry_t *oahttslf_getent(uint64_t key, bool *vacant);

/*
 * oahttslf_insert()
 *
 * Insert a key/value pair into the hashtable, or update the value
 * for existing key.
 *
 * Parameters:
 *    key   - key to insert
 *    value - value to insert for the key
 *    thread_id - our thread identifer (0,1,2,.. NOT pthread_t) used
 *                 only for instrumentation
 *
 * Return value:
 *    Value for the key prior to the new insertion (OAHTTSLF_EMPTY_VALUE
 *    for a new key)
 */
uint64_t oahttslf_insert(uint64_t key, uint64_t value, int thread_id)
{
  static const char *funcname = "oahttslf_insert";
  volatile oahttslf_entry_t *ent;
  uint64_t entkey, oldvalue;
  bool vacant;

  assert(key != OAHTTSLF_EMPTY_KEY);
  assert(value != OAHTTSLF_EMPTY_VALUE);

  ent = oahttslf_getent(key, &vacant);
  if (!ent)
  {
    bpa_fatal_error(funcname, "hash table full\n"); /* TODO expand table */
  }
  oldvalue = ent->value;
  entkey = ent->key;
  if (vacant)
  {
    if (CAS64(&ent->key, OAHTTSLF_EMPTY_KEY, key) != OAHTTSLF_EMPTY_KEY) {
#ifdef USE_CONTENTION_INSTRUMENT
      retry_count[thread_id]++;
#endif
      return oahttslf_insert(key, value, thread_id);  /* tail-recursive call to retry */
    }

    /* in lookup we test for OAHHTSLF_EMPTY_VALUE so if someone looks up 
       in anothe thread before the value is set here, we return key not found.
       NB depends on 64-bit atomic writes */
/*    if (ent->value != OAHTTSLF_EMPTY_VALUE)
      bpa_error_msg(funcname, "value = %llX\n", ent->value);     
      */
/*  __asm__ __volatile__ ("" : : : "memory");  */

    entkey = key ;

#ifdef USE_INSTRUMENT
    key_count[thread_id]++;  /* count new keys only, to get total in table */
#endif
  }
#ifdef DEBUG
  /*assert(key == ent->key);*/
  if (key != entkey)
  {
          fprintf(stderr, "OAHTTSLF ASSERTION FAILURE: key=%llX entkey=%llX\n",  key, entkey);
          exit(1);
  }
#endif
#ifdef ALLOW_UPDATE
  if (oldvalue == value)  /* shortcut to avoid expense of CAS instruction */
    return oldvalue;
  if (CAS64(&ent->value, oldvalue, value) != oldvalue)
    return oahttslf_insert(key, value, thread_id);  /* tail-recursive call to retry */
#else
  ent->value = value;
#endif
  return oldvalue;
}



/*
 * oahttslf_lookup()
 *
 * Get the value for a key from the hashtable
 *
 * Parameters:
 *     key -  key to look up
 *     value - (output) value for key,ony set if TRUE returned.
 * * Return value: *      TRUE if key found, FALSE otherwise.
 */
bool oahttslf_lookup(uint64_t key, uint64_t *value)
{
  volatile oahttslf_entry_t *ent;
  bool vacant;
  uint64_t val;

  ent = oahttslf_getent(key, &vacant);
  if (ent)
  {
    val = ent->value;
    if (!vacant && val != OAHTTSLF_EMPTY_VALUE)
    {
      *value = val;
      return TRUE;
    }
  }
  return FALSE;
}




/*
 * oahttslf_getent()
 *
 * Get the entry for a key from the hashtable
 *
 * Parameters:
 *     key -  key to look up
 *     vacant - (OUT) TRUE if reutrn pointer to entry for key
 *                    is not currently occupied by key
 *  Return value:
 *     pointer to entry with key, or for key (but currently empty) in hashtable
 *     or NULL if hashtable is full
 */
volatile oahttslf_entry_t *oahttslf_getent(uint64_t key, bool *vacant)
{
  unsigned int h;
  volatile oahttslf_entry_t *ent;
  int probes = 0;
  uint64_t entkey;

  h = hash_function(key);
  ent = &hashtable[h];
  entkey = ent->key;
  while (probes < OAHTTSLF_SIZE - 1 && entkey != key && entkey != OAHTTSLF_EMPTY_KEY)
  {
    ++probes;
    h = (h + OAHTTSLF_PROBE_STEP) & (OAHTTSLF_SIZE - 1); /*SIZE must be 2^n*/
    ent = &hashtable[h];
    entkey = ent->key;
  }
  if (probes >= OAHTTSLF_SIZE - 1)
    return NULL;
  else if (entkey == OAHTTSLF_EMPTY_KEY)
    *vacant = TRUE;
  else
    *vacant = FALSE;
  return ent;
}




/*
 * oahttslf_validate()
 *
 * Test for duplicate keys  -this should not happen
 *
 * Parameters:
 *    None
 *
 * Return value:
 *    0 if duplicate keys found else 1
 */
int oahttslf_validate(void)
{

  int i,j;

  for (i = 0; i < OAHTTSLF_SIZE; i++)
    if (hashtable[i].key != OAHTTSLF_EMPTY_KEY)
      for (j = i + 1; j < OAHTTSLF_SIZE; j++)
        if (hashtable[j].key == hashtable[i].key)
          return 0;
    
  return 1;
}

/* TODO function to free hash table entries */


/*
 * oahttslf_printstats()
 *
 *   Compute and print statistics about the hash table to stdout
 *
 *   Parameters: None
 *   Return value: None
 */
void oahttslf_printstats(void)
{
  unsigned int num_items=0;
  int i;

  for (i = 0; i < OAHTTSLF_SIZE; i++)
  {

    if (hashtable[i].key != OAHTTSLF_EMPTY_KEY)
      num_items++;
  }
  printf("num items       : %u (%f%% full)\n", num_items,
         100.0*(float)num_items/OAHTTSLF_SIZE);
}

/*
 * oahttslf_reset()
 * 
 * reset all the table entries to empty
 *
 * Parameters: None
 * Return value: None
 *
 */
void oahttslf_reset(void)
{
#if defined(USE_INSTRUMENT) || defined(USE_CONTENTION_INSTRUMENT)
  int i;
#endif
  assert(0 == OAHTTSLF_EMPTY_KEY);
  assert(0 == OAHTTSLF_EMPTY_VALUE);
  memset(hashtable, 0, sizeof(hashtable));
#ifdef USE_INSTRUMENT
  for (i = 0; i < MAX_NUM_THREADS; i++)
    key_count[i] = 0;
#endif
#ifdef USE_CONTENTION_INSTRUMENT
  for (i = 0; i < MAX_NUM_THREADS; i++)
    retry_count[i] = 0;
#endif
}



/*
 * oahttslf_num_entries()
 *
 *   count number of keys in the hash table
 *   WARNING: may be very slow - iterates thriough whole table; we do
 *   not have a counter.
 *
 *   Parameters: None
 *   Return value: None
 */
unsigned int oahttslf_num_entries()
{
  unsigned int num_items=0;
  int i;

  for (i = 0; i < OAHTTSLF_SIZE; i++)
  {

    if (hashtable[i].key != OAHTTSLF_EMPTY_KEY)
      num_items++;
  }
  return num_items;
}

#ifdef USE_INSTRUMENT
/*
 *  add up the per-thread key counters and return total
 * Parameters: None
 * Return value: Total number of keys in the hash table
 */
unsigned int oahttslf_total_key_count()
{
  unsigned int num_items = 0;
  int i;
  for (i = 0; i < MAX_NUM_THREADS; i++)
    num_items += key_count[i];
  return num_items;
}
#endif


#ifdef USE_CONTENTION_INSTRUMENT
/*
 *  add up the per-thread retry counters and return total
 * Parameters: None
 * Return value: Total number of times an insertion had to be retried
 */
unsigned int oahttslf_total_retry_count()
{
  unsigned int total_retries = 0;
  int i;
  for (i = 0; i < MAX_NUM_THREADS; i++)
    total_retries += retry_count[i];
  return total_retries;
}
#endif

/*
 * oahttslf_insert_double()
 *
 * Insert a key/value pair into the hashtable, or update the value
 * for existing key.
 *
 * Parameters:
 *    key   - key to insert
 *    value - value to insert for the key
 *    thread_id - id (0,...n, not pthread id) of this thread
 *
 * Return value:
 *    Value for the key prior to the new insertion (OAHTTSLF_EMPTY_VALUE
 *    for a new key)
 */
double oahttslf_insert_double(uint64_t key, double value, int thread_id)
{
  static const char *funcname = "oahttslf_insert_double";
  volatile oahttslf_entry_t *ent;
  uint64_t entkey, oldvalue;
  bool vacant;
  uint64_t newvalue;
  double d_oldvalue;

  assert(key != OAHTTSLF_EMPTY_KEY);
  assert(value != OAHTTSLF_EMPTY_VALUE);

  ent = oahttslf_getent(key, &vacant);
  if (!ent)
  {
    bpa_fatal_error(funcname, "hash table full\n"); /* TODO expand table */
  }
  oldvalue = ent->value;
  entkey = ent->key;
  if (vacant)
  {
    if (CAS64(&ent->key, OAHTTSLF_EMPTY_KEY, key) != OAHTTSLF_EMPTY_KEY)
      return oahttslf_insert_double(key, value, thread_id);  /* tail-recursive call to retry */

    /* in lookup we test for OAHHTSLF_EMPTY_VALUE so if someone looks up 
       in anothe thread before the value is set here, we return key not found.
       NB depends on 64-bit atomic writes */
/*    if (ent->value != OAHTTSLF_EMPTY_VALUE)
      bpa_error_msg(funcname, "value = %llX\n", ent->value);     
      */
/*  __asm__ __volatile__ ("" : : : "memory");  */

    entkey = key ;
  }
#ifdef DEBUG
  /*assert(key == ent->key);*/
  if (key != entkey)
  {
          fprintf(stderr, "OAHTTSLF ASSERTION FAILURE: key=%llX entkey=%llX\n",  key, entkey);
          exit(1);
  }
#endif
  memcpy(&newvalue, &value, sizeof(double));
#ifdef ALLOW_UPDATE
  if (oldvalue == newvalue)  /* shortcut to avoid expense of CAS instruction */
    return newvalue;
  if (CAS64(&ent->value, oldvalue, newvalue) != oldvalue)
    return oahttslf_insert_double(key, value, thread_id);  /* tail-recursive call to retry */
#else
  ent->value = newvalue;
#endif
  memcpy(&d_oldvalue, &oldvalue, sizeof(double));
  return d_oldvalue;
}



/*
 * oahttslf_lookup_double()
 *
 * Get the value for a key from the hashtable
 *
 * Parameters:
 *     key -  key to look up
 *     value - (output) value for key,ony set if TRUE returned.
 * * Return value: *      TRUE if key found, FALSE otherwise.
 */
bool oahttslf_lookup_double(uint64_t key, double *value)
{
  volatile oahttslf_entry_t *ent;
  bool vacant;
  double val;

  ent = oahttslf_getent(key, &vacant);
  if (ent)
  {
    if (!vacant && ent->value != OAHTTSLF_EMPTY_VALUE)
    {
      memcpy(&val, &ent->value, sizeof(double));
      *value = val;
      return TRUE;
    }
  }
  return FALSE;
}




