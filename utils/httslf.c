/*****************************************************************************
 * 
 * File:    httslf.c
 * Author:  Alex Stivala
 * Created: April 2009
 *
 * Separate chaining thread-safe lock-free hash table.
 * 
 * gcc version 4.1.0 or greater is required, in order to use the
 * __sync_val_compare_and_swap()
 * and __sync_val_fetch_and_add() builtins
 * (this module was developed on Linux 2.6.22 (x86) with gcc 4.1.3),
 * except on Solaris, where we use SUNWSpro compiler and
 * cas32()/caslong()/cas64()/casptr() in 
 * /usr/include/sys/atomic.h)
 *
 * $Id: httslf.c 2324 2009-05-06 02:44:04Z astivala $
 *
 *****************************************************************************/

#include <stdio.h>
#include <string.h>

#include "bpautils.h"
#include "httslf.h"
#include "cellpool.h"
#include "atomicdefs.h"

/* do not use cell-pool allocator on Linux (malloc() is faster) ,
  but cell-pool allocator is faster on Solaris (SPARC) */
#ifdef SOLARIS
#define USE_CP_ALLOC 
#endif

/*****************************************************************************
 *
 * static data
 *
 *****************************************************************************/


/* The hash table */
/* Each entry is head of chain pointer and
   is serialized with CAS logic in httslf_insert() */
/* TODO: change to allocate dynamically so can have multiple */
static httslf_entry_t *hashtable[HTTSLF_SIZE];

/* user data sizes and callback functions set by ht_initialize() */
static size_t key_size;             /* size of key data */
static size_t value_size;           /* size of value data */
hash_function_t hash_function;      /* hash function */
keymatch_function_t keymatch_function;  /* key match function (compare keys) */
copy_function_t keycopy_function;   /* copy key data (use memcpy if NULL) */
copy_function_t valuecopy_function; /* copy value data (use memcpy if NULL) */

 /* counters for hash collisions: serialize with __sync_fetch_and_add() */
static unsigned int insert_collision_count = 0; /* TODO implement this */
static unsigned int lookup_collision_count = 0; /* TODO implement this */


/*****************************************************************************
 *
 * external functions
 *
 *****************************************************************************/

/*
 * httslf_insert()
 *
 * Insert a key/value pair into the hashtable
 * NB This only allows insertion of a NEW key - if the key already
 * exists, we do nothing (and do NOT update the value if it is different).
 * This is for use in dynamic programming with multiple threads
 * simple case (no bounding) where a particular key once its value is set
 * is the optimal - any other thread can only ever compute the same value
 * anyway. The case where we have bounds and values can change is more
 * difficult, and not handled by this function.
 *
 * Parameters:
 *    key   - ptr to key to insert
 *    value - value to insert for the key
 *
 * Return value:
 *    Pointer to entry inserted.
 */
httslf_entry_t *httslf_insert(void *key, void *value)
{
  static const char *funcname = "httslf_insert";
  unsigned int h;
  httslf_entry_t *ent, *oldent, *newent = NULL;
  httslf_entry_t *inserted_entry = NULL;

  h = hash_function(key);
  do
  {
    ent = hashtable[h];
    oldent = ent; /* old value for CAS logic */
    while (ent && !keymatch_function(key, &ent->data))
      ent = ent->next;
    if (!ent)
    {
      /* bpa_malloc() is just malloc(). When compiling & linking
         with -pthread (wtih gcc at least), malloc() is threadsafe but
         it is not lock-free so there is a cost to using this, a
         lock-free malloc() would be nice.  There is one in nbds
         (http://code.google.com/p/nbds/) [could also just use that
         library instead of using this at all I suppose anyway -
         using separate chaining like this is inefficient in many ways,
         the closed hash table as in nbds does no malloc() at all
         on insertions so this problem goes away completely,
         and it is faster and works better with cache]; see also
         Michael 2004 "Scalable Lock-Free Dynamic Memory Allocation"
         PLDI'04 [note reference to CAS first appearing in S/370 POP!]).
         Although
         http://developers.sun.com/solaris/articles/multiproc/multiproc.html
         seems to show that the GNU libc malloc (ptmalloc) performs
         pretty well in multithreaded code anyway so maybe this is fine
         as it is (at least when using gcc on Solaris).
         But we can't really claim the hashtable is lockfree if the
         malloc() it uses is not...
         so cp_alloc() uses a trivial cell pool allocator to just
         get a new cell from a preallocated pool, completely lock-free.
      */
      if (!newent)
      {
#ifdef USE_CP_ALLOC
        newent = (httslf_entry_t *)cellpool_alloc();
        if (!newent)
          bpa_fatal_error(funcname, "cellpool_alloc() failed\n");
#else
        newent = (httslf_entry_t *)bpa_malloc(sizeof(httslf_entry_t) +
                                              key_size + value_size);
#endif
        if (keycopy_function)
          keycopy_function(&newent->data, key);
        else
          memcpy(&newent->data, key, key_size);
        if (valuecopy_function)
          valuecopy_function(&newent->data + key_size, value);
        else
          memcpy(&newent->data + key_size, value, value_size);

        /* insert at head of list using CAS instruction - if we lose
           the race (another thread is inserting this key also) then
           we re-try (in do while loop).
        */
      }
      newent->next = oldent;
      inserted_entry = newent;
    }
    else
    {
      /* key already exists, just ignore the new one
         NB we do NOT update the value here, see header comment */
      inserted_entry = ent;
      break;
    }
  }
  while (CASPTR(&hashtable[h], oldent, newent) != oldent);
  
  return inserted_entry;
}



/*
 * httslf_lookup()
 *
 * Get the value for a key from the hashtable
 *
 * Parameters:
 *     key - ptr to key to look up
 * 
 * Return value:
 *      pointer to value or NULL if not present
 */
void *httslf_lookup(void *key)
{
  unsigned int h;
  httslf_entry_t *ent;
  h = hash_function(key);
  ent = hashtable[h];
  while (ent && !keymatch_function(key, &ent->data))
    ent = ent->next;
  if (ent)
    return &ent->data + key_size;
  else
    return NULL;
}





/* 
 * httslf_initialize()
 *
 * setup hashtable key and value types and functions, and initialize
 * the cell pool allocator for the size of the entries.
 *
 * Parameters:
 *    keysize     - size of key data
 *    valuesize   - size of value data
 *    hashfunc    - hash function: take ptr to key and return hash value
 *    keycopy     - key copy function: copy key data from second param ptr
 *                  to first param ptr, return first param ptr
 *                  May be NULL, then  memcpy() is used.
 *    keymatch    - key match function: given two key data ptrs, return
 *                  nonzero iff keys match (are equal).
 *    valuecopy   - value copy function, similar to keycopy
 *                  May be NULL, then memcpy() is used.
 *
 * Return value:
 *    None.
 *    
 */
void httslf_initialize(size_t keysize, size_t valuesize, hash_function_t hashfunc,
                   copy_function_t keycopy, keymatch_function_t keymatch,
                   copy_function_t valuecopy)
{
  static const char *funcname = "httslf_initialize";
  size_t cell_size;

  key_size = keysize;
  value_size = valuesize;
  hash_function = hashfunc;
  keycopy_function = keycopy;
  keymatch_function = keymatch;
  valuecopy_function = valuecopy;
#ifdef USE_CP_ALLOC
  cell_size = sizeof(httslf_entry_t) + key_size + value_size;
  /* FIXME need to pass number of cells as a parameter in this function */
  if (!(cellpool_initialize(cell_size, (size_t)1*1024*1024*1024 / cell_size)))
      bpa_fatal_error(funcname, "cellpool_initialize() failed\n");
#endif
}



/*
 * httslf_validate()
 *
 * Test for duplicate keys in the lists -this should not happen
 *
 * Parameters:
 *    None
 *
 * Return value:
 *    0 if duplicate keys found else 1
 */
int httslf_validate(void)
{
  httslf_entry_t *ent1, *ent2;
  int i;

  for (i = 0; i < HTTSLF_SIZE; i++)
      for (ent1 = hashtable[i]; ent1 != NULL; ent1 = ent1->next)
          for (ent2 = ent1->next; ent2 != NULL; ent2 = ent2->next)
              if (keymatch_function(&ent1->data, &ent2->data))
              {
                  return 0;
              }
  return 1;
}

/* TODO function to free hash table entries */


/*
 * httslf_printstats()
 *
 *   Compute and print statistics about the hash table to stdout
 *
 *   Parameters: None
 *   Return value: None
 */
void httslf_printstats(void)
{
  unsigned int num_items=0, num_entries=0;
  unsigned int chain_length,max_chain_length=0,sum_chain_length=0;
  float avg_chain_length;
  int i;
  httslf_entry_t *ent;

  for (i = 0; i < HTTSLF_SIZE; i++)
  {
    chain_length = 0;
    if ((ent = hashtable[i]) != NULL)
      num_entries++;
    while (ent)
    {
      num_items++;
      chain_length++;
      ent = ent->next;
    }
    sum_chain_length += chain_length;
    if (chain_length > max_chain_length)
      max_chain_length = chain_length;
  }
  avg_chain_length = (float)sum_chain_length / num_entries;
  printf("num slots used  : %u\n", num_entries);
  printf("num items       : %u (%f%% full)\n", num_items,
         100.0*(float)num_items/HTTSLF_SIZE);
  printf("max chain length: %u\n", max_chain_length);
  printf("avg chain length: %f\n", avg_chain_length);
}
