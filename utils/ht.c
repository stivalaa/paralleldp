/*****************************************************************************
 * 
 * File:    ht.c
 * Author:  Alex Stivala
 * Created: April 2009
 *
 * Separate chaining hash table.
 * 
 *
 * $Id: ht.c 3084 2009-12-20 03:39:46Z alexs $
 *
 *****************************************************************************/

#include <stdio.h>
#include <string.h>

#include "bpautils.h"
#include "ht.h"

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
/* TODO: change to allocate dynamically so can have multiple */
static ht_entry_t *hashtable[HT_SIZE];

/* user data sizes and callback functions set by ht_initialize() */
static size_t key_size;             /* size of key data */
static size_t value_size;           /* size of value data */
hash_function_t hash_function;      /* hash function */
keymatch_function_t keymatch_function;  /* key match function (compare keys) */
copy_function_t keycopy_function;   /* copy key data (use memcpy if NULL) */
copy_function_t valuecopy_function; /* copy value data (use memcpy if NULL) */


/*****************************************************************************
 *
 * external functions
 *
 *****************************************************************************/

/*
 * ht_insert()
 *
 * Insert a key/value pair into the hashtable
 * NB This only allows insertion of a NEW key - if the key already
 * exists, it is considered a fatal error.
 * This is for use in dynamic programming 
 * simple case (no bounding) where a particular key once its value is set
 * is the optimal - we should never try to set for the same key more than once.
 *
 *
 * Parameters:
 *    key   - ptr to key to insert
 *    value - value to insert for the key
 *
 * Return value:
 *    None.
 */
void ht_insert(void *key, void *value)
{
  unsigned int h;
  ht_entry_t *ent;

  h = hash_function(key);
  ent = hashtable[h];
  while (ent && !keymatch_function(key, &ent->data))
    ent = ent->next;
  if (!ent)
  {
#ifdef USE_CP_ALLOC
    ent = (ht_entry_t *)cellpool_alloc();
#else
    ent = (ht_entry_t *)bpa_malloc(sizeof(ht_entry_t) + key_size + value_size);
#endif
    if (keycopy_function)
      keycopy_function(&ent->data, key);
    else
      memcpy(&ent->data, key, key_size);
    if (valuecopy_function)
      valuecopy_function(&ent->data + key_size, value);
    else
      memcpy(&ent->data + key_size, value, value_size);
    /* insert at head of list */
    ent->next = hashtable[h];
    hashtable[h] = ent;
  }
  else
  {
    bpa_fatal_error("ht_insert", "key already set");
  }
}



/*
 * ht_lookup()
 *
 * Get the value for a key from the hashtable
 *
 * Parameters:
 *     key - ptr to key to look up
 * 
 * Return value:
 *     pointer to value for the key, NULL if not found.
 */
void *ht_lookup(void *key)
{
  unsigned int h;
  ht_entry_t *ent;
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
 * ht_initialize()
 *
 * setup hashtable key and value types and functions 
 *
 * Parameters:
 *    keysize     - size of key data
 *    valuesize   - size of value data
 *    hashfunc    - hash function: take ptr to key and return hash value
 *    keycopy     - key copy function: copy key data from second param ptr
 *                  to first param ptr, return first param ptr
 *                  May be NULL, then  memcpy() is used.
 *    keymatch    - key match function: given two key data ptrs, return
 *                  0 iff keys match (are equal).
 *    valuecopy   - value copy function, similar to keycopy
 *                  May be NULL, then memcpy() is used.
 *
 * Return value:
 *    None.
 *    
 */
void ht_initialize(size_t keysize, size_t valuesize, hash_function_t hashfunc,
                   copy_function_t keycopy, keymatch_function_t keymatch,
                   copy_function_t valuecopy)
{
  static const char *funcname = "ht_initialize";
  size_t cell_size;

  key_size = keysize;
  value_size = valuesize;
  hash_function = hashfunc;
  keycopy_function = keycopy;
  keymatch_function = keymatch;
  valuecopy_function = valuecopy;
#ifdef USE_CP_ALLOC
  cell_size = sizeof(ht_entry_t) + key_size + value_size;
  /* FIXME need to pass number of cells as a parameter in this function */
  if (!(cellpool_initialize(cell_size, (size_t)1*1024*1024*1024 / cell_size)))
    bpa_fatal_error(funcname, "cellpool_initialize() failed\n");
#endif
}




/*
 * ht_validate()
 *
 * Test for duplicate keys in the lists -this should not happen
 *
 * Parameters:
 *    None
 *
 * Return value:
 *    0 if duplicate keys found else 1
 */
int ht_validate(void)
{
  ht_entry_t *ent1, *ent2;
  int i;

  for (i = 0; i < HT_SIZE; i++)
      for (ent1 = hashtable[i]; ent1 != NULL; ent1 = ent1->next)
          for (ent2 = ent1->next; ent2 != NULL; ent2 = ent2->next)
              if (keymatch_function(&ent1->data, &ent2->data)) {
                  return 0;
              }
  return 1;
}


/*
 * ht_printstats()
 *
 *   Compute and print statistics about the hash table to stdout
 *
 *   Parameters: None
 *   Return value: None
 */
void ht_printstats(void)
{
  unsigned int num_items=0, num_entries=0;
  unsigned int chain_length,max_chain_length=0,sum_chain_length=0;
  float avg_chain_length;
  int i;
  ht_entry_t *ent;

  for (i = 0; i < HT_SIZE; i++)
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
         100.0*(float)num_items/HT_SIZE);
  printf("max chain length: %u\n", max_chain_length);
  printf("avg chain length: %f\n", avg_chain_length);
}
