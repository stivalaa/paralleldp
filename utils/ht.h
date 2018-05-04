#ifndef HT_H
#define HT_H
/*****************************************************************************
 * 
 * File:    ht.h
 * Author:  Alex Stivala
 * Created: April 2009
 *
 * Declarations for separate chaining hash table.
 * 
 *
 * $Id: ht.h 3084 2009-12-20 03:39:46Z alexs $
 *
 *****************************************************************************/

#include "bpautils.h"




/* size of the hash table (must be a power of 2)*/
/*#define HT_SIZE  134217728 */   /* 2^27 */ /* too large */
#define HT_SIZE  67108864   /* 2^26 */


typedef struct ht_entry_s
{
    struct ht_entry_s *next;
    char data[1];  /* overlaid with key followed by value, 
                      sizes defined by user. 
                   FIXME: of course this is dodgy because of alignment padding
                   etc., but it seems to work... */
} ht_entry_t;

/* hash function type */
typedef unsigned int (*hash_function_t)(const void *key);

/* key/value copy function type */
typedef void *(*copy_function_t)(void *dest, const void *src);

/* key match function type */
typedef int (*keymatch_function_t)(const void *k1, const void *k2);

/* setup hashtable key and value types and functions */
void ht_initialize(size_t keysize, size_t valuesize, hash_function_t hashfunc,
                   copy_function_t keycopy, keymatch_function_t keymatch,
                   copy_function_t valuecopy);

/* insert into hashtable */
void ht_insert(void *key, void *value);

/* lookup in hashtable */
void *ht_lookup(void *key);

/* test for invalid structure */
int ht_validate(void);

/* compute and print stats about hash table */
void ht_printstats(void);

#endif /* HT_H */
