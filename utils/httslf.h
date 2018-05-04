#ifndef HTTSLF_H
#define HTTSLF_H
/*****************************************************************************
 * 
 * File:    httslf.h
 * Author:  Alex Stivala
 * Created: April 2009
 *
 * Declarations for separate chaining thread-safe lock-free hash table.
 * 
 *
 * $Id: httslf.h 2519 2009-06-14 04:12:41Z astivala $
 *
 *****************************************************************************/

#include "bpautils.h"


/* size of the hash table (must be a power of 2)*/
/*#define HTTSLF_SIZE  134217728 */   /* 2^27 */ /* too large */
#define HTTSLF_SIZE  67108864   /* 2^26 */

typedef struct httslf_entry_s
{
    struct httslf_entry_s *next;
    char data[1];  /* overlaid with key followed by value, 
                      sizes defined by user. 
                   FIXME: of course this is dodgy because of alignment padding
                   etc., but it seems to work... */
} httslf_entry_t;

/* hash function type */
typedef unsigned int (*hash_function_t)(const void *key);

/* key/value copy function type */
typedef void *(*copy_function_t)(void *dest, const void *src);

/* key match function type */
typedef int (*keymatch_function_t)(const void *k1, const void *k2);

/* setup hashtable key and value types and functions */
void httslf_initialize(size_t keysize, size_t valuesize, hash_function_t hashfunc,
                   copy_function_t keycopy, keymatch_function_t keymatch,
                   copy_function_t valuecopy);

/* insert into hashtable */
httslf_entry_t *httslf_insert(void *key, void *value);

/* lookup in hashtable */
void *httslf_lookup(void *key);

/* test for invalid structure */
int httslf_validate(void);

/* compute and print stats about hash table */
void httslf_printstats(void);

#endif /* HTTSLF_H */
