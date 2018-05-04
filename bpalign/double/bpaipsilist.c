/*****************************************************************************
 * 
 * File:    bpaipsilist.c
 * Author:  Alex Stivala
 * Created: June 2007
 *
 *
 * $Id: bpaipsilist.c 2584 2009-06-25 07:54:41Z astivala $
 *
 * ipsilist implementation.
 *
 * An indexed psi list (ipsilist) is a table indexed by i
 * where ipsilist[i] is a list of (j,psi) tuples where
 * psi = log(P[i,j] / Pmin)
 * where P[i,j] is the prob at i,j in the list and Pmin is the minimum
 * probability deemed signficant (global constant PMIN).
 * For values i,j,psi with no probability or probability less than PMIN
 * there is no entry in the list at i.
 *
 *****************************************************************************/

#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "bpaglobals.h"
#include "bpautils.h"
#include "bpaipsilist.h"

/*****************************************************************************
 *
 * external functions
 *
 *****************************************************************************/

#define INITIAL_LIST_LEN 1 /* allocate this many ipsi elements initially */


/*
 * add_ipsi_element()
 *
 * add a new ipsi element on the end of the list
 * Note we memcpy the ipsi element onto the array (allocating or
 * reallocating if necessary),
 * this is not a linked list but an array - the element is not added to 
 * the list itself, it is copied.
 * Using memcpy assumes of course no pointers in ipsi element. 
 *
 * Parameters:
 *     ipsilist - ipsilist to add to
 *     ipsi     - an ipsi element to and add on to end of list
 *
 * Return value:
 *     None
 */
void bpa_add_ipsi_element(ipsi_list_t *ipsilist, const ipsi_element_t *ipsi)
{
  assert(ipsilist->num_elements <= ipsilist->allocated_size);

  if (ipsilist->allocated_size == 0)
  {
    /* nothing allocated yet, make initial allocation */
    ipsilist->ipsi = (ipsi_element_t*)bpa_calloc(INITIAL_LIST_LEN, sizeof(ipsi_element_t));
    ipsilist->allocated_size = INITIAL_LIST_LEN;
  }
  else if (ipsilist->num_elements == ipsilist->allocated_size)
  {
    /* realloc for one more element  */
    ipsilist->ipsi = (ipsi_element_t*)bpa_realloc(ipsilist->ipsi, 
                                 ++ipsilist->allocated_size *
                                 sizeof(ipsi_element_t));
  }
  /* else num_elements < allocated_size (see assertion above) */

  memcpy(&ipsilist->ipsi[ipsilist->num_elements++], ipsi, sizeof(ipsi_element_t));
}


/*
 * bpa_pairlist_to_ipsilist()
 *
 * Convert the base pairing probability list
 * i.e.  the list of tuples (i, j, p)
 * where i,j are (zero-based) co-ordinates in the sequence
 * and p is the the base pair probabillity at i,j
 * (assumed to be 0 at i,j when i,j not in list) into a
 * an indexed psi list, as described above.
 *
 * Note it is a precondition that the list passed in already has elements
 * with p < PMIN filtered out. bpa_read_basepairs() does this.
 *
 * It is required that the ipsilist[i]s are sorted by right co-ord
 * ascending. This is assured by the input pairlists generated
 * by rnafold2list.py being sorted by i ascending and within i by j ascending.
 *
 * Parameters:
 *    pairlist      - pairlist (array of basepair_t) from bpa_read_basepairs()
 *    pairlist_len  - number of elements in pairlist
 *    seq_len       - length of sequence (this will be length of array returned)
 *
 * Return value:
 *    Newly allocated array of ipsi_list_s structs, each
 *    of which has a pointer to a newly allocated array of
 *    of ipsi_list_elements, being the (j,psi) values for that i.
 * 
 */
ipsi_list_t *bpa_pairlist_to_ipsilist(basepair_t *pairlist, 
                                      int pairlist_len,
                                      int seq_len)
{
  const double   inv_log_invpmin = 1 / log(1 / PMIN); /* normalization factor */
  ipsi_list_t    *ipsilist       = NULL;              /* return value */
  ipsi_element_t ipsi_element; /* on stack as copied to end of list not linked*/
  int            i;
  
  /* allocate ipsilist header and actual list with one element for each base */
  ipsilist = (ipsi_list_t*)bpa_calloc(seq_len, sizeof(ipsi_list_t)); /* all fields zeroed */
  

  for (i = 0; i < pairlist_len; i++)
  {
    assert(pairlist[i].left < pairlist[i].right);
    assert(pairlist[i].prob >= PMIN);
    
    ipsi_element.right = pairlist[i].right;
    ipsi_element.psi = log(pairlist[i].prob / PMIN) * inv_log_invpmin;
    ipsi_element.arclen_diff = -1; /* only used if useordering in bpadynprog*/
    assert(ipsi_element.psi >=  0);  /* can be 0 if prob == PMIN */
    bpa_add_ipsi_element(&ipsilist[pairlist[i].left], &ipsi_element);
    
  }
  return ipsilist;
}


/*
 * bpa_free_ipsilist()
 *
 * free ipsilist allocated by bpa_pairlist_to_ipsilist 
 * Frees the ipsilist and all memory allocted under it
 *
 * Parameters:
 *     ipsilist - ipsilist to free
 *     list_len - length of the array ipsilist
 *
 * Return value:
 *     None.
 */
void bpa_free_ipsilist(ipsi_list_t *ipsilist, int list_len)
{
  int i;

  for (i = 0; i < list_len; i++)
    free(ipsilist[i].ipsi);
  free(ipsilist);
}

/*
 * bpa_dump_ipsilist()
 *
 * debugging function to dump ipsilist to stder 
 *
 * Parameters:
 *    ipsilist - ipsilist (array of ipsislist_t structs) to print
 *    list_len - number of elements in ipsilist array
 */
void bpa_dump_ipsilist(const ipsi_list_t *ipsilist, int list_len)
{
  int i,j;
  
  for (i = 0; i < list_len; i++)
  {
    fprintf(stderr, "%d ", i);
    if (ipsilist[i].num_elements > 0)
    {
      for (j = 0; j <ipsilist[i].num_elements; j++)
      {
        fprintf(stderr, "(%d, %g) ", ipsilist[i].ipsi[j].right,
                ipsilist[i].ipsi[j].psi);
      }
      fprintf(stderr, "\n");
    }
  }
}

/*
 * bpa_serialize_ipsilist()
 *
 * convert the ipsilist to a format contiguous in memory with no pointers
 * This way it can be sent to the GPU with CUDA.
 * We will not worry about memory efficiency, but make it efficient
 * to access without having to deserialize it: it will simply be a
 * 2D array of ipisi_elements, with the leading dimension the maximum
 * number of elements in any row. Unused elements have all fields 0.
 *
 * Probably should do something better like use packed column
 * (Harwell-Boeing) sparse matrix format for psi matrices. (TODO)
 *
 * Parameters:
 *    ipsilist - ipsilist (array of ipsislist_t structs) to print
 *    list_len - number of elements in ipsilist array
 *    n        - (output) leading dimension of returned array (number of
 *               entries in each of the list_len rows)
 *
 * Return value:
 *    Pointer to serialzied ipsilist, newly allocated.
 */
ipsi_element_t *bpa_serialize_ipsilist(const ipsi_list_t *ipsilist, int list_len,
                                    int *n)
{
  int i,j,max_num_elements;
  ipsi_element_t *serial_ipsilist;
  
  max_num_elements = 0;
  for (i = 0; i < list_len; i++)
  {
    if (ipsilist[i].num_elements > max_num_elements)
      max_num_elements = ipsilist[i].num_elements;
  }        
  serial_ipsilist = (ipsi_element_t *)bpa_calloc(list_len * max_num_elements,
                                                 sizeof(ipsi_element_t));
  for (i = 0; i < list_len; i++)
  {
    for (j = 0; j <ipsilist[i].num_elements; j++)
    {
      serial_ipsilist[i * max_num_elements + j].right = ipsilist[i].ipsi[j].right;
      serial_ipsilist[i * max_num_elements + j].psi = ipsilist[i].ipsi[j].psi;
      serial_ipsilist[i * max_num_elements + j].arclen_diff = ipsilist[i].ipsi[j].arclen_diff;
      
    }
  }  
  *n = max_num_elements;
  return serial_ipsilist;
}

/*
 * bpa_dump_seripsilist()
 *
 * debugging function to dump serialized ipsilist to stderr
 *
 * Parameters:
 *    seripsilist - serizliaed ipsilist to print
 *    list_len - number of rows in serizliaed ipsilist
 *    n        - leading dimension of seripsilist array (number of
 *               entries in each of the list_len rows)
 */
void bpa_dump_seripsilist(const ipsi_element_t *seripsilist, int list_len,
        int n)
{
  int i,j;
  
  for (i = 0; i < list_len; i++)
  {
    fprintf(stderr, "%d ", i);
    for (j = 0; j < n; j++)
    {
      if (seripsilist[i*n+j].right == 0)
        break;
      fprintf(stderr, "(%d, %g) ", seripsilist[i*n + j].right,
              seripsilist[i*n + j].psi);
    }
    fprintf(stderr, "\n");
  }
}

