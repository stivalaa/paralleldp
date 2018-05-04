#ifndef BPAIPSILIST_H
#define BPAIPSILIST_H
/*****************************************************************************
 * 
 * File:    bpaipsilist.h
 * Author:  Alex Stivala
 * Created: June 2007
 *
 *
 * $Id: bpaipsilist.h 2584 2009-06-25 07:54:41Z astivala $
 *
 * Declarations of structures and functions for the ipsilist.
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

#include "bpaparse.h"

/*
 * ipsi_element_t is a single element in the ipsilist.
 */
typedef struct ipsi_element_s
{
    int    right;   /* co-ord of 2nd base in pairing (1st is index in list) */
    double psi;     /* psi value of this pairing */
    int    arclen_diff; /* difference in arc len used by -o (ordering) opt */
} ipsi_element_t;


/* 
 * the overall ipsilist is an array of these ipsi_list_t lists (one
 * for each i i.e. each position in sequence
 */
typedef struct ipsi_list_s
{
    int allocated_size;   /* number of elements currently allocated for */
    int num_elements;     /* number of elements actually in list */
    ipsi_element_t *ipsi; /* array of ipsi_element_t structs */
} ipsi_list_t;


/* Convert basepair list from bpaparse to ipsilist */
ipsi_list_t *bpa_pairlist_to_ipsilist(basepair_t *pairlist, 
                                      int pairlist_len, int seq_len);

/* free ipsilist allocated by bpa_pairlist_to_ipsilist */
void bpa_free_ipsilist(ipsi_list_t *ipsilist, int list_len);

/* add a new ipsi element onto the end of the list */
void bpa_add_ipsi_element(ipsi_list_t *ipsilist, const ipsi_element_t *ipsi);

/* debugging function to dump ipsilist to stder */
void bpa_dump_ipsilist(const ipsi_list_t *ipsilist, int list_len);

/* serialize the ipsilist into contiguous memory block */
ipsi_element_t *bpa_serialize_ipsilist(const ipsi_list_t *ipsilist, int list_len,
                                    int *n);

void bpa_dump_seripsilist(const ipsi_element_t *seripsilist, int list_len,
        int n);

#endif /* BPAIPSILIST_H */
