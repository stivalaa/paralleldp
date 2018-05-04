#ifndef BPAPARSE_H
#define BPAPARSE_H
/*****************************************************************************
 * 
 * File:    bpaparse.h
 * Author:  Alex Stivala
 * Created: June 2007
 *
 * Declarations of structures functions to parse the input files for bpalign.
 * 
 *
 * $Id: bpaparse.h 2584 2009-06-25 07:54:41Z astivala $
 *
 *****************************************************************************/

/*
 * The basepair_t type is a struct representing a single base pair (probability)
 *  as parsed from the input (each line results in one of these structs)
 */

typedef struct basepair_s
{
    int    left;    /* co-ord of first base in pairing */
    int    right;   /* co-ord of second base in pairing */
    double prob;    /* probability of this base pairing */
} basepair_t;

/* read specified input filename and create array of basepair stucts */
basepair_t *bpa_read_basepairs(const char *filename, double pmin,
                               char **sequence,
                               int *bplist_len);

/* debugging function to dump basepairlist to stderr */
void bpa_dump_bp_list(int n, const basepair_t *bplist, const char *sequence);

#endif /* BPAPARSE_H */
