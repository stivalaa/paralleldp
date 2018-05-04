/*****************************************************************************
 * 
 * File:    bpaparse.c
 * Author:  Alex Stivala
 * Created: June 2007
 *
 * $Id: bpaparse.c 2584 2009-06-25 07:54:41Z astivala $
 *
 * Functions to parse the input files for bpalign. 
 * Instead of having to write C code to parse the PostScript dotplot output
 * for Vienna RNA package RNAfold -p, we instead parse a very simple format
 * that is generated from the former by the rnafold2list.py script (Python
 * being a much more appropriate language for parsing this than C, and
 * conversely as it turns out C being much more appropriate for the computation
 * than Python).
 * 
 *
 *   The first line is two integers separeated by whitespace
 *
 *   n m
 *
 *   where n is the length of the sequence and m is the number of base pairs
 *   listed.
 *
 *   The second line is the sequence as a string of A,U,C,G charcters
 *   (terminated by newline). Regardless of sequence length it is all on 
 *   one line (of length n)
 *
 *   Each subsequent line (of which there are m) is of the form:
 *
 *   i j p
 *
 *   where i and j are integers denoting two base pairs (starting from 0)
 *   that are paired with probability p. p is floating point number in printf
 *   %d format.
 *
 *   Lines starting with '#' in column 1 are comments and are ignored.
 *
 * See also rnafold2list.py and readrnafold.py.
 *
 *****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "bpaglobals.h"
#include "bpautils.h"
#include "bpaparse.h"

#define COMMENT_CHAR '#'   /* comments in input file have this in column 1 */
#define MAX_LINE_LEN 8192  /* maximum input line length NB does NOT apply to the
                              line with the sequence on it! */


/*
 * bpa_read_basepairs()
 *
 * read specified input filename and create array of basepair stucts 
 * and also the sequence. Discard those with probability < pmin.
 *
 * Parameters:
 *  filename      - input filename (containing output of rnafold2list.py)
 *  pmin          - minimum base pairing probability; discard those with
 *                  probability < pmin.
 *  sequence      - (OUT) newly allocated string containg RNA sequence
 *  bplist_len    - (OUT) number of elements in bplist
 *
 * Return value:
 *   Newly allocated array of basepair_t structs, one for each of the 
 *   i,j,p lines in the input, or
 *   NULL on error.
 *
 */
basepair_t *bpa_read_basepairs(const char *filename, double pmin,
                               char **sequence,
                               int *bplist_len)
{
  static const char *funcname = "bpa_read_basepairs";
  char linebuf[MAX_LINE_LEN];
  basepair_t *bplist; /* return value: newly allocated list of basepairs */
  FILE *fp;
  int bp_index;
  int seqlen;         /* length of sequence */
  int num_basepairs;  /* number of basepairs */
  int lineno;         /* current line number in input */

  if ((fp = fopen(filename, "r")) == NULL)
  {
    bpa_error_msg(funcname, "Can't open input file %s for reading (%s)\n", 
                  filename,
                  strerror(errno));
    return NULL;
  }
  lineno = 1;
  while (fgets(linebuf, MAX_LINE_LEN, fp) && linebuf[0] == COMMENT_CHAR)
    lineno++; /* skip comment lines */

  /* first line has two ints, length of seq and length of bp list */

  if (sscanf(linebuf, "%d %d\n", &seqlen, &num_basepairs) != 2)
  {
    bpa_error_msg(funcname, 
                  "Parse error in file %s at line %d: expecting 2 ints\n",
                  filename, lineno);
    return NULL;
  }
  if (seqlen < 0 || num_basepairs < 0)
  {
    bpa_error_msg(funcname, 
                  "Parse error in file %s at line %d: expecting 2 +ve ints\n",
                  filename, lineno);
    return NULL;
  }
  lineno++;

  /* TODO: there is a bug here that can't have comment between two ints line
     and sequence, probably should fix it */

  /* second line is the sequence, all on one line */
  *sequence = (char*)bpa_malloc((seqlen+2) * sizeof(char)); /* add 2 to allow newline*/
  if (fgets(*sequence, seqlen+2, fp) == NULL) /* fgets includes the newline */
  {
    bpa_error_msg(funcname, 
                  "Parse error in file %s at line %d\n", filename, lineno);
    free(*sequence);
    return NULL;
  }
  (*sequence)[strlen(*sequence)-1] = '\0'; /* overwrite the newline with nul */
  
  if (strlen(*sequence) != (unsigned)seqlen)
  {
    bpa_error_msg(funcname, 
    "Parse error in file %s at line %d: expecting sequence length %d got %d\n",
                  filename, lineno, seqlen, strlen(*sequence));
    free(*sequence);
    return NULL;
  }
  lineno++;

  /* 
   * subsequent lines are 
   *  i j p
   * base pair probabilities
   */
  bplist = (basepair_t*)bpa_malloc(num_basepairs * sizeof(basepair_t));
  bp_index = 0;
  while (fgets(linebuf, MAX_LINE_LEN, fp))
  {
    if (linebuf[0] == COMMENT_CHAR)
    {
      lineno++; 
      continue; /* skip comment */
    }
    if (bp_index >= num_basepairs)
    {
      bpa_error_msg(funcname, "warning: not reading %s after line %d\n",
                    filename, lineno);
      break;
    }
    if (sscanf(linebuf, "%d %d %lf\n", &bplist[bp_index].left, 
               &bplist[bp_index].right,
               &bplist[bp_index].prob) != 3)
    {
      bpa_error_msg(funcname, 
                    "Parse error in file %s at line %d\n",
                    filename, lineno);
      free(*sequence);
      free(bplist);
      return NULL;
    }
    lineno++;
    if (bplist[bp_index].prob >= pmin)
      bp_index++; /* only keep entry if prob >= pmin */
  }
  *bplist_len = bp_index;
  return bplist;
}


/*
 * bpa_dump_bp_list()
 *
 * Debugging function to dump basepair list and sequence to stderr.
 *
 * Parmeters:
 *   n      - length of bplist
 *   bplist - list of basepair structs
 *   sequence - also set by bpa_read_basepairs()
 */
void bpa_dump_bp_list(int n, const basepair_t *bplist, const char *sequence)
{
  static const char *funcname = "dump_bp_list";
  int i;

  if (!sequence)
    fprintf(stderr, "%s null sequence\n", funcname);
  else
    fprintf(stderr, "%s\n", sequence);
        
  if (!bplist)
    fprintf(stderr, "%s null bplist\n", funcname);
  else
  {
    for (i = 0; i < n; i++)
      fprintf(stderr, "%d %d %f\n", bplist[i].left, 
              bplist[i].right, bplist[i].prob);
  }
}
