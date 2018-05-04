/*****************************************************************************
 * 
 * File:    numcores.c
 * Author:  Alex Stivala
 * Created: April 2009
 *
 * Write the number of processors online to stdout.
 * 
 *
 * $Id: numcores.c 2316 2009-05-05 08:38:59Z astivala $
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>

int main(int argc, char *argv[])
{
  long ncores;
  if (argc != 1) 
  {
    fprintf(stderr, "usage: %s\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  if ((ncores =  sysconf(_SC_NPROCESSORS_ONLN)) < 0)
  {
    fprintf(stderr, "sysconf() failed (%d)\n", errno);
    exit(EXIT_FAILURE);
  }
  printf("%ld\n", ncores);
  exit(EXIT_SUCCESS);
}
    
