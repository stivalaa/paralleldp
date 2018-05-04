/*************************************************************************
 *
 * timeguard - run a process with an elapsed time limit
 *
 * File:    timeguard.
 * Author:  Alex Stivala
 * Created: December 2009
 *
 * Usage:
 *     timeguard timeout comand_and_args
 *
 *     timeout is timeout value in seconds
 *
 *     command_and_args is the command and its arugment list to run.
 *
 * Run the specified command, and terminate it with SIGALRM if it has
 * not already completed after timeout elapsed time.  This script is
 * needed because ulimit can only limit CPU time not elapsed time, so
 * cannot guard against processes sleeping or in some sort of deadlock
 * or I/O wait (ps state D) etc.
 *
 * $Id: timeguard.c 3046 2009-12-15 06:32:13Z alexs $
 * 
 *************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>

static void usage(const char *progname)
{
  fprintf(stderr, "Usage: %s timeout command_and_args\n",progname);
  exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
  unsigned int timeout_seconds;
  char **new_argv;

  if (argc < 2)
    usage(argv[0]);

  timeout_seconds = (unsigned)atoi(argv[1]);
  if (timeout_seconds < 1)
    usage(argv[0]);
  
  new_argv = (char **)malloc(argc * sizeof(char *));
  new_argv[0] = argv[2];
  memcpy(&new_argv[1], &argv[3], (argc - 1) * sizeof(char *));

  alarm(timeout_seconds);

  if (execvp(new_argv[0], new_argv) < 0)
    perror("execvp failed");

  return 0; /* just to avoid warning: we never reach here because of execvp() */
}

