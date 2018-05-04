###############################################################################
#
# File      : Makefile
# Author    : Alex Stivala (astivala)
# Created   : June 2007
#
# $Id: Makefile 2584 2009-06-25 07:54:41Z astivala $
#
# Makefile for building the the C 'bpalign' program.
#
#
# The default target "all" will compile and link the program.
# The "clean" target will remove object and executable files.
#
#
# set MODE=DEBUG to build with debugging and profiling on. Otherwise
# default is to build with optimizations on and no debug or profile.
#
# This Makefile uses implicit rules to compile .c files into .o files
# and some features specific to GNU Makee.
#
###############################################################################


# include the common.mk file from parent diretory, and override it with
# local.mk if it exists

include ../../common.mk
-include ../../local.mk

INCDIRS =  -I../../utils
NBDSINC =  -I../../nbds.0.4.3/include

COMMONSRCS  = bpaglobals.c bpaparse.c bpaipsilist.c parbpamain.c \
              bpadynprog_single.c
OAHTTSLFSRCS  = bpadynprog_oahttslf.c bpadynprog_threadcall.c
NBDSSRCS    = bpadynprog_nbds.c
RANDSRCS    = bpadynprog_rand_oahttslf.c

SRCS = $(COMMONSRCS) $(OAHTTSLFSRCS) $(NBDSSRCS)

COMMONOBJS    = $(COMMONSRCS:.c=.o)
OAHTTSLFOBJS  = $(OAHTTSLFSRCS:.c=.o)
NBDSOBJS      = $(NBDSSRCS:.c=.o)
RANDOBJS      = $(RANDSRCS:.c=.o)
NBDSLIBS      =  ../../nbds.0.4.3/libnbds.a

OBJS = $(COMMONOBJS) $(OAHTTSLFOBJS) $(NBDSOBJS) $(RANDOBJS)

LIBS = ../../utils/libbpautils_thread.a

CFLAGS += $(INCDIRS)
CFLAGS += $(PTHREAD_CFLAGS)
LDFLAGS += $(PTHREAD_LDFLAGS)

all: parbpalign parbpalign_rand #  parbpalign_nbds


parbpalign: $(COMMONOBJS) $(OAHTTSLFOBJS) ../../utils/ht.o ../../utils/oahttslf.o
	$(LD) -o $@ $^  $(LIBS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)

parbpalign_nbds:  $(COMMONOBJS) $(NBDSOBJS) ../../utils/ht.o
	$(LD) -o $@ $^ $(LIBS) $(NBDSLIBS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)

parbpalign_rand: $(COMMONOBJS)  $(RANDOBJS) ../../utils/ht.o ../../utils/oahttslf.o
	$(LD) -o $@ $^  $(LIBS) $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)

depend: $(SRCS) $(HDRS)
	$(MAKEDEPEND) $(INCDIRS) $(NBDSINC) $(SRCS) $(HDRS) > $(DEPENDFILE)


# special rule for users of nbds - need to compile with GNU C99 
C99FLAGS = -std=gnu99 -Wall -pthread $(INCDIRS) $(NBDSINC)
ifeq ($(MODE),DEBUG)
    C99FLAGS     += $(CDEBUG)
else
    C99FLAGS     += $(COPTIMIZE)
endif
bpadynprog_nbds.o: bpadynprog_nbds.c
	$(CC) $(C99FLAGS) -c -o $@ $<


# make the TAGS file for EMACS
tags:
	$(MAKETAGS) $(SRCS) $(HDRS)


clean:
	$(RM) $(OBJS) httest.o
	$(RM) parbpalign parbpalign_nbds parbpalign_rand

include $(DEPENDFILE)


