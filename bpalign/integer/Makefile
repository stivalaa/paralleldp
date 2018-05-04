###############################################################################
#
# File      : Makefile
# Author    : Alex Stivala (astivala)
# Created   : June 2007
#
# $Id: Makefile 3146 2009-12-27 01:49:26Z alexs $
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

R       = R --vanilla --slave

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

GRAPHS = mundara.eps \
          mundara-hashcountincrease.eps mundara-hascount.eps \
          mundara-hascountoverreuse.eps \
		  mundara_array.eps \
          mundara_array-hashcountincrease.eps mundara_array-hascount.eps \
          mundara_array-hascountoverreuse.eps \
		  mungera.eps \
		  mungera_array.eps \
		  tango.eps tango_array.eps

RTABS = mundara.rtab mundara_array.rtab \
        mungera.rtab mungera_array.rtab \
		mungera_norandomization_instrument.rtab \
		mundara_norandomization_instrument.rtab \
		mundara_norandomization.rtab \
		mungera_bottomup.rtab mundara_bottomup.rtab \
		mungera_instrument.rtab \
		tango.rtab tango_norandomization.rtab \
		tango_bottomup.rtab tango_bottomup.rtab \

BASETIMES = mundara.basetime mungera.basetime tango.basetime

all: parbpalign parbpalign_rand #  parbpalign_nbds

times : $(RTABS)
graphs: mundara.eps mundara_array.eps mungera.eps mungera_array.eps
	$(R) -f plotalgoeffic.r

basetimes: $(BASETIMES)


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

mundara_norandomization_instrument.rtab: mundara_norandomization.rtab
	cp $< $@

%_array.eps: %_array.rtab %_bottomup.rtab plotarrayspeedup.r
	$(R) -f plotarrayspeedup.r --args $*_array.rtab $*_bottomup.rtab

%-times-tabsize2.eps: %_tablesize26.rtab %_tablesize24.rtab ../../mosp/plottimestabsize2.r
	$(R) -f ../../mosp/plottimestabsize2.r --args $< 

%-tabsize2.eps: %_tablesize26.rtab %_tablesize24.rtab ../../mosp/plottimestabsize2.r
	$(R) -f ../../mosp/plotspeeduptabsize2.r --args $< 


%.eps: %.rtab ../../mosp/plotspeedup.r plotalgoeffic.r
	$(R) -f ../../mosp/plotspeedup.r --args $*.rtab

%.rtab:
	results2rtab.sh $*_output > $@


%.basetime: %_output
	grep -v '^#' $*_output/bpalign_0.out | awk 'BEGIN {min=9999;max=0;} {sum += $$5 ; if ($$5 < min) min = $$5 ; if ($$5 > max) max = $$5;} END {printf("total %.0f\nmin  %.0f\nmax  %.0f\n", sum/1000,min/1000,max/1000)}' > $@

# make the TAGS file for EMACS
tags:
	$(MAKETAGS) $(SRCS) $(HDRS)


clean:
	$(RM) $(OBJS) httest.o
	$(RM) parbpalign parbpalign_nbds parbpalign_rand
	$(RM) $(GRAPHS)

realclean: clean
	$(RM) $(RTABS)

include $(DEPENDFILE)


