###############################################################################
#
# File      : Makefile
# Author    : Alex Stivala (astivala)
# Created   : June 2007
#
# $Id: Makefile 3184 2010-01-01 04:10:13Z alexs $
#
# Makefile for building the the C utils library containing lock-free hash
# table etc.
#
#
# The default target "all" will compile and build the library.
# The "test" target will build test harness executables.
# The "runtests" target will run the test harness executables
# The "clean" target will remove object and executable files and libraries.
# The "times" target will make table of times for tests
# The "graphs" target will make speedup graph from times, using R.
#
# set MODE=DEBUG to build with debugging and profiling on. Otherwise
# default is to build with optimizations on and no debug or profile.
#
# This Makefile uses implicit rules to compile .c files into .o files
# and some features specific to GNU Make.
#
#
# In order to build the test that use the Intel Threading Building Blocks
# for Open Source,
# these need to first be installed, download from
#   http://www.threadingbuildingblocks.org/
# then set env var TBB21_INSTALL_DIR to the install directory and
# source the tbbvars.sh (or tbbvars.csh) script in the appropriate 
# architecture/OS subdirctory to setup environment variables for building
# with the Intel Threading Building Blocks library.
# Developed with version 2.1 (20080605 Open Source Linux version).
# 
###############################################################################


# include the common.mk file from parent diretory, and override it with
# local.mk if it exists

include ../common.mk
-include ../local.mk

INCDIRS =  
LIB_THREAD_SRCS = bpautils.c httslf.c cellpool.c
LIB_NOTHREAD_SRCS = bpautils.c ht.c cellpool.c

TEST_SRCS =  httest.c httslftest.c oahttslftest.c
OTHER_SRCS = oahttslf.c
SRCS = $(LIB_THREAD_SRCS) $(LIB_NOTHREAD_SRCS) $(TEST_SRCS) $(OTHER_SRCS)

LIB_THREAD_OBJS  = $(LIB_THREAD_SRCS:.c=.o)
LIB_NOTHREAD_OBJS = $(LIB_NOTHREAD_SRCS:.c=.nothread.o)

TEST_OBJS = $(TEST_SRCS:.c=.o)
OTHER_OBJS = $(OTHER_SRCS:.c=.o)
OBJS = $(LIB_NOTHREAD_OBJS) $(LIB_THREAD_OBJS) $(TEST_OBJS) $(OTHER_OBJS)


CFLAGS += $(INCDIRS) 

TEST_EXES = httest httslftest oahttslftest
LIBS = libbpautils_thread.a libbpautils_nothread.a

HOSTNAME = ${shell hostname | cut -d. -f1}

GRAPHS = charikar.eps mundara.eps mungera.eps node00.eps \
         tbbhashmap.charikar.eps tbbhashmap.node00.eps \
		 charikar.testspeedupmulti.eps node00.testspeedupmulti.eps \
         mundara.testspeedupmulti2.eps \
         mungera.testspeedupmulti2.eps \
         tango.eps 

R       = R --vanilla --slave

all: libbpautils_thread.a libbpautils_nothread.a tests gprof-helper.so \
     numcores timeguard

tests: $(TEST_EXES)

runtests: tests
	for i in $(TEST_EXES) ; do  \
          $$i ; \
	  if [ $$? -ne 0  ]; then  \
	    echo $$i  FAILED ; \
          else  \
            echo $$i  PASSED ; \
          fi ; \
        done 

times: $(HOSTNAME).rtab $(HOSTNAME).oahttslftest.rtab
graphs: $(GRAPHS)

libbpautils_thread.a: $(LIB_THREAD_OBJS)
	$(AR) rv $@ $^

libbpautils_nothread.a: $(LIB_NOTHREAD_OBJS)
	$(AR) rv $@ $^


tbbhashmaptest: tbbhashmaptest.o tbbhashmap.o bpautils.o
	$(CXX) $(CXXFLAGS) $(PTHREAD_CFLAGS) -o $@ $^ -ltbbmalloc -ltbb 

%.testspeedupmulti.eps: plottestspeedupmulti.r
	$(R) -f plottestspeedupmulti.r

%.testspeedupmulti2.eps: plottestspeedupmulti2.r
	$(R) -f plottestspeedupmulti2.r

%.eps: %.rtab plottestspeedup.r 
	$(R) -f plottestspeedup.r --args $<

%.oahttslftest.rtab: #timetests.sh $(TEST_EXES) numcores
	numcores=`numcores` ; \
	timetests.sh httest oahttslftest $$numcores > $@

%.rtab: #timetests.sh $(TEST_EXES) numcores
	numcores=`numcores` ; \
	timetests.sh httest httslftest $$numcores > $@


%.nothread.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(INCS) -c -o $@ $<

oahttslf.o: oahttslf.c $(INLINE_ASM)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(INCS) $(PTHREAD_CFLAGS) $(INLINE_ASM) -c -o $@ $<


%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(PTHREAD_CFLAGS) -c -o $@ $<

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(INCS) $(PTHREAD_CFLAGS) -c -o $@ $<

httest: httest.o libbpautils_nothread.a
	$(LD) -o $@ $^ libbpautils_nothread.a $(LDFLAGS) $(LDLIBPATH) $(LDLIBS)

httslftest: httslftest.o libbpautils_thread.a
	$(LD) -o $@ $^  libbpautils_thread.a $(LDFLAGS) $(PTHREAD_LDFLAGS) $(LDLIBPATH) $(LDLIBS)

oahttslftest: oahttslftest.o oahttslf.o bpautils.o
	$(LD) -o $@ $^ $(LDFLAGS) $(PTHREAD_LDFLAGS) $(LDLIBPATH) $(LDLIBS)

nomemorytest: nomemorytest.o oahttslf.o bpautils.o
	$(LD) -o $@ $^ $(LDFLAGS) $(PTHREAD_LDFLAGS) $(LDLIBPATH) $(LDLIBS)

simpletest: simpletest.o bpautils.o
	$(LD) -o $@ $^ $(LDFLAGS) $(PTHREAD_LDFLAGS) $(LDLIBPATH) $(LDLIBS)

numcores: numcores.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $<

timeguard: timeguard.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $<

gprof-helper.so: gprof-helper.c
	gcc -shared -fPIC $(LDFLAGS) gprof-helper.c -o gprof-helper.so -lpthread -ldl

depend: $(SRCS) $(HDRS)
	$(MAKEDEPEND) $(INCDIRS) $(NBDSINC) $(SRCS) $(HDRS) > $(DEPENDFILE)


# make the TAGS file for EMACS
tags: TAGS
TAGS:
	$(MAKETAGS) $(SRCS) $(HDRS)


clean:
	$(RM) $(OBJS)
	$(RM) $(LIBS) $(TEST_EXES)
	$(RM) gprof-helper.so numcores timeguard
	$(RM) tbbhashmaptest.o tbbhashmap.o tbbhashmaptest

realclean: clean
	$(RM) $(HOSTNAME).rtab $(HOSTNAME).eps

include $(DEPENDFILE)


