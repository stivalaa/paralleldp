###############################################################################
#
# File      : common.mk
# Author    : Alex Stivala (astivala)
# Created   : July 2008
#
# $Id: Makefile 1690 2008-07-16 04:35:47Z astivala $
#
# Definitions of compilers and compile options for all Makefiles.
# Use GNU make.
#
# The _XOPEN_SOURCE and _XOPEN_SOURCE_EXTENDED defines are set to 1 to
# to allow certain prototypes to be included properly and indicate
# that the program is to work with the X/Open XPG4v2 (SUS)
# standard. Specifically, this wasy getopt() will compile and work
# nicely with no warnings on Solaris, and also allow the use of
# strcasecmp() and strdup() with no warnings. 
#
# set MODE=DEBUG to build with debugging and profiling on. Otherwise
# default is to build with optimizations on and no debug or profile.
#
#
###############################################################################


CC         = gcc
CXX        = g++
CPPFLAGS   = -D_XOPEN_SOURCE=1 -D_XOPEN_SOURCE_EXTENDED=1
CPPFLAGS  += -DUSE_INSTRUMENT
CDEBUG     = -g  -O0 -pg -DDEBUG
COPTIMIZE  = -O3 
CFLAGS     =  -Wall
CXXFLAGS   = 
ifeq ($(MODE),DEBUG)
    CFLAGS     += $(CDEBUG)
    CXXFLAGS   += $(CDEBUG)
else
    CFLAGS     += $(COPTIMIZE)
    CXXFLAGS   += $(COPTIMIZE)
endif
#              the following warnings are not implied by -Wall
CFLAGS     += -Wextra -Wfloat-equal  \
              -Wdeclaration-after-statement -Wundef -Wshadow \
              -Wpointer-arith -Wbad-function-cast -Wcast-qual -Wcast-align\
              -Wwrite-strings -Wmissing-prototypes \
              -Wmissing-declarations -Wunreachable-code

PTHREAD_CLFAGS = -pthread -DUSE_THREADING
LD         = gcc
LDFLAGS    = 
ifeq ($(MODE),DEBUG)
    LDFLAGS    += -g 
    LDFLAGS    += -pg  # for profiler gprof
endif
LDLIBPATH  = 
LDLIBS     = -lm 
STREAMFLOW_DIR=$(HOME)/phd/paralleldp/streamflow
PTHREAD_LDFLAGS = -pthread -L$(STREAMFLOW_DIR) -lstreamflow

MAKEDEPEND = gcc -MM $(CPPFLAGS)
DEPENDFILE = .depend

# Program to build TAGS file for EMACS
MAKETAGS   = etags

