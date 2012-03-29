# Emacs options: -*- makefile -*-
#
# Makefile pruned from cellMD
# Written by Sriram Swaminarayan
# 
##############################################################
# Description of magic variables:
#
# DO_SINGLE=yes|no  controls whether single precision or double precision
#
# DO_GZIP=yes|no    controls whether gzip support is included
#
# TARGET:  the name of the executable
#
# OPT:            Optimizatin flag, default -O5, set to -g to debug
#
# USER_FLAGS:     The compilation flags
#
# USER_LIBS:      This is additional libraries needed for compiling.
#                 THis needs more work
#
# DIR_SRC:        The source directory
#
##################################################################

TARGET = cruft

#
# flags that control behavior
#
# default optimization
OPT ?= -O5

# GZIP support for input
DO_GZIP ?= yes

# enable single precision
DO_SINGLE ?= no

DIR_SRC = src
SRCS = pmd.c \
       cmdLineParser.c \
       domains.c \
       eam.c \
       ljForce.c \
       mycommand.c \
       read.c \
       write.c \
       utility.c \
       timestep.c

##########################
# Set user defined C FLAGS
##########################
USER_CFLAGS = $(OPT)
ifeq ($(DO_SINGLE),yes)
  USER_CFLAGS += -DSINGLE
else
  USER_CFLAGS += -DDOUBLE
endif

USER_CFLAGS += -I$(DIR_TOP)/include
USER_CFLAGS += -I$(DIR_TOP)/$(DIR_SRC)/include
USER_CFLAGS += -I/usr/include

USER_LIBS = -lm
ifeq ($(DO_GZIP),yes)
  USER_CFLAGS += -DGZIPSUPPORT
  USER_LIBS += -lz
endif

ifeq ($(DO_NOSELF),yes)
  USER_CFLAGS += -DNOSELF=1
endif

ifeq ($(DO_NOLIST),yes)
  USER_CFLAGS += -DNOLIST=1
endif

#OBJS = $(patsubst %.c,$(DIR_OBJ)/%.o,$(SRCS))
OBJS = $(patsubst %.c,%.o,$(SRCS))

# set up paths

DIR_TOP ?= $(shell pwd)
DIR_OBJ = $(DIR_TOP)/obj

VPATH = $(DIR_TOP)
VPATH += :$(DIR_TOP)/$(DIR_SRC)
DIR_MYOBJS = $(DIR_OBJ)
CC ?= gcc

CFLAGS = $(USER_CFLAGS)


##########
# Targets
##########
%.o : %.c
	cd $(DIR_OBJ); $(CC) $(CFLAGS) -c $< -o $@

flat: dirs
	$(MAKE) cpu DIR_SRC=src-flat

cpu: dirs
	cd $(DIR_OBJ); $(MAKE) -f $(DIR_TOP)/Makefile cpu0


cpu0: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(DIR_TOP)/$(TARGET) $(USER_LIBS)

test:
	@echo ''
	@echo it lives in $(DIR_TOP)' !'
	@echo $(OBJS)

clean: rmdirs
	touch DUMMY~; rm -f *~ */*.o */*.d */*~ */*/*~ lib*.a $(TARGET)


dirs:
	@if [ ! -e $(DIR_OBJ) ] ; then mkdir $(DIR_OBJ) ; fi

rmdirs:
	@if [ -e $(DIR_OBJ) ] ; then /bin/rm -rf $(DIR_OBJ) ; fi

#########--end--#########
export 
