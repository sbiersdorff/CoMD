#All the paths need to be adjusted to your system. This file provides the 
#base SRCS that you need.

STREAMSDK =  /home/cuda/sdk4.1.28/OpenCL/common/
INCS = -I$(STREAMSDK)/inc 
LIBS = -L/usr/lib64/ -lOpenCL -lm


HEAD = -I. -I../src-flat/include
SRCS = pmdOCL.c helpers.c cl_utils.c docompute.c\
	../src-flat/cmdLineParser.c \
           ../src-flat/domains.c \
           ../src-flat/eam.c \
           ../src-flat/ljForce.c \
           ../src-flat/ic_fcc.c \
           ../src-flat/mycommand.c \
           ../src-flat/read.c \
           ../src-flat/write.c \
           ../src-flat/utility.c \
           ../src-flat/timestep.c

OBJS = $(SRCS:.c=.o)

DEPS = *.h
CC = gcc

all: CoMD

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(INCS) $(HEAD)

cruft: $(OBJS)
	$(CC) -o ../CoMDOCL -O3 $(OBJS) $(HEAD)  $(INCS) $(LIBS)

clean:
	rm *.o ../src-flat/*.o

