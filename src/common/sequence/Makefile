CC=gcc
# with profiling & debugging:
# CFLAGS= -m64 -O0 -g -Wall -Wshadow -pedantic   -pg -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 #-mcmodel=medium
# without profiling & debugging:
#CFLAGS= -g -Wall -Wshadow -pedantic -DSTRMAT -DSTATS
#CFLAGS= -O3 -DSTRMAT -DSTATS

#32 bit version
#when the option VER32 is set, the source will be compiled as the 32 bit version          
# preprocessor options: VER32 -> 32-bit version, DEBUG -> debug version, UNIX -> unix version (as opposed to windows version); m64 for unix/mac, exclude for linux!
CFLAGS= -O3 -Wall -Wshadow -pedantic -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -std=c99 -DVER32 -DUNIX #-m64  

#64 bit version
# option VER32 is not set
CFLAGS64= -O3 -Wall -Wshadow -pedantic -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -std=c99 -DUNIX #-m64  

################################################
#include directories
INCLUDE_COMMON_DIR = ..
INCLUDE_UTIL_DIR = ../util
# INCLUDE_DIR32 = ../../../Include/DeepShallow
# INCLUDE_DIR64 = ../../../Include/DeepShallow64
INCLUDE_DIR32 = ../../../externSrc/deepShallow/src
INCLUDE_DIR64 = ../../../externSrc/deepShallow64/src

INCLUDES = -I$(INCLUDE_COMMON_DIR) -I$(INCLUDE_UTIL_DIR) 
INCLUDES32 = $(INCLUDES) -I$(INCLUDE_DIR32)
INCLUDES64 = $(INCLUDES) -I$(INCLUDE_DIR64)

#library directories
LIB_DIR = ../../../lib
LIB_DIR32 = ../../../lib/deepShallow/
LIB_DIR64 = ../../../lib/deepShallow64/

LIBS32 = -lm -L$(LIB_DIR) -lutil
LIBS64 = -lm -L$(LIB_DIR) -lutil64

################################################
#
# The source, object and library files
#
SRCFILES = expectedShulen.c sequenceData.c sequenceUnion.c shulen.c sequenceAnalysis.c

OBJFILES = $(addsuffix 32.o, $(basename $(SRCFILES)))
OBJFILES64 = $(addsuffix 64.o, $(basename $(SRCFILES)))

SEQ_LIB = libsequence.a 
SEQ64_LIB = libsequence64.a
################################################

VERSION = 2.0.0
#
# The dependencies for each of the *.o and *64.o files.
#
.PHONY : all
all: $(SEQ_LIB) $(SEQ64_LIB)

$(SEQ_LIB) : $(OBJFILES)
	ar rcs $(SEQ_LIB) $(OBJFILES)

$(SEQ64_LIB) : $(OBJFILES64)
	ar rcs $(SEQ64_LIB) $(OBJFILES64)

#$(OBJFILES): $(SRCFILES) # this doesn't work!
%32.o : %.c
	$(CC) $(CFLAGS) $(INCLUDES32) -c -o $*32.o $*.c 
#$(LIBS32)

%64.o : %.c
	$(CC) $(CFLAGS64) $(INCLUDES64) -c -o $*64.o $*.c 
#$(LIBS64)

###############################################

#
# Other Standard make rules
#
lint : 
	lint -x $(SRCFILES) | more

clean:
	rm -f *.o *.a

remove: 
	make clean


