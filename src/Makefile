CC = gcc
version = $(shell bash ../scripts/getVersion.sh)
date = $(shell bash ../scripts/getDate.sh)
sys = $(shell uname -s)

CFLAGS = -O3 -std=c99 -Wall -Wshadow -pedantic -DVERSION="\"$(version)\"" -DDATE="\"$(date)\"" -I/opt/homebrew/include
src = analysis.c annotation.c auxiliaryLcpTree.c eprintf.c expectedShulen.c fileListUnix.c interface.c \
      interval.c intervalStack.c lcpTree.c mainAlfy.c queryBTNode.c queryInterval.c ran.c sequenceData.c \
      sequenceUnion.c shulen.c stringUtil.c subjectNode.c
obj = analysis.o annotation.o auxiliaryLcpTree.o eprintf.o expectedShulen.o fileListUnix.o interface.o \
      interval.o intervalStack.o lcpTree.o mainAlfy.o queryBTNode.o queryInterval.o ran.o sequenceData.o \
      sequenceUnion.o shulen.o stringUtil.o subjectNode.o
lib = -lm -ldivsufsort -L/opt/homebrew/lib
exe = alfy

all: $(exe)

$(exe): $(obj)
	$(CC) $(CFLAGS) -o $(exe) $(obj) $(lib)

mac: MakefileMac
	rm -f alfy
	make -f MakefileMac
MakefileMac: Makefile configMac.awk
	awk -f configMac.awk Makefile > MakefileMac

clean:
	rm -f *.o alfy
test: $(exe)
	bash test.sh
