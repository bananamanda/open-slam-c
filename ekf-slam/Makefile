DEBUG=0
OMP=-fopenmp -DOMP
CC 	:= gcc
CCFLAGS	:= -g -O3 -Wall -DDEBUG=$(DEBUG)
LDFLAGS	:= -lm
DDIR = ./test
UPATH = ../util/

STYLE	:= llvm

CFILES = simulator.c $(UPATH)util.c $(UPATH)cycletimer.c $(UPATH)instrument.c $(UPATH)particles.c
HFILES = $(UPATH)util.h $(UPATH)cycletimer.h $(UPATH)instrument.h $(UPATH)particles.h

.PHONY: all clean

all: slam-seq slam-omp

slam-seq: $(CFILES) $(HFILES) $(UPATH)fake_omp.h $(UPATH)fake_omp.c
	clang-format -style=$(STYLE) -i $(CFILES)
	$(CC) $(CCFLAGS) -o slam-seq $(CFILES) $(UPATH)fake_omp.c $(LDFLAGS)

slam-omp: $(CFILES) $(HFILES)
	clang-format -style=$(STYLE) -i $(CFILES)
	$(CC) $(CCFLAGS) $(OMP) -o slam-omp $(CFILES) $(LDFLAGS)

clean:
	rm -f $(TARGETS) $(OBJ)
