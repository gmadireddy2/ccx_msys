DIR=../../../
LIBDIR=../../../../../../../../msys64/mingw64/lib

CFLAGS = -Wall -O3 -DUSE_MT -pthread -I $(DIR)/ARPACK -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DNETWORKOUT
FFLAGS = -Wall -O2

CC=cc
FC=gfortran

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

include Makefile.inc

SCCXMAIN = ccx_2.21.c

OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)

DIR2=../../../../../

LIBS = \
	   $(DIR)/SPOOLES.2.2/MT/src/spoolesMT.a \
	   $(DIR)/SPOOLES.2.2/spooles.a \
	   $(DIR)/ARPACK/libarpack_INTEL.a \
	   $(LIBDIR)/libpthread.a \
	   $(LIBDIR)/libm.a \
	   $(DIR2)/Charlie/abqsub.dll

ccx_2.21.exe: $(OCCXMAIN) ccx_2.21.a $(LIBS)
	./date.pl; $(CC) $(CFLAGS) -c ccx_2.21.c; $(FC) -Wall -O2 -o $@ $(OCCXMAIN) ccx_2.21.a $(LIBS) -fopenmp

ccx_2.21.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?
