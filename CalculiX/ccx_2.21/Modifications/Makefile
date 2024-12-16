DIR=/cygdrive/c/Users/gmadireddy/CalculiX/Windows/ccx_prool/SPOOLES.2.2

CFLAGS = -Wall -O2  -I $(DIR) -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DNETWORKOUT 
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



LIBS = \
       $(DIR)/MT/src/spoolesMT.a \
	   $(DIR)/spooles.a \
	/cygdrive/c/Users/gmadireddy/CalculiX/Windows/ccx_prool/ARPACK/libarpack_INTEL.a \
       -lpthread -lm -lc

ccx_2.21.exe: $(OCCXMAIN) ccx_2.21.a  $(LIBS)
	./date.pl; $(CC) $(CFLAGS) -c ccx_2.21.c; $(FC)  -Wall -O2 -o $@ $(OCCXMAIN) ccx_2.21.a $(LIBS) -fopenmp


ccx_2.21.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?
