FC := ifort

.SUFFIXES: .o .f90

MKLROOT=/opt/intel/mkl/

#/opt/intel/compilers_and_libraries_2018/mac/mkl/interfaces/include/intel64/lp64/ /opt/intel/mkl/lib/intel64_lin/

FC_OPTIONS =  -fpp -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include -g
LINKER_OPTIONS =  -g ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

%.o: %.f90
	$(FC) $(FC_OPTIONS) -c $<

all: gslib2 sgsim

OBJS= kinds_kit.o geometry.o variography.o sorting.o searching.o solvers.o simulation.o

gslib2: $(OBJS)
#	$(FC) -g -shared -fPIC -o gslib2.so kinds_kit.o geometry.o variography.o simulation.o

sgsim: gslib2 sgsim.o
	$(FC) $(OBJS) sgsim.o -o sgsim $(LINKER_OPTIONS)

sgsim.o: sgsim.for
	$(FC) -c sgsim.for

clean:
	rm -fr *.mod *.o
