FC := ifort

.SUFFIXES: .o .f90

MKLROOT=/opt/intel/compilers_and_libraries_2018/mac/mkl

#/opt/intel/compilers_and_libraries_2018/mac/mkl/interfaces/include/intel64/lp64/

FC_OPTIONS = -I${MKLROOT}/interfaces/include/intel64/lp64 -I${MKLROOT}/interfaces/include -g
LINKER_OPTIONS =  -g ${MKLROOT}/lib/libmkl_lapack95_lp64.a ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -lpthread -lm -ldl

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
