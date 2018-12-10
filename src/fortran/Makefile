FC := ifort

.SUFFIXES: .o .f90

MKLROOT=/opt/intel/mkl/

MKL_LIBS=${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
MKL_INCS=-I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include

#/opt/intel/compilers_and_libraries_2018/mac/mkl/interfaces/include/intel64/lp64/ /opt/intel/mkl/lib/intel64_lin/
TRACE_OPTIONS = -DNO_TRACE 
OPTIMIZATION_FLAGS = -O4

FC_OPTIONS =  -fPIC $(TRACE_OPTIONS) -fpp  -g $(OPTIMIZATION_FLAGS) -qopenmp # -init=snan,arrays -traceback
LINKER_OPTIONS =  -g -qopenmp

%.o: %.f90
	$(FC) $(FC_OPTIONS) -c $<

%.o: %.for
	$(FC) $(FC_OPTIONS) -c $<

all: gslib2 sgsim

OBJS= kinds_kit.o geometry.o interpolation.o sorting.o searching.o gaussian_utils.o variography.o solvers.o kriging.o sgsimulation.o 

gslib2: $(OBJS)
	$(FC) -g -shared -fPIC -o gslib2.so $(OBJS)

sgsim: gslib2 sgsim.o
	$(FC) $(OBJS) sgsim.o -o sgsim $(LINKER_OPTIONS)

kt3d: gslib2 kt3d.o
	$(FC) $(OBJS) kt3d.o -o kt3d $(LINKER_OPTIONS) 

clean:
	rm -fr *.mod *.o
