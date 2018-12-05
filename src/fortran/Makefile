FC := gfortran

all: gslib2

gslib2: kinds_kit.o geometry.o variography.o searching.o sgsim.o sorting.o
	$(FC) -g -shared -fPIC -o gslib2.so kinds_kit.o geometry.o variography.o sgsim.o

geometry.o: geometry.f90
	$(FC) -c geometry.f90

searching.o: searching.f90
	$(FC) -c searching.f90

sorting.o: sorting.f90
	$(FC) -c sorting.f90

variography.o: variography.f90
	$(FC) -c variography.f90

kinds_kit.o: kinds_kit.f90
	$(FC) -c kinds_kit.f90

sgsim.o: sgsim.f90
		$(FC) -c sgsim.f90

clean:
	rm -fr *.mod *.o
