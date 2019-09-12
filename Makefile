F77     =      ifort
#F77     =      gfortran
FFLAGS =  -ipo -O3 -no-prec-div -fp-model fast=2 -xHost # = fast w/o -static
#FFLAGS  =      -O2 -g -traceback -i8 -mcmodel=medium
#FFLAGS  =       -g -O2 -Wall -ffpe-summary='none' -fcray-pointer -mcmodel=medium -finteger-4-integer-8 #-fdefault-integer-8 -ffpe-trap=denormal 
LIBS    = -L${MKLROOT} -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lm -ldl
#FFLAGS  =       -C -Wall
#LIBS    =       -L/opt/f77/lib:/home/hmh/lib -leispack -lsu2
#LIBS    =       -L/home/hmh/lib -leispack -lsu2 -llapack -lblas
#LIBS    =       -L/home/hmh/lib -L/opt/intel_mkl_61/lib 
#LIBS    =       -L/home/hmh/lib -leispack -lsu2  -lblas

.f.o:
	$(F77) $(FFLAGS) -c -o $*.o $<

clean:
	rm -f *.o core

qual: qual.o
	$(F77) $(FFLAGS) -o qual.exe qual.o $(LIBS)

luise: luise.o
	$(F77) $(FFLAGS) -o luise.exe luise.o $(LIBS)

obem: obem.o
	$(F77) $(FFLAGS) -o obem.exe obem.o $(LIBS)

obelma: obelma.o
	$(F77) $(FFLAGS) -o obelma.exe obelma.o $(LIBS)

enemb: enemb.o
	$(F77) $(FFLAGS) -o enemb.exe enemb.o $(LIBS)
