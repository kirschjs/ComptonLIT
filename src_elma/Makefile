F77     =      gfortran
FFLAGS  =       -Ofast -finteger-4-integer-8 -funroll-loops -funroll-all-loops  
# -g -pg  -ffpe-summary='none' -fcray-pointer  -mcmodel=medium
#-fno-align-commons#-fdefault-integer-8 -ffpe-trap=denormal  -fno-align-commons 

LIBS    = -L. -llapack -lblas

#FFLAGS  =       -Ofast -g -traceback -i8 -mcmodel=medium
#FFLAGS =  -ipo -O3 -no-prec-div -fp-model fast=2 -xHost # = fast w/o -static
#FFLAGS  =      -O2 -g -traceback -i8 
#LIBS    = -L${MKLROOT} -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lm -ldl

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

enem_joki: enem_joki.o
	$(F77) $(FFLAGS) -o enem_joki.exe enem_joki.o $(LIBS)
