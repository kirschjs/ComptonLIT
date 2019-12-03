F77     =       gfortran
FFLAGS  =       -Ofast -funroll-loops -ffpe-summary='none' -fcray-pointer -mcmodel=large 
#-ffpe-trap=denormal -pg -funroll-loops -g 

## orig settings
#F77     =        ifort
#FFLAGS  =       -O1 -g -traceback -i8 -mcmodel=medium
#FFLAGS  =       -g -O2 -mcmodel=medium #-ffpe-trap=denormal 
LIBS    = -L. -llapack -lblas
#LIBS    = -L${MKLROOT} -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5# -lpthread -lm -ldl
# ohne -g und -traceback erzeugt spoleinv runtime error in smatest,rrandb subroutinen wenn auf die variable nbess zugegriffen wird;
# mit den optionen tritt der fehler nicht auf;
#FFLAGS  =       -C -Wall
#LIBS    =       -L/opt/f77/lib:/home/hmh/lib -leispack -lsu2
#LIBS    =       -L/home/hmh/lib -leispack -lsu2 -llapack -lblas
#LIBS    =       -L/home/hmh/lib -L/opt/intel_mkl_61/lib -leispack -lsu2 -lmkl -lmkl_lapack -lmkl_lapack32 -lguide -lblas
#LIBS    =       -L/opt/intel/mkl/lib/ia32 
#LIBS    =       -L/home/hmh/lib -leispack -lsu2  -lblas
#LIBS    =       -L/home/jkirscher/zweistein/lib -lmkl_lapack -lguide -lpthread

.f.o:
	$(F77) $(FFLAGS) -c -o $*.o $<


clean:
	rm -f *.o core

Tmat: Tmat.o 
	$(F77) $(FFLAGS) -o Tmat.exe Tmat.o $(LIBS)

best: best.o 
	$(F77) $(FFLAGS) -o best.exe best.o $(LIBS)

QUAFL_N: QUAFL_N.o par/QUAF
	$(F77) $(FFLAGS) -o QUAFL_N.exe QUAFL_N.o 

LUDW_CN: LUDW_CN.o par/LUCN
	$(F77) $(FFLAGS) -o LUDW_CN.exe LUDW_CN.o $(LIBS)

DRQUA_AK_N: DRQUA_AK_N.o par/DRQUA_AK
	$(F77) $(FFLAGS) -o DRQUA_AK_N.exe DRQUA_AK_N.o $(LIBS)

DRQUA_I: DRQUA_I.o par/DRQUA_I
	$(F77) $(FFLAGS) -o DRQUA_I.exe DRQUA_I.o $(LIBS)

DRQUA_I_2: DRQUA_I_2.o par/DRQUA_I_2
	$(F77) $(FFLAGS) -o DRQUA_I_2.exe DRQUA_I_2.o $(LIBS)

TDR2ENDM_AK_N: TDR2ENDM_AK_N.o par/DR2END_AK
	$(F77) $(FFLAGS) -o TDR2ENDM_AK_N.exe TDR2ENDM_AK_N.o $(LIBS)

TDR2END_AK: TDR2END_AK.o par/DR2END_AK
	$(F77) $(FFLAGS) -o TDR2END_AK.exe TDR2END_AK.o $(LIBS)

TDR2END_I_2: TDR2END_I_2.o par/DR2END_AK_I_2
	$(F77) $(FFLAGS) -o TDR2END_I_2.exe TDR2END_I_2.o $(LIBS)

DR2END_AK_exp: DR2END_AK_exp.o par/DR2END_AK
	$(F77) $(FFLAGS) -o DR2END_AK_exp.exe DR2END_AK_exp.o $(LIBS)

DR2END_AK: DR2END_AK.o par/DR2END_AK
	$(F77) $(FFLAGS) -o DR2END_AK.exe DR2END_AK.o $(LIBS)

DR2END_AK_I: DR2END_AK_I.o par/DR2END_AK_I
	$(F77) $(FFLAGS) -o DR2END_AK_I.exe DR2END_AK_I.o $(LIBS)

DR2END_AK_I_2: DR2END_AK_I_2.o par/DR2END_AK_I_2
	$(F77) $(FFLAGS) -o DR2END_AK_I_2.exe DR2END_AK_I_2.o $(LIBS)

TDR2ENDM: TDR2ENDM.o par/DR2END
	$(F77) $(FFLAGS) -o TDR2ENDM.exe TDR2ENDM.o $(LIBS)

TDRQUA: TDRQUA.o par/DRQUA
	$(F77) $(FFLAGS) -o TDRQUA.exe TDRQUA.o $(LIBS)

DR2ENDM: DR2ENDM.o par/DR2END
	$(F77) $(FFLAGS) -o DR2ENDM.exe DR2ENDM.o $(LIBS)

NOHAWRI: NOHAWRI.o par/BENDM
	$(F77) $(FFLAGS) -o NOHAWRI.exe NOHAWRI.o $(LIBS)

BENDM: BENDM.o par/BENDM
	$(F77) $(FFLAGS) -o BENDM.exe BENDM.o $(LIBS)

DRQUA: DRQUA.o par/DRQUA
	$(F77) $(FFLAGS) -o DRQUA.exe DRQUA.o $(LIBS)

COBER: COBER.o par/COBER
	$(F77) $(FFLAGS) -o COBER.exe COBER.o $(LIBS)

DRLUD: DRLUD.o par/DRLUD
	$(F77) $(FFLAGS) -o DRLUD.exe DRLUD.o $(LIBS)

DRLUD_I: DRLUD_I.o par/DRLUD_I
	$(F77) $(FFLAGS) -o DRLUD_I.exe DRLUD_I.o $(LIBS)

DRLUD_I_2: DRLUD_I_2.o par/DRLUD_I_2
	$(F77) $(FFLAGS) -o DRLUD_I_2.exe DRLUD_I_2.o $(LIBS)

DROBER: DROBER.o par/DROBER
	$(F77) $(FFLAGS) -o DROBER.exe DROBER.o $(LIBS)

DROBER_I: DROBER_I.o par/DROBER_I
	$(F77) $(FFLAGS) -o DROBER_I.exe DROBER_I.o $(LIBS)

DROBER_I_2: DROBER_I_2.o par/DROBER_I_2
	$(F77) $(FFLAGS) -o DROBER_I_2.exe DROBER_I_2.o $(LIBS)

KOBER: KOBER.o par/KOBER
	$(F77) $(FFLAGS) -o KOBER.exe KOBER.o $(LIBS)

umnum: umnum.o 
	$(F77) $(FFLAGS) -o umnum.exe umnum.o $(LIBS)

samel: samel.o 
	$(F77) $(FFLAGS) -o samel.exe samel.o $(LIBS)

nana_exp: nana_exp.o 
	$(F77) $(FFLAGS) -o nana_exp.exe nana_exp.o $(LIBS)

nana_c: nana_c.o 
	$(F77) $(FFLAGS) -o nana_c.exe nana_c.o $(LIBS)

nana_R: nana_R.o 
	$(F77) $(FFLAGS) -o nana_R.exe nana_R.o $(LIBS)

streimat: streimat.o par.endm
	$(F77) $(FFLAGS) -o streimat.exe streimat.o $(LIBS)

S-POLE_RC: S-POLE_RC.o par/verma 
	$(F77) $(FFLAGS) -o S-POLE_RC.exe S-POLE_RC.o $(LIBS)

S-POLE_INV: S-POLE_INV.o par/verma coul.o
	$(F77) $(FFLAGS) -o S-POLE_INV.exe S-POLE_INV.o coul.o $(LIBS)

S-POLE_INV-old: S-POLE_INV-old.o par/verma coul.o
	$(F77) $(FFLAGS) -o S-POLE_INV-old.exe S-POLE_INV-old.o coul.o $(LIBS)

S-POLE_Rmat: S-POLE_Rmat.o par/verma coul.o
	$(F77) $(FFLAGS) -o S-POLE_Rmat.exe S-POLE_Rmat.o coul.o $(LIBS)

S-POLE_cmplx: S-POLE_cmplx.o par/verma coul.o
	$(F77) $(FFLAGS) -o S-POLE_cmplx.exe S-POLE_cmplx.o coul.o $(LIBS)

S-POLE_kohn: S-POLE_kohn.o par/verma coul.o
	$(F77) $(FFLAGS) -o S-POLE_kohn.exe S-POLE_kohn.o coul.o $(LIBS)

S-POLE_zget: S-POLE_zget.o par/verma coul.o
	$(F77) $(FFLAGS) -o S-POLE_zget.exe S-POLE_zget.o coul.o $(LIBS)

S-POLE_zget_PWBA: S-POLE_zget_PWBA.o par/verma coul.o
	$(F77) $(FFLAGS) -o S-POLE_zget_PWBA.exe S-POLE_zget_PWBA.o coul.o $(LIBS)

