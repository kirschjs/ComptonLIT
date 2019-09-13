      REAL FUNCTION CLG(J1,J2,J3,M1,M2)
C
C     CLG BERECHNET DIE CLEBSCH-GORDAN-KOEFFIZIENTEN
C     (J1/2,M1/2;J2/2,M2/2|J3/2,(M1+M2)/2) NACH
C     EDMONDS 'ANGULAR MOMENTUM IN QUANTUM MECHANICS',
C     PRINCETON, 1960 GLEICHUNGEN (3.10.60), (3.7.3)
C     UND TABELLE 2 (1. GLEICHUNG)
C
C     BENUTZT COMMON /COMY/ MIT DEN LOGRITHMEN DER
C     FAKULTAETEN
C
C     M. UNKELBACH 1989
C     LETZTE AENDERUNG: 06.02.89
C
C
      INTEGER JW1, JW2, JW3, MW1, MW2, MW3, JSUM, JSUM1,
     *        JDIF1, JDIF2, JDIF3, JMSUM1, JMSUM2, JMSUM3,
     *        JMDIF1, JMDIF2, JMDIF3, JJM1, JJM2, IMAX, IMIN,
     *        I, J1, J2, J3, M1, M2
C
      REAL FAKLN, CLGH
C
      COMMON /COMY/ FAKLN(0:99)
C     FAKLN(I) = LOG(I!)
C
C
C
C
      JW1=J1
      JW2=J2
      JW3=J3
      MW1=M1
      MW2=M2
C
C     CHECK, OB CLG = 0
      CLG=0.
      IF (JW1.LT.IABS(MW1)) RETURN
      IF (JW2.LT.IABS(MW2)) RETURN
      IF (JW3.GT.JW1+JW2.OR.JW3.LT.IABS(JW1-JW2)) RETURN
      MW3=MW1+MW2
      IF (JW3.LT.IABS(MW3)) RETURN
      JMSUM1=JW1+MW1
      JMSUM2=JW2+MW2
      JMSUM3=JW3+MW3
      IF (MOD(JMSUM1,2).EQ.1) RETURN
      IF (MOD(JMSUM2,2).EQ.1) RETURN
      IF (MOD(JMSUM3,2).EQ.1) RETURN
C
C
      JSUM=(JW1+JW2+JW3)/2
      JSUM1=JSUM+1
      JDIF1=JSUM-JW1
      JDIF2=JSUM-JW2
      JDIF3=JSUM-JW3
C
      IF (IABS(MW1)+IABS(MW2).EQ.0) GOTO 100
C
C     NORMALE CLEBSCH-GORDAN-KOEFFIZIENTEN
      JMSUM1=JMSUM1/2
      JMDIF1=JMSUM1-MW1
      JMSUM2=JMSUM2/2
      JMDIF2=JMSUM2-MW2
      JMSUM3=JMSUM3/2
      JMDIF3=JMSUM3-MW3
      JJM1=JDIF1+JMDIF1
      JJM2=JDIF3-JMDIF1
      IMIN=MAX0(0,-JJM2)
      IMAX=MIN0(JMDIF1,JMDIF3)
C
      CLGH=0.
      DO 50, I=IMIN, IMAX
       CLGH=CLGH+DBLE(1-2*MOD(I,2))*
     *     EXP(FAKLN(JMSUM1+I)+FAKLN(JJM1-I)-FAKLN(I)-FAKLN(JMDIF1-I)-
     *         FAKLN(JMDIF3-I)-FAKLN(JJM2+I))
50    CONTINUE
C
      IF (IMIN.GT.IMAX) CLGH=1.
      CLGH=CLGH*EXP((FAKLN(JDIF3)+FAKLN(JMDIF1)+FAKLN(JMDIF2)+
     *             FAKLN(JMDIF3)+FAKLN(JMSUM3)-FAKLN(JSUM1)-
     *             FAKLN(JDIF1)-FAKLN(JDIF2)-FAKLN(JMSUM1)-
     *             FAKLN(JMSUM2)+FAKLN(JW3+1)-FAKLN(JW3))*.5)
      CLG=CLGH*DBLE(1-2*MOD(JMDIF1,2))
C
C     ENDE DER BERECHNUNG FUER NORMALE CLEBSCH-GORDAN-KOEFFIZIENTEN
      RETURN
C
C
C
100   CONTINUE
C     PARITAETSCLEBSCH
C
      IF (MOD(JSUM,2).EQ.1) RETURN
C
      CLGH=EXP((FAKLN(JDIF1)+FAKLN(JDIF2)+FAKLN(JDIF3)-FAKLN(JSUM1)+
     *         FAKLN(JW3+1)-FAKLN(JW3))*.5+
     *        FAKLN(JSUM/2)-FAKLN(JDIF1/2)-FAKLN(JDIF2/2)-
     *        FAKLN(JDIF3/2))
      CLG=CLGH*DBLE(1-2*MOD((JSUM+JW1-JW2)/2,2))
C
C
C     ENDE DER RECHNUNG FUER PARITAETSCLEBSCH
      RETURN
      END
C
      REAL FUNCTION P2L(A,N) 
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C    P2L=PLM WITH M=2 (EDMONDS) 
      COMMON/PLM/ XLM,Y(20) 
       IF(N.GE.2)  GOTO 2 
      P2L=0.
      RETURN
2      IF(ABS(A-XLM).GE.1.E-10) F=PLEGEN(A,N) 
       F=N
       P2L=-(A*(F-1.)*P1L(A,N)-(F+1.)*P1L(A,N-1))/SQRT(1.-A*A)
      RETURN                                                                 458
      END                                                                    459
      REAL FUNCTION P1L(A,N) 
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     P1L=PLM WITH M=1 (EDMONDS)
      COMMON/PLM/XLM,Y(20)
      IF(N.GE.1) GOTO 2 
      P1L=0.
      RETURN
2     IF(ABS(A-XLM).GE.1.E-10)  F=PLEGEN(A,N) 
      F=N 
      P1L=-F*(A*Y(N+1)-Y(N))/SQRT(1.-A*A) 
      RETURN                                                                 686
      END                                                                    687
      REAL FUNCTION Z1(A,B,C,D,E,F)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     Z FAKTOR AUS WELTON 
      RNULL=0.
      Z1=(2.*A+1.)*(2.*B+1.)*(2.*C+1.)*(2.*D+1.) 
      Z1=SQRT(Z1)*YG(A,C,F,RNULL,RNULL)*Y(A,B,C,D,E,F)
      RETURN                    
      END                                                                    692
      FUNCTION IAD(J,I,L,N) 
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C    J=2JD, L=LD , I=NUMMER KANALSPIN, N=1 EINGANG,N=2 AUSGANG
C    IAD=SUM(IZ=1 TO I) MIN0(2JD,2S(IZ)) +LD -(JD+S(I)) 
      COMMON SP(2,3)
      IAD=0                                                                  713
      IZ=1                                                                   714
1     IP=2*SP(N,IZ) 
      IAD=IAD+MIN0(J,IP)                                                     716
      IZ=IZ+1                                                                717
      IF(IZ-I)1,1,2                                                          718
2     IAD=IAD+L-(J+IP)/2 +IZ-1
      RETURN                                                                 720
      END                                                                    721
      REAL FUNCTION PLEGEN(X,N)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     PLEGEN(X,N) COMPUTES LEGENDREPOLYNOMIALS OF ORDER N FOR ARGUMENT X     723
      COMMON/PLM/XLM,Y(20)
C     Y(N+1)=PL(XLM,N) NACH EDMONDS 
      Y(1)=1.                                                                725
      IF(N.EQ.0) GOTO 1 
      Y(2)=X
      IF(N.EQ.1) GOTO 1 
      DO 4  I=2,N 
      G=X*Y(I)                                                               730
4     Y(I+1)=2.*G-Y(I-1)-(G-Y(I-1))/FLOAT(I)                                 731
1     PLEGEN=Y(N+1) 
      XLM=X 
      RETURN                                                                 733
      END                                                                    734
      REAL FUNCTION YACTLN(N,K)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /COMY/ C(100)                
      DIMENSION   N(100)                 
      B=0.                              
      DO 10 I=1,K                      
      IF(IABS(N(I)).GT.100) GO TO 30  
      IF(N(I).LT.0) GO TO 20         
      J=N(I)+1                      
      B=C(J)+B                     
      GO TO 10                    
   20 J=IABS(N(I))+1             
      B=B-C(J)                  
   10 CONTINUE                 
      YACTLN=B                
      RETURN                 
   30 CONTINUE              
      RETURN               
      END                                                                    477
      REAL FUNCTION YACTRT(N,K)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION N(100)                    
      YACTRT=EXP(0.5*YACTLN(N,K))        
      RETURN                            
      END                              
      REAL FUNCTION YACTOR(N,K)   
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION N(100)                       
      YACTOR=EXP(YACTLN(N,K))               
      RETURN                               
      END                                 
      REAL FUNCTION Y(A,B,C,D,E,F)   
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  NARG(100)                                                   489
C     RACAH COEFFICIENT FUNCTION                                             490
CRACAH     RACAH COEFFICIENT FUNCTION                                        491
C                                                                            492
C     USES THE CLOSED EXPRESSION OF WIGNER.                                  493
C     %SEE ROSE, ELEMENTARY THEORY OF ANGULAR MOMENTUM,%6.7<,P.110<          494
C                                                                            495
      S=0.0                                                                  496
      P=YRI(A,B,E)*YRI(C,D,E)*YRI(A,C,F)*YRI(B,D,F)                          497
      IF(P) 1,1,2                                                            498
    1 Q=0.0                                                                  499
      GO TO 7                                                                500
    2 N1=A+B+E+.1                                                            501
      N2=C+D+E+.1                                                            502
      N3=A+C+F+.1                                                            503
      N4=B+D+F+.1                                                            504
      NMIN=1+MAX0(N1,N2,N3,N4)                                               505
      IF(NMIN)1,1,3                                                          506
    3 M1=A+B+C+D+.1                                                          507
      M2=A+D+E+F+.1                                                          508
      M3=B+C+E+F+.1                                                          509
      NMAX=1+MIN0(M1,M2,M3)                                                  510
      IF(NMAX)1,1,4                                                          511
    4 IF(NMAX-NMIN)1,5,5                                                     512
    5 DO 6 NN=NMIN,NMAX                                                      513
      N=NN-1                                                                 514
      M=M1+N                                                                 515
      T=1-2*MOD(M,2)                                                         516
      NARG(1)=N+1                                                            517
      NARG(2)=N1-N                                                           518
      NARG(3)=N2-N                                                           519
      NARG(4)=N3-N                                                           520
      NARG(5)=N4-N                                                           521
      NARG(6)=N-M1                                                           522
      NARG(7)=N-M2                                                           523
      NARG(8)=N-M3                                                           524
      KARG=8                                                                 525
      S=S+T*YACTOR(NARG,KARG)                                                526
    6 CONTINUE                                                               527
      Q=S*P                                                                  528
    7 Y=Q                                                                    529
      RETURN                                                                 530
      END                                                                    531
      REAL FUNCTION YR(A,B,C) 
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                                                                            533
C     TESTS FOR TRIANGULAR RELATIONSHIP BETWEEN A,B,AND C,AND FOR            534
C     A&B&C AN INTEGER.                  
C                                       
      K=A+B+C+0.1                      
      IF(K)1,5,5                      
    5 IF(MIN(A+B-C,B+C-A,C+A-B)+0.1) 1,2,2  
    1 YR=0.0                               
      GO TO 4                             
    2 M=A+B+C+0.1                        
      N=A+B+C+0.6                       
      IF(N-M)1,3,1                     
    3 YR=1.0                          
    4 RETURN                         
      END                           
      REAL FUNCTION YG(AJ1,AJ2,AJ3,AM1,AM2) 
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  NARG(100)                             
C     CLEBSCH-GORDAN COEFFICIENT  ...  CG%AJ1,AJ2,AJ3,AM1,AM2<               550
CCGCOEF     CLEBSCH-GORDAN COEFFICIENT PROGRAM                               551
C                                                                            552
C                                                                            553
C     USES THE CLOSED EXPRESSION DERIVED BY RACAH                            554
C     %SEE ROSE, ELEMENTARY THEORY OF ANGULAR MOMENTUM,%3.19<,P.40<          555
C
      S=0
      IF(AJ1.LT.ABS(AM1)) GOTO 60
      IF(AJ2.LT.ABS(AM2)) GOTO 60
      T=SQRT(2.0*AJ3+1.0)*YRI(AJ1,AJ2,AJ3)
      IF(T.GT.0.) GOTO 90
   60 CONTINUE
   70 YG=0
      GO TO 410
   90 IF(ABS(AM1)+ABS(AM2).EQ.0.) GOTO 300
      AM3=AM1+AM2
       IF(AJ3.LT.ABS(AM3)) GOTO 60
      K1=AJ1+AM1+.1
      K2=AJ2+AM2+.1
      K3=AJ3+AM3+.1
      L1=AJ1-AM1+.1
      L2=AJ2-AM2+.1
      L3=AJ3-AM3+.1
      I1=L1-L3
      I2=K2-K3
      I3=AJ1+AJ2-AJ3+.1
  200 MIN=MAX0(0,I1,I2)                                                      576
  210 MAX=MIN0(L1,K2,I3)                                                     577
  220 IF(MAX) 60,230,230                                                     578
  230 IF(MAX-MIN) 60,240,240                                                 579
  240 DO 270 N=MIN,MAX                                                       580
  250 A=1-2*MOD(N,2)                                                         581
      NARG(1)=L1                                                             582
      NARG(2)=K2                                                             583
      NARG(3)=K3                                                             584
      NARG(4)=N-L1                                                           585
      NARG(5)=N-K2                                                           586
      NARG(6)=N-I3                                                           587
      NARG(7)=I1-N                                                           588
      NARG(8)=I2-N                                                           589
      NARG(9)=-N                                                             590
      KARG=9                                                                 591
      S=S+A*YACTOR(NARG,KARG)                                                592
  270 CONTINUE                                                               593
      NARG(1)=K1                                                             594
      NARG(2)=-K2                                                            595
      NARG(3)=-K3                                                            596
      NARG(4)=-L1                                                            597
      NARG(5)=L2                                                             598
      NARG(6)=L3                                                             599
      KARG=6                                                                 600
      YG=S*T*YACTRT(NARG,KARG)                                               601
  290 GO TO 410                                                              602
  300 AJ=AJ1+AJ2+AJ3+.1                                                      603
  310 K1=AJ                                                                  604
  320 IF(MOD(K1,2)) 70,330,70                                                605
  330 AJ=AJ*.5                                                               606
  340 K1=AJ                                                                  607
  350 K2=AJ3+AJ                                                              608
  360 L1=AJ1-AJ                                                              609
  370 L2=AJ2-AJ                                                              610
  380 L3=AJ3-AJ                                                              611
  390 A=1-2*MOD(K2,2)                                                        612
      NARG(1)=K1                                                             613
      NARG(2)=L1                                                             614
      NARG(3)=L2                                                             615
      NARG(4)=L3                                                             616
      KARG=4                                                                 617
      YG=A*T*YACTOR(NARG,KARG)                                               618
  410 RETURN                                                                 619
      END                                                                    620
      REAL FUNCTION YRI(A,B,C)                                                    621
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  NARG(100)                                                   622
C     TRIANGLE COEFFICIENT FUNCTION                                          623
CTRIF      TRIANGLE COEFFICIENT FUNCTION                                     624
C                                                                            625
C     %SEE ROSE, ELEMENTARY THEORY OF ANGULAR MOMENTUM, P.111<               626
C                                                                            627
      IF(YR(A,B,C))1,1,2                                                     628
    1 T=0.0                                                                  629
      GO TO 6                                                                630
    2 I1=A+B-C+.1                                                            631
    3 I2=A-B+C+.1                                                            632
    4 I3=-A+B+C+.1                                                           633
      I4=A+B+C+1.1                                                           634
    5 NARG(1)=I1                                                             635
      NARG(2)=I2                                                             636
      NARG(3)=I3                                                             637
      NARG(4)=-I4                                                            638
      KARG=4                                                                 639
      T=YACTRT(NARG,KARG)                                                    640
    6 YRI=T                                                                  641
      RETURN                                                                 642
      END                                                                    643
      REAL FUNCTION X(A1,B1,C1,
     1           A2,B2,C2,                
     2           A3,B3,C3)                                                   646
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     X COEFFICIENT FUNCTION                                                 647
CXCOEFF     X COEFFICIENT FUNCTION    %9-J SYMBOL<                           648
C     USES EXPRESSION OF ROSE INVOLVING SUM OVER GROUPS OF 3 RACAH COEFS     649
C     %SEE ROSE, ELEMENTARY THEORY OF ANGULAR MOMENTUM,P.112<                650
      S=0.0                                                                  651
      P=YR(A1,B1,C1)*YR(A2,B2,C2)*YR(A3,B3,C3)*YR(A1,A2,A3)*YR(B1,B2,B3)     652
     1*YR(C1,C2,C3)                                                          653
      IF(P) 1,1,2                                                            654
    1 Q=0.0                                                                  655
      GO TO 7                                                                656
    2 D=B1+A2+.1                                                             657
      I1=D                                                                   658
      I2=D+.5                                                                659
      IF(I1-I2) 9,8,1                                                        660
    8 D=1.0                                                                  661
      GO TO 10                                                               662
    9 D=.5                                                                   663
   10 M=A1+A2+A3+B1+B2+B3+C1+C2+C3+0.1                                       664
      NMIN=MAX(ABS(B1-A2),ABS(C1-A3),ABS(C2-B3))+D+.1   
      IF(NMIN) 1,1,3                    
    3 NMAX=MIN(B1+A2,C1+A3,C2+B3)+D+.1   
      IF(NMAX) 1,1,4                                                         668
    4 IF(NMAX-NMIN) 1,5,5                                                    669
    5 DO 6 N=NMIN,NMAX                                                       670
      AN=FLOAT(N)-D                                                          671
      S=S+(2.0*AN+1.0)*Y(B1,A2,C1,A3,AN,A1)*Y(A2,B1,C2,B3,AN,B2)*Y(A3,C1     672
     1,B3,C2,AN,C3)                                                          673
    6 CONTINUE                                                               674
      A=1-2*MOD(M,2)                                                         675
      Q=S*A                                                                  676
    7 X=Q                                                                    677
      RETURN                                                                 678
      END                                                                    679

      FUNCTION F6J(JD1,JD2,JD3,LD1,LD2,LD3)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     VERSION I F6J FUNCTION CALLS S6J  FORTRAN IV
C     VEREINFACHT 27.9.95 H.M.H
      J1=JD1
      J2=JD2
      J3=JD3
      L1=LD1
      L2=LD2
      L3=LD3
C     ANGULAR MOMENTUM COUPLING TESTS FOR 6J COEFFICIENT
      F6J=0.0
      IF(J1.LT.0 .OR. J2.LT.0 .OR. J3.LT.0) RETURN
      IF(L1.LT.0 .OR. L2.LT.0 .OR. L3.LT.0) RETURN
      IF(MOD(J1+J2+J3,2).NE.0) RETURN
      IF(J3.GT.J1+J2 .OR. J3.LT.ABS(J1-J2)) RETURN
      IF(MOD(J1+L2+L3,2).NE.0) RETURN
      IF(L3.GT.J1+L2 .OR. L3.LT.ABS(J1-L2)) RETURN
      IF(MOD(L1+J2+L3,2).NE.0) RETURN
      IF(L3.GT.L1+J2 .OR. L3.LT.ABS(L1-J2)) RETURN
      IF(MOD(L1+L2+J3,2).NE.0) RETURN
      IF(J3.GT.L1+L2 .OR. J3.LT.ABS(L1-L2)) RETURN
      F6J=S6J(J1,J2,J3,L1,L2,L3)
      RETURN
      END
      FUNCTION F9J(JD1,JD2,JD3,JD4,JD5,JD6,JD7,JD8,JD9)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     F9J VERSION I  CALLS S6J  FORTRAN IV
C     VEREINFACHT 27.9.95 H.M.H.
      DIMENSION KN(6),KX(6),NN(6)
      J1=JD1   
      J2=JD2  
      J3=JD3 
      J4=JD4
      J5=JD5
      J6=JD6
      J7=JD7
      J8=JD8
      J9=JD9
      F9J= 0.0
C     ANGULAR MOMENTUM COUPLING TESTS FOR 9J COEFFICIENT 
      IF(MOD(J1+J2+J3,2).NE.0) RETURN
      IF(J3.GT.J1+J2 .OR. J3.LT.ABS(J1-J2)) RETURN
      IF(MOD(J4+J5+J6,2).NE.0) RETURN
      IF(J6.GT.J4+J5 .OR. J6.LT.ABS(J4-J5)) RETURN
      IF(MOD(J7+J8+J9,2).NE.0) RETURN
      IF(J9.GT.J7+J8 .OR. J9.LT.ABS(J7-J8)) RETURN
      IF(MOD(J1+J4+J7,2).NE.0) RETURN
      IF(J7.GT.J1+J4 .OR. J7.LT.ABS(J1-J4)) RETURN
      IF(MOD(J2+J5+J8,2).NE.0) RETURN
      IF(J8.GT.J2+J5 .OR. J8.LT.ABS(J2-J5)) RETURN
      IF(MOD(J3+J6+J9,2).NE.0) RETURN
      IF(J9.GT.J3+J6 .OR. J9.LT.ABS(J3-J6)) RETURN
      KN(1)=MAX0(IABS(J2-J6),IABS(J1-J9),IABS(J4-J8))
      KN(2)=MAX0(IABS(J2-J7),IABS(J5-J9),IABS(J4-J3))
      KN(3)=MAX0(IABS(J6-J7),IABS(J5-J1),IABS(J8-J3))
      KN(4)=MAX0(IABS(J6-J1),IABS(J2-J9),IABS(J5-J7))
      KN(5)=MAX0(IABS(J2-J4),IABS(J3-J7),IABS(J6-J8))
      KN(6)=MAX0(IABS(J3-J5),IABS(J1-J8),IABS(J4-J9))
      KX(1)=MIN0(J2+J6,J1+J9,J4+J8)
      KX(2)=MIN0(J2+J7,J5+J9,J4+J3)
      KX(3)=MIN0(J6+J7,J5+J1,J8+J3)
      KX(4)=MIN0(J1+J6,J2+J9,J5+J7)
      KX(5)=MIN0(J2+J4,J3+J7,J6+J8)
      KX(6)=MIN0(J3+J5,J1+J8,J4+J9)
      DO 35 K=1,6
   35 NN(K)=KX(K)-KN(K)
      KSIGN=1
      I=MIN0(NN(1),NN(2),NN(3),NN(4),NN(5),NN(6))
      DO 40 K=1,6
      IF(I-NN(K))40,50,40
   40 CONTINUE
   50 KMIN=KN(K)+1
      KMAX=KX(K)+1
      GO TO(130,52,53,54,55,56),K
   52 J=J1
      J1=J5
      J5=J
      J=J3
      J3=J8
      J8=J
      J=J6
      J6=J7
      J7=J
      GO TO 130
   53 J=J2
      J2=J7
      J7=J
      J=J3
      J3=J4
      J4=J
      J=J5
      J5=J9
      J9=J
      GO TO 130
   54 J=J1
      J1=J2
      J2=J
      J=J4
      J4=J5
      J5=J
      J=J7
      J7=J8
      J8=J
      GO TO 120
   55 J=J1
      J1=J3
      J3=J
      J=J4
      J4=J6
      J6=J
      J=J7
      J7=J9
      J9=J
      GO TO 120
   56 J=J2
      J2=J3
      J3=J
      J=J5
      J5=J6
      J6=J
      J=J8
      J8=J9
      J9=J
  120 KSIGN=(1-MOD(J1+J2+J3+J4+J5+J6+J7+J8+J9,4))
C     SUMMATION OF SERIES OF EQUATION (2)  
  130 SUM=0.0                             
      SIG=(-1)**(KMIN-1)*KSIGN
      FLK=KMIN                           
      DO 200 K=KMIN,KMAX,2              
      TERM=FLK*S6J(J1,J4,J7,J8,J9,K-1)*S6J(J2,J5,J8,J4,K-1,J6)
     1*S6J(J3,J6,J9,K-1,J1,J2)
      FLK=FLK+2.0                    
  200 SUM=SUM+TERM                    
      F9J=SUM*SIG
      RETURN                       
      END                         
      FUNCTION S6J(JD1,JD2,JD3,LD1,LD2,LD3)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     VERSION I  FORTRAN IV
      DIMENSION MA(4),MB(3),MED(12)
      COMMON /FACT/FL(322),NCALL
      DATA NCALL/1/
      J1=JD1
      J2=JD2
      J3=JD3
      L1=LD1
      L2=LD2
      L3=LD3
C     DETERMINE WHETHER TO CALCULATE FL(N) S
      IF(NCALL.EQ.0) GOTO 15
      NCALL=0
C     CALCULATE FL(N) S
      FL(1)=0.0
      FL(2)=0.0
      DO 50 N= 3,322
      FN=N-1
   50 FL(N)=FL(N-1)+LOG(FN)
   15 MED(1)=(-J1+J2+J3)/2
      MED(2)=(+J1-J2+J3)/2
      MED(3)=(+J1+J2-J3)/2
      MED(4)=(-J1+L2+L3)/2
      MED(5)=(+J1-L2+L3)/2
      MED(6)=(+J1+L2-L3)/2
      MED(7)=(-L1+J2+L3)/2
      MED(8)=(+L1-J2+L3)/2
      MED(9)=(+L1+J2-L3)/2
      MED(10)=(-L1+L2+J3)/2
      MED(11)=(+L1-L2+J3)/2
      MED(12)=(+L1+L2-J3)/2
      MA(1)=MED(1)+MED(2)+MED(3)
      MA(2)=MED(4)+MED(5)+MED(6)
      MA(3)=MED(7)+MED(8)+MED(9)
      MA(4)=MED(10)+MED(11)+MED(12)
      MB(1)=MA(1)+MED(12)
      MB(2)=MA(1)+MED(4)
      MB(3)=MA(1)+MED(8)
C     DETERMINE MAXIMUM OF (J1+J2+J3),(J1+L2+L3),(L1+J2+L3),(L1+L2+J3)
      MAX=MA(1)
      DO 30 N=2,4
      IF (MAX-MA(N)) 20,30,30
   20 MAX=MA(N)
   30 CONTINUE
C     DETERMINE MINIMUM OF (J1+J2+L1+L2), (J2+J3+L2+L3),(J3+J1+L3+L1)
      MIN=MB(1)
      DO 51 N=2,3
      IF (MIN-MB(N)) 51,51,40
   40 MIN=MB(N)
   51 CONTINUE
      MINH=MIN
      KMAX=MIN-MAX
      MINP1=MIN+1
      MINI  =MINP1-MA(1)
      MIN2=MINP1-MA(2)
      MIN3=MINP1-MA(3)
      MIN4=MINP1-MA(4)
      MIN5=MINP1+1
      MIN6=MB(1)-MIN
      MIN7=MB(2)-MIN
      MIN8=MB(3)-MIN
      UK=1.E-15
      S=1.0E-15
      IF (KMAX) 65,65,55
   55 DO 60 K=1,KMAX
      UK=-UK*DBLE(MINI-K)*DBLE(MIN2-K)*DBLE(MIN3-K)*DBLE(MIN4-K)/
     1 (DBLE(MIN5-K)*DBLE(MIN6+K)*DBLE(MIN7+K)*DBLE(MIN8+K))
C     CUT OFF SERIES AT 1.0D-25
      IF(ABS(UK)-1.E-25) 65,60,60
   60 S=S+UK
   65 S=S*1.0E+15
C     CALCULATE DELTA FUNCTIONS
      DELOG=0.0
      DO 70 N=1,12
      NUM=MED(N)
   70 DELOG=DELOG+FL(NUM+1)
      NUM1=MA(1)+2
      NUM2=MA(2)+2
      NUM3=MA(3)+2
      NUM4=MA(4)+2
      DELOG=DELOG-FL(NUM1)-FL(NUM2)-FL(NUM3)-FL(NUM4)
      DELOG=0.5*DELOG
      ULOG=FL(MIN5)-FL(MINI)-FL(MIN2)-FL(MIN3)-FL(MIN4)-FL(MIN6+1)-
     1   FL(MIN7+1)-FL(MIN8+1)
      PLOG=DELOG+ULOG
      P=EXP (PLOG)
      S6J=P*S
      IF(MOD(MINH,2).NE.0)  S6J=-S6J
      RETURN
      END      