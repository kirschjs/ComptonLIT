      PROGRAM OBELMA
C    OBERPROGRAMM FUER EINTEILCHENOPERATOREN MIT EINTEILCHEN PUNKTIERTEN DC
C     ELEKTROMAGNETISCHE UEBERGAENGE
C
C
C                       T.M. 1985
C     14.1.93 AUF PARAMETER UMGESTELLT H.M.H
C
C EINGABE                                                                      6
C                                                                              7
C JEDER INDEX I) BEZEICHNET EINE KARTE                                         8
C INTEGERS IM FORMAT 10I3,REALS IM FORMAT E12.4                                9
C                                                                             10
C    1) NBAND2,KBAND            NBAND2 WIRD AN QUAFOR                         11
C          ACHTUNG: NBAND3 = 4 IST EIN WEITERES ZWISCHENBAND]                 12
C    2) STEUERZAHLEN NREG(K),K=1,4
C       K=1,4  NREG(K)=0  OPERATOR K WIRD NICHT GERECHNET
C                     =1  OPERATOR K WIRD GERECHNET                           15
C    3) ZAHL DER MAXIMAL ZU PERMUTIERENDEN TEILCHEN ODER GRUPPEN     KZP      17
C    4) ZAHL DER ZERLEGUNGEN=NZF,ZAHL DER TEILCHEN =NZT          NZF NZT      18
C                                                                             19
C                                                                             20
C DER FOLGENDE DATENSATZ MUSS FUER JEDE ZERLEGUNG WIEDERHOLT WERDEN           21
C    A) ZAHL DER CLUSTER,ZAHL DER SPIN-ISOSPIN-          NZC,NZG,MZG,NOL      22
C       PRODUKTE,ZAHL DER KOMBINATIONEN DER SPIN-ISOSPIN-PRODUKTE,            23
C       NUMMER DER AUSGEZEICHNETEN RELATIVFUNKTION                            24
C    B) ZAHL DER TEILCHEN IM CLUSTER N, N=1,NZC              NGRV(2,N, )      25
C   C1) 1.SPIN-ISOSPIN-PRODUKT                                       NZH      26
C   C2) 2.SPIN-ISOSPIN-PRODUKT                                       NZH      27
C ....                                                               NZH      28
C CNZG) NZG.SPIN-ISOSPIN-PRODUKT                                     NZH      29
C                                                                             30
C     MIT DEN FOLGENDEN COF KOENNEN FESTGELEGT ERDEN                          31
C     S,S3 UND T,T3                                                           32
C                                                                             33
C   D1) KOEFFIZIENT DES  1.  SPIN-ISOSPIN-PROD. ZUR  1. KOMBINATION  COF      34
C       KOEFFIZIENT DES  1.  SPIN-ISOSPIN-PROD. ZUR  2. KOMBINATION  COF      35
C                                                   ...                       36
C       KOEFFIZIENT DES  1.  SPIN-ISOSPIN-PROD. ZUR MZG.KOMBINATION  COF      37
C   D2) KOEFFIZIENT DES  2.  SPIN-ISOSPIN-PROD. ZUR  1. KOMBINATION  COF      38
C       KOEFFIZIENT DES  2.  SPIN-ISOSPIN-PROD. ZUR  2. KOMBINATION  COF      39
C                                                   ...                       40
C       KOEFFIZIENT DES  2.  SPIN-ISOSPIN-PROD. ZUR MZG.KOMBINATION  COF      41
C ....                                                                        42
C DNZG) KOEFFIZIENT DES NZG. SPIN-ISOSPIN-PROD. ZUR  1. KOMBINATION  COF      43
C       KOEFFIZIENT DES NZG. SPIN-ISOSPIN-PROD. ZUR  2. KOMBINATION  COF      44
C                                                   ...                       45
C       KOEFFIZIENT DES NZG. SPIN-ISOSPIN-PROD. ZUR MZG.KOMBINATION  COF      46
C     PROTON  SPIN AUF  = 1                                                   47
C     PROTON  SPIN AB   = 2                                                   48
C     NEUTRON SPIN AUF  = 3                                                   49
C     NEUTRON SPIN AB   = 4                                                   50
C     BEI JEDER FUNKTION MUSS DIE Z-KOMPONENTE DES SPINS MAXIMAL SEIN         51
C     NOL GIBT DIE ZAHL DER CLUSTER IN DER ERSTEN GRUUPE AN                   52
      PARAMETER (NZOPER=4, NZFMAX=10, NZTMAX=12, MZGMAX=9,
     *           NZGMAX=16, NZCMAX=4, NZCMAF=24, NDIM1=3850,
     *           NDIM2=3000, NDIM3=300, NDIM4=6850)
C     NZOPER: ANZAHL DER OPERATOREN
C     NZFMAX: MAXIMALE ANZAHL DER ZERLEGUNGEN
C     NZTMAX:    "       "     "  TEILCHEN
C     MZGMAX:    "       "     "  GEKOPPELTEN   SPIN-ISOSPIN-FUNKTIONEN
C     NZGMAX:    "       "     "  UNGEKOPPELTEN    "      "         "
C     NZCMAX:    "       "     "  CLUSTER
C     NZCMAF: NZCMAX!
C
      PARAMETER (NZTMA1=NZTMAX-1)
C
C                                                                             53
C                                                                             54
      COMMON /PLATZ/ V(NDIM4)
      COMMON/COMY/D(100)                                                      57
      COMMON /CSPIN/ NALG(NZTMAX,3)
C
      COMMON /CUM/ NPERM(NZCMAX+1,NZCMAF),NH(7,4),NKOR(NZTMAX,7),NP(5)
C
      COMMON /CKO/ NBAND2,VEC(NZTMAX,NZTMA1,NZFMAX)
C
      COMMON /CPER/ NZD(NZTMAX),MAXLIQ,MAXLIS,MAXLIT,MZGS,NBAND3,NBAND4
C
      COMMON /CEL/ KBAND,MAXLI,MERKS(NZGMAX,NZGMAX),MG(NZCMAX),
     *             NDEL,NFAK(NZCMAX+1),NFL,NFR,MELAUS
C
      COMMON NGRU(2,NZCMAX,2),NOUT,NZT,NZV,I1,I2,I3,MGW(NZCMAX+1),
     *       MMOEG(NDIM3),MSR(NZCMAX,NZCMAX),MVV(35,NZCMAX),ANSP(NDIM1)
C
      DIMENSION NV(NZCMAX),NREG(NZOPER)
      DIMENSION NS(NZGMAX,NZFMAX),MS(MZGMAX,NZFMAX)
      DIMENSION MB(NZFMAX),MC(NZFMAX),MZG(NZFMAX),NOL(NZFMAX)
      DIMENSION NZC(NZFMAX),NZG(NZFMAX),NSS(2,NZGMAX,NZFMAX)
      DIMENSION NMASSE(2,NZGMAX,NZFMAX),NLAD(2,NZGMAX,NZFMAX)
      DIMENSION MSS(2,MZGMAX,NZFMAX),MMASSE(2,MZGMAX,NZFMAX)
      DIMENSION MLAD(2,MZGMAX,NZFMAX)
      DIMENSION UU(MZGMAX,MZGMAX)
      DIMENSION  ANTP(NDIM1), U(NDIM1), UV(NDIM2)
      DIMENSION COF(MZGMAX,NZGMAX,NZFMAX),NZH(NZTMAX,NZGMAX,NZFMAX)
      DIMENSION NGRV(2,NZCMAX,NZFMAX),NCOF1(MZGMAX,NZGMAX,NZFMAX),
     *         NCOF2(MZGMAX,NZGMAX,NZFMAX)
C
      EQUIVALENCE (V(1),UV(1)), (V(NDIM2+1),U(1))
      EQUIVALENCE (V(1), NCOF1(1,1,1)),
     *             (V(MZGMAX*NZGMAX*NZFMAX+1), NCOF2(1,1,1))
      EQUIVALENCE (V(1),ANTP(1))
C      OPEN(UNIT=3,STATUS='SCRATCH',FORM='UNFORMATTED',
      OPEN(UNIT=3,STATUS='SCRATCH',FORM='UNFORMATTED')
      OPEN(UNIT=4,STATUS='SCRATCH',FORM='UNFORMATTED')
      OPEN(UNIT=5,FILE='INOB',STATUS='OLD')
      OPEN(UNIT=6,FILE='OUTPUT')
      OPEN(UNIT=7,STATUS='SCRATCH',FORM='UNFORMATTED')
      OPEN(UNIT=8,FILE='OBOUT',STATUS='UNKNOWN',FORM='UNFORMATTED')
C
      NIM1 = NDIM1
      NIM2 = NDIM2
      NIM3 = NDIM3
      IF (NDIM4.LT.MAX0(MZGMAX*NZGMAX*NZFMAX*2,NDIM1+NDIM2,
     *    2*NZTMAX*NZTMA1+NZTMA1*(2*NZCMAX-1))) STOP 1
      NOUT=6                                                                  88
      INPUT=5                                                                 89
      NBAND3 = 4                                                              90
  612 FORMAT("0 PROTON-OPERATOREN")
  613 FORMAT("0 NEUTRON-OPERATOREN")
614   FORMAT("0 PROTON-SPIN-OPERATOREN")
615   FORMAT("0 NEUTRON-SPIN-OPERATOREN")
 1000 FORMAT(24I3)                                                            99
 1001 FORMAT(10I3)                                                           100
 1006 FORMAT (41H0 LISTE DER MATRIXELEMENTE ZWISCHEN DER  ,I3,               102
     1          19H TEN STRUKTUR DER  ,I3,20H CLUSTEREINTEILUNG  /           103
     2        30X, 9H UND DER ,I4,19H TEN STRUKTUR DER  ,I3,                 104
     3          20H CLUSTEREINTEILUNG  )                                     105
 1010 FORMAT (1H1)                                                           106
 1011 FORMAT (16H0 DEFINITION DES   ,I3,4H TEN,
     127H SPIN-ISOSPIN-PRODUKTES ZUR    ,I3,4H TEN,
     2                                     18H CLUSTEREINTEILUNG   )         109
 1012 FORMAT(9H  CLUSTER     ,I3,5X,4I3)
 1013 FORMAT  (45H0 DEFINITION DER ZUSAMMENGESETZTEN STRUKTUREN  )           111
 1014 FORMAT (3H0 /,I3,2H ,,I3,4H ) =  /)
 1015 FORMAT(5(F10.4,4H  // ,I3,2H ,,I3,4H )   ))
 1020 FORMAT(22H SPEICHERPLATZ ZU ENG ,3I10)                                 114
2011  FORMAT   (19H0 BERECHNET WERDEN )                                      115
 1090 FORMAT(1H0)                                                            116
 5000 FORMAT (21H0SPEICHERPLATZBEDARF:,3I10)                                 117
C                                                                            118
      D(1)=.0                                                                122
      D(2)=.0                                                                123
      DO 199 I=2,99                                                          124
  199 D(I+1) = ALOG(FLOAT(I))+D(I)                                           125
      READ (INPUT,1000) NBAND2,KBAND,MAIAUS,MELAUS
      READ(INPUT,1000) (NREG(K),K=1,4)
      READ(INPUT,1000) KZP                                                   169
      NFAK(1)=1                                                              170
      DO 300 K=1,KZP                                                         171
  300 NFAK(K+1)=K*NFAK(K)                                                    172
C         NFAK(I) = (I-1)!
C         DEFINITION VON NPERM, NPERM(1,K) =SIGNUM VON PERMUTATION K
C         NPERM(I,K) ENTHAELT DIE PERMUTATION
      NPERM(1,1)=1                                                           174
      NPERM(2,1)=1                                                           175
      DO702 NZP=2,KZP                                                        176
      NBL=NFAK(NZP)                                                          177
      NBH=NZP-1                                                              178
      DO 702 L=1,NZP                                                         179
      NZVL=L-1                                                               180
      DO702 M=1,NBL                                                          181
      ML=(L-1)*NBL+M                                                         182
      NPERM(1,ML)=NPERM(1,M)*(-1)**(L+1)                                     183
      DO 704 NZL=1,NBH                                                       184
      IF(NZL+L.GT.NZP) GOTO 706
      NVL=NZL
      GO TO 704                                                              187
  706 NVL=NZL+1                                                              188
  704 NPERM(NVL+1,ML)=NPERM(NZL+1,M)                                         189
      NZVS=NZP+1-NZVL                                                        190
      NPERM(NZVS,ML)=NZP                                                     191
  702 CONTINUE                                                               192
      READ(INPUT,1000)NZF,NZT                                                194
  714 FORMAT(7H0NZT = ,I3,30HIST ZU GROSS, MAX 12 MOEGLICH    )              195
      IF (NZT.LE.12) GOTO 710
        WRITE(NOUT,714) NZT
      STOP                                                                   199
  710 CONTINUE                                                               200
      NZV = NZT-1                                                            201
C        BLOCK EINGABE
      DO 1 K=1,NZF                                                           202
      READ (INPUT,1001) NZC(K),NZG(K),MZG(K),NOL(K)
      M=NZC(K)                                                               204
      READ (INPUT,1001) (NGRV(2,L,K), L=1,M)                                 205
      NGRV(1,1,K) =1                                                         206
      DO 2 L=2,M                                                             207
    2 NGRV(1,L,K) =NGRV(1,L-1,K)+NGRV(2,L-1,K)                               208
C  NGRV(2,L,.)=ZAHL DER TEILCHEN IM CLUSTER L
C  NGRV(1,L,.)=NUMMER DES 1.TEILCHENS IM CLUSTER L
      KK = NOL(K) + 1                                                        209
      KK = NGRV(1,KK,K)                                                      210
      LM=NZG(K)                                                              211
      DO 3 L=1,LM                                                            212
    3 READ (INPUT,1000) (NZH(N,L,K),N=1,NZT)                                 213
C       EINLESEN DER ELEMENTAREN SPINFUNKTIONEN
      DO 940 L=1,LM                                                          214
C       BESTIMMUNG DER GESAMTSPINS NS, FRAGMENTSPINS NSS
C        FRAGMENTMASSEN NMASSE, FRAGMENTLADUNGEN NLAD
      NS(L,K)=0                                                              215
      NSS(1,L,K) = 0                                                         216
      NSS(2,L,K) = 0                                                         217
      NMASSE(1,L,K) = KK - 1                                                 218
      NMASSE(2,L,K) = NZT - KK + 1                                           219
      NLAD(1,L,K) = 0                                                        220
      NLAD(2,L,K) = 0                                                        221
      DO 940 N=1,NZT                                                         222
      NNS =           2*(NZH(N,L,K)-2*(NZH(N,L,K)/2))-1                      223
      IF(N.GE.KK) GOTO 942
        KZ = 1
      GO TO 943                                                              226
  942 KZ = 2                                                                 227
  943 NLAD(KZ,L,K) = NLAD(KZ,L,K) -(NZH(N,L,K)-1)/2 + 1                      228
      NSS(KZ,L,K) = NSS(KZ,L,K) + NNS                                        229
  940 NS(L,K)=NS(L,K)+NNS                                                    230
      N=MZG(K)                                                               231
      DO 7   L = 1,LM                                                        232
    7 READ (INPUT,1000) (NCOF1(LL,L,K),NCOF2(LL,L,K),LL=1,N)                 233
C      CHECK EINGABE
      DO 950 LL=1,N                                                          234
      I=0                                                                    235
      DO 951 L=1,LM                                                          236
      IF (NCOF1(LL,L,K).EQ.0) GOTO 951
       IF(I.GT.0) GOTO 954
       I=1
C       SETZEN DER WERTE
      MS(LL,K)=NS(L,K)                                                       240
      MMASSE(1,LL,K) = NMASSE(1,L,K)                                         241
      MMASSE(2,LL,K) = NMASSE(2,L,K)                                         242
      MLAD(1,LL,K) = NLAD(1,L,K)                                             243
      MLAD(2,LL,K) = NLAD(2,L,K)                                             244
      MSS(1,LL,K) = IABS(NSS(1,L,K))                                         245
      MSS(2,LL,K) = IABS(NSS(2,L,K))                                         246
      GO TO 951                                                              247
C         CHECK DER WERTE
  954 IF((MS(LL,K)-NS(L,K))+(MMASSE(1,LL,K)- NMASSE(1,L,K))+
     1     (MMASSE(2,LL,K)-NMASSE(2,L,K))+(MLAD(1,LL,K)-NLAD(1,L,K))+
     2     (MLAD(2,LL,K)-NLAD(2,L,K)).EQ.0) GOTO 956
       WRITE(NOUT,1040) L,LL,K
 1040 FORMAT(7H FEHLER,3I3)                                                  252
      STOP                                                                   253
  956 MSS(1,LL,K) = MAX0(MSS(1,LL,K),IABS(NSS(1,L,K)))                       254
      MSS(2,LL,K) = MAX0(MSS(2,LL,K),IABS(NSS(2,L,K)))
951   CONTINUE
  950 CONTINUE                                                               256
C      BERECHNUNG DER COEFFIZIENTEN
      DO 1   L = 1,LM                                                        257
      DO 1   LL = 1,N                                                        258
      COF(LL,L,K) = .0
      IF(NCOF1(LL,L,K).EQ.0) GOTO1
      AB1 = NCOF1(LL,L,K)
      AB2 = NCOF2(LL,L,K)                                                    263
      COF(LL,L,K) = (AB1/ABS(AB1)) * SQRT(ABS(AB1)/AB2)
    1 CONTINUE                                                               266
C      HIERNACH WERDEN NCOF1,NCOF2,NMASSE,NLAD,NSS NICHT MEHR VERWENDET
C      NEUER BLOCK
C     AUSDRUCK DER FUNKTIONSEIGENSCHAFTEN                                    267
      WRITE (NOUT,1010)                                                      268
      DO 6 K=1,NZF                                                           269
      LM=NZG(K)                                                              270
      DO 5 L=1,LM                                                            271
      WRITE (NOUT,1011) L,K                                                  272
      WRITE (NOUT,1090)                                                      273
      MM=NZC(K)                                                              274
      DO 5 M=1,MM                                                            275
      N1=NGRV(1,M,K)                                                         276
      N2=NGRV(2,M,K)+N1-1                                                    277
    5 WRITE (NOUT,1012) M,(NZH(N,L,K),N=N1,N2)                               278
      WRITE (NOUT,1013)                                                      279
      WRITE (NOUT,1090)                                                      280
      LM=MZG(K)                                                              281
      DO 6 L=1,LM                                                            282
      WRITE (NOUT,1014) L,K                                                  283
      MM=NZG(K)                                                              284
      WRITE (NOUT,1015) (COF(L,M,K),M,K,M=1,MM)                              285
    6 WRITE(NOUT,1016)  MMASSE(1,L,K),MMASSE(2,L,K),MLAD(1,L,K),             286
     1        MLAD(2,L,K),MSS(1,L,K),MSS(2,L,K),MS(L,K)                      287
 1016 FORMAT(/30H ZERLEGUNG IN DIE GRUPPEN VON   ,I2,4H UND,I2,              288
     1     10H NUKLEONEN  /4H MIT,I2,4H UND,I2,10H LADUNGEN   /              289
     2     12H SPINS SIND  ,I2,7H /2 UND,I2,3H /2  /                         290
     3   15H GESAMTSPIN IST   ,I2,3H /2    /)                                291
      WRITE (NOUT,1090)                                                      292
      WRITE (NOUT,2011)                                                      293
      IF(NREG(1).NE.0) WRITE(NOUT,612)
      IF(0.NE.NREG( 2)) WRITE (NOUT,613)
      IF(NREG(3).NE.0) WRITE(NOUT,614)
      IF(NREG(4).NE.0) WRITE(NOUT,615)
C       ENDE BLOCK AUSSCHREIBEN
C     NOTIEREN DER FUNKTIONSEIGENSCHAFTEN                                    315
      REWIND NBAND2                                                          316
      WRITE(NBAND2) NZF,NZT,NZV,(NREG(K),K=1,4)
      DO 4 K=1,NZF                                                           318
      M=NZC(K)                                                               319
      WRITE(NBAND2) NZC(K),MZG(K),NOL(K)                                     320
      LL=MZG(K)                                                              321
      WRITE(NBAND2) ((MMASSE(N,L,K),MLAD(N,L,K),MSS(N,L,K),N=1,2),           322
     1     MS(L,K),L=1,LL)                                                   323
    4 WRITE (NBAND2) ((NGRV(N,L,K),L=1,M),N=1,2)                             324
C      HIERNACH WERDEN MLAD,MMASSE,MSS NICHT MEHR VERWENDET
C     MATRIXELEMENTE ZWISCHEN EINFACHEN ALGEBRAISCHEN STRUKTUREN             326
      DO 10 MFL=1,NZF                                                        327
      I1=NZC(MFL)                                                            328
      I3=I1-1                                                                329
      I5=I3-1                                                                330
      DO 12 M=1,I1                                                           332
      NGRU(1,M,1) =NGRV(1,M,MFL)                                             333
   12 NGRU(2,M,1)=NGRV(2,M,MFL)                                              334
      CALL SJACK(I1,NOL(MFL),MFL)                                            335
      MFRR=MFL                                                               336
      DO 10 MFR=1,MFRR                                                       337
      I2=NZC(MFR)                                                            338
      I4=I2-1                                                                339
      WRITE(NBAND2)  ((VEC (M,K,MFR),M=1,NZT),K=1,I4)                        340
      DO 13 M=1,I2                                                           341
      NGRU(1,M,2)=NGRV(1,M,MFR)                                              342
   13 NGRU(2,M,2)=NGRV(2,M,MFR)                                              343
C     ADRESSENSUCHREGISTER
C      HIER WERDEN DIE MOEGLICHEN TEILCHENAUSTAEUSCHE VORBEREITET
      DO 22 K=1,I1                                                           345
      I=0                                                                    346
      NZ=NGRU(2,K,1)                                                         347
      NZZ=(NZ+1)**I2-(NZ+1)**I4
C      NZZ ZAEHLT DIE MOEGLICHKEITEN AB NZ TEILCHEN AUF I2 CLUSTER ZU
C     VERTEILEN, OHNE TEILCHENZAHLERHALTUNG,DIESE WIRD ERST NACH LOOP 24
C     IN ORDNUNG GEBRACHT
      DO 23 L=1,NZZ                                                          350
      J=0                                                                    351
      KK=L                                                                   352
      DO 24 M=1,I2                                                           353
      MM=I2-M+1                                                              354
      NV(MM)=KK/(NZ+1)**(MM-1)                                               355
      KK=KK-NV(MM)*((NZ+1)**(MM-1))                                          356
   24 J=J+NV(MM)                                                             357
      IF(J.NE.NZ) GOTO 23
      DO 26 M=1,I2
      IF(NV(M).GT.NGRU(2,M,2)) GOTO 23
   26 CONTINUE                                                               361
      I=I+1                                                                  362
      MVV(I,K) = L
C      IN MVV(K,L) STEHT DIE NUMMER NZZ,DIE TEILCHEN VON CLUSTER L LINKS
C      AUF DIE CLUSTER RECHTS ZU VERTEILEN. ES GIBT MG(L) MOEGLICHKEITEN
C      DAFUER
   23 CONTINUE                                                               364
   22 MG(K)=I
C      MG(K) ENTHAELT DIE ANZAHL DER TEILCHENAUSTAEUSCHE DES K-TEN
C      CLUSTERS LINKS IN DIE CLUSTER RECHTS
      MGW(1)=1                                                               366
      DO 27 K=1,I1                                                           367
   27 MGW(K+1)=MGW(K)*MG(K)                                                  368
C        UNTERBLOCK ENDE
      MAXLIQ=MGW(I1)
C      MGW(I1) ENTHAELT EINE OBERE ABSCHAETZUNG DER ANZAHL DER DOPPEL-
C      NEBENKLASSEN SYMBOLE, NAEMLICH DIE PRODUKTE DER ZEILENAUSTAEUSCHE
      N1 = 1                                                                 370
       N2 = MAXLIQ                                                           371
      N3 = NIM3                                                              372
      WRITE (NOUT,5000) N1,N2,N3                                             373
      IF(N2.LE.N3) GOTO 910
 1021 WRITE(NOUT,1020) N1,N2,N3                                              376
      STOP                                                                   377
  910 CONTINUE                                                               378
      DO 48 K=1,MAXLIQ                                                       379
   48 MMOEG(K)=0                                                             380
      I=0                                                                    381
      KK=1                                                                   382
   28 DO 29 K=KK,I1                                                          383
   29 NV(K)=1                                                                384
      DO 30 K=1,I1                                                           385
      M=NV(K)                                                                386
   30 MB(K)=MVV(M,K)                                                         387
      DO 31 L=1,I2                                                           388
      L1=I2-L+1                                                              389
      M=0                                                                    390
      DO 32 K=1,I1                                                           391
      MC(K)=MB(K)/((NGRU(2,K,1)+1)**(I2-L))                                  392
      MB(K)=MB(K)-MC(K)*((NGRU(2,K,1)+1)**(I2-L))                            393
   32 M=M+MC(K)                                                              394
      IF(M.NE.NGRU(2,L1,2)) GOTO 40
   31 CONTINUE                                                               396
      I=I+1                                                                  397
      J=NV(1)                                                                398
      IF(I5.LE.0) GOTO 1401
      DO 34    K = 1,I5
   34 J=J+(NV(K+1)-1)*MGW(K+1)                                               401
 1401 MMOEG(J) = I                                                           402
   40 DO 41 K=1,I1                                                           404
      KK=I1-K+1                                                              405
      IF(NV(KK).LT.MG(KK)) GOTO 42
   41 CONTINUE                                                               407
      GO TO 44                                                               408
   42 NV(KK)=NV(KK)+1                                                        409
      KK=KK+1                                                                410
      GO TO 28
C      UNTERBLOCK ENDE
   44 MAXLIT=I
C      MAXLIT IST DIE ANZAHL DER DOPPELNEBENKLASSEN
C    MMOEG(J)=N BEDEUTET,DAS DC-SYMBOL N HAT DIE ADRESSE J,WOBEI
C     J FESTLEGT WELCHE SPALTENAUSTAEUSCHE VORZUNEHMEN SIND VON LINKS
C     NACH RECHTS. J=NV(1)+SUMME (NV(K+1)-1)*MGW(K+1),K=1,NZCL-2
      J=1                                                                    413
      DO 46 K=1,I1                                                           414
      DO 46 L=1,I2                                                           415
      MSR(L,K)=J                                                             416
C  ADDRESSEN DER DC-SYMBOLMATRIXELEMENTE ZUM KENNZEICHNEN DER WECHSELW
   46 J=J+1
      M=I1*I2
      MAXLIS=M
C      MAXLIS =ANZAHL DER MOEGLICHEN ET-WW IN EINEM DC
      MAXLI=MAXLIT*MAXLIS
C     MAXLI= ANZAHL DER PUNKTIERTEN DC
      N3 = NIM1                                                              428
      N1 = 2                                                                 429
      WRITE (NOUT,5000) N1,MAXLI,N3                                          430
      IF(MAXLI.LE.N3) GOTO 1026
      N2=MAXLI                                                               433
      GO TO 1021                                                             434
1026  NZGL=NZG(MFL)
      NZGR=NZG(MFR)                                                          437
      MZGL=MZG(MFL)                                                          438
      MZGR=MZG(MFR)                                                          439
      MZGS = MZGR*MZGL                                                       440
      DO 50   MKC = 1,4
C        LOOP UEBER OPERATOREN
      IF(0.EQ.NREG(MKC))GO TO 50                                             442
      DO 52   K= 1,MAXLI                                                     443
52    ANTP(K)=0.
      GO TO(602,603,604,605),MKC
  602 WRITE(NOUT,612)                                                        446
      GO TO 4022                                                             447
  603 WRITE(NOUT,613)                                                        448
      GOTO 4022
604   WRITE(NOUT,614)
      GOTO 4022
605   WRITE(NOUT,615)
 4022 CONTINUE                                                               459
      REWIND KBAND                                                           460
      DO 60 NFL=1,NZGL                                                       462
C         LOOP ELEMENTARE SPINFUNKTIONEN LINKS
      DO 64 K=1,NZT                                                          463
   64 NALG(K,1)=NZH(K,NFL,MFL)                                               464
      DO 60 NFR=1,NZGR                                                       465
C         LOOP ELEMENTARE SPINFUNKTIONEN RECHTS
      DO 65 K=1,NZT                                                          466
   65 NALG(K,2)=NZH(K,NFR,MFR)                                               467
      NDEL = NS(NFL,MFL) - NS(NFR,MFR)                                       468
      NDEL = NDEL/2                                                          469
      CALL ELEML(MKC)                                                        470
60    CONTINUE                                                               471
C         ELEMENTARE LISTEN ERSTELLT
      REWIND NBAND3                                                          472
      WRITE (NBAND3) (ANTP(K),K=1,MAXLI)
      CALL PERME (MKC,1,NNN)
C      MIT DIESEM AUFRUF WERDEN DIE NNN REPRESENTANTEN FUER WW INDC AUF
C      NBAND3 GESCHRIEBEN
      N1 = 3                                                                 476
      N2 = NNN*MZGL*MZGR                                                     477
      N3 = NIM2                                                              478
      WRITE (NOUT,5000) N1,N2,N3                                             479
      IF (N2.GT.N3) GOTO 1021
      REWIND NBAND3                                                          482
      WRITE(NBAND2)    NNN,NNN,NNN                                           483
      IF(NNN.EQ.0) GOTO 700
C     MATRIXELEMENTE ZWISCHEN ZUSAMMENGESETZTEN STRUKTUREN                   486
      DO 130 KFL=1,MZGL                                                      487
      KFRR=MZGR                                                              488
      DO 130 KFR=1,KFRR                                                      489
      NDEL=MS(KFL,MFL)-MS(KFR,MFR)                                           490
      NDEL = NDEL/2                                                          491
C                                                                            492
C                                                                            493
      SL=MS(KFL,MFL)                                                         494
      SML=.5*SL                                                              495
      SR=MS(KFR,MFR)                                                         496
      SMR=.5*SR                                                              497
      GOTO(419,419,420,420),MKC
419   X=0.
      F=1.                                                                   501
      IIX=0
      IIF=1      
      GOTO 422
420   X=1.
      F=1.
      IIX=1
      IIF=1
422   AMQ=NDEL                                                               512
      FAKTOR=0.
C      Y=YG(SMR,X,SML,SMR,AMQ)                                                513
C      write(nout,*)'clgtest:',SMR,X,SML,SMR,AMQ,Y
      Y=CLG(MS(KFR,MFR),2*IIX,MS(KFL,MFL),
     1       MS(KFR,MFR),MS(KFL,MFL)-MS(KFR,MFR))
c      write(nout,*)'clgtest:',MS(KFR,MFR),2*IIX,MS(KFL,MFL),
c     1       MS(KFR,MFR),MS(KFL,MFL)-MS(KFR,MFR),Y2
      IF(Y.EQ.0.) GOTO 8913
      F = F * SQRT(2.*SML+1.)                                                518
      FAKTOR=F/Y                                                             519
C                                                                            521
8913  IF(MAIAUS.EQ.0)GO TO 121
      WRITE(NOUT,1006) KFL,MFL,KFR,MFR                                       524
       WRITE (NOUT,1090)                                                     525
121    NBAND4 = (KFL - 1)*MZGR + KFR
      REWIND KBAND                                                           528
      DO 132 K=1,MAXLI                                                       529
132   U(K)=.0                                                                530
      DO 134 NFL=1,NZGL                                                      531
      DO 134 NFR=1,NZGR                                                      532
      IF(MERKS(NFL,NFR).EQ.0) GOTO 134
       READ(KBAND) (ANSP(KKX),KKX=1,MAXLI)
      DO 140 K=1,MAXLI
      AKN = ANSP(K)
      U(K) = U(K) + COF(KFL,NFL,MFL)*COF(KFR,NFR,MFR)* AKN*FAKTOR
  140 CONTINUE                                                               540
134   CONTINUE                                                               541
      REWIND NBAND3                                                          542
      READ (NBAND3) (ANSP(K),K=1,MAXLI)
C     NOTIEREN DER LISTEN ZUSAMMENGESETZTER STRUKTUREN                       544
      IF (MAIAUS.EQ.0) GOTO 126
      CALL PERME(MKC,0,NNK)
126   CALL PERME(MKC,-2,NNK)
130   CONTINUE                                                               554
      I = 0                                                                  555
      DO 810    K = 1,NNN                                                    556
      IK = 0                                                                 557
      DO 800   KFL = 1,MZGL                                                  558
      II = IK + I                                                            559
      IK = IK + MZGR                                                         560
      DO 800   KFR = 1,MZGR                                                  561
      JK = II + KFR                                                          562
      UU(KFL,KFR) = UV(JK)
  800 CONTINUE                                                               565
      I = I + MZGS                                                           566
      READ (NBAND3) N1, (NZD(I6),I6=1,NZT)
      WRITE(NBAND2) N1, (NZD(I6),I6=1,NZT) ,
     1    ((UU(KFL,KFR),KFL=1,MZGL),KFR=1,MZGR)
      IF(MAIAUS.EQ.0)GOTO 810
      WRITE(NOUT,2012)N1,(NZD(I6),I6=1,NZT)
2012  FORMAT(1X,20I3)
      WRITE(NOUT,2013)((UU(KFL,KFR),KFL=1,MZGL),KFR=1,MZGR)
2013  FORMAT(1X,1P10E12.4)
810   CONTINUE
  700 CONTINUE                                                               570
C      ENDE LOOP OPERATOREN
50    CONTINUE
   10 CONTINUE                                                               572
      WRITE(NOUT,3010)                                                       573
 3010 FORMAT(//18H ENDE DER RECHNUNG )                                       574
      STOP                                                                   575
      END                                                                    576
      SUBROUTINE ELEML(MKC)                                                  577
C     SUBROUTINE ZUR AUFSTELLUNG DER ELEMENTAREN LISTEN                      578
C     WIRD FUER JEDE ZERL. LI. U. RE., JEDEN OPERATOR, JEDES ELEMENTARE
C     SPIN-ISOSPIN-PRODUKT AUFGERUFEN
      PARAMETER (NZFMAX=10, NZTMAX=12, NZGMAX=16, NZCMAX=4,
     *           NZCMAF=24, NDIM1=3850, NDIM3=300, NDIM4=6850)
C
      COMMON /PLATZ/ V(NDIM4)
C
      COMMON /CSPIN/ NALG(NZTMAX,3)
C
      COMMON /CUM/ NPERM(NZCMAX+1,NZCMAF),NH(7,4),NKOR(NZTMAX,7),NP(5)
C
      COMMON /CEL/ KBAND,MAXLI,MERKS(NZGMAX,NZGMAX),MG(NZCMAX),
     *             NDEL,NFAK(NZCMAX+1),NFL,NFR,MELAUS
C
      COMMON NGRU(2,NZCMAX,2),NOUT,NZT,NZV,I1,I2,I3,MGW(NZCMAX+1),
     *       MMOEG(NDIM3),MSR(NZCMAX,NZCMAX),MVV(35,NZCMAX),ANSP(NDIM1)
C
      DIMENSION NV(NZCMAX),NSG(NZCMAX,NZCMAX)
      DIMENSION  ANTP(NDIM1)
      DIMENSION NW(NZCMAX)
      EQUIVALENCE (V(1),ANTP(1))
      DO 320 K=1,4                                                           600
      I=0                                                                    601
      DO 321 L=1,NZT                                                         602
      IF(NALG(L,1).NE.K) GOTO 321
          I=I+1
      NH(I+2,K)=L                                                            605
  321 NH(1,K)=NFAK(I+1)                                                      606
      NH(2,K)=I                                                              607
  320 CONTINUE                                                               608
C      NH(2,K) =  I =ANZAHL DER TEILCHEN VOM TYP K,  NH(1,K)= I!
C       NH(2+...I,K)= NUMMMER DES TEILCHEN VOM TYP K
      GOTO(301,301,305,305),MKC
301   LKC=1
      GOTO 20
305   LKC=3
20    DO 1 MLKC=1,LKC
      IF(MLKC.GT.1)GOTO 22
      DO 21 K=1,MAXLI
21    ANSP(K) = 0.
      ISPZ1 = 0                                                              627
22    CONTINUE
C      CHECK OB MATRIXELEMENT FUER SPINDIFFERENZ MOEGLICH IST
      IF(MKC.GT.2)GOTO 200
      FAK=1.
      IF(NDEL.NE.0) GOTO 90
      GOTO 220
200   IF(IABS(NDEL).GT.1) GOTO 90
      IF(NDEL)207,208,209
207   IF(MLKC-1)90,220,90
208   IF(MLKC-2)90,220,90
209   IF(MLKC-3)90,220,90
220   ISPZ1 = ISPZ1+1
         LT = 1                                                              665
      IF(MKC.LT.3) GOTO 253
      IF(MLKC-2) 250,251,252
250   FAK=1./SQRT(2.)
      GOTO 253
251   FAK=.5
      GOTO 253
252   FAK=-1./SQRT(2.)
253   CONTINUE
510   DO 29 K=1,NZT
      NKOR(K,2)=K                                                            671
   29 NALG(K,3)=NALG(K,2)                                                    672
C      AN DIESER STELLE IST NALG(.,3) DIE SPINFUNKTION DER RECHTEN SEITE
C     JETZT ANWENDUNG DER SPINOPERATOREN
      IFT=1
      IF(MKC.LT.3)GOTO 100
      IF(MLKC-2)46,47,48
46    CALL SPINDO(LT,I)
      IF(I) 161,161,100
47    CALL SPINWE(LT,I)
      IFT=I
      GOTO 100
48    CALL SPINUP(LT,I)
      IF(I.EQ.0) GOTO 161
C      NALG(.,3) IST DIE SPINFUNKTION DER RECHTEN SEITE NACH ANWENDUNG
C       DER SPINOPERATOREN
C     AUFSUCHEN DER PERMUTATION P0                                           724
100   DO 101 K=1,NZT
      NKOR(K,1)=K                                                            726
  101 NKOR(K,3)=K                                                            727
      I=1                                                                    728
      DO 104 K=1,NZT                                                         729
      IF(NALG(K,1).EQ.NALG(K,3)) GOTO 104
       K1=K+1
      DO 106 L=K1,NZT                                                        732
      IF(NALG(K,1).EQ.NALG(L,3)) GOTO 107
  106 CONTINUE                                                               734
      GO TO 161
C     SPINMATRIXELEMENT =0, NEUES PAAR WW-TEILCHEN
  107 NALG(L,3)=NALG(K,3)                                                    736
      NALG(K,3)=NALG(K,1)                                                    737
      J=NKOR(L,3)                                                            738
      NKOR(L,3)=NKOR(K,3)                                                    739
      NKOR(K,3)=J                                                            740
      I=-I                                                                   741
  104 CONTINUE                                                               742
      NP(1)=I                                                                743
C      NP(1) ENTHAELT VORZEICHEN DER PERMUTATION P0,
C      NKOR(.,3) ENTHAELT DIE PERMUTATION P0, ANGEWANDT AUF DIE FUNKTIONEN
      DO 108 K=1,NZT                                                         744
      L=NKOR(K,3)                                                            745
  108 NKOR(L,4) =K
C      NKOR(.,4) INVERSE PERMUTATION ZU NKOR(.,3), P0**-1
      DO 110 K=1,NZT                                                         747
      L=NKOR(K,2)                                                            748
  110 NKOR(K,5)=NKOR(L,4)
C    NKOR(.,5) PERMUTATION P0**-1 * NKOR(.,2)
C    NKOR(.,2) IST HIER DIE IDENTITAET
      LT1=NKOR(LT,4)                                                         750
C      CHECK OB WW-TEILCHEN PROTON ODER NEUTRON
      GOTO(114,115,114,115),MKC
114   IF(NALG(LT1,1)-3) 117,161,161
115   IF(NALG(LT1,1)-3) 161,117,117
  117 MG1=NH(1,1)                                                            755
      DO 120 NF1=1,MG1                                                       756
      CALL UM(NF1,1)                                                         757
      MG2=NH(1,2)                                                            758
      DO 120 NF2=1,MG2                                                       759
      CALL UM(NF2,2)                                                         760
      MG3=NH(1,3)                                                            761
      DO 120 NF3=1,MG3                                                       762
      CALL UM(NF3,3)                                                         763
      MG4=NH(1,4)                                                            764
      DO 120 NF4=1,MG4                                                       765
      CALL UM(NF4,4)
C     DURCH DIE AUFRUFE VON UM WERDEN ALLE PERMUTAIONEN P' IN NKOR(.,1)
C     ERZEUGT
      DO 124 K=1,NZT                                                         767
      DO 124 L=1,2                                                           768
      M=NKOR(K,L+3)                                                          769
      N=NKOR(M,1)                                                            770
  124 NKOR(K,L+5)=N                                                          771
      LT2=NKOR(LT1,1)                                                        772
C      NKOR(.,6)=P'*(P0**-1),NKOR(.,6) WIRD IM FOLGENDEN NICHT VERWENDET
C      NKOR(.,7)=P'*(P0**-1)*NKOR(.,2)
C      HIER IST NKOR(.,6)=NKOR(.,7)
C      NP(5) IST VORZEICHEN DER GESAMTPERMUTATION
C     KONSTRUKTION DES PUNKTIERTEN DC-SYMBOLS
      DO 130 K=1,I1
      K1=NGRU(1,K,1)                                                         776
      K2=NGRU(2,K,1)+K1-1                                                    777
      DO 130 L=1,I2                                                          778
      L1=NGRU(1,L,2)                                                         779
      L2=NGRU(2,L,2)+L1-1                                                    780
      NSG(L,K)=0
      DO 132 M=K1,K2
      DO 1132 N=L1,L2
      IF(NKOR(N,7).NE.M)GOTO 135
      NSG(L,K)=NSG(L,K) +1
C      NSG(L,K) ENTHAELT DAS DC-SYMBOL FUER DAS L-TE CLUSTER RECHTS
C      UND K-TE CLUSTER LINKS
  135 IF(NKOR(N,7).EQ.LT2) MS21=L
1132   CONTINUE
      IF(M.EQ.LT2) MS11=K
  132 CONTINUE                                                               795
  130 CONTINUE                                                               796
C      ADRESSRECHNUNG
      DO 150 K=1,I1                                                          853
      NW(K)=0                                                                854
      DO 151 L=1,I2                                                          855
  151 NW(K)=NW(K)+NSG(L,K)*(NGRU(2,K,1)+1)**(L-1)                          856
      MM=MG(K)                                                               857
      DO 154 M=1,MM                                                          858
      IF(NW(K).EQ.MVV(M,K)) GOTO 150
  154 CONTINUE                                                               860
  150 NV(K)=M
C     NV ENTHAELT DIE JEWEILIGE NR DES AUSTAUSCHES
      J=1                                                                    862
      DO 152 K=1,I3                                                          863
  152 J=J+(NV(K)-1)*MGW(K)
C     J = ADRESSE DES DC-SYMBOLS,WIE IN HAUPT BEI 1401 BERECHNET
      J=MMOEG(J)
C     J=NUMMER DES DC-SYMBOLS
      J=I1*I2*(J-1)
      M=MSR(MS21,MS11)
C    KENNZEICHNUNG DER WECHSELWIRKUNG IM DC
C      DOPPELNEBENKLASSENINDEX VON TEILCHEN 1
      J=J+M
C      J=ADRESSE DES PUNKTIERTEN ET-DC,SO DASS ALLE MOEGLICHEN WW FUER EIN
C      DC NACHEINANDER ERFOLGEN
      ANSP(J) = ANSP(J) + FAK*FLOAT(IFT*NP(5))
      IF(MELAUS.LT.2) GOTO 120
      WRITE(NOUT,126) MLKC,J,M,IFT,ANSP(J),NP(5),LT,LT1,LT2,
     $(NKOR(I,7),I=1,NZT),(NKOR(I,4),I=1,NZT),(NKOR(I,1),I=1,NZT)
126   FORMAT(1X,2I6,I3,I6,F8.4,I4,'  ',3I3,'   ',36I2)
  120 CONTINUE                                                               871
  161 LT=LT+1                                                                877
      IF(LT.LE.NZT) GOTO 510
   90 CONTINUE                                                               880
C     ADDITION DER LISTEN FUER GLEICHE SPINDIFFERENZEN                       884
      IF(MLKC.NE.LKC) GOTO 1
      MERKS(NFL,NFR) = 0
      IF (ISPZ1.EQ.0) GOTO 1
      AJJ = 0
      DO 53    I = 1,MAXLI
C      LOOP UEBER ALLE PUNKTIERTEN DC
      AIW3 = ABS(ANSP(I))
      AJJ = AJJ + AIW3
      ANTP(I) = ANTP(I) + AIW3
   53 CONTINUE                                                               894
      IF(AJJ.EQ.0) GOTO 1
      MERKS(NFL,NFR) = 1
      WRITE(KBAND) (ANSP(I),I=1,MAXLI)
      IF (MELAUS.EQ.0) GOTO 1
      WRITE(NOUT,600)  NFL,NFR
      CALL PERME (MKC,-1,NNN)
  600 FORMAT(20H ELEMENTARE LISTE            ,4I3/)
    1 CONTINUE                                                               903
      RETURN                                                                 904
      END                                                                    905
      SUBROUTINE PERME(MKC,KENN,NZAHL)
C      PERME WIRD GERUFEN MIT KENN=-1 VON ELEML ZUM AUSDRUCK FUER JEDE
C      ELEMENTARE SPINFUNKTION
C      PERME WIRD GERUFEN MIT KENN=1 VON HAUPT IMMER,SCHREIBT REPRESENTANTEN
C      FUER PUNKTIERTES DC AUF NBAND3
C      PERME WIRD GERUFEN MIT KENN=0 VON HAUPT FALLS AUSDRUCK FUER
C      GEKOPPELTE SPINFUNKTIONEN
C      PERME WIRD GERUFEN MIT KENN =-2 VON HAUPT IMMER,SPEICHERT DIE
C      MATRIXELEMENTE FUER GEKOPPELTE SPINFUNKTIONEN IN  UV
      PARAMETER (NZTMAX=12, NZCMAX=4, NDIM1=3850,
     *           NDIM2=3000, NDIM3=300, NDIM4=6850)
C
      COMMON /PLATZ/ V(NDIM4)
C
      COMMON /CPER/ NZD(NZTMAX),MAXLIQ,MAXLIS,MAXLIT,MZGS,NBAND3,NBAND4
C      MAXLIT=ANZAHL DER DC
C      MAXLIQ=ANZAHL DER ZEILENAUSTAUSCHE >= MAXLIT
C      MAXLIS=ANZAHL DER MOEGLICHEN 2-TEILCHEN WW IN EINEM DC
C
      COMMON /CADD/ NSH(4,NZCMAX,NZCMAX)
C
      COMMON NGRU(2,NZCMAX,2),NOUT,NZT,NZV,I1,I2,I3,MGW(NZCMAX+1),
     *       MMOEG(NDIM3),MSR(NZCMAX,NZCMAX),MVV(35,NZCMAX),ANSP(NDIM1)
C
      DIMENSION ANTP(NDIM1), U(NDIM1), UV(NDIM2)
C
      EQUIVALENCE (V(1),UV(1)), (V(NDIM2+1),U(1))
      EQUIVALENCE (V(1),ANTP(1))
      NZAHL=0                                                                927
      DO 1 K=1,MAXLIT
C     LOOP UEBER DC
      NAD=(K-1)*MAXLIS                                                       929
      A=.0                                                                   930
      DO 2 L=1,MAXLIS
C      LOOP UEBER WW-TERME
      N=NAD+L                                                                932
      IF(KENN) 5,3,8                                                         933
8     AKNM = ANTP(N)
      A=A+ AKNM
      GO TO 2                                                                936
    3 A = A +  ABS(U(N))                                                     937
      GO TO2                                                                 938
5     A=A+ ANSP(N)**2
    2 CONTINUE                                                               940
      IF(A.LT.0.00001) GOTO 1
C     FALLS KEIN BEITRAG FUER DIESES DC
      CALL DCADD(K)
      DO 30 K1=1,I1                                                          979
      DO 30 L1=1,I2                                                          980
      M1=MSR(L1,K1)
C     DC-INDEX TEILCHEN 1
       NADD=NAD+M1
C      NADD= GESAMTADRESSE
      N1=NSH(1,L1,K1)
C     N1 REPRESENTANT FUER WW TEILCHEN 1
      IF(KENN) 50,51,80
80    IF(ANTP(NADD).EQ.0) GOTO 30
       NZAHL=NZAHL+1
      WRITE(NBAND3) N1, (NZD(I6),I6=1,NZT)
      GO TO 30
51    IF(U(NADD).EQ.0.) GOTO 30
       WRITE (NOUT,100) U(NADD),MKC,N1,(NZD(I),I=1,NZT)
      GO TO 30                                                              1001
  100 FORMAT(1PE12.4,4H (/O  ,I1,2H ( ,I2,3H )/ ,20I3,2H )  /)
140   WRITE(NOUT,101) ANSP(NADD),MKC,N1,(NZD(I),I=1,NZT)
      GO TO 30
50    IF(ANSP(NADD).EQ.0) GOTO 30
             IF((KENN+1).EQ.0) GO TO 140
      IF( ABS(U(NADD)).GT.0.00001)     GO TO 144                            1006
      U(NADD)=.0                                                            1007
  144 UV(NBAND4) = U(NADD)
      NBAND4 = NBAND4 + MZGS                                                1009
101   FORMAT(F12.4  ,4H (/O  ,I1,2H ( ,I2,3H )/ ,20I3,2H )  /)
   30 CONTINUE                                                              1011
    1 CONTINUE                                                              1012
      RETURN                                                                1013
      END                                                                   1014
      SUBROUTINE DCADD(K)
C     DCADD BERECHNET AUS DER ADRESSE DEN ENTSPRECHENDEN DOPPELNEBEN
C     KLASSENREPRESENTANTEN ZURUECK
C     DIESER STEHT IN DC-FORM IN NSH,BZW ALS PERMUTATION IN NZD
C
      PARAMETER (NZTMAX=12, NZCMAX=4, NDIM3=300)
C
      COMMON /CPER/ NZD(NZTMAX),MAXLIQ,MAXLIS,MAXLIT,MZGS,NBAND3,NBAND4
C
      COMMON /CADD/ NSH(4,NZCMAX,NZCMAX)
C
      COMMON NGRU(2,NZCMAX,2),NOUT,NZT,NZV,I1,I2,I3,MGW(NZCMAX+1),
     *       MMOEG(NDIM3),MSR(NZCMAX,NZCMAX),MVV(35,NZCMAX)
C
      DIMENSION NV(NZCMAX),NSG(NZCMAX,NZCMAX)
C
C       KONSTRUKTION EINES REPRESENTANTEN DES PUNKTIERTEN DC-SYMBOLS
       DO 6 L=1,MAXLIQ
      IF(K.EQ.MMOEG(L)) GOTO 10
    6 CONTINUE
   10 M=L
C     M=ADRESSE DES K-TEN DC-SYMBOLS
C     KONSTRUKTION DES DC-SYMBOLS
      DO 12 N=2,I1                                                           946
      NN=I1-N+1                                                              947
      NV(NN)=(M-1)/MGW(NN)+1                                                 948
      M=M-(NV(NN)-1)*MGW(NN)                                                 949
      I=NV(NN)                                                               950
      II=MVV(I,NN)                                                           951
      DO 16 I=1,I2                                                           952
      IJ=I2-I+1                                                              953
      NSG(IJ,NN)=II/(NGRU(2,NN,1)+1)**(IJ-1)                               954
   16 II=II-NSG(IJ,NN)*(NGRU(2,NN,1)+1)**(IJ-1)                            955
   12 CONTINUE                                                               956
      DO 18 M=1,I2                                                           957
      NSG(M,I1)=NGRU(2,M,2)                                                958
      DO 18 N=1,I3                                                           959
   18 NSG(M,I1)=NSG(M,I1)-NSG(M,N)
C     NSG(L,K) ENHAELT DAS DC,L-TES CLUSTER RECHTS,K-TES CLUSTER LINKS
C     BESTIMMEN DES REPRESENTANTEN DES DC-SYMBOLS
      I=1                                                                    961
      DO 20 J=1,I1                                                           962
      DO 20 L=1,I2                                                           963
      M=NSG(L,J)                                                           964
      IF(M.LE.0) GOTO 20
       DO 21   N = 1,M
      NSH(N,L,J)=I                                                           967
   21 I=I+1                                                                  968
   20 CONTINUE                                                               969
      I=1                                                                    970
      DO 24 L=1,I2                                                           971
      DO 24 J=1,I1                                                           972
      M=NSG(L,J)                                                           973
      IF(M.LE.0) GOTO 24
       DO 25   N = 1,M
      NZD(I)=NSH(N,L,J)                                                      976
   25 I=I+1                                                                  977
   24 CONTINUE
C     REPRESENTANT BESTIMMT
      RETURN
       END
      SUBROUTINE SJACK(K1,LL,MFL)                                           1015
C     DIESE SUBROUTINE BERECHNET DIE JACOBI KOORDINATEN UND SCHREIBT
C     SIE AUF BAND
C     EINGABE !  K1 = ZAHL DER CLUSTER, MFL= NR DER ZERLEGUNG
C     LL = ZAHL DER CLUSTER IM ERSTEN FRAGMENT ,WIRD UM 1 ERHOEHT,ABER
C     NOL WIRD ANSCHLIESEND NICHT MEHR VERWENDET
C     SVEC ENTHAELT DIE JACOBIKOORDINATEN ALS FUNKTION DER EINTEILCHEKOORDINATEN
C     DIE ORTHOGONALE TRANSFORMATION WIRD OHNE!!! DEN SCHWERPUNKT AUSGEFUEHRT
C     RVEC ENTHAELT DIE DAZU TRANSPONIERTE MATRIX
C     S(I)=SUMME UEBER J (RVEC(I,J)*R(J))
C     C(.,.,K) =1, KENNZEICHNET DIE INNEREN JACOBIKOOR. VON CLUSTER K
C      C(.,.,K) =1 FUER K=NZC CHARAKTERISIERT DIE RELATIV KOORDINATEN
      PARAMETER (NZFMAX=10, NZTMAX=12, NZCMAX=4, NDIM4=6850)
      PARAMETER (NZTMA1=NZTMAX-1)
C
      COMMON /PLATZ/ V(NDIM4)
C
      COMMON /CKO/ NBAND2,VEC(NZTMAX,NZTMA1,NZFMAX)
C
      COMMON NGRU(2,NZCMAX,2),NOUT,NZT,NZV
C
      DIMENSION RVEC(NZTMA1,NZTMAX),SVEC(NZTMAX,NZTMA1),
     *          C(NZTMA1,2*NZCMAX-1)
      EQUIVALENCE (V(1),RVEC(1,1)),(V(NZTMAX*NZTMA1+1),SVEC(1,1)),
     *             (V(2*NZTMAX*NZTMA1+1),C(1,1))
C
      LL = LL + 1                                                           1032
      L5 = NGRU(1,LL,1)                                                     1033
      K2=K1-1                                                               1034
      K4=K1 + K2                                                            1035
      DO 2124 K=1,NZV
      DO 2124 L=1,NZT
2124  SVEC(L,K)=0.
      DO 2125 L=1,NZV
      DO 2125 I=1,K4
2125  C(L,I)=0.
      I = 0                                                                 1043
      DO 10   K = 1,K1                                                      1044
      L4 = NGRU(2,K,1) - 1                                                  1045
      IF(L4.EQ.0) GOTO 10
       DO 12   L = 1,L4
      I = I + 1                                                             1048
      DO 13   M = 1,L                                                       1049
      L1 =      NGRU(1,K,1) + M - 1                                         1050
   13 SVEC(L1,I) = -1./ FLOAT(L)                                            1051
      L1=NGRU(1,K,1) + L                                                    1052
   12 SVEC(L1,I) = 1.                                                       1053
   10 CONTINUE                                                              1054
      LL2 = LL - 2                                                          1055
      IF(LL2.LE.0) GOTO 14
       DO 16   K = 1,LL2
      I = I + 1                                                             1058
      L1 = NGRU(1,K+1,1) - 1                                                1059
      L2 = L1 + 1                                                           1060
      L3 = NGRU(2,K+1,1)                                                    1061
      L4 = L1 + L3                                                          1062
      DO 17   M = 1,L1                                                      1063
   17 SVEC(M,I) = -1./ FLOAT(L1)                                            1064
      DO 16   M = L2,L4                                                     1065
   16 SVEC(M,I) =  1./ FLOAT(L3)                                            1066
   14 IF(K1.LE.LL) GOTO 24
       LL1 = LL + 1
      DO 21   K = LL1,K1                                                    1069
      I = I + 1                                                             1070
      L1 = NGRU(1,K,1) - 1                                                  1071
      L2 = L1 + 1                                                           1072
      L3 = NGRU(2,K,1)                                                      1073
      L4 = L1 + L3                                                          1074
      DO 22   M = L5,L1                                                     1075
   22 SVEC(M,I) =-1./ FLOAT(L1-L5+1)                                        1076
      DO 23   M = L2,L4                                                     1077
   23 SVEC(M,I) = 1./ FLOAT(L3)                                             1078
   21 CONTINUE                                                              1079
   24  L4 = L5 - 1                                                          1080
      DO 26   M  = 1,L4                                                     1081
      SVEC(M,NZV) = -1./ FLOAT(L4)                                          1082
   26 SVEC(M,NZT) = 1./ FLOAT(NZT)                                          1083
      DO 27   M = L5,NZT                                                    1084
      SVEC(M,NZV) = 1./ FLOAT(NZT-L4)                                       1085
   27 SVEC(M,NZT) = 1./ FLOAT(NZT)                                          1086
      DO 60   K = 1,NZT                                                     1087
      A = .0                                                                1088
      DO 61   L = 1,NZT                                                     1089
   61 A = A + SVEC(L,K)**2                                                  1090
      A = 1./ SQRT(A)                                                       1091
      DO 62  L = 1,NZT                                                      1092
   62 SVEC(L,K) = A*SVEC(L,K)                                               1093
   60 CONTINUE                                                              1094
      DO 30   K = 1,NZT                                                     1095
      DO 30   L = 1,NZT                                                     1096
      RVEC(L,K) = SVEC(K,L)                                                 1097
   30 CONTINUE                                                              1098
      WRITE(NOUT,40)                                                        1099
       DO 34   K = 1,NZV                                                    1100
      L1 = NZV - K2 + K                                                     1101
      WRITE(NOUT,41)    K,(SVEC(M,K ),M,M=1,NZT)                            1102
      IF(K.GE.K1) GOTO 34
        DO 70    M = 1,NZT
   70 VEC(M,K,MFL) = SVEC(M,L1)
C      VEC ENTHAELT DIE CLUSTERRELATIVKOORDINATEN
C     DIE K-TE RELATIVKOORD. IST SUMME UEBER M (VEC(M,K)*R(M))
   34 CONTINUE                                                              1106
   40 FORMAT(//28H DEFINITION DER KOORDINATEN      //)                      1107
   41 FORMAT(3H S(,I2,3H) =,10(F6.3,3H R(,I2,1H)  ) / 7X,                   1108
     1         2(F6.3,3H R(,I2,1H)))                                        1109
   42 FORMAT(3H R(,I2,3H) =,10(F6.3,3H S(,I2,1H)  ) / 7X,                   1110
     1         2(F6.3,3H S(,I2,1H)))                                        1111
   43 FORMAT(//)                                                            1112
      WRITE(NOUT,43)                                                        1113
      DO 36   K = 1,NZT                                                     1114
   36 WRITE(NOUT,42)   K,(RVEC(M,K),M,M=1,NZV)                              1115
      WRITE(NOUT,43)
      DO 37 K= 1,K2
37    WRITE (NOUT,45) K,(VEC(M,K,MFL),M,M=1,NZT)
45    FORMAT(" VEC(",I2,")=",10(F6.3," R(",I2,")"))
      WRITE(NOUT,43)
      IIZ = 0                                                               1117
      DO 2130 K=1,K1                                                        1118
      NNZ   = NGRU(2,K,1)  - 1                                              1119
      IF(NNZ.EQ.0) GOTO 2130
       DO 2132    M = 1,NNZ
      IIZ = IIZ + 1                                                         1122
 2132 C(IIZ,K) = 1.
       WRITE(NOUT,44) (C(M,K),M=1,NZV)
   44 FORMAT(10F10.4)                                                       1125
 2130 CONTINUE                                                              1126
      WRITE(NOUT,43)
      DO 50   K = 1,K2                                                      1127
      KK=NZV-K2+K                                                           1128
      KI = K1 + K                                                           1129
   50 C(KK ,KI) = 1.
      DO 51 K=1,K4
51      WRITE (NOUT,44) (C(M,K),M=1,NZV)
      WRITE(NBAND2) ((RVEC(M,N),M=1,NZV),N=1,NZT)                           1131
      WRITE(NBAND2) ((SVEC(N,M),M=1,NZV),N=1,NZT)                           1132
      WRITE(NBAND2) ((C(N,K),N=1,NZV),K=1,K4)
C      DIE UEBERGABE VON C KANN AUF C(M,M,K) EINGESCHRAENKT WERDEN
      RETURN                                                                1134
      END                                                                   1135
      SUBROUTINE SPINWE(K,I)                                                1136
      PARAMETER (NZTMAX=12)
C
      COMMON /CSPIN/ NALG(NZTMAX,3)
      M=-2*(NALG(K,3)/2)+NALG(K,3)                                          1146
      I=2*M-1                                                               1147
      RETURN                                                                1148
      END                                                                   1149
      SUBROUTINE SPINUP(K,I)                                                1150
      PARAMETER (NZTMAX=12)
C
      COMMON /CSPIN/ NALG(NZTMAX,3)
C
C
      M=-2*(NALG(K,3)/2)+NALG(K,3)                                          1160
      IF (M)   1,1,2                                                        1161
    1 NALG(K,3) =NALG(K,3)-1                                                1162
      I=1                                                                   1163
      RETURN                                                                1164
    2 I=0                                                                   1165
      RETURN                                                                1166
      END                                                                   1167
      SUBROUTINE SPINDO(K,I)                                                1168
      PARAMETER (NZTMAX=12)
C
      COMMON /CSPIN/ NALG(NZTMAX,3)
C
      M=-2*(NALG(K,3)/2)+NALG(K,3)                                          1178
      IF(M)  2,2,1                                                          1179
    1 NALG(K,3)=NALG(K,3)+1                                                 1180
      I=1                                                                   1181
      RETURN                                                                1182
    2 I=0                                                                   1183
      RETURN                                                                1184
      END                                                                   1185
      SUBROUTINE UM(NF,MF)                                                  1186
C
      PARAMETER (NZTMAX=12, NZCMAX=4, NZCMAF=24)
C
      COMMON /CUM/ NPERM(NZCMAX+1,NZCMAF),NH(7,4),NKOR(NZTMAX,7),NP(5)
C
      KK=NH(2,MF)                                                           1196
      IF(KK.EQ.0) GOTO 3
      DO 1   K = 1,KK
      M=NH(K+2,MF)                                                          1199
      N=NPERM(K+1,NF)                                                       1200
      NN=NH(N+2,MF)                                                         1201
    1 NKOR(M,1)=NN                                                          1202
3      NP(MF+1)=NPERM(1,NF)*NP(MF)
C        NP(5) ENTHAELT GESAMTVORZEICHEN DER PERMUTATION P
C        NKOR(.,1) ENTHAELT PERMUTATION P
      RETURN                                                                1205
      END                                                                   1206

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