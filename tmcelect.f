!!! HADI MADANIAN: THERE IS 'MAKEFILE' NOW IN THE FOLDER TO RUN ALL THE CODES
      
      PROGRAM TMCELECT
      USE MODULE01
      USE MODULE02
      USE MODULE03
      USE MODULE04
      IMPLICIT NONE
      INTEGER I,J,II,JJ, I2, N, IP, IRES
      REAL DENSITY2NDDER(NGRIDST)      
      REAL T, T2, TD, T3
      REAL SD, SSTEP, GYROPERIOD, DTSAMPLE
      REAL VTOT, VPERP, VPARA, U  
      REAL SVDA(2), SVDA2(2), SOLD, INIANG, INIANGDEG, BRATIO
      REAL SMIN, SMAX, DELTAS, B01, DSS, DS(2)
      REAL PROD(SGRID), PEPROD1(IPMAX) 
      REAL RKM(IPMAX), NSDENS(NNSPEC, SGRID)
      REAL HSUBS, HSTEP, TOTALPRRATE
      REAL PEEXT(SGRID,JBINS), PR_AREA(SGRID,JBINS)
      INTEGER ALTREF, ALTIND
      REAL PRATZ(JBINS), ZCUMPROB(JBINS), PRATE(SGRID), ECUMPROB(SGRID)
      REAL PEPERTOT(JBINS), INIELDATA(ELENUM,4), SIGANN(NNSPEC)
      REAL VTOT_EV, T4, P, E2, R1, R2
      REAL SIGAFF(NNSPEC+1), NEUTDENS
      INTEGER IND2, IND1, IND3, COID,COID2, COUNT2
      REAL FIONZ, DELTA, TEMP, CURE, ET, EET, TSIGNE, DELTAE
      REAL TE(SGRID), SNE(SGRID)
      INTEGER OFAG, CLTT
      REAL PRIMEVPERP(8), PRIMESVDA2(2,8), PRIMET(8), CCOUNTER(8)
      REAL V2ND_CMS
!----------- SOLAR WIND INITIAL BOUNDARY PARAMETERS -----------------      
      REAL PHIINF(JBINS)
      REAL PHIINF_CUMSUM(JBINS)
      REAL EEG2(JBINS)        
      
      REAL TESTTSIGNE(SGRID,100)
      
!-----------------OPENING OUTPUT FILES ------------------------------
      CALL READFILENAMES(FLNM,FLNM5,FLNM3,FLNM4,FLNM13,FLNM16,FLNM17)
      OPEN(UNIT=13,FILE=FLNM13,STATUS='UNKNOWN')
      OPEN(UNIT=14,FILE=FLNM4,STATUS='UNKNOWN')
      OPEN(UNIT=16,FILE=FLNM16,STATUS='UNKNOWN')
      
!--------------------------------------------------------------------
!------------ CALL MODULE SUBROUTINES TO READ IN DATA----------------
      CALL READINBFIELD(RAD1, SS1, BGRID, B2NDDER)
      !!!!! RAD1() IS THE RADIAL DISTANCE TO COMET IN CM
      !!!!! SS1(SGRID) IS THE DISTANCE ALONG FIELD LINE IN CM
      !!!!! BGRID IS THE MAGNITUDE OF THE MAGNETIC FIELD IN NT
      !!!!! B2NDDER IS THE SECOND DERIVATIVE OF THE FIELD WRT S
      
      CALL READXSECTION(ESIGA, BSPROBE, BSPROBI, ICURRE,ISIGA, IINDEX, IFSSIGA, SEC, TSIGA)
      !      CALL EIMPCOEFF(IONTHR,AK,AKB,AJ,AJB,TS,TA,TB,GAMS,GAMB,NIONS)
!-------------------------------------------------------------------------------------
      DO I=1,JBINS
      EEG(I)=0.05 + (I-1)*0.1
      ENDDO
      DELTAE=EEG(5)-EEG(4)      
      !-------------------------------------------------------------------------------------
C~    CALL READNEUTDENS(RADIUS2ST,N2DENSITY)    ! FOR TITAN CASE      
      
      FIONZ = 1.0E-6/DHEL/DHEL
      DO N = 1, SGRID
      TEMP = 100.
      TE(N) = TEMP
      DELTA = 1.0E-6*(300./TEMP)**0.5
      DO COUNT2 = 1, NNSPEC
      CALL NEUTCOMET(COUNT2,RAD1(N),NEUTDENS)
      NSDENS(COUNT2,N) = NEUTDENS
      IF(COUNT2.EQ.1) THEN
      SNE(N) = (FIONZ*NEUTDENS/DELTA)**0.5
      ENDIF
      ENDDO
      ENDDO
!C~       DO I=1,SGRID
!C~       WRITE(17,*) I, TE(I), SNE(I), NSDENS(1,I), RAD1(I)
!C~       ENDDO
!C~       CLOSE(17)
!C~       STOP
!-------------------------------------------------------------------------------------
!CREATE SAMPLING BIN LOCATIONS & WRITE BIN LOCATION WITH B-FIELD VALUES TO FILE
!-------------------------------------------------------------------------------------
      
      SMIN=MINVAL(SS1)
      SMAX=MAXVAL(SS1)
      SSTEP = (SS1(10)-SS1(9))*10  
      DELTAS=ABS(SMIN - SMAX)/(BINNUM - 1)
      WRITE(13,649) SWOK, ELENUM, SGRID, BINNUM
      WRITE(13,*), 'INITIAL BIN SPECS. BEFORE TRACKING:'
      WRITE(13,*), 'BIN_IND, BIN_DIST, BIN_RAD, B-FIELD, NEUT_DENS'
      DO I = 1, BINNUM
      BINLOC(I) = SMIN + (I-1)*DELTAS
      CALL BFIELDINTERP(BINLOC(I), B01, DSS)
      IP=MINLOC(ABS(SS1-BINLOC(I)),1)
      WRITE(13,650) I,BINLOC(I), RAD1(IP), B01, NSDENS(1, IP)
      ENDDO
      
649   FORMAT(1X,I3,1X,I8,1X,2(1X, I5),'  !!SOLARWIND, #ELECTRON, #SGRID, #BIN')
650   FORMAT(1X, I4, 4(1X, (ES12.4E2)))
      
      
      
!!!!!!!---------------------------------------------------------
!------------------------READ NEUTRAL DENSITY PROFILE
      
C~       CALL SPLINE(RADIUS2ST, N2DENSITY, NGRIDST, DENSITY2NDDER)
C~       DO I=1, SGRID
C~       PRINT*, RADIUS2ST(I), N2DENSITY(I)
C~       CALL INTERPSPLINE(RADIUS2ST, N2DENSITY, DENSITY2NDDER, NGRIDST, BINLOC(I), NDENEXTRAP(I), TEMP)
CC~       PRINT*, NDENEXTRAP(I)
C~       ENDDO
C~       PRINT*, NDENEXTRAP(1), NDENEXTRAP(SGRID)
C~       IRES=1      ! FOR DEFINITION OF THIS SEE THE SUBROUTINE
C~       IPMAX=200
C~       DO I=1,SGRID
C~             CALL EXTRAP(IPMAX,RADIUS2ST,NDEN(I,1), CORESRAD(I),NDENEXTRAP(I),1,IRES)
C~             IRES=0
C~       ENDDO
      
      
!-------------------------------------------------------------------------------------
!READ PHOTOELECTRON PRODUCTION DATA
!-------------------------------------------------------------------------------------
      CALL READINEPHOT(PHELECPROD)
      
C~         HSUBS=8.000E+05 
C~         HSTEP=2.000E+05
C~       RKM(1) = HSUBS
C~       DO I = 2, IPMAX       ! IPMAX IS DIFFERENT THAN THE BFIELD GRID POINTS
C~         RKM(I) = RKM(I-1) + HSTEP
C~       ENDDO
      
      !--------------------------------------------------------------------
      !-- EXTRAPOLATION OF EPHOT PRODUCTIONS TO RAD1 AND HIGHER ALTITUDES
      !--------------------------------------------------------------------      
C~       DO J=JBINS,1,-1  ! LOOP OVER ENERGIES

C~           ! EXTRAPOLATING PRODUCTIONS FROM EPHOT GRID TO 
C~           ! CURRENT DISTANCE ALONG FIELD GRID:
C~           IRES=1
C~           DO I=1,SGRID
C~             CALL EXTRAP(IPMAX,RKM,PHELECPROD(:,J),RAD1(I),PROD(I),1,IRES)
C~             IF(RAD1(I).GT.RKM(IPMAX)) THEN
C~                PROD(I) = PHELECPROD(IPMAX,J)*(RKM(IPMAX)/RAD1(I))**2
C~             ENDIF 
C~             PEEXT(I,J)=PROD(I)
C~             IRES=0
C~           ENDDO
C~           PR_AREA(:,J)=(PEEXT(:,J)/BGRID(:))* BGRID(1)
C~           PEPERTOT(J)=SUM(PEEXT(:,J)/BGRID(:))* BGRID(1)*(SS1(2)-SS1(1))
C~       ENDDO
      
!------------ WRITE OUT TOTAL PRODICTION RATE AT EACH ENERGY -----------
!------------ TO OUTPUT FILE FOR SCALING FACTOR LATER IN PROCESSING-----
      WRITE(16,*),'ENERGY   FLUX'
      DO J=1,JBINS
      WRITE(16,*) EEG(J), PEPERTOT(J)
      ENDDO
      CLOSE(16)

! Monoenergetic beam of electrons: 
      IF(SWOK .EQ. 1) THEN
      DO I=1,JBINS
      EEG2(I)=0
      ENDDO
      EEG2(EBEAM)=EEG(EBEAM)
      
      PHIINF=3.3456E7*NES*(EEG2-EJAMD)*EXP(-(EEG2-EJAMD)/TES)/TES**1.5+3.3456E7*NHALO*(EEG2-EJAMD)*EXP(-(EEG2-EJAMD)/THALO)/THALO**1.5
      PHIINF_CUMSUM=CUMSUM_R(PHIINF/SUM(PHIINF))
      DO I=1,ELENUM
      CALL INIT_RANDOM_SEED()
      CALL RANDOM_NUMBER(R1)
      IND1 = IFIRSTLOC(R1 < PHIINF_CUMSUM)
      VTOT=EV2CMS(EEG(IND1))
      CALL RANDIRECTION(SVDA, VPERP, VTOT)
      
      
      
      
      IF (SVDA(2).GT.0.0) THEN
      INIELDATA(I,1)=SS1(1)
      ELSE
      INIELDATA(I,1)=SS1(SIZE(SS1))
      ENDIF
      INIELDATA(I,2)=SVDA(2)         ! V_PARALLEL 
      INIELDATA(I,3)=VPERP              ! V_PERPENDICULAR 
      INIELDATA(I,4)=VTOT               ! TOTAL ENERGY
      ENDDO
      PRINT*, 'SOLAR WIND SELECTED'
      ELSE
      ALTREF=10
      PRATZ=PR_AREA(ALTREF,:)/SUM(PR_AREA(ALTREF,:))
      ZCUMPROB=CUMSUM_R(PRATZ)
      DO I=1,ELENUM
      CALL INIT_RANDOM_SEED()
      CALL RANDOM_NUMBER(R1)   
      CALL RANDOM_NUMBER(R2)   
      IND1 = IFIRSTLOC(R1 < ZCUMPROB)
      PRATE=PR_AREA(:,IND1)/SUM(PR_AREA(:,IND1))
      ECUMPROB=CUMSUM_R(PRATE)
      ALTIND = IFIRSTLOC(R2 < ECUMPROB)
      VTOT=EV2CMS(EEG(IND1))
      CALL RANDIRECTION(SVDA, VPERP, VTOT)
      INIELDATA(I,1)=SS1(ALTIND) 
      INIELDATA(I,2)=SVDA(2)         
      INIELDATA(I,3)=VPERP         
      INIELDATA(I,4)=VTOT          
      ENDDO
      ENDIF
      WRITE(14,*), 'TESTPAR# EFLAG GRIDIND E-IND COLL_COUNT COLLISIONTYPE POSITION  VPARALLEL   VPERP   MAG_MOM   TIME   TAU'
      
! LOOP OVER ELECTRONS STARTS HERE
      DO I = 1, ELENUM
      PRINT*, I
      T = 0.0
      T2 = 0.0
      T3 = 0.0
      SVDA(1)=INIELDATA(I,1)       
      SVDA(2)=INIELDATA(I,2)        
      VPERP=INIELDATA(I,3)           
      VTOT=INIELDATA(I,4)               
      IND3=1
      E2=0.0
      COID=0
      COID2=0
660   CONTINUE
      SOLD = SVDA(1)		
      CALL BFIELDINTERP(SVDA(1), B01, DSS)
      VTOT = SQRT(VPERP**2 + SVDA(2)**2)
      INIANG=ACOS(SVDA(2)/VTOT)
      INIANGDEG=INIANG*180.0/PI
      BRATIO=SQRT(MINVAL(BGRID)/MAXVAL(BGRID))
      !-------------------------------------------------------------------------------------
      !FIND MAGNETIC MOMENT 'U' ELECTRON KEEPS THIS VALUE UNTIL A COLLISION OCCURS
      !-------------------------------------------------------------------------------------
      U=(.5*ME*VPERP**2)/B01
      !-------------------------------------------------------------------------------------
      !TRAJECTORY BEGINS HERE WHICH STARTS AN ELECTRON OFF WITH THE ABOVE CONDITIONS ALONG
      !A GIVEN B-FIELD LINE
      !-------------------------------------------------------------------------------------
      CALL INIT_RANDOM_SEED()
      DO WHILE (T.LE.TMAX)
      OFAG=IND3
      CALL RK4(SVDA2, SVDA, T, DELTAT, U)
      CALL BFIELDINTERP(SVDA2(1), B01, DSS)
C~                   CALL VPER(VPERP, U, SS)
      VPERP = SQRT(ABS(2*U*B01/ME))
      VTOT = SQRT(VPERP**2 + SVDA2(2)**2)
      SIGAFF(:)=0.0
      VTOT_EV=.5*ME*(VTOT/100)**2/QE
      IND2 = MINLOC(ABS(SS1 - SVDA2(1)), 1)
      IND1 = MINLOC(ABS(EEG - VTOT_EV), 1)
      SIGANN=NSDENS(:,IND2)*TSIGA(:,IND1)
      SIGAFF(1:NNSPEC)=SIGANN
      CURE=EEG(IND1)
      ET = 8.618E-5*TE(IND2)
      EET= CURE - ET
      IF (EET .GT. 0.0) THEN
      TSIGNE =((3.37E-12*SNE(IND2)**0.97)/(CURE**0.94)) * ((EET)/(CURE-(0.53*ET)))**2.36
C~                   TSIGNE=TSIGNE*(RAVMU)/SINDIP/DAG
      ELSE
      TSIGNE=0.0
      ENDIF
      SIGAFF(NNSPEC+1)=TSIGNE/DELTAE
      T4=1/(VTOT*SUM(SIGAFF))
      TD = T - T2
      CLTT=0
      COID2=COID
      CALL RANDOM_NUMBER(R1)
      P=1-EXP(-TD/T4)
      IF (R1 < P) THEN
      COID=COID+1
      T2 = T     
      CALL COLLISION(SVDA2, VPERP, SIGAFF, IND2, IND1,IND3,E2,CLTT)
      U=(.5*ME*VPERP**2)/B01
      ENDIF
      SD = ABS(SVDA2(1) - SOLD)
      IF ((COID2.NE.COID) .OR. (SD.GT.SSTEP)) THEN
      CALL BINS(IND2, SVDA2, VPERP, I, U, T, IND1,OFAG, COID, T4, CLTT)
      SOLD = SVDA2(1)
      ENDIF
      IF (IND3.NE.OFAG) THEN
      PRIMEVPERP(OFAG)=VPERP
      PRIMESVDA2(:,OFAG)=SVDA2
      PRIMET(OFAG)=T
      CCOUNTER(OFAG)= COID
      INIANG=ACOS(SVDA2(2)/SQRT((SVDA2(2)**2+VPERP**2)))
      V2ND_CMS=EV2CMS(E2)
      VPERP=V2ND_CMS*SIN(INIANG)
      SVDA2(2)=V2ND_CMS*COS(INIANG)
      SVDA(1)=SVDA2(1)
      SVDA(2)=SVDA2(2)
      T=0.0
      T2=0.0
      COID=0
      GOTO 660
      ENDIF
      T = T + DELTAT
      SVDA(1) = SVDA2(1)
      SVDA(2) = SVDA2(2)
      IF ((SVDA(1).GT.SMAX).OR.(SVDA(1).LT.SMIN)) THEN
      PRINT*, 'PARTICLE EXITS THE BOX AT:'
      PRINT*, SVDA(1)/1E5,SVDA(2)/1E5
      T = TMAX + DELTAT
      ENDIF
      VTOT = SQRT(VPERP**2 + SVDA2(2)**2)
      IF (VTOT.LT.EV2CMS(EEG(1))) T=TMAX+DELTAT
      ENDDO
!-------------------------------------------------------------------------------------
!END OF TRAJECTORY FOR THIS ELECTRON
!-------------------------------------------------------------------------------------
      IND3=IND3-1
      IF (IND3 .GE. 1) THEN
      SVDA=PRIMESVDA2(:,IND3)
      VPERP=PRIMEVPERP(IND3)
      T=PRIMET(IND3)
      T2=T
      COID=CCOUNTER(IND3)
      GOTO 660
      ENDIF
      ENDDO      
      CLOSE(13)
      CLOSE(14)
      END PROGRAM
      

!-------------------------------------------------------------------------------------
!THIS IS THE RUNGE-KUTTA 4TH ORDER ODE SOLVE; THE ARGUMENTS AND DESCRIPTIONS ARE:
!Y: (INPUT) IS A TWO ELEMENT ARRAY:
!Y(1) = S: (LOCATION ON THE FIELD LINE)
!Y(2): (PERPENDICULAR LINE VELOCITY)
!DYDX: (INPUT) TWO ELEMENT ARRAY WITH THE DERIVATIVES OF Y
!N: (INPUT) NUMBER OF VARIABLES, IN THIS CASE 2
!SGRID1: # OF GRID POINTS
!X: (INPUT) TIME
!H: (INPUT) TIME STEP
!YOUT: (OUTPUT) TWO ELEMENT ARRAY OF S AND VPERP SOLVED
!DERIVS: A SUBROUTINE CALLED WITHIN THIS SUBROUTINE TO CALCULATE THE DERIVATIVES
!			OF THE Y(1) AND Y(2)
!U2: (INPUT) THE MAGNETIC MOMENT IS PASSED TO DERIVS
!XXO, YYO, ZZO, SS1, B2NDDER, EE1O, NN1O, XX1, YY1, ZZ1, BGRID, EE1, NN1:
!	(ALL INPUTS) ARE ALL 50 ELEMENT ARRAYS CONTAINING SPLINE COEFFICIENTS
!	WHICH ARE SENT TO DERIVS. THERE, THEY ARE USED BY INTERPSPLINE (AS DESCRIBED
!	IN THE MAIN PROGAM) TO FIND FIELDLINE VALUES
!-------------------------------------------------------------------------------------
      SUBROUTINE RK4(YOUT, Y, X, H, U2)
      USE MODULE01
      IMPLICIT NONE
      INTEGER N, I
      INTEGER, PARAMETER :: NMAX = 2
      REAL H, X, Y(NMAX), YOUT(NMAX), VPER, U2
      REAL H6, HH, XH, DYDX(NMAX), DYM(NMAX), DYT(NMAX), YT(NMAX)
      HH = H*0.5
      H6 = H/6
      XH = X + HH
      CALL DERIVS(DYDX, Y, U2)
      YT = Y + HH*DYDX
      CALL DERIVS(DYT, YT, U2)
      YT = Y + HH*DYT
      CALL DERIVS(DYM, YT, U2)
      YT = Y + H*DYM
      DYM = DYT + DYM
      CALL DERIVS(DYT, YT, U2)
      YOUT = Y + H6*(DYDX + DYT + 2.*DYM)
      RETURN
      END
!-------------------------------------------------------------------------------------
!THIS SUBROUTINE CALCULATES THE DERIVATIVES OF THE PARALLEL VELOCITY AND LOCATION:
!X: (INPUT) TIME
!SGRID1: NUMBER OF GRID POINTS
!Y1: (INPUT) TWO ELEMENT ARRAY
!Y1(1): LOCATION S
!Y1(2): PARALLEL VELOCITY
!DYDX: (OUTPUT)
!	DYDX(1): DERIVATIVE OF Y1(1)
!	DYDX(2): DERIVATIVE OF Y1(2)
!U3: (INPUT) MAGNETIC MOMENT
!XXO, YYO, ZZO, SS1, B2NDDER, EE1O, NN1O, XX1, YY1, ZZ1, BGRID, EE1, NN1:
!(INPUT) SPLINE COEFFICIENTS USED BY SUBROUTINE ZBFLD2
!-------------------------------------------------------------------------------------
C~       EXTERNAL UMAG
C~       EXTERNAL INTERPSPLINE
!       DERIVS(DS,SVAARRAY,U)
      SUBROUTINE DERIVS(DYDX, Y1, U3)
      USE MODULE01
      IMPLICIT NONE
      REAL Y1(2), DYDX(2), SL, BFIELD, U3
      REAL ST
      ST = Y1(1)
!       CALL BFLD2(Y1(1), BFIELD, SL, SS1, BGRID, B2NDDER)
      CALL BFIELDINTERP(ST, BFIELD, SL)
!       CALL INTERPSPLINE(SS1, EE1, EE1O, N2, ST, EDEN, DEDS)
!       CALL INTERPSPLINE(SS1, NN1, NN1O, N2, ST, NUMDEN, DNDS)
!       EFIELD = DEDS/NUMDEN*2.587E-4
      DYDX(1) = Y1(2)
      DYDX(2) = -(U3*SL)/ME - QE*EFIELD/ME
      RETURN
      END
      SUBROUTINE BFIELDINTERP(X,B01, B1STDER)
      USE MODULE01, ONLY : SS1, BGRID, B2NDDER, SGRID
      IMPLICIT NONE
      REAL X, B01, B1STDER
      CALL INTERPSPLINE(SS1, BGRID, B2NDDER, SGRID, X, B01, B1STDER)
      END SUBROUTINE BFIELDINTERP
      SUBROUTINE INTERPSPLINE(XA, YA, Y2A, N, X, Y, DYDX)
      
      !INPUT: XA, YA, Y2A, N, X
      !OUTPUT: Y, DYDS
      !-------------------------------------------------------------------------------------
      !THIS SUBROUTINE USES THE COEFFICIENTS CALCULATED IN SUBROUTINE SPLINE. CODE USED 
      !IS FROM PRESS, DESCRIPTION FROM PRESS IS LISTED BELOW
      !-------------------------------------------------------------------------------------
      !-------------------------------------------------------------------------------------
      !GIVEN THE ARRAYS XA(1:N) AND YA(1:N) OF LENGTH N, WHICH TABULATE A FUNCTION (WITH
      !XAI'S IN ORDER), AND GIVEN THE ARRAY Y2A(1:N), WHICH IS THE OUTPUT FROM SPLINE
      !ABOVE, AND GIVEN A VALUE OF X, THIS ROUTINE RETURNS A CUBIC-SPLINE INTERPOLATED
      !VALUE Y
      !-------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N
      REAL X,Y, DYDX, XA(N), Y2A(N), YA(N)
      INTEGER K, KHI, KLO
      REAL A, B, H
      KLO = 1
      KHI = N
1     IF(KHI-KLO.GT.1) THEN
      K = (KHI+KLO)/2
      IF(XA(K).GT.X) THEN
      KHI = K
      ELSE
      KLO = K
      ENDIF
      GOTO 1
      ENDIF
      
      H = XA(KHI) - XA(KLO)
      IF(H.EQ.0.) THEN
      PRINT*, 'BAD XA INPUT IN INTERPSPLINE'
      STOP
      ENDIF
      A = (XA(KHI) - X)/H
      B = (X - XA(KLO))/H
      Y = A*YA(KLO) + B*YA(KHI) + ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6
      DYDX = (YA(KHI) - YA(KLO))/H - ((3*A**2-1)*H*Y2A(KLO))/6 + (3*B**2-1)*H*Y2A(KHI)/6
      RETURN
      END SUBROUTINE
      SUBROUTINE SPLINE(X, Y, N, Y2)
      IMPLICIT NONE
      INTEGER N, NMAX
      REAL YP1, YPN, X(N), Y(N), Y2(N)
      PARAMETER (NMAX = 500)
!-------------------------------------------------------------------------------------
!GIVEN ARRAYS X(1:N) AND Y(1:N) CONTAINING A TABULARED FUNCTION I.E. YI = F(XI),
!WITH X1 < X2 < ::: < XN, AND GIVEN VALUES YP1, AND YPN FOR THE RST DERIVATIVE OF THE
!INTERPOLATING FUNCTION AT POINTS 1, AND N, RESPECTIVELY, THIS ROUTINE RETURNS AN ARRAY
!CY2(1:N) OF LENGTH N WHICH CONTAINS THE SECOND DERIVATIVES OF THE INTERPOLATING FUNCTION
!AT THE TABULATED POINTS XI. IF YP1 AND/OR YPN ARE EQUAL TO 1E10 30 OR LARGER, THE ROUTINE
!IS SIGNALED TO SET THE CORRESPONDING BOUNDARY CONDITION FOR A NATURAL APLINE, WITH
!ZERO SECOND DERIVATIVE ON THAT BOUNDARY
!-------------------------------------------------------------------------------------
      INTEGER I, K
      REAL P, QN, SIG, UN, U(NMAX)
      YP1=0.0
      YPN=0.0
      IF (YP1.GT.99E30) THEN
      Y2(1) = 0
      U(1) = 0
      ELSE 	
      Y2(1) = -0.5
      U(1) = (.3/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO I = 2, N-1
      SIG = (X(I)-X(I-1))/(X(I+1)-X(I-1))
      P = SIG*Y2(I-1)+2
      Y2(I) = (SIG-1.)/P
      U(I) = (6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))/(X(I) - X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      ENDDO
      IF (YPN.GT.99E30) THEN 		
      QN = 0
      UN = 0
      ELSE						
      QN = 0.5
      UN = (3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N) = (UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO K = N-1, 1, -1 			
      Y2(K) = Y2(K)*Y2(K+1)+U(K)
      ENDDO      
      RETURN
      END SUBROUTINE
      
      
!       SUBROUTINE VTOTMAG(VTOT)
!-------------------------------------------------------------------------------------
!THIS SUBROUTINE CREATES A VELOCITY FROM A MAXWELLIAN IN MONTE CARLO FASHION
!VTOT: (OUTPUT) TOTAL VELOCITY
!RAN2: EXTERNAL FUNCTION
!SEED1: (INPUT) A NEGATIVE INTEGER USED IN FUNCTION RAN2
!-------------------------------------------------------------------------------------
C~       USE MODULE01, ONLY : EINIT, QE, ME
C~       IMPLICIT NONE
C~       REAL VTOT
C~       VTOT = SQRT(2*EINIT*QE/ME)*10
CC~      CALL INIT_RANDOM_SEED()
CC~      CALL RANDOM_NUMBER(R1)
CC~       VTOT = VTOT*SQRT(-LOG(1-R1))
C~       RETURN
C~       END SUBROUTINE
      
      SUBROUTINE RANDIRECTION(SVDA, VPERP, VTOT)
      IMPLICIT NONE
      REAL SVDA(2), VPERP, VTOT
      REAL R, R1, R2
      CALL INIT_RANDOM_SEED()
      CALL RANDOM_NUMBER(R)
      R1=2*R-1	  
      SVDA(2) = VTOT*R1
      VPERP = VTOT*(SQRT(1.0-R1**2))      
      RETURN
      END
      
      SUBROUTINE BINS(LOCIND, SVDA2, VPERP, PARNUM, U, T, IND1, IND3, COUNT1, T4, CLTT)
      USE MODULE01, ONLY : BINNUM, BINLOC
      IMPLICIT NONE
      INTEGER PARNUM, III, IND1, IND3, COUNT1, CLTT,LOCIND
      REAL SOLD, SVDA2(2), VPERP, U, T, T4
      III=LOCIND
      DO III = 1, BINNUM
      IF (((SOLD.LE.BINLOC(III)).AND.(SVDA2(1).GE.BINLOC(III))).OR.((SOLD.GE.BINLOC(III)).AND.(SVDA2(1).LE.BINLOC(III)))) THEN
      WRITE(14,701) PARNUM, IND3, III, IND1, COUNT1,CLTT, SVDA2(1), SVDA2(2), VPERP, U, T, T4
      ENDIF
      ENDDO
701   FORMAT((1X, I6),3(2X, I5), 2(3X, I5), 6(2X, (ES12.4E2)))
      RETURN
      END SUBROUTINE
      
      SUBROUTINE INIT_RANDOM_SEED()
      USE ISO_FORTRAN_ENV, ONLY: INT64
      IMPLICIT NONE
      INTEGER, ALLOCATABLE :: SEED(:)
      INTEGER :: I, N, UN, ISTAT, DT(8), PID
      INTEGER(INT64) :: T
      CALL RANDOM_SEED(SIZE = N)
      ALLOCATE(SEED(N))
      OPEN(NEWUNIT=UN, FILE="/DEV/URANDOM", ACCESS="STREAM", FORM="UNFORMATTED", ACTION="READ", STATUS="OLD", IOSTAT=ISTAT)
      IF (ISTAT == 0) THEN
      READ(UN) SEED
      CLOSE(UN)
      ELSE
      CALL SYSTEM_CLOCK(T)
      IF (T == 0) THEN
      CALL DATE_AND_TIME(VALUES=DT)
      T = (DT(1) - 1970) * 365_INT64 * 24 * 60 * 60 * 1000 
     &                  + DT(2) * 31_INT64 * 24 * 60 * 60 * 1000 
     &                  + DT(3) * 24_INT64 * 60 * 60 * 1000 
     &                  + DT(5) * 60 * 60 * 1000 
     &                  + DT(6) * 60 * 1000 + DT(7) * 1000 
     &                  + DT(8)
      END IF
      PID = GETPID()
      T = IEOR(T, INT(PID, KIND(T)))
      DO I = 1, N
      SEED(I) = LCG(T)
      END DO
      END IF
      CALL RANDOM_SEED(PUT=SEED)
      CONTAINS
      FUNCTION LCG(S)
      INTEGER :: LCG
      INTEGER(INT64) :: S
      IF (S == 0) THEN
      S = 104729
      ELSE
      S = MOD(S, 4294967296_INT64)
      END IF
      S = MOD(S * 279470273_INT64, 4294967291_INT64)
      LCG = INT(MOD(S, INT(HUGE(0), INT64)), KIND(0))
      END FUNCTION LCG
      END SUBROUTINE INIT_RANDOM_SEED
      
C~       SUBROUTINE ENERGYLOSE(S5, VPE, VT, EL)
C~       REAL S5(2), VPE, VT, EL, ME, E, VF
C~       ME = 9.1093897E-31
C~       VF = SQRT(VT**2 - (2*EL*1.602177E-19)/ME)
C~       VPE = VPE*(VF/VT)
C~       S5(2) = S5(2)*(VF/VT)
C~       RETURN
C~       END
      
      
      SUBROUTINE NEUTCOMET(FLAG,RR,DENS)
      USE MODULE05
      USE MODULE01, ONLY : QZ, DHEL
      IMPLICIT NONE
      REAL UN, TAU1, TAU2, TAU3, LAMBDA1, LAMBDA2
      REAL LAMBDA3, DENSTOT, DENS
      REAL RR
      INTEGER FLAG
      UN = 1.0E+5
      TAU1 = 1.0E+6*DHEL*DHEL
      TAU2 = TAU1
      TAU3 = TAU1
      LAMBDA1 = TAU1*UN
      LAMBDA2 = TAU2*UN
      LAMBDA3 = TAU3*UN
      DENSTOT = QZ/(RR*RR*UN*4*PI)
      IF(FLAG.EQ.1) THEN
      DENS = 0.85*DENSTOT*EXP(-RR/LAMBDA1)
      ELSEIF(FLAG.EQ.2) THEN
      DENS = 0.08*DENSTOT*EXP(-RR/LAMBDA2)
      ELSE
      DENS = 0.07*DENSTOT*EXP(-RR/LAMBDA3)
      ENDIF
      
      RETURN
      END SUBROUTINE NEUTCOMET
      
      SUBROUTINE EXTRAP(N,X,Y,XP,YP,IEXP,IRESET)
      INTEGER N, IEXP, IRESET, IFLAG
      REAL X(N),Y(N)
      IF (IRESET.EQ.1) J=1
      CALL FINDJ(N,X,XP,J,IFLAG)
      IF (IEXP.EQ.0) THEN
      CALL LINEAR(N,X,Y,XP,J,IFLAG,YP)
      ELSE IF (IEXP.EQ.1) THEN
      CALL EXPO(N,X,Y,XP,J,IFLAG,YP)
      ENDIF
      RETURN
      END
      
      SUBROUTINE LINEAR(N,X,Y,XP,J,IFLAG,YP)
      REAL X(N),Y(N)
      IF (IFLAG .EQ. 1) THEN
      YP=Y(J)
      ELSE
      IF (IFLAG .EQ. -1) THEN
      JL=1
      JR=2
      ELSE IF (IFLAG .EQ. -2) THEN
      JL=N-1
      JR=N
      ELSE IF (IFLAG .EQ. 0) THEN
      JL=J
      JR=J+1
      ELSE 
      PRINT *,'IFLAG NOT RETURNED RIGHT'
      ENDIF
      YP=(Y(JR)-Y(JL))*(XP-X(JL))/(X(JR)-X(JL))+Y(JL)
      ENDIF
      
      RETURN
      END
      
      SUBROUTINE EXPO(N,X,Y,XP,J,IFLAG,YP)
      REAL X(N),Y(N),YR,YL,XR,XL,YP,XPL
      INTEGER JL,JR
      IF (IFLAG .EQ. 1) THEN
      YP=Y(J)
      ELSE 
      IF (IFLAG .EQ. -1) THEN
      JL=1
      JR=2
      ELSE IF (IFLAG .EQ. -2) THEN
      JL=N-1
      JR=N
      ELSE IF (IFLAG .EQ. 0) THEN
      JL=J
      JR=J+1
      ELSE 
      PRINT *,'IFLAG NOT RETURNED RIGHT'
      ENDIF
      
      IF ((Y(JR).GT.1.E-30) .AND. (Y(JL).GT.1.E-30) 
     &	.AND. (X(JR).GT.1.E-30) .AND. (X(JL).GT.1.E-30)) THEN
      YR=ALOG(Y(JR))
      YL=ALOG(Y(JL))
      XR=ALOG(X(JR))
      XL=ALOG(X(JL))
      XPL=ALOG(XP)
      YP=(YR-YL)*(XPL-XL)/(XR-XL)+YL
      YP=EXP(YP)
      ELSE
      YP=0.0
      ENDIF
      ENDIF
      
      RETURN
      END
      
      SUBROUTINE FINDJ(NPT,XX,X,J,IFLAG)
      INTEGER IFLAG, NPT
      INTEGER, PARAMETER :: NMX=500
      REAL XX(NPT)
      IF (X .LT. XX(1)) GOTO 1
      IF (X .GT. XX(NPT)) GOTO 2
      J=1
13    IF (X .EQ. XX(J)) GOTO 14
      IF ((X .GT. XX(J)) .AND. (X .LT. XX(J+1))) GOTO 15
      J=J+1
      IF (J .LE. NPT) GOTO 13
      PRINT *,' ## INFINITE LOOP IN FINDJ, X=',X
      STOP 
1     IFLAG=-1
      RETURN
2     IFLAG=-2
      RETURN
14    IFLAG=1
      RETURN
15    IFLAG=0
      RETURN
      END
      
      SUBROUTINE PARA2(X0,S,R)
      REAL X0, S, R, RP
      REAL X1, RFLANK
      REAL THETA, ALPHA
      X1=X0*0.5
      CALL PARASR(X1,S,RP)
      THETA=2.0*ACOS(SQRT(X1/RP))
      R=RP**2+X1**2+RP*X0*COS(THETA)
      R=SQRT(R)
      ALPHA=ASIN(RP/R*SIN(THETA))
      ALPHA=ALPHA*180.0/3.1415926
      RFLANK=SQRT(2.0)*X0
      IF (R.GT.RFLANK) ALPHA=180.0-ALPHA
      
      RETURN
      END
      
      SUBROUTINE PARASR(X0,S,R)
      REAL X0, S, R
      REAL S1, R1, S2, R2
      IF (S .LE. (X0*0.3)) GOTO 404
      R=X0+S
      S1=((R-X0)*R)**0.5+X0*ALOG(((R-X0)**0.5+R**0.5)/X0**0.5)
      R1=(S-X0*ALOG(((R-X0)**0.5+R**0.5)/X0**0.5))**2/R+X0
      DO I=1,50
      S2=((R1-X0)*R1)**0.5+X0*ALOG(((R1-X0)**0.5+R1**0.5)/X0**0.5)
      IF (ABS((S2-S)/S) .LT. 1.E-4) GOTO 402
      IF (ABS((S2-S1)/S) .LE. 1.E-4) GOTO 402
      R2=(R1-R)*(S-S1)/(S2-S1)+R
      IF (R2 .LE. X0) R2=X0
      R=R1
      S1=S2
      R1=R2
      ENDDO
402   CONTINUE
      R=R1
      GOTO 405
404   R=0.2492*S**2/X0+X0
405   CONTINUE
      RETURN
      END
      
      
      
      
      
      
      
