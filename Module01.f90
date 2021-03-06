      
      
MODULE MODULE01
      IMPLICIT NONE
      INTEGER, PARAMETER :: ELENUM=2		! NUMBER OF ELECTRONS
      REAL, PARAMETER :: ME = 9.1093897E-31, QE = 1.6022E-19
      INTEGER, PARAMETER:: SGRID=2000       ! ACTUAL GRID
      INTEGER, PARAMETER:: BINNUM=100      ! SAMPLING GRID
      REAL, PARAMETER :: DELTAT = 0.0001		!RK4 TIME STEP, (IN SECONDS)
      REAL, PARAMETER :: DELTASP = 1000.0		!DISTANCE SAMPLE STEP FOR CHECKING BINS
      REAL, PARAMETER :: TMAX = 10.0 		!MAXIMUM TIME TO FOLLOW ELECTRON (IN SECONDS)
      INTEGER, PARAMETER :: SWOK=1	
      REAL, PARAMETER :: NES=7.0, TES=10.0, NHALO=1.0, THALO=100.0
      INTEGER, PARAMETER :: EBEAM=400        ! MONOENERGETIC ENERGY BEAM
      REAL, PARAMETER :: EJAMD=0.0
      REAL, PARAMETER :: QZ=5.0E27			! OUTGASSING RATE OF COMET
      REAL, PARAMETER :: DHEL=1.2			! HELIOCENTRIC DISTANCE
      REAL, PARAMETER :: EFIELD = 0.00  
      !!! NEUTRAL DENSITY FILE
      INTEGER, PARAMETER :: NNEUT=10, NGRIDST=200   
      REAL NDEN(NGRIDST,NNEUT)
      REAL N2DENSITY(NGRIDST)
      REAL RADIUS2ST(NGRIDST)
      !!! B-FIELD FILE PARAMETERS
      !!! THESE ARE USED IN DERIVS AND INTERPSPLINE
      REAL RAD1(SGRID), SS1(SGRID), BGRID(SGRID), B2NDDER(SGRID)
      INTEGER, PARAMETER :: IPMAX=200            
      INTEGER, PARAMETER :: JBINS=5000          
      REAL PHELECPROD(IPMAX,JBINS)
      INTEGER, PARAMETER :: NNSPEC=3    
      REAL NEUTSPECDENSITY(NNSPEC, SGRID)
      
      !!!  MONTE CARLO (AND EPHOT) ENERGY GRID:
      REAL EEG(JBINS)                 ! MC ENERGY GRID
      REAL BINLOC(BINNUM)
      
      CONTAINS
      SUBROUTINE READNEUTDENS(RADIUS2ST,N2DENSITY)
      USE MODULE03, ONLY : FLNM
      IMPLICIT NONE
      REAL NDEN(NGRIDST,NNEUT)
      REAL RADIUS2ST(NGRIDST), N2DENSITY(NGRIDST)
      INTEGER I, J
      OPEN(UNIT=12,FILE=FLNM,STATUS='OLD')
      DO  I=1,4
      READ(12,*)
      ENDDO
      DO  I=1,NGRIDST
      READ(12,615) RADIUS2ST(I), (NDEN(I,J),J=1,10)
      RADIUS2ST(I)=(RADIUS2ST(I)+2575.0)*1E3  
      N2DENSITY(I)=NDEN(I,1)
      ENDDO
615   FORMAT(1X,F7.0,2X,1P10E11.3)
      CLOSE(12)
      RETURN
      END SUBROUTINE READNEUTDENS
      
      SUBROUTINE READINBFIELD(RAD1, SS1, BGRID, B2NDDER)
      USE MODULE03, ONLY : FLNM5
      IMPLICIT NONE
      REAL RAD1(SGRID), SS1(SGRID), BGRID(SGRID), B2NDDER(SGRID)
      INTEGER SGRID1, I
      OPEN(UNIT=21,FILE=FLNM5,STATUS='OLD') ! IN NANO-TESLA
      READ(21,*) SGRID1
      DO  I=1,SGRID
      READ(21,*) RAD1(I), SS1(I), BGRID(I), B2NDDER(I)
      ENDDO
      CLOSE(21)
      RETURN
      END SUBROUTINE READINBFIELD
      
      SUBROUTINE READINEPHOT(PHELECPROD)  
      USE MODULE03, ONLY : FLNM3
      IMPLICIT NONE
      REAL PHELECPROD(200,5000)
      INTEGER IPMAX, JBINS, I, J, J1
      OPEN(UNIT=15,FILE=FLNM3,STATUS='OLD')
      READ(15,*)
      READ(15,*)
      READ(15,620)IPMAX,JBINS
      DO J=1,JBINS
      J1=JBINS-J+1
      READ(15,621)(PHELECPROD(I,J1),I=1,IPMAX)
      ENDDO
620   FORMAT(2I6)
621   FORMAT(1X,1P10E11.3)
      
      CLOSE(15)
      RETURN
      END SUBROUTINE READINEPHOT
      
      REAL FUNCTION EV2CMS(E)
      IMPLICIT NONE
      REAL E
      EV2CMS=SQRT(2*E*QE/ME)*100
      RETURN
      END
      
      
      END MODULE MODULE01
      
      
