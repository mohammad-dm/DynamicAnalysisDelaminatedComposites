MODULE MAS_SLAVE
USE CARP,ONLY: DGEMMD,DGEMVD
IMPLICIT NONE
PUBLIC CONSTRAINT
CONTAINS
SUBROUTINE CONSTRAINT(NDOF,ME,SSA,TA,MSAY)
DOUBLE PRECISION,DIMENSION(:,:),INTENT(OUT)::TA
INTEGER,INTENT(IN)::NDOF,ME
INTEGER,INTENT(OUT)::SSA
INTEGER,DIMENSION(:),INTENT(IN)::MSAY
DOUBLE PRECISION::Ym,Ys
DOUBLE PRECISION ,DIMENSION(:,:),ALLOCATABLE::TAC !,TAB
DOUBLE PRECISION::dy
DOUBLE PRECISION ,DIMENSION(:,:),ALLOCATABLE::TR
INTEGER F,SL,SC,SD,SE,SF,I,J,K,m,HA,SER,SEK,a,b,ssb,ssb_pr
CHARACTER SOR
INTEGER,DIMENSION(:),ALLOCATABLE::BEWL
INTEGER,DIMENSION(:),ALLOCATABLE::NOL

!Note that everynode is a master first
!The design of the subroutine is such that for each MPC application a T matrix is generated in an orderly fashion
!A master node can be turned into a slave for the next MPC application however a slave cannot be a master afterwards


OPEN(13,FILE="INPUT_CHECK.DAC", FORM="FORMATTED",STATUS="UNKNOWN")
Print*,'MPC CONSTRAINT SUBROUTINE RUNNING'

ALLOCATE(BEWL(NDOF),TR(ME,ME),NOL(NDOF))
ALLOCATE(TAC(ME,ME))
SSA=ME
ssb=ME

DO I=1,SSA
DO J=1,SSA
TAC(I,J)=0.0
END DO
TAC(I,I)=1.0
END DO

!******************
OPEN(1125,FILE="CONBAG.TXT", FORM="FORMATTED",STATUS="UNKNOWN")


J=0
NoL=0
DO m=1,NDOF
IF(MSAY(m)==0) THEN
 J=J+1
NOL(m)=J
END IF
END DO


 WRITE(13,*)' '
 WRITE(13,*)'MPC CONSTRAINTS'
 WRITE(13,*)'   MASTER NOD ','   CONNECTED DOF ','   SLAVE NOD'
i=0

DO !This is a loop for each y

READ(1125,16)SOR


IF ((SOR=="Y").or. (SOR=="y")) THEN


!DO a=1,ME
!DO b=1,ME
!TR(a,b)=0
!END DO
!END DO
TR=0

!DO a=1,NDOF
!BEWL(a)=0
!END DO
BEWL=0

ssb_pr=ssb
!*****************************************************
!*****************************************************
!********


!********
!*****************************************************
!*****************************************************

i=i+1

READ(1125,17)F  !MASTER

!Enter the contraint type
! 1 X DIRECTON'
! 2 Y DIRECTON'
! 3 Z Rot '

! There has to be more types to impose rigid flange type conditions
! For instance 4 Z rot + nodes stay on the same line

ssb=ssb_pr-1	
	
SSA=SSB

 
READ(1125,17)SER						 !Constraint degree of freedom type

if( SER<=3) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SEK=SER
dy=0
ELSE if ( SER==4) then
SEK=1
READ(1125,*)Ym
READ(1125,*)Ys
dy=Ys-Ym
END IF

READ(1125,17)SL								  !SLAVE nodes
If (NOL(3*F-3+SEK)== 0) then
BEWL( NOL(3*SL-3+SEK) )=-1
else
BEWL( NOL(3*SL-3+SEK) )=NOL(3*F-3+SEK)
      !find new dof of slave  and relate it with new dof of master
end if

WRITE(13,*)F,'  ',SER,'   ',SL


!********
NOL(3*SL-3+SEK)=0
do m=3*SL-3+SEK+1,NDOF
if (NOL(m)/=0) then
NOL(m)=NOL(m)-1
end if
end do

!********
J=0
DO k=1,ssb_pr
!
IF( BEWL(k)>0 ) THEN
 TR(k,NOL(3*F-3+SEK))=1.0
 		 if ( SER==4) then
		 TR(k,NOL(3*F-3+3))=dy
		  end if
ELSE if ( BEWL(k)==-1 ) then
         if ( SER==4) then
		 TR(k,NOL(3*F-3+3))=dy
		 end if
ELSE 
  J=J+1
 TR(k,j)=1.0
END IF
!
END DO
!
!********
!
If (i==1) then
DO a=1,SSB_pr
DO b=1,SSB
TAC(a,b)=TR(a,b)
END DO
END DO

SC=SSB

else if (i>1) then

!DO a=1,ME
!DO b=1,SC
!TA(a,b)=TAC(a,b)
!END DO
!END DO
TA=TAC


SD=SSB
!TAC=0
CALL DGEMMD(ME, SD, SC, TA, TR ,TAC )
SC=SSB

end if


ELSE
!************************
EXIT
END IF
END DO


   DO a=1,ME
   DO b=1,SSA
   TA(a,b)=TAC(a,b)
   END DO
   END DO
!NEWL  related the new degree of freedom numbering with the previous one
!For example NEWL(4,2)=2 the fourth degree of freedom becomes second after an MPC condition is applied
!Note that first MPC is the support conditions
!SEK(n) knows what type of MPC is applied for the nth application

16 FORMAT(A1)
17 FORMAT(I5)
58 CLOSE(1125)

!******************************************************************************

Print*,'MPC CONSTRAINT SUBROUTINE ENDED'


END SUBROUTINE CONSTRAINT


	
 END MODULE
 

