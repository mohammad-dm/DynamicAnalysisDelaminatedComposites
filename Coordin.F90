MODULE COORD
IMPLICIT NONE
PUBLIC COORDINATES
CONTAINS
SUBROUTINE COORDINATES (N,X,Cc,Lc,Ie,Je,B_type)
INTEGER,INTENT(out)::N,Cc
DOUBLE PRECISION,DIMENSION(:), allocatable,INTENT(OUT)::X,Lc
INTEGER,DIMENSION(:), allocatable,INTENT(OUT)::Ie,Je ,B_type ! I VE J NODES OF ELEMENTS
INTEGER:: S	,M
!OPEN(3,FILE="COOR.TXT", FORM="FORMATTED",STATUS="UNKNOWN")
OPEN(1,FILE="GEO.TXT", FORM="FORMATTED",STATUS="UNKNOWN")
!
!
    read(1,*)N			!Read NJ from BarStruc.inp
    read(1,*)Cc;		!Read NE from BarStruc.inp
    !
    write(13,*)'Number of Nodes    =',N    !
    write(13,*)'Number of Elements =',Cc	!
    write(13,*)'----------------------------'
    write(13,*)' '
!
        allocate (X(N),Lc(Cc))		    !Coor is of dimension NJ
        allocate (Ie(Cc), Je(Cc),B_type(Cc) )
!
DO S=1,N
! 'ENTER X,Y,Z COORDINATES OF THE NODE S'
READ(1,*)X(S)  
END DO
!150 CLOSE(3)
!
!OPEN(4,FILE="BAG.TXT", FORM="FORMATTED",STATUS="UNKNOWN")
!OPEN(4,FILE="ELEM_NODES.TXT", FORM="FORMATTED",STATUS="UNKNOWN")
DO M=1,Cc
! PRINT WHICH NODE IS ATTACHED TO THE ELEMENT M, FIRST I THEN J
READ(1,*)Ie(M),Je(M)   
END DO
!150 CLOSE(1)
!
!
DO M=1,Cc
	Lc(M)=SQRT((X(Je(M))-X(Ie(M)))**2)
END DO
!
OPEN(13,FILE="INPUT_CHECK.DAC",FORM="FORMATTED",STATUS="UNKNOWN")					 !************
!WRITE(13,*)'DUGUM NOKTASI KOORDINATLARI' 
WRITE(13,*)' '
WRITE(13,*)'COORDINATES OF NODES'
WRITE(13,*)'        NOD','      X(mm)    ','         Y(mm)    ' 
!
DO S=1,N
WRITE(13,9)S,X(S)
END DO
9 FORMAT(I10,F16.7)
!
WRITE(13,*)' '
!WRITE(13,*)'      CUBUK','       I UCU','          J UCU','       BOYU (m)'
WRITE(13,*)'      ELEMENT','       I END','          J END','       LENGTH (mm)'
!
DO M=1,Cc
WRITE(13,10)M,Ie(M),Je(M),Lc(M)
END DO
10 FORMAT(I10,2I15,F18.8)
WRITE(13,*)
!
!
DO M=1,Cc
READ(1,*)B_type(M)
!
write (13,*)'For beam number ',M 
    !
    If (B_type(M)==1) Then
        write (13,*)'Beam type was selected as Euler-Bernoulli'
    Else If (B_type(M)==2) Then
        write (13,*)'Beam type was selected as Timoshenko'
    Else If (B_type(M)==3) Then
        write (13,*)'Beam type was selected as Jun et al Composite'
    Else If (B_type(M)==4) Then
         write (13,*)'Beam type was selected as Gopalakrishnan et al Composite'
    Else If (B_type(M)==5) Then
         write (13,*)'Interface Spring'
    End If
    !
    Write(13,*)'-----------------------'
    Write(13,*)' '
    !
END DO
!
!
END SUBROUTINE
END MODULE COORD

