Module Spectral_Stiffness
!
USE PolynomialRoots
USE COMPLEX_INVERSE 
!
IMPLICIT NONE
!
Public:: PROPERTIES,COMPOSITE_PROPERTIES,INTERFACE_PROPERTIES,EULER_BERNOULLI_ELEMENT_SPECTRAL_STIFFNESS,INTERFACE_SPRING_SPECTRAL_STIFFNESS,INTERFACE_SPRING_CONVENTIONAL_STIFFNESS,GOPALA_COMPOSITE_ELEMENT_SPECTRAL_STIFFNESS,ASSEMB,solver,JUN_COMPOSITE_ELEMENT_SPECTRAL_STIFFNESS,TIMOSHENKO_ELEMENT_SPECTRAL_STIFFNESS
!
!
contains 
!
!-----------------------------------------------------------------
!
    SUBROUTINE PROPERTIES (IELEM,A_,B_,D_,I_)

    !DEFINING VARIABLES-------------------------------------------
    Double precision, dimension(:,:,:), intent (inout) :: A_,B_,D_	
    Double precision, dimension(:,:), intent (inout) ::I_       
    Double precision::E,G,A,II,Den
    !
    integer, intent (in) ::IELEM 
    !
    !-----------------------------------------------------------------
    OPEN(2,FILE="Isotropic_Beam.txt", FORM="FORMATTED",STATUS="UNKNOWN")     
    !Program reads from file BarStruc.inp
    OPEN(13,FILE="INPUT_CHECK.DAC",FORM="FORMATTED",STATUS="UNKNOWN")
    !
    write(13,*)'----------------------------'
    write(13,*)' '
    write(13,*)' '
    write(13,*)'****************************'
    write(13,*)' '
    !
    !


    !***********************************************************************************/
				    !****	Material and Cross-sectional properties	****/
!   
    !do  i=1,NE  !for each member
    !    
    read(2,*)E				    !read the modulus of elasticity
    read(2,*)G				    !read the shear modulus 
    read(2,*)A				    !read the cross-sectinal area
    read(2,*)Den			    !read density
    read(2,*)II					!read the moment of inertia
    write(13,*)'  Elasticity Modulus of member   ',IELEM,'  = ',E	
    !write the modulus of elasticity for each member
    write(13,*)'  Shear Modulus of member   ',IELEM,'  = ',G	
    !write the modulus of elasticity for each member
    write(13,*)'  Cross-sectional Area of member ',IELEM,'  = ',A	
    !write the he cross-sectinal area for each member
    write(13,*)'  Density of member ',IELEM,'  = ',Den	
    !write the he cross-sectinal area for each member
    write(13,*)'  Moment of inertia ',IELEM,'  = ',II	
    !write the moment of inertia for each member
    write(13,*)' '
    !
    A_(1,1,IELEM)=E*A
    !A_(2,2,IELEM)=G*A
    !changed
    A_(3,3,IELEM)=G*A
    D_(1,1,IELEM)=E*II
    I_(1,IELEM)=Den*A
    I_(2,IELEM)=Den*A
        !
    !end do
    !
    write(13,*)'----------------------------'
    write(13,*)' '
    !
    !
    END SUBROUTINE PROPERTIES
!
!--------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------
!
    SUBROUTINE COMPOSITE_PROPERTIES (IELEM,A_,B_,D_,I_)  

    !DEFINING VARIABLES-------------------------------------------
    Double precision, dimension(:,:,:), intent (inout) :: A_,B_,D_	
    Double precision, dimension(:,:), intent (inout) ::I_ 
    integer, intent (in) ::IELEM

    
!    
!        
    double precision,dimension(:),allocatable::E1,E2,nu12,nu21,G12,G13,G23,fi,Z,Den,TA
    double precision::Q11,Q12,Q22,Q66
    double precision::Q_11,Q_12,Q_16,Q_22,Q_26,Q_66,t
    integer M,k,N_lay,udo,i,j
    !integer fj,NDOF
    integer::Nlay
    !Added by Mohammad=========================
    ! This block shares I_ between ELEMENT_SPECTRAL_STIFFNESS and COMPOSITE_ELEMENT_SPECTRAL_STIFFNESS Subroutines
    !The shared variable between two subroutines can not be allocatable, so a large value, 10000, was selected for I_ dimension
    double precision::Q44,Q55,Q_55
    double precision,dimension(:),allocatable::WIDTH
    !Added by Mohammad=========================
!************************************************************************
!    
    !-----------------------------------------------------------------
    OPEN(2,FILE="Composite_Beam.txt", FORM="FORMATTED",STATUS="UNKNOWN")     
    !Program reads from file BarStruc.inp
    OPEN(13,FILE="INPUT_CHECK.DAC",FORM="FORMATTED",STATUS="UNKNOWN")

				    !****	Material and Cross-sectional properties	****/
    N_lay=20
    allocate(E1(N_lay),E2(N_lay),nu12(N_lay),nu21(N_lay),G12(N_lay),G13(N_lay),G23(N_lay),fi(N_lay)) 
    allocate(Z(N_lay+1),Den(N_lay ) )
    !ADDED BY mOHAMMAD
    ALLOCATE ( WIDTH (N_lay) )
    !-----------------
!
!    
!        
read(2,*)Nlay
!
Do k=1,Nlay
!
READ(2,*)E1(k)
READ(2,*)E2(k)
READ(2,*)nu12(k)
READ(2,*)nu21(k)
READ(2,*)G12(k)
READ(2,*)G13(k)
READ(2,*)G23(k)
READ(2,*)fi(k)
READ(2,*)Z(k),Z(k + 1)
!
read(2,*)Den(k)
read(2,*)WIDTH(k)
!
!
WRITE(13,*)'Material Properties of Layers'
WRITE(13,*)'Lay','   E1    ','       E2    ','   nu12  ','  nu21     ',' G12  ','    G13    ','     G23 ','    fi' 	 													
WRITE(13,*)'   ','  N/mm2  ','  N/mm2   ','        ','            ',' N/mm2 ',' N/mm2   ',' N/mm2  ','  rad '
WRITE(13,8)k,' ',E1(k) ,'  ',E2(k),' ',nu12(k),' ',nu21(k),' ',G12(k),' ',G13(k),' ',G23(k),' ',fi(k)
WRITE(13,*)' ' 
!
t=Z(k + 1)-Z(k)
!
write(13,*)'Layer tickness',t
write(13,*)'Layer Width',WIDTH(k)
WRITE(13,*)' ' 
WRITE(13,*)'--------------------------------------------- ' 
!
!
End do     
!
WRITE(13,*)' ' 
WRITE(13,*)'--------------------------------------------- ' 
!

!--------------------------------------------    

!
!
!-----------------------------------------------------
!
i=IELEM
!
Do k=1,Nlay
!
! 
 Q11=E1(k)/(1-nu12(k)*nu21(k))
! 
 Q12=nu12(k)*E2(k)/(1-nu12(k)*nu21(k))
!
 Q22=E2(k)/(1-nu12(k)*nu21(k))
!
 Q66 = G12(k)
 Q44 = G23(k)
 Q55 = G13(k)
!
  !*******************************
!
  Q_11= Q11*(cos(fi(k)))**4+2*(Q12+2*Q66)*(sin(fi(k)))**2*(cos(fi(k)))**2+Q22*(sin(fi(k)))**4
! 
  Q_12=(Q11+Q22-4*Q66)*(sin(fi(k)))**2*(cos(fi(k)))**2 + Q12*(sin(fi(k))**4+cos(fi(k))**4)
!
  Q_22= Q11*(sin(fi(k)))**4+2*(Q12+2*Q66)*(sin(fi(k)))**2*(cos(fi(k)))**2+Q22*(cos(fi(k)))**4
!
  Q_16=(Q11-Q12-2*Q66)*(sin(fi(k)))*(cos(fi(k)))**3 + (Q12-Q22+2*Q66)*(sin(fi(k))**3*cos(fi(k)))
!
  Q_26=(Q11-Q12-2*Q66)*(sin(fi(k)))**3*(cos(fi(k))) + (Q12-Q22+2*Q66)*(sin(fi(k))*cos(fi(k))**3)
!
  Q_66=(Q11+Q22-2*Q12-2*Q66)*(sin(fi(k)))**2*(cos(fi(k)))**2 + Q66*(sin(fi(k))**4+cos(fi(k))**4)
!
    !Added by Mohammad
  Q_11 = Q11*(cos(fi(k)))**4+2*(Q12+2*Q66)*(sin(fi(k)))**2*(cos(fi(k)))**2+Q22*(sin(fi(k)))**4
  Q_55 = Q55*(cos(fi(k)))**2 + Q44*(sin(fi(k)))**2
!==========================
!   
   A_(1,1,i) = A_(1,1,i) + Q_11*( Z(k + 1)-Z(k) ) * width(k)
   A_(1,2,i) = A_(1,2,i) + Q_12*( Z(k + 1)-Z(k) ) * width(k)
   A_(2,2,i) = A_(2,2,i) + Q_22*( Z(k + 1)-Z(k) ) * width(k)
   A_(1,3,i) = A_(1,3,i) + Q_16*( Z(k + 1)-Z(k) ) * width(k)
   A_(2,3,i) = A_(2,3,i) + Q_26*( Z(k + 1)-Z(k) ) * width(k)
   !A_(3,3,i) = A_(3,3,i) + Q_66*( Z(k + 1)-Z(k) ) * width(k)
   A_(2,1,i) = A_(1,2,i)        ! Symmetry of the matrix A
   A_(3,1,i) = A_(1,3,i)        ! Symmetry of the matrix A
   A_(3,2,i) = A_(2,3,i)        ! Symmetry of the matrix A
!   
      !Edited by Mohammad: Changed Q_66 that is not needed to Q_55---
   A_(3,3,i) = A_(3,3,i)+Q_55*( Z(k + 1)-Z(k) )
!
   B_(1,1,i)= B_(1,1,i) + 1.0/2.0 * Q_11 * ( Z(k + 1)**2-Z(k)**2 ) * width(k)
   B_(1,2,i)= B_(1,2,i) + 1.0/2.0 * Q_12 * ( Z(k + 1)**2-Z(k)**2 ) * width(k)
   B_(2,2,i)= B_(2,2,i) + 1.0/2.0 * Q_22 * ( Z(k + 1)**2-Z(k)**2 ) * width(k)
   B_(1,3,i)= B_(1,3,i) + 1.0/2.0 * Q_16 * ( Z(k + 1)**2-Z(k)**2 ) * width(k)
   B_(2,3,i)= B_(2,3,i) + 1.0/2.0 * Q_26 * ( Z(k + 1)**2-Z(k)**2 ) * width(k)
   !B_(3,3,i)= B_(3,3,i) + 1.0/2.0 * Q_66 * ( Z(k + 1)**2-Z(k)**2 ) * width(k)
   B_(2,1,i)= B_(1,2,i)         ! Symmetry of the matrix B
   B_(3,1,i)= B_(1,3,i)         ! Symmetry of the matrix B
   B_(3,2,i)= B_(2,3,i)         ! Symmetry of the matrix B
!   
      !Edited by Mohammad: Changed Q_66 that is not needed to Q_55---
   B_(3,3,i)= B_(3,3,i)+1.0/2.0*Q_55*( Z(k + 1)**2-Z(k)**2 ) * width(k)
!
   D_(1,1,i)= D_(1,1,i)+1.0/3.0*Q_11*( Z(k + 1)**3-Z(k)**3 ) * width(k)
   D_(1,2,i)= D_(1,2,i)+1.0/3.0*Q_12*( Z(k + 1)**3-Z(k)**3 ) * width(k)
   D_(2,2,i)= D_(2,2,i)+1.0/3.0*Q_22*( Z(k + 1)**3-Z(k)**3 ) * width(k)
   D_(1,3,i)= D_(1,3,i)+1.0/3.0*Q_16*( Z(k + 1)**3-Z(k)**3 ) * width(k)
   D_(2,3,i)= D_(2,3,i)+1.0/3.0*Q_26*( Z(k + 1)**3-Z(k)**3 ) * width(k)
   !D_(3,3,i)= D_(3,3,i)+1.0/3.0*Q_66*( Z(k + 1)**3-Z(k)**3 ) * width(k)
   D_(2,1,i)= D_(1,2,i)     ! Symmetry of the matrix D
   D_(3,1,i)= D_(1,3,i)     ! Symmetry of the matrix D
   D_(3,2,i)= D_(2,3,i)     ! Symmetry of the matrix D
!   
      !Edited by Mohammad: Changed Q_66 that is not needed to Q_55---
   D_(3,3,i)= D_(3,3,i)+1.0/3.0*Q_55*( Z(k + 1)**3-Z(k)**3 ) * width(k)
!
   I_ (1,i) = I_(1,i) + Den(k) * WIDTH(K) * ( Z(k + 1)-Z(k) )
   I_ (2,i) = I_(2,i) + Den(k) * WIDTH(K) / 2 * ( Z(k + 1)**2-Z(k)**2 )
   I_ (3,i) = I_(3,i) + Den(k) * WIDTH(K) / 3 * ( Z(k + 1)**3-Z(k)**3 )
!
!
!
!t=Z(k + 1,i)-Z(k,i)
!!
!!
!densA(i)=densA(i)+Den(k)*t
!
End do
!
WRITE(13,*)'||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||' 
!
If (B_(1,1,i)<A_(1,1,i) *1.0E-08) then
B_(1,1,i)=0
End If
!
If (abs(I_ (2,i))<abs(I_ (1,i)) *1.0E-08 )then
I_ (2,i)=0
End If
!
WRITE(13,*)'Sectional Properties of Beam ',i
WRITE(13,*)'A_(1,1) ',' A_(2,2) ',' A_(3,3) ',' B_(1,1) ',' B_(2,2) ',' B_(3,3) ',' D_(1,1) ',' D_(2,2)','  D_(3,3)'  	 													
WRITE(13,*)'   N     ','      N  ','          N  ','       Nmm ','    Nmm  ','   Nmm  ','   Nmm2 ',' Nmm2 ',' Nmm2'
WRITE(13,9)A_(1,1,i) ,'  ',A_(2,2,i),' ',A_(3,3,i),' ',B_(1,1,i) ,'  ',B_(2,2,i),' ',B_(3,3,i),' ',D_(1,1,i) ,'  ',D_(2,2,i),' ',D_(3,3,i)
WRITE(13,*)' ' 
!
!
!
!
!
8 FORMAT(I2,A3,E6.1,A3,E6.1,A2,E6.1,A3,E6.1,A3,E6.1,A3,E6.1,A2,E6.1,A2,E6.1,A3,E8.2,E6.1)
9 FORMAT(E6.1,A3,E6.1,A2,E6.1,A3,E6.1,A3,E6.1,A3,E6.1,A2,E6.1,A2,E6.1,A3,E6.1)
!
!
!      
    END SUBROUTINE COMPOSITE_PROPERTIES
!--------------------------------------------------------------------------------------------  
!    
!
!-----------------------------------------------------------------
!
    SUBROUTINE INTERFACE_PROPERTIES (IELEM,A_)

    !DEFINING VARIABLES-------------------------------------------
    Double precision, dimension(:,:,:), intent (inout) :: A_
   Double precision::Cx,Cy
    !
    integer, intent (in) ::IELEM 
    !
    !-----------------------------------------------------------------
    OPEN(3,FILE="Interface_Spring.txt", FORM="FORMATTED",STATUS="UNKNOWN")     
    !Program reads from file BarStruc.inp
    OPEN(13,FILE="INPUT_CHECK.DAC",FORM="FORMATTED",STATUS="UNKNOWN")
    !
    write(13,*)'----------------------------'
    write(13,*)' '
    write(13,*)' '
    write(13,*)'****************************'
    write(13,*)' '
    !
    !
    !***********************************************************************************/
				    !****	Material and Cross-sectional properties	****/
!   
    !do  i=1,NE  !for each member
    !    
    read(3,*)Cx				    !read the modulus of elasticity
    read(3,*)Cy			    !read the cross-sectinal area

    write(13,*)'  Cohesivity of Spring in Horizontal  ',IELEM,'  = ',Cx	
    !write the modulus of elasticity for each member
    write(13,*)'  Cohesivity of Spring in Vertical ',IELEM,'  = ',Cy
    !write the he cross-sectinal area for each member
    !
    A_(1,1,IELEM)=Cx
    A_(2,2,IELEM)=Cy
    !
    !end do
    !
    write(13,*)'----------------------------'
    write(13,*)' '
    !
    !
    END SUBROUTINE INTERFACE_PROPERTIES
!
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------    
    SUBROUTINE  INTERFACE_SPRING_SPECTRAL_STIFFNESS (A_,K_SPECTRAL_ELEMENT,OM,IELEM)
    !DEFINING VARIABLES-------------------
    Double precision, dimension(:,:,:), intent (in) :: A_
    Double  complex, DIMENSION(6,6) , intent (inout) :: K_SPECTRAL_ELEMENT
    Double precision , intent (in) :: OM
    INTEGER , intent (in) ::IELEM
    Double precision KL,K1,GAMA,GAMABAR,ALFA,ALFABAR,BETA,BETABAR,DET,STF
    Double precision ak,exp_1,exp_2
    integer::i,j
    INTEGER :: P,W
    !-----------------------------------------------
    !
    !
    K_SPECTRAL_ELEMENT=0
    !
    !BAR ELEMENT*************************************
    !
    ak=om       ! Uncoupled Wave number 
    exp_1 = exp( (0.0,1) * (-1) * ak  )
    exp_2 = exp( (0.0,1) * ak  )  
    K_SPECTRAL_ELEMENT (1,1) = A_(1,1,IELEM) * (ak/dsin(ak))* dcos(ak)
    K_SPECTRAL_ELEMENT (1,4) = -1*A_(1,1,IELEM) * (ak/dsin(ak))		             
    K_SPECTRAL_ELEMENT (4,1) = K_SPECTRAL_ELEMENT (1,4)
    K_SPECTRAL_ELEMENT (4,4) = K_SPECTRAL_ELEMENT (1,1)
    !
    !BAR ELEMENT******************************
    !
    ak=om         ! Uncoupled Wave number 
    exp_1 = exp( (0.0,1) * (-1) * ak  )
    exp_2 = exp( (0.0,1) * ak  )  
    K_SPECTRAL_ELEMENT (2,2) = A_(2,2,IELEM) * (ak/dsin(ak))* dcos(ak)
    K_SPECTRAL_ELEMENT (2,5) = -1*A_(2,2,IELEM) * (ak/dsin(ak))		             
    K_SPECTRAL_ELEMENT (5,2) = K_SPECTRAL_ELEMENT (1,4)
    K_SPECTRAL_ELEMENT (5,5) = K_SPECTRAL_ELEMENT (1,1)
    !
    !
    !-----------------------------------------------
    !   
    END SUBROUTINE INTERFACE_SPRING_SPECTRAL_STIFFNESS 
!
!---------------------------------Delpasand 21/09/23---------------
!--------------------------------------------------------------------------------------------    
    SUBROUTINE  INTERFACE_SPRING_CONVENTIONAL_STIFFNESS (A_,K_SPECTRAL_ELEMENT,OM,IELEM)
    !DEFINING VARIABLES-------------------
    Double precision, dimension(:,:,:), intent (in) :: A_
    Double  complex, DIMENSION(6,6) , intent (inout) :: K_SPECTRAL_ELEMENT
    Double precision , intent (in) :: OM
    INTEGER , intent (in) ::IELEM
    Double precision KL,K1,GAMA,GAMABAR,ALFA,ALFABAR,BETA,BETABAR,DET,STF
    Double precision ak,exp_1,exp_2
    integer::i,j
    INTEGER :: P,W
    !-----------------------------------------------
    !
    K_SPECTRAL_ELEMENT=0
    !
    !BAR ELEMENT X direcrtion*************************************

    K_SPECTRAL_ELEMENT (1,1) = A_(1,1,IELEM) 
    K_SPECTRAL_ELEMENT (1,4) = -1*A_(1,1,IELEM)              
    K_SPECTRAL_ELEMENT (4,1) = K_SPECTRAL_ELEMENT (1,4)
    K_SPECTRAL_ELEMENT (4,4) = K_SPECTRAL_ELEMENT (1,1)
    !
    !BAR ELEMENT Y direction******************************
    K_SPECTRAL_ELEMENT (2,2) = A_(2,2,IELEM) 
    K_SPECTRAL_ELEMENT (2,5) = -1*A_(2,2,IELEM) 
    K_SPECTRAL_ELEMENT (5,2) = K_SPECTRAL_ELEMENT (1,4)
    K_SPECTRAL_ELEMENT (5,5) = K_SPECTRAL_ELEMENT (1,1)

    END SUBROUTINE INTERFACE_SPRING_CONVENTIONAL_STIFFNESS 
!----------------------------------Finished delpasand 21/09/23
!--------------------------------------------------------------------------------------------    
    SUBROUTINE EULER_BERNOULLI_ELEMENT_SPECTRAL_STIFFNESS   (A_,B_,D_,I_,L,K_SPECTRAL_ELEMENT,OM,IELEM)
    !DEFINING VARIABLES-------------------
    Double precision, dimension(:,:,:), intent (in) :: A_,B_,D_
    Double precision, dimension(:,:), intent (in) :: I_
    Double precision, dimension(:), intent (in) :: L	
    Double  complex, DIMENSION(6,6) , intent (inout) :: K_SPECTRAL_ELEMENT
    Double precision , intent (in) :: OM
    INTEGER , intent (in) ::IELEM
    Double precision KL,K1,GAMA,GAMABAR,ALFA,ALFABAR,BETA,BETABAR,DET,STF
    Double precision ak,exp_1,exp_2
    integer::i,j
    INTEGER :: P,W
    !-----------------------------------------------
    !!!!DOUBLE precision::ca,cb,cs,ka,kb,kr,r,s1,s2,c,d
    !!!!DOUBLE precision::omega,i2,i1,i0
    !!!!DOUBLE precision::D11,B11,a55,a11,B55
    !!!!double complex ::a,b,alpha,k_6,k_5,k_4,k_3,k_2,k_1
    !!!!INTEGER,PARAMETER:: DP=KIND(1.0D0)
    !!!!REAL(dp),DIMENSION(0:4):: p
    !!!!COMPLEX(dp),DIMENSION(4):: z
    !--------------------------------------
    !
    !
    !BEAM ELEMENT******************************
    K_SPECTRAL_ELEMENT=0
    !
    K1= ( I_(1,IELEM)*(OM**2)/( D_(1,1,IELEM) ) )**0.25    ! Uncoupled Wave number
    KL=K1*L(IELEM)  !**2 
    DET=1-dCOS(KL)*dCOSH(KL)
    STF=D_(1,1,IELEM)/(L(IELEM)**3)
    ALFA=(dCOS(KL)*dSINH(KL)+dSIN(KL)*dCOSH(KL))*(KL**3)/DET
    ALFABAR=(dSIN(KL)+dSINH(KL))*(KL**3)/DET    
    BETA=((-1)*dCOS(KL)*dSINH(KL)+dSIN(KL)*dCOSH(KL))*(KL)/DET
    BETABAR=((-1)*dSIN(KL)+dSINH(KL))*(KL)/DET
    GAMA=((-1)*dCOS(KL)+dCOSH(KL))*(KL**2)/DET
    GAMABAR=dSIN(KL)*dSINH(KL)*(KL**2)/DET
    !
    K_SPECTRAL_ELEMENT(2,2) = STF*ALFA
    K_SPECTRAL_ELEMENT(2,3) = STF*GAMABAR*L(IELEM)
    K_SPECTRAL_ELEMENT(2,5) = STF*(-1)*ALFABAR
    K_SPECTRAL_ELEMENT(2,6) = STF*GAMA*L(IELEM)
    K_SPECTRAL_ELEMENT(3,3) = STF*BETA*(L(IELEM)**2)
    K_SPECTRAL_ELEMENT(3,5) = STF*(-1)*GAMA*L(IELEM)                !!!!!! It was written gammabar in doyles paper which is wrong
    K_SPECTRAL_ELEMENT(3,6) = STF*BETABAR*(L(IELEM)**2)
    K_SPECTRAL_ELEMENT(5,5) = STF*ALFA
    K_SPECTRAL_ELEMENT(5,6) = STF*(-1)*GAMABAR*L(IELEM) 
    K_SPECTRAL_ELEMENT(6,6) = STF*BETA*(L(IELEM)**2)
    !
    K_SPECTRAL_ELEMENT(3,2) = K_SPECTRAL_ELEMENT(2,3)
    K_SPECTRAL_ELEMENT(5,2) = K_SPECTRAL_ELEMENT(2,5)
    K_SPECTRAL_ELEMENT(6,2) = K_SPECTRAL_ELEMENT(2,6)
    K_SPECTRAL_ELEMENT(5,3) = K_SPECTRAL_ELEMENT(3,5)
    K_SPECTRAL_ELEMENT(6,3) = K_SPECTRAL_ELEMENT(3,6)
    K_SPECTRAL_ELEMENT(6,5) = K_SPECTRAL_ELEMENT(5,6)
    !
    !
    !BAR ELEMENT*************************************
    ak=om*dsqrt(I_(1,ielem)/( A_(1,1,IELEM) ) )         ! Uncoupled Wave number 
    exp_1 = exp( (0.0,1) * (-1) * ak * L(ielem) )
    exp_2 = exp( (0.0,1) * ak * L(ielem) )  
    K_SPECTRAL_ELEMENT (1,1) = A_(1,1,IELEM)/L(IELEM) * (ak*L(IELEM)/dsin(ak*L(IELEM)))* dcos(ak*L(IELEM))
    K_SPECTRAL_ELEMENT (1,4) = -1*A_(1,1,IELEM)/L(IELEM) * (ak*L(IELEM)/dsin(ak*L(IELEM)))		             
    K_SPECTRAL_ELEMENT (4,1) = K_SPECTRAL_ELEMENT (1,4)
    K_SPECTRAL_ELEMENT (4,4) = K_SPECTRAL_ELEMENT (1,1)
    
      
    !-----------------------------------------------
    !   
    END SUBROUTINE EULER_BERNOULLI_ELEMENT_SPECTRAL_STIFFNESS
!
!--------------------------------------------------------------------------------------------
    Subroutine TIMOSHENKO_ELEMENT_SPECTRAL_STIFFNESS (A_,B_,D_,I_,L,K_SPECTRAL_ELEMENT,OM,IELEM)
    !
    Double  complex, DIMENSION(6,6) , intent (inout) :: K_SPECTRAL_ELEMENT
    !
    Double precision, dimension(:,:,:), intent (in) :: A_,B_,D_
    Double precision, dimension(:,:), intent (in) :: I_
    Double precision, dimension(:), intent (in) :: L	
    Double precision, intent (in) :: OM
    !
    INTEGER , intent (in) ::IELEM
    !
    Double precision :: GAMA,ALFA,BETA,GAMMA,DELTA,EPS,KSI,LAMDA,R,S,Q,T,PSI_A,LAMDA_A,LAMDA_B,S1,S2,S3,C1,C2,C3,ZETA,ETHA,LAMDA_C,SIGMA_A,CHI,K
    Double precision :: ak,exp_1,exp_2,EA,MU,EI,kGA
    integer::i,j,P,W
    !
    !
    !
    EA = A_(1,1,iELEM) 
    EI = D_(1,1,iELEM)
    MU = I_(1,iELEM)
    !kGA= A_(2,2,iELEM)
    !changed
    kGA= A_(3,3,iELEM)
    !
    !!!
    !
    K = EI / L(IELEM)
    LAMDA = ( ( MU*(L(IELEM)**4)*(OM**2) ) / EI )**0.25
    R = EI / ( EA * L(IELEM)**2 )
    S = EI / ( kGA* L(IELEM)**2 )
    !ENFORCING EULER-BRNOULI
    !R=0
    !S=0
    !
    Q = 1 - LAMDA**4 * R**2 * S**2
    T = 1
    PSI_A = ( (R**2)*(LAMDA*2) + (S**2)*(LAMDA**2) ) / (2*T)
    !
    LAMDA_A = LAMDA * SQRT ( ( PSI_A + SQRT ( PSI_A**2 + Q / T ) ) ) 
    !
    !
    S1 = SIN (LAMDA_A)
    C1 = COS (LAMDA_A)
    !
    IF ( Q > 0 ) THEN 
        J = 1
        LAMDA_B = LAMDA * SQRT ( (J*(-1*PSI_A + SQRT(PSI_A**2 + Q/T) ) ) )
        S2 = SINH(LAMDA_B)
        C2 = COSH(LAMDA_B)
    ENDIF
    !
    IF ( Q <= 0 )THEN
        J = -1
        LAMDA_B = LAMDA * SQRT ( (J*(-1*PSI_A + SQRT(PSI_A**2 + Q/T) ) ) )
        S2 = SIN(LAMDA_B)
        C2 = COS(LAMDA_B)
    ENDIF
    !
    !
    !
    ZETA = T - S**2 * LAMDA**4 / (LAMDA_A**2)
    ETHA = ZETA * (LAMDA_A/LAMDA_B) / (J*T + S**2 * LAMDA**4 / (LAMDA_B**2) )
    LAMDA_C = 0.5 * ( ETHA * LAMDA_A + LAMDA_B )
    SIGMA_A = ( ETHA * ( 1 - C1*C2 ) + 0.5* ( 1 - J*ETHA**2 )*S1*S2)/ LAMDA_C
    CHI = ( LAMDA_B - ETHA*LAMDA_A ) / (LAMDA_A + J*ETHA*LAMDA_B)
    !
    ALFA = (S1*C2 - J*ETHA*C1*S2) / SIGMA_A 
    BETA = ( J*ETHA*S2 - S1) / SIGMA_A
    GAMMA = ( ( C1*S2 + ETHA*S1*C2) * (LAMDA**4/(LAMDA_A*LAMDA_B )) ) /SIGMA_A
    DELTA = (C2 - C1)*ZETA*LAMDA_A / SIGMA_A
    EPS = (S2 + ETHA*S1) * ( LAMDA**4/(LAMDA_A*LAMDA_B ) ) / SIGMA_A
    KSI = 0.5 * (S1*S2 - CHI*(1-C1*C2) ) * ZETA * (LAMDA_A + J*ETHA*LAMDA_B)*(LAMDA_A/LAMDA_C)/SIGMA_A
    !
    K_SPECTRAL_ELEMENT(2,2) = GAMMA * K / ( L(IELEM)**2 )
    K_SPECTRAL_ELEMENT(2,3) = KSI * K / L(IELEM)
    K_SPECTRAL_ELEMENT(2,5) = -1*EPS * K / (L(IELEM)**2)
    K_SPECTRAL_ELEMENT(2,6) = DELTA * K / L(IELEM)
    K_SPECTRAL_ELEMENT(3,3) = ALFA * K
    K_SPECTRAL_ELEMENT(3,5) = -1 * DELTA * K / L(IELEM)
    K_SPECTRAL_ELEMENT(3,6) = BETA * K
    K_SPECTRAL_ELEMENT(5,5) = GAMMA * K / ( L(IELEM)**2 )
    K_SPECTRAL_ELEMENT(5,6) = -1*KSI * K / L(IELEM)  
    K_SPECTRAL_ELEMENT(6,6) = ALFA * K
    !
    K_SPECTRAL_ELEMENT(3,2) = K_SPECTRAL_ELEMENT(2,3)
    K_SPECTRAL_ELEMENT(5,2) = K_SPECTRAL_ELEMENT(2,5)
    K_SPECTRAL_ELEMENT(6,2) = K_SPECTRAL_ELEMENT(2,6)
    K_SPECTRAL_ELEMENT(5,3) = K_SPECTRAL_ELEMENT(3,5)
    K_SPECTRAL_ELEMENT(6,3) = K_SPECTRAL_ELEMENT(3,6)
    K_SPECTRAL_ELEMENT(6,5) = K_SPECTRAL_ELEMENT(5,6)
    !
    !BAR
    ak=om*dsqrt(I_(1,ielem)/( A_(1,1,IELEM) ) )         ! Uncoupled Wave number 
    exp_1 = exp( (0.0,1) * (-1) * ak * L(ielem) )
    exp_2 = exp( (0.0,1) * ak * L(ielem) )  
    K_SPECTRAL_ELEMENT (1,1) =    A_(1,1,IELEM)/L(IELEM) * (ak*L(IELEM)/dsin(ak*L(IELEM)))* dcos(ak*L(IELEM))
    K_SPECTRAL_ELEMENT (1,4) = -1*A_(1,1,IELEM)/L(IELEM) * (ak*L(IELEM)/dsin(ak*L(IELEM)))		             
    K_SPECTRAL_ELEMENT (4,1) = K_SPECTRAL_ELEMENT (1,4)
    K_SPECTRAL_ELEMENT (4,4) = K_SPECTRAL_ELEMENT (1,1)
    !
    !
   
    END SUBROUTINE TIMOSHENKO_ELEMENT_SPECTRAL_STIFFNESS
!
!-----------------------------------------------------
!*********************ISOTROPIC ENDS******************
!-----------------------------------------------------
!
!-----------------------------------------------------
!********************COMPOSITE BEGINS******************
!-----------------------------------------------------
!
!-----------------------------------------------------
SUBROUTINE JUN_COMPOSITE_ELEMENT_SPECTRAL_STIFFNESS ( A_, B_, D_, I_, L,  K_SPECTRAL_ELEMENT, OM, IELEM  )
!
    Double complex, DIMENSION(6,6), intent (inout) :: K_SPECTRAL_ELEMENT
    !
    Double precision, dimension(:,:,:), intent (in) :: A_,B_,D_	
    Double precision, dimension(:,:), intent (in) ::I_
    Double precision, dimension(:), intent (in) :: L
    Double precision, intent (in) :: OM
    !
    integer, intent (in) :: IELEM 
    !
    !
    DOUBLE PRECISION, PARAMETER :: PI = 3.1415927
    DOUBLE precision::omega,i3,i2,i1
    DOUBLE precision::D11,B11,a55,a11,B55
    DOUBLE precision::ETHA_0,ETHA_1,ETHA_2,ETHA_3,A1,A2,A3,Q,V
    DOUBLE COMPLEX :: CHI_1,CHI_2,CHI_3
    !trying complex just to be sure NOTHING IS ACTUALLY COMPLEX!
    DOUBLE complex, DIMENSION(6,6) :: R,H
    DOUBLE complex, DIMENSION(3) :: KAPPA,T,T_BAR,T_HAT,T_DOUBLE_HAT,T_DOUBLE_BAR
    !
    !
    Integer,DIMENSION(6):: IPVT
    integer::i,j,info
    integer :: p,w
    !
    !
    a11 = A_ (1,1,ielem)
    a55 = A_ (3,3,ielem)
    b11 = B_ (1,1,ielem)
    b55 = B_ (3,3,ielem)
    d11 = D_ (1,1,ielem)
    omega = OM
    i1 = I_ (1,ielem)
    i2 = I_ (2,ielem)
    i3 = I_ (3,ielem)
    !
    If ( B_ (1,1,ielem) /=0 .or. I_ (2,ielem) /=0 ) then             !----------------Coupled-Axial-Bending------------
    !
    ETHA_3 = -1 * (B11**2 - A11*D11)*(A55) 
    ETHA_2 = -1 * B11**2 * I1 * OMEGA**2 + ( A11*(D11*I1 + I3*A55) + (D11*I1 - 2*B11*I2)*A55 )* OMEGA**2
    ETHA_1 = -1 * A55 * I1 * A11 * OMEGA**2 + (D11*I1**2 - 2*B11*I1*I2 + A11*I1*I3 - (I2**2 - I1*I3)*A55 ) * OMEGA**4
    ETHA_0 = I1*OMEGA**4 * (-1*A55*I1 + (-1*I2**2 + I1*I3)*OMEGA**2)
    !
    A1 = ETHA_2/ETHA_3
    A2 = ETHA_1/ETHA_3
    A3 = ETHA_0/ETHA_3
    !
    Q = -1*A2 + (A1**2)/3 
    V = ACOS ( -27*A3 + 9*A1*A2 - 2*A1**3 / (2 * SQRT( (A1**2-3*A2)**3 ) ) )
    !
    CHI_1 = 2*SQRT(Q/3) * COS( V/3 ) - A1/3
    CHI_2 = 2*SQRT(Q/3) * COS( (2*PI+V)/3 ) - A1/3
    CHI_3 = 2*SQRT(Q/3) * COS( (4*PI+V)/3 ) - A1/3
    !
    KAPPA(1) = SQRT (CHI_1)
    KAPPA(2) = SQRT (CHI_2)
    KAPPA(3) = SQRT (CHI_3)
    !
    DO J=1,3
        T(J) = -1 * ( B11*KAPPA(J)**2 + I2*OMEGA**2 ) / (A11*KAPPA(J)**2 + I1*OMEGA**2)
        T_BAR(J) = -1*A55*KAPPA(J) / ( A55*KAPPA(J)**2 + I1*OMEGA**2 )
        !T(J) = -1 * ( B11*KAPPA(J)**2 + I2*OMEGA**2 ) / (A55*KAPPA(J)**2 + I1*OMEGA**2)
        !T_BAR(J) = -1*A55*KAPPA(J) / ( A11*KAPPA(J)**2 + I1*OMEGA**2 )
    ENDDO
    !
    DO J=1,3
        T_HAT(J) = A11*KAPPA(J)*T(J) + B11*KAPPA(J)
        T_DOUBLE_HAT(J) = -1 * ( B11*KAPPA(J)*T(J) + D11*KAPPA(J) )
        T_DOUBLE_BAR(J) = -1 * ( A55*KAPPA(J)*T_BAR(J) + A55 )
    ENDDO
    !
    !R MATRIX==============================================
    !
    DO J=1,3
        R(1,J) = T(J)
        R(2,J) = T_BAR(J)
        R(3,J) = 1
        R(4,J) = T(J) * EXP ( KAPPA(J) * L(IELEM) )
        R(5,J) = T_BAR(J) * EXP ( KAPPA(J) * L(IELEM) )
        R(6,J) = EXP ( KAPPA(J) * L(IELEM) )
    ENDDO
    !
    DO J=1,3
        R(1,J+3) = T(J)
        R(2,J+3) = -1*T_BAR(J)
        R(3,J+3) = 1
        R(4,J+3) = T(J) * EXP ( -1 * KAPPA(J) * L(IELEM) )
        R(5,J+3) = -1*T_BAR(J) * EXP ( -1 * KAPPA(J) * L(IELEM) )
        R(6,J+3) = EXP ( -1 * KAPPA(J) * L(IELEM) )
    ENDDO
    !
    !
    !H MATRIX===================================================
    DO J=1,3
        H(1,J) = -1*T_HAT(J)
        H(2,J) = T_DOUBLE_BAR(J)
        H(3,J) = T_DOUBLE_HAT(J)
        H(4,J) = T_HAT(J) * EXP ( KAPPA(J) * L(IELEM) )
        H(5,J) = -1 * T_DOUBLE_BAR(J) * EXP ( KAPPA(J) * L(IELEM) )
        H(6,J) = -1 * T_DOUBLE_HAT(J) * EXP ( KAPPA(J) * L(IELEM) )
    ENDDO
    !
    DO J=1,3
        H(1,J+3) = T_HAT(J)
        H(2,J+3) = T_DOUBLE_BAR(J)
        H(3,J+3) = -1 * T_DOUBLE_HAT(J)
        H(4,J+3) = -1 * T_HAT(J) * EXP ( -1 * KAPPA(J) * L(IELEM) )
        H(5,J+3) = -1 * T_DOUBLE_BAR(J) * EXP ( -1 * KAPPA(J) * L(IELEM) )
        H(6,J+3) = T_DOUBLE_HAT(J) * EXP ( -1 * KAPPA(J) * L(IELEM) )
    ENDDO
    !
    !
    !complex inverse
    CALL zgetrf(R,6,6,ipvt,info)
    CALL zgetri(R,6,6,ipvt,info)
    !
        IF ( info/=0 ) THEN
        print*,'Inverse operation is not successful'
        print*,'info=',info
        write (22,*)'Inverse operation is not successful'
        write (22,*)'info=',info
        print*,'STOP'
        pause
        END IF
    !
        !print*,'     '
        !print*,'===Jun et al Stiffness Matrix is being used as selected==='
        !print*,'     '
        !!
        !write(13,*),'     '
        !write(13,*),'===Jun et al Stiffness Matrix is being used as selected==='
        !write(13,*),'     '
    !
    K_SPECTRAL_ELEMENT = MATMUL (H,R)
    !
    Else               !----------------Un-Coupled-Axial-Bending------------
    !  
        !print*,'***********************************************************'
        !print*,'!!!Switched to Timoshenko Stiffness Matrix due to Symmetry!!!'
        !print*,'***********************************************************'
        !!
        !write(13,*),'***********************************************************'
        !write(13,*),'!!!Switched to Timoshenko Stiffness Matrix due to Symmetry!!!'
        !write(13,*),'***********************************************************'
    !    
    CALL TIMOSHENKO_ELEMENT_SPECTRAL_STIFFNESS (A_,B_,D_,I_,L,K_SPECTRAL_ELEMENT,OM,IELEM) 
    !CALL EULER_BERNOULLI_ELEMENT_SPECTRAL_STIFFNESS (A_,B_,D_,I_,L,K_SPECTRAL_ELEMENT,OM,IELEM)
    !    
    End If              !-------------------------------------
    !
    !
    !
    END SUBROUTINE JUN_COMPOSITE_ELEMENT_SPECTRAL_STIFFNESS
!  
!    
!-----------------------------------------------------
!
    SUBROUTINE GOPALA_COMPOSITE_ELEMENT_SPECTRAL_STIFFNESS ( A_, B_, D_, I_, L,  K_SPECTRAL_ELEMENT, OM, IELEM  )
    !DEFINING VARIABLES-------------------
    !
    Double complex, DIMENSION(6,6), intent (inout) :: K_SPECTRAL_ELEMENT
    !
    Double precision, dimension(:,:,:), intent (in) :: A_,B_,D_	
    Double precision, dimension(:,:), intent (in) ::I_
    Double precision, dimension(:), intent (in) :: L
    Double precision, intent (in) :: OM
    !
    integer, intent (in) :: IELEM 
    !
    !--------------------------------------
    !
    DOUBLE COMPLEX, DIMENSION (3,6) :: r_matrix,T1_0,T1_L,k_1,k_2,q0_r_10,q1_r_10,q0_r_1L,q1_r_1L,q1_r,q0_r
    DOUBLE COMPLEX, DIMENSION (6,6) :: lamda_0_0,T2
    DOUBLE COMPLEX, DIMENSION (6,6) :: lamda_0_L
    DOUBLE COMPLEX, DIMENSION (6,6) :: lamda_1_0
    DOUBLE COMPLEX, DIMENSION (6,6) :: lamda_1_L 
    DOUBLE complex, DIMENSION (6)   :: ks
    double complex ::a,b,alpha,k6,k5,k4,k3,k2,k1,kj
    !
    !
    DOUBLE precision, DIMENSION (6,6) ::T2_inv
    DOUBLE precision, DIMENSION (3,3) :: q0
    DOUBLE precision, DIMENSION (3,3) :: q1
    Double precision::det(2)
    Double precision::KL,GAMA,GAMABAR,ALFA,ALFABAR,BETA,BETABAR,STF
    Double precision::ak,exp_1,exp_2
    DOUBLE precision::ca,cb,cs,ka,kb,kr,r,s1,s2,c,d
    DOUBLE precision::omega,i2,i1,i0
    DOUBLE precision::D11,B11,a55,a11,B55
    !    !
    Integer,DIMENSION(6):: IPVT
    integer::i,j,info
    !
    !Added by Mohammad=========================
    ! This block shares I_ between ELEMENT_SPECTRAL_STIFFNESS and COMPOSITE_ELEMENT_SPECTRAL_STIFFNESS Subroutines
    !The shared variable between two subroutines can not be allocatable, so a large value, 10000, was selected for I_ dimension
    !----------------------------------------------------------------------------------------------------------------------
    !Variables to be used in the POLYROOTS MODULE
    INTEGER,PARAMETER:: DP=KIND(1.0D0)
    REAL(dp),DIMENSION(0:4):: p
    COMPLEX(dp),DIMENSION(4):: z,sqz
    !variables for exp_compkex
    double complex exp_compkex
    !Added by Mohammad=========================
    !
!--------------------------------------
!--------------------------------------
    !
    a11 = A_ (1,1,ielem)
    a55 = A_ (3,3,ielem)
    b11 = B_ (1,1,ielem)
    b55 = B_ (3,3,ielem)
    d11 = D_ (1,1,ielem)
    omega = OM
    i0 = I_ (1,ielem)
    i1 = I_ (2,ielem)
    i2 = I_ (3,ielem)
    !
    If ( B_ (1,1,ielem) /=0 .or. I_ (2,ielem) /=0 ) then             !----------------Coupled-Axial-Bending------------
    !
!Parameters 
!(Eq 12)
    ca = ( a11 / i0 ) ** 0.5
    cb = ( d11 * omega**2 / i0 ) ** 0.25
    cs = ( a55 / i0 ) ** 0.5 
!!(Eq 12)
    ka = omega / ca
    kb = omega / cb
    kr = omega / cs
!!(Eq 13)
    r = ( b11**2 / ( a11 * d11 ) ) ** 0.5
    s1 =  omega * ( i2 / a55 ) ** 0.5
    s2 = ( i1**2 / ( i0 * i2 ) ) ** 0.5
!!(Eq 27,28,29 )
    a = 1 - r**2
    b = ( 2 * r * s1 * s2 * ka * kb**2 ) / kr - (1 - r**2) * kr**2 - ( s1**2 * kb**4 ) / ( kr**2 ) - ka**2
    c = ka**2 * kr**2 - 2 * r * s1 * s2 * ka * kr * kb**2  - ( 1 - s1**2)*kb**4 + s1**2 * ( 1 - s2**2 ) * ka**2 *kb**4 / (kr**2)
    d = ( 1 - s1**2 * ( 1 - s2**2 ) ) * ka**2 * kb**4
    !  
    !solving the characteristic equation
    p(0) = d
    p(1) = c
    p(2) = b
    p(3) = a
    p(4) = 0
    !
    call QuarticSolver (p, z)
    !
    !write (*,*) 'Roots of aX^3 + bX^2 + cX +d = 0 are ', ( z(i), i=1,3 )
    !write (*,*) 'First pair of roots ak^6 + bk^4 + ck^2 +d = 0', sqrt(z(1)) , -1*sqrt(z(1))
    !write (*,*) 'Second pair of roots ak^6 + bk^4 + ck^2 +d = 0', sqrt(z(2)), -1*sqrt(z(2))
    !write (*,*) 'Third pair of roots ak^6 + bk^4 + ck^2 +d = 0', sqrt(z(3)) , -1*sqrt(z(3))
    alpha = z(1)
    sqz=sqrt ( z )
    
    k1 = sqrt ( z(1) )
    k2 = -1 * sqrt ( z(1) )

    k3 =  sqrt ( -1*( (b+alpha*a)/(2*a) )     + sqrt( ( (b+alpha*a)/(2*a) )**2 - (alpha*(b+alpha*a)+c) / a ) )
    k4 =  -1 *sqrt ( -1*( (b+alpha*a)/(2*a) ) + sqrt( ( (b+alpha*a)/(2*a) )**2 - (alpha*(b+alpha*a)+c) / a ) )
    
    k5 =  -1 * sqrt( -1*( (b+alpha*a)/(2*a) ) - sqrt( ( (b+alpha*a)/(2*a) )**2 - (alpha*(b+alpha*a)+c) / a ) )
    k6 =       sqrt( -1*( (b+alpha*a)/(2*a) ) - sqrt( ( (b+alpha*a)/(2*a) )**2 - (alpha*(b+alpha*a)+c) / a ) )
    

!ks = [k1,k2,k3,k4,k5,k6]                ! An array that containts the roots (Wavenumbers)
    ks(1) =k1
    ks(2) =k2
    ks(3) =k3
    ks(4) =k4
    ks(5) =k5
    ks(5) =k5
    ks(6) =k6
   
!!
!!
!!Lamda0 matrix at X=0 . Eq 34   
!!L is Element Length and x is the location we want to calculate
!!could be zero or L
    lamda_0_0 (1,1) = 1
    lamda_0_0 (2,2) = exp (  (0,1) * k2 * ( L(ielem) )  )
    lamda_0_0 (3,3) = 1
    lamda_0_0 (4,4) = exp (  (0,1) * k4 * ( L(ielem) )  )
    lamda_0_0 (5,5) = 1
    lamda_0_0 (6,6) = exp (  (0,1) * k6 * ( L(ielem) )  )
!!!
!!!Lamda0 matrix at X=L . Eq 34   
!!!L is Element Length and x is the location we want to calculate
!!!could be zero or L
    lamda_0_L (1,1) = exp (  (0,-1) * k1 * ( L(ielem) )  )
    lamda_0_L (2,2) = 1
    lamda_0_L (3,3) = exp (  (0,-1) * k3 * ( L(ielem) )  )
    lamda_0_L (4,4) = 1
    lamda_0_L (5,5) = exp (  (0,-1) * k5 * ( L(ielem) )  )
    lamda_0_L (6,6) = 1
!!!
!!!Lamda1_0 matrix at X=0 . Eq 47   
    lamda_1_0 (1,1) = (0,-1) * k1 
    lamda_1_0 (2,2) = (0,-1) * k2 * exp (  (0,1) * k2 * ( L(ielem) )  )
    lamda_1_0 (3,3) = (0,-1) * k3 
    lamda_1_0 (4,4) = (0,-1) * k4 * exp (  (0,1) * k4 * ( L(ielem) )  )
    lamda_1_0 (5,5) = (0,-1) * k5 
    lamda_1_0 (6,6) = (0,-1) * k6 * exp (  (0,1) * k6 * ( L(ielem) )  )
!!!
!!!Lamda1_L matrix at X=L . Eq 47  
    lamda_1_L (1,1) = (0,-1) * k1 * exp (  (0,-1) * k1 * ( L(ielem) )  )
    lamda_1_L (2,2) = (0,-1) * k2 
    lamda_1_L (3,3) = (0,-1) * k3 * exp (  (0,-1) * k3 * ( L(ielem) )  )
    lamda_1_L (4,4) = (0,-1) * k4
    lamda_1_L (5,5) = (0,-1) * k5 * exp (  (0,-1) * k5 * ( L(ielem) )  )
    lamda_1_L (6,6) = (0,-1) * k6
!!!
!!
!!R matrix - Equations 35,36,37,38
    r_matrix = 0.0
    do j=1,6
       kj = ks (j) !Take the jth root
       !
       if ( j==1 .or. j==2 ) then 
        r_matrix (1,j) = 1
        r_matrix (2,j) = ( (0,1) * kj * ( (r*kr**2*kj**2) / (ka*kb**2) - s1*s2*kr ) ) / ( (kj**2 - kr**2 ) * (1 - s1**2 + ( kr**2*kj**2 ) / ( kb**4 ) ) - kj**2 )
        r_matrix (3,j) = ( (kj**2 - kr**2 ) * ( (r*kr**2*kj**2) / (ka*kb**2) - s1*s2*kr ) ) / ( (kj**2 - kr**2 ) * (1 - s1**2 + ( kr**2*kj**2 ) / ( kb**4 ) ) - kj**2 )
       end if
       !
       if ( j >=3 .and. j<=6 ) then
        r_matrix (1,j) = ( (0,1)*(kj**2 - kr**2) * ( (s1*s2*ka**2)/kr - (r*ka*kj**2)/(kb**2) )) / ( kj*(kj**2 - ka**2) )
        r_matrix (2,j) = 1
        r_matrix (3,j) = (0,1) * (kr**2 - kj**2) / kj
       endif
      !
    end do
     
!    ! 
 !T1 matrix Equation 39
 !PLease use the multiplication module to get T1_0 and T1_L.   
!!!!!T1_0 = r_matrix * lamda_0_0
T1_0  =  MATMUL ( r_matrix , lamda_0_0 )
T1_L  =  MATMUL ( r_matrix , lamda_0_L )
!!
!! Eq 40
!! If you put T1_0 and T1_L on top of each other you'll get a 6*6 matrix called T2
Do i=1,3
    Do j=1,6
    T2(i,j)   = T1_0(i,j)
    T2(i+3,j) = T1_L(i,j)
    End Do
End Do
!!
!!
!! Eq 45
!! Matrix Q0 and Q1
!!
!!
    q0 = 0.0
    q0(2,3) = -1 * a55
!!
    q1 = 0.0
    q1(1,1) = a11
    q1(1,3) = -1 * b11
    q1(2,2) = a55
    q1(3,1) = b11
    q1(3,3) = -1 * d11
!
!! The K matrix is calculated according to Eq 46
!! @ is mutrix multplication
!! Please use frotran multiplication module to get k_1 and k_2
!
!
!!!!!!!k_1 = -1 * ( q0 @ r @ lamda_0_0 + q1 @ r @ lamda_1_0)
!!!!!!!k_2 = q0 @ r @ lamda_0_l + q1 @ r @ lamda_1_l
!
    q0_r     =  MATMUL ( q0 , r_matrix )
    q0_r_10  =  MATMUL ( q0_r , lamda_0_0 )
    !
    q1_r     =  MATMUL ( q1 , r_matrix )
    q1_r_10  =  MATMUL ( q1_r , lamda_1_0 )
    !
    k_1 =   -1*q0_r_10  -   1*q1_r_10
    !
    q0_r_1L  =  MATMUL ( q0_r , lamda_0_L )
    q1_r_1L  =  MATMUL ( q1_r , lamda_1_L )
    !
    k_2 =   q0_r_1L     +   q1_r_1L
          
    Do i=1,3
    Do j=1,6
    K_SPECTRAL_ELEMENT (i,j)   = k_1(i,j)
    K_SPECTRAL_ELEMENT (i+3,j) = k_2(i,j)
    End Do
    End Do
!
    !T2_inv=real(T2)
    !CALL DGEFA(T2,6,6,ipvt,info)	        ! Complex Inverse	
    !CALL dgedi(T2,6,6,ipvt,11,det)
    !
    CALL zgetrf(T2,6,6,ipvt,info)
    CALL zgetri(T2,6,6,ipvt,info)
    !
    IF ( info/=0 ) THEN
        print*,'Inverse operation is not successful'
        print*,'info=',info
        write (22,*)'Inverse operation is not successful'
        write (22,*)'info=',info
        print*,'STOP'
        pause
    END IF
    !
    !
    !print*,'  '
    !print*,'===Gopalakrishnan et al Stiffness Matrix is being used as selected==='
    !print*,'  '
    !!
    !write(13,*)'  '
    !write(13,*)'===Gopalakrishnan et al Stiffness Matrix is being used as selected==='
    !write(13,*)'  '
    !
    K_SPECTRAL_ELEMENT  =   MATMUL (  K_SPECTRAL_ELEMENT , T2 )
    !
    Else               !----------------Un-Coupled-Axial-Bending------------
    !
        !print*,'***********************************************************'
        !print*,'!!!Switched to Timoshenko Stiffness Matrix due to Symmetry!!!'
        !print*,'***********************************************************'
        !!
        !write(13,*),'***********************************************************'
        !write(13,*),'!!!Switched to Timoshenko Stiffness Matrix due to Symmetry!!!'
        !write(13,*),'***********************************************************'
    !    
    CALL TIMOSHENKO_ELEMENT_SPECTRAL_STIFFNESS (A_,B_,D_,I_,L,K_SPECTRAL_ELEMENT,OM,IELEM)       
    ! 
    !CALL EULER_BERNOULLI_ELEMENT_SPECTRAL_STIFFNESS (A_,B_,D_,I_,L,K_SPECTRAL_ELEMENT,OM,IELEM)
    !
    End If              !-------------------------------------
    !
    !
    END SUBROUTINE GOPALA_COMPOSITE_ELEMENT_SPECTRAL_STIFFNESS
!
!--------------------------------------------------------------------------------------------
!
    SUBROUTINE ASSEMB ( IELEM,ELEMENT_MATRIX,K_SPECTRAL_GLOBAL,BE,EN )
        !DEFINING VARIABLES---------------------------
        Double complex ,DIMENSION (6,6), intent (in)  :: ELEMENT_MATRIX
        Double complex, dimension (:,:), intent (inout) :: K_SPECTRAL_GLOBAL
        !
        integer, dimension(:), intent (in) :: Be, En
        INTEGER, intent (in) :: IELEM
        INTEGER NDOF,JDOFN
        integer INODE,NODEI,IDOFN,NROWS,NROWE,JNODE,NODEJ,NCOLS,NCOLE,NDOFN,NNODE
        !---------------------------------------------
        NDOFN = 3
        NNODE = 2
        !
             DO  INODE=1,NNODE
                 !
                  IF ( INODE == 1 ) THEN
                      NODEI = BE (IELEM)
                  ENDIF
                  !
                  IF ( INODE == 2 ) THEN
                      NODEI = EN (IELEM)
                  ENDIF
                  !
                  !NODEI=LNODS(IELEM,INODE)                                          
                  DO  IDOFN=1,NDOFN 
                      !
                      NROWS=(NODEI-1)*NDOFN+IDOFN                                     
                      NROWE=(INODE-1)*NDOFN+IDOFN                                     
    !                  ASLOD(NROWS)=ASLOD(NROWS)+ELOAD(IELEM,NROWE)                    
    !                                                                    
    !     ASSEMBLE THE ELEMENT STIFFNESS MATRICES
    !                                                                       
                      DO  JNODE =1,NNODE
                          !
                          IF ( JNODE == 1 ) THEN
                              NODEJ = BE (IELEM)
                          ENDIF
                          !
                          IF ( JNODE == 2 ) THEN
                              NODEJ = EN (IELEM)
                          ENDIF 
                          !
                          !NODEJ=LNODS(IELEM,JNODE)                                          
                          DO  JDOFN =1,NDOFN                                              
                              NCOLS=(NODEJ-1)*NDOFN+JDOFN                                     
                              NCOLE=(JNODE-1)*NDOFN+JDOFN                                     
                              K_SPECTRAL_GLOBAL(NROWS,NCOLS)= K_SPECTRAL_GLOBAL(NROWS,NCOLS) + ELEMENT_MATRIX(NROWE,NCOLE)  !real(ELEMENT_MATRIX(NROWE,NCOLE))
                          END DO
                          !
                      END DO
                  END DO
             END DO 
!             
    END SUBROUTINE ASSEMB
!
!--------------------------------------------------------------------------------------------
!    
    SUBROUTINE SOLVER (SFR,SFI,SKS,NF,Ur,Ui)
        !
        !defining variables
        integer, intent (in) :: NF
        Double complex, dimension(:), intent (out) :: Ur,Ui
        
        !
        !///////////////
		Double complex, dimension(:,:), intent (in):: SKS		            
		Double complex, dimension (:), intent (in):: SFr,SFi	
        Double complex, dimension(:,:), allocatable:: AG
        Double complex:: c,SCr,SCi
        integer m,k,i,j,p
        !-----------------------------------
	!		                        
							!************************************/
    !
	! Matrix to eliminate degrees of freedom corresponding to fixed nodes
    allocate (AG(NF,NF+2) ) 

        !Gauss elimination procedure
   

		do  i=1,NF
			do  j=1,NF
				AG(i,j)=SKS(i,j)
			end do
		end do


		do  i=1,NF
			AG(i,NF+1)=SFr(i)
			AG(i,NF+2)=SFi(i)
        end do
        
        
!!!        CALL CONSTRAINT(X,Y,Nn,ME,SSA,TA,MSAY)

        
        !Upper triangulization

		!double precision c

		do   k=1,NF-1				   
			do j=k+1,NF
			 c=AG(j,k)/AG(k,k)
				do  p=1,NF+2
				AG(j,p)=AG(j,p)-c*AG(k,p);
				end do
			end do
        end do
        
        if(AG(NF,NF)==0) then
		!write(3,*)"No Solution"
        print*,'No Solution'
		else 
		Ur(NF)=AG(NF,NF+1)/AG(NF,NF)
		Ui(NF)=AG(NF,NF+2)/AG(NF,NF)
        end if
        
        !Back Substitution
		   do i=NF-1,1,-1 
			SCr=0
			SCi=0
			 do  j=i+1,NF
				SCr=SCr+AG(i,j)*Ur(j);
				SCi=SCi+AG(i,j)*Ui(j);
			 end do
			 Ur(i)=1/AG(i,i)*(AG(i,NF+1)-SCr)
			 Ui(i)=1/AG(i,i)*(AG(i,NF+2)-SCi)

           end do
           
        !***********************************************************************************/
		!matrix vector multiplications Uj=  SUP x U

    END SUBROUTINE SOLVER
!
!--------------------------------------------------------------------------------------------
!   
    End Module Spectral_Stiffness
!
!
!----------------------------------------
!
!
!!Eq 34
!!You need to define a function that calculates the exponantial of a
!!complex number : e^(b+ic) = (e^b)(e^(ic)) = (e^b) * ((cos c) + i(sin c)) 
!double complex function exp_complex(i) result(j)
!  complex :: i ! input
!  b = real (i)
!  c = aimag (i)
!  j = exp (b) * ( cos(c) + (0,1) * sin(c) )
!    end function
!
!----------------------------------------
!
double complex function sqrt_complex (x)
    double precision :: x
    if ( x < 0 ) then
        sqrt_complex = sqrt (-1*x) * (0,1)
    endif
    if ( x >= 0 )then
        sqrt_complex = sqrt (x)
    endif
    end function