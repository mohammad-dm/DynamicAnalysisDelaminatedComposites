Module Wittrick_Williams
USE Spectral_Stiffness
USE CUSH
USE TRANSIENT
!
!USE INVERSE 
IMPLICIT NONE
!
Public:: Wittrick_Williams_Eigenvalue, SELECT_STIFFNESS, JJ_general, UPPER_TRIANGULAR, J_0_ISOTROPIC, J_0_ADAPTIVE,SEL_STIFF,BOUNDRY_IMPOSER
    contains
!    
    SUBROUTINE Wittrick_Williams_Eigenvalue (LAMD, NE, A_, B_, D_,I_, L, BE, EN, NF, SI, NJ, NUM_FREQ, B_Type,Me,RES,RES_tr,Ca,Ca_tr)  
        !CALCULATING NATURAL FREQUENCIES USING Wittrick Williams Method
        !DEFINING VARIABLES 
        !------------------------------------------------------------------
        !OUTPUT
        double precision,DIMENSION(:),ALLOCATABLE, INTENT (INOUT) :: LAMD  
        !------------------------------------------------------------------
        !INPUT
        integer, INTENT (IN) :: NE,NF,NJ,NUM_FREQ,Me
        integer, dimension(:), INTENT (IN) :: Be, En, SI  ,B_Type
        Double precision, dimension(:),INTENT (IN) :: L
        Double precision, dimension(:,:),INTENT (IN) ::I_,RES,RES_tr,Ca,Ca_tr
        Double precision, dimension(:,:,:), INTENT (IN) :: A_,B_,D_	
        !---------------------------------------------------------------------
        INTEGER:: NUM_ITERATION,j,k,i
        integer,DIMENSION(:),ALLOCATABLE::MN,SIGNCOUNT_J0
        DOUBLE PRECISION:: OM,OM_L, OM_U,error,eps,step
        Double complex ,DIMENSION(:,:),ALLOCATABLE :: K_G   !STIFNESS MATRIX AFTER BOUNDARY CONDITIONS ARE IMPOSED
        !---------------------------------------------------------------------
        OPEN(103,FILE="Natural_Frequencies.DAC", FORM="FORMATTED",STATUS="UNKNOWN") 
        !
        !
        WRITE (103,*)' '
        WRITE (103,*)'Number of Elements used for the Natural Frequency Analysis is  ',NE
        WRITE (103,*)' '
        WRITE (103,*)' '
        !
    Do i=1,NE
            !
    If (B_type(i)==1) Then
        write (103,*)'Beam type was selected as Euler-Bernoulli'
    Else If (B_type(i)==2) Then
        write (103,*)'Beam type was selected as Timoshenko'
    Else If (B_type(i)==3) Then
        write (103,*)'Beam type was selected as Jun et al Composite'
    Else If (B_type(i)==4) Then
         write (103,*)'Beam type was selected as Gopalakrishnan et al Composite'
    Else If (B_type(i)==5) then                    
         write (103,*) 'Interface element was used'
    End If
        !
        If ( B_ (1,1,i) ==0 .and. I_ (2,i) ==0 ) then             !----------------UnCoupled-Axial-Bending------------
        !
        print*,'******************Warning for Natural Frequency Calculations***************'
        print*,'Element ',i,' !!!Switched to Timoshenko Stiffness Matrix due to Symmetry!!!'
        print*,'***********************************************************'
        !
        write(103,*),'********************Warning********************************'
        write(103,*),'Element ',i,' !!!Switched to Timoshenko Stiffness Matrix due to Symmetry!!!'
        write(103,*),'***********************************************************'
        !
        Else
                !
        print*,'    '
        print*,'Element ',i,' ---calculated as NON-symmetric composite---'
        print*,'    '
        !
        write(103,*)'    '
        write(103,*),'Element ',i,'---calculated as NON-symmetric composite---'
        write(103,*)'    '
        !

        End If
        !
    END Do
        !
        ALLOCATE ( K_G (NF,NF),MN(NE),SIGNCOUNT_J0(NE) )
        !WRITE (*,*)"ENTER NUMBER OF MODES YOU WISH TO USE"
        !READ (*,*)NUM_FREQ
        allocate ( lamd (num_freq) )
        !NUM_ITERATION = 100
        eps=1.0E-04
        step=10
        OM_U=0
        NUM_ITERATION=0
                MN=2
                SIGNCOUNT_J0=1
        !
        DO I = 1,NUM_FREQ   !,  -1
            K=0
            OM=OM_U
            NUM_ITERATION=0
            
          DO while ( NUM_ITERATION==0 ) 
                K=K+1
                    OM_L = OM
                    OM_U = OM_U + step
               error=1
               J=0
            !                
            DO while (error>eps.or.J>I )
                !
                !
                if (NUM_ITERATION==0) then
                 OM =OM_U
                Else
                 OM = 0.5 * ( OM_U + OM_L )
                End if
                !
                CALL SELECT_STIFFNESS ( OM, NE, A_,B_,D_,I_,L, BE, EN, NF, NJ,K_G, SI,B_Type,Me,RES,RES_tr,Ca,Ca_tr )  !
                !
                !
                CALL JJ_general ( K_G, NF, A_,B_,D_,I_, L, NE, OM, J,MN,B_type,SIGNCOUNT_J0 )
                !
                !
                 K_G = 0.0
                IF ( J >= I ) THEN
                 OM_U = OM
                 NUM_ITERATION=NUM_ITERATION+1
                Else IF ( J < I ) THEN
                 OM_L = OM
                END IF 
                error=OM_U-OM_L
                !
            END DO
            !   
            !                
            !           
          END DO 
            !
                    LAMD (I) = OM
                    WRITE (103,*)"OMEGA FOR MODE", I, "=", OM
                    OM_L = OM_U
            !          
        END DO
            !
        WRITE (*,*)"Cyclic Natural Frequencies are stored in Natural Frequencies.DAC"
            !
    END SUBROUTINE Wittrick_Williams_Eigenvalue
!
!----------------------------------------------------------------------------------------
!
    SUBROUTINE JJ_general ( K_G, NF, A_,B_,D_,I_,L, NE, OM, J,MN,B_type,SIGNCOUNT_J0 )
        !INPUT--------------------------------------------------------------------------- 
        DOUBLE complex, DIMENSION (:,:), INTENT (IN) :: K_G
        Double precision, dimension(:,:,:), INTENT (IN) :: A_,B_,D_
        Double precision, dimension(:,:), INTENT (IN) :: I_
        Double precision, dimension(:), INTENT (IN) :: L
        DOUBLE PRECISION, INTENT (IN) :: OM
        !
        INTEGER, dimension(:), INTENT (INOUT) ::MN,SIGNCOUNT_J0
        INTEGER, dimension(:), INTENT (IN) ::B_type
        INTEGER, INTENT (IN) :: NF, NE
        INTEGER, INTENT (INOUT) :: J
        !
        DOUBLE PRECISION::c
        INTEGER:: SIGNCOUNT, J0
        !-------------------------------------------------------------------------------- 
        !-------------------------------------------------------------------------------- 
        CALL UPPER_TRIANGULAR (K_G, NF, SIGNCOUNT)
        !
        CALL J_0_ISOTROPIC (A_,B_,D_,I_,L,NE,OM,J0, B_type)
        !CALL J_0_ADAPTIVE (A_,B_,D_,I_,L,NE,OM,J0,MN,B_type,SIGNCOUNT_J0)
        !
        J = J0 + SIGNCOUNT
        !
            OPEN(352,FILE="TRY_SIGNCOUNT.DAC", FORM="FORMATTED",STATUS="UNKNOWN")
            WRITE(352,*)'OMEGA  ',OM
            WRITE(352,*)'SIGN_COUNT  ',SIGNCOUNT
        !  
    END SUBROUTINE JJ_general 
!
!----------------------------------------------------------------------------------------
!
    SUBROUTINE UPPER_TRIANGULAR (AG, NF, SIGN)
        !
        DOUBLE complex, DIMENSION (:,:), INTENT (IN) :: AG
        INTEGER, INTENT (out) :: SIGN
        !
        INTEGER NF,k,j,p ,i,info    
        INTEGER,DIMENSION(:),allocatable::ipvt
        DOUBLE complex, dimension (:,:), allocatable:: TRIANGULAR_AG
        DOUBLE complex, dimension (:), allocatable:: COC
        DOUBLE complex::c, SCON_PR,SCON,SCOC_PR,SCOC
        !
        ALLOCATE ( TRIANGULAR_AG (NF,NF), ipvt (NF), COC(NF) )
        !
        TRIANGULAR_AG = AG
        !Upper triangulization
        !DO I=1,NF
        !    DO J=1,NF
        !        TRIANGULAR_AG(I,J) = dreal(AG (I,J))            !*****Using only real part
        !    ENDDO
        !ENDDO
		do   k=1,NF-1				   
			do j=k+1,NF
			    c=TRIANGULAR_AG(j,k)/TRIANGULAR_AG(k,k)
				do  p=1,NF
				TRIANGULAR_AG(j,p)=TRIANGULAR_AG(j,p)-c*TRIANGULAR_AG(k,p);
				end do
			end do
        end do
        !
        ! 
                 SIGN = 0
       !
        !
        SCON=1
        !SCOC=1
        !------------------------------------|||||||||||||||||||||||||||||||------------------------------------
        DO I=1,NF
            !
            SCON_PR=SCON
            !SCOC_PR=SCOC
            !
            COC(I) = dcmplx ( real(TRIANGULAR_AG (I,I)), -1*imag(TRIANGULAR_AG (I,I)) )
            !
            !SCOC=SCOC*COC(I)
            !
            SCON=SCON_PR*TRIANGULAR_AG (I,I)
            !
            !IF ( REAL(SCON_PR)*REAL(SCON)<=0 .and. IMAG(SCON_PR)*IMAG(SCON)<=0 ) THEN !          !.and. REAL(SCOC_PR*SCOC)<0  .and. REAL(SCOC_PR*SCON)<0 .and. REAL(SCOC_PR*SCON)<0
            !    SIGN = SIGN + 1
            !ENDIF
        !    
            !IF (  real(TRIANGULAR_AG (I,I)) <= 0.0  .and. imag(TRIANGULAR_AG (I,I)) <= 0.0   ) THEN         ! Existance of Imaginary is questionable!!!
            !    SIGN = SIGN + 1
            !ENDIF
        !
            IF (   real(TRIANGULAR_AG (I,I)) <= 0.0    ) THEN
                SIGN = SIGN + 1
            ENDIF
        !
        !
        ENDDO
        !------------------------------------|||||||||||||||||||||||||||||||------------------------------------
        !       
        !
        !
        !
        !DO I=1,NF
        !    IF (  (real(TRIANGULAR_AG (I,I)) <= 0.0 ) ) THEN !.and.imag(TRIANGULAR_AG (I,I)) <= 0.0 )  ) THEN
        !        SIGN = SIGN + 1
        !    ENDIF
        !ENDDO
        !        
            ! if (SIGN>=1) then
            !    print*,'S=',SIGN
            !end if
        !
        !
    END SUBROUTINE UPPER_TRIANGULAR
!
!----------------------------------------------------------------------------------------
!
    SUBROUTINE J_0_ISOTROPIC (A_,B_,D_,I_,L,NE,OM,J0,B_type)
        !DEFINING VARIABLES---------- 
        INTEGER,intent(in):: NE
        INTEGER,intent(out):: J0
        !
        Double precision,intent(in):: OM
        Double precision, dimension(:,:,:),intent(in):: A_,B_,D_
        Double precision, dimension(:,:),intent(in)::I_
        Double precision, dimension(:),intent(in):: L
        !
        INTEGER IELEM,JA,JB,I, SIGN
        Double precision R,BEAM_WAVENUMBER, BAR_WAVENUMBER
        DOUBLE PRECISION, PARAMETER :: PI = 3.14159265359
        INTEGER, dimension(:), INTENT (IN) ::B_type
        !
        !----------------------------
        !Willams and Wittrick (1970) - An automatic computational procedure for calculating natural frequencies of skeletal structures
        !International Journal of Mechanical Sciences page 788
        !----------------------------
        j0 = 0
        DO IELEM=1,NE
            BAR_WAVENUMBER = om*dsqrt(I_(1,ielem)/(A_(1,1,ielem))) * L(IELEM) 
            !BEAM_WAVENUMBER = dsqrt(dSQRT( (DENSITY(IELEM)*A(IELEM)*(OM**2) ) / (E(IELEM) * II(IELEM)) )) * L(IELEM) 
            BEAM_WAVENUMBER =L(IELEM) * ( I_(1,IELEM)*(OM**2)/(D_(1,1,ielem)) )**0.25 !
            !CALCULATING JA
            JA = int(BAR_WAVENUMBER / PI)
            !CALCULATING JB
            I = int(BEAM_WAVENUMBER / PI)
            R = 1 - COS(BEAM_WAVENUMBER)*COSH(BEAM_WAVENUMBER)
            IF ( R >= 0.0 ) THEN
                SIGN = 1
            ELSE
                SIGN = -1
            ENDIF
            JB = I - ( 1 - SIGN* (-1)**I  )/2
            
            !Delpasand 21/-9/23
            IF ( B_TYPE (IELEM) == 5 .OR. B_TYPE (IELEM) == 6 )THEN
                JA = 0
                JB = 0
            ENDIF
            !--------------------
            
            J0 = J0 + JA + JB
            
        END DO
        
            !if (J0>=1) then
            !    print*,'At OM=',OM
            !    print*,'For IELEM=',IELEM
            !    print*,'J0=',J0
            !end if
            !
    END SUBROUTINE J_0_ISOTROPIC
!
!----------------------------------------------------------------------------------------
!
    SUBROUTINE J_0_ADAPTIVE (A_,B_,D_,I_,L,NE,OM,J0,MN,B_type,SIGNCOUNT_J0)
        !DEFINING VARIABLES---------- 
        Double precision, dimension(:,:,:),intent(in):: A_,B_,D_
        Double precision, dimension(:,:),intent(in):: I_
        Double precision, dimension(:),intent(in):: L
        Double precision,intent(in):: OM
        !
        INTEGER, dimension(:),intent(inout) ::MN,SIGNCOUNT_J0
        INTEGER, dimension(:),intent(in) ::B_type
        INTEGER,intent(in) :: NE
        INTEGER,intent(out):: J0
        !
        INTEGER:: IELEM, SIGNCOUNT_
        !
        !-------------------------------------------------------------------
        !
        Double precision, dimension(:,:,:), allocatable::A__,B__,D__
        Double precision, dimension(:,:), allocatable::I__
        Double precision, dimension(:), allocatable:: COOR_	,LE_
        Double precision:: Lc
        Double  complex ,DIMENSION(:,:), ALLOCATABLE :: K_SPECTRAL_GLOBAL_,K_G_
        Double  complex ,DIMENSION(6,6) :: K_SPECTRAL_ELEMENT_
        integer, dimension(:), allocatable:: Be_, En_, SI_,B__type
        integer::NJ_,IELEM_, WN,MNC,MNF,NF_,SIGNCOUNT_COARSE,SIGNCOUNT_FINE,SIGNCOUNT_PR
        !
        OPEN(13,FILE="INPUT_CHECK.DAC",FORM="FORMATTED",STATUS="UNKNOWN")
        !
        !
        WN=2      ! Mesh refinement factor  
        !----------------------------
        !
        !   FOR EACH CLAMPED-CLAMPED MEMBER
            j0 = 0
            DO IELEM=1,NE
            !
                SIGNCOUNT_PR=SIGNCOUNT_J0(IELEM)
            !   START A MESH NUMBER
            !
            !   PRODUCE A COARSE MESH ACCORDING TO MESH NUMBER 
            !        
            !If (B_type==1) then         ! Isotropic beam
            !    
            Lc=L(IELEM)
            MNC=MN(IELEM)
            !
            CALL PROPERTY_GENERATOR (NJ_,COOR_,BE_,EN_,NF_,Le_,SI_,MNC,Lc)           ! GENERATE MESH AND REFINED MESH
            !
            ALLOCATE ( K_SPECTRAL_GLOBAL_ (3*NJ_,3*NJ_) ,K_G_ (NF_,NF_),A__(3,3,MNC),B__(3,3,MNC),D__(3,3,MNC),I__(3,MNC) ,B__Type(MNC) )
            !
            A__(1,1,:)=A_(1,1,IELEM)
            B__(1,1,:)=B_(1,1,IELEM)
            D__(1,1,:)=D_(1,1,IELEM)
            I__(1,:)=I_(1,IELEM)
            !
            A__(3,3,:)=A_(3,3,IELEM)
            B__(3,3,:)=B_(3,3,IELEM)
            I__(2,:)=I_(2,IELEM)
            I__(3,:)=I_(3,IELEM)
            !
            B__Type(:)=B_Type(IELEM)
            !
                CALL SEL_STIFF ( OM, MNC, A__,B__,D__,I__, Le_, BE_, EN_, NF_, NJ_,K_G_, SI_,B__Type)           
            !
                    CALL UPPER_TRIANGULAR (K_G_, NF_, SIGNCOUNT_COARSE)
            !
                    DEALLOCATE ( K_SPECTRAL_GLOBAL_  ,K_G_ , Le_, BE_, EN_, SI_, A__, B__, D__,I__,B__Type )
            !
            !   PRODUCE A FINE MESH ACCORDING TO MESH NUMBER AND THE MESH FACTOR
            !        
            ! 
            Lc=L(IELEM)
            MNF=MN(IELEM)*WN
            !
            CALL PROPERTY_GENERATOR (NJ_,COOR_,BE_,EN_,NF_,Le_,SI_,MNF,Lc)           ! GENERATE MESH AND REFINED MESH
            !
            !----------------------------------------------------------------------------------------------
            !
            ALLOCATE ( K_SPECTRAL_GLOBAL_ (3*NJ_,3*NJ_) ,K_G_ (NF_,NF_),A__(3,3,MNF),B__(3,3,MNF),D__(3,3,MNF),I__(3,MNF),B__Type(MNF) )
            !
            A__(1,1,:)=A_(1,1,IELEM)
            B__(1,1,:)=B_(1,1,IELEM)
            D__(1,1,:)=D_(1,1,IELEM)
            I__(1,:)=I_(1,IELEM)
            !
            A__(3,3,:)=A_(3,3,IELEM)
            B__(3,3,:)=B_(3,3,IELEM)
            I__(2,:)=I_(2,IELEM)
            I__(3,:)=I_(3,IELEM)
            !
            B__Type(:)=B_Type(IELEM)
            !
                CALL SEL_STIFF ( OM, MNF, A__,B__,D__,I__, Le_, BE_, EN_, NF_, NJ_,K_G_, SI_,B__Type )          
            !
                    CALL UPPER_TRIANGULAR (K_G_  , NF_,   SIGNCOUNT_FINE)
            !
                    DEALLOCATE ( K_SPECTRAL_GLOBAL_  ,K_G_ , Le_, BE_, EN_, SI_, A__, B__, D__ ,I__,B__Type)
            !
            !----------------------------------------------------------                    
            !        
                    IF (SIGNCOUNT_COARSE==SIGNCOUNT_FINE) THEN
                            SIGNCOUNT_=SIGNCOUNT_COARSE
                    ELSE IF (SIGNCOUNT_COARSE<SIGNCOUNT_FINE) THEN
                            SIGNCOUNT_=SIGNCOUNT_FINE
            !                UPDATE MESH NUMBER USING THE MESH FACTOR
                             MN(IELEM)=MN(IELEM)*WN
            !
                     WRITE (13,*)'------------------------------------- '
                     WRITE (13,*)' AT OMEGA_TRIAL ',OM
                     WRITE (13,*)' SIGNCOUNT_COARSE < SIGNCOUNT_FINE '
                     WRITE (13,*)' SIGNCOUNT_COARSE ',SIGNCOUNT_COARSE, ' < ',SIGNCOUNT_FINE, ' SIGNCOUNT_FINE'

                     WRITE (13,*)' Number of Elements for J0 count is updated in element',IELEM
                     WRITE (13,*)' J0 count is now based on ',MN(IELEM),'elements'
                     WRITE (13,*)'------------------------------------- '
            !
                    ELSE IF (SIGNCOUNT_COARSE>SIGNCOUNT_FINE) THEN
            !
            !        UPDATE MESH NUMBER USING THE MESH FACTOR
                             MN(IELEM)=MN(IELEM)*WN
            !
                     WRITE (13,*)'------------------------------------- '
                     WRITE (13,*)' AT OMEGA_TRIAL ',OM
                     WRITE (13,*)' SIGNCOUNT_COARSE',SIGNCOUNT_COARSE, ' > ',SIGNCOUNT_FINE, ' SIGNCOUNT_FINE'
                     WRITE (13,*)' **Something went wrong** '
            !
                    END IF
            !
                    j0=j0+SIGNCOUNT_
            !-----------------------------------------------------------------
                    IF (SIGNCOUNT_>SIGNCOUNT_PR) THEN
            !        UPDATE MESH NUMBER USING THE MESH FACTOR
                             MN(IELEM)=MN(IELEM)*WN
            !
                     WRITE (13,*)' '
                     WRITE (13,*)'------------------------------------- '
                     WRITE (13,*)' AT OMEGA_TRIAL ',OM
                     WRITE (13,*)' SIGNCOUNT ',SIGNCOUNT_  ,' > ', SIGNCOUNT_PR, ' PREVIOUS SIGNCOUNT'
                     WRITE (13,*)' Number of Elements for J0 count is updated in element',IELEM
                     WRITE (13,*)' J0 count is now based on ',MN(IELEM),'elements'
                     WRITE (13,*)'------------------------------------- '
                     WRITE (13,*)' '
            !
                    END IF
                    !
                     SIGNCOUNT_J0(IELEM)=SIGNCOUNT_
            !        
            !if (J0>=1) then
            !    print*,'At OM=',OM
            !    print*,'For IELEM=',IELEM
            !    print*,'SIGNCOUNT=',SIGNCOUNT_J0
            !end if
            !
            END DO       !  IELEM
            !
            !  
            !----------------------------
            !
    END SUBROUTINE J_0_ADAPTIVE
!
!----------------------------------------------------------------------------------------
!
    SUBROUTINE SELECT_STIFFNESS ( OM, NE, A_,B_,D_,I_,L, BE, EN, NF, NJ,K_G, SI,B_Type ,Me,RES,RES_tr,Ca,Ca_tr) 
        !----------------------------------------------------------
        !INPUT
        Double precision, dimension(:,:,:), INTENT (IN) :: A_,B_,D_	
        Double precision, dimension(:,:), INTENT (IN) :: I_,RES,RES_tr,Ca,Ca_tr
        Double precision, dimension(:), INTENT (IN) :: L
        DOUBLE PRECISION, INTENT (IN) :: OM 
        !
        integer, dimension(:), INTENT (IN) :: Be, En, SI,B_Type 
        integer, INTENT (IN) :: NE ,NF,NJ ,Me
        !----------------------------------------------------------
        !OUTPUT
        Double complex ,DIMENSION(:,:), INTENT (INOUT) :: K_G
        !----------------------------------------------------------
        Double complex ,DIMENSION(:,:),ALLOCATABLE :: K_SPECTRAL_GLOBAL,CZC,ZCS,ZC
        Double  complex , DIMENSION(6,6) :: K_SPECTRAL_ELEMENT
        integer::IELEM,I,J
        !----------------------------------------------------------
        !
        ALLOCATE ( K_SPECTRAL_GLOBAL (3*NJ,3*NJ) )
        K_SPECTRAL_GLOBAL = 0.0
        !       
OPEN(22,FILE="K_matrix.txt", FORM="FORMATTED",STATUS="UNKNOWN") 
write (22,*) '  '
write (22,*) '///////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\'  
write (22,*) '----STIFFNESS MATRIX FOR NATURAL FREQUENCY ANALYSIS----'
write (22,*) '            ----Column-by-column---- '
write (22,*) '  '

!        
        DO IELEM = 1,NE
        K_SPECTRAL_ELEMENT = 0.0
        !
            If (B_type(IELEM)==1) then
                CALL EULER_BERNOULLI_ELEMENT_SPECTRAL_STIFFNESS ( A_, B_, D_, I_, L, K_SPECTRAL_ELEMENT, OM, IELEM )
                write (22,*) '----EULER-BERNOULLI BEAM FOR OMEGA----',OM
            Else If (B_type(IELEM)==2) then    
                CALL TIMOSHENKO_ELEMENT_SPECTRAL_STIFFNESS      ( A_, B_, D_, I_, L, K_SPECTRAL_ELEMENT, OM, IELEM )  
                write (22,*) '----TIMOSHENKO BEAM FOR OMEGA----',OM
            Else If (B_type(IELEM)==3) then    
                CALL JUN_COMPOSITE_ELEMENT_SPECTRAL_STIFFNESS   ( A_, B_, D_, I_, L, K_SPECTRAL_ELEMENT, OM, IELEM )
                write (22,*) '----JUN COMPOSITE BEAM FOR OMEGA----',OM
            Else If (B_type(IELEM)==4) then                    
                CALL GOPALA_COMPOSITE_ELEMENT_SPECTRAL_STIFFNESS( A_, B_, D_, I_, L, K_SPECTRAL_ELEMENT, OM, IELEM )  
                write (22,*) '----GOPALAKRISHNAN COMPOSITE BEAM FOR OMEGA----',OM
            Else If (B_type(IELEM)==5) then                    
                !CALL INTERFACE_SPRING_SPECTRAL_STIFFNESS    (A_,K_SPECTRAL_ELEMENT,OM,IELEM)     
                CALL INTERFACE_SPRING_CONVENTIONAL_STIFFNESS    (A_,K_SPECTRAL_ELEMENT,OM,IELEM)            
                write (22,*) 'INTERFACE ELEMENT FOR OMEGA:',OM
            End If 
        !
write (22,*) 'Element Number:',IELEM
write (22,*) '  '
    do i=1,6
        do j=1,6
            write (22,*) K_SPECTRAL_ELEMENT(i,j)
        enddo
        write (22,*)
    enddo
        !
    write (22,*) '-------------------------------------------------'
    write (22,*) '-------------------NATURAL------------------'
    write (22,*) '-------------------------------------------------'
    write (22,*) ' '     
        !
            CALL ASSEMB ( IELEM, K_SPECTRAL_ELEMENT, K_SPECTRAL_GLOBAL, BE, EN)  
        ENDDO
        !
        
!!********************************************14/7/2021************************************************
ALLOCATE(ZC(ME,3*NJ),ZCS(ME,ME),CZC(NF,ME))
!
!CALL DGEMMD( ME, 3*NJ, 3*NJ, RES_tr,K_SPECTRAL_GLOBAL ,ZC )  
ZC=MATMUL(RES_tr,K_SPECTRAL_GLOBAL)
!
!CALL DGEMMD( ME, ME, 3*NJ, ZC,RES ,ZCS ) 
ZCS=MATMUL(ZC,RES)
!
!CALL DGEMMD( NF, ME, ME, Ca_tr,ZCS ,CZC )  
CZC=MATMUL(Ca_tr,ZCS)
!
!CALL DGEMMD( NF, NF, ME, CZC,Ca ,K_G ) 
K_G=MATMUL(CZC,Ca)
!
DEALLOCATE(ZC,ZCS,CZC)

!!********************************************14/7/2021************************************************
!
!
    END SUBROUTINE SELECT_STIFFNESS
!!--------------------------------------------------
     SUBROUTINE SEL_STIFF ( OM, NE, A_,B_,D_,I_,L, BE, EN, NF, NJ,K_G, SI,B_Type ) 
        !----------------------------------------------------------
        !INPUT
        Double precision, dimension(:,:,:), INTENT (IN) :: A_,B_,D_	
        Double precision, dimension(:,:), INTENT (IN) :: I_
        Double precision, dimension(:), INTENT (IN) :: L
        DOUBLE PRECISION, INTENT (IN) :: OM 
        !
        integer, dimension(:), INTENT (IN) :: Be, En, SI ,B_Type 
        integer, INTENT (IN) :: NE ,NF,NJ
        !----------------------------------------------------------
        !OUTPUT
        Double complex ,DIMENSION(:,:), INTENT (INOUT) :: K_G
        !----------------------------------------------------------
        Double complex ,DIMENSION(:,:),ALLOCATABLE :: K_SPECTRAL_GLOBAL,CZC,ZCS,ZC
        Double  complex , DIMENSION(6,6) :: K_SPECTRAL_ELEMENT
        integer::IELEM,I,J
        !----------------------------------------------------------
        !
        ALLOCATE ( K_SPECTRAL_GLOBAL (3*NJ,3*NJ) )
        K_SPECTRAL_GLOBAL = 0.0
        !       
OPEN(22,FILE="K_matrix.txt", FORM="FORMATTED",STATUS="UNKNOWN") 
write (22,*) '  '
write (22,*) '///////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\'  
write (22,*) '----STIFFNESS MATRIX FOR NATURAL FREQUENCY ANALYSIS----'
write (22,*) '            ----Column-by-column---- '
write (22,*) '  '

!        
        DO IELEM = 1,NE
        K_SPECTRAL_ELEMENT = 0.0
        !
            If (B_type(IELEM)==1) then
                CALL EULER_BERNOULLI_ELEMENT_SPECTRAL_STIFFNESS ( A_, B_, D_, I_, L, K_SPECTRAL_ELEMENT, OM, IELEM )
                write (22,*) '----EULER-BERNOULLI BEAM FOR OMEGA----',OM
            Else If (B_type(IELEM)==2) then    
                CALL TIMOSHENKO_ELEMENT_SPECTRAL_STIFFNESS      ( A_, B_, D_, I_, L, K_SPECTRAL_ELEMENT, OM, IELEM )  
                write (22,*) '----TIMOSHENKO BEAM FOR OMEGA----',OM
            Else If (B_type(IELEM)==3) then    
                CALL JUN_COMPOSITE_ELEMENT_SPECTRAL_STIFFNESS   ( A_, B_, D_, I_, L, K_SPECTRAL_ELEMENT, OM, IELEM )
                write (22,*) '----JUN COMPOSITE BEAM FOR OMEGA----',OM
            Else If (B_type(IELEM)==4) then                    
                CALL GOPALA_COMPOSITE_ELEMENT_SPECTRAL_STIFFNESS( A_, B_, D_, I_, L, K_SPECTRAL_ELEMENT, OM, IELEM )  
                write (22,*) '----GOPALAKRISHNAN COMPOSITE BEAM FOR OMEGA----',OM
            Else If (B_type(IELEM)==5) then                    
                !CALL INTERFACE_SPRING_SPECTRAL_STIFFNESS    (A_,K_SPECTRAL_ELEMENT,OM,IELEM)  
                CALL INTERFACE_SPRING_CONVENTIONAL_STIFFNESS    (A_,K_SPECTRAL_ELEMENT,OM,IELEM)            
                write (22,*) 'INTERFACE ELEMENT FOR OMEGA:',OM
            End If 
        !
write (22,*) 'Element Number:',IELEM
write (22,*) '  '
    do i=1,6
        do j=1,6
            write (22,*) K_SPECTRAL_ELEMENT(i,j)
        enddo
        write (22,*)
    enddo
        !
    write (22,*) '-------------------------------------------------'
    write (22,*) '-------------------NATURAL------------------'
    write (22,*) '-------------------------------------------------'
    write (22,*) ' '     
        !
            CALL ASSEMB ( IELEM, K_SPECTRAL_ELEMENT, K_SPECTRAL_GLOBAL, BE, EN)  
        ENDDO
        !
        CALL BOUNDRY_IMPOSER( K_SPECTRAL_GLOBAL, NJ, NF, SI, K_G)

        !
    END SUBROUTINE SEL_STIFF
!
!----------------------------------------------------------------------------------------
!
    SUBROUTINE PROPERTY_GENERATOR (NJ,COOR,BE,EN,NF,Le,SI,NE,Lc)
    !DEFINING VARIABLES-------------------------------------------
    !
    Double precision, dimension(:), allocatable, intent (out) :: Le, Coor	
    Double precision, intent (in) :: Lc	
    !
    integer, dimension(:), allocatable, intent (out) :: Be,En,SI
    integer, intent (out) :: NF,NJ
    integer, intent (in) ::  NE
    !
    integer :: i
    !
    !-----------------------------------------------------------------
    !
    !
    !*********************************************************************************/
    !
    NJ=NE+1
    !
    !***********************************************************************************/
			    !****	Joint Coordinates (Geometry of the strucutre)  ****/
    allocate (Coor(NJ))		    !Coor is of dimension NJ
    !
    do i=1,NJ
    Coor(i)	= Real((i-1))/NE*Lc		!Read Coordinates of the joints  
    end do
    !
    !only one coordinate along the longitudinal axis of the bar is needed for every joint
    !***********************************************************************************/
            		!****	Connectivity Matrices	****/
    allocate (Be(NE), En(NE) )
    !
    do  i=1,NE          !for each member
    Be(i)=i	        	!read the joint number connected to the first end of the member
    En(i)=i+1		        !read the joint number connected to the secondend of the member
    end do
    !
    !***********************************************************************************/
				    !****	Material and Cross-sectional properties	****/
    allocate( Le(NE) )
    !
    !   
    !
    !-----------------------------------------------------    
    ! Length of the member
    do i=1,NE			!for each member
    Le(i)=dsqrt( (Coor(En(i))-Coor(Be(i)))*(Coor(En(i))-Coor(Be(i))) ) 
    ! Length of the members is the difference between the second and first coordinates 
    end do
    !
    !
    !
    allocate (SI (3*NJ) )
    !
    do i=1,3*NJ
	SI(i)=0		                            
    end do
    !
    SI(1)=1		                            
    SI(2)=1	
    SI(3)=1
    !
    SI(3*NJ-2)=1	
    SI(3*NJ-1)=1		                            
    SI(3*NJ  )=1	
    !
    NF=3*NJ-6
    END SUBROUTINE PROPERTY_GENERATOR
!------------------------------------------------------------------------
    SUBROUTINE BOUNDRY_IMPOSER ( MATRIX, NJ , NF, SI, SKS)
        !DEFINING VARIABLES--------------------
        Double complex ,DIMENSION(:,:), INTENT(IN) :: MATRIX
        Double complex, dimension(:,:), intent(out):: SKS
        integer, dimension(:), INTENT(in):: SI
        INTEGER,INTENT(IN):: NF,NJ
       !
	    Double precision, dimension(:,:), allocatable:: SUP,SUPtr
        Double complex, dimension(:,:), allocatable::SK
        integer m,i,j,k            
	    !--------------------------------------  
        allocate ( SUP(3*NJ,NF) ) 
        allocate ( SK(NF, 3*NJ) )	
	    allocate ( SUPtr(NF,3*NJ) )
        !--------------------------------------
        SKS = 0.0
        do i=1,3*NJ
	        do j=1,NF
	    	    SUP(i,j)=0
	        end do
        end do
        m=0
	
        do i=1,3*NJ	
		    if (SI(i)==0) then	
		    m=m+1
		    SUP(i,m)=1
		    end if
        end do
    
        do i=1,3*NJ
			    do j=1,NF
				    SUPtr(j,i)=SUP(i,j)
			    end do
        end do

			    do  i=1,NF
			    do  k=1,3*NJ
			    SK(i,k)=0
			    end do
			    end do

		    do i=1,NF
			    do  k=1,3*NJ 
				    do j=1,3*NJ
				    SK(i,k)=SK(i,k)+SUPtr(i,j)*MATRIX(j,k)
				    end do
		    end do
	      end do
	    do   i=1,NF
		    do   k=1,NF
			    SKS(i,k)=0
			    end do
	    end do
		    do  i=1,NF
			    do  k=1,NF
				    do j=1,3*NJ
				    SKS(i,k)=SKS(i,k)+SK(i,j)*SUP(j,k)
				    end do
			    end do
		    end do
    END SUBROUTINE BOUNDRY_IMPOSER
   !*******************************************************************************************       
!----------------------------------------------------------------------------------------
!
    SUBROUTINE Wittrick_Williams_Eigenvector ( V, NE, A_,B_,D_,I_, L, BE, EN, NF, SI, NJ, LAMD, NUM_FREQ,B_Type,Me,RES,RES_tr,Ca,Ca_tr )
        !INPUT---------------------------------------------------------------------------
        Double precision, dimension(:,:,:), INTENT (IN) :: A_,B_,D_	
        Double precision, dimension(:,:), INTENT (IN) :: I_,RES,RES_tr,Ca,Ca_tr
        Double precision, dimension(:), INTENT (IN) :: L
        !
        integer, dimension(:), INTENT (IN) :: B_Type,SI,Be, En 
        integer, INTENT (IN) :: NF,NJ,NE,Me
        !------------------------------------------------------------------------------
        double precision,DIMENSION(:,:), ALLOCATABLE :: V
        Double precision, dimension(:),allocatable :: LAMD
        !------------------------------------------------------------------------------
        INTEGER NUM_FREQ,idof,i
        DOUBLE PRECISION OM
        DOUBLE complex, DIMENSION (:,:), ALLOCATABLE :: K_G
        DOUBLE complex, DIMENSION (:,:), ALLOCATABLE :: K_G_R
        DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: F, D, REAL_K_G,    REAL_K_G_R
        INTEGER IROW,ICOL,JROW,JCOL,index,ifreq,flag
        double precision rcond, rcond_original_k, rcond_reduced_k
        !------------------------------------------------------------------------------------
        !OPEN(100,FILE="EigenFrequencies.DAC", FORM="FORMATTED",STATUS="UNKNOWN") 
        !OPEN(101,FILE="EigenVectors.DAC", FORM="FORMATTED",STATUS="UNKNOWN") 
        allocate ( v (nf,num_freq) )
        ALLOCATE ( K_G (NF,NF) ,     K_G_R (NF-1,NF-1) )
        ALLOCATE ( REAL_K_G(NF,NF),   REAL_K_G_R(NF-1,NF-1) )
        ALLOCATE ( F(NF,NF), D(NF,NF) )
        v = 0.0
        DO IFREQ = 1, NUM_FREQ
            !WRITE (101,*)"MODE NUMBER:", IFREQ
            D = 0.0
            F = 0.0
            OM = LAMD (IFREQ)
            CALL SELECT_STIFFNESS ( OM, NE, A_,B_,D_,I_, L, BE, EN, NF, NJ,K_G, SI,B_Type,Me,RES,RES_tr,Ca,Ca_tr )
            !
            REAL_K_G=dreal(K_G)     !*****Using only real part
            CALL DSYSVX( NF , 1, REAL_K_G, NF, f, d, rcond )
            rcond_original_k = rcond
            IDOF = 1
            FLAG = 0
            DO WHILE ( IDOF <= NF .AND. FLAG == 0 )
                JROW = 0
                JCOL = 0
                DO IROW=1,NF
                    DO ICOL=1,NF
                        IF ( IROW < IDOF .AND. ICOL < IDOF )THEN
                            JROW = IROW
                            JCOL = ICOL
                        ENDIF
                        IF ( IROW < IDOF .AND. ICOL > IDOF )THEN
                            JROW = IROW
                            JCOL = ICOL - 1
                        ENDIF
                        IF ( IROW > IDOF .AND. ICOL < IDOF )THEN
                            JROW = IROW - 1
                            JCOL = ICOL
                        ENDIF
                        IF ( IROW > IDOF .AND. ICOL > IDOF )THEN
                            JROW = IROW - 1  
                            JCOL = ICOL - 1 
                        ENDIF
                        IF ( IROW /= IDOF .AND. ICOL /= IDOF  )THEN
                            K_G_R ( JROW, JCOL ) = K_G ( IROW, ICOL )
                        ENDIF
                    ENDDO
                ENDDO
                !WRITE (*,*)
                DO IROW=1,NF
                    IF ( IROW < IDOF ) THEN
                        INDEX = IROW
                    ENDIF
                    IF ( IROW > IDOF ) THEN
                        INDEX = IROW - 1
                    ENDIF
                    IF (IROW .NE. IDOF )THEN
                        F (INDEX,1) = -1 * K_G (IROW,IDOF)
                    ENDIF
                ENDDO
                !
                REAL_K_G_R=dreal(K_G_R)     !*****Using only real part
                CALL DSYSVX( NF-1 , 1, REAL_K_G_R, NF-1, f, d, rcond )
                rcond_reduced_k = rcond
                if (  rcond_reduced_k / rcond_original_k  > 10**4 ) then
                    IROW = 1 
                    JROW = 1
                    FLAG = 1
                    DO WHILE ( JROW <= NF )
                        IF ( JROW == IDOF )THEN
                            V(JROW,IFREQ) = 1.0
                            JROW = JROW + 1
                        ELSE
                            V(JROW,IFREQ) = D (IROW,1) 
                            IROW = IROW + 1
                            JROW = JROW + 1
                        ENDIF                        
                    ENDDO
                endif
                IDOF = IDOF + 1
            ENDDO
        ENDDO
        !WRITE (*,*)
             !     
    END SUBROUTINE Wittrick_Williams_Eigenvector
!----------------------------------------------------------------------------------------  
    
   END Module Wittrick_Williams 