Module Steady_State
USE Spectral_Stiffness
USE RAN
!IMPLICIT NONE
!
Public:: Steady_State_D 
    contains
!    
    SUBROUTINE Steady_State_D ( N, NF, SF, lj, NJ, NE, A_,B_,D_,I_, L, BE, EN, ND1 , X, B_type,Me,RES,RES_tr,Ca,Ca_tr )
    !DEFINING VARIABLES---------------------
    !SAMPLING FREQUENCY AND LOAD VARIABLES
    !
    Double complex ,DIMENSION(:), intent (in)::X
    Double precision SF,OM
    integer, intent(in):: ND1, lj,NJ,NE,NF,Me
    integer, intent(inout):: N
   !STRUCTURE PROPERTIES INPUT
    Double precision, dimension(:,:,:), intent (in):: A_,B_,D_
    Double precision, dimension(:,:), intent (in)::I_,RES,RES_tr,Ca,Ca_tr
    Double precision, dimension(:), intent (in):: L
    integer, dimension(:), intent (in) ::  B_type,Be, En
    !STIFFNESS MATRICES 
    Double complex ,DIMENSION (:,:), ALLOCATABLE ::  K_SPECTRAL_GLOBAL,CZC,ZC,ZCS,K_G
    DOUBLE COMPLEX ,DIMENSION (6,6) :: K_SPECTRAL_ELEMENT
    DOUBLE COMPLEX, DIMENSION (:), ALLOCATABLE :: U_HAT
    !
    DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: ALL_REAL, ALL_IMAG
    DOUBLE complex, DIMENSION (:), ALLOCATABLE :: FR,FI,UI,UR,F_r ,SFr,F_i ,SFi
    Double precision, dimension(:), allocatable:: Urj,Uij
    Double precision::Pi=3.14159265359
    !---------------------------------------
    !variables for separatingthe steady state response of a selected degree of freedom from all the degrees 
    !of freedom written in file number 350 and write it into file 39 
    Double precision INPUT
    INTEGER DOF_START_INDEX,DOF_END_INDEX,i,j
    !
    !
    OPEN(350,FILE="STEADY_STATE_ALL_DOFS_ALL_FREQUENCIES.DAC", FORM="FORMATTED",STATUS="UNKNOWN")
    OPEN(39,FILE="STEADY_STATE_RESPONSE_OF_SELECTED_DOF.DAC", FORM="FORMATTED",STATUS="UNKNOWN")
    !
    !
        ALLOCATE ( FR(3*NJ), FI(3*NJ) )
        allocate ( Ur(NF), Ui(NF) )
        ALLOCATE ( ALL_REAL (NF,N), ALL_IMAG(NF,N) )
        allocate ( Urj (3*NJ), Uij (3*NJ) )
        ALLOCATE ( K_SPECTRAL_GLOBAL (3*NJ,3*NJ) )
        ALLOCATE(ZC(ME,3*NJ),ZCS(ME,ME),K_G(NF,NF),CZC(NF,ME),F_r(ME),SFr(NF),F_i(ME),SFi(NF))

        FR=0
        FI=0
        K_SPECTRAL_GLOBAL=0
        !
        DO ik=1,N       !Each frequency loop
            !
            K_SPECTRAL_GLOBAL = 0.0
            Fr(lj)=dreal(X(ik))
            Fi(lj)=dimag(X(ik))
            om=sf/N*ik*2*Pi     !This is angular frequency    !23/11/2020
            !
OPEN(22,FILE="K_matrix.txt", FORM="FORMATTED",STATUS="UNKNOWN") 
write (22,*) '  '
write (22,*) '***********************************************'  
write (22,*) 'STIFFNESS MATRIX FOR STEADY-STATE ANALYSIS'
write (22,*) '         cOLUMN-BY-cOLUMN '
write (22,*) '  '
            !
            DO IELEM = 1,NE
            !    
            If (B_type(IELEM)==1) then
                CALL EULER_BERNOULLI_ELEMENT_SPECTRAL_STIFFNESS     (A_,B_,D_,I_,L,K_SPECTRAL_ELEMENT,OM,IELEM)
write (22,*) 'EULER-BERNOULLI BEAM FOR OMEGA:',OM
            Else If (B_type(IELEM)==2) then    
                CALL TIMOSHENKO_ELEMENT_SPECTRAL_STIFFNESS          (A_,B_,D_,I_,L,K_SPECTRAL_ELEMENT,OM,IELEM)   
write (22,*) 'TIMOSHENKO BEAM FOR OMEGA:',OM
            Else If (B_type(IELEM)==3) then    
                CALL JUN_COMPOSITE_ELEMENT_SPECTRAL_STIFFNESS       (A_,B_,D_,I_,L,K_SPECTRAL_ELEMENT,OM,IELEM)
write (22,*) 'JUN COMPOSITE BEAM FOR OMEGA:',OM
            Else If (B_type(IELEM)==4) then                    
                CALL GOPALA_COMPOSITE_ELEMENT_SPECTRAL_STIFFNESS    (A_,B_,D_,I_,L,K_SPECTRAL_ELEMENT,OM,IELEM)            
write (22,*) 'GOPALAKRISHNAN COMPOSITE BEAM FOR OMEGA:',OM
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
    write (22,*) '-------------------STEADY---------------------'
    write (22,*) '-------------------------------------------------'
    write (22,*) ' '               
            !           
                CALL ASSEMB ( IELEM,K_SPECTRAL_ELEMENT, K_SPECTRAL_GLOBAL,BE,EN ) 
            !    
            ENDDO   !IELEM
            !
write (22,*) 'oooooooooooooooooooooooooooooooooooooooooooooooo'  
!
!********************************************14/7/2021************************************************
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
!CALL DGEMVD (ME,3*NJ , RES_tr, Fr,F_r)
F_R=MATMUL(RES_tr,FR)
!
!CALL DGEMVD( NF, ME, Ca_tr,F_r ,SFr )  
SFR=MATMUL(Ca_tr,F_R)
!
!CALL DGEMVD (ME,3*NJ , RES_tr, Fi,F_i)
F_i=MATMUL(RES_tr,Fi)
!
!CALL DGEMVD( NF, ME, Ca_tr,F_i ,SFi )  
SFi=MATMUL(Ca_tr,F_i)

!!
DEALLOCATE(ZC,ZCS,CZC,F_r,F_i)

!********************************************14/7/2021************************************************
!
CALL SOLVER(SFR,SFI,K_G,NF,Ur,Ui)
!
           !
            !STORING FOR GETTING THE TRANSIENT LATER
            DO K=1,NF       !Number of Free degrees-of-freedom
                ALL_REAL(K,IK)= UR(K)
                ALL_IMAG(K,IK)= UI(K)
            END DO
            ! 
        END DO      !Frequency Loop ik
        
     Do i=1,NE
       
    If (B_type(i)==1) Then
        write (39,*)'Beam type was selected as Euler-Bernoulli'
    Else If (B_type(i)==2) Then
        write (39,*)'Beam type was selected as Timoshenko'
    Else If (B_type(i)==3) Then
        write (39,*)'Beam type was selected as Jun et al Composite'
    Else If (B_type(i)==4) Then
         write (39,*)'Beam type was selected as Gopalakrishnan et al Composite'
    Else If (B_type(i)==5) then                    
         write (39,*) 'Interface element was used'      
    End If
        !
        !
        If ( B_ (1,1,i) ==0 .and. I_ (2,i) ==0 ) then             !----------------UnCoupled-Axial-Bending------------
        !
        print*,'**********************Warning for Steady State Analysis********************'
        print*,'Element ',i,' !!!Switched to Timoshenko Stiffness Matrix due to Symmetry!!!'
        print*,'***********************************************************'
        !
        write(39,*),'********************Warning********************************'
        write(39,*),'Element ',i,' !!!Switched to Timoshenko Stiffness Matrix due to Symmetry!!!'
        write(39,*),'***********************************************************'
        !
        Else
                !
        print*,'    '
        print*,'Element ',i,' ---calculated as NON-symmetric composite---'
        print*,'    '
        !
        write(39,*)'    '
        write(39,*),'Element ',i,'---calculated as NON-symmetric composite---'
        write(39,*)'    '
        !

        End If
        !
    END Do

        
   ! 
   !////////////////////////////////////////////Bringing the Steady response back to time domain\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        Nyquist= N/2+1
        ALLOCATE ( U_hat(N+ND1) ) 
        U_HAT = 0.0
        DO I=1,NF
        !    
            DO IK=1,Nyquist
                U_HAT (IK) = dcmplx ( ALL_REAL (I,IK), ALL_IMAG (I,IK) )  ! inverse Fourier Transform for U_hat
            END DO 
        !
            DO IK=1,Nyquist
                U_HAT (IK-1+Nyquist) = dcmplx ( ALL_REAL (I,Nyquist+1-IK), -1*ALL_IMAG (I,Nyquist+1-IK) )  ! U_hat includes complex conjugat beyond Nyquist because the signal was real
            END DO 
        !
      ! I have obtained the solution by using the complex conjugate of U_hat after Nyquist.
      ! This makes sense as real signals in time are symmetric around origin and fNyquist.       
      ! Reconstruction of U_HAT page 37 Doyle's book Make sure U_HAT real only consider Nyquist and Complex conjugate. Lecture 14 Notes
            Call FFT(N,ND1,1,U_hat)
        !
            DO IK=1,N
            WRITE (350,*) real( U_HAT(IK) )/N    
            END DO
      !
        END DO
      !
        
                !WRITE (*,*)"PLEASE ENTER DEGREE OF FREEDOM NUMBER YOU WOULD LIKE TO SEE THE DISPLACEMENT OF"
        !WRITE (*,*)"PLEASE NOTE THAT THE NUMBER IS FOR THE FREE DOFS NOT THE GLOBAL ONES "
        !READ(*,*)USER_DOF
        REWIND 350
        !
        DO I=1,NF*N
            READ (350,*) INPUT
            DOF_START_INDEX = ( USER_DOF - 1 ) * N + 1
            DOF_END_INDEX =  ( USER_DOF ) * N 
            IF ( I >=  DOF_START_INDEX .AND. I <= DOF_END_INDEX )THEN
                WRITE (39,*) INPUT
            ENDIF
        ENDDO
        
       WRITE (*,*)"The steady state response is stored in STEADYSTATE RESPONSE OF SELECTED DOF.DAC"
!       
!
      write (22,*) ' '  
      write (22,*) '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||'  
      write (22,*) ' |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| '  
      write (22,*) '  ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||  '  
      write (22,*) '   ||||||||||||||||||||||||||||||||||||||||||||||||||||||||   '  
      write (22,*) '    ||||||||||||||||||||||||||||||||||||||||||||||||||||||    '  
      write (22,*) '     ||||||||||||||||||||||||||||||||||||||||||||||||||||     ' 
      write (22,*) '      ||||||||||||||||||||||||||||||||||||||||||||||||||      '  
      write (22,*) '       ||||||||||||||||||||||||||||||||||||||||||||||||       '  
      write (22,*) '        ||||||||||||||||||||||||||||||||||||||||||||||        '  
      write (22,*) '         ||||||||||||||||||||||||||||||||||||||||||||         '  
      write (22,*) '          ||||||||||||||||||||||||||||||||||||||||||          '  
      write (22,*) '           ||||||||||||||||||||||||||||||||||||||||           '  
      write (22,*) '            ||||||||||||||||||||||||||||||||||||||            '  
      write (22,*) '             ||||||||||||||||||||||||||||||||||||             '  
      write (22,*) '              ||||||||||||||||||||||||||||||||||              '  
      write (22,*) '               ||||||||||||||||||||||||||||||||               '  
      write (22,*) '                ||||||||||||||||||||||||||||||                '  
      write (22,*) '                 ||||||||||||||||||||||||||||                 ' 
      write (22,*) '                  ||||||||||||||||||||||||||                  '  
      write (22,*) '                   ||||||||||||||||||||||||                   '  
      write (22,*) '                    ||||||||||||||||||||||                    '  
      write (22,*) '                     ||||||||||||||||||||                     '  
      write (22,*) '                      ||||||||||||||||||                      '  
      write (22,*) '                       ||||||||||||||||                       '  
      write (22,*) '                        ||||||||||||||                        '  
      write (22,*) '                         ||||||||||||                         '  
      write (22,*) '                          ||||||||||                          '  
      write (22,*) '                           ||||||||                           '  
      write (22,*) '                            ||||||                            '  
      write (22,*) '                             ||||                             '  
      write (22,*) '                              ||                              '  
     write (22,*) ' '  
    END SUBROUTINE Steady_State_D
    !
    !
   END Module Steady_State 