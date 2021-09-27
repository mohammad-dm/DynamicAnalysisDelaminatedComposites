Module TRANSIENT
!
!IMPLICIT NONE
!
Public:: TRANSIENT_D,mass_matrix,SPECTRAL_MASS_MATRIX,GAUSSIAN_POINTS,ASSEMB_FE,ALFA,TRANSIENT_RESPONSE,INITIAL_DISP_FROM_STEADY_RESPONSE,INITIAL_VELOCITY_FROM_STEADY_RESPONSE
    contains
    !
    !*******************************************************************************************   
    subroutine mass_matrix ( ne, L, Be, En, MASS_CONSISTENT_IMPOSED,NJ, SI, NF, OM, A_,B_,D_,I_,Me,RES,RES_tr,Ca,Ca_tr)
        !defining variables 
        Double complex, dimension (:,:), intent(out) :: MASS_CONSISTENT_IMPOSED       !IN AND OUT
        Double precision, dimension(:,:,:), intent (in) :: A_,B_,D_                     !input	
        Double precision, dimension(:,:), intent (in) ::I_,RES,RES_tr,Ca,Ca_tr
        Double precision, dimension(:), intent(in) :: L                        !input
        Double precision, intent(in) :: OM
        integer, dimension(:), intent (in) :: Be, En, SI                                !input
        integer, intent (in):: NJ,NF,Me
        !
        double complex ,DIMENSION(:,:),ALLOCATABLE :: MASS_CONSISTENT_GLOBAL          !IN AND OUT
        double complex ,DIMENSION(:,:),ALLOCATABLE ::ZC,ZCS,CZC
        DOUBLE complex, dimension(6,6) :: M_SPECTRAL_ELEMENT
        double complex, dimension(6,6) :: M_CONSISTENT_ELEMENT
        double complex, dimension(6,6) :: M_LUMPED_ELEMENT
        !
        !---------------------
        ALLOCATE ( MASS_CONSISTENT_GLOBAL (3*NJ,3*NJ) )
        !WRITE (*,*)
        MASS_CONSISTENT_IMPOSED = 0.0
        MASS_CONSISTENT_GLOBAL = 0.0
        do ielem=1,ne
            CALL LUMPED_MASS_MATRIX (  L, IELEM,A_,B_,D_,I_, M_LUMPED_ELEMENT, OM)
            CALL ASSEMB_FE ( IELEM,M_LUMPED_ELEMENT,MASS_CONSISTENT_GLOBAL,BE,EN,NJ )
            !CALL SPECTRAL_MASS_MATRIX ( L, IELEM,A_,B_,D_,I_, M_SPECTRAL_ELEMENT, OM)
            !CALL ASSEMB_FE ( IELEM,M_spectral_ELEMENT,MASS_CONSISTENT_GLOBAL,BE,EN,NJ )
        ENDDO
        !********************************************14/7/2021************************************************
!
ALLOCATE(ZC(ME,3*NJ),ZCS(ME,ME),CZC(NF,ME))
!
!CALL DGEMMD( ME, 3*NJ, 3*NJ, RES_tr,MASS_CONSISTENT_GLOBAL ,ZC )  
ZC=MATMUL(RES_tr,MASS_CONSISTENT_GLOBAL)

!
!CALL DGEMMD( ME, ME, 3*NJ, ZC,RES ,ZCS ) 
ZCS=MATMUL(ZC,RES)

!
!CALL DGEMMD( NF, ME, ME, Ca_tr,ZCS ,CZC ) 
CZC=MATMUL(Ca_tr,ZCS)
!
!CALL DGEMMD( NF, NF, ME, CZC,Ca ,MASS_CONSISTENT_IMPOSED ) 
MASS_CONSISTENT_IMPOSED=MATMUL(CZC,Ca)
!
DEALLOCATE(ZC,ZCS,CZC)

!********************************************14/7/2021************************************************

        WRITE(*,*)
    end subroutine mass_matrix 
   !*******************************************************************************************    
    SUBROUTINE LUMPED_MASS_MATRIX (  L, IELEM,A_,B_,D_,I_, M_LUMPED_ELEMENT, OM)
    
        DOUBLE complex, dimension (6,6), intent (OUT) :: M_LUMPED_ELEMENT
        DOUBLE PRECISION, dimension(:,:,:), intent (in) :: A_,B_,D_	
        Double precision, dimension(:,:), intent (in) :: I_
        DOUBLE PRECISION, dimension(:), intent (in) :: L
        DOUBLE PRECISION, intent (in) :: OM
        INTEGER, intent (in) :: IELEM
        DOUBLE PRECISION MM 
        
        M_LUMPED_ELEMENT = 0.0
        
		MM = I_(1,ielem) * L(ielem)
        
		M_LUMPED_ELEMENT (1,1) = MM / 2
        M_LUMPED_ELEMENT (2,2) = MM / 2
        M_LUMPED_ELEMENT (4,4) = MM / 2
        M_LUMPED_ELEMENT (5,5) = MM / 2
		
    END SUBROUTINE LUMPED_MASS_MATRIX 
    
    SUBROUTINE SPECTRAL_MASS_MATRIX (  L, IELEM,A_,B_,D_,I_, M_SPECTRAL_ELEMENT, OM)
        !DEFINING VARIABLES
        DOUBLE complex, dimension (6,6), intent (OUT) :: M_SPECTRAL_ELEMENT
        !
        DOUBLE PRECISION, dimension(:,:,:), intent (in) :: A_,B_,D_	
        Double precision, dimension(:,:), intent (in) :: I_
        DOUBLE PRECISION, dimension(:), intent (in) :: L
        DOUBLE PRECISION, intent (in) :: OM
        INTEGER, intent (in) :: IELEM
        !
        DOUBLE PRECISION, allocatable, DIMENSION (:) :: ksiQ, WeiQ 
        DOUBLE PRECISION, DIMENSION (1,4) :: N_BEAM
        DOUBLE PRECISION, DIMENSION (4,1) :: N_BEAM_T
        DOUBLE PRECISION, DIMENSION (4,4) :: N_BEAM_T_N_BEAM
        DOUBLE PRECISION, DIMENSION (1,2) :: N_BAR
        DOUBLE PRECISION, DIMENSION (2,1) :: N_BAR_T
        DOUBLE PRECISION, DIMENSION (2,2) :: N_BAR_T_N_BAR
        DOUBLE PRECISION, DIMENSION (2,2) :: N_BAR_T_N_BAR_X
        DOUBLE PRECISION, DIMENSION (4,4) :: N_BEAM_T_N_BEAM_X
        DOUBLE PRECISION X
        DOUBLE PRECISION KL
        DOUBLE PRECISION KF, L_BAR
        INTEGER NUM_GAUSSPOINTS
        INTEGER NUM_INTERVALS
        !------------------
        M_SPECTRAL_ELEMENT = 0.0
        NUM_GAUSSPOINTS = 31
        ALLOCATE ( ksiQ( NUM_GAUSSPOINTS ), WeiQ ( NUM_GAUSSPOINTS ) )
        CALL GAUSSIAN_POINTS ( NUM_GAUSSPOINTS , ksiQ , WeiQ )
        !---------------------BAR-------------------------------------------------------
        KL = OM * SQRT ( I_ (1,IELEM) / ( A_(1,1,IELEM)  ) )
                !GAUSSIAN QUDRATURE

        !Direct integration
        N_BAR = 0.0
        N_BAR_T = 0.0
        N_BAR_T_N_BAR = 0.0
        N_BAR_T_N_BAR_X = 0.0
        iteration = 10000
        do i=1,iteration
            delta_x =  l(ielem)/iteration
            x=i*l(ielem)/iteration
            N_BAR (1,1) =  ( 1 / SIN ( KL * L(IELEM) ) ) * SIN (KL * ( L(IELEM) - X) ) 
            N_BAR (1,2) =  ( 1 / SIN ( KL * L(IELEM) ) ) * SIN ( KL * X ) 
            N_BAR_T = TRANSPOSE (N_BAR)
            N_BAR_T_N_BAR_X = MATMUL (N_BAR_T,N_BAR)
            M_SPECTRAL_ELEMENT (1,1) =  M_SPECTRAL_ELEMENT (1,1) + N_BAR_T_N_BAR_X(1,1) * delta_x 
            M_SPECTRAL_ELEMENT (1,4) =  M_SPECTRAL_ELEMENT (1,2) + N_BAR_T_N_BAR_X(1,2) * delta_x 
            M_SPECTRAL_ELEMENT (4,1) =  M_SPECTRAL_ELEMENT (2,1) + N_BAR_T_N_BAR_X(2,1) * delta_x 
            M_SPECTRAL_ELEMENT (4,4) =  M_SPECTRAL_ELEMENT (2,2) + N_BAR_T_N_BAR_X(2,2) * delta_x 
        enddo
        M_SPECTRAL_ELEMENT (1,1) = M_SPECTRAL_ELEMENT (1,1) * I_(1,ielem)
        M_SPECTRAL_ELEMENT (1,4) = M_SPECTRAL_ELEMENT (1,4) * I_(1,ielem)
        M_SPECTRAL_ELEMENT (4,1) = M_SPECTRAL_ELEMENT (4,1) * I_(1,ielem)
        M_SPECTRAL_ELEMENT (4,4) = M_SPECTRAL_ELEMENT (4,4) * I_(1,ielem)
        !---------------------------------------BEAM-------------------
        KF = SQRT (OM) * ( ( I_(1,IELEM) ) / ( D_(1,1,IELEM) ) )**0.25
        L_BAR = KF * L (IELEM)
        ETA = 2 * KF * ( 1 - COS(L_BAR)*COSH(L_BAR) )

                        
        !Direct integration
        N_BEAM = 0.0
        N_BEAM_T = 0.0
        N_BEAM_T_N_BEAM = 0.0
        N_BEAM_T_N_BEAM_X = 0.0
         do i=1,iteration
            delta_x =  l(ielem)/iteration
            x=i*l(ielem)/iteration
            X_BAR = KF * X
            N_BEAM (1,1) =  (1/ETA) * KF * ( COS(X_BAR) - COS(L_BAR-X_BAR)*COSH(L_BAR) - COS(L_BAR)*COSH(L_BAR-X_BAR) + COSH(X_BAR) + SIN(L_BAR-X_BAR)*SINH(L_BAR) - SIN(L_BAR)*SINH(L_BAR-X_BAR) )
            N_BEAM (1,2) =  (1/ETA) * ( -1*COSH(L_BAR-X_BAR)*SIN(L_BAR) +  COSH(L_BAR)*SIN(L_BAR-X_BAR) + SIN(X_BAR) - COSH(L_BAR-X_BAR)*SINH(L_BAR) +COS(L_BAR)*SINH(L_BAR-X_BAR) + SINH(X_BAR)  )
            N_BEAM (1,3) = (1/ETA) * KF * ( COS(L_BAR-X_BAR) - COS(X_BAR)*COSH(L_BAR) - COS(L_BAR)*COSH(X_BAR) + COSH(L_BAR-X_BAR) + SIN(X_BAR)*SINH(L_BAR) - SIN(L_BAR)*SINH(X_BAR)  )
            N_BEAM (1,4) = -1 * (1/ETA) * (-1*COSH(X_BAR)*SIN(L_BAR) + COSH(L_BAR)*SIN(X_BAR) + SIN(L_BAR-X_BAR) - COS(X_BAR)*SINH(L_BAR) + COS(L_BAR)*SINH(X_BAR) + SINH(L_BAR-X_BAR) )
            N_BEAM_T = TRANSPOSE (N_BEAM)
            N_BEAM_T_N_BEAM_X = MATMUL (N_BEAM_T,N_BEAM)
            M_SPECTRAL_ELEMENT (2,2) =  M_SPECTRAL_ELEMENT (2,1) + N_beam_T_N_beam_X(1,1) * delta_x 
            M_SPECTRAL_ELEMENT (2,3) =  M_SPECTRAL_ELEMENT (2,3) + N_beam_T_N_beam_X(1,2) * delta_x 
            M_SPECTRAL_ELEMENT (2,5) =  M_SPECTRAL_ELEMENT (2,5) + N_beam_T_N_beam_X(1,3) * delta_x 
            M_SPECTRAL_ELEMENT (2,6) =  M_SPECTRAL_ELEMENT (2,6) + N_beam_T_N_beam_X(1,4) * delta_x 
            
            M_SPECTRAL_ELEMENT (3,2) =  M_SPECTRAL_ELEMENT (3,2) + N_beam_T_N_beam_X(1,1) * delta_x 
            M_SPECTRAL_ELEMENT (3,3) =  M_SPECTRAL_ELEMENT (3,3) + N_beam_T_N_beam_X(1,2) * delta_x 
            M_SPECTRAL_ELEMENT (3,5) =  M_SPECTRAL_ELEMENT (3,5) + N_beam_T_N_beam_X(1,3) * delta_x 
            M_SPECTRAL_ELEMENT (3,6) =  M_SPECTRAL_ELEMENT (3,6) + N_beam_T_N_beam_X(1,4) * delta_x 
            
            M_SPECTRAL_ELEMENT (5,2) =  M_SPECTRAL_ELEMENT (5,2) + N_beam_T_N_beam_X(1,1) * delta_x 
            M_SPECTRAL_ELEMENT (5,3) =  M_SPECTRAL_ELEMENT (5,3) + N_beam_T_N_beam_X(1,2) * delta_x 
            M_SPECTRAL_ELEMENT (5,5) =  M_SPECTRAL_ELEMENT (5,5) + N_beam_T_N_beam_X(1,3) * delta_x 
            M_SPECTRAL_ELEMENT (5,6) =  M_SPECTRAL_ELEMENT (5,6) + N_beam_T_N_beam_X(1,4) * delta_x 
            
            M_SPECTRAL_ELEMENT (6,2) =  M_SPECTRAL_ELEMENT (6,2) + N_beam_T_N_beam_X(1,1) * delta_x 
            M_SPECTRAL_ELEMENT (6,3) =  M_SPECTRAL_ELEMENT (6,3) + N_beam_T_N_beam_X(1,2) * delta_x 
            M_SPECTRAL_ELEMENT (6,5) =  M_SPECTRAL_ELEMENT (6,5) + N_beam_T_N_beam_X(1,3) * delta_x 
            M_SPECTRAL_ELEMENT (6,6) =  M_SPECTRAL_ELEMENT (6,6) + N_beam_T_N_beam_X(1,4) * delta_x 
         enddo
         M_SPECTRAL_ELEMENT (2,2) = M_SPECTRAL_ELEMENT (2,2) * I_(1,ielem)
         M_SPECTRAL_ELEMENT (2,3) = M_SPECTRAL_ELEMENT (2,3) * I_(1,ielem)
         M_SPECTRAL_ELEMENT (2,5) = M_SPECTRAL_ELEMENT (2,5) * I_(1,ielem)
         M_SPECTRAL_ELEMENT (2,6) = M_SPECTRAL_ELEMENT (2,6) * I_(1,ielem)
         
         M_SPECTRAL_ELEMENT (3,2) = M_SPECTRAL_ELEMENT (3,2) * I_(1,ielem)
         M_SPECTRAL_ELEMENT (3,3) = M_SPECTRAL_ELEMENT (3,3) * I_(1,ielem)
         M_SPECTRAL_ELEMENT (3,5) = M_SPECTRAL_ELEMENT (3,5) * I_(1,ielem)
         M_SPECTRAL_ELEMENT (3,6) = M_SPECTRAL_ELEMENT (3,6) * I_(1,ielem) 
         
         M_SPECTRAL_ELEMENT (5,2) = M_SPECTRAL_ELEMENT (5,2) * I_(1,ielem)
         M_SPECTRAL_ELEMENT (5,3) = M_SPECTRAL_ELEMENT (5,3) * I_(1,ielem)
         M_SPECTRAL_ELEMENT (5,5) = M_SPECTRAL_ELEMENT (5,5) * I_(1,ielem)
         M_SPECTRAL_ELEMENT (5,6) = M_SPECTRAL_ELEMENT (5,6) * I_(1,ielem)
         
         M_SPECTRAL_ELEMENT (6,2) = M_SPECTRAL_ELEMENT (6,2) * I_(1,ielem)
         M_SPECTRAL_ELEMENT (6,3) = M_SPECTRAL_ELEMENT (6,3) * I_(1,ielem)
         M_SPECTRAL_ELEMENT (6,5) = M_SPECTRAL_ELEMENT (6,5) * I_(1,ielem)
         M_SPECTRAL_ELEMENT (6,6) = M_SPECTRAL_ELEMENT (6,6) * I_(1,ielem)
         
    END SUBROUTINE SPECTRAL_MASS_MATRIX 
   !*******************************************************************************************    
    SUBROUTINE GAUSSIAN_POINTS ( PIE , ksiQ , WeiQ)
    !
        DOUBLE PRECISION, DIMENSION (:),INTENT (OUT) :: ksiQ, WeiQ 
        INTEGER ,INTENT (IN) :: PIE
        !----------------------------
        !
        if (pie==1) then

             ksiQ(1)=0 
            !---------------------------
             WeiQ(1)=2.0

        else if (pie==2) then

             ksiQ(1)=-0.57735026919           
             ksiQ(2)= 0.57735026919
            !---------------------------
              WeiQ(1)=1.0
              WeiQ(2)=1.0

        else if (pie==3) then   

             ksiQ(1)=-0.77459666924         
             ksiQ(2)=0
             ksiQ(3)= 0.77459666924 
            !---------------------------
              WeiQ(1)=5.0/9.0
              WeiQ(2)=8.0/9.0
              WeiQ(3)=5.0/9.0

        else if (pie==4) then   

             ksiQ(1)=-0.8611363116
             ksiQ(2)=-0.339981436
             ksiQ(3)= 0.339981436
             ksiQ(4)= 0.8611363116
            !---------------------------
              WeiQ(1)=0.3478548451
              WeiQ(2)=0.6521451549
              WeiQ(3)=0.6521451549
              WeiQ(4)=0.3478548451

        else  if (pie==5) then

             ksiQ(1)=-0.906179845938664
             ksiQ(2)=-0.538469310105683
             ksiQ(3)= 0.0
             ksiQ(4)= 0.538469310105683
             ksiQ(5)= 0.906179845938664
            !---------------------------
              WeiQ(1)=0.236926885056189
              WeiQ(2)=0.478628670499366
              WeiQ(3)=0.568888888888889
              WeiQ(4)=0.478628670499366
              WeiQ(5)=0.236926885056189

        else  if (pie==6) then

             ksiQ(1)=-0.9324695142031521
             ksiQ(2)=-0.6612093864662645
             ksiQ(3)=-0.2386191860831969
             ksiQ(4)= 0.2386191860831969
             ksiQ(5)= 0.6612093864662645
             ksiQ(6)= 0.9324695142031521
            !---------------------------
              WeiQ(1)=0.1713244923791704
              WeiQ(2)=0.3607615730481386
              WeiQ(3)=0.4679139345726910
              WeiQ(4)=0.4679139345726910
              WeiQ(5)=0.3607615730481386
              WeiQ(6)=0.1713244923791704

        else if (pie==7) then   

             ksiQ(1)=-0.9491079123427585
             ksiQ(2)=-0.7415311855993945
             ksiQ(3)=-0.4058451513773972
             ksiQ(4)=0
             ksiQ(5)=0.4058451513773972
             ksiQ(6)=0.7415311855993945
             ksiQ(7)=0.9491079123427585
            !---------------------------
              WeiQ(1)=0.1294849661688697
              WeiQ(2)=0.2797053914892766
              WeiQ(3)=0.3818300505051189
              WeiQ(4)=0.4179591836734694 
              WeiQ(5)=0.3818300505051189 
              WeiQ(6)=0.2797053914892766 
              WeiQ(7)=0.1294849661688697

        else if (pie==8) then   

             ksiQ(1)=-0.9602898564975363
             ksiQ(2)=-0.7966664774136267
             ksiQ(3)=-0.5255324099163290
             ksiQ(4)=-0.1834346424956498
             ksiQ(5)=0.1834346424956498
             ksiQ(6)=0.5255324099163290
             ksiQ(7)=0.7966664774136267
             ksiQ(8)=0.9602898564975363
            !---------------------------
              WeiQ(1)=0.1012285362903763
              WeiQ(2)=0.2223810344533745
              WeiQ(3)=0.3137066458778873
              WeiQ(4)=0.3626837833783620
              WeiQ(5)=0.3626837833783620
              WeiQ(6)=0.3137066458778873
              WeiQ(7)=0.2223810344533745
              WeiQ(8)=0.1012285362903763

          else  if (pie==9) then

             ksiQ(1)=-0.9681602395076261
             ksiQ(2)=-0.8360311073266358
             ksiQ(3)=-0.6133714327005904
             ksiQ(4)=-0.3242534234038089
             ksiQ(5)=0
             ksiQ(6)=0.3242534234038089
             ksiQ(7)=0.6133714327005904
             ksiQ(8)=0.8360311073266358
             ksiQ(9)=0.9681602395076261
            !---------------------------       
             WeiQ(1)=0.0812743883615744
             WeiQ(2)=0.1806481606948574
             WeiQ(3)=0.2606106964029354
             WeiQ(4)=0.3123470770400029
             WeiQ(5)=0.3302393550012598
             WeiQ(6)=0.3123470770400029
             WeiQ(7)=0.2606106964029354
             WeiQ(8)=0.1806481606948574
             WeiQ(9)=0.0812743883615744

         else  if (pie==10) then

             ksiQ(1)=-0.9739065285171717
             ksiQ(2)=-0.8650633666889845
             ksiQ(3)=-0.6794095682990244
             ksiQ(4)=-0.4333953941292472
             ksiQ(5)=-0.1488743389816312
             ksiQ(6)=0.1488743389816312
             ksiQ(7)=0.4333953941292472
             ksiQ(8)=0.6794095682990244
             ksiQ(9)=0.8650633666889845
             ksiQ(10)=0.9739065285171717
           !---------------------------
             WeiQ(1)=0.0666713443086881
             WeiQ(2)=0.1494513491505806
             WeiQ(3)=0.2190863625159820
             WeiQ(4)=0.2692667193099963
             WeiQ(5)=0.2955242247147529
             WeiQ(6)=0.2955242247147529
             WeiQ(7)=0.2692667193099963
             WeiQ(8)=0.2190863625159820
             WeiQ(9)=0.1494513491505806
             WeiQ(10)=0.0666713443086881
             
         else if (pie == 31) then
             
               ksiQ(1) = 0.0
               ksiQ(2) =  -0.099555312
               ksiQ(3) =0.099555312
               ksiQ(4) = -0.198121199
               ksiQ(5) = 0.198121199
               ksiQ(6) = -0.29471807
               ksiQ(7) = 0.29471807
               ksiQ(8) = -0.388385902
               ksiQ(9) = 0.388385902
               ksiQ(10) = -0.478193782
               ksiQ(11) = 0.478193782
               ksiQ(12) = -0.563249161
               ksiQ(13) = 0.563249161
               ksiQ(14) = -0.642706723
               ksiQ(15) = 0.642706723
               ksiQ(16) = -0.715776785
               ksiQ(17) = 0.715776785
               ksiQ(18) = -0.781733148
               ksiQ(19) = 0.781733148
               ksiQ(20) = -0.83992032
               ksiQ(21) = 0.83992032
               ksiQ(22) = -0.88976003
               ksiQ(23) = 0.88976003
               ksiQ(24) = -0.930756998
               ksiQ(25) = 0.930756998
               ksiQ(26) = -0.962503925
               ksiQ(27) = 0.962503925
               ksiQ(28) = -0.98468591
               ksiQ(29) = 0.98468591
               ksiQ(30) = -0.997087482
               ksiQ(31) = 0.997087482 
              !---------------------------
               WeiQ(1) = 0.099720545
               WeiQ(2) =0.099225011
               WeiQ(3) =0.099225011
               WeiQ(4) =0.097743335
               WeiQ(5) =0.097743335
               WeiQ(6) =0.095290243
               WeiQ(7) =0.095290243
               WeiQ(8) =0.091890114
               WeiQ(9) =0.091890114
               WeiQ(10) =0.087576741
               WeiQ(11) =0.087576741
               WeiQ(12) =0.082392992
               WeiQ(13) =0.082392992
               WeiQ(14) =0.076390387
               WeiQ(15) =0.076390387
               WeiQ(16) =0.069628583
               WeiQ(17) =0.069628583
               WeiQ(18) =0.062174787
               WeiQ(19) =0.062174787
               WeiQ(20) =0.054103082
               WeiQ(21) =0.054103082
               WeiQ(22) =0.045493708
               WeiQ(23) =0.045493708
               WeiQ(24) =0.036432274
               WeiQ(25) =0.036432274
               WeiQ(26) =0.027009019
               WeiQ(27) =0.027009019
               WeiQ(28) =0.017318621
               WeiQ(29) =0.017318621
               WeiQ(30) =0.007470832
               WeiQ(31) =0.007470832
        end if
        
    END SUBROUTINE GAUSSIAN_POINTS
  !*******************************************************************************************  
    SUBROUTINE ASSEMB_FE ( IELEM, ELEMENT_MATRIX, GLOBAL_MATRIX,BE,EN,NJ )
        !DEFINING VARIABLES---------------------------
        DOUBLE complex  ,DIMENSION(:,:), intent (in)  :: ELEMENT_MATRIX
        DOUBLE complex  ,DIMENSION(:,:), intent (out) :: GLOBAL_MATRIX
        INTEGER, dimension(:),  intent (in) :: Be, En
        INTEGER, intent (in) :: IELEM,NJ
        !global matrix should be changed to a allocatable matrix
        INTEGER NDOF,INODE,NODEI,IDOFN,NROWS,NROWE,JNODE,NODEJ,NCOLS,NCOLE,NDOFN,NNODE
        !---------------------------------------------
        NDOFN = 3
        NNODE = 2
        !ALLOCATE ( GLOBAL_MATRIX(3*NJ,3*NJ) )
        DO  INODE=1,NNODE
            IF ( INODE == 1 ) THEN
                NODEI = BE (IELEM)
            ENDIF
            IF ( INODE == 2 ) THEN
                NODEI = EN (IELEM)
            ENDIF
            !NODEI=LNODS(IELEM,INODE)                                          
            DO  IDOFN=1,NDOFN                                               
                NROWS=(NODEI-1)*NDOFN+IDOFN                                     
                NROWE=(INODE-1)*NDOFN+IDOFN                                     
!                  ASLOD(NROWS)=ASLOD(NROWS)+ELOAD(IELEM,NROWE)                    
!                                                                    
!     ASSEMBLE THE ELEMENT STIFFNESS MATRICES
!                                                                       
                DO  JNODE =1,NNODE
                    IF ( JNODE == 1 ) THEN
                        NODEJ = BE (IELEM)
                    ENDIF
                    IF ( JNODE == 2 ) THEN
                        NODEJ = EN (IELEM)
                    ENDIF 
                    !NODEJ=LNODS(IELEM,JNODE)                                          
                    DO  JDOFN =1,NDOFN                                              
                        NCOLS=(NODEJ-1)*NDOFN+JDOFN                                     
                        NCOLE=(JNODE-1)*NDOFN+JDOFN                                     
                        GLOBAL_MATRIX(NROWS,NCOLS)=GLOBAL_MATRIX(NROWS,NCOLS)+ELEMENT_MATRIX(NROWE,NCOLE) 
                          
                    END DO
                END DO
            END DO
        END DO 
    END SUBROUTINE ASSEMB_FE
  !*******************************************************************************************  
 
    SUBROUTINE ALFA (V, X_0, Y_0, X_DOT_0, Y_DOT_0, NF, ALFA_0, ALFA_DOT_0, NUM_FREQ, NE,Be, En,L, MASS , NJ, SI,A_,B_,D_,I_, LAMD,Me,RES,RES_tr,Ca,Ca_tr )
        !DEFINING VARIABLES
        DOUBLE complex, dimension (:,:),  intent(out) :: MASS
        !
        DOUBLE PRECISION, dimension(:,:,:), intent (in) :: A_,B_,D_ !input	
        DOUBLE PRECISION, DIMENSION (:,:),  intent (in) :: V,I_,RES,RES_tr,Ca,Ca_tr
        DOUBLE PRECISION, dimension (:), intent(in)  :: L !input
        DOUBLE PRECISION, DIMENSION (:), intent(in)  :: Y_0, Y_DOT_0
        DOUBLE PRECISION, DIMENSION (:), intent(out) :: X_0, X_DOT_0
        DOUBLE PRECISION, DIMENSION (:), intent(out) :: ALFA_0, ALFA_DOT_0
        INTEGER, dimension(:), intent (in) :: Be, En !input
        INTEGER, intent(in) :: NF, NE,Me
       !
        DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: PHI, PHI_T
        DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: DIFF,DIFF_DOT
        DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: PHI_T_M
        DOUBLE PRECISION, DIMENSION (1,1) :: B,C,D
        integer R,NUM_FREQ
        !
        !-----FOR SPECTRAL MASS
        DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE :: LAMD 
        DOUBLE PRECISION OM
        integer, dimension(:), allocatable:: SI
        !-----------------
        R=1
        ALLOCATE ( PHI(NF,R), PHI_T(R,NF) )
        ALLOCATE ( DIFF(NF,R), DIFF_DOT(NF,R) )
        ALLOCATE ( PHI_T_M (R,NF) )
        ALFA_0 = 0.0
        ALFA_DOT_0 = 0.0
        B = 0.0
        C=0.0
        D=0.0
        PHI = 0.0
        DIFF = 0.0
        DIFF_DOT = 0.0
        X_DOT_0 = 0.0
        x_0 = 0.0
        DO I=1,NUM_FREQ
            OM = LAMD (I)
            CALL mass_matrix ( ne, L, Be, En, MASS, NJ, SI, NF,  OM, A_,B_,D_,I_,Me,RES,RES_tr,Ca,Ca_tr)
            DO J=1,NF
                PHI(J,R) = V ( J,I ) 
                DIFF (J,R) = X_0 (J) - Y_0 (J)
                DIFF_DOT (J,R) = X_DOT_0 (J) - Y_DOT_0 (J)
            ENDDO
            PHI_T = TRANSPOSE (PHI)    
            PHI_T_M = MATMUL ( PHI_T, MASS )   
            B  =  MATMUL ( PHI_T_M , PHI )
            C  =  MATMUL ( PHI_T_M , DIFF )
            D  =  MATMUL ( PHI_T_M , DIFF_DOT )
            ALFA_0 (I) = C(1,1) / B(1,1)
            ALFA_DOT_0 (I) = D(1,1) / B(1,1)
        ENDDO
    END SUBROUTINE ALFA  
   !*******************************************************************************************     
    SUBROUTINE TRANSIENT_RESPONSE (ALFA_0, ALFA_DOT_0, V, NF, LAMD, N, SF, NUM_FREQ, USER_DOF )
        !DEFINING VARIABLES---------------------
        DOUBLE PRECISION,DIMENSION(:,:), INTENT(IN) :: V
        DOUBLE PRECISION, DIMENSION (:), INTENT(IN) :: ALFA_0, ALFA_DOT_0,LAMD
        DOUBLE PRECISION, INTENT(IN):: SF
        INTEGER, INTENT(IN):: N,NF,NUM_FREQ,USER_DOF
        !
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: OMEGA
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PHI
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TR
        DOUBLE PRECISION, DIMENSION (1) :: ALFA_T
        DOUBLE PRECISION T, INPUT
        !---------------------------------------
        ALLOCATE ( OMEGA (NUM_FREQ) , PHI (NF), TR(NF) )
        OPEN(36,FILE="TRANSIENT_ALLDOFS_ALLTIMESTEPS.DAC", FORM="FORMATTED",STATUS="UNKNOWN")
        OPEN(37,FILE="TRANSIENT_RESPONSE_OF_SELECTED_DOF.DAC", FORM="FORMATTED",STATUS="UNKNOWN")
        !OPEN(41,FILE="ALFA_T.DAC", FORM="FORMATTED",STATUS="UNKNOWN")
        !
        DO I=1,NUM_FREQ
            OMEGA(I) = LAMD(I)
        ENDDO
        !WRITE (*,*)
        TS = 1/SF
        DO I = 1, N !EACH TIME STEP
            TR = 0.0
            T = (I-1)*TS
            DO J=1,NUM_FREQ !SUMMING UP THE CONTRIBUTION OF ALL MODES 
                DO K=1,NF
                    PHI (K) = V (K,J)
                ENDDO
                ALFA_T (1) = ALFA_0 (J) * COS ( OMEGA(J)*T ) + ALFA_DOT_0 (J) * SIN ( OMEGA(J) * T  ) / OMEGA (J)
                DO K=1,NF
                    TR(K) = TR(K) + PHI(K) * ALFA_T(1)
                ENDDO
            ENDDO
            !WRITE (36,*) T , "--------"
            DO J=1,NF 
                WRITE (36,*) TR (J)
            ENDDO
        ENDDO
        REWIND 36
        DO I=1,NF*N
            READ (36,*) INPUT 
            IF ( MOD (I,NF) == USER_DOF )THEN
                WRITE (37,*)INPUT
            ENDIF
        ENDDO
        WRITE (*,*)"The transient response is stored in TRANSIENT RESPONSE OF SELECTED DOF.DAC"
        
    
    END SUBROUTINE TRANSIENT_RESPONSE
    !*******************************************************************************************    
!
        SUBROUTINE INITIAL_DISP_FROM_STEADY_RESPONSE (Y_0, N, NF,USER_DOF)
        !DEFINING VARIABLES----------
        double precision, DIMENSION (:), intent (out) :: Y_0
        INTEGER, intent (in) :: N, NF,USER_DOF
        !
        double precision INPUT 
        !----------------------------
        OPEN(35,FILE="STEADY_STATE_ALL_DOFS_ALL_FREQUENCIES.DAC", FORM="FORMATTED",STATUS="UNKNOWN")
        REWIND 35
        K=1
        DO I=1,NF*N
            READ (35,*) INPUT
            IF ( MOD(I,N) == 1 )THEN
                Y_0(K) = INPUT
                K = K + 1
            ENDIF
        ENDDO
        !

    END SUBROUTINE INITIAL_DISP_FROM_STEADY_RESPONSE
    !
    SUBROUTINE INITIAL_VELOCITY_FROM_STEADY_RESPONSE ( Y_0, Y_DOT_0 , N, NF, SF)
        !DEFINING VARIABLES----------
        DOUBLE PRECISION, DIMENSION (:), intent (in) :: Y_0
        DOUBLE PRECISION, DIMENSION (:), intent (out) ::Y_DOT_0 
        DOUBLE PRECISION, intent (in) :: sf
        INTEGER, intent (in) :: N, NF
        !
        double precision, DIMENSION (:), ALLOCATABLE :: Y_1
        Double precision ts,INPUT
        !---------------------------
       ALLOCATE ( Y_1 (NF) )
       TS = 1/SF
       REWIND 35
        K=1
        DO I=1,NF*N
            READ  (35,*) INPUT
            IF ( MOD(I,N) == 2 )THEN
                Y_1(K) = INPUT
                K = K + 1
            ENDIF
        ENDDO
    
        DO I=1,NF 
             Y_DOT_0(I) = ( Y_1(I) - Y_0(I) ) / TS
        ENDDO
        
    END SUBROUTINE INITIAL_VELOCITY_FROM_STEADY_RESPONSE
    !
END Module TRANSIENT 