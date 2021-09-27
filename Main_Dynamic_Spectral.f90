PROGRAM DYNAMIC_SPECTRAL_ELEMENT 
!
Use Spectral_Stiffness
Use Sampling
Use RAN
Use Steady_State
Use Transient
Use Wittrick_Williams
Use COORD
USE SUPPO
USE MAS_SLAVE
!
!---------------Overall Picture----------------
!This Program performs Spectral Euler-Bernoulli Beam-Bar analysis
!-----A given time-domain signal is converted to Frequency Domain
!       +FFT module is called
!           The sampling time needs to be entered. The sampling frequency is the inverse
!           Length of the time signal needs to be entered. The Length of the frequency signal is the nearest power of it.
!           When the time signal is real only. The frequency signal is symmetric around Nyquist and zero. So only half of the frequency length (content) is used.
!               The sampling time determines the frequency range of the transformed signal
!               The Length of signal (which is the period in time domain) determines the spacing of the discrete frequencies of the transformed signal 
!-----Important to note that the FFT provides the steady state response as it is the integral up to infinity. However, we want the integral up to the time t
!----The integral from infinity to t is the transient part (see Boyd).
!---Even if there is no initial displacement or velocity imposed there is transient response due to excitation only. 
!--The Transient response is composed of natural frequency content and the Steady state response is composed of excitation frequency content.
!-In order to determine the Transient contribution modal analysis is used. The amplitude of each natural mode is determined by using mode participation of initial conditions (see alfa in Veletsos)
!-Natural Frequencies are calculated using the Wittrick-Willams algorithm as the Spectral stiffness matrix is trancendental
!-An adaptive mesh refinement approach is implemented in the Wittrick-Willams algorithm to calculate the clamped=clamped frequencies
!----------------------------------------------
!
Implicit none
    !
    Double precision OM
    double complex , DIMENSION(:,:), ALLOCATABLE :: MASS
    integer:: i,j,k,ik,IELEM, NUM_FREQ
    !----------------------------------
    !Defining Sampling_freq variables
    Double complex ,DIMENSION(:),allocatable::X
    Double precision  sf
    integer N,lj,NS,ND1 
    INTEGER,DIMENSION(:),allocatable:: B_type
    !----------------------------------
    !Defining STRUCT_PROPERTIES variables
    integer NJ,NE,NF ,Me
    Double precision, dimension(:,:,:),allocatable:: A_,B_,D_ 
    Double precision, dimension(:,:), allocatable::I_ ,RES,RES_tr,Ca,Ca_tr,TA
    Double precision, dimension(:), allocatable:: Coor,L
    integer, dimension(:), allocatable:: Be, En, SI
    !-----------------------------------------
    !Defining ELEMENT_SPECTRAL_STIFFNESS variables
    !------------------------------
    !------------------------------
    !DEFINING EIGSOLVE VARIABLES
    double precision,DIMENSION(:,:), ALLOCATABLE :: V
    double precision,DIMENSION(:), ALLOCATABLE :: LAMD  
    integer :: ierr
    !------------------------------
    !DEFINING VELETSOS METHOD PARAMETERS
    Double precision, DIMENSION (:), ALLOCATABLE :: X_0, X_DOT_0, Y_0, Y_DOT_0
    Double precision, DIMENSION (:), ALLOCATABLE :: ALFA_0, ALFA_DOT_0
    Double precision::A0(3,3),A1(3,3),A01(3,3),P_0(3,3),P_0T(3,3),A_0(3,3)
    INTEGER USER_DOF,USED_DOF
    !    
    !
    CALL Coordinates (NJ,COOR,Ne,L,BE,EN,B_type)
    ALLOCATE(SI(3*NJ))
    CALL Supports (NJ,NS,SI)
    ME=3*NJ-NS
    ALLOCATE(RES(3*NJ,ME),RES_tr(ME,3*NJ),TA(ME,ME))
    CALL MINDG(RES,RES_tr,SI,ME,NJ)
    CALL CONSTRAINT(3*NJ,ME,NF,TA,SI)
    !    NF 		!Final degree of freedom after boundary conditions and MPCs are imposed
    ALLOCATE(Ca_tr(NF,ME),Ca(Me,NF) )
    Ca=0
    Ca_tr=0
    DO I=1,ME
    DO J=1,NF
    Ca_tr(J,I)=TA(I,J)
    Ca(I,J)=TA(I,J)
    END DO
    END DO
    !
    !
    allocate(  A_(3,3,NE),B_(3,3,NE),D_(3,3,NE),I_(3,NE) )
    !
!    
A_=0
B_=0
D_=0
I_=0
!
    DO IELEM = 1,NE
    !** Element Properties are collected. Later they are called for Spectral Element Calculations at Different Stages
    If      (B_type(IELEM)==1.or.B_type(IELEM)==2) then                                                                           ! Isotropic beam
    CALL PROPERTIES           (IELEM,A_,B_,D_,I_)    ! Information about Element Sectional Properties 
    Else If (B_type(IELEM)==3.or.B_type(IELEM)==4) then                                                                           ! Composite beam
    CALL COMPOSITE_PROPERTIES (IELEM,A_,B_,D_,I_)    ! Information about Element Sectional Properties 
    Else If (B_type(IELEM)==5.or.B_type(IELEM)==6) then                                                                           ! Interface spring
    CALL INTERFACE_PROPERTIES (IELEM,A_)             ! Information about Element Sectional Properties 
    End If
    !
    End Do
    !
    !NJ: Number of Nodes
    !NE: Number of Elements
    !NS: Number of Restraints
    !NF: Number of Free-Degrees-of-Freedom
    !SI: Support information
    !BE: First node of an element    
!------------------------------------------------------------   
!!!!!!!!!!!!!BEGINNING OF STEADY STATE ANALYSIS!!!!!!!!!!!!!!!
!------------------------------------------------------------   
    !** Collects info about FFT sampling and also where to apply excitation and collect output
    CALL SAMPLING_INFO (sf,N,ND1,X,lj,NUM_FREQ,USED_DOF,NJ,SI,USER_DOF)     
    !
    !sf : Sampling frequency
    !N  : Length of the Time Signal
    !ND1: Number of Excitation Data  ! If N>ND1 then those X's are zero. The difference might be due to N being power of 2
    !lj : Loaded Degree-Of-Freedom
    !
    !** Performs Fast-Fourier-Transform  - Here it is the Forward Transform (e^-i) of the Excitation Data
    Call FFT(N,ND1,-1,X)
    !
    !
    !EN: Second node of an element 
    !
    !** At each frequency content the response is obtained by using the frequency dependent dynamic stiffness matrix and combined by using Inverse FFT
    Call Steady_State_D ( N, NF, SF, lj, NJ, NE, A_, B_, D_, I_, L, BE, EN, ND1, X, B_type, Me, RES, RES_tr, Ca, Ca_tr )
    !
    Print*,' '
    Print*,'------------------------------------'
    !    
!------------------------------------------------------------   
!!!!!!!!!!!!!BEGINNING OF TRANSIENT ANALYSIS!!!!!!!!!!!!!!!
!------------------------------------------------------------
    !** Determines Natural Frequencies and Mode Shapes
    CALL Wittrick_Williams_Eigenvalue  ( LAMD, NE, A_, B_, D_, I_, L, BE, EN, NF, SI, NJ,  NUM_FREQ, B_Type,Me,RES,RES_tr,Ca,Ca_tr )
    CALL Wittrick_Williams_Eigenvector ( V, NE, A_, B_, D_, I_, L, BE, EN, NF, SI, NJ, LAMD, NUM_FREQ, B_Type,Me,RES,RES_tr,Ca,Ca_tr )
    Print*,' '
    Print*,'------------------------------------'
    !
    ! Veletsos Method----------------------------------------------------------------
    !**Determines the Initial Deflection and Velocity based on the Steady-State Response
    ALLOCATE ( Y_0(NF), Y_DOT_0(NF), X_0(NF), X_DOT_0(NF) )
    ALLOCATE ( ALFA_0 (NF), ALFA_DOT_0 (NF) )
    CALL INITIAL_DISP_FROM_STEADY_RESPONSE ( Y_0, N, NF, USER_DOF)              ! Identifies the Initial Displacement at the Excited DOF
    CALL INITIAL_VELOCITY_FROM_STEADY_RESPONSE ( Y_0, Y_DOT_0 , N, NF, SF)      ! Identifies the Initial Velocity at the Excited DOF
    Print*,' '
    Print*,'------------------------------------'
    !
    ALLOCATE ( MASS (NF,NF) )
    !
    !**Determines the Modal Contributions using Natural Frequencies - ALFA_O for Initial Displacement and ALFA_DOT_0 for Initial Velocity
    CALL ALFA (V, X_0, Y_0, X_DOT_0, Y_DOT_0, NF, ALFA_0, ALFA_DOT_0, NUM_FREQ, NE, Be, En, L, MASS, NJ, SI, A_, B_, D_, I_, LAMD, Me, RES, RES_tr, Ca, Ca_tr )
    CALL TRANSIENT_RESPONSE (ALFA_0, ALFA_DOT_0, V, NF, LAMD, N, SF, NUM_FREQ, USER_DOF)  ! Each mode contributions are combined using that they are trigonometric time functions 
    Print*,' '
    Print*,'------------------------------------'
!
!-------------------------------------------------------------------------------
!!------------------------------------------------------------   
!!!!!!!!!!!!!!!!END OF TRANSIENT ANALYSIS!!!!!!!!!!!!!!!
!!------------------------------------------------------------   
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!
!    Print*,'Press Enter to close the screen'
!    Pause' '
!!
END PROGRAM DYNAMIC_SPECTRAL_ELEMENT
    