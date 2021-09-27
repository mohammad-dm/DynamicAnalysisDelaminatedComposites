Module Sampling
IMPLICIT NONE
!
Public:: SAMPLING_INFO
contains 
    SUBROUTINE SAMPLING_INFO(SF,N,ND1,X,lj,N_MO,USED_DOF,NJ,SI,USER_DOF) !
    !DEFINING VARIABLES-------------------
    Double precision, intent (out) :: sf
    integer , intent (out) :: N,N_MO,USED_DOF,USER_DOF,lj,ND1
    integer , intent (in) :: NJ
    integer , dimension(:),intent (in) ::SI
    Double  complex ,DIMENSION(:) ,allocatable, intent (out) ::X
    !
    Double precision om,TS
    integer n_l,fj,i,udo,a
    double precision p,t
    integer nearest_power_of_2
    !-----------------------------------------
    !--------------------
    !
    !
    OPEN(13,FILE="INPUT_CHECK.DAC",FORM="FORMATTED",STATUS="UNKNOWN")
    !OPEN(141,FILE="timesteps.DAC", FORM="FORMATTED",STATUS="UNKNOWN")
    !
    OPEN(25,FILE="FFT_INFO.TXT", FORM="FORMATTED",STATUS="UNKNOWN")
    OPEN(7,FILE="EXCITATION_DATA.TXT", FORM="FORMATTED",STATUS="UNKNOWN")
    !
    !
    Print*,'What is the sampling time? - Delta_T'
    Read(25,*),TS!
    sf=1/TS
    Print*,'The sampling frequency (Delta_f) is ',sf
    !
    Write(13,*)'The sampling time (Delta_T) is ',TS
    Write(13,*)'--------------------------------'
    Write(13,*)' '
    Write(13,*)'The sampling frequency (Delta_f) is ',sf
    Write(13,*)'--------------------------------'
    Write(13,*)' '
    !
    !
    !
    Print*,'What is the length of the time signal? - n_L '
    Read(25,*),n_L
    Write(13,*)'Length of Frequency signal is',n_L
    !
   !
    do i=0,n_L
    if (n_L <= 2**i) exit 
    end do 
    !
    nearest_power_of_2 = int(2**i)      !Next power of 2 from length of x
    N=nearest_power_of_2           !TOTAL NUMBER OF TRANSFORMED VALUES-FFT contains N many number of frequency component
    !
    Print*,'Next power of 2 from the signal length - N = ',N
    !  
    Write(13,*)'Next power of 2 from the signal length - N = ',N
    !
    Write(13,*)'--------------------------------'
    Write(13,*)' '
    !
    !    
    do i=0,N/2                     !FFT produces symmetric results and only half of it is required
    om=sf/N*i
    end do
    !        
    WRITE(13,*)'DATA'
    ND1=n_L
    WRITE(13,*)'NUMBER OF X (DATA) ND1=',ND1
    !N=NF
    WRITE(13,*)'TOTAL NUMBER OF DATA AND TRANSFORMED VALUES N=',N 
    Write(13,*)'--------------------------------'
    Write(13,*)' '
    !
    !
    !
    Read(25,*),N_MO !
    write (13,*)'Number of Modes to be determined in the Wittrick-Williams',N_MO
    Write(13,*)'--------------------------------'
    Write(13,*)' '
    
    Read(7,*),USED_DOF !
    write (13,*)'Selected Output DOF',USED_DOF
    Write(13,*)'--------------------------------'
    Write(13,*)' '

    !
    ! Creating data_x file
    !   do i=1,n_l
    !    t = (i-1) * ts
    !    write (141,*) t
    !    p = 100*dsin(3.14159265359*2*1*t) !excitation function
    !    write (7,*) p
    !   end do
    !
    Write(13,*)'-----**************************-----'
    Write(13,*)' '
    !
    !
    !
    read(7,*)lj                 ! Which degree-of-freedom is loaded
    !
    write (13,*)'The degree-of-freedom that the excitation is applied ',lj
    Write(13,*)'--------------------------------'
    Write(13,*)' '
    !  
    ALLOCATE(X(N+ND1))
    X=0.0
    !
    DO I=1,ND1
    READ(7,*)X(I)               ! Excitation values corresponding discrete times    
    END DO
    !
    !
    write (13,*)'Excitation values corresponding to discrete times'
    DO I=1,ND1
    write(13,*)X(I)             ! Excitation values corresponding discrete times    
    END DO
    !
    Write(13,*)'--------------------------------'
    Write(13,*)' '
    !
    !
        Print*,' '
        Print*,'------------------------------------'
    !
    udo=0
    do i=1,3*NJ
        a=a+1
            If (a==3) then
            a=0		                    ! read the fixed Dof
            End If   
    fj=3*SI(i)-3+a    
    If (fj<USED_DOF) then
    udo=udo+1
    end if
    end do
    !
    USER_DOF=USED_DOF-udo
    !
    END SUBROUTINE SAMPLING_INFO
    !
    END Module Sampling  
    !*************************************************************************


        
    