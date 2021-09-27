	  MODULE CUSH
	  USE CARP,ONLY:DGEMMD,DGEMVD,XERBLA
	  PUBLIC DSYSVX 
	  CONTAINS
	 SUBROUTINE DSYSVX( N,NRHS,A,LDA,BAK,DP,rcond)

      !SUBROUTINE DSYSVX( FACT,UPLO,N,NRHS,A,LDA,AF,LDAF,IPIV,B,LDB,X,LDX,RCOND,FERR,BERR,WORK,LWORK,IWORK,INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          FACT, UPLO
      INTEGER            INFO, LDA, LDAF, LDB, LDX, LWORK, N, NRHS	,I,J
      DOUBLE PRECISION ,  intent (out) :: RCOND
!     ..
!     .. Array Arguments ..
    !  INTEGER            IPIV( * ), IWORK( * )
    !  DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ),BERR( * ), FERR( * ), WORK( * ), X( LDX, * )
	INTEGER,DIMENSION(:),ALLOCATABLE::IPIV,IWORK
	DOUBLE PRECISION,DIMENSION(:,:),INTENT(IN)::A
	DOUBLE PRECISION,DIMENSION(:,:),INTENT(IN):: BAK
	DOUBLE PRECISION,DIMENSION(:,:),INTENT(OUT)::DP
	DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::  AF	 
	DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::BERR, FERR  ,WORK
	DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::X,B
!     ..
!
!  Purpose
!  =======
!
!  DSYSVX uses the diagonal pivoting factorization to compute the
!  solution to a real system of linear equations A * X = B,
!  where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
!  matrices.
!
!  Error bounds on the solution and a condition estimate are also
!  provided.
!
!  Description
!  ===========
!
!  The following steps are performed:
!
!  1. If FACT = 'N', the diagonal pivoting method is used to factor A.
!     The form of the factorization is
!        A = U * D * U**T,  if UPLO = 'U', or
!        A = L * D * L**T,  if UPLO = 'L',		  ! BKZ BATHE SAYFA 440
!     where U (or L) is a product of permutation and unit upper (lower)
!     triangular matrices, and D is symmetric and block diagonal with
!     1-by-1 and 2-by-2 diagonal blocks.
!
!  2. If some D(i,i)=0, so that D is exactly singular, then the routine
!     returns with INFO = i. Otherwise, the factored form of A is used
!     to estimate the condition number of the matrix A.  If the
!     reciprocal of the condition number is less than machine precision,
!     INFO = N+1 is returned as a warning, but the routine still goes on
!     to solve for X and compute error bounds as described below.
!
!  3. The system of equations is solved for X using the factored form
!     of A.
!
!  4. Iterative refinement is applied to improve the computed solution
!     matrix and calculate error bounds and backward error estimates
!     for it.
!
!  Arguments
!  =========
!
!  FACT    (input) CHARACTER*1
!          Specifies whether or not the factored form of A has been
!          supplied on entry.
!          = 'F':  On entry, AF and IPIV contain the factored form of
!                  A.  AF and IPIV will not be modified.
!          = 'N':  The matrix A will be copied to AF and factored.
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The number of linear equations, i.e., the order of the
!          matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrices B and X.  NRHS >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The symmetric matrix A.  If UPLO = 'U', the leading N-by-N
!          upper triangular part of A contains the upper triangular part
!          of the matrix A, and the strictly lower triangular part of A
!          is not referenced.  If UPLO = 'L', the leading N-by-N lower
!          triangular part of A contains the lower triangular part of
!          the matrix A, and the strictly upper triangular part of A is
!          not referenced.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  AF      (input or output) DOUBLE PRECISION array, dimension (LDAF,N)
!          If FACT = 'F', then AF is an input argument and on entry
!          contains the block diagonal matrix D and the multipliers used
!          to obtain the factor U or L from the factorization
!          A = U*D*U**T or A = L*D*L**T as computed by DSYTRF.
!
!          If FACT = 'N', then AF is an output argument and on exit
!          returns the block diagonal matrix D and the multipliers used
!          to obtain the factor U or L from the factorization
!          A = U*D*U**T or A = L*D*L**T.
!
!  LDAF    (input) INTEGER
!          The leading dimension of the array AF.  LDAF >= max(1,N).
!
!  IPIV    (input or output) INTEGER array, dimension (N)
!          If FACT = 'F', then IPIV is an input argument and on entry
!          contains details of the interchanges and the block structure
!          of D, as determined by DSYTRF.
!          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!          interchanged and D(k,k) is a 1-by-1 diagonal block.
!          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!
!          If FACT = 'N', then IPIV is an output argument and on exit
!          contains details of the interchanges and the block structure
!          of D, as determined by DSYTRF.
!
!  B       (input) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          The N-by-NRHS right hand side matrix B.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  X       (output) DOUBLE PRECISION array, dimension (LDX,NRHS)
!          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X.
!
!  LDX     (input) INTEGER
!          The leading dimension of the array X.  LDX >= max(1,N).
!
!  RCOND   (output) DOUBLE PRECISION
!          The estimate of the reciprocal condition number of the matrix
!          A.  If RCOND is less than the machine precision (in
!          particular, if RCOND = 0), the matrix is singular to working
!          precision.  This condition is indicated by a return code of
!          INFO > 0.
!
!  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)
!          The estimated forward error bound for each solution vector
!          X(j) (the j-th column of the solution matrix X).
!          If XTRUE is the true solution corresponding to X(j), FERR(j)
!          is an estimated upper bound for the magnitude of the largest
!          element in (X(j) - XTRUE) divided by the magnitude of the
!          largest element in X(j).  The estimate is as reliable as
!          the estimate for RCOND, and is almost always a slight
!          overestimate of the true error.
!
!  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)
!          The componentwise relative backward error of each solution
!          vector X(j) (i.e., the smallest relative change in
!          any element of A or B that makes X(j) an exact solution).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The length of WORK.  LWORK >= 3*N, and for best performance
!          LWORK >= N*NB, where NB is the optimal blocksize for
!          DSYTRF.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  IWORK   (workspace) INTEGER array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, and i is
!                <= N:  D(i,i) is exactly zero.  The factorization
!                       has been completed but the factor D is exactly
!                       singular, so the solution and error bounds could
!                       not be computed. RCOND = 0 is returned.
!                = N+1: D is nonsingular, but RCOND is less than machine
!                       precision, meaning that the matrix is singular
!                       to working precision.  Nevertheless, the
!                       solution and error bounds are computed because
!                       there are a number of situations where the
!                       computed solution can be more accurate than the
!                       value of RCOND would suggest.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, NOFACT
      INTEGER            LWKOPT, NB
      DOUBLE PRECISION   ANORM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANSY
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLACPY, DSYCON, DSYRFS, DSYTRF, DSYTRS, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.



!	********************************************

integer sf1,sf2,sf3

FACT="N"
UPLO="L"
LDAF=N
LDB=N
LDX=N
LWORK=65*N !N*100 nb=64 tahmin edildi

ALLOCATE(IPIV(N),IWORK(N),STAT=sf1)
ALLOCATE(AF(LDAF,N),BERR(NRHS),FERR(NRHS),WORK(LWORK), STAT=sf2)!
ALLOCATE(X(LDX,NRHS),B(LDB,NRHS), STAT=sf3)
 !Print*,'sf',sf1,sf2,sf3
IF(sf1.or.sf2.or.sf3/=0) THEN
  STOP
  End if
!PRINT*,' allocation is done'





DO J=1,NRHS
DO I=1,N
B(I,J)=BAK(I,J)
END DO
End do


!**********************************************
    


      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.NOFACT .AND. .NOT.LSAME( FACT, 'F' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) )THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDAF.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -13
      ELSE IF( LWORK.LT.MAX( 1, 3*N ) .AND. .NOT.LQUERY ) THEN
         INFO = -18
      END IF
!
      IF( INFO.EQ.0 ) THEN
         NB = ILAENV( 1, 'DSYTRF', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
	   WORK( 1 ) = LWKOPT
       END IF
!
      IF( INFO.NE.0 ) THEN
		 CALL XERBLA( 'DSYSVX', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
      IF( NOFACT ) THEN
!
!        Compute the factorization A = U*D*U' or A = L*D*L'.
!
         CALL DLACPY( UPLO, N, N, A, LDA, AF, LDAF )
         CALL DSYTRF( UPLO, N, AF, LDAF, IPIV, WORK, LWORK, INFO )
!
!        Return if INFO is non-zero.
!
         IF( INFO.NE.0 ) THEN
            IF( INFO.GT.0 )RCOND = ZERO
			 DEALLOCATE(IPIV,IWORK,AF,BERR,FERR,WORK,X,B)
			 print*,'deallocated at line 310'
            RETURN
         END IF
      END IF
!
!     Compute the norm of the matrix A.
!
      ANORM = DLANSY( 'I', UPLO, N, A, LDA, WORK )
!
!     Compute the reciprocal of the condition number of A.
!
      CALL DSYCON( UPLO, N, AF, LDAF, IPIV, ANORM, RCOND, WORK, IWORK,INFO )
!
!     Set INFO = N+1 if the matrix is singular to working precision.
!
      IF( RCOND.LT.DLAMCH( 'Epsilon' ) )INFO = N + 1
!
!     Compute the solution vectors X.
!
      CALL DLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL DSYTRS( UPLO, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO )
!
!     Use iterative refinement to improve the computed solutions and
!     compute error bounds and backward error estimates for them.
!
      CALL DSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )
!
      
!
!     End of DSYSVX
!
!************
DO J=1,NRHS
DO I=1,N
DP(I,J)=X(I,J)
END DO
End do
!***********

!Print*,'EqSolver DSYSVX   INFO=',INFO,'       NB=',NB ,'       UPLO =',UPLO
!Print*,'RCOND=',RCOND,'                FERR=',FERR,'                     BERR=',BERR


 
  DEALLOCATE(IPIV,IWORK,AF,BERR,FERR,WORK,X,B,stat=sf1)
	IF(sf1/=0) THEN
PRINT*,'BELLEK TAHSISINDE HATA VAR-DSYSVX2!'
STOP
END IF
!Print*,'Deallocation is complete for n=',n ,sf1
  
   RETURN
       END SUBROUTINE DSYSVX



  END MODULE CUSH
	  !*********************************************************

	  !*************************************************
	        SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DLACPY copies all or part of a two-dimensional matrix A to another
!  matrix B.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be copied to B.
!          = 'U':      Upper triangular part
!          = 'L':      Lower triangular part
!          Otherwise:  All of the matrix A
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The m by n matrix A.  If UPLO = 'U', only the upper triangle
!          or trapezoid is accessed; if UPLO = 'L', only the lower
!          triangle or trapezoid is accessed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  B       (output) DOUBLE PRECISION array, dimension (LDB,N)
!          On exit, B = A in the locations specified by UPLO.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,M).
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
      RETURN
!
!     End of DLACPY
!
      END SUBROUTINE DLACPY


	  !************************************************
         SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DSYTRF computes the factorization of a real symmetric matrix A using
!  the Bunch-Kaufman diagonal pivoting method.  The form of the
!  factorization is
!
!     A = U*D*U**T  or  A = L*D*L**T
!
!  where U (or L) is a product of permutation and unit upper (lower)
!  triangular matrices, and D is symmetric and block diagonal with
!  1-by-1 and 2-by-2 diagonal blocks.
!
!  This is the blocked version of the algorithm, calling Level 3 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, the block diagonal matrix D and the multipliers used
!          to obtain the factor U or L (see below for further details).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (output) INTEGER array, dimension (N)
!          Details of the interchanges and the block structure of D.
!          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!          interchanged and D(k,k) is a 1-by-1 diagonal block.
!          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The length of WORK.  LWORK >=1.  For best performance
!          LWORK >= N*NB, where NB is the block size returned by ILAENV.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
!                has been completed, but the block diagonal matrix D is
!                exactly singular, and division by zero will occur if it
!                is used to solve a system of equations.
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', then A = U*D*U', where
!     U = P(n)*U(n)* ... *P(k)U(k)* ...,
!  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
!  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!
!             (   I    v    0   )   k-s
!     U(k) =  (   0    I    0   )   s
!             (   0    0    I   )   n-k
!                k-s   s   n-k
!
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
!  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
!  and A(k,k), and v overwrites A(1:k-2,k-1:k).
!
!  If UPLO = 'L', then A = L*D*L', where
!     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
!  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
!  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!
!             (   I    0     0   )  k-1
!     L(k) =  (   0    I     0   )  s
!             (   0    v     I   )  n-k-s+1
!                k-1   s  n-k-s+1
!
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
!  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
!  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            IINFO, IWS, J, K, KB, LDWORK, LWKOPT, NB, NBMIN
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASYF, DSYTF2, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Determine the block size
!
         NB = ILAENV( 1, 'DSYTRF', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = LDWORK*NB
         IF( LWORK.LT.IWS ) THEN
            NB = MAX( LWORK / LDWORK, 1 )
            NBMIN = MAX( 2, ILAENV( 2, 'DSYTRF', UPLO, N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = 1
      END IF
      IF( NB.LT.NBMIN )NB = N
!
      IF( UPPER ) THEN
!
!        Factorize A as U*D*U' using the upper triangle of A
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        KB, where KB is the number of columns factorized by DLASYF;
!        KB is either NB or NB-1, or K for the last block
!
         K = N
   10    CONTINUE
!
!        If K < 1, exit from loop
!
         IF( K.LT.1 ) GOTO 40
!
         IF( K.GT.NB ) THEN
!
!           Factorize columns k-kb+1:k of A and use blocked code to
!           update columns 1:k-kb
!
            CALL DLASYF( UPLO, K, NB, KB, A, LDA, IPIV, WORK, LDWORK,IINFO )
         ELSE
!
!           Use unblocked code to factorize columns 1:k of A
!
            CALL DSYTF2( UPLO, K, A, LDA, IPIV, IINFO )
            KB = K
         END IF
!
!        Set INFO on the first occurrence of a zero pivot
!
         IF( INFO.EQ.0 .AND. IINFO.GT.0 )INFO = IINFO
!
!        Decrease K and return to the start of the main loop
!
         K = K - KB
         goto 10
!
      ELSE
!
!        Factorize A as L*D*L' using the lower triangle of A
!
!        K is the main loop index, increasing from 1 to N in steps of
!        KB, where KB is the number of columns factorized by DLASYF;
!        KB is either NB or NB-1, or N-K+1 for the last block
!
         K = 1
   20    CONTINUE
!
!        If K > N, exit from loop
!
         IF( K.GT.N )GOTO 40
!
         IF( K.LE.N-NB ) THEN
!
!           Factorize columns k:k+kb-1 of A and use blocked code to
!           update columns k+kb:n
!
            CALL DLASYF( UPLO, N-K+1, NB, KB, A( K, K ), LDA, IPIV( K ), WORK, LDWORK, IINFO )
         ELSE
!
!           Use unblocked code to factorize columns k:n of A
!
            CALL DSYTF2( UPLO, N-K+1, A( K, K ), LDA, IPIV( K ), IINFO )
            KB = N - K + 1
         END IF
!
!        Set INFO on the first occurrence of a zero pivot
!
         IF( INFO.EQ.0 .AND. IINFO.GT.0 )INFO = IINFO + K - 1
!
!        Adjust IPIV
!
         DO 30 J = K, K + KB - 1
            IF( IPIV( J ).GT.0 ) THEN
               IPIV( J ) = IPIV( J ) + K - 1
            ELSE
               IPIV( J ) = IPIV( J ) - K + 1
            END IF
   30    CONTINUE
!
!        Increase K and return to the start of the main loop
!
         K = K + KB
         GOTO 20
!
      END IF
!
   40 CONTINUE
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of DSYTRF
!
      END SUBROUTINE DSYTRF
 
! ****************************************
      SUBROUTINE DSYCON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK,IWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   ANORM, RCOND
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * ), IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DSYCON estimates the reciprocal of the condition number (in the
!  1-norm) of a real symmetric matrix A using the factorization
!  A = U*D*U**T or A = L*D*L**T computed by DSYTRF.
!
!  An estimate is obtained for norm(inv(A)), and the reciprocal of the
!  condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the details of the factorization are stored
!          as an upper or lower triangular matrix.
!          = 'U':  Upper triangular, form is A = U*D*U**T;
!          = 'L':  Lower triangular, form is A = L*D*L**T.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The block diagonal matrix D and the multipliers used to
!          obtain the factor U or L as computed by DSYTRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          Details of the interchanges and the block structure of D
!          as determined by DSYTRF.
!
!  ANORM   (input) DOUBLE PRECISION
!          The 1-norm of the original matrix A.
!
!  RCOND   (output) DOUBLE PRECISION
!          The reciprocal of the condition number of the matrix A,
!          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
!          estimate of the 1-norm of inv(A) computed in this routine.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)
!
!  IWORK    (workspace) INTEGER array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, KASE
      DOUBLE PRECISION   AINVNM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLACON, DSYTRS, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYCON', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.LE.ZERO ) THEN
         RETURN
      END IF
!
!     Check that the diagonal matrix D is nonsingular.
!
      IF( UPPER ) THEN
!
!        Upper triangular storage: examine D from bottom to top
!
         DO 10 I = N, 1, -1
            IF( IPIV( I ).GT.0 .AND. A( I, I ).EQ.ZERO ) RETURN
   10    CONTINUE
      ELSE
!
!        Lower triangular storage: examine D from top to bottom.
!
         DO 20 I = 1, N
            IF( IPIV( I ).GT.0 .AND. A( I, I ).EQ.ZERO ) RETURN
   20    CONTINUE
      END IF
!
!     Estimate the 1-norm of the inverse.
!
      KASE = 0
   30 CONTINUE
      CALL DLACON( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE )
      IF( KASE.NE.0 ) THEN
!
!        Multiply by inv(L*D*L') or inv(U*D*U').
!
         CALL DSYTRS( UPLO, N, 1, A, LDA, IPIV, WORK, N, INFO )
         goto 30
      END IF
!
!     Compute the estimate of the reciprocal condition number.
!
      IF( AINVNM.NE.ZERO )RCOND = ( ONE / AINVNM ) / ANORM
!
      RETURN
!
!     End of DSYCON
!
      END SUBROUTINE DSYCON
!********************************************************
	        SUBROUTINE DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DSYTRS solves a system of linear equations A*X = B with a real
!  symmetric matrix A using the factorization A = U*D*U**T or
!  A = L*D*L**T computed by DSYTRF.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the details of the factorization are stored
!          as an upper or lower triangular matrix.
!          = 'U':  Upper triangular, form is A = U*D*U**T;
!          = 'L':  Lower triangular, form is A = L*D*L**T.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The block diagonal matrix D and the multipliers used to
!          obtain the factor U or L as computed by DSYTRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          Details of the interchanges and the block structure of D
!          as determined by DSYTRF.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, K, KP
      DOUBLE PRECISION   AK, AKM1, AKM1K, BK, BKM1, DENOM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 )RETURN
!
      IF( UPPER ) THEN
!
!        Solve A*X = B, where A = U*D*U'.
!
!        First solve U*D*X = B, overwriting B with X.
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = N
   10    CONTINUE
!
!        If K < 1, exit from loop.
!
         IF( K.LT.1 )GOTO 30
!
         IF( IPIV( K ).GT.0 ) THEN
!
!           1 x 1 diagonal block
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            IF( KP.NE.K )CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in column K of A.
!
            CALL DGER( K-1, NRHS, -ONE, A( 1, K ), 1, B( K, 1 ), LDB,B( 1, 1 ), LDB )
!
!           Multiply by the inverse of the diagonal block.
!
            CALL DSCAL( NRHS, ONE / A( K, K ), B( K, 1 ), LDB )
            K = K - 1
         ELSE
!
!           2 x 2 diagonal block
!
!           Interchange rows K-1 and -IPIV(K).
!
            KP = -IPIV( K )
            IF( KP.NE.K-1 )CALL DSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in columns K-1 and K of A.
!
            CALL DGER( K-2, NRHS, -ONE, A( 1, K ), 1, B( K, 1 ), LDB,B( 1, 1 ), LDB )
            CALL DGER( K-2, NRHS, -ONE, A( 1, K-1 ), 1, B( K-1, 1 ),LDB, B( 1, 1 ), LDB )
!
!           Multiply by the inverse of the diagonal block.
!
            AKM1K = A( K-1, K )
            AKM1 = A( K-1, K-1 ) / AKM1K
            AK = A( K, K ) / AKM1K
            DENOM = AKM1*AK - ONE
            DO 20 J = 1, NRHS
               BKM1 = B( K-1, J ) / AKM1K
               BK = B( K, J ) / AKM1K
               B( K-1, J ) = ( AK*BKM1-BK ) / DENOM
               B( K, J ) = ( AKM1*BK-BKM1 ) / DENOM
   20       CONTINUE
            K = K - 2
         END IF
!
         goto 10
   30    CONTINUE
!
!        Next solve U'*X = B, overwriting B with X.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = 1
   40    CONTINUE
!
!        If K > N, exit from loop.
!
         IF( K.GT.N )GOTO 50
!
         IF( IPIV( K ).GT.0 ) THEN
!
!           1 x 1 diagonal block
!
!           Multiply by inv(U'(K)), where U(K) is the transformation
!           stored in column K of A.
!
            CALL DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB, A( 1, K ),1, ONE, B( K, 1 ), LDB )
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            IF( KP.NE.K )CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K + 1
         ELSE
!
!           2 x 2 diagonal block
!
!           Multiply by inv(U'(K+1)), where U(K+1) is the transformation
!           stored in columns K and K+1 of A.
!
            CALL DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB, A( 1, K ),1, ONE, B( K, 1 ), LDB )
            CALL DGEMV( 'Transpose', K-1, NRHS, -ONE, B, LDB,A( 1, K+1 ), 1, ONE, B( K+1, 1 ), LDB )
!
!           Interchange rows K and -IPIV(K).
!
            KP = -IPIV( K )
            IF( KP.NE.K )CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K + 2
         END IF
!
         goto 40
   50    CONTINUE
!
      ELSE
!
!        Solve A*X = B, where A = L*D*L'.
!
!        First solve L*D*X = B, overwriting B with X.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = 1
   60    CONTINUE
!
!        If K > N, exit from loop.
!
         IF( K.GT.N )GOTO 80
!
         IF( IPIV( K ).GT.0 ) THEN
!
!           1 x 1 diagonal block
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            IF( KP.NE.K ) CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in column K of A.
!
            IF( K.LT.N )CALL DGER( N-K, NRHS, -ONE, A( K+1, K ), 1, B( K, 1 ),LDB, B( K+1, 1 ), LDB )
!
!           Multiply by the inverse of the diagonal block.
!
            CALL DSCAL( NRHS, ONE / A( K, K ), B( K, 1 ), LDB )
            K = K + 1
         ELSE
!
!           2 x 2 diagonal block
!
!           Interchange rows K+1 and -IPIV(K).
!
            KP = -IPIV( K )
            IF( KP.NE.K+1 )CALL DSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in columns K and K+1 of A.
!
            IF( K.LT.N-1 ) THEN
               CALL DGER( N-K-1, NRHS, -ONE, A( K+2, K ), 1, B( K, 1 ),LDB, B( K+2, 1 ), LDB )
               CALL DGER( N-K-1, NRHS, -ONE, A( K+2, K+1 ), 1,B( K+1, 1 ), LDB, B( K+2, 1 ), LDB )
            END IF
!
!           Multiply by the inverse of the diagonal block.
!
            AKM1K = A( K+1, K )
            AKM1 = A( K, K ) / AKM1K
            AK = A( K+1, K+1 ) / AKM1K
            DENOM = AKM1*AK - ONE
            DO 70 J = 1, NRHS
               BKM1 = B( K, J ) / AKM1K
               BK = B( K+1, J ) / AKM1K
               B( K, J ) = ( AK*BKM1-BK ) / DENOM
               B( K+1, J ) = ( AKM1*BK-BKM1 ) / DENOM
   70       CONTINUE
            K = K + 2
         END IF
!
         goto 60
   80    CONTINUE
!
!        Next solve L'*X = B, overwriting B with X.
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
         K = N
   90    CONTINUE
!
!        If K < 1, exit from loop.
!
         IF( K.LT.1 )GOTO 100
!
         IF( IPIV( K ).GT.0 ) THEN
!
!           1 x 1 diagonal block
!
!           Multiply by inv(L'(K)), where L(K) is the transformation
!           stored in column K of A.
!
            IF( K.LT.N )CALL DGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ),LDB, A( K+1, K ), 1, ONE, B( K, 1 ), LDB )
!
!           Interchange rows K and IPIV(K).
!
            KP = IPIV( K )
            IF( KP.NE.K )CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K - 1
         ELSE
!
!           2 x 2 diagonal block
!
!           Multiply by inv(L'(K-1)), where L(K-1) is the transformation
!           stored in columns K-1 and K of A.
!
            IF( K.LT.N ) THEN
               CALL DGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ),LDB, A( K+1, K ), 1, ONE, B( K, 1 ), LDB )
               CALL DGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ),LDB, A( K+1, K-1 ), 1, ONE, B( K-1, 1 ),LDB )
            END IF
!
!           Interchange rows K and -IPIV(K).
!
            KP = -IPIV( K )
            IF( KP.NE.K )CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K = K - 2
         END IF
!
         goto 90
  100    CONTINUE
      END IF
!
      RETURN
!
!     End of DSYTRS
!
      END SUBROUTINE DSYTRS
!****************************************************
      SUBROUTINE DSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB,X, LDX, FERR, BERR, WORK, IWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * ), IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ),BERR( * ), FERR( * ), WORK( * ), X( LDX, * )
!     ..
!
!  Purpose
!  =======
!
!  DSYRFS improves the computed solution to a system of linear
!  equations when the coefficient matrix is symmetric indefinite, and
!  provides error bounds and backward error estimates for the solution.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrices B and X.  NRHS >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The symmetric matrix A.  If UPLO = 'U', the leading N-by-N
!          upper triangular part of A contains the upper triangular part
!          of the matrix A, and the strictly lower triangular part of A
!          is not referenced.  If UPLO = 'L', the leading N-by-N lower
!          triangular part of A contains the lower triangular part of
!          the matrix A, and the strictly upper triangular part of A is
!          not referenced.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  AF      (input) DOUBLE PRECISION array, dimension (LDAF,N)
!          The factored form of the matrix A.  AF contains the block
!          diagonal matrix D and the multipliers used to obtain the
!          factor U or L from the factorization A = U*D*U**T or
!          A = L*D*L**T as computed by DSYTRF.
!
!  LDAF    (input) INTEGER
!          The leading dimension of the array AF.  LDAF >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          Details of the interchanges and the block structure of D
!          as determined by DSYTRF.
!
!  B       (input) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          The right hand side matrix B.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  X       (input/output) DOUBLE PRECISION array, dimension (LDX,NRHS)
!          On entry, the solution matrix X, as computed by DSYTRS.
!          On exit, the improved solution matrix X.
!
!  LDX     (input) INTEGER
!          The leading dimension of the array X.  LDX >= max(1,N).
!
!  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)
!          The estimated forward error bound for each solution vector
!          X(j) (the j-th column of the solution matrix X).
!          If XTRUE is the true solution corresponding to X(j), FERR(j)
!          is an estimated upper bound for the magnitude of the largest
!          element in (X(j) - XTRUE) divided by the magnitude of the
!          largest element in X(j).  The estimate is as reliable as
!          the estimate for RCOND, and is almost always a slight
!          overestimate of the true error.
!
!  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)
!          The componentwise relative backward error of each solution
!          vector X(j) (i.e., the smallest relative change in
!          any element of A or B that makes X(j) an exact solution).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
!
!  IWORK   (workspace) INTEGER array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Internal Parameters
!  ===================
!
!  ITMAX is the maximum number of steps of iterative refinement.
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D+0 )
      DOUBLE PRECISION   THREE
      PARAMETER          ( THREE = 3.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            COUNT, I, J, K, KASE, NZ
      DOUBLE PRECISION   EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK
!     ..
!     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DLACON, DSYMV, DSYTRS, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDAF.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYRFS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) THEN
         DO 10 J = 1, NRHS
            FERR( J ) = ZERO
            BERR( J ) = ZERO
   10    CONTINUE
         RETURN
      END IF
!
!     NZ = maximum number of nonzero elements in each row of A, plus 1
!
      NZ = N + 1
      EPS = DLAMCH( 'Epsilon' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS
!
!     Do for each right hand side
!
      DO 140 J = 1, NRHS
!
         COUNT = 1
         LSTRES = THREE
   20    CONTINUE
!
!        Loop until stopping criterion is satisfied.
!
!        Compute residual R = B - A * X
!
         CALL DCOPY( N, B( 1, J ), 1, WORK( N+1 ), 1 )
         CALL DSYMV( UPLO, N, -ONE, A, LDA, X( 1, J ), 1, ONE,WORK( N+1 ), 1 )
!
!        Compute componentwise relative backward error from formula
!
!        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )
!
!        where abs(Z) is the componentwise absolute value of the matrix
!        or vector Z.  If the i-th component of the denominator is less
!        than SAFE2, then SAFE1 is added to the i-th components of the
!        numerator and denominator before dividing.
!
         DO 30 I = 1, N
            WORK( I ) = ABS( B( I, J ) )
   30    CONTINUE
!
!        Compute abs(A)*abs(X) + abs(B).
!
         IF( UPPER ) THEN
            DO 50 K = 1, N
               S = ZERO
               XK = ABS( X( K, J ) )
               DO 40 I = 1, K - 1
                  WORK( I ) = WORK( I ) + ABS( A( I, K ) )*XK
                  S = S + ABS( A( I, K ) )*ABS( X( I, J ) )
   40          CONTINUE
               WORK( K ) = WORK( K ) + ABS( A( K, K ) )*XK + S
   50       CONTINUE
         ELSE
            DO 70 K = 1, N
               S = ZERO
               XK = ABS( X( K, J ) )
               WORK( K ) = WORK( K ) + ABS( A( K, K ) )*XK
               DO 60 I = K + 1, N
                  WORK( I ) = WORK( I ) + ABS( A( I, K ) )*XK
                  S = S + ABS( A( I, K ) )*ABS( X( I, J ) )
   60          CONTINUE
               WORK( K ) = WORK( K ) + S
   70       CONTINUE
         END IF
         S = ZERO
         DO 80 I = 1, N
            IF( WORK( I ).GT.SAFE2 ) THEN
               S = MAX( S, ABS( WORK( N+I ) ) / WORK( I ) )
            ELSE
               S = MAX( S, ( ABS( WORK( N+I ) )+SAFE1 ) /( WORK( I )+SAFE1 ) )
            END IF
   80    CONTINUE
         BERR( J ) = S
!
!        Test stopping criterion. Continue iterating if
!           1) The residual BERR(J) is larger than machine epsilon, and
!           2) BERR(J) decreased by at least a factor of 2 during the
!              last iteration, and
!           3) At most ITMAX iterations tried.
!
         IF( BERR( J ).GT.EPS .AND. TWO*BERR( J ).LE.LSTRES .AND.COUNT.LE.ITMAX ) THEN
!
!           Update solution and try again.
!
            CALL DSYTRS( UPLO, N, 1, AF, LDAF, IPIV, WORK( N+1 ), N,INFO )
            CALL DAXPY( N, ONE, WORK( N+1 ), 1, X( 1, J ), 1 )
            LSTRES = BERR( J )
            COUNT = COUNT + 1
            goto 20
         END IF
!
!        Bound error from formula
!
!        norm(X - XTRUE) / norm(X) .le. FERR =
!        norm( abs(inv(A))*
!           ( abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) / norm(X)
!
!        where
!          norm(Z) is the magnitude of the largest component of Z
!          inv(A) is the inverse of A
!          abs(Z) is the componentwise absolute value of the matrix or
!             vector Z
!          NZ is the maximum number of nonzeros in any row of A, plus 1
!          EPS is machine epsilon
!
!        The i-th component of abs(R)+NZ*EPS*(abs(A)*abs(X)+abs(B))
!        is incremented by SAFE1 if the i-th component of
!        abs(A)*abs(X) + abs(B) is less than SAFE2.
!
!        Use DLACON to estimate the infinity-norm of the matrix
!           inv(A) * diag(W),
!        where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) )))
!
         DO 90 I = 1, N
            IF( WORK( I ).GT.SAFE2 ) THEN
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I )
            ELSE
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I ) + SAFE1
            END IF
   90    CONTINUE
!
         KASE = 0
  100    CONTINUE
         CALL DLACON( N, WORK( 2*N+1 ), WORK( N+1 ), IWORK, FERR( J ),KASE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN
!
!              Multiply by diag(W)*inv(A').
!
               CALL DSYTRS( UPLO, N, 1, AF, LDAF, IPIV, WORK( N+1 ), N,INFO )
               DO 110 I = 1, N
                  WORK( N+I ) = WORK( I )*WORK( N+I )
  110          CONTINUE
            ELSE IF( KASE.EQ.2 ) THEN
!
!              Multiply by inv(A)*diag(W).
!
               DO 120 I = 1, N
                  WORK( N+I ) = WORK( I )*WORK( N+I )
  120          CONTINUE
               CALL DSYTRS( UPLO, N, 1, AF, LDAF, IPIV, WORK( N+1 ), N,INFO )
            END IF
            goto 100
         END IF
!
!        Normalize error.
!
         LSTRES = ZERO
         DO 130 I = 1, N
            LSTRES = MAX( LSTRES, ABS( X( I, J ) ) )
  130    CONTINUE
         IF( LSTRES.NE.ZERO )FERR( J ) = FERR( J ) / LSTRES
!
  140 CONTINUE
!
      RETURN
!
!     End of DSYRFS
!
      END SUBROUTINE DSYRFS

!**********************************************************************
      SUBROUTINE DLASYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KB, LDA, LDW, N, NB
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), W( LDW, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASYF computes a partial factorization of a real symmetric matrix A
!  using the Bunch-Kaufman diagonal pivoting method. The partial
!  factorization has the form:
!
!  A  =  ( I  U12 ) ( A11  0  ) (  I    0   )  if UPLO = 'U', or:
!        ( 0  U22 ) (  0   D  ) ( U12' U22' )
!
!  A  =  ( L11  0 ) (  D   0  ) ( L11' L21' )  if UPLO = 'L'
!        ( L21  I ) (  0  A22 ) (  0    I   )
!
!  where the order of D is at most NB. The actual order is returned in
!  the argument KB, and is either NB or NB-1, or N if N <= NB.
!
!  DLASYF is an auxiliary routine called by DSYTRF. It uses blocked code
!  (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
!  A22 (if UPLO = 'L').
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NB      (input) INTEGER
!          The maximum number of columns of the matrix A that should be
!          factored.  NB should be at least 2 to allow for 2-by-2 pivot
!          blocks.
!
!  KB      (output) INTEGER
!          The number of columns of A that were actually factored.
!          KB is either NB-1 or NB, or N if N <= NB.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!          On exit, A contains details of the partial factorization.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (output) INTEGER array, dimension (N)
!          Details of the interchanges and the block structure of D.
!          If UPLO = 'U', only the last KB elements of IPIV are set;
!          if UPLO = 'L', only the first KB elements are set.
!
!          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!          interchanged and D(k,k) is a 1-by-1 diagonal block.
!          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!
!  W       (workspace) DOUBLE PRECISION array, dimension (LDW,NB)
!
!  LDW     (input) INTEGER
!          The leading dimension of the array W.  LDW >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
!               has been completed, but the block diagonal matrix D is
!               exactly singular.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            IMAX, J, JB, JJ, JMAX, JP, K, KK, KKW, KP,KSTEP, KW
      DOUBLE PRECISION   ABSAKK, ALPHA, COLMAX, D11, D21, D22, R1,ROWMAX, T
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      EXTERNAL           LSAME, IDAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DGEMV, DSCAL, DSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      INFO = 0
!
!     Initialize ALPHA for use in choosing pivot block size.
!
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
!
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Factorize the trailing columns of A using the upper triangle
!        of A and working backwards, and compute the matrix W = U12*D
!        for use in updating A11
!
!        K is the main loop index, decreasing from N in steps of 1 or 2
!
!        KW is the column of W which corresponds to column K of A
!
         K = N
   10    CONTINUE
         KW = NB + K - N
!
!        Exit from loop
!
         IF( ( K.LE.N-NB+1 .AND. NB.LT.N ) .OR. K.LT.1 )GOTO 30
!
!        Copy column K of A to column KW of W and update it
!
         CALL DCOPY( K, A( 1, K ), 1, W( 1, KW ), 1 )
         IF( K.LT.N )CALL DGEMV( 'No transpose', K, N-K, -ONE, A( 1, K+1 ), LDA,W( K, KW+1 ), LDW, ONE, W( 1, KW ), 1 )
!
         KSTEP = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
         ABSAKK = ABS( W( K, KW ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!
         IF( K.GT.1 ) THEN
            IMAX = IDAMAX( K-1, W( 1, KW ), 1 )
            COLMAX = ABS( W( IMAX, KW ) )
         ELSE
            COLMAX = ZERO
         END IF
!
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
!
!           Column K is zero: set INFO and continue
!
            IF( INFO.EQ.0 )INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
               KP = K
            ELSE
!
!              Copy column IMAX to column KW-1 of W and update it
!
               CALL DCOPY( IMAX, A( 1, IMAX ), 1, W( 1, KW-1 ), 1 )
               CALL DCOPY( K-IMAX, A( IMAX, IMAX+1 ), LDA,W( IMAX+1, KW-1 ), 1 )
               IF( K.LT.N ) CALL DGEMV( 'No transpose', K, N-K, -ONE, A( 1, K+1 ),LDA, W( IMAX, KW+1 ), LDW, ONE,W( 1, KW-1 ), 1 )
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
               JMAX = IMAX + IDAMAX( K-IMAX, W( IMAX+1, KW-1 ), 1 )
               ROWMAX = ABS( W( JMAX, KW-1 ) )
               IF( IMAX.GT.1 ) THEN
                  JMAX = IDAMAX( IMAX-1, W( 1, KW-1 ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( W( JMAX, KW-1 ) ) )
               END IF
!
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
                  KP = K
               ELSE IF( ABS( W( IMAX, KW-1 ) ).GE.ALPHA*ROWMAX ) THEN
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
                  KP = IMAX
!
!                 copy column KW-1 of W to column KW
!
                  CALL DCOPY( K, W( 1, KW-1 ), 1, W( 1, KW ), 1 )
               ELSE
!
!                 interchange rows and columns K-1 and IMAX, use 2-by-2
!                 pivot block
!
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
!
            KK = K - KSTEP + 1
            KKW = NB + KK - N
!
!           Updated column KP is already stored in column KKW of W
!
            IF( KP.NE.KK ) THEN
!
!              Copy non-updated column KK to column KP
!
               A( KP, K ) = A( KK, K )
               CALL DCOPY( K-1-KP, A( KP+1, KK ), 1, A( KP, KP+1 ),LDA )
               CALL DCOPY( KP, A( 1, KK ), 1, A( 1, KP ), 1 )
!
!              Interchange rows KK and KP in last KK columns of A and W
!
               CALL DSWAP( N-KK+1, A( KK, KK ), LDA, A( KP, KK ), LDA )
               CALL DSWAP( N-KK+1, W( KK, KKW ), LDW, W( KP, KKW ),LDW )
            END IF
!
            IF( KSTEP.EQ.1 ) THEN
!
!              1-by-1 pivot block D(k): column KW of W now holds
!
!              W(k) = U(k)*D(k)
!
!              where U(k) is the k-th column of U
!
!              Store U(k) in column k of A
!
               CALL DCOPY( K, W( 1, KW ), 1, A( 1, K ), 1 )
               R1 = ONE / A( K, K )
               CALL DSCAL( K-1, R1, A( 1, K ), 1 )
            ELSE
!
!              2-by-2 pivot block D(k): columns KW and KW-1 of W now
!              hold
!
!              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
!
!              where U(k) and U(k-1) are the k-th and (k-1)-th columns
!              of U
!
               IF( K.GT.2 ) THEN
!
!                 Store U(k) and U(k-1) in columns k and k-1 of A
!
                  D21 = W( K-1, KW )
                  D11 = W( K, KW ) / D21
                  D22 = W( K-1, KW-1 ) / D21
                  T = ONE / ( D11*D22-ONE )
                  D21 = T / D21
                  DO 20 J = 1, K - 2
                     A( J, K-1 ) = D21*( D11*W( J, KW-1 )-W( J, KW ) )
                     A( J, K ) = D21*( D22*W( J, KW )-W( J, KW-1 ) )
   20             CONTINUE
               END IF
!
!              Copy D(k) to A
!
               A( K-1, K-1 ) = W( K-1, KW-1 )
               A( K-1, K ) = W( K-1, KW )
               A( K, K ) = W( K, KW )
            END IF
         END IF
!
!        Store details of the interchanges in IPIV
!
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         END IF
!
!        Decrease K and return to the start of the main loop
!
         K = K - KSTEP
         goto 10
!
   30    CONTINUE
!
!        Update the upper triangle of A11 (= A(1:k,1:k)) as
!
!        A11 := A11 - U12*D*U12' = A11 - U12*W'
!
!        computing blocks of NB columns at a time
!
         DO 50 J = ( ( K-1 ) / NB )*NB + 1, 1, -NB
            JB = MIN( NB, K-J+1 )
!
!           Update the upper triangle of the diagonal block
!
            DO 40 JJ = J, J + JB - 1
               CALL DGEMV( 'No transpose', JJ-J+1, N-K, -ONE,A( J, K+1 ), LDA, W( JJ, KW+1 ), LDW, ONE,A( J, JJ ), 1 )
   40       CONTINUE
!
!           Update the rectangular superdiagonal block
!
            CALL DGEMM( 'No transpose', 'Transpose', J-1, JB, N-K, -ONE,A( 1, K+1 ), LDA, W( J, KW+1 ), LDW, ONE,A( 1, J ), LDA )
   50    CONTINUE
!
!        Put U12 in standard form by partially undoing the interchanges
!        in columns k+1:n
!
         J = K + 1
   60    CONTINUE
         JJ = J
         JP = IPIV( J )
         IF( JP.LT.0 ) THEN
            JP = -JP
            J = J + 1
         END IF
         J = J + 1
         IF( JP.NE.JJ .AND. J.LE.N )CALL DSWAP( N-J+1, A( JP, J ), LDA, A( JJ, J ), LDA )
         IF( J.LE.N )GOTO 60
!
!        Set KB to the number of columns factorized
!
         KB = N - K
!
      ELSE
!
!        Factorize the leading columns of A using the lower triangle
!        of A and working forwards, and compute the matrix W = L21*D
!        for use in updating A22
!
!        K is the main loop index, increasing from 1 in steps of 1 or 2
!
         K = 1
   70    CONTINUE
!
!        Exit from loop
!
         IF( ( K.GE.NB .AND. NB.LT.N ) .OR. K.GT.N )GOTO 90
!
!        Copy column K of A to column K of W and update it
!
         CALL DCOPY( N-K+1, A( K, K ), 1, W( K, K ), 1 )
         CALL DGEMV( 'No transpose', N-K+1, K-1, -ONE, A( K, 1 ), LDA,W( K, 1 ), LDW, ONE, W( K, K ), 1 )
!
         KSTEP = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
         ABSAKK = ABS( W( K, K ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!
         IF( K.LT.N ) THEN
            IMAX = K + IDAMAX( N-K, W( K+1, K ), 1 )
            COLMAX = ABS( W( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
!
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
!
!           Column K is zero: set INFO and continue
!
            IF( INFO.EQ.0 )INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
               KP = K
            ELSE
!
!              Copy column IMAX to column K+1 of W and update it
!
               CALL DCOPY( IMAX-K, A( IMAX, K ), LDA, W( K, K+1 ), 1 )
               CALL DCOPY( N-IMAX+1, A( IMAX, IMAX ), 1, W( IMAX, K+1 ), 1 )
               CALL DGEMV( 'No transpose', N-K+1, K-1, -ONE, A( K, 1 ),LDA, W( IMAX, 1 ), LDW, ONE, W( K, K+1 ), 1 )
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
               JMAX = K - 1 + IDAMAX( IMAX-K, W( K, K+1 ), 1 )
               ROWMAX = ABS( W( JMAX, K+1 ) )
               IF( IMAX.LT.N ) THEN
                  JMAX = IMAX + IDAMAX( N-IMAX, W( IMAX+1, K+1 ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( W( JMAX, K+1 ) ) )
               END IF
!
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
                  KP = K
               ELSE IF( ABS( W( IMAX, K+1 ) ).GE.ALPHA*ROWMAX ) THEN
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
                  KP = IMAX
!
!                 copy column K+1 of W to column K
!
                  CALL DCOPY( N-K+1, W( K, K+1 ), 1, W( K, K ), 1 )
               ELSE
!
!                 interchange rows and columns K+1 and IMAX, use 2-by-2
!                 pivot block
!
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
!
            KK = K + KSTEP - 1
!
!           Updated column KP is already stored in column KK of W
!
            IF( KP.NE.KK ) THEN
!
!              Copy non-updated column KK to column KP
!
               A( KP, K ) = A( KK, K )
               CALL DCOPY( KP-K-1, A( K+1, KK ), 1, A( KP, K+1 ), LDA )
               CALL DCOPY( N-KP+1, A( KP, KK ), 1, A( KP, KP ), 1 )
!
!              Interchange rows KK and KP in first KK columns of A and W
!
               CALL DSWAP( KK, A( KK, 1 ), LDA, A( KP, 1 ), LDA )
               CALL DSWAP( KK, W( KK, 1 ), LDW, W( KP, 1 ), LDW )
            END IF
!
            IF( KSTEP.EQ.1 ) THEN
!
!              1-by-1 pivot block D(k): column k of W now holds
!
!              W(k) = L(k)*D(k)
!
!              where L(k) is the k-th column of L
!
!              Store L(k) in column k of A
!
               CALL DCOPY( N-K+1, W( K, K ), 1, A( K, K ), 1 )
               IF( K.LT.N ) THEN
                  R1 = ONE / A( K, K )
                  CALL DSCAL( N-K, R1, A( K+1, K ), 1 )
               END IF
            ELSE
!
!              2-by-2 pivot block D(k): columns k and k+1 of W now hold
!
!              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
!
!              where L(k) and L(k+1) are the k-th and (k+1)-th columns
!              of L
!
               IF( K.LT.N-1 ) THEN
!
!                 Store L(k) and L(k+1) in columns k and k+1 of A
!
                  D21 = W( K+1, K )
                  D11 = W( K+1, K+1 ) / D21
                  D22 = W( K, K ) / D21
                  T = ONE / ( D11*D22-ONE )
                  D21 = T / D21
                  DO 80 J = K + 2, N
                     A( J, K ) = D21*( D11*W( J, K )-W( J, K+1 ) )
                     A( J, K+1 ) = D21*( D22*W( J, K+1 )-W( J, K ) )
   80             CONTINUE
               END IF
!
!              Copy D(k) to A
!
               A( K, K ) = W( K, K )
               A( K+1, K ) = W( K+1, K )
               A( K+1, K+1 ) = W( K+1, K+1 )
            END IF
         END IF
!
!        Store details of the interchanges in IPIV
!
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         END IF
!
!        Increase K and return to the start of the main loop
!
         K = K + KSTEP
         goto 70
!
   90    CONTINUE
!
!        Update the lower triangle of A22 (= A(k:n,k:n)) as
!
!        A22 := A22 - L21*D*L21' = A22 - L21*W'
!
!        computing blocks of NB columns at a time
!
         DO 110 J = K, N, NB
            JB = MIN( NB, N-J+1 )
!
!           Update the lower triangle of the diagonal block
!
            DO 100 JJ = J, J + JB - 1
               CALL DGEMV( 'No transpose', J+JB-JJ, K-1, -ONE,A( JJ, 1 ), LDA, W( JJ, 1 ), LDW, ONE,A( JJ, JJ ), 1 )
  100       CONTINUE
!
!           Update the rectangular subdiagonal block
!
            IF( J+JB.LE.N ) CALL DGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB,K-1, -ONE, A( J+JB, 1 ), LDA, W( J, 1 ), LDW,ONE, A( J+JB, J ), LDA )
  110    CONTINUE
!
!        Put L21 in standard form by partially undoing the interchanges
!        in columns 1:k-1
!
         J = K - 1
  120    CONTINUE
         JJ = J
         JP = IPIV( J )
         IF( JP.LT.0 ) THEN
            JP = -JP
            J = J - 1
         END IF
         J = J - 1
         IF( JP.NE.JJ .AND. J.GE.1 )CALL DSWAP( J, A( JP, 1 ), LDA, A( JJ, 1 ), LDA )
         IF( J.GE.1 )GOTO 120
!
!        Set KB to the number of columns factorized
!
         KB = K - 1
!
      END IF
      RETURN
!
!     End of DLASYF
!
      END SUBROUTINE DLASYF
!**************************************************************
	      SUBROUTINE DSYTF2( UPLO, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DSYTF2 computes the factorization of a real symmetric matrix A using
!  the Bunch-Kaufman diagonal pivoting method:
!
!     A = U*D*U'  or  A = L*D*L'
!
!  where U (or L) is a product of permutation and unit upper (lower)
!  triangular matrices, U' is the transpose of U, and D is symmetric and
!  block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, the block diagonal matrix D and the multipliers used
!          to obtain the factor U or L (see below for further details).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (output) INTEGER array, dimension (N)
!          Details of the interchanges and the block structure of D.
!          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!          interchanged and D(k,k) is a 1-by-1 diagonal block.
!          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
!               has been completed, but the block diagonal matrix D is
!               exactly singular, and division by zero will occur if it
!               is used to solve a system of equations.
!
!  Further Details
!  ===============
!
!  1-96 - Based on modifications by J. Lewis, Boeing Computer Services
!         Company
!
!  If UPLO = 'U', then A = U*D*U', where
!     U = P(n)*U(n)* ... *P(k)U(k)* ...,
!  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
!  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!
!             (   I    v    0   )   k-s
!     U(k) =  (   0    I    0   )   s
!             (   0    0    I   )   n-k
!                k-s   s   n-k
!
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
!  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
!  and A(k,k), and v overwrites A(1:k-2,k-1:k).
!
!  If UPLO = 'L', then A = L*D*L', where
!     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
! i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
! n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
! and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
! defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!
!              (   I    0     0   )  k-1
!      L(k) =  (   0    I     0   )  s
!             (   0    v     I   )  n-k-s+1
!                k-1   s  n-k-s+1
!
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
!  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
!  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IMAX, J, JMAX, K, KK, KP, KSTEP
      DOUBLE PRECISION   ABSAKK, ALPHA, COLMAX, D11, D12, D21, D22, R1,ROWMAX, T, WK, WKM1, WKP1
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      EXTERNAL           LSAME, IDAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSCAL, DSWAP, DSYR, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTF2', -INFO )
         RETURN
      END IF
!
!     Initialize ALPHA for use in choosing pivot block size.
!
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
!
      IF( UPPER ) THEN
!
!        Factorize A as U*D*U' using the upper triangle of A
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2
!
         K = N
   10    CONTINUE
!
!        If K < 1, exit from loop
!
         IF( K.LT.1 )GOTO 70
         KSTEP = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
         ABSAKK = ABS( A( K, K ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!
         IF( K.GT.1 ) THEN
            IMAX = IDAMAX( K-1, A( 1, K ), 1 )
            COLMAX = ABS( A( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
!
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
!
!           Column K is zero: set INFO and continue
!
            IF( INFO.EQ.0 )INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
               KP = K
            ELSE
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
               JMAX = IMAX + IDAMAX( K-IMAX, A( IMAX, IMAX+1 ), LDA )
               ROWMAX = ABS( A( IMAX, JMAX ) )
               IF( IMAX.GT.1 ) THEN
                  JMAX = IDAMAX( IMAX-1, A( 1, IMAX ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( A( JMAX, IMAX ) ) )
               END IF
!
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
                  KP = K
               ELSE IF( ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX ) THEN
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
                  KP = IMAX
               ELSE
!
!                 interchange rows and columns K-1 and IMAX, use 2-by-2
!                 pivot block
!
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
!
            KK = K - KSTEP + 1
            IF( KP.NE.KK ) THEN
!
!              Interchange rows and columns KK and KP in the leading
!              submatrix A(1:k,1:k)
!
               CALL DSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 )
               CALL DSWAP( KK-KP-1, A( KP+1, KK ), 1, A( KP, KP+1 ),LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               IF( KSTEP.EQ.2 ) THEN
                  T = A( K-1, K )
                  A( K-1, K ) = A( KP, K )
                  A( KP, K ) = T
               END IF
            END IF
!
!           Update the leading submatrix
!
            IF( KSTEP.EQ.1 ) THEN
!
!              1-by-1 pivot block D(k): column k now holds
!
!              W(k) = U(k)*D(k)
!
!              where U(k) is the k-th column of U
!
!              Perform a rank-1 update of A(1:k-1,1:k-1) as
!
!              A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)'
!
               R1 = ONE / A( K, K )
               CALL DSYR( UPLO, K-1, -R1, A( 1, K ), 1, A, LDA )
!
!              Store U(k) in column k
!
               CALL DSCAL( K-1, R1, A( 1, K ), 1 )
            ELSE
!
!              2-by-2 pivot block D(k): columns k and k-1 now hold
!
!              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
!
!              where U(k) and U(k-1) are the k-th and (k-1)-th columns
!              of U
!
!              Perform a rank-2 update of A(1:k-2,1:k-2) as
!
!              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )'
!                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )'
!
               IF( K.GT.2 ) THEN
!
                  D12 = A( K-1, K )
                  D22 = A( K-1, K-1 ) / D12
                  D11 = A( K, K ) / D12
                  T = ONE / ( D11*D22-ONE )
                  D12 = T / D12
!
                  DO 30 J = K - 2, 1, -1
                     WKM1 = D12*( D11*A( J, K-1 )-A( J, K ) )
                     WK = D12*( D22*A( J, K )-A( J, K-1 ) )
                     DO 20 I = J, 1, -1
                        A( I, J ) = A( I, J ) - A( I, K )*WK - A( I, K-1 )*WKM1
   20                CONTINUE
                     A( J, K ) = WK
                     A( J, K-1 ) = WKM1
   30             CONTINUE
!
               END IF
!
            END IF
         END IF
!
!        Store details of the interchanges in IPIV
!
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         END IF
!
!        Decrease K and return to the start of the main loop
!
         K = K - KSTEP
         goto 10
!
      ELSE
!
!        Factorize A as L*D*L' using the lower triangle of A
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2
!
         K = 1
   40    CONTINUE
!
!        If K > N, exit from loop
!
         IF( K.GT.N )GOTO 70
         KSTEP = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
         ABSAKK = ABS( A( K, K ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!
         IF( K.LT.N ) THEN
            IMAX = K + IDAMAX( N-K, A( K+1, K ), 1 )
            COLMAX = ABS( A( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
!
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
!
!           Column K is zero: set INFO and continue
!
            IF( INFO.EQ.0 )INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
               KP = K
            ELSE
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
               JMAX = K - 1 + IDAMAX( IMAX-K, A( IMAX, K ), LDA )
               ROWMAX = ABS( A( IMAX, JMAX ) )
               IF( IMAX.LT.N ) THEN
                  JMAX = IMAX + IDAMAX( N-IMAX, A( IMAX+1, IMAX ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( A( JMAX, IMAX ) ) )
               END IF
!
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
                  KP = K
               ELSE IF( ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX ) THEN
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
                  KP = IMAX
               ELSE
!
!                 interchange rows and columns K+1 and IMAX, use 2-by-2
!                 pivot block
!
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
!
            KK = K + KSTEP - 1
            IF( KP.NE.KK ) THEN
!
!              Interchange rows and columns KK and KP in the trailing
!              submatrix A(k:n,k:n)
!
               IF( KP.LT.N )CALL DSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 )
               CALL DSWAP( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ),LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               IF( KSTEP.EQ.2 ) THEN
                  T = A( K+1, K )
                  A( K+1, K ) = A( KP, K )
                  A( KP, K ) = T
               END IF
            END IF
!
!           Update the trailing submatrix
!
            IF( KSTEP.EQ.1 ) THEN
!
!              1-by-1 pivot block D(k): column k now holds
!
!              W(k) = L(k)*D(k)
!
!              where L(k) is the k-th column of L
!
               IF( K.LT.N ) THEN
!
!                 Perform a rank-1 update of A(k+1:n,k+1:n) as
!
!                 A := A - L(k)*D(k)*L(k)' = A - W(k)*(1/D(k))*W(k)'
!
                  D11 = ONE / A( K, K )
                  CALL DSYR( UPLO, N-K, -D11, A( K+1, K ), 1,A( K+1, K+1 ), LDA )
!
!                 Store L(k) in column K
!
                  CALL DSCAL( N-K, D11, A( K+1, K ), 1 )
               END IF
            ELSE
!
!              2-by-2 pivot block D(k)
!
               IF( K.LT.N-1 ) THEN
!
!                 Perform a rank-2 update of A(k+2:n,k+2:n) as
!
!                 A := A - ( (A(k) A(k+1))*D(k)**(-1) ) * (A(k) A(k+1))'
!
!                 where L(k) and L(k+1) are the k-th and (k+1)-th
!                 columns of L
!
                  D21 = A( K+1, K )
                  D11 = A( K+1, K+1 ) / D21
                  D22 = A( K, K ) / D21
                  T = ONE / ( D11*D22-ONE )
                  D21 = T / D21
!
                  DO 60 J = K + 2, N
!
                     WK = D21*( D11*A( J, K )-A( J, K+1 ) )
                     WKP1 = D21*( D22*A( J, K+1 )-A( J, K ) )
!
                     DO 50 I = J, N
                        A( I, J ) = A( I, J ) - A( I, K )*WK -A( I, K+1 )*WKP1
   50                CONTINUE
!
                     A( J, K ) = WK
                     A( J, K+1 ) = WKP1
!
   60             CONTINUE
               END IF
            END IF
         END IF
!
!        Store details of the interchanges in IPIV
!
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         END IF
!
!        Increase K and return to the start of the main loop
!
         K = K + KSTEP
         goto 40
!
      END IF
!
   70 CONTINUE
!
      RETURN
!
!     End of DSYTF2
!
      END SUBROUTINE DSYTF2
 !******************************************
       SUBROUTINE DLACON( N, V, X, ISGN, EST, KASE )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      INTEGER            KASE, N
      DOUBLE PRECISION   EST
!     ..
!     .. Array Arguments ..
      INTEGER            ISGN( * )
      DOUBLE PRECISION   V( * ), X( * )
!     ..
!
!  Purpose
!  =======
!
!  DLACON estimates the 1-norm of a square, real matrix A.
!  Reverse communication is used for evaluating matrix-vector products.
!
!  Arguments
!  =========
!
!  N      (input) INTEGER
!         The order of the matrix.  N >= 1.
!
!  V      (workspace) DOUBLE PRECISION array, dimension (N)
!         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
!         (W is not returned).
!
!  X      (input/output) DOUBLE PRECISION array, dimension (N)
!         On an intermediate return, X should be overwritten by
!               A * X,   if KASE=1,
!               A' * X,  if KASE=2,
!         and DLACON must be re-called with all the other parameters
!         unchanged.
!
!  ISGN   (workspace) INTEGER array, dimension (N)
!
!  EST    (output) DOUBLE PRECISION
!         An estimate (a lower bound) for norm(A).
!
!  KASE   (input/output) INTEGER
!         On the initial call to DLACON, KASE should be 0.
!         On an intermediate return, KASE will be 1 or 2, indicating
!         whether X should be overwritten by A * X  or A' * X.
!         On the final return from DLACON, KASE will again be 0.
!
!  Further Details
!  ======= =======
!
!  Contributed by Nick Higham, University of Manchester.
!  Originally named SONEST, dated March 16, 1988.
!
!  Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of
!  a real or complex matrix, with applications to condition estimation",
!  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, ITER, J, JLAST, JUMP
      DOUBLE PRECISION   ALTSGN, ESTOLD, TEMP
!     ..
!     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DASUM
      EXTERNAL           IDAMAX, DASUM
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, NINT, SIGN
!     ..
!     .. Save statement ..
      SAVE
!     ..
!     .. Executable Statements ..
!
      IF( KASE.EQ.0 ) THEN
         DO 10 I = 1, N
            X( I ) = ONE / DBLE( N )
   10    CONTINUE
         KASE = 1
         JUMP = 1
         RETURN
      END IF
!
      goto ( 20, 40, 70, 110, 140 )JUMP
!
!     ................ ENTRY   (JUMP = 1)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
!
   20 CONTINUE
      IF( N.EQ.1 ) THEN
         V( 1 ) = X( 1 )
         EST = ABS( V( 1 ) )
!        ... QUIT
         goto 150
      END IF
      EST = DASUM( N, X, 1 )
!
      DO 30 I = 1, N
         X( I ) = SIGN( ONE, X( I ) )
         ISGN( I ) = NINT( X( I ) )
   30 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
!
!     ................ ENTRY   (JUMP = 2)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X.
!
   40 CONTINUE
      J = IDAMAX( N, X, 1 )
      ITER = 2
!
!     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
!
   50 CONTINUE
      DO 60 I = 1, N
         X( I ) = ZERO
   60 CONTINUE
      X( J ) = ONE
      KASE = 1
      JUMP = 3
      RETURN
!
!     ................ ENTRY   (JUMP = 3)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
   70 CONTINUE
      CALL DCOPY( N, X, 1, V, 1 )
      ESTOLD = EST
      EST = DASUM( N, V, 1 )
      DO 80 I = 1, N
         IF( NINT( SIGN( ONE, X( I ) ) ).NE.ISGN( I ) )GOTO 90
   80 CONTINUE
!     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GOTO 120
!
   90 CONTINUE
!     TEST FOR CYCLING.
      IF( EST.LE.ESTOLD )GOTO 120
!
      DO 100 I = 1, N
         X( I ) = SIGN( ONE, X( I ) )
         ISGN( I ) = NINT( X( I ) )
  100 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
!
!     ................ ENTRY   (JUMP = 4)
!     X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X.
!
  110 CONTINUE
      JLAST = J
      J = IDAMAX( N, X, 1 )
      IF( ( X( JLAST ).NE.ABS( X( J ) ) ) .AND. ( ITER.LT.ITMAX ) ) THEN
         ITER = ITER + 1
         goto 50
      END IF
!
!     ITERATION COMPLETE.  FINAL STAGE.
!
  120 CONTINUE
      ALTSGN = ONE
      DO 130 I = 1, N
         X( I ) = ALTSGN*( ONE+DBLE( I-1 ) / DBLE( N-1 ) )
         ALTSGN = -ALTSGN
  130 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
!
!     ................ ENTRY   (JUMP = 5)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
  140 CONTINUE
      TEMP = TWO*( DASUM( N, X, 1 ) / DBLE( 3*N ) )
      IF( TEMP.GT.EST ) THEN
         CALL DCOPY( N, X, 1, V, 1 )
         EST = TEMP
      END IF
!
  150 CONTINUE
      KASE = 0
      RETURN
!
!     End of DLACON
!
      END SUBROUTINE DLACON

!  ***************************************************
       subroutine  dswap (n,dx,incx,dy,incy)
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)goto 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
   20 m = mod(n,3)
      if( m .eq. 0 ) goto 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end  SUBROUTINE DSWAP
	  !****************************************
	        SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  DGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Parameters
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
!     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
!
      RETURN
!
!     End of DGER  .
!
      END SUBROUTINE DGER

!********************************************
      subroutine  dscal(n,da,dx,incx)
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
!
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)goto 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,5)
      if( m .eq. 0 ) goto 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end SUBROUTINE DSCAL
!*************************************************
      SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,BETA, Y, INCY )
!     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )

!  Purpose
!  =======
!
!  DGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.		   N
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.		 BOYUT
!           M must be at least zero.										SATIR A
!           Unchanged on exit.												  
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.	 BOYUT
!           N must be at least zero.										SUTUN A
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.						   1
!           Unchanged on exit.										 		  
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared		  M
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta. When BETA is		   0
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!

!     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!  

 

      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND..NOT.LSAME( TRANS, 'T' ).AND. .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMV ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )   RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO ) RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
!
!        Form  y := alpha*A*x + y.
!
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
!
!        Form  y := alpha*A'*x + y.
!
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of DGEMV .
!
      END  SUBROUTINE DGEMV
!*******************************************
      subroutine  dcopy(n,dx,incx,dy,incy)
!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)goto 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,7)
      if( m .eq. 0 ) goto 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end SUBROUTINE DCOPY
!**************************************************

      SUBROUTINE DSYMV ( UPLO, N, ALPHA, A, LDA, X, INCX,BETA, Y, INCY )
!     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER*1        UPLO
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  DSYMV  performs the matrix-vector  operation
!
!     y := alpha*A*x + beta*y,
!
!  where alpha and beta are scalars, x and y are n element vectors and
!  A is an n by n symmetric matrix.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y. On exit, Y is overwritten by the updated
!           vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND..NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 5
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 7
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYMV ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )RETURN
!
!     Set up the start points in  X  and  Y.
!
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
!     First form  y := beta*y.
!
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, N
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, N
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )RETURN
      IF( LSAME( UPLO, 'U' ) )THEN
!
!        Form  y  when A is stored in upper triangle.
!
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               DO 50, I = 1, J - 1
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2  = TEMP2  + A( I, J )*X( I )
   50          CONTINUE
               Y( J ) = Y( J ) + TEMP1*A( J, J ) + ALPHA*TEMP2
   60       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 80, J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               IX    = KX
               IY    = KY
               DO 70, I = 1, J - 1
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2   = TEMP2   + A( I, J )*X( IX )
                  IX      = IX      + INCX
                  IY      = IY      + INCY
   70          CONTINUE
               Y( JY ) = Y( JY ) + TEMP1*A( J, J ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
   80       CONTINUE
         END IF
      ELSE
!
!        Form  y  when A is stored in lower triangle.
!
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 100, J = 1, N
               TEMP1  = ALPHA*X( J )
               TEMP2  = ZERO
               Y( J ) = Y( J )       + TEMP1*A( J, J )
               DO 90, I = J + 1, N
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2  = TEMP2  + A( I, J )*X( I )
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP2
  100       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 120, J = 1, N
               TEMP1   = ALPHA*X( JX )
               TEMP2   = ZERO
               Y( JY ) = Y( JY )       + TEMP1*A( J, J )
               IX      = JX
               IY      = JY
               DO 110, I = J + 1, N
                  IX      = IX      + INCX
                  IY      = IY      + INCY
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2   = TEMP2   + A( I, J )*X( IX )
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of DSYMV .
!
      END  SUBROUTINE DSYMV
!*********************************************************
      subroutine daxpy(n,da,dx,incx,dy,incy)
!
!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)goto 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,4)
      if( m .eq. 0 ) goto 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end SUBROUTINE DAXPY
!********************************************
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,BETA, C, LDC )
!     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  DGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.		 N
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.		   N
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix	 BOYUT 
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.	 SATIR A
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix    BOYUT 
!           op( B ) and the number of columns of the matrix C. N must be	SUTUN B
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix	BOYUT 
!           op( A ) and the number of rows of the matrix op( B ). K must	SUTUN A
!           be at least  zero.												SATIR B
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.				            1
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared		BOYUT
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then		 M
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.			  
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n		
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared	  BOYUT
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then		K
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is	  0
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared    BETA 0
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
!     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Executable Statements ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.



!
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.( .NOT.LSAME( TRANSA, 'C' ) ).AND.( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.( .NOT.LSAME( TRANSB, 'C' ) ).AND.( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )  RETURN
!
!     And if  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
!
!     Start the operations.
!
      IF( NOTB )THEN
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B + beta*C.
!
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B + beta*C
!
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B' + beta*C
!
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B' + beta*C
!
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of DGEMM .
!
      END SUBROUTINE DGEMM
!*************************************************
      SUBROUTINE DSYR  ( UPLO, N, ALPHA, X, INCX, A, LDA )
!     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, LDA, N
      CHARACTER*1        UPLO
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
!     ..
!
!  Purpose
!  =======
!
!  DSYR   performs the symmetric rank 1 operation
!
!     A := alpha*x*x' + A,
!
!  where alpha is a real scalar, x is an n element vector and A is an
!  n by n symmetric matrix.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced. On exit, the
!           upper triangular part of the array A is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced. On exit, the
!           lower triangular part of the array A is overwritten by the
!           lower triangular part of the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
!     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND..NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYR  ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )RETURN
!
!     Set the start point in X if the increment is not unity.
!
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
      IF( LSAME( UPLO, 'U' ) )THEN
!
!        Form  A  when A is stored in upper triangle.
!
         IF( INCX.EQ.1 )THEN
            DO 20, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO 10, I = 1, J
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE
            JX = KX
            DO 40, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IX   = KX
                  DO 30, I = 1, J
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
   30             CONTINUE
               END IF
               JX = JX + INCX
   40       CONTINUE
         END IF
      ELSE
!
!        Form  A  when A is stored in lower triangle.
!
         IF( INCX.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( J ).NE.ZERO )THEN
                  TEMP = ALPHA*X( J )
                  DO 50, I = J, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
            JX = KX
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IX   = JX
                  DO 70, I = J, N
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of DSYR  .
!
      END SUBROUTINE DSYR
!****************************************************
!************************************************FUNCTIONS*************************
      LOGICAL          FUNCTION LSAME( CA, CB )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          CA, CB
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      LSAME = CA.EQ.CB
      IF( LSAME )RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
!
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.INTA.GE.145 .AND. INTA.LE.153 .OR.INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.INTB.GE.145 .AND. INTB.LE.153 .OR.INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
!
!     RETURN
!
!     End of LSAME
	   END FUNCTION LSAME
!***************************************************
      
      INTEGER    FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,N4 )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR and QZ methods
!               for nonsymmetric eigenvalue problems.
!          = 9: maximum size of the subproblems at the bottom of the
!               computation tree in the divide-and-conquer algorithm
!               (used by xGELSD and xGESDD)
!          =10: ieee NaN arithmetic can be trusted not to trap
!          =11: infinity arithmetic can be trusted not to trap
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. External Functions ..
      INTEGER            IEEECK
      EXTERNAL           IEEECK
!     ..
!     .. Executable Statements ..
!
      goto ( 100, 100, 100, 400, 500, 600, 700, 800, 900, 1000,1100 ) ISPEC
!
!     Invalid value for ISPEC
!
      ILAENV = -1
      RETURN
!
  100 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!        ASCII character set
!
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 ) SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!        EBCDIC character set
!
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.( IC.GE.145 .AND. IC.LE.153 ) .OR.( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.( IC.GE.145 .AND. IC.LE.153 ) .OR. ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!        Prime machines:  ASCII+128
!
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 ) SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
!
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
!
      goto ( 110, 200, 300 ) ISPEC
!
  110 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
!
  200 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
!
  300 CONTINUE
!
!     ISPEC = 3:  crossover point
!
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
!
  400 CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
      ILAENV = 6
      RETURN
!
  500 CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      ILAENV = 2
      RETURN
!
  600 CONTINUE 
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
!
  700 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
      ILAENV = 1
      RETURN
!
  800 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ILAENV = 50
      RETURN
!
  900 CONTINUE
!
!     ISPEC = 9:  maximum size of the subproblems at the bottom of the
!                 computation tree in the divide-and-conquer algorithm
!                 (used by xGELSD and xGESDD)
!
      ILAENV = 25
      RETURN
!
 1000 CONTINUE
!
!     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
!
!     ILAENV = 0	???????????????????????????? DIKKAT
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 ) 
      END IF
      RETURN
!
 1100 CONTINUE
!
!     ISPEC = 11: infinity arithmetic can be trusted not to trap
!
!     ILAENV = 0  ????????????????????????????	DIKKAT
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 ) 
      END IF
      RETURN
!
!     End of ILAENV
!
      END FUNCTION ILAENV
!**********************************************************

      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          CMACH
!     ..
!
!  Purpose
!  =======
!
!  DLAMCH determines double precision machine parameters.
!
!  Arguments
!  =========
!
!  CMACH   (input) CHARACTER*1
!          Specifies the value to be returned by DLAMCH:
!          = 'E' or 'e',   DLAMCH := eps
!          = 'S' or 's ,   DLAMCH := sfmin
!          = 'B' or 'b',   DLAMCH := base
!          = 'P' or 'p',   DLAMCH := eps*base
!          = 'N' or 'n',   DLAMCH := t
!          = 'R' or 'r',   DLAMCH := rnd
!          = 'M' or 'm',   DLAMCH := emin
!          = 'U' or 'u',   DLAMCH := rmin
!          = 'L' or 'l',   DLAMCH := emax
!          = 'O' or 'o',   DLAMCH := rmax
!
!          where
!
!          eps   = relative machine precision
!          sfmin = safe minimum, such that 1/sfmin does not overflow
!          base  = base of the machine
!          prec  = eps*base
!          t     = number of (base) digits in the mantissa
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!          emin  = minimum exponent before (gradual) underflow
!          rmin  = underflow threshold - base**(emin-1)
!          emax  = largest exponent before overflow
!          rmax  = overflow threshold  - (base**emax)*(1-eps)
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,RND, SFMIN, SMALL, T
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAMC2
!     ..
!     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,EMAX, RMAX, PREC
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
!
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
!
      DLAMCH = RMACH
      RETURN
!
!     End of DLAMCH
!
      END  FUNCTION DLAMCH
!
!***********************************************************************
!
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
!     ..
!
!  Purpose
!  =======
!
!  DLAMC1 determines the machine parameters given by BETA, T, RND, and
!  IEEE1.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  IEEE1   (output) LOGICAL
!          Specifies whether rounding appears to be done in the IEEE
!          'round to nearest' style.
!
!  Further Details
!  ===============
!
!  The routine is based on the routine  ENVRON  by Malcolm and
!  incorporates suggestions by Gentleman and Marovich. See
!
!     Malcolm M. A. (1972) Algorithms to reveal properties of
!        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
!
!     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
!        that reveal properties of floating point arithmetic units.
!        Comms. of the ACM, 17, 276-277.
!
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
!     ..
!     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
!
!        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
!        IEEE1, T and RND.
!
!        Throughout this routine  we use the function  DLAMC3  to ensure
!        that relevant values are  stored and not held in registers,  or
!        are not affected by optimizers.
!
!        Compute  a = 2.0**m  with the  smallest positive integer m such
!        that
!
!           fl( a + 1.0 ) = a.
!
         A = 1
         C = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            goto 10
         END IF
!+       END WHILE
!
!        Now compute  b = 2.0**m  with the smallest positive integer m
!        such that
!
!           fl( a + b ) .gt. a.
!
         B = 1
         C = DLAMC3( A, B )
!
!+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            goto 20
         END IF
!+       END WHILE
!
!        Now compute the base.  a and c  are neighbouring floating point
!        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
!        their difference is beta. Adding 0.25 to c is to ensure that it
!        is truncated to beta and not ( beta - 1 ).
!
         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = C + QTR
!
!        Now determine whether rounding or chopping occurs,  by adding a
!        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
!
         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )LRND = .FALSE.
!
!        Try and decide whether rounding is done in the  IEEE  'round to
!        nearest' style. B/2 is half a unit in the last place of the two
!        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
!        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
!        A, but adding B/2 to SAVEC should change SAVEC.
!
         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
!
!        Now find  the  mantissa, t.  It should  be the  integer part of
!        log to the base beta of a,  however it is safer to determine  t
!        by powering.  So we find t as the smallest positive integer for
!        which
!
!           fl( beta**t + 1.0 ) = 1.0.
!
         LT = 0
         A = 1
         C = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            goto 30
         END IF
!+       END WHILE
!
      END IF
!
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
!
!     End of DLAMC1
!
      END SUBROUTINE DLAMC1
!
!************************************************************************
!
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
!     ..
!
!  Purpose
!  =======
!
!  DLAMC2 determines the machine parameters specified in its argument
!  list.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  EPS     (output) DOUBLE PRECISION
!          The smallest positive number such that
!
!             fl( 1.0 - EPS ) .LT. 1.0,
!
!          where fl denotes the computed value.
!
!  EMIN    (output) INTEGER
!          The minimum exponent before (gradual) underflow occurs.
!
!  RMIN    (output) DOUBLE PRECISION
!          The smallest normalized number for the machine, given by
!          BASE**( EMIN - 1 ), where  BASE  is the floating point value
!          of BETA.
!
!  EMAX    (output) INTEGER
!          The maximum exponent before overflow occurs.
!
!  RMAX    (output) DOUBLE PRECISION
!          The largest positive number for the machine, given by
!          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
!          value of BETA.
!
!  Further Details
!  ===============
!
!  The computation of  EPS  is based on a routine PARANOIA by
!  W. Kahan of the University of California at Berkeley.
!
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,SIXTH, SMALL, THIRD, TWO, ZERO
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,LRMIN, LT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
!
!        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
!        BETA, T, RND, EPS, EMIN and RMIN.
!
!        Throughout this routine  we use the function  DLAMC3  to ensure
!        that relevant values are stored  and not held in registers,  or
!        are not affected by optimizers.
!
!        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
!
         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
!
!        Start to find EPS.
!
         B = LBETA
         A = B**( -LT )
         LEPS = A
!
!        Try some tricks to see whether or not this is the correct  EPS.
!
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )B = LEPS
!
         LEPS = 1
!
!+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            goto 10
         END IF
!+       END WHILE
!
         IF( A.LT.LEPS ) LEPS = A
!
!        Computation of EPS complete.
!
!        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
!        Keep dividing  A by BETA until (gradual) underflow occurs. This
!        is detected when we cannot recover the previous A.
!
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
!
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
!            ( Non twos-complement machines, no gradual underflow;
!              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
!            ( Non twos-complement machines, with gradual underflow;
!              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
!            ( Twos-complement machines, no gradual underflow;
!              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
!            ( Twos-complement machines with gradual underflow;
!              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
!         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
!!!
! Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
!!!
!
!        Assume IEEE arithmetic if we found denormalised  numbers above,
!        or if arithmetic seems to round in the  IEEE style,  determined
!        in routine DLAMC1. A true IEEE machine should have both  things
!        true; however, faulty machines may have one or the other.
!
         IEEE = IEEE .OR. LIEEE1
!
!        Compute  RMIN by successive division by  BETA. We could compute
!        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
!        this computation.
!
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
!
!        Finally, call DLAMC5 to compute EMAX and RMAX.
!
         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
!
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
!
      RETURN
!
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-','  EMIN = ', I8, /' If, after inspection, the value EMIN looks',' acceptable please comment out ',/ ' the IF block as marked within the code of routine',' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
!
!     End of DLAMC2
!
      END SUBROUTINE DLAMC2
!
!***********************************************************************
!
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
!     ..
!
!  Purpose
!  =======
!
!  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
!  the addition of  A  and  B ,  for use in situations where optimizers
!  might hold one of these in a register.
!
!  Arguments
!  =========
!
!  A, B    (input) DOUBLE PRECISION
!          The values A and B.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      DLAMC3 = A + B
!
      RETURN
!
!     End of DLAMC3
!
      END FUNCTION DLAMC3
!
!***********************************************************************
!
      SUBROUTINE DLAMC4( EMIN, START, BASE )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
!     ..
!
!  Purpose
!  =======
!
!  DLAMC4 is a service routine for DLAMC2.
!
!  Arguments
!  =========
!
!  EMIN    (output) EMIN
!          The minimum exponent before (gradual) underflow, computed by
!          setting A = START and dividing by BASE until the previous A
!          can not be recovered.
!
!  START   (input) DOUBLE PRECISION
!          The starting point for determining EMIN.
!
!  BASE    (input) INTEGER
!          The base of the machine.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
!     ..
!     .. Executable Statements ..
!
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
!+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
!    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         goto 10
      END IF
!+    END WHILE
!
      RETURN
!
!     End of DLAMC4
!
      END SUBROUTINE DLAMC4
!
!***********************************************************************
!
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
!     ..
!
!  Purpose
!  =======
!
!  DLAMC5 attempts to compute RMAX, the largest machine floating-point
!  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
!  approximately to a power of 2.  It will fail on machines where this
!  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
!  EMAX = 28718).  It will also fail if the value supplied for EMIN is
!  too large (i.e. too close to zero), probably with overflow.
!
!  Arguments
!  =========
!
!  BETA    (input) INTEGER
!          The base of floating-point arithmetic.
!
!  P       (input) INTEGER
!          The number of base BETA digits in the mantissa of a
!          floating-point value.
!
!  EMIN    (input) INTEGER
!          The minimum exponent before (gradual) underflow.
!
!  IEEE    (input) LOGICAL
!          A logical flag specifying whether or not the arithmetic
!          system is thought to comply with the IEEE standard.
!
!  EMAX    (output) INTEGER
!          The largest exponent before overflow
!
!  RMAX    (output) DOUBLE PRECISION
!          The largest machine floating-point number.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
!     ..
!     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..
!
!     First compute LEXP and UEXP, two powers of 2 that bound
!     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
!     approximately to the bound that is closest to abs(EMIN).
!     (EMAX is the exponent of the required number RMAX).
!
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         goto 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
!
!     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
!     than or equal to EMIN. EXBITS is the number of bits needed to
!     store the exponent.
!
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
!
!     EXPSUM is the exponent range, approximately equal to
!     EMAX - EMIN + 1 .
!
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
!
!     NBITS is the total number of bits needed to store a
!     floating-point number.
!
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
!
!        Either there are an odd number of bits used to store a
!        floating-point number, which is unlikely, or some bits are
!        not used in the representation of numbers, which is possible,
!        (e.g. Cray machines) or the mantissa has an implicit bit,
!        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
!        most likely. We have to assume the last alternative.
!        If this is true, then we need to reduce EMAX by one because
!        there must be some way of representing zero in an implicit-bit
!        system. On machines like Cray, we are reducing EMAX by one
!        unnecessarily.
!
         EMAX = EMAX - 1
      END IF
!
      IF( IEEE ) THEN
!
!        Assume we are on an IEEE machine which reserves one exponent
!        for infinity and NaN.
!
         EMAX = EMAX - 1
      END IF
!
!     Now create RMAX, the largest machine number, which should
!     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
!
!     First compute 1.0 - BETA**(-P), being careful that the
!     result is less than 1.0 .
!
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )OLDY = Y
         Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE ) Y = OLDY
!
!     Now multiply by BETA**EMAX to get RMAX.
!
      DO 30 I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE
!
      RMAX = Y
      RETURN
!
!     End of DLAMC5
!
      END  SUBROUTINE DLAMC5
!*********************************************

        DOUBLE PRECISION FUNCTION DLANSY( NORM, UPLO, N, A, LDA, WORK )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLANSY  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real symmetric matrix A.
!
!  Description
!  ===========
!
!  DLANSY returns the value
!
!     DLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in DLANSY as described
!          above.
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is to be referenced.
!          = 'U':  Upper triangular part of A is referenced
!          = 'L':  Lower triangular part of A is referenced
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, DLANSY is
!          set to zero.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The symmetric matrix A.  If UPLO = 'U', the leading n by n
!          upper triangular part of A contains the upper triangular part
!          of the matrix A, and the strictly lower triangular part of A
!          is not referenced.  If UPLO = 'L', the leading n by n lower
!          triangular part of A contains the lower triangular part of
!          the matrix A, and the strictly upper triangular part of A is
!          not referenced.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(N,1).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
!          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
!          WORK is not referenced.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASSQ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 J = 1, N
               DO 10 I = 1, J
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40 J = 1, N
               DO 30 I = J, N
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR.( NORM.EQ.'1' ) ) THEN
!
!        Find normI(A) ( = norm1(A), since A is symmetric).
!
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 60 J = 1, N
               SUM = ZERO
               DO 50 I = 1, J - 1
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   50          CONTINUE
               WORK( J ) = SUM + ABS( A( J, J ) )
   60       CONTINUE
            DO 70 I = 1, N
               VALUE = MAX( VALUE, WORK( I ) )
   70       CONTINUE
         ELSE
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( A( J, J ) )
               DO 90 I = J + 1, N
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   90          CONTINUE
               VALUE = MAX( VALUE, SUM )
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 2, N
               CALL DLASSQ( J-1, A( 1, J ), 1, SCALE, SUM )
  110       CONTINUE
         ELSE
            DO 120 J = 1, N - 1
               CALL DLASSQ( N-J, A( J+1, J ), 1, SCALE, SUM )
  120       CONTINUE
         END IF
         SUM = 2*SUM
         CALL DLASSQ( N, A, LDA+1, SCALE, SUM )
         VALUE = SCALE*SQRT( SUM )
      END IF
!
      DLANSY = VALUE
      RETURN
!
!     End of DLANSY
!
      END FUNCTION DLANSY
!******************************************************
      integer function idamax(n,dx,incx)
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dmax
      integer i,incx,ix,n
!
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)goto 20
!
!        code for increment not equal to 1
!
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) goto 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) goto 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end FUNCTION IDAMAX
!***********************************************
      double precision function dasum(n,dx,incx)
!
!     takes the sum of the absolute values.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dtemp
      integer i,incx,m,mp1,n,nincx
!
      dasum = 0.0d0
      dtemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)goto 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
   10 continue
      dasum = dtemp
      return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,6)
      if( m .eq. 0 ) goto 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) goto 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))+ dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue
   60 dasum = dtemp
      return
      end FUNCTION DASUM
!*******************************************************
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1998
!
!     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
!     ..
!
!  Purpose
!  =======
!
!  IEEECK is called from the ILAENV to verify that Infinity and
!  possibly NaN arithmetic is safe (i.e. will not trap).
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies whether to test just for inifinity arithmetic
!          or whether to test for infinity and NaN arithmetic.
!          = 0: Verify infinity arithmetic only.
!          = 1: Verify infinity and NaN arithmetic.
!
!  ZERO    (input) REAL
!          Must contain the value 0.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  ONE     (input) REAL
!          Must contain the value 1.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  RETURN VALUE:  INTEGER
!          = 0:  Arithmetic failed to produce the correct answers
!          = 1:  Arithmetic produced the correct answers
!
!     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,NEGZRO, NEWZRO, POSINF
!     ..
!     .. Executable Statements ..
      IEEECK = 1
!
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
!
!
!
!     Return if we were only asked to check infinity arithmetic
!
      IF( ISPEC.EQ.0 )RETURN
!
      NAN1 = POSINF + NEGINF
!
      NAN2 = POSINF / NEGINF
!
      NAN3 = POSINF / POSINF
!
      NAN4 = POSINF*ZERO
!
      NAN5 = NEGINF*NEGZRO
!
      NAN6 = NAN5*0.0
!
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      RETURN
      END FUNCTION IEEECK
!****************************************
      SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASSQ  returns the values  scl  and  smsq  such that
!
!     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
!
!  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
!  assumed to be non-negative and  scl  returns the value
!
!     scl = max( scale, abs( x( i ) ) ).
!
!  scale and sumsq must be supplied in SCALE and SUMSQ and
!  scl and smsq are overwritten on SCALE and SUMSQ respectively.
!
!  The routine makes only one pass through the vector x.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements to be used from the vector X.
!
!  X       (input) DOUBLE PRECISION array, dimension (N)
!          The vector for which a scaled sum of squares is computed.
!             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
!
!  INCX    (input) INTEGER
!          The increment between successive values of the vector X.
!          INCX > 0.
!
!  SCALE   (input/output) DOUBLE PRECISION
!          On entry, the value  scale  in the equation above.
!          On exit, SCALE is overwritten with  scl , the scaling factor
!          for the sum of squares.
!
!  SUMSQ   (input/output) DOUBLE PRECISION
!          On entry, the value  sumsq  in the equation above.
!          On exit, SUMSQ is overwritten with  smsq , the basic sum of
!          squares from which  scl  has been factored out.
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   ABSXI
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( X( IX ).NE.ZERO ) THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
!
!     End of DLASSQ
!
      END SUBROUTINE DLASSQ




