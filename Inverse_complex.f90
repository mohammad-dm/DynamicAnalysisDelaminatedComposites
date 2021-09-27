MODULE COMPLEX_INVERSE 
    
!USE    CUSH
!USE    ZGET
!
    implicit none
    PUBLIC:: ZGETRF,zgetri !, zgemv, zswap, ztrsm, zgemm, ztrtri, ztrmm, ztrti2, zscal, ztrmv
    
    CONTAINS 
!
     SUBROUTINE ZGETRF (A, LDA,  N, IPIV, INFO )
 !    
 !ZGETRF computes an LU factorization of a general M-by-N matrix A
 !using partial pivoting with row interchanges.
 !
 !The factorization has the form
 !   A = P * L * U
 !where P is a permutation matrix, L is lower triangular with unit
 !diagonal elements (lower trapezoidal if m > n), and U is upper
 !triangular (upper trapezoidal if m < n).
 !
 !This is the right-looking Level 3 BLAS version of the algorithm.
 !
![in]	M	
!          M is INTEGER
!          The number of rows of the matrix A.  M >= 0.
![in]	N	
!          N is INTEGER
!          The number of columns of the matrix A.  N >= 0.
![in,out]	A	
!          A is COMPLEX*16 array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
![in]	LDA	
!          LDA is INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
![out]	IPIV	
!          IPIV is INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
![out]	INFO	
!          INFO is INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.

 !  -- LAPACK computational routine --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !
 !     .. Scalar Arguments ..
       INTEGER            LDA, M, N
       INTEGER,INTENT(OUT)::INFO
 !     ..
 !     .. Array Arguments ..
       INTEGER,DIMENSION(:), intent(out)::IPIV
       COMPLEX*16,DIMENSION(:,:), intent(inout):: A


 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16         ONE
       parameter( one = ( 1.0d+0, 0.0d+0 ) )
 !     ..
 !     .. Local Scalars ..
       INTEGER            I, IINFO, J, JB, NB
 !     ..
 !     .. External Subroutines ..
       EXTERNAL           xerbla, zgemm, zgetrf2, zlaswp, ztrsm
 !     ..
 !     .. External Functions ..
       INTEGER            ILAENV
       EXTERNAL           ilaenv
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          max, min
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input parameters.
       M=LDA
 !
       info = 0
       IF( m.LT.0 ) THEN
          info = -1
       ELSE IF( n.LT.0 ) THEN
          info = -2
       ELSE IF( lda.LT.max( 1, m ) ) THEN
          info = -4
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'ZGETRF', -info )
          RETURN
       END IF
 !
 !     Quick return if possible
 !
       IF( m.EQ.0 .OR. n.EQ.0 ) RETURN
 !
 !     Determine the block size for this environment.
 !
       nb = ilaenv( 1, 'ZGETRF', ' ', m, n, -1, -1 )
       IF( nb.LE.1 .OR. nb.GE.min( m, n ) ) THEN
 !
 !        Use unblocked code.
 !
          CALL zgetrf2( m, n, a, lda, ipiv, info )
       ELSE
 !
 !        Use blocked code.
 !
          DO 20 j = 1, min( m, n ), nb
             jb = min( min( m, n )-j+1, nb )
 !
 !           Factor diagonal and subdiagonal blocks and test for exact
 !           singularity.
 !
             CALL zgetrf2( m-j+1, jb, a( j, j ), lda, ipiv( j ), iinfo )
 !
 !           Adjust INFO and the pivot indices.
 !
             IF( info.EQ.0 .AND. iinfo.GT.0 ) info = iinfo + j - 1
             DO 10 i = j, min( m, j+jb-1 )
                ipiv( i ) = j - 1 + ipiv( i )
    10       CONTINUE
 !
 !           Apply interchanges to columns 1:J-1.
 !
             CALL zlaswp( j-1, a, lda, j, j+jb-1, ipiv, 1 )
 !
             IF( j+jb.LE.n ) THEN
 !
 !              Apply interchanges to columns J+JB:N.
 !
                CALL zlaswp( n-j-jb+1, a( 1, j+jb ), lda, j, j+jb-1,ipiv, 1 )
 !
 !              Compute block row of U.
 !
                CALL ztrsm( 'Left', 'Lower', 'No transpose', 'Unit', jb,n-j-jb+1, one, a( j, j ), lda, a( j, j+jb ),lda )
                IF( j+jb.LE.m ) THEN
 !
 !                 Update trailing submatrix.
 !
                   CALL zgemm( 'No transpose', 'No transpose', m-j-jb+1,n-j-jb+1, jb, -one, a( j+jb, j ), lda, a( j, j+jb ), lda, one, a( j+jb, j+jb ),lda )
                END IF
             END IF
    20    CONTINUE
       END IF
       RETURN
 !
 !     End of ZGETRF
       
       End SUBROUTINE  ZGETRF
 ! ------------------------------------------------------------------------------------     
      SUBROUTINE zgetri(A, LDA,  N, IPIV, INFO )
 !
 !ZGETRI computes the inverse of a matrix using the LU factorization
 !computed by ZGETRF.
 !
 !This method inverts U and then computes inv(A) by solving the system
 !inv(A)*L = inv(U) for inv(A).
![in]	N	
!          N is INTEGER
!          The order of the matrix A.  N >= 0.
![in,out]	A	
!          A is COMPLEX*16 array, dimension (LDA,N)
!          On entry, the factors L and U from the factorization
!          A = P*L*U as computed by ZGETRF.
!          On exit, if INFO = 0, the inverse of the original matrix A.
![in]	LDA	
!          LDA is INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
![in]	IPIV	
!          IPIV is INTEGER array, dimension (N)
!          The pivot indices from ZGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
![out]	WORK	
!          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
![in]	LWORK	
!          LWORK is INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,N).
!          For optimal performance LWORK >= N*NB, where NB is
!          the optimal blocksize returned by ILAENV.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
![out]	INFO	
!          INFO is INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
!                singular and its inverse could not be computed.
!Author
!Univ. of Tennessee
!Univ. of California Berkeley
!Univ. of Colorado Denver
!NAG Ltd.
 !  -- LAPACK computational routine --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !
 !     .. Scalar Arguments ..
       INTEGER            INFO, LDA, LWORK, N
 !     ..
 !     .. Array Arguments ..
       INTEGER            IPIV( * )
       COMPLEX*16,DIMENSION(:,:), intent(inout):: A
       COMPLEX*16,DIMENSION(:),ALLOCATABLE:: WORK
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16         ZERO, ONE
       parameter( zero = ( 0.0d+0, 0.0d+0 ),one = ( 1.0d+0, 0.0d+0 ) )
 !     ..
 !     .. Local Scalars ..
       LOGICAL            LQUERY
       INTEGER            I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB,NBMIN, NN
 !     ..
 !     .. External Functions ..
       INTEGER            ILAENV
       EXTERNAL           ilaenv
 !     ..
 !     .. External Subroutines ..
       EXTERNAL           xerbla, zgemm, zgemv, zswap, ztrsm, ztrtri
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          max, min
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input parameters.
 !      
       LWORK=65*N
      
       ALLOCATE(WORK(LWORK))

       info = 0
       nb = ilaenv( 1, 'ZGETRI', ' ', n, -1, -1, -1 )
       lwkopt = n   !nb
       work( 1 ) = lwkopt
       lquery = ( lwork.EQ.-1 )
       IF( n.LT.0 ) THEN
          info = -1
       ELSE IF( lda.LT.max( 1, n ) ) THEN
          info = -3
       ELSE IF( lwork.LT.max( 1, n ) .AND. .NOT.lquery ) THEN
          info = -6
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'ZGETRI', -info )
          RETURN
       ELSE IF( lquery ) THEN
          RETURN
       END IF
 !
 !     Quick return if possible
 !
       IF( n.EQ.0 )RETURN
 !
 !     Form inv(U).  If INFO > 0 from ZTRTRI, then U is singular,
 !     and the inverse is not computed.
 !
       CALL ztrtri( 'Upper', 'Non-unit', n, a, lda, info )
       IF( info.GT.0 )RETURN
 !
       nbmin = 2
       ldwork = n
       IF( nb.GT.1 .AND. nb.LT.n ) THEN
          iws = max( ldwork*nb, 1 )
          IF( lwork.LT.iws ) THEN
             nb = lwork / ldwork
             nbmin = max( 2, ilaenv( 2, 'ZGETRI', ' ', n, -1, -1, -1 ) )
          END IF
       ELSE
          iws = n
       END IF
 !
 !     Solve the equation inv(A)!L = inv(U) for inv(A).
 !
       IF( nb.LT.nbmin .OR. nb.GE.n ) THEN
 !
 !        Use unblocked code.
 !
          DO 20 j = n, 1, -1
 !
 !           Copy current column of L to WORK and replace with zeros.
 !
             DO 10 i = j + 1, n
                work( i ) = a( i, j )
                a( i, j ) = zero
    10       CONTINUE
 !
 !           Compute current column of inv(A).
 !
             IF( j.LT.n ) CALL zgemv( 'No transpose', n, n-j, -one, a( 1, j+1 ),lda, work( j+1 ), 1, one, a( 1, j ), 1 )
    20    CONTINUE
       ELSE
 !
 !        Use blocked code.
 !
          nn = ( ( n-1 ) / nb )!nb + 1
          DO 50 j = nn, 1, -nb
             jb = min( nb, n-j+1 )
 !
 !           Copy current block column of L to WORK and replace with
 !           zeros.
 !
             DO 40 jj = j, j + jb - 1
                DO 30 i = jj + 1, n
                   work( i+( jj-j )*ldwork ) = a( i, jj )
                   a( i, jj ) = zero
    30          CONTINUE
    40       CONTINUE
 !
 !           Compute current block column of inv(A).
 !
             IF( j+jb.LE.n )CALL zgemm( 'No transpose', 'No transpose', n, jb,n-j-jb+1, -one, a( 1, j+jb ), lda,work( j+jb ), ldwork, one, a( 1, j ), lda )
             CALL ztrsm( 'Right', 'Lower', 'No transpose', 'Unit', n, jb,one, work( j ), ldwork, a( 1, j ), lda )
    50    CONTINUE
       END IF
 !
 !     Apply column interchanges.
 !
       DO 60 j = n - 1, 1, -1
          jp = ipiv( j )
          IF( jp.NE.j ) CALL zswap( n, a( 1, j ), 1, a( 1, jp ), 1 )
    60 CONTINUE
 !
       work( 1 ) = iws
       RETURN
 !
 !     End of ZGETRI
 !
      END SUBROUTINE
   ! ------------------------------------------------------------------------------------     
    
    END MODULE COMPLEX_INVERSE 
!******************************************************************************************
!------------------------------------------------------------------------------------------
!******************************************************************************************

      !------------------------------- // ZGEMV \\ -------------------------------------------
      
          SUBROUTINE zgemv(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
 !
 !  -- Reference BLAS level2 routine --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !
 !     .. Scalar Arguments ..
       COMPLEX*16 ALPHA,BETA
       INTEGER INCX,INCY,LDA,M,N
       CHARACTER TRANS
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 A(LDA,*),X(*),Y(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16 ONE
       parameter(one= (1.0d+0,0.0d+0))
       COMPLEX*16 ZERO
       parameter(zero= (0.0d+0,0.0d+0))
 !     ..
 !     .. Local Scalars ..
       COMPLEX*16 TEMP
       INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
       LOGICAL NOCONJ
 !     ..
 !     .. External Functions ..
       LOGICAL LSAME
       EXTERNAL lsame
 !     ..
 !     .. External Subroutines ..
       EXTERNAL xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC dconjg,max
 !     ..
 !
 !     Test the input parameters.
 !
       info = 0
       IF (.NOT.lsame(trans,'N') .AND. .NOT.lsame(trans,'T') .AND. .NOT.lsame(trans,'C')) THEN
           info = 1
       ELSE IF (m.LT.0) THEN
           info = 2
       ELSE IF (n.LT.0) THEN
           info = 3
       ELSE IF (lda.LT.max(1,m)) THEN
           info = 6
       ELSE IF (incx.EQ.0) THEN
           info = 8
       ELSE IF (incy.EQ.0) THEN
           info = 11
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('ZGEMV ',info)
           RETURN
       END IF
 !
 !     Quick return if possible.
 !
       IF ((m.EQ.0) .OR. (n.EQ.0) .OR. ((alpha.EQ.zero).AND. (beta.EQ.one))) RETURN
 !
       noconj = lsame(trans,'T')
 !
 !     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
 !     up the start points in  X  and  Y.
 !
       IF (lsame(trans,'N')) THEN
           lenx = n
           leny = m
       ELSE
           lenx = m
           leny = n
       END IF
       IF (incx.GT.0) THEN
           kx = 1
       ELSE
           kx = 1 - (lenx-1)*incx
       END IF
       IF (incy.GT.0) THEN
           ky = 1
       ELSE
           ky = 1 - (leny-1)*incy
       END IF
 !
 !     Start the operations. In this version the elements of A are
 !     accessed sequentially with one pass through A.
 !
 !     First form  y := beta!y.
 !
       IF (beta.NE.one) THEN
           IF (incy.EQ.1) THEN
               IF (beta.EQ.zero) THEN
                   DO 10 i = 1,leny
                       y(i) = zero
    10             CONTINUE
               ELSE
                   DO 20 i = 1,leny
                       y(i) = beta*y(i)
    20             CONTINUE
               END IF
           ELSE
               iy = ky
               IF (beta.EQ.zero) THEN
                   DO 30 i = 1,leny
                       y(iy) = zero
                       iy = iy + incy
    30             CONTINUE
               ELSE
                   DO 40 i = 1,leny
                       y(iy) = beta*y(iy)
                       iy = iy + incy
    40             CONTINUE
               END IF
           END IF
       END IF
       IF (alpha.EQ.zero) RETURN
       IF (lsame(trans,'N')) THEN
 !
 !        Form  y := alpha*A*x + y.
 !
           jx = kx
           IF (incy.EQ.1) THEN
               DO 60 j = 1,n
                   temp = alpha*x(jx)
                   DO 50 i = 1,m
                       y(i) = y(i) + temp*a(i,j)
    50             CONTINUE
                   jx = jx + incx
    60         CONTINUE
           ELSE
               DO 80 j = 1,n
                   temp = alpha*x(jx)
                   iy = ky
                   DO 70 i = 1,m
                       y(iy) = y(iy) + temp*a(i,j)
                       iy = iy + incy
    70             CONTINUE
                   jx = jx + incx
    80         CONTINUE
           END IF
       ELSE
 !
 !        Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.
 !
           jy = ky
           IF (incx.EQ.1) THEN
               DO 110 j = 1,n
                   temp = zero
                   IF (noconj) THEN
                       DO 90 i = 1,m
                           temp = temp + a(i,j)*x(i)
    90                 CONTINUE
                   ELSE
                       DO 100 i = 1,m
                           temp = temp + dconjg(a(i,j))*x(i)
   100                 CONTINUE
                   END IF
                   y(jy) = y(jy) + alpha*temp
                   jy = jy + incy
   110         CONTINUE
           ELSE
               DO 140 j = 1,n
                   temp = zero
                   ix = kx
                   IF (noconj) THEN
                       DO 120 i = 1,m
                           temp = temp + a(i,j)*x(ix)
                           ix = ix + incx
   120                 CONTINUE
                   ELSE
                       DO 130 i = 1,m
                           temp = temp + dconjg(a(i,j))*x(ix)
                           ix = ix + incx
   130                 CONTINUE
                   END IF
                   y(jy) = y(jy) + alpha*temp
                   jy = jy + incy
   140         CONTINUE
           END IF
       END IF
 !
       RETURN
 !
 !     End of ZGEMV .
 !
    END SUBROUTINE zgemv
    
    !====--=--=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    
          SUBROUTINE zswap(N,ZX,INCX,ZY,INCY)
 !
 !  -- Reference BLAS level1 routine --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !
 !     .. Scalar Arguments ..
       INTEGER INCX,INCY,N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 ZX(*),ZY(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Local Scalars ..
       COMPLEX*16 ZTEMP
       INTEGER I,IX,IY
 !     ..
       IF (n.LE.0) RETURN
       IF (incx.EQ.1 .AND. incy.EQ.1) THEN
 !
 !       code for both increments equal to 1
          DO i = 1,n
             ztemp = zx(i)
             zx(i) = zy(i)
             zy(i) = ztemp
          END DO
       ELSE
 !
 !       code for unequal increments or equal increments not equal
 !         to 1
 !
          ix = 1
          iy = 1
          IF (incx.LT.0) ix = (-n+1)*incx + 1
          IF (incy.LT.0) iy = (-n+1)*incy + 1
          DO i = 1,n
             ztemp = zx(ix)
             zx(ix) = zy(iy)
             zy(iy) = ztemp
             ix = ix + incx
             iy = iy + incy
          END DO
       END IF
       RETURN
          END subroutine zswap
          
    !-------------------------------------------------------------------------------
          
    SUBROUTINE ztrsm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
 !
 !  -- Reference BLAS level3 routine --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !
 !     .. Scalar Arguments ..
       COMPLEX*16 ALPHA
       INTEGER LDA,LDB,M,N
       CHARACTER DIAG,SIDE,TRANSA,UPLO
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 A(LDA,*),B(LDB,*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. External Functions ..
       LOGICAL LSAME
       EXTERNAL lsame
 !     ..
 !     .. External Subroutines ..
       EXTERNAL xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC dconjg,max
 !     ..
 !     .. Local Scalars ..
       COMPLEX*16 TEMP
       INTEGER I,INFO,J,K,NROWA
       LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
 !     ..
 !     .. Parameters ..
       COMPLEX*16 ONE
       parameter(one= (1.0d+0,0.0d+0))
       COMPLEX*16 ZERO
       parameter(zero= (0.0d+0,0.0d+0))
 !     ..
 !
 !     Test the input parameters.
 !
       lside = lsame(side,'L')
       IF (lside) THEN
           nrowa = m
       ELSE
           nrowa = n
       END IF
       noconj = lsame(transa,'T')
       nounit = lsame(diag,'N')
       upper = lsame(uplo,'U')
 !
       info = 0
       IF ((.NOT.lside) .AND. (.NOT.lsame(side,'R'))) THEN
           info = 1
       ELSE IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
           info = 2
       ELSE IF ((.NOT.lsame(transa,'N')) .AND. (.NOT.lsame(transa,'T')) .AND. (.NOT.lsame(transa,'C'))) THEN
           info = 3
       ELSE IF ((.NOT.lsame(diag,'U')) .AND. (.NOT.lsame(diag,'N'))) THEN
           info = 4
       ELSE IF (m.LT.0) THEN
           info = 5
       ELSE IF (n.LT.0) THEN
           info = 6
       ELSE IF (lda.LT.max(1,nrowa)) THEN
           info = 9
       ELSE IF (ldb.LT.max(1,m)) THEN
           info = 11
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('ZTRSM ',info)
           RETURN
       END IF
 !
 !     Quick return if possible.
 !
       IF (m.EQ.0 .OR. n.EQ.0) RETURN
 !
 !     And when  alpha.eq.zero.
 !
       IF (alpha.EQ.zero) THEN
           DO 20 j = 1,n
               DO 10 i = 1,m
                   b(i,j) = zero
    10         CONTINUE
    20     CONTINUE
           RETURN
       END IF
 !
 !     Start the operations.
 !
       IF (lside) THEN
           IF (lsame(transa,'N')) THEN
 !
 !           Form  B := alpha*inv( A )!B.
 !
               IF (upper) THEN
                   DO 60 j = 1,n
                       IF (alpha.NE.one) THEN
                           DO 30 i = 1,m
                               b(i,j) = alpha*b(i,j)
    30                     CONTINUE
                       END IF
                       DO 50 k = m,1,-1
                           IF (b(k,j).NE.zero) THEN
                               IF (nounit) b(k,j) = b(k,j)/a(k,k)
                               DO 40 i = 1,k - 1
                                   b(i,j) = b(i,j) - b(k,j)*a(i,k)
    40                         CONTINUE
                           END IF
    50                 CONTINUE
    60             CONTINUE
               ELSE
                   DO 100 j = 1,n
                       IF (alpha.NE.one) THEN
                           DO 70 i = 1,m
                               b(i,j) = alpha*b(i,j)
    70                     CONTINUE
                       END IF
                       DO 90 k = 1,m
                           IF (b(k,j).NE.zero) THEN
                               IF (nounit) b(k,j) = b(k,j)/a(k,k)
                               DO 80 i = k + 1,m
                                   b(i,j) = b(i,j) - b(k,j)*a(i,k)
    80                         CONTINUE
                           END IF
    90                 CONTINUE
   100             CONTINUE
               END IF
           ELSE
 !
 !           Form  B := alpha!inv( A!!T )!B
 !           or    B := alpha!inv( A!!H )!B.
 !
               IF (upper) THEN
                   DO 140 j = 1,n
                       DO 130 i = 1,m
                           temp = alpha*b(i,j)
                           IF (noconj) THEN
                               DO 110 k = 1,i - 1
                                   temp = temp - a(k,i)*b(k,j)
   110                         CONTINUE
                               IF (nounit) temp = temp/a(i,i)
                           ELSE
                               DO 120 k = 1,i - 1
                                   temp = temp - dconjg(a(k,i))*b(k,j)
   120                         CONTINUE
                               IF (nounit) temp = temp/dconjg(a(i,i))
                           END IF
                           b(i,j) = temp
   130                 CONTINUE
   140             CONTINUE
               ELSE
                   DO 180 j = 1,n
                       DO 170 i = m,1,-1
                           temp = alpha*b(i,j)
                           IF (noconj) THEN
                               DO 150 k = i + 1,m
                                   temp = temp - a(k,i)*b(k,j)
   150                         CONTINUE
                               IF (nounit) temp = temp/a(i,i)
                           ELSE
                               DO 160 k = i + 1,m
                                   temp = temp - dconjg(a(k,i))*b(k,j)
   160                         CONTINUE
                               IF (nounit) temp = temp/dconjg(a(i,i))
                           END IF
                           b(i,j) = temp
   170                 CONTINUE
   180             CONTINUE
               END IF
           END IF
       ELSE
           IF (lsame(transa,'N')) THEN
 !
 !           Form  B := alpha!B!inv( A ).
 !
               IF (upper) THEN
                   DO 230 j = 1,n
                       IF (alpha.NE.one) THEN
                           DO 190 i = 1,m
                               b(i,j) = alpha*b(i,j)
   190                     CONTINUE
                       END IF
                       DO 210 k = 1,j - 1
                           IF (a(k,j).NE.zero) THEN
                               DO 200 i = 1,m
                                   b(i,j) = b(i,j) - a(k,j)*b(i,k)
   200                         CONTINUE
                           END IF
   210                 CONTINUE
                       IF (nounit) THEN
                           temp = one/a(j,j)
                           DO 220 i = 1,m
                               b(i,j) = temp*b(i,j)
   220                     CONTINUE
                       END IF
   230             CONTINUE
               ELSE
                   DO 280 j = n,1,-1
                       IF (alpha.NE.one) THEN
                           DO 240 i = 1,m
                               b(i,j) = alpha*b(i,j)
   240                     CONTINUE
                       END IF
                       DO 260 k = j + 1,n
                           IF (a(k,j).NE.zero) THEN
                               DO 250 i = 1,m
                                   b(i,j) = b(i,j) - a(k,j)*b(i,k)
   250                         CONTINUE
                           END IF
   260                 CONTINUE
                       IF (nounit) THEN
                           temp = one/a(j,j)
                           DO 270 i = 1,m
                               b(i,j) = temp*b(i,j)
   270                     CONTINUE
                       END IF
   280             CONTINUE
               END IF
           ELSE
 !
 !           Form  B := alpha*B*inv( A**T )
 !           or    B := alpha*B*inv( A**H ).
 !
               IF (upper) THEN
                   DO 330 k = n,1,-1
                       IF (nounit) THEN
                           IF (noconj) THEN
                               temp = one/a(k,k)
                           ELSE
                               temp = one/dconjg(a(k,k))
                           END IF
                           DO 290 i = 1,m
                               b(i,k) = temp*b(i,k)
   290                     CONTINUE
                       END IF
                       DO 310 j = 1,k - 1
                           IF (a(j,k).NE.zero) THEN
                               IF (noconj) THEN
                                   temp = a(j,k)
                               ELSE
                                   temp = dconjg(a(j,k))
                               END IF
                               DO 300 i = 1,m
                                   b(i,j) = b(i,j) - temp*b(i,k)
   300                         CONTINUE
                           END IF
   310                 CONTINUE
                       IF (alpha.NE.one) THEN
                           DO 320 i = 1,m
                               b(i,k) = alpha*b(i,k)
   320                     CONTINUE
                       END IF
   330             CONTINUE
               ELSE
                   DO 380 k = 1,n
                       IF (nounit) THEN
                           IF (noconj) THEN
                               temp = one/a(k,k)
                           ELSE
                               temp = one/dconjg(a(k,k))
                           END IF
                           DO 340 i = 1,m
                               b(i,k) = temp*b(i,k)
   340                     CONTINUE
                       END IF
                       DO 360 j = k + 1,n
                           IF (a(j,k).NE.zero) THEN
                               IF (noconj) THEN
                                   temp = a(j,k)
                               ELSE
                                   temp = dconjg(a(j,k))
                               END IF
                               DO 350 i = 1,m
                                   b(i,j) = b(i,j) - temp*b(i,k)
   350                         CONTINUE
                           END IF
   360                 CONTINUE
                       IF (alpha.NE.one) THEN
                           DO 370 i = 1,m
                               b(i,k) = alpha*b(i,k)
   370                     CONTINUE
                       END IF
   380             CONTINUE
               END IF
           END IF
       END IF
 !
       RETURN
 !
 !     End of ZTRSM .
 !
    END subroutine ztrsm
    
! ------------------------------------------------------------------------------------     

          SUBROUTINE zgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
 !
 !  -- Reference BLAS level3 routine --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !
 !     .. Scalar Arguments ..
       COMPLEX*16 ALPHA,BETA
       INTEGER K,LDA,LDB,LDC,M,N
       CHARACTER TRANSA,TRANSB
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. External Functions ..
       LOGICAL LSAME
       EXTERNAL lsame
 !     ..
 !     .. External Subroutines ..
       EXTERNAL xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC dconjg,max
 !     ..
 !     .. Local Scalars ..
       COMPLEX*16 TEMP
       INTEGER I,INFO,J,L,NROWA,NROWB
       LOGICAL CONJA,CONJB,NOTA,NOTB
 !     ..
 !     .. Parameters ..
       COMPLEX*16 ONE
       parameter(one= (1.0d+0,0.0d+0))
       COMPLEX*16 ZERO
       parameter(zero= (0.0d+0,0.0d+0))
 !     ..
 !
 !     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
 !     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
 !     B  respectively are to be  transposed but  not conjugated  and set
 !     NROWA and NROWB  as the number of rows  of  A  and  B  respectively.
 !
       nota = lsame(transa,'N')
       notb = lsame(transb,'N')
       conja = lsame(transa,'C')
       conjb = lsame(transb,'C')
       IF (nota) THEN
           nrowa = m
       ELSE
           nrowa = k
       END IF
       IF (notb) THEN
           nrowb = k
       ELSE
           nrowb = n
       END IF
 !
 !     Test the input parameters.
 !
       info = 0
       IF ((.NOT.nota) .AND. (.NOT.conja) .AND.(.NOT.lsame(transa,'T'))) THEN
           info = 1
       ELSE IF ((.NOT.notb) .AND. (.NOT.conjb) .AND.(.NOT.lsame(transb,'T'))) THEN
           info = 2
       ELSE IF (m.LT.0) THEN
           info = 3
       ELSE IF (n.LT.0) THEN
           info = 4
       ELSE IF (k.LT.0) THEN
           info = 5
       ELSE IF (lda.LT.max(1,nrowa)) THEN
           info = 8
       ELSE IF (ldb.LT.max(1,nrowb)) THEN
           info = 10
       ELSE IF (ldc.LT.max(1,m)) THEN
           info = 13
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('ZGEMM ',info)
           RETURN
       END IF
 !
 !     Quick return if possible.
 !
       IF ((m.EQ.0) .OR. (n.EQ.0) .OR.(((alpha.EQ.zero).OR. (k.EQ.0)).AND. (beta.EQ.one))) RETURN
 !
 !     And when  alpha.eq.zero.
 !
       IF (alpha.EQ.zero) THEN
           IF (beta.EQ.zero) THEN
               DO 20 j = 1,n
                   DO 10 i = 1,m
                       c(i,j) = zero
    10             CONTINUE
    20         CONTINUE
           ELSE
               DO 40 j = 1,n
                   DO 30 i = 1,m
                       c(i,j) = beta*c(i,j)
    30             CONTINUE
    40         CONTINUE
           END IF
           RETURN
       END IF
 !
 !     Start the operations.
 !
       IF (notb) THEN
           IF (nota) THEN
 !
 !           Form  C := alpha!A!B + beta!C.
 !
               DO 90 j = 1,n
                   IF (beta.EQ.zero) THEN
                       DO 50 i = 1,m
                           c(i,j) = zero
    50                 CONTINUE
                   ELSE IF (beta.NE.one) THEN
                       DO 60 i = 1,m
                           c(i,j) = beta*c(i,j)
    60                 CONTINUE
                   END IF
                   DO 80 l = 1,k
                       temp = alpha*b(l,j)
                       DO 70 i = 1,m
                           c(i,j) = c(i,j) + temp*a(i,l)
    70                 CONTINUE
    80             CONTINUE
    90         CONTINUE
           ELSE IF (conja) THEN
 !
 !           Form  C := alpha*A**H*B + beta*C.
 !
               DO 120 j = 1,n
                   DO 110 i = 1,m
                       temp = zero
                       DO 100 l = 1,k
                           temp = temp + dconjg(a(l,i))*b(l,j)
   100                 CONTINUE
                       IF (beta.EQ.zero) THEN
                           c(i,j) = alpha*temp
                       ELSE
                           c(i,j) = alpha*temp + beta*c(i,j)
                       END IF
   110             CONTINUE
   120         CONTINUE
           ELSE
 !
 !           Form  C := alpha*A**T*B + beta*C
 !
               DO 150 j = 1,n
                   DO 140 i = 1,m
                       temp = zero
                       DO 130 l = 1,k
                           temp = temp + a(l,i)*b(l,j)
   130                 CONTINUE
                       IF (beta.EQ.zero) THEN
                           c(i,j) = alpha*temp
                       ELSE
                           c(i,j) = alpha*temp + beta*c(i,j)
                       END IF
   140             CONTINUE
   150         CONTINUE
           END IF
       ELSE IF (nota) THEN
           IF (conjb) THEN
 !
 !           Form  C := alpha*A*B**H + beta*C.
 !
               DO 200 j = 1,n
                   IF (beta.EQ.zero) THEN
                       DO 160 i = 1,m
                           c(i,j) = zero
   160                 CONTINUE
                   ELSE IF (beta.NE.one) THEN
                       DO 170 i = 1,m
                           c(i,j) = beta*c(i,j)
   170                 CONTINUE
                   END IF
                   DO 190 l = 1,k
                       temp = alpha*dconjg(b(j,l))
                       DO 180 i = 1,m
                           c(i,j) = c(i,j) + temp*a(i,l)
   180                 CONTINUE
   190             CONTINUE
   200         CONTINUE
           ELSE
 !
 !           Form  C := alpha!A!B!!T + beta!C
 !
               DO 250 j = 1,n
                   IF (beta.EQ.zero) THEN
                       DO 210 i = 1,m
                           c(i,j) = zero
   210                 CONTINUE
                   ELSE IF (beta.NE.one) THEN
                       DO 220 i = 1,m
                           c(i,j) = beta*c(i,j)
   220                 CONTINUE
                   END IF
                   DO 240 l = 1,k
                       temp = alpha*b(j,l)
                       DO 230 i = 1,m
                           c(i,j) = c(i,j) + temp*a(i,l)
   230                 CONTINUE
   240             CONTINUE
   250         CONTINUE
           END IF
       ELSE IF (conja) THEN
           IF (conjb) THEN
 !
 !           Form  C := alpha!A!!H!B!!H + beta!C.
 !
               DO 280 j = 1,n
                   DO 270 i = 1,m
                       temp = zero
                       DO 260 l = 1,k
                           temp = temp + dconjg(a(l,i))*dconjg(b(j,l))
   260                 CONTINUE
                       IF (beta.EQ.zero) THEN
                           c(i,j) = alpha*temp
                       ELSE
                           c(i,j) = alpha*temp + beta*c(i,j)
                       END IF
   270             CONTINUE
   280         CONTINUE
           ELSE
 !
 !           Form  C := alpha!A!!H!B!!T + beta!C
 !
               DO 310 j = 1,n
                   DO 300 i = 1,m
                       temp = zero
                       DO 290 l = 1,k
                           temp = temp + dconjg(a(l,i))*b(j,l)
   290                 CONTINUE
                       IF (beta.EQ.zero) THEN
                           c(i,j) = alpha*temp
                       ELSE
                           c(i,j) = alpha*temp + beta*c(i,j)
                       END IF
   300             CONTINUE
   310         CONTINUE
           END IF
       ELSE
           IF (conjb) THEN
 !
 !           Form  C := alpha!A!!T!B!!H + beta!C
 !
               DO 340 j = 1,n
                   DO 330 i = 1,m
                       temp = zero
                       DO 320 l = 1,k
                           temp = temp + a(l,i)*dconjg(b(j,l))
   320                 CONTINUE
                       IF (beta.EQ.zero) THEN
                           c(i,j) = alpha*temp
                       ELSE
                           c(i,j) = alpha*temp + beta*c(i,j)
                       END IF
   330             CONTINUE
   340         CONTINUE
           ELSE
 !
 !           Form  C := alpha!A!!T!B!!T + beta!C
 !
               DO 370 j = 1,n
                   DO 360 i = 1,m
                       temp = zero
                       DO 350 l = 1,k
                           temp = temp + a(l,i)*b(j,l)
   350                 CONTINUE
                       IF (beta.EQ.zero) THEN
                           c(i,j) = alpha*temp
                       ELSE
                           c(i,j) = alpha*temp + beta*c(i,j)
                       END IF
   360             CONTINUE
   370         CONTINUE
           END IF
       END IF
 !
       RETURN
 !
 !     End of ZGEMM .
 !
       END subroutine zgemm

 ! ------------------------------------------------------------------------------------     

      SUBROUTINE ztrtri( UPLO, DIAG, N, A, LDA, INFO )
 !
 !      ZTRTRI computes the inverse of a complex upper or lower triangular
 !matrix A.
 !
 !This is the Level 3 BLAS version of the algorithm.
 !  -- LAPACK computational routine --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !
 !     .. Scalar Arguments ..
       CHARACTER          DIAG, UPLO
       INTEGER            INFO, LDA, N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16         A( LDA, * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16         ONE, ZERO
       parameter( one = ( 1.0d+0, 0.0d+0 ),zero = ( 0.0d+0, 0.0d+0 ) )
 !     ..
 !     .. Local Scalars ..
       LOGICAL            NOUNIT, UPPER
       INTEGER            J, JB, NB, NN
 !     ..
 !     .. External Functions ..
       LOGICAL            LSAME
       INTEGER            ILAENV
       EXTERNAL           lsame, ilaenv
 !     ..
 !     .. External Subroutines ..
       EXTERNAL           xerbla, ztrmm, ztrsm, ztrti2
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          max, min
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input parameters.
 !
       info = 0
       upper = lsame( uplo, 'U' )
       nounit = lsame( diag, 'N' )
       IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
          info = -1
       ELSE IF( .NOT.nounit .AND. .NOT.lsame( diag, 'U' ) ) THEN
          info = -2
       ELSE IF( n.LT.0 ) THEN
          info = -3
       ELSE IF( lda.LT.max( 1, n ) ) THEN
          info = -5
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'ZTRTRI', -info )
          RETURN
       END IF
 !
 !     Quick return if possible
 !
       IF( n.EQ.0 )RETURN
 !
 !     Check for singularity if non-unit.
 !
       IF( nounit ) THEN
          DO 10 info = 1, n
             IF( a( info, info ).EQ.zero ) RETURN
    10    CONTINUE
          info = 0
       END IF
 !
 !     Determine the block size for this environment.
 !
       nb = ilaenv( 1, 'ZTRTRI', uplo // diag, n, -1, -1, -1 )
       IF( nb.LE.1 .OR. nb.GE.n ) THEN
 !
 !        Use unblocked code
 !
          CALL ztrti2( uplo, diag, n, a, lda, info )
       ELSE
 !
 !        Use blocked code
 !
          IF( upper ) THEN
 !
 !           Compute inverse of upper triangular matrix
 !
             DO 20 j = 1, n, nb
                jb = min( nb, n-j+1 )
 !
 !              Compute rows 1:j-1 of current block column
 !
                CALL ztrmm( 'Left', 'Upper', 'No transpose', diag, j-1,jb, one, a, lda, a( 1, j ), lda )
                CALL ztrsm( 'Right', 'Upper', 'No transpose', diag, j-1,jb, -one, a( j, j ), lda, a( 1, j ), lda )
 !
 !              Compute inverse of current diagonal block
 !
                CALL ztrti2( 'Upper', diag, jb, a( j, j ), lda, info )
    20       CONTINUE
          ELSE
 !
 !           Compute inverse of lower triangular matrix
 !
             nn = ( ( n-1 ) / nb )!nb + 1
             DO 30 j = nn, 1, -nb
                jb = min( nb, n-j+1 )
                IF( j+jb.LE.n ) THEN
 !
 !                 Compute rows j+jb:n of current block column
 !
                   CALL ztrmm( 'Left', 'Lower', 'No transpose', diag,n-j-jb+1, jb, one, a( j+jb, j+jb ), lda,a( j+jb, j ), lda )
                   CALL ztrsm( 'Right', 'Lower', 'No transpose', diag,n-j-jb+1, jb, -one, a( j, j ), lda,a( j+jb, j ), lda )
                END IF
 !
 !              Compute inverse of current diagonal block
 !
                CALL ztrti2( 'Lower', diag, jb, a( j, j ), lda, info )
    30       CONTINUE
          END IF
       END IF
 !
       RETURN
 !
 !     End of ZTRTRI
 !
    END subroutine ztrtri
    
 ! ------------------------------------------------------------------------------------     
      
 SUBROUTINE ztrmm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
 !
 !  -- Reference BLAS level3 routine --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !
 !     .. Scalar Arguments ..
       COMPLEX*16 ALPHA
       INTEGER LDA,LDB,M,N
       CHARACTER DIAG,SIDE,TRANSA,UPLO
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 A(LDA,*),B(LDB,*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. External Functions ..
       LOGICAL LSAME
       EXTERNAL lsame
 !     ..
 !     .. External Subroutines ..
       EXTERNAL xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC dconjg,max
 !     ..
 !     .. Local Scalars ..
       COMPLEX*16 TEMP
       INTEGER I,INFO,J,K,NROWA
       LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
 !     ..
 !     .. Parameters ..
       COMPLEX*16 ONE
       parameter(one= (1.0d+0,0.0d+0))
       COMPLEX*16 ZERO
       parameter(zero= (0.0d+0,0.0d+0))
 !     ..
 !
 !     Test the input parameters.
 !
       lside = lsame(side,'L')
       IF (lside) THEN
           nrowa = m
       ELSE
           nrowa = n
       END IF
       noconj = lsame(transa,'T')
       nounit = lsame(diag,'N')
       upper = lsame(uplo,'U')
 !
       info = 0
       IF ((.NOT.lside) .AND. (.NOT.lsame(side,'R'))) THEN
           info = 1
       ELSE IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
           info = 2
       ELSE IF ((.NOT.lsame(transa,'N')) .AND.(.NOT.lsame(transa,'T')) .AND.(.NOT.lsame(transa,'C'))) THEN
           info = 3
       ELSE IF ((.NOT.lsame(diag,'U')) .AND. (.NOT.lsame(diag,'N'))) THEN
           info = 4
       ELSE IF (m.LT.0) THEN
           info = 5
       ELSE IF (n.LT.0) THEN
           info = 6
       ELSE IF (lda.LT.max(1,nrowa)) THEN
           info = 9
       ELSE IF (ldb.LT.max(1,m)) THEN
           info = 11
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('ZTRMM ',info)
           RETURN
       END IF
 !
 !     Quick return if possible.
 !
       IF (m.EQ.0 .OR. n.EQ.0) RETURN
 !
 !     And when  alpha.eq.zero.
 !
       IF (alpha.EQ.zero) THEN
           DO 20 j = 1,n
               DO 10 i = 1,m
                   b(i,j) = zero
    10         CONTINUE
    20     CONTINUE
           RETURN
       END IF
 !
 !     Start the operations.
 !
       IF (lside) THEN
           IF (lsame(transa,'N')) THEN
 !
 !           Form  B := alpha*A*B.
 !
               IF (upper) THEN
                   DO 50 j = 1,n
                       DO 40 k = 1,m
                           IF (b(k,j).NE.zero) THEN
                               temp = alpha*b(k,j)
                               DO 30 i = 1,k - 1
                                   b(i,j) = b(i,j) + temp*a(i,k)
    30                         CONTINUE
                               IF (nounit) temp = temp*a(k,k)
                               b(k,j) = temp
                           END IF
    40                 CONTINUE
    50             CONTINUE
               ELSE
                   DO 80 j = 1,n
                       DO 70 k = m,1,-1
                           IF (b(k,j).NE.zero) THEN
                               temp = alpha*b(k,j)
                               b(k,j) = temp
                               IF (nounit) b(k,j) = b(k,j)*a(k,k)
                               DO 60 i = k + 1,m
                                   b(i,j) = b(i,j) + temp*a(i,k)
    60                         CONTINUE
                           END IF
    70                 CONTINUE
    80             CONTINUE
               END IF
           ELSE
 !
 !           Form  B := alpha*A**T*B   or   B := alpha*A**H*B.
 !
               IF (upper) THEN
                   DO 120 j = 1,n
                       DO 110 i = m,1,-1
                           temp = b(i,j)
                           IF (noconj) THEN
                               IF (nounit) temp = temp*a(i,i)
                               DO 90 k = 1,i - 1
                                   temp = temp + a(k,i)*b(k,j)
    90                         CONTINUE
                           ELSE
                               IF (nounit) temp = temp*dconjg(a(i,i))
                               DO 100 k = 1,i - 1
                                   temp = temp + dconjg(a(k,i))*b(k,j)
   100                         CONTINUE
                           END IF
                           b(i,j) = alpha*temp
   110                 CONTINUE
   120             CONTINUE
               ELSE
                   DO 160 j = 1,n
                       DO 150 i = 1,m
                           temp = b(i,j)
                           IF (noconj) THEN
                               IF (nounit) temp = temp*a(i,i)
                               DO 130 k = i + 1,m
                                   temp = temp + a(k,i)*b(k,j)
   130                         CONTINUE
                           ELSE
                               IF (nounit) temp = temp*dconjg(a(i,i))
                               DO 140 k = i + 1,m
                                   temp = temp + dconjg(a(k,i))*b(k,j)
   140                         CONTINUE
                           END IF
                           b(i,j) = alpha*temp
   150                 CONTINUE
   160             CONTINUE
               END IF
           END IF
       ELSE
           IF (lsame(transa,'N')) THEN
 !
 !           Form  B := alpha*B*A.
 !
               IF (upper) THEN
                   DO 200 j = n,1,-1
                       temp = alpha
                       IF (nounit) temp = temp*a(j,j)
                       DO 170 i = 1,m
                           b(i,j) = temp*b(i,j)
   170                 CONTINUE
                       DO 190 k = 1,j - 1
                           IF (a(k,j).NE.zero) THEN
                               temp = alpha*a(k,j)
                               DO 180 i = 1,m
                                   b(i,j) = b(i,j) + temp*b(i,k)
   180                         CONTINUE
                           END IF
   190                 CONTINUE
   200             CONTINUE
               ELSE
                   DO 240 j = 1,n
                       temp = alpha
                       IF (nounit) temp = temp*a(j,j)
                       DO 210 i = 1,m
                           b(i,j) = temp*b(i,j)
   210                 CONTINUE
                       DO 230 k = j + 1,n
                           IF (a(k,j).NE.zero) THEN
                               temp = alpha*a(k,j)
                               DO 220 i = 1,m
                                   b(i,j) = b(i,j) + temp*b(i,k)
   220                         CONTINUE
                           END IF
   230                 CONTINUE
   240             CONTINUE
               END IF
           ELSE
 !
 !           Form  B := alpha*B*A**T   or   B := alpha*B*A**H.
 !
               IF (upper) THEN
                   DO 280 k = 1,n
                       DO 260 j = 1,k - 1
                           IF (a(j,k).NE.zero) THEN
                               IF (noconj) THEN
                                   temp = alpha*a(j,k)
                               ELSE
                                   temp = alpha*dconjg(a(j,k))
                               END IF
                               DO 250 i = 1,m
                                   b(i,j) = b(i,j) + temp*b(i,k)
   250                         CONTINUE
                           END IF
   260                 CONTINUE
                       temp = alpha
                       IF (nounit) THEN
                           IF (noconj) THEN
                               temp = temp*a(k,k)
                           ELSE
                               temp = temp*dconjg(a(k,k))
                           END IF
                       END IF
                       IF (temp.NE.one) THEN
                           DO 270 i = 1,m
                               b(i,k) = temp*b(i,k)
   270                     CONTINUE
                       END IF
   280             CONTINUE
               ELSE
                   DO 320 k = n,1,-1
                       DO 300 j = k + 1,n
                           IF (a(j,k).NE.zero) THEN
                               IF (noconj) THEN
                                   temp = alpha*a(j,k)
                               ELSE
                                   temp = alpha*dconjg(a(j,k))
                               END IF
                               DO 290 i = 1,m
                                   b(i,j) = b(i,j) + temp*b(i,k)
   290                         CONTINUE
                           END IF
   300                 CONTINUE
                       temp = alpha
                       IF (nounit) THEN
                           IF (noconj) THEN
                               temp = temp*a(k,k)
                           ELSE
                               temp = temp*dconjg(a(k,k))
                           END IF
                       END IF
                       IF (temp.NE.one) THEN
                           DO 310 i = 1,m
                               b(i,k) = temp*b(i,k)
   310                     CONTINUE
                       END IF
   320             CONTINUE
               END IF
           END IF
       END IF
 !
       RETURN
 !
 !     End of ZTRMM .
 !
    END SUBROUTINE
    
 ! ------------------------------------------------------------------------------------     
    
       SUBROUTINE ztrti2( UPLO, DIAG, N, A, LDA, INFO )
 !
 !  -- LAPACK computational routine --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !
 !     .. Scalar Arguments ..
       CHARACTER          DIAG, UPLO
       INTEGER            INFO, LDA, N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16         A( LDA, * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16         ONE
       parameter( one = ( 1.0d+0, 0.0d+0 ) )
 !     ..
 !     .. Local Scalars ..
       LOGICAL            NOUNIT, UPPER
       INTEGER            J
       COMPLEX*16         AJJ
 !     ..
 !     .. External Functions ..
       LOGICAL            LSAME
       EXTERNAL           lsame
 !     ..
 !     .. External Subroutines ..
       EXTERNAL           xerbla, zscal, ztrmv
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          max
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input parameters.
 !
       info = 0
       upper = lsame( uplo, 'U' )
       nounit = lsame( diag, 'N' )
       IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
          info = -1
       ELSE IF( .NOT.nounit .AND. .NOT.lsame( diag, 'U' ) ) THEN
          info = -2
       ELSE IF( n.LT.0 ) THEN
          info = -3
       ELSE IF( lda.LT.max( 1, n ) ) THEN
          info = -5
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'ZTRTI2', -info )
          RETURN
       END IF
 !
       IF( upper ) THEN
 !
 !        Compute inverse of upper triangular matrix.
 !
          DO 10 j = 1, n
             IF( nounit ) THEN
                a( j, j ) = one / a( j, j )
                ajj = -a( j, j )
             ELSE
                ajj = -one
             END IF
 !
 !           Compute elements 1:j-1 of j-th column.
 !
             CALL ztrmv( 'Upper', 'No transpose', diag, j-1, a, lda,a( 1, j ), 1 )
             CALL zscal( j-1, ajj, a( 1, j ), 1 )
    10    CONTINUE
       ELSE
 !
 !        Compute inverse of lower triangular matrix.
 !
          DO 20 j = n, 1, -1
             IF( nounit ) THEN
                a( j, j ) = one / a( j, j )
                ajj = -a( j, j )
             ELSE
                ajj = -one
             END IF
             IF( j.LT.n ) THEN
 !
 !              Compute elements j+1:n of j-th column.
 !
                CALL ztrmv( 'Lower', 'No transpose', diag, n-j,a( j+1, j+1 ), lda, a( j+1, j ), 1 )
                CALL zscal( n-j, ajj, a( j+1, j ), 1 )
             END IF
    20    CONTINUE
       END IF
 !
       RETURN
 !
 !     End of ZTRTI2
 !
    END SUBROUTINE

! ------------------------------------------------------------------------------------     

       SUBROUTINE zscal(N,ZA,ZX,INCX)
 !
 !  -- Reference BLAS level1 routine --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !
 !     .. Scalar Arguments ..
       COMPLEX*16 ZA
       INTEGER INCX,N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 ZX(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Local Scalars ..
       INTEGER I,NINCX
 !     ..
       IF (n.LE.0 .OR. incx.LE.0) RETURN
       IF (incx.EQ.1) THEN
 !
 !        code for increment equal to 1
 !
          DO i = 1,n
             zx(i) = za*zx(i)
          END DO
       ELSE
 !
 !        code for increment not equal to 1
 !
          nincx = n*incx
          DO i = 1,nincx,incx
             zx(i) = za*zx(i)
          END DO
       END IF
       RETURN
    END SUBROUTINE

! ------------------------------------------------------------------------------------     

       SUBROUTINE ztrmv(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
 !
 !  -- Reference BLAS level2 routine --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !
 !     .. Scalar Arguments ..
       INTEGER INCX,LDA,N
       CHARACTER DIAG,TRANS,UPLO
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 A(LDA,*),X(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16 ZERO
       parameter(zero= (0.0d+0,0.0d+0))
 !     ..
 !     .. Local Scalars ..
       COMPLEX*16 TEMP
       INTEGER I,INFO,IX,J,JX,KX
       LOGICAL NOCONJ,NOUNIT
 !     ..
 !     .. External Functions ..
       LOGICAL LSAME
       EXTERNAL lsame
 !     ..
 !     .. External Subroutines ..
       EXTERNAL xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC dconjg,max
 !     ..
 !
 !     Test the input parameters.
 !
       info = 0
       IF (.NOT.lsame(uplo,'U') .AND. .NOT.lsame(uplo,'L')) THEN
           info = 1
       ELSE IF (.NOT.lsame(trans,'N') .AND. .NOT.lsame(trans,'T') .AND. .NOT.lsame(trans,'C')) THEN
           info = 2
       ELSE IF (.NOT.lsame(diag,'U') .AND. .NOT.lsame(diag,'N')) THEN
           info = 3
       ELSE IF (n.LT.0) THEN
           info = 4
       ELSE IF (lda.LT.max(1,n)) THEN
           info = 6
       ELSE IF (incx.EQ.0) THEN
           info = 8
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('ZTRMV ',info)
           RETURN
       END IF
 !
 !     Quick return if possible.
 !
       IF (n.EQ.0) RETURN
 !
       noconj = lsame(trans,'T')
       nounit = lsame(diag,'N')
 !
 !     Set up the start point in X if the increment is not unity. This
 !     will be  ( N - 1 )!INCX  too small for descending loops.
 !
       IF (incx.LE.0) THEN
           kx = 1 - (n-1)*incx
       ELSE IF (incx.NE.1) THEN
           kx = 1
       END IF
 !
 !     Start the operations. In this version the elements of A are
 !     accessed sequentially with one pass through A.
 !
       IF (lsame(trans,'N')) THEN
 !
 !        Form  x := A!x.
 !
           IF (lsame(uplo,'U')) THEN
               IF (incx.EQ.1) THEN
                   DO 20 j = 1,n
                       IF (x(j).NE.zero) THEN
                           temp = x(j)
                           DO 10 i = 1,j - 1
                               x(i) = x(i) + temp*a(i,j)
    10                     CONTINUE
                           IF (nounit) x(j) = x(j)*a(j,j)
                       END IF
    20             CONTINUE
               ELSE
                   jx = kx
                   DO 40 j = 1,n
                       IF (x(jx).NE.zero) THEN
                           temp = x(jx)
                           ix = kx
                           DO 30 i = 1,j - 1
                               x(ix) = x(ix) + temp*a(i,j)
                               ix = ix + incx
    30                     CONTINUE
                           IF (nounit) x(jx) = x(jx)*a(j,j)
                       END IF
                       jx = jx + incx
    40             CONTINUE
               END IF
           ELSE
               IF (incx.EQ.1) THEN
                   DO 60 j = n,1,-1
                       IF (x(j).NE.zero) THEN
                           temp = x(j)
                           DO 50 i = n,j + 1,-1
                               x(i) = x(i) + temp*a(i,j)
    50                     CONTINUE
                           IF (nounit) x(j) = x(j)*a(j,j)
                       END IF
    60             CONTINUE
               ELSE
                   kx = kx + (n-1)*incx
                   jx = kx
                   DO 80 j = n,1,-1
                       IF (x(jx).NE.zero) THEN
                           temp = x(jx)
                           ix = kx
                           DO 70 i = n,j + 1,-1
                               x(ix) = x(ix) + temp*a(i,j)
                               ix = ix - incx
    70                     CONTINUE
                           IF (nounit) x(jx) = x(jx)*a(j,j)
                       END IF
                       jx = jx - incx
    80             CONTINUE
               END IF
           END IF
       ELSE
 !
 !        Form  x := A!!T!x  or  x := A!!H!x.
 !
           IF (lsame(uplo,'U')) THEN
               IF (incx.EQ.1) THEN
                   DO 110 j = n,1,-1
                       temp = x(j)
                       IF (noconj) THEN
                           IF (nounit) temp = temp*a(j,j)
                           DO 90 i = j - 1,1,-1
                               temp = temp + a(i,j)*x(i)
    90                     CONTINUE
                       ELSE
                           IF (nounit) temp = temp*dconjg(a(j,j))
                           DO 100 i = j - 1,1,-1
                               temp = temp + dconjg(a(i,j))*x(i)
   100                     CONTINUE
                       END IF
                       x(j) = temp
   110             CONTINUE
               ELSE
                   jx = kx + (n-1)*incx
                   DO 140 j = n,1,-1
                       temp = x(jx)
                       ix = jx
                       IF (noconj) THEN
                           IF (nounit) temp = temp*a(j,j)
                           DO 120 i = j - 1,1,-1
                               ix = ix - incx
                               temp = temp + a(i,j)*x(ix)
   120                     CONTINUE
                       ELSE
                           IF (nounit) temp = temp*dconjg(a(j,j))
                           DO 130 i = j - 1,1,-1
                               ix = ix - incx
                               temp = temp + dconjg(a(i,j))*x(ix)
   130                     CONTINUE
                       END IF
                       x(jx) = temp
                       jx = jx - incx
   140             CONTINUE
               END IF
           ELSE
               IF (incx.EQ.1) THEN
                   DO 170 j = 1,n
                       temp = x(j)
                       IF (noconj) THEN
                           IF (nounit) temp = temp*a(j,j)
                           DO 150 i = j + 1,n
                               temp = temp + a(i,j)*x(i)
   150                     CONTINUE
                       ELSE
                           IF (nounit) temp = temp*dconjg(a(j,j))
                           DO 160 i = j + 1,n
                               temp = temp + dconjg(a(i,j))*x(i)
   160                     CONTINUE
                       END IF
                       x(j) = temp
   170             CONTINUE
               ELSE
                   jx = kx
                   DO 200 j = 1,n
                       temp = x(jx)
                       ix = jx
                       IF (noconj) THEN
                           IF (nounit) temp = temp*a(j,j)
                           DO 180 i = j + 1,n
                               ix = ix + incx
                               temp = temp + a(i,j)*x(ix)
   180                     CONTINUE
                       ELSE
                           IF (nounit) temp = temp*dconjg(a(j,j))
                           DO 190 i = j + 1,n
                               ix = ix + incx
                               temp = temp + dconjg(a(i,j))*x(ix)
   190                     CONTINUE
                       END IF
                       x(jx) = temp
                       jx = jx + incx
   200             CONTINUE
               END IF
           END IF
       END IF
 !
       RETURN
 !
 !     End of ZTRMV .
 !
    END SUBROUTINE
    
 ! ------------------------------------------------------------------------------------     
    
    subroutine zgetrf2	(	M,N,A,LDA,IPIV,INFO )	  !( m, n1, a, lda, ipiv, iinfo )	
    
 ! ZGETRF2 computes an LU factorization of a general M-by-N matrix A
 !using partial pivoting with row interchanges.
 !
 !The factorization has the form
 !   A = P * L * U
 !where P is a permutation matrix, L is lower triangular with unit
 !diagonal elements (lower trapezoidal if m > n), and U is upper
 !triangular (upper trapezoidal if m < n).
 !
 !This is the recursive version of the algorithm. It divides
 !the matrix into four submatrices:
 !
 !       [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
 !   A = [ -----|----- ]  with n1 = min(m,n)/2
 !       [  A21 | A22  ]       n2 = n-n1
 !
 !                                      [ A11 ]
 !The subroutine calls itself to factor [ --- ],
 !                                      [ A12 ]
 !                [ A12 ]
 !do the swaps on [ --- ], solve A12, update A22,
 !                [ A22 ]
 !
 !then calls itself to factor A22 and do the swaps on A21.   
 !
 !     .. Scalar Arguments ..
       INTEGER            INFO, LDA, M, N
 !     ..
 !     .. Array Arguments ..
       INTEGER            IPIV( * )
       COMPLEX*16         A( LDA, * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16         ONE, ZERO
       parameter( one = ( 1.0d+0, 0.0d+0 ), zero = ( 0.0d+0, 0.0d+0 ) )
 !     ..
 !     .. Local Scalars ..
       DOUBLE PRECISION   SFMIN
       COMPLEX*16         TEMP
       INTEGER            I, IINFO, N1, N2
 !     ..
 !     .. External Functions ..
       DOUBLE PRECISION   DLAMCH
       INTEGER            IZAMAX
       EXTERNAL           dlamch, izamax
 !     ..
 !     .. External Subroutines ..
       EXTERNAL           zgemm, zscal, zlaswp, ztrsm, xerbla,ZGETRF3
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          max, min
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input parameters
 !
       info = 0
       IF( m.LT.0 ) THEN
          info = -1
       ELSE IF( n.LT.0 ) THEN
          info = -2
       ELSE IF( lda.LT.max( 1, m ) ) THEN
          info = -4
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'ZGETRF2', -info )
          RETURN
       END IF
 !
 !     Quick return if possible
 !
       IF( m.EQ.0 .OR. n.EQ.0 ) RETURN
  
       IF ( m.EQ.1 ) THEN
 !
 !        Use unblocked code for one row case
 !        Just need to handle IPIV and INFO
 !
          ipiv( 1 ) = 1
          IF ( a(1,1).EQ.zero )info = 1
 !
       ELSE IF( n.EQ.1 ) THEN
 !
 !        Use unblocked code for one column case
 !
 !
 !        Compute machine safe minimum
 !
          sfmin = dlamch('S')
 !
 !        Find pivot and test for singularity
 !
          i = izamax( m, a( 1, 1 ), 1 )
          ipiv( 1 ) = i
          IF( a( i, 1 ).NE.zero ) THEN
 !
 !           Apply the interchange
 !
             IF( i.NE.1 ) THEN
                temp = a( 1, 1 )
                a( 1, 1 ) = a( i, 1 )
                a( i, 1 ) = temp
             END IF
 !
 !           Compute elements 2:M of the column
 !
             IF( abs(a( 1, 1 )) .GE. sfmin ) THEN
                CALL zscal( m-1, one / a( 1, 1 ), a( 2, 1 ), 1 )
             ELSE
                DO 10 i = 1, m-1
                   a( 1+i, 1 ) = a( 1+i, 1 ) / a( 1, 1 )
    10          CONTINUE
             END IF
 !
          ELSE
             info = 1
          END IF
  
       ELSE
 !
 !        Use recursive code
 !
          n1 = min( m, n ) / 2
          n2 = n-n1
 !
 !               [ A11 ]
 !        Factor [ --- ]
 !               [ A21 ]
 !
          CALL zgetrf3( m, n1, a, lda, ipiv, iinfo )
  
          IF ( info.EQ.0 .AND. iinfo.GT.0 )   info = iinfo
 !
 !                              [ A12 ]
 !        Apply interchanges to [ --- ]
 !                              [ A22 ]
 !
          CALL zlaswp( n2, a( 1, n1+1 ), lda, 1, n1, ipiv, 1 )
 !
 !        Solve A12
 !
          CALL ztrsm( 'L', 'L', 'N', 'U', n1, n2, one, a, lda,  a( 1, n1+1 ), lda )
 !
 !        Update A22
 !
          CALL zgemm( 'N', 'N', m-n1, n2, n1, -one, a( n1+1, 1 ), lda, a( 1, n1+1 ), lda, one, a( n1+1, n1+1 ), lda )
 !
 !        Factor A22
 !
          CALL zgetrf3( m-n1, n2, a( n1+1, n1+1 ), lda, ipiv( n1+1 ),   iinfo )
 !
 !        Adjust INFO and the pivot indices
 !
          IF ( info.EQ.0 .AND. iinfo.GT.0 ) info = iinfo + n1
          DO 20 i = n1+1, min( m, n )
             ipiv( i ) = ipiv( i ) + n1
    20    CONTINUE
 !
 !        Apply interchanges to A21
 !
          CALL zlaswp( n1, a( 1, 1 ), lda, n1+1, min( m, n), ipiv, 1 )
 !
       END IF
       RETURN
 !
 !     End of ZGETRF2
 !   
    End SUBROUTINE ZGETRF2

! ------------------------------------------------------------------------------------     
!!
subroutine zlaswp	(N, A, LDA, K1, K2, IPIV, INCX )	
    !  -- LAPACK auxiliary routine --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !
 !     .. Scalar Arguments ..
       INTEGER            INCX, K1, K2, LDA, N
 !     ..
 !     .. Array Arguments ..
       INTEGER            IPIV( * )
       COMPLEX*16         A( LDA, * )
 !     ..
 !
 ! =====================================================================
 !
 !     .. Local Scalars ..
       INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
       COMPLEX*16         TEMP
 !     ..
 !     .. Executable Statements ..
 !
 !     Interchange row I with row IPIV(K1+(I-K1)*abs(INCX)) for each of rows
 !     K1 through K2.
 !
       IF( incx.GT.0 ) THEN
          ix0 = k1
          i1 = k1
          i2 = k2
          inc = 1
       ELSE IF( incx.LT.0 ) THEN
          ix0 = k1 + ( k1-k2 )*incx
          i1 = k2
          i2 = k1
          inc = -1
       ELSE
          RETURN
       END IF
 !
       n32 = ( n / 32 )*32
       IF( n32.NE.0 ) THEN
          DO 30 j = 1, n32, 32
             ix = ix0
             DO 20 i = i1, i2, inc
                ip = ipiv( ix )
                IF( ip.NE.i ) THEN
                   DO 10 k = j, j + 31
                      temp = a( i, k )
                      a( i, k ) = a( ip, k )
                      a( ip, k ) = temp
    10             CONTINUE
                END IF
                ix = ix + incx
    20       CONTINUE
    30    CONTINUE
       END IF
       IF( n32.NE.n ) THEN
          n32 = n32 + 1
          ix = ix0
          DO 50 i = i1, i2, inc
             ip = ipiv( ix )
             IF( ip.NE.i ) THEN
                DO 40 k = n32, n
                   temp = a( i, k )
                   a( i, k ) = a( ip, k )
                   a( ip, k ) = temp
    40          CONTINUE
             END IF
             ix = ix + incx
    50    CONTINUE
       END IF
 !
       RETURN
 !
 !     End of ZLASWP
 
    End Subroutine
    !
 ! ------------------------------------------------------------------------------------     
    
integer function izamax	(N,ZX,INCX )
!     .. Scalar Arguments ..
       INTEGER INCX,N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 ZX(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Local Scalars ..
       DOUBLE PRECISION DMAX
       INTEGER I,IX
 !     ..
 !     .. External Functions ..
       DOUBLE PRECISION DCABS1
       EXTERNAL dcabs1
 !     ..
       izamax = 0
       IF (n.LT.1 .OR. incx.LE.0) RETURN
       izamax = 1
       IF (n.EQ.1) RETURN
       IF (incx.EQ.1) THEN
 !
 !        code for increment equal to 1
 !
          dmax = dcabs1(zx(1))
          DO i = 2,n
             IF (dcabs1(zx(i)).GT.dmax) THEN
                izamax = i
                dmax = dcabs1(zx(i))
             END IF
          END DO
       ELSE
 !
 !        code for increment not equal to 1
 !
          ix = 1
          dmax = dcabs1(zx(1))
          ix = ix + incx
          DO i = 2,n
             IF (dcabs1(zx(ix)).GT.dmax) THEN
                izamax = i
                dmax = dcabs1(zx(ix))
             END IF
             ix = ix + incx
          END DO
       END IF
       RETURN
    End Function IZAMAX
    
! ------------------------------------------------------------------------------------     
    
          DOUBLE PRECISION FUNCTION dcabs1(Z)
 !
 !  -- Reference BLAS level1 routine --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !
 !     .. Scalar Arguments ..
       COMPLEX*16 z
 !     ..
 !     ..
 !  =====================================================================
 !
 !     .. Intrinsic Functions ..
       INTRINSIC abs,dble,dimag
 !
       dcabs1 = abs(dble(z)) + abs(dimag(z))
       RETURN
    End Function
   
    subroutine zgetrf3	(	M,N,A,LDA,IPIV,INFO )	  !( m, n1, a, lda, ipiv, iinfo )	
       INTEGER            INFO, LDA, M, N
       INTEGER            IPIV( * )
       COMPLEX*16         A( LDA, * )
       COMPLEX*16         ONE, ZERO
       parameter( one = ( 1.0d+0, 0.0d+0 ), zero = ( 0.0d+0, 0.0d+0 ) )
       DOUBLE PRECISION   SFMIN
       COMPLEX*16         TEMP
       INTEGER            I, IINFO, N1, N2
       DOUBLE PRECISION   DLAMCH
       INTEGER            IZAMAX
       EXTERNAL           dlamch, izamax
       EXTERNAL           zgemm, zscal, zlaswp, ztrsm, xerbla
       INTRINSIC          max, min
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input parameters
 !
       info = 0
       IF( m.LT.0 ) THEN
          info = -1
       ELSE IF( n.LT.0 ) THEN
          info = -2
       ELSE IF( lda.LT.max( 1, m ) ) THEN
          info = -4
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'ZGETRF2', -info )
          RETURN
       END IF
 !
 !     Quick return if possible
 !
       IF( m.EQ.0 .OR. n.EQ.0 ) RETURN
  
       IF ( m.EQ.1 ) THEN
 !
 !        Use unblocked code for one row case
 !        Just need to handle IPIV and INFO
 !
          ipiv( 1 ) = 1
          IF ( a(1,1).EQ.zero )info = 1
 !
       ELSE IF( n.EQ.1 ) THEN
 !
 !        Use unblocked code for one column case
 !
 !
 !        Compute machine safe minimum
 !
          sfmin = dlamch('S')
 !
 !        Find pivot and test for singularity
 !
          i = izamax( m, a( 1, 1 ), 1 )
          ipiv( 1 ) = i
          IF( a( i, 1 ).NE.zero ) THEN
 !
 !           Apply the interchange
 !
             IF( i.NE.1 ) THEN
                temp = a( 1, 1 )
                a( 1, 1 ) = a( i, 1 )
                a( i, 1 ) = temp
             END IF
 !
 !           Compute elements 2:M of the column
 !
             IF( abs(a( 1, 1 )) .GE. sfmin ) THEN
                CALL zscal( m-1, one / a( 1, 1 ), a( 2, 1 ), 1 )
             ELSE
                DO 10 i = 1, m-1
                   a( 1+i, 1 ) = a( 1+i, 1 ) / a( 1, 1 )
    10          CONTINUE
             END IF
 !
          ELSE
             info = 1
          END IF
  
       ELSE
 !
 !        Use recursive code
 !
          n1 = min( m, n ) / 2
          n2 = n-n1
 !
 !               [ A11 ]
 !        Factor [ --- ]
 !               [ A21 ]
 !
          CALL zgetrf2( m, n1, a, lda, ipiv, iinfo )
  
          IF ( info.EQ.0 .AND. iinfo.GT.0 )   info = iinfo
 !
 !                              [ A12 ]
 !        Apply interchanges to [ --- ]
 !                              [ A22 ]
 !
          CALL zlaswp( n2, a( 1, n1+1 ), lda, 1, n1, ipiv, 1 )
 !
 !        Solve A12
 !
          CALL ztrsm( 'L', 'L', 'N', 'U', n1, n2, one, a, lda,  a( 1, n1+1 ), lda )
 !
 !        Update A22
 !
          CALL zgemm( 'N', 'N', m-n1, n2, n1, -one, a( n1+1, 1 ), lda, a( 1, n1+1 ), lda, one, a( n1+1, n1+1 ), lda )
 !
 !        Factor A22
 !
          CALL zgetrf2( m-n1, n2, a( n1+1, n1+1 ), lda, ipiv( n1+1 ),   iinfo )
 !
 !        Adjust INFO and the pivot indices
 !
          IF ( info.EQ.0 .AND. iinfo.GT.0 ) info = iinfo + n1
          DO 20 i = n1+1, min( m, n )
             ipiv( i ) = ipiv( i ) + n1
    20    CONTINUE
 !
 !        Apply interchanges to A21
 !
          CALL zlaswp( n1, a( 1, 1 ), lda, n1+1, min( m, n), ipiv, 1 )
 !
       END IF
       RETURN
 !
 !    
 !   
    End SUBROUTINE ZGETRF3
 ! ------------------------------------------------------------------------------------     
 ! ------------------------------------------------------------------------------------     
 ! ------------------------------------------------------------------------------------     
    
    
!       SUBROUTINE xerbla( SRNAME, INFO )
! !
! !  -- Reference BLAS level1 routine --
! !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
! !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
! !
! !     .. Scalar Arguments ..
!       CHARACTER*(*)      SRNAME
!       INTEGER            INFO
! !     ..
! !
! ! =====================================================================
! !
! !     .. Intrinsic Functions ..
!       INTRINSIC          len_trim
! !     ..
! !     .. Executable Statements ..
! !
!       WRITE( *, fmt = 9999 )srname( 1:len_trim( srname ) ), info
! !
!       stop
! !
!9999 FORMAT( ' ** On entry to ', a, ' parameter number ', i2, ' had ','an illegal value' )
! !
! !     End of XERBLA
! !
!       END SUBROUTINE

 !----------------------------------------------------     
!INTEGER  FUNCTION  ilaenv( ISPEC, NAME, OPTS, N1, N2, N3,N4 )
! !
! !  -- LAPACK test routine --
! !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
! !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
! !
! !     .. Scalar Arguments ..
!       CHARACTER*( * )    name, opts
!       INTEGER            ispec, n1, n2, n3, n4
! !     ..
! !
! !  =====================================================================
! !
! !     .. Intrinsic Functions ..
!       INTRINSIC          int, min, real
! !     ..
! !     .. External Functions ..
!       INTEGER            ieeeck
!       EXTERNAL           ieeeck
! !     ..
! !     .. Arrays in Common ..
!       INTEGER            iparms( 100 )
! !     ..
! !     .. Common blocks ..
!       COMMON             / claenv / iparms
! !     ..
! !     .. Save statement ..
!       SAVE               / claenv /
! !     ..
! !     .. Executable Statements ..
! !
!       IF( ispec.GE.1 .AND. ispec.LE.5 ) THEN
! !
! !        Return a value from the common block.
! !
!          IF ( name(2:6).EQ.'GEQR ' ) THEN
!             IF (n3.EQ.2) THEN
!                ilaenv = iparms( 2 )
!             ELSE
!                ilaenv = iparms( 1 )
!             END IF
!          ELSE IF ( name(2:6).EQ.'GELQ ' ) THEN
!             IF (n3.EQ.2) THEN
!                ilaenv = iparms( 2 )
!             ELSE
!                ilaenv = iparms( 1 )
!             END IF
!          ELSE
!             ilaenv = iparms( ispec )
!          END IF
! !
!       ELSE IF( ispec.EQ.6 ) THEN
! !
! !        Compute SVD crossover point.
! !
!          ilaenv = int( real( min( n1, n2 ) )*1.6e0 )
! !
!       ELSE IF( ispec.GE.7 .AND. ispec.LE.9 ) THEN
! !
! !        Return a value from the common block.
! !
!          ilaenv = iparms( ispec )
! !
!       ELSE IF( ispec.EQ.10 ) THEN
! !
! !        IEEE NaN arithmetic can be trusted not to trap
! !
! !        ILAENV = 0
!          ilaenv = 1
!          IF( ilaenv.EQ.1 ) THEN
!             ilaenv = ieeeck( 1, 0.0, 1.0 )
!          END IF
! !
!       ELSE IF( ispec.EQ.11 ) THEN
! !
! !        Infinity arithmetic can be trusted not to trap
! !
!!        ILAENV = 0
!          ilaenv = 1
!          IF( ilaenv.EQ.1 ) THEN
!             ilaenv = ieeeck( 0, 0.0, 1.0 )
!          END IF
! !
!       ELSE
! !
! !        Invalid value for ISPEC
! !
!          ilaenv = -1
!       END IF
! !
!       RETURN
! !
! !     End of ILAENV
! !
!END FUNCTION  ilaenv

!LOGICAL FUNCTION lsame(CA,CB)
! !
! !  -- Reference BLAS level1 routine --
! !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
! !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
! !
! !     .. Scalar Arguments ..
!       CHARACTER ca,cb
! !     ..
! !
! ! =====================================================================
! !
! !     .. Intrinsic Functions ..
!       INTRINSIC ichar
! !     ..
! !     .. Local Scalars ..
!       INTEGER inta,intb,zcode
! !     ..
! !
! !     Test if the characters are equal
! !
!       lsame = ca .EQ. cb
!       IF (lsame) RETURN
! !
! !     Now test for equivalence if both characters are alphabetic.
! !
!       zcode = ichar('Z')
! !
! !     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
! !     machines, on which ICHAR returns a value with bit 8 set.
! !     ICHAR('A') on Prime machines returns 193 which is the same as
! !     ICHAR('A') on an EBCDIC machine.
! !
!       inta = ichar(ca)
!       intb = ichar(cb)
! !
!       IF (zcode.EQ.90 .OR. zcode.EQ.122) THEN
! !
! !        ASCII is assumed - ZCODE is the ASCII code of either lower or
! !        upper case 'Z'.
! !
!           IF (inta.GE.97 .AND. inta.LE.122) inta = inta - 32
!           IF (intb.GE.97 .AND. intb.LE.122) intb = intb - 32
! !
!       ELSE IF (zcode.EQ.233 .OR. zcode.EQ.169) THEN
! !
! !        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
! !        upper case 'Z'.
! !
!           IF (inta.GE.129 .AND. inta.LE.137 .OR.inta.GE.145 .AND. inta.LE.153 .OR.inta.GE.162 .AND. inta.LE.169) inta = inta + 64
!           IF (intb.GE.129 .AND. intb.LE.137 .OR.intb.GE.145 .AND. intb.LE.153 .OR.intb.GE.162 .AND. intb.LE.169) intb = intb + 64
! !
!       ELSE IF (zcode.EQ.218 .OR. zcode.EQ.250) THEN
! !
! !        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
! !        plus 128 of either lower or upper case 'Z'.
! !
!           IF (inta.GE.225 .AND. inta.LE.250) inta = inta - 32
!           IF (intb.GE.225 .AND. intb.LE.250) intb = intb - 32
!       END IF
!       lsame = inta .EQ. intb
! !
! !     RETURN
! !
! !     End of LSAME
! !
!       END FUNCTION lsame
    
!INTEGER FUNCTION ieeeck( ISPEC, ZERO, ONE )
! !
! !  -- LAPACK auxiliary routine --
! !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
! !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
! !
! !     .. Scalar Arguments ..
!       INTEGER            ispec
!       REAL               one, zero
! !     ..
! !
! !  =====================================================================
! !
! !     .. Local Scalars ..
!       REAL               nan1, nan2, nan3, nan4, nan5, nan6, neginf,negzro, newzro, posinf
! !     ..
! !     .. Executable Statements ..
!       ieeeck = 1
! !
!       posinf = one / zero
!       IF( posinf.LE.one ) THEN
!          ieeeck = 0
!          RETURN
!       END IF
! !
!       neginf = -one / zero
!       IF( neginf.GE.zero ) THEN
!          ieeeck = 0
!          RETURN
!       END IF
! !
!       negzro = one / ( neginf+one )
!       IF( negzro.NE.zero ) THEN
!          ieeeck = 0
!          RETURN
!       END IF
! !
!       neginf = one / negzro
!       IF( neginf.GE.zero ) THEN
!          ieeeck = 0
!          RETURN
!       END IF
! !
!       newzro = negzro + zero
!       IF( newzro.NE.zero ) THEN
!          ieeeck = 0
!          RETURN
!       END IF
! !
!       posinf = one / newzro
!       IF( posinf.LE.one ) THEN
!          ieeeck = 0
!          RETURN
!       END IF
! !
!       neginf = neginf*posinf
!       IF( neginf.GE.zero ) THEN
!          ieeeck = 0
!          RETURN
!       END IF
! !
!       posinf = posinf*posinf
!       IF( posinf.LE.one ) THEN
!          ieeeck = 0
!          RETURN
!       END IF
! !
! !
! !
! !
! !     Return if we were only asked to check infinity arithmetic
! !
!       IF( ispec.EQ.0 )RETURN
! !
!       nan1 = posinf + neginf
! !
!       nan2 = posinf / neginf
! !
!       nan3 = posinf / posinf
! !
!       nan4 = posinf*zero
! !
!       nan5 = neginf*negzro
! !
!       nan6 = nan5*zero
! !
!       IF( nan1.EQ.nan1 ) THEN
!          ieeeck = 0
!          RETURN
!       END IF
! !
!       IF( nan2.EQ.nan2 ) THEN
!          ieeeck = 0
!          RETURN
!       END IF
! !
!       IF( nan3.EQ.nan3 ) THEN
!          ieeeck = 0
!          RETURN
!       END IF
! !
!       IF( nan4.EQ.nan4 ) THEN
!          ieeeck = 0
!          RETURN
!       END IF
! !
!       IF( nan5.EQ.nan5 ) THEN
!          ieeeck = 0
!          RETURN
!       END IF
! !
!       IF( nan6.EQ.nan6 ) THEN
!          ieeeck = 0
!          RETURN
!       END IF
! !
!       RETURN
!END FUNCTION ieeeck


    