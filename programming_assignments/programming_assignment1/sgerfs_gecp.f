*> \brief \b SGERFS_GECP
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download SGERFS + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgerfs.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgerfs.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgerfs.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGERFS_GECP( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, JPIV, B, LDB,
*                          X, LDX, FERR, BERR, WORK, IWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          TRANS
*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * ), JPIV( * ), IWORK( * )
*       REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
*      $                   BERR( * ), FERR( * ), WORK( * ), X( LDX, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGERFS_GECP improves the computed solution to a system of linear
*> equations and provides error bounds and backward error estimates for
*> the solution.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          Specifies the form of the system of equations:
*>          = 'N':  A * X = B     (No transpose)
*>          = 'T':  A**T * X = B  (Transpose)
*>          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of columns
*>          of the matrices B and X.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          The original N-by-N matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] AF
*> \verbatim
*>          AF is REAL array, dimension (LDAF,N)
*>          The factors L and U from the factorization A = P*L*U*Q
*>          as computed by SGETC2.
*> \endverbatim
*>
*> \param[in] LDAF
*> \verbatim
*>          LDAF is INTEGER
*>          The leading dimension of the array AF.  LDAF >= max(1,N).
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          The pivot indices from SGETC2; for 1<=i<=N, row i of the
*>          matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[in] JPIV
*> \verbatim
*>          JPIV is INTEGER array, dimension (N)
*>          The pivot indices from SGETC2; for 1<=i<=N, column j of the
*>          matrix was interchanged with column JPIV(i).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is REAL array, dimension (LDB,NRHS)
*>          The right hand side matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is REAL array, dimension (LDX,NRHS)
*>          On entry, the solution matrix X, as computed by SGETRS_GECP.
*>          On exit, the improved solution matrix X.
*> \endverbatim
*>
*> \param[in] LDX
*> \verbatim
*>          LDX is INTEGER
*>          The leading dimension of the array X.  LDX >= max(1,N).
*> \endverbatim
*>
*> \param[out] FERR
*> \verbatim
*>          FERR is REAL array, dimension (NRHS)
*>          The estimated forward error bound for each solution vector
*>          X(j) (the j-th column of the solution matrix X).
*>          If XTRUE is the true solution corresponding to X(j), FERR(j)
*>          is an estimated upper bound for the magnitude of the largest
*>          element in (X(j) - XTRUE) divided by the magnitude of the
*>          largest element in X(j).  The estimate is as reliable as
*>          the estimate for RCOND, and is almost always a slight
*>          overestimate of the true error.
*> \endverbatim
*>
*> \param[out] BERR
*> \verbatim
*>          BERR is REAL array, dimension (NRHS)
*>          The componentwise relative backward error of each solution
*>          vector X(j) (i.e., the smallest relative change in
*>          any element of A or B that makes X(j) an exact solution).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (3*N)
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*> \par Internal Parameters:
*  =========================
*>
*> \verbatim
*>  ITMAX is the maximum number of steps of iterative refinement.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup gerfs
*
*  =====================================================================

      SUBROUTINE sgerfs_gecp( TRANS, N, NRHS, A, LDA, AF, LDAF,
     $                   IPIV, JPIV, B, LDB, X, LDX, FERR, BERR,
     $                   WORK, IWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), JPIV( * ), IWORK( * )
      REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
     $                   berr( * ), ferr( * ), work( * ), x( ldx, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      REAL               ZERO
      parameter( zero = 0.0e+0 )
      REAL               ONE
      parameter( one = 1.0e+0 )
      REAL               TWO
      parameter( two = 2.0e+0 )
      REAL               THREE
      parameter( three = 3.0e+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
      CHARACTER          TRANST
      INTEGER            COUNT, I, J, K, KASE, NZ
      REAL               EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           saxpy, scopy, sgemv, sgetrs_gecp, slacn2,
     $                   xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH
      EXTERNAL           lsame, slamch
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      info = 0
      notran = lsame( trans, 'N' )
      IF( .NOT.notran .AND. .NOT.lsame( trans, 'T' ) .AND. .NOT.
     $    lsame( trans, 'C' ) ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( nrhs.LT.0 ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( ldaf.LT.max( 1, n ) ) THEN
         info = -7
      ELSE IF( ldb.LT.max( 1, n ) ) THEN
         info = -10
      ELSE IF( ldx.LT.max( 1, n ) ) THEN
         info = -12
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'SGERFS', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 .OR. nrhs.EQ.0 ) THEN
         DO 10 j = 1, nrhs
            ferr( j ) = zero
            berr( j ) = zero
   10    CONTINUE
         RETURN
      END IF
*
      IF( notran ) THEN
         transt = 'T'
      ELSE
         transt = 'N'
      END IF
*
*     NZ = maximum number of nonzero elements in each row of A, plus 1
*
      nz = n + 1
      eps = slamch( 'Epsilon' )
      safmin = slamch( 'Safe minimum' )
      safe1 = real( nz )*safmin
      safe2 = safe1 / eps
*
*     Do for each right hand side
*
      DO 140 j = 1, nrhs
*
         count = 1
         lstres = three
   20    CONTINUE
*
*        Loop until stopping criterion is satisfied.
*
*        Compute residual R = B - op(A) * X,
*        where op(A) = A, A**T, or A**H, depending on TRANS.
*
         CALL scopy( n, b( 1, j ), 1, work( n+1 ), 1 )
         CALL sgemv( trans, n, n, -one, a, lda, x( 1, j ), 1, one,
     $               work( n+1 ), 1 )
*
*        Compute componentwise relative backward error from formula
*
*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )
*
*        where abs(Z) is the componentwise absolute value of the matrix
*        or vector Z.  If the i-th component of the denominator is less
*        than SAFE2, then SAFE1 is added to the i-th components of the
*        numerator and denominator before dividing.
*
         DO 30 i = 1, n
            work( i ) = abs( b( i, j ) )
   30    CONTINUE
*
*        Compute abs(op(A))*abs(X) + abs(B).
*
         IF( notran ) THEN
            DO 50 k = 1, n
               xk = abs( x( k, j ) )
               DO 40 i = 1, n
                  work( i ) = work( i ) + abs( a( i, k ) )*xk
   40          CONTINUE
   50       CONTINUE
         ELSE
            DO 70 k = 1, n
               s = zero
               DO 60 i = 1, n
                  s = s + abs( a( i, k ) )*abs( x( i, j ) )
   60          CONTINUE
               work( k ) = work( k ) + s
   70       CONTINUE
         END IF
         s = zero
         DO 80 i = 1, n
            IF( work( i ).GT.safe2 ) THEN
               s = max( s, abs( work( n+i ) ) / work( i ) )
            ELSE
               s = max( s, ( abs( work( n+i ) )+safe1 ) /
     $             ( work( i )+safe1 ) )
            END IF
   80    CONTINUE
         berr( j ) = s
*
*        Test stopping criterion. Continue iterating if
*           1) The residual BERR(J) is larger than machine epsilon, and
*           2) BERR(J) decreased by at least a factor of 2 during the
*              last iteration, and
*           3) At most ITMAX iterations tried.
*
         IF( berr( j ).GT.eps .AND. two*berr( j ).LE.lstres .AND.
     $       count.LE.itmax ) THEN
*
*           Update solution and try again.
*
            CALL sgetrs_gecp( trans, n, 1, af, ldaf, ipiv, jpiv,
     $                   work( n+1 ), n, info )
            CALL saxpy( n, one, work( n+1 ), 1, x( 1, j ), 1 )
            lstres = berr( j )
            count = count + 1
            GO TO 20
         END IF
*
*        Bound error from formula
*
*        norm(X - XTRUE) / norm(X) .le. FERR =
*        norm( abs(inv(op(A)))*
*           ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X)
*
*        where
*          norm(Z) is the magnitude of the largest component of Z
*          inv(op(A)) is the inverse of op(A)
*          abs(Z) is the componentwise absolute value of the matrix or
*             vector Z
*          NZ is the maximum number of nonzeros in any row of A, plus 1
*          EPS is machine epsilon
*
*        The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B))
*        is incremented by SAFE1 if the i-th component of
*        abs(op(A))*abs(X) + abs(B) is less than SAFE2.
*
*        Use SLACN2 to estimate the infinity-norm of the matrix
*           inv(op(A)) * diag(W),
*        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))
*
         DO 90 i = 1, n
            IF( work( i ).GT.safe2 ) THEN
               work( i ) = abs( work( n+i ) ) + real( nz )*eps*work( i )
            ELSE
               work( i ) = abs( work( n+i ) ) + real( nz )*eps*work( i )
     $                     + safe1
            END IF
   90    CONTINUE
*
         kase = 0
  100    CONTINUE
         CALL slacn2( n, work( 2*n+1 ), work( n+1 ), iwork,
     $                ferr( j ),
     $                kase, isave )
         IF( kase.NE.0 ) THEN
            IF( kase.EQ.1 ) THEN
*
*              Multiply by diag(W)*inv(op(A)**T).
*
               CALL sgetrs_gecp( transt, n, 1, af, ldaf, ipiv, jpiv,
     $                      work( n+1 ), n, info )
               DO 110 i = 1, n
                  work( n+i ) = work( i )*work( n+i )
  110          CONTINUE
            ELSE
*
*              Multiply by inv(op(A))*diag(W).
*
               DO 120 i = 1, n
                  work( n+i ) = work( i )*work( n+i )
  120          CONTINUE
               CALL sgetrs_gecp( trans, n, 1, af, ldaf, ipiv, jpiv,
     $                      work( n+1 ), n, info )
            END IF
            GO TO 100
         END IF
*
*        Normalize error.
*
         lstres = zero
         DO 130 i = 1, n
            lstres = max( lstres, abs( x( i, j ) ) )
  130    CONTINUE
         IF( lstres.NE.zero )
     $      ferr( j ) = ferr( j ) / lstres
*
  140 CONTINUE
*
      RETURN
*
*     End of SGERFS_GECP
*

      END
