*> \brief <b> SGESVX_GECP computes the solution to system of linear equations A * X = B for GE matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download SGESVX + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgesvx.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgesvx.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgesvx.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGESVX_GECP( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, JPIV
*                          EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR,
*                          WORK, IWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          EQUED, FACT, TRANS
*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
*       REAL               RCOND
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * ), JPIV( * ), IWORK( * )
*       REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
*      $                   BERR( * ), C( * ), FERR( * ), R( * ),
*      $                   WORK( * ), X( LDX, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGESVX_GECP uses the LU factorization to compute the solution to a real
*> system of linear equations
*>    A * X = B,
*> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
*>
*> Error bounds on the solution and a condition estimate are also
*> provided.
*> \endverbatim
*
*> \par Description:
*  =================
*>
*> \verbatim
*>
*> The following steps are performed:
*>
*> 1. If FACT = 'E', real scaling factors are computed to equilibrate
*>    the system:
*>       TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B
*>       TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B
*>       TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B
*>    Whether or not the system will be equilibrated depends on the
*>    scaling of the matrix A, but if equilibration is used, A is
*>    overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')
*>    or diag(C)*B (if TRANS = 'T' or 'C').
*>
*> 2. If FACT = 'N' or 'E', the LU decomposition is used to factor the
*>    matrix A (after equilibration if FACT = 'E') as
*>       A = P * L * U,
*>    where P is a permutation matrix, L is a unit lower triangular
*>    matrix, and U is upper triangular.
*>
*> 3. If some U(i,i)=0, so that U is exactly singular, then the routine
*>    returns with INFO = i. Otherwise, the factored form of A is used
*>    to estimate the condition number of the matrix A.  If the
*>    reciprocal of the condition number is less than machine precision,
*>    INFO = N+1 is returned as a warning, but the routine still goes on
*>    to solve for X and compute error bounds as described below.
*>
*> 4. The system of equations is solved for X using the factored form
*>    of A.
*>
*> 5. Iterative refinement is applied to improve the computed solution
*>    matrix and calculate error bounds and backward error estimates
*>    for it.
*>
*> 6. If equilibration was used, the matrix X is premultiplied by
*>    diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so
*>    that it solves the original system before equilibration.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] FACT
*> \verbatim
*>          FACT is CHARACTER*1
*>          Specifies whether or not the factored form of the matrix A is
*>          supplied on entry, and if not, whether the matrix A should be
*>          equilibrated before it is factored.
*>          = 'F':  On entry, AF, IPIV and JPIV contain the factored form of A.
*>                  If EQUED is not 'N', the matrix A has been
*>                  equilibrated with scaling factors given by R and C.
*>                  A, AF, IPIV and JPIV are not modified.
*>          = 'N':  The matrix A will be copied to AF and factored.
*>          = 'E':  The matrix A will be equilibrated if necessary, then
*>                  copied to AF and factored.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          Specifies the form of the system of equations:
*>          = 'N':  A * X = B     (No transpose)
*>          = 'T':  A**T * X = B  (Transpose)
*>          = 'C':  A**H * X = B  (Transpose)
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of linear equations, i.e., the order of the
*>          matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of columns
*>          of the matrices B and X.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is
*>          not 'N', then A must have been equilibrated by the scaling
*>          factors in R and/or C.  A is not modified if FACT = 'F' or
*>          'N', or if FACT = 'E' and EQUED = 'N' on exit.
*>
*>          On exit, if EQUED .ne. 'N', A is scaled as follows:
*>          EQUED = 'R':  A := diag(R) * A
*>          EQUED = 'C':  A := A * diag(C)
*>          EQUED = 'B':  A := diag(R) * A * diag(C).
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] AF
*> \verbatim
*>          AF is REAL array, dimension (LDAF,N)
*>          If FACT = 'F', then AF is an input argument and on entry
*>          contains the factors L and U from the factorization
*>          A = P*L*U*Q as computed by SGETC2.  If EQUED .ne. 'N', then
*>          AF is the factored form of the equilibrated matrix A.
*>
*>          If FACT = 'N', then AF is an output argument and on exit
*>          returns the factors L and U from the factorization A = P*L*U
*>          of the original matrix A.
*>
*>          If FACT = 'E', then AF is an output argument and on exit
*>          returns the factors L and U from the factorization A = P*L*U
*>          of the equilibrated matrix A (see the description of A for
*>          the form of the equilibrated matrix).
*> \endverbatim
*>
*> \param[in] LDAF
*> \verbatim
*>          LDAF is INTEGER
*>          The leading dimension of the array AF.  LDAF >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          If FACT = 'F', then IPIV is an input argument and on entry
*>          contains the pivot indices from the factorization A = P*L*U*Q
*>          as computed by SGETC2; row i of the matrix was interchanged
*>          with row IPIV(i).
*>
*>          If FACT = 'N', then IPIV is an output argument and on exit
*>          contains the pivot indices from the factorization A = P*L*U*Q
*>          of the original matrix A.
*>
*>          If FACT = 'E', then IPIV is an output argument and on exit
*>          contains the pivot indices from the factorization A = P*L*U*Q
*>          of the equilibrated matrix A.
*> \endverbatim
*>
*> \param[in,out] JPIV
*> \verbatim
*>          JPIV is INTEGER array, dimension (N)
*>          If FACT = 'F', then JPIV is an input argument and on entry
*>          contains the pivot indices from the factorization A = P*L*U*Q
*>          as computed by SGETC2; column i of the matrix was interchanged
*>          with column JPIV(i).
*>
*>          If FACT = 'N', then JPIV is an output argument and on exit
*>          contains the pivot indices from the factorization A = P*L*U*Q
*>          of the original matrix A.
*>
*>          If FACT = 'E', then JPIV is an output argument and on exit
*>          contains the pivot indices from the factorization A = P*L*U*Q
*>          of the equilibrated matrix A.
*> \endverbatim
*>
*> \param[in,out] EQUED
*> \verbatim
*>          EQUED is CHARACTER*1
*>          Specifies the form of equilibration that was done.
*>          = 'N':  No equilibration (always true if FACT = 'N').
*>          = 'R':  Row equilibration, i.e., A has been premultiplied by
*>                  diag(R).
*>          = 'C':  Column equilibration, i.e., A has been postmultiplied
*>                  by diag(C).
*>          = 'B':  Both row and column equilibration, i.e., A has been
*>                  replaced by diag(R) * A * diag(C).
*>          EQUED is an input argument if FACT = 'F'; otherwise, it is an
*>          output argument.
*> \endverbatim
*>
*> \param[in,out] R
*> \verbatim
*>          R is REAL array, dimension (N)
*>          The row scale factors for A.  If EQUED = 'R' or 'B', A is
*>          multiplied on the left by diag(R); if EQUED = 'N' or 'C', R
*>          is not accessed.  R is an input argument if FACT = 'F';
*>          otherwise, R is an output argument.  If FACT = 'F' and
*>          EQUED = 'R' or 'B', each element of R must be positive.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is REAL array, dimension (N)
*>          The column scale factors for A.  If EQUED = 'C' or 'B', A is
*>          multiplied on the right by diag(C); if EQUED = 'N' or 'R', C
*>          is not accessed.  C is an input argument if FACT = 'F';
*>          otherwise, C is an output argument.  If FACT = 'F' and
*>          EQUED = 'C' or 'B', each element of C must be positive.
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is REAL array, dimension (LDB,NRHS)
*>          On entry, the N-by-NRHS right hand side matrix B.
*>          On exit,
*>          if EQUED = 'N', B is not modified;
*>          if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by
*>          diag(R)*B;
*>          if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is
*>          overwritten by diag(C)*B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is REAL array, dimension (LDX,NRHS)
*>          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X
*>          to the original system of equations.  Note that A and B are
*>          modified on exit if EQUED .ne. 'N', and the solution to the
*>          equilibrated system is inv(diag(C))*X if TRANS = 'N' and
*>          EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C'
*>          and EQUED = 'R' or 'B'.
*> \endverbatim
*>
*> \param[in] LDX
*> \verbatim
*>          LDX is INTEGER
*>          The leading dimension of the array X.  LDX >= max(1,N).
*> \endverbatim
*>
*> \param[out] RCOND
*> \verbatim
*>          RCOND is REAL
*>          The estimate of the reciprocal condition number of the matrix
*>          A after equilibration (if done).  If RCOND is less than the
*>          machine precision (in particular, if RCOND = 0), the matrix
*>          is singular to working precision.  This condition is
*>          indicated by a return code of INFO > 0.
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
*>          WORK is REAL array, dimension (MAX(1,4*N))
*>          On exit, WORK(1) contains the reciprocal pivot growth
*>          factor norm(A)/norm(U). The "max absolute element" norm is
*>          used. If WORK(1) is much less than 1, then the stability
*>          of the LU factorization of the (equilibrated) matrix A
*>          could be poor. This also means that the solution X, condition
*>          estimator RCOND, and forward error bound FERR could be
*>          unreliable. If factorization fails with 0<INFO<=N, then
*>          WORK(1) contains the reciprocal pivot growth factor for the
*>          leading INFO columns of A.
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
*>          > 0:  if INFO = i, and i is
*>                <= N:  U(i,i) is exactly zero.  The factorization has
*>                       been completed, but the factor U is exactly
*>                       singular, so the solution and error bounds
*>                       could not be computed. RCOND = 0 is returned.
*>                = N+1: U is nonsingular, but RCOND is less than machine
*>                       precision, meaning that the matrix is singular
*>                       to working precision.  Nevertheless, the
*>                       solution and error bounds are computed because
*>                       there are a number of situations where the
*>                       computed solution can be more accurate than the
*>                       value of RCOND would suggest.
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
*> \ingroup gesvx
*
*  =====================================================================

      SUBROUTINE sgesvx_gecp( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF,
     $                   IPIV, JPIV,
     $                   EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR,
     $                   WORK, IWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED, FACT, TRANS
      INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
      REAL               RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), JPIV( * ), IWORK( * )
      REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
     $                   BERR( * ), C( * ), FERR( * ), R( * ),
     $                   work( * ), x( ldx, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0e+0, one = 1.0e+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU
      CHARACTER          NORM
      INTEGER            I, INFEQU, J
      REAL               AMAX, ANORM, BIGNUM, COLCND, RCMAX, RCMIN,
     $                   rowcnd, rpvgrw, smlnum
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH, SLANGE, SLANTR
      EXTERNAL           LSAME, SLAMCH, SLANGE, SLANTR
*     ..
*     .. External Subroutines ..
      EXTERNAL           sgecon_gecp, sgeequ, sgerfs_gecp, sgetc2,
     $                   sgetrs_gecp, slacpy, slaqge, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
      info = 0
      nofact = lsame( fact, 'N' )
      equil = lsame( fact, 'E' )
      notran = lsame( trans, 'N' )
      IF( nofact .OR. equil ) THEN
         equed = 'N'
         rowequ = .false.
         colequ = .false.
      ELSE
         rowequ = lsame( equed, 'R' ) .OR. lsame( equed, 'B' )
         colequ = lsame( equed, 'C' ) .OR. lsame( equed, 'B' )
         smlnum = slamch( 'Safe minimum' )
         bignum = one / smlnum
      END IF
*
*     Test the input parameters.
*
      IF( .NOT.nofact .AND.
     $    .NOT.equil .AND.
     $    .NOT.lsame( fact, 'F' ) )
     $     THEN
         info = -1
      ELSE IF( .NOT.notran .AND. .NOT.lsame( trans, 'T' ) .AND. .NOT.
     $         lsame( trans, 'C' ) ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( nrhs.LT.0 ) THEN
         info = -4
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -6
      ELSE IF( ldaf.LT.max( 1, n ) ) THEN
         info = -8
      ELSE IF( lsame( fact, 'F' ) .AND. .NOT.
     $         ( rowequ .OR. colequ .OR. lsame( equed, 'N' ) ) ) THEN
         info = -10
      ELSE
         IF( rowequ ) THEN
            rcmin = bignum
            rcmax = zero
            DO 10 j = 1, n
               rcmin = min( rcmin, r( j ) )
               rcmax = max( rcmax, r( j ) )
   10       CONTINUE
            IF( rcmin.LE.zero ) THEN
               info = -11
            ELSE IF( n.GT.0 ) THEN
               rowcnd = max( rcmin, smlnum ) / min( rcmax, bignum )
            ELSE
               rowcnd = one
            END IF
         END IF
         IF( colequ .AND. info.EQ.0 ) THEN
            rcmin = bignum
            rcmax = zero
            DO 20 j = 1, n
               rcmin = min( rcmin, c( j ) )
               rcmax = max( rcmax, c( j ) )
   20       CONTINUE
            IF( rcmin.LE.zero ) THEN
               info = -12
            ELSE IF( n.GT.0 ) THEN
               colcnd = max( rcmin, smlnum ) / min( rcmax, bignum )
            ELSE
               colcnd = one
            END IF
         END IF
         IF( info.EQ.0 ) THEN
            IF( ldb.LT.max( 1, n ) ) THEN
               info = -14
            ELSE IF( ldx.LT.max( 1, n ) ) THEN
               info = -16
            END IF
         END IF
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'SGESVX', -info )
         RETURN
      END IF
*
      IF( equil ) THEN
*
*        Compute row and column scalings to equilibrate the matrix A.
*
         CALL sgeequ( n, n, a, lda, r, c, rowcnd, colcnd, amax,
     $                infequ )
         IF( infequ.EQ.0 ) THEN
*
*           Equilibrate the matrix.
*
            CALL slaqge( n, n, a, lda, r, c, rowcnd, colcnd, amax,
     $                   equed )
            rowequ = lsame( equed, 'R' ) .OR. lsame( equed, 'B' )
            colequ = lsame( equed, 'C' ) .OR. lsame( equed, 'B' )
         END IF
      END IF
*
*     Scale the right hand side.
*
      IF( notran ) THEN
         IF( rowequ ) THEN
            DO 40 j = 1, nrhs
               DO 30 i = 1, n
                  b( i, j ) = r( i )*b( i, j )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( colequ ) THEN
         DO 60 j = 1, nrhs
            DO 50 i = 1, n
               b( i, j ) = c( i )*b( i, j )
   50       CONTINUE
   60    CONTINUE
      END IF
*
      IF( nofact .OR. equil ) THEN
*
*        Compute the LU factorization of A.
*
         CALL slacpy( 'Full', n, n, a, lda, af, ldaf )
         CALL sgetc2( n, af, ldaf, ipiv, jpiv, info )
*
*        Note: SGETC2 returns INFO > 0 if U(INFO,INFO) is small and may
*        cause overflow. Unlike SGETRF, this is not a hard failure.
*        We continue and let SGECON determine if the matrix is singular.
*        Reset INFO to 0 to avoid interfering with later error checks.
*
         info = 0
      END IF
*
*     Compute the norm of the matrix A and the
*     reciprocal pivot growth factor RPVGRW.
*
      IF( notran ) THEN
         norm = '1'
      ELSE
         norm = 'I'
      END IF
      anorm = slange( norm, n, n, a, lda, work )
      rpvgrw = slantr( 'M', 'U', 'N', n, n, af, ldaf, work )
      IF( rpvgrw.EQ.zero ) THEN
         rpvgrw = one
      ELSE
         rpvgrw = slange( 'M', n, n, a, lda, work ) / rpvgrw
      END IF
*
*     Compute the reciprocal of the condition number of A.
*
      CALL sgecon_gecp( norm, n, af, ldaf, anorm, rcond, work,
     $                  iwork, ipiv, jpiv, info )
*
*     Compute the solution matrix X.
*
      CALL slacpy( 'Full', n, nrhs, b, ldb, x, ldx )
      CALL sgetrs_gecp( trans, n, nrhs, af, ldaf, ipiv, jpiv, x,
     $                  ldx, info )
*
*     Use iterative refinement to improve the computed solution and
*     compute error bounds and backward error estimates for it.
*
      CALL sgerfs_gecp( trans, n, nrhs, a, lda, af, ldaf, ipiv,
     $             jpiv, b, ldb, x, ldx, ferr, berr, work, iwork,
     $             info )
*
*     Transform the solution matrix X to a solution of the original
*     system.
*
      IF( notran ) THEN
         IF( colequ ) THEN
            DO 80 j = 1, nrhs
               DO 70 i = 1, n
                  x( i, j ) = c( i )*x( i, j )
   70          CONTINUE
   80       CONTINUE
            DO 90 j = 1, nrhs
               ferr( j ) = ferr( j ) / colcnd
   90       CONTINUE
         END IF
      ELSE IF( rowequ ) THEN
         DO 110 j = 1, nrhs
            DO 100 i = 1, n
               x( i, j ) = r( i )*x( i, j )
  100       CONTINUE
  110    CONTINUE
         DO 120 j = 1, nrhs
            ferr( j ) = ferr( j ) / rowcnd
  120    CONTINUE
      END IF
*
*     Set INFO = N+1 if the matrix is singular to working precision.
*
      IF( rcond.LT.slamch( 'Epsilon' ) )
     $   info = n + 1
*
      work( 1 ) = rpvgrw
      RETURN
*
*     End of SGESVX_GECP
*

      END
