      SUBROUTINE sgecon_gecp( NORM, N, A, LDA, ANORM, RCOND, WORK,
     $                        IWORK, IPIV, JPIV, INFO )
*
*  -- LAPACK variant for complete pivoting --
*  Based on LAPACK routine SGECON but modified to handle complete pivoting
*
*  Purpose
*  =======
*
*  SGECON_GECP estimates the reciprocal of the condition number of a
*  general real matrix A with complete pivoting LU factorization, in
*  either the 1-norm or the infinity-norm, using the LU factorization
*  computed by SGETC2.
*
*  An estimate is obtained for norm(inv(A)), and the reciprocal of the
*  condition number is computed as RCOND = 1 / (norm(A) * norm(inv(A))).
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies whether the 1-norm condition number or the
*          infinity-norm condition number is required:
*          = '1' or 'O':  1-norm;
*          = 'I':         Infinity-norm.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input) REAL array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U*Q
*          as computed by SGETC2.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  ANORM   (input) REAL
*          If NORM = '1' or 'O', the 1-norm of the original matrix A.
*          If NORM = 'I', the infinity-norm of the original matrix A.
*
*  RCOND   (output) REAL
*          The reciprocal of the condition number of the matrix A,
*          computed as RCOND = 1/(norm(A) * norm(inv(A))).
*
*  WORK    (workspace) REAL array, dimension (4*N)
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The row pivot indices from SGETC2.
*
*  JPIV    (input) INTEGER array, dimension (N)
*          The column pivot indices from SGETC2.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            N, LDA, INFO
      REAL               ANORM, RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * ), IPIV( * ), JPIV( * )
      REAL               A( LDA, * ), WORK( * )
*     ..
*     .. Local Scalars ..
      LOGICAL            ONENRM
      CHARACTER          NORMIN
      INTEGER            KASE, KASE1
      REAL               AINVNM, SCALE, SL, SMLNUM, SU
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH
      EXTERNAL           LSAME, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLACN2, SGETRS_GECP, XERBLA, SSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      ONENRM = NORM.EQ.'1' .OR. LSAME( NORM, 'O' )
      IF( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGECON_GECP', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.EQ.ZERO ) THEN
         RETURN
      END IF
*
      SMLNUM = SLAMCH( 'Safe minimum' )
*
*     Estimate the norm of inv(A) using complete pivoting LU factors.
*
      AINVNM = ZERO
      NORMIN = 'N'
      IF( ONENRM ) THEN
         KASE1 = 1
      ELSE
         KASE1 = 2
      END IF
      KASE = 0
   10 CONTINUE
      CALL SLACN2( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.KASE1 ) THEN
*
*           Multiply by inv(A): solve A * x = v
*           SGETRS_GECP handles all permutations internally
            CALL SGETRS_GECP( 'No transpose', N, 1, A, LDA, IPIV, JPIV,
     $                        WORK, N, INFO )
         ELSE
*
*           Multiply by inv(A^T): solve A^T * x = v
*           SGETRS_GECP handles all permutations internally
            CALL SGETRS_GECP( 'Transpose', N, 1, A, LDA, IPIV, JPIV,
     $                        WORK, N, INFO )
         END IF
*
         GO TO 10
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM.NE.ZERO )
     $   RCOND = ( ONE / AINVNM ) / ANORM
*
      RETURN
      END
*
*     =====================================================================
*     Helper subroutine: Apply column permutation Q to a vector
*     =====================================================================
*
      SUBROUTINE slapiv_col( N, X, JPIV )
      INTEGER            N
      INTEGER            JPIV( N )
      REAL               X( N )
*
*     Apply column permutation: X_out(i) = X_in(JPIV(i))
*
      INTEGER            I, J
      REAL               TEMP( 1000 )
*
      DO 10 I = 1, N
         TEMP( I ) = X( I )
   10 CONTINUE
*
      DO 20 I = 1, N
         X( I ) = TEMP( JPIV( I ) )
   20 CONTINUE
*
      RETURN
      END
*
*     =====================================================================
*     Helper subroutine: Apply inverse column permutation Q^T to a vector
*     =====================================================================
*
      SUBROUTINE slapiv_col_t( N, X, JPIV )
      INTEGER            N
      INTEGER            JPIV( N )
      REAL               X( N )
*
*     Apply inverse column permutation: X_out(JPIV(i)) = X_in(i)
*
      INTEGER            I
      REAL               TEMP( 1000 )
*
      DO 10 I = 1, N
         TEMP( I ) = X( I )
   10 CONTINUE
*
      DO 20 I = 1, N
         X( JPIV( I ) ) = TEMP( I )
   20 CONTINUE
*
      RETURN
      END
