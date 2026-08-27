*> \brief \b SGETRS_GECP
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download SGETRS + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgetrs.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgetrs.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgetrs.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGETRS_GECP( TRANS, N, NRHS, A, LDA, IPIV, JPIV, B, LDB, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          TRANS
*       INTEGER            INFO, LDA, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * ), JPIV( * )
*       REAL               A( LDA, * ), B( LDB, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGETRS_GECP solves a system of linear equations
*>    A * X = B  or  A**T * X = B
*> with a general N-by-N matrix A using the LU factorization computed
*> by SGETC2.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          Specifies the form of the system of equations:
*>          = 'N':  A * X = B  (No transpose)
*>          = 'T':  A**T* X = B  (Transpose)
*>          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
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
*>          of the matrix B.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          The factors L and U from the factorization A = P*L*U
*>          as computed by SGETC2.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
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
*> \param[in,out] B
*> \verbatim
*>          B is REAL array, dimension (LDB,NRHS)
*>          On entry, the right hand side matrix B.
*>          On exit, the solution matrix X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
*> \ingroup getrs
*
*  =====================================================================

      SUBROUTINE sgetrs_gecp( TRANS, N, NRHS, A, LDA, IPIV, JPIV,
     $                        B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), JPIV( * )
      REAL               A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      parameter( one = 1.0e+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           slaswp, strsm, xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max
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
      ELSE IF( ldb.LT.max( 1, n ) ) THEN
         info = -8
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'SGETRS', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 .OR. nrhs.EQ.0 )
     $   RETURN
*
      IF( notran ) THEN
*
*        Solve A * X = B.
*
*        Apply row interchanges to the right hand sides.
*
         CALL slaswp( nrhs, b, ldb, 1, n, ipiv, 1 )
*
*        Solve L*X = B, overwriting B with X.
*
         CALL strsm( 'Left', 'Lower', 'No transpose', 'Unit', n,
     $               nrhs,
     $               one, a, lda, b, ldb )
*
*        Solve U*X = B, overwriting B with X.
*
         CALL strsm( 'Left', 'Upper', 'No transpose', 'Non-unit', n,
     $               nrhs, one, a, lda, b, ldb )
*
*        Apply column interchanges to the right hand sides.
*
         CALL slaswp( nrhs, b, ldb, 1, n, jpiv, -1 )  
      ELSE
*
*        Solve A**T * X = B.
*
*        Apply column interchanges to the solution vectors.
*
         CALL slaswp( nrhs, b, ldb, 1, n, jpiv, 1 )
*
*        Solve U**T *X = B, overwriting B with X.
*
         CALL strsm( 'Left', 'Upper', 'Transpose', 'Non-unit', n,
     $               nrhs,
     $               one, a, lda, b, ldb )
*
*        Solve L**T *X = B, overwriting B with X.
*
         CALL strsm( 'Left', 'Lower', 'Transpose', 'Unit', n, nrhs,
     $               one,
     $               a, lda, b, ldb )
*
*        Apply row interchanges to the solution vectors.
*
         CALL slaswp( nrhs, b, ldb, 1, n, ipiv, -1 )
      END IF
*
      RETURN
*
*     End of SGETRS_GECP
*

      END
