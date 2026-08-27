      subroutine gepp_wrapper(A, b, n, x, rcond, ferr, berr, work1)
c
c     Wrapper for LAPACK SGESVX expert driver to solve Ax = b
c
c     Inputs:
c       A(n,n) - coefficient matrix (float32, Fortran order)
c       b(n)   - right-hand side vector (float32)
c       n      - dimension of the system
c
c     Outputs:
c       x(n)   - solution vector
c       rcond  - reciprocal condition number estimate
c       ferr   - forward error bound
c       berr   - backward error bound
c       work1  - reciprocal pivot growth factor (from WORK(1))
c
      implicit none

c     Input/output parameters
      integer, intent(in) :: n
      real, intent(in) :: A(n,n), b(n)
      real, intent(out) :: x(n), rcond, ferr, berr, work1

c     Local variables for SGESVX
      real :: AF(n,n), Acopy(n,n), bcopy(n)
      real :: R(n), C(n)
      real :: WORK(4*n)
      real :: FERR_ARR(1), BERR_ARR(1)
      integer :: IPIV(n), IWORK(n)
      character :: FACT, TRANS, EQUED
      integer :: INFO, LDA, LDAF, LDB, LDX, NRHS

c     Set SGESVX parameters
      FACT = 'N'      ! No factorization provided, compute fresh
      TRANS = 'N'     ! Solve A*X = B (no transpose)
      EQUED = 'N'     ! No equilibration
      NRHS = 1        ! Single right-hand side
      LDA = n         ! Leading dimension of A
      LDAF = n        ! Leading dimension of AF
      LDB = n         ! Leading dimension of B
      LDX = n         ! Leading dimension of X

c     Copy input arrays (SGESVX may modify A and b)
      Acopy = A
      bcopy = b

c     Call LAPACK SGESVX expert driver
c     SGESVX performs:
c       - LU factorization with partial pivoting
c       - Solve the system
c       - Estimate condition number
c       - Compute error bounds
c       - Return pivot growth factor
      call SGESVX(FACT, TRANS, n, NRHS, Acopy, LDA, AF, LDAF, IPIV,
     &            EQUED, R, C, bcopy, LDB, x, LDX, rcond,
     &            FERR_ARR, BERR_ARR, WORK, IWORK, INFO)

c     Extract scalar outputs from arrays
      ferr = FERR_ARR(1)
      berr = BERR_ARR(1)
      work1 = WORK(1)  ! Reciprocal pivot growth factor

c     Check for errors
c     Note: INFO = N+1 means matrix is singular to working precision
c     but the solution was still computed. This is a warning, not an error.
c     We suppress this message as it's expected for ill-conditioned systems.
      if (INFO .lt. 0) then
         print *, 'SGESVX: illegal value in argument', -INFO
      else if (INFO .gt. 0 .and. INFO .le. n) then
         print *, 'SGESVX: U(',INFO,',',INFO,') is exactly zero'
         print *, 'The factorization has been completed, but the'
         print *, 'factor U is exactly singular, so the solution'
         print *, 'and error bounds could not be computed.'
      endif

      return
      end subroutine gepp_wrapper
