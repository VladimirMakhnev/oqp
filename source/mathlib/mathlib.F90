module mathlib

  use precision, only: dp
  use oqp_linalg

  implicit none

  private
  public orb_to_dens
  public traceprod_sym_packed
  public solve_linear_equations
  public symmetrize_matrix
  public antisymmetrize_matrix
  public triangular_to_full
  public orthogonal_transform_sym
  public orthogonal_transform
  public orthogonal_transform2
  public matrix_invsqrt

  interface pack_matrix
    module procedure :: PACK_F90, PACK_F77
  end interface pack_matrix
  interface unpack_matrix
    module procedure :: UNPACK_F90, UNPACK_F77
  end interface unpack_matrix
  public :: pack_matrix, unpack_matrix
  public :: pack_f90, unpack_f90
  public :: jacobi_rotate_mo

contains

!> @brief Compute density matrix from a set of orbitals and respective occupation numbers
!> @detail Compute the transformation: D = V * X * V^T
!> @param[out]  d   density matrix
!> @param[in]   v   matrix of orbitals
!> @param[in]   x   vector of occupation numbers
!> @param[in]   m   number of columns in `V`
!> @param[in]   n   dimension of `D`
!> @param[in]   ldv leading dimension of `V`
  subroutine orb_to_dens(d, v, x, m, n, ldv)

    use precision, only: dp

    implicit none

    real(kind=dp), intent(out) :: d(*)
    real(kind=dp), intent(in) :: v(ldv,*), x(*)
    integer, intent(in) :: m, n, ldv

    real(kind=dp), allocatable :: d2(:,:), tmp(:,:)
    integer :: n2, i

    allocate(d2(n,n), tmp(n,m))
    do i = 1, m
      tmp(:,i) = v(:,i)*x(m)
    end do

    call dsyr2k('u', 'n', n, m, 0.5_dp, v, ldv, tmp, n, 0.0_dp, d2, n)

    n2 = (n*n+n)/2
    call pack_matrix(d2, d(:n2))

  end subroutine

!> @brief  Compute the trace of the product of two symmetric matrices in packed format
!> @detail The trace is actually an inner product of matrices, assuming they are vectors
!> @param[in]  a  first matrix
!> @param[in]  b  second matrix
!> @param[in]  n  dimension of matrices `A` and `B`
  function traceprod_sym_packed(a, b, n) result(res)
    use precision, only: dp

    implicit none

    integer :: i, k, n, n2
    real(kind=dp) :: res, a(*), b(*)

    n2 = (n*n+n)/2
    res = 2*dot_product(a(:n2), b(:n2))

    ! subtract the product of the diagonal elements (it is counted twice above)
    k = 0
    do i = 1, n
      k = k+i
      res = res - a(k)*b(k)
    end do

  end function

!> @brief Compute the solution to a real system of linear equations A * X = B
!> @detai Wrapper for DSYSV from Lapack
!> @param[in,out] A      general matrix, destroyed on exit.
!> @param[in,out] B      RHS on entry. The solution on exit.
!> @param[in]     n      size of the problem.
!> @param[in]     nrhs   number of RHS vectors
!> @param[in]     lda    leading dimension of matrix A
!> @param[out]    ierr   Error flag, ierr=0 if no errors. Read DSYEV manual for details
  subroutine solve_linear_equations(a, b, n, nrhs, lda, ierr)

    use precision, only: dp
    use messages, only: show_message
    implicit none

    integer, intent(in) :: lda, n, nrhs
    real(dp) :: a(*)
    real(dp) :: b(*)
    integer, intent(INOUT) :: ierr

    real(dp), allocatable :: work(:)
    integer, allocatable :: ipvt(:)
    integer :: lwork
    real(dp) :: rwork(1)
    integer, external :: ilaenv

    call dsysv('U', n, nrhs, a, lda, ipvt, b, n, rwork, -1, ierr)
    lwork = int(rwork(1))

    allocate(work(lwork), ipvt(n))
    call dsysv('U', n, nrhs, a, lda, ipvt, b, n, work, lwork, ierr)
    deallocate(work, ipvt)

    if (ierr /= 0) then
      call show_message('DSYSV FAILED')
    end if
  end subroutine

!> @brief Compute `A = A + A^T` of a square matrix
!> @param[inout] a  square NxN matrix
!> @param[in]    n  matrix dimension
  subroutine symmetrize_matrix(a,n)
    use precision, only: dp
    real(kind=dp), intent(inout) :: a(n,*)
    integer, intent(in) :: n
    integer :: i
    do i = 1, n
      a(i:n,i) = a(i:n,i) + a(i,i:n)
      a(1:i-1,i) = a(i,1:i-1)
    end do
  end subroutine symmetrize_matrix

!> @brief Compute `A = A - A^T` of a square matrix
!> @param[inout] a  square NxN matrix
!> @param[in]    n  matrix dimension
  subroutine antisymmetrize_matrix(a,n)
    use precision, only: dp
    real(kind=dp), intent(inout) :: a(n,*)
    integer, intent(in) :: n
    integer :: i
    do i = 1, n
      a(i:n,i) = a(i:n,i) - a(i,i:n)
      a(1:i-1,i) = -a(i,1:i-1)
    end do
  end subroutine


!> @brief Fill the upper/lower triangle of the symmetric matrix
!         in triangular form
!> @param[inout] a     square NxN matrix in triangular form
!> @param[in]    n     matrix dimension
!> @param[in]    uplo  U if `A` is upper triangular, L if lower triangular
  subroutine triangular_to_full(a,n,uplo)
    use precision, only: dp
    use messages, only: show_message, with_abort

    implicit none

    real(kind=dp), intent(inout) :: a(n,*)
    integer, intent(in) :: n
    character(len=1), intent(in) :: uplo

    integer :: i

    if (uplo=='u' .or. uplo=='U') then
      do i = 1, n
        a(i+1:n,i) = a(i,i+1:n)
      end do
    else if (uplo=='l' .or. uplo=='L') then
      do i = 1, n
        a(i,i+1:n) = a(i+1:n,i)
      end do
    else
      call show_message('Invalid parameter UPLO='//uplo// &
              ' in `triangular_to_full`. Use either `L` or `U`.', with_abort)
    end if
  end subroutine triangular_to_full

!> @brief Compute orthogonal transformation of a symmetric marix A
!>        in packed format:
!>        B = U^T * A * U
!> @param[in]    a      Matrix to transform
!> @param[in]    u      Orthogonal matrix U(ldu,m)
!> @param[in]    n      dimension of matrix A
!> @param[in]    m      dimension of matrix B
!> @param[in]    ldu    leading dimension of matrix U
!> @param[out]   b      Result
!> @param[inout] wrk    Scratch space
!> @author Vladimir Mironov
  subroutine orthogonal_transform_sym(n, m, a, u, ldu, b)
    use messages, only: show_message, with_abort

    implicit none

    integer, intent(in) :: n, m, ldu
    real(kind=8), intent(in) :: a(*)
    real(kind=8), intent(in) :: u(n,*)
    real(kind=8), intent(out) :: b(*)

    real(kind=8), allocatable :: tmp(:,:), a2(:,:)
    integer :: info

    ! Allocate workspace array
    allocate(a2(n,n), tmp(n,m), stat=info)
    if (info /= 0) then
      call show_message('Cannot allocate memory', WITH_ABORT)
    end if

    ! Unpack the symmetric matrix A
    call dtpttr('u', n, a, a2, n, info)
    if (info /= 0) then
      call show_message('(A,I8)', &
              'Error: DTPTTR returned info =', info, WITH_ABORT)
    end if

    ! Compute A * U
    call dsymm('l', 'u', n, m, 1.0d0, a2, n, u, ldu, 0.0d0, tmp, n)

    ! Compute symmetric matrix B = U^T * (A * U)
    call dgemm('t', 'n', m, m, n, 1.0d0, u, ldu, tmp, n, 0.0d0, a2, n)

    ! Pack the symmetric matrix B
    call dtrttp('u', m, a2, n, b, info)
    if (info /= 0) then
      call show_message('(A,I8)', &
              'Error: DTRTTP returned info =', info, WITH_ABORT)
    end if

  end subroutine

!> @brief Compute orthogonal transformation of a square marix
!> @param[in]    trans  If trans='n' compute B = U^T * A * U
!>                      If trans='t' compute B = U * A * U^T
!> @param[in]    ld     Dimension of matrices
!> @param[in]    u      Square orthogonal matrix
!> @param[inout] a      Matrix to transform, optionally output matrix
!> @param[out]   b      Result, can be absent for in-place transform of matrix A
!> @param[inout] wrk    Scratch space, optional
!> @author Vladimir Mironov
  subroutine orthogonal_transform(trans, ld, u, a, b, wrk)
    use messages, only: show_message, with_abort
    implicit none
    character(len=1), intent(in) :: trans
    integer :: ld
    real(kind=dp), intent(in)    :: u(*)
    real(kind=dp), intent(in)    :: a(*)
    real(kind=dp), optional, intent(out)   :: b(*)
    real(kind=dp), optional, target, intent(inout) :: wrk(*)
    real(kind=dp), pointer :: pwrk(:)
    real(kind=dp), allocatable, target :: wrk_internal(:)
    if (present(wrk)) then
      pwrk(1:ld*ld) => wrk(1:ld*ld)
    else
      allocate(wrk_internal(ld*ld))
      pwrk(1:ld*ld) => wrk_internal(1:ld*ld)
    end if
    select case (trans)
    case ('n', 'N')
        call dgemm('n', 'n', ld, ld, ld, &
                   1.0_dp, a,   ld, &
                           u,   ld, &
                   0.0_dp, pwrk, ld)

        if (present(b)) then
        call dgemm('t', 'n', ld, ld, ld, &
                   1.0_dp, u,   ld, &
                           pwrk, ld, &
                   0.0_dp, b,   ld)
       else
        call dgemm('t', 'n', ld, ld, ld, &
                   1.0_dp, u,   ld, &
                           pwrk, ld, &
                   0.0_dp, a,   ld)
       end if
    case ('t', 'T')
        call dgemm('n', 'n', ld, ld, ld, &
                   1.0_dp, u,   ld, &
                           a,   ld, &
                   0.0_dp, pwrk, ld)

        if (present(b)) then
        call dgemm('n', 't', ld, ld, ld, &
                   1.0_dp, pwrk, ld, &
                           u,   ld, &
                   0.0_dp, b,   ld)
       else
        call dgemm('n', 't', ld, ld, ld, &
                   1.0_dp, pwrk, ld, &
                           u,   ld, &
                   0.0_dp, a,   ld)
       end if
    case default
      call show_message('Invalid parameter TRANS='//trans// &
              ' in `orthogonal_transform`', with_abort)
    end select
  end subroutine
!> @brief Compute orthogonal transformation of a square marix
!> @param[in]    trans  If trans='n' compute B = U^T * A * U
!>                      If trans='t' compute B = U * A * U^T
!> @param[in]    ld     Dimension of matrices
!> @param[in]    u      Square orthogonal matrix
!> @param[in]    a      Matrix to transform
!> @param[out]   b      Result
!> @param[inout] wrk    Scratch space
!> @author Vladimir Mironov
  subroutine orthogonal_transform2(trans, m, n, u, ldu, a, lda, b, ldb, wrk)
    use messages, only: show_message, with_abort
    implicit none
    character(len=1), intent(in) :: trans
    integer :: m, n, ldu, lda, ldb
    real(kind=dp), intent(in)    :: u(*)
    real(kind=dp), intent(in)    :: a(*)
    real(kind=dp), intent(out)   :: b(*)
    real(kind=dp), intent(inout) :: wrk(*)
    select case (trans)
    case ('n', 'N')
        call dgemm('n', 'n', m, n, m, &
                   1.0_dp, a,   lda, &
                           u,   ldu, &
                   0.0_dp, wrk, n)

        call dgemm('t', 'n', n, n, m, &
                   1.0_dp, u,   ldu, &
                           wrk, n, &
                   0.0_dp, b,   ldb)
    case ('t', 'T')
        call dgemm('n', 'n', m, n, n, &
                   1.0_dp, u,   ldu, &
                           a,   lda, &
                   0.0_dp, wrk, n)

        call dgemm('n', 't', m, m, n, &
                   1.0_dp, wrk, n, &
                           u,   ldu, &
                   0.0_dp, b,   ldb)
    case default
      call show_message('Invalid parameter TRANS='//trans// &
              ' in `orthogonal_transform`', with_abort)
    end select
  end subroutine

!> @brief Compute matrix inverse square root using SVD and removing
!>  linear dependency
!
!> @detail This subroutine is used to obtain set of `canonical orbitals`
!>   by diagonalization of the basis set overlap matrix
!>   Q = S^{-1/2}, Q^T * S * Q = I
!
!> @param[in]   s    Overlap matrix, symmetric packed format
!> @param[out]  q    Matrix inverse square root, square matrix
!> @param[in]   nbf  Dimeension of matrices S and Q, basis set size
!> @param[out]  qrnk Rank of matrix Q
!> @param[in]   tol  optional, tolerance to remove linear dependency,
!>                   default = 1.0e-8
  subroutine matrix_invsqrt(s, q, nbf, qrnk, tol)
    use messages,  only: show_message, with_abort
    use eigen,     only: diag_symm_packed
    implicit none

    real(kind=dp), intent(in) :: s(*)
    real(kind=dp), intent(out) :: q(nbf,*)
    integer, intent(in) :: nbf
    real(kind=dp), optional, intent(in) :: tol
    integer, optional, intent(out) :: qrnk

    real(kind=dp), parameter :: deftol = 1.0d-08

    real(kind=dp), allocatable :: tmp(:), eig(:)
    real(kind=dp) :: rtol
    integer :: nbf2, ok, i, j

    rtol = deftol
    if (present(tol)) rtol = tol

    nbf2 = nbf*(nbf+1)/2

    allocate(tmp(nbf2), &
             eig(nbf), &
             stat=ok)
    if (ok/=0) call show_message('Cannot allocate memory', WITH_ABORT)

    tmp(:) = s(1:nbf2)

!   Compute SVD
    call diag_symm_packed(1, nbf, nbf, nbf, tmp, eig, q, ok)

!   Compute Q = S^{-1/2}, eliminating eigenvectors corresponding
!   to small eigenvalues
    j  = 0
    do i = 1, nbf
      if (eig(i) >= rtol) then
        j = j+1
        q(:,j) = q(:,i) / sqrt(eig(i))
      end if
    end do

    q(:,j+1:nbf) = 0

    if (present(qrnk)) qrnk = j

  end subroutine matrix_invsqrt

  !> @brief Fortran-90 routine for packing symmetric matrix to 1D array
  !>
  !> @date      6 October 2021   - Initial release -
  !> @author    Igor S. Gerasimov
  !>
  !> @param[in]          A    - matrix for packing (N x N)
  !> @param[out]         AP   - packed matrix ( N*(N+1)/2 )
  !> @param[in,optional] UPLO - format of packed matrix, `U` for upper and `L` for lower
  subroutine PACK_F90(A, AP, UPLO)
    real(dp), intent(IN) :: A(:, :)
    real(dp), intent(OUT) :: AP(:)
    character(len=1), intent(in), optional :: UPLO
    character(len=1) :: DUPLO
    if (.not. present(UPLO)) then
      DUPLO = 'U'
    else
      DUPLO = UPLO
    end if
    call PACK_F77(A, size(A, 1), AP, DUPLO)
  end subroutine PACK_F90
  !> @brief Fortran-77 routine for packing symmetric matrix to 1D array
  !>
  !> @date      6 October 2021   - Initial release -
  !> @author    Igor S. Gerasimov
  !>
  !> @param[in]          A    - matrix for packing (N x N)
  !> @param[in]          N    - shape of matrix A
  !> @param[out]         AP   - packed matrix ( N*(N+1)/2 )
  !> @param[in]          UPLO - format of packed matrix, `U` for upper and `L` for lower
  subroutine PACK_F77(A, N, AP, UPLO) bind(C, name="MTX_PACK")
    use messages, only: show_message, WITH_ABORT
    real(dp), intent(in) :: A(N, *)
    integer, intent(in) :: N
    real(dp), intent(out) :: AP(*)
    character(len=1), intent(in) :: UPLO
    integer :: INFO
    call dtrttp(uplo, n, a, n, ap, info)
    if (info /= 0) then
      call show_message("error in pack procedure. please, check arguments", with_abort)
    end if
  end subroutine PACK_F77
  !> @brief Fortran-90 routine for unpacking 1D array to symmetric matrix
  !>
  !> @details   LAPACK returns only upper or lower filling of matrix
  !<              so then matrix is symmetrised by couple of cycles
  !>
  !> @date      6 October 2021   - Initial release -
  !> @author    Igor S. Gerasimov
  !>
  !> @param[in]          AP   - packed matrix (N x N)
  !> @param[out]         A    - matrix for unpacking ( N*(N+1)/2 )
  !> @param[in,optional] UPLO - format of packed matrix, `U` for upper and `L` for lower
  subroutine UNPACK_F90(AP, A, UPLO)
    real(dp), intent(IN) :: AP(*)
    real(dp), intent(OUT) :: A(:, :)
    character(len=1), intent(in), optional :: UPLO
    character(len=1) :: DUPLO
    if (.not. present(UPLO)) then
      DUPLO = 'U'
    else
      DUPLO = UPLO
    end if
    call UNPACK_F77(AP, A, size(A, 1), DUPLO)
  end subroutine UNPACK_F90
  !> @brief Fortran-77 routine for unpacking 1D array to symmetric matrix
  !>
  !> @details   LAPACK returns only upper or lower filling of matrix
  !<              so then matrix is symmetrised by couple of cycles
  !>
  !> @date      6 October 2021   - Initial release -
  !> @author    Igor S. Gerasimov
  !>
  !> @param[in]          AP   - packed matrix
  !> @param[out]         A    - matrix for unpacking (N x N)
  !> @param[in]          N    - shape of matrix A
  !> @param[in]          UPLO - format of packed matrix, `U` for upper and `L` for lower
  subroutine UNPACK_F77(AP, A, N, UPLO) bind(C, name="MTX_UNPACK")
    use messages, only: show_message, WITH_ABORT
    real(dp), intent(in) :: AP(*)
    real(dp), intent(out) :: A(N, *)
    integer, intent(in) :: N
    character(len=1), intent(in) :: UPLO
    integer :: info, i
    call dtpttr(uplo, n, ap, a, n, info)
    if (INFO /= 0) then
      call show_message("Error in PACK procedure. Please, check arguments", WITH_ABORT)
    end if
    if (UPLO == 'u' .or. UPLO == 'U') then
      do i = 1, N
        A(i+1:N, i) = A(i, i+1:N)
      end do
    else if (UPLO == 'l' .or. UPLO == 'L') then
      do i = 1, N
        A(i, i+1:N) = A(i+1:N, i)
      end do
    else
      call show_message("UNPACK_F77: UPLO can have only `l`, `L`, `u` or `U` value", WITH_ABORT)
    end if
  end subroutine UNPACK_F77

!> @brief Jacobi pair-rotations of MO based on off-diagonal elements of S_mo (overlap matrix)
!> @details Performs Jacobi rotations to align alpha and beta MO orbitals for UHF reference.
!>          Called twice with isegm=0 (occupied-virtual alignment) and isegm=1 (virtual orbitals).
!> @param[inout] mo_a       Alpha MO coefficients
!> @param[inout] mo_b       Beta MO coefficients
!> @param[in]    smat_full  AO overlap matrix (full, not packed)
!> @param[in]    nocca      Number of alpha occupied orbitals
!> @param[in]    nbf        Number of basis functions
!> @param[inout] work1      Work array (nbf x nbf)
!> @param[inout] work2      Work array (nbf x nbf)
!> @param[in]    debug_mode Print debug information
!> @author Vladimir Yu. Makhnev
  subroutine jacobi_rotate_mo(mo_a, mo_b, smat_full, nocca, nbf, work1, work2, debug_mode)
    use precision, only: dp
    use io_constants, only: iw
    use iso_c_binding, only: c_bool

    implicit none

    real(kind=dp), intent(inout), dimension(:,:) :: mo_a, mo_b
    real(kind=dp), intent(in), dimension(:,:) :: smat_full
    integer, intent(in) :: nocca, nbf
    real(kind=dp), intent(inout), dimension(:,:) :: work1, work2
    logical(c_bool), intent(in) :: debug_mode

    real(kind=dp), allocatable :: s_mo(:,:)
    integer :: nmo, ok

    nmo = size(mo_a, 2)

    allocate(s_mo(nbf, nbf), stat=ok)
    if (ok /= 0) return

    ! Perform Jacobi rotations for occupied-virtual alignment (isegm=0)
    call jacobi_rotate_segment(mo_a, mo_b, smat_full, s_mo, nocca, nbf, nmo, &
                               work1, 0, debug_mode)

    ! Perform Jacobi rotations for virtual orbitals (isegm=1)
    call jacobi_rotate_segment(mo_a, mo_b, smat_full, s_mo, nocca, nbf, nmo, &
                               work1, 1, debug_mode)

    deallocate(s_mo)

  end subroutine jacobi_rotate_mo

!> @brief Jacobi rotations for a specific MO segment
!> @param[inout] mo_a       Alpha MO coefficients
!> @param[inout] mo_b       Beta MO coefficients
!> @param[in]    smat_full  AO overlap matrix
!> @param[inout] s_mo       MO overlap matrix (work array)
!> @param[in]    nocca      Number of alpha occupied orbitals
!> @param[in]    nbf        Number of basis functions
!> @param[in]    nmo        Number of MOs
!> @param[inout] work       Work array
!> @param[in]    isegm      Segment: 0=occ-vir, 1=virtual
!> @param[in]    dgprint    Debug print flag
  subroutine jacobi_rotate_segment(mo_a, mo_b, smat_full, s_mo, nocca, nbf, nmo, &
                                   work, isegm, dgprint)
    use precision, only: dp
    use io_constants, only: iw
    use iso_c_binding, only: c_bool

    implicit none

    real(kind=dp), intent(inout), dimension(:,:) :: mo_a, mo_b
    real(kind=dp), intent(in), dimension(:,:) :: smat_full
    real(kind=dp), intent(inout), dimension(:,:) :: s_mo
    integer, intent(in) :: nocca, nbf, nmo, isegm
    real(kind=dp), intent(inout), dimension(:,:) :: work
    logical(c_bool), intent(in) :: dgprint

    integer :: i, p, q, iterj
    integer :: p_start, p_end, q_start
    integer :: i_max, j_max, max_iter
    real(kind=dp) :: max_off, thresh
    logical :: if_conv
    real(kind=dp), parameter :: go2ev = 27.211386245988d+00

    if_conv = .false.
    thresh = 1d-3

    write(iw,'(A)') '                    ++++++++++++++++++++++++++++++++++++++++'
    write(iw,'(A)') '                       MODULE: HF_DFT_Energy'
    write(iw,'(A)') '                       Rotation MO orbitals (Jacobi)'
    write(iw,'(A)') '                    ++++++++++++++++++++++++++++++++++++++++'
    write(iw,'(A)') ''

    ! Calculate overlap between alpha and beta MOs: s_mo = mo_a^T * smat_full * mo_b
    call dgemm('t', 'n', nbf, nbf, nbf, 1.0_dp, mo_a, nbf, smat_full, nbf, 0.0_dp, work, nbf)
    call dgemm('n', 'n', nbf, nbf, nbf, 1.0_dp, work, nbf, mo_b, nbf, 0.0_dp, s_mo, nbf)

    ! Normalize columns
    do i = 1, nbf
      s_mo(:,i) = s_mo(:,i) / max(norm2(s_mo(:,i)), 1.0e-10_dp)
    end do

    if (dgprint) then
      write(iw,'(A)') '-----------------------------------------'
      write(iw,'(A)') 'Diagonal elements of overlap matrix (BEFORE ROTATIONS)'
      write(iw,'(A)') '-----------------------------------------'
      write(iw,'(A)') '# orb.      Overlap'
      write(iw,'(A)') '-----------------------------------------'
      do i = 1, nmo
         write(iw,'(I5,F12.6)') i, s_mo(i,i)
      enddo
      write(iw,'(A)') '-----------------------------------------'
      write(iw,'(A)') ''
    endif

    ! Set iteration parameters based on segment
    if (isegm == 0) then
        p_start = nocca-1
        p_end   = 2
        q_start = 1
        max_iter = 10000
    else if (isegm == 1) then
        p_start = nmo
        p_end   = nocca+1
        q_start = nocca
        max_iter = 10000
    else
        write(iw, *) "WRONG ISEGM"
        return
    endif

    ! Main iteration loop
    do iterj = 1, max_iter

        max_off = 0.0D0
        i_max = -1
        j_max = -1

        ! Find maximum off-diagonal element
        do p = p_start, p_end, -1
            do q = q_start, p-1
                if (dabs(s_mo(p,q)) > max_off) then
                    max_off = dabs(s_mo(p,q))
                    i_max = p
                    j_max = q
                endif
                if (dabs(s_mo(q,p)) > max_off) then
                    max_off = dabs(s_mo(q,p))
                    i_max = q
                    j_max = p
                endif
            end do
        end do

        if (dgprint) write(iw, *) "max_off", max_off, i_max, j_max

        if (max_off <= thresh) then
           write(iw,'("segment ",I1," converged at iter ",I0)') isegm, iterj
           call flush(iw)
           exit
        else if (if_conv) then
            write(iw,'("segment ",I1," reached the min theta at iter ",I0)') isegm, iterj
            call flush(iw)
            exit
        end if

        call jacobi_rotate_pair(mo_a, mo_b, s_mo, nmo, nbf, isegm, i_max, j_max, if_conv)
    enddo

    if (dgprint) then
      write(iw,'(A)') '-----------------------------------------'
      write(iw,'(A)') 'Diagonal elements of overlap matrix (AFTER ROTATIONS)'
      write(iw,'(A)') '-----------------------------------------'
      write(iw,'(A)') '# orb.      Overlap'
      write(iw,'(A)') '-----------------------------------------'
      do i = 1, nmo
         write(iw,'(I5,F12.6)') i, s_mo(i,i)
      enddo
      write(iw,'(A)') '-----------------------------------------'
    endif

    ! Fix signs of orbitals
    call jacobi_check_sign(mo_a, mo_b, smat_full, s_mo, nmo, nbf)

    if (dgprint) then
     write(iw,'(A)') '-----------------------------------------'
     write(iw,'(A)') 'Diagonal elements of overlap matrix (FINAL/SIGN FIXED)'
     write(iw,'(A)') '-----------------------------------------'
     write(iw,'(A)') '# orb.      Overlap'
     write(iw,'(A)') '-----------------------------------------'
     do i = 1, nmo
        write(iw,'(I5,F12.6)') i, s_mo(i,i)
     enddo
     write(iw,'(A)') '-----------------------------------------'
    endif

  end subroutine jacobi_rotate_segment

!> @brief Rotate a pair of MOs using Jacobi rotation
  subroutine jacobi_rotate_pair(mo_a, mo_b, s_mo, nbf, norb, isegm, i_idx, j_idx, if_conv)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(inout), dimension(:,:) :: mo_a, mo_b
    real(kind=dp), intent(inout), dimension(:,:) :: s_mo
    integer, intent(in) :: nbf, norb, isegm, i_idx, j_idx
    logical, intent(inout) :: if_conv

    real(kind=dp) :: tht, aa, bb, cc, dd, att, btt, cth, sth, tmp
    integer :: k

    aa = s_mo(i_idx, i_idx)
    bb = s_mo(j_idx, j_idx)
    cc = s_mo(i_idx, j_idx)
    dd = s_mo(j_idx, i_idx)

    att = 0.5d0 * (aa*aa + bb*bb - cc*cc - dd*dd)

    if (isegm == 0) then
        btt = aa*dd - bb*cc
    else if (isegm == 1) then
        btt = aa*cc - bb*dd
    end if

    tht = 0.5d0 * datan2(btt, att)

    if (abs(tht) < 1.0d-4) then
        if_conv = .true.
        return
    end if

    cth = dcos(tht)
    sth = dsin(tht)

    if (isegm == 0) then
        ! Rotate columns of mo_a: Givens rotation inline
        do k = 1, norb
            tmp = cth * mo_a(k, i_idx) + sth * mo_a(k, j_idx)
            mo_a(k, j_idx) = cth * mo_a(k, j_idx) - sth * mo_a(k, i_idx)
            mo_a(k, i_idx) = tmp
        end do
        ! Rotate rows of s_mo
        do k = 1, nbf
            tmp = cth * s_mo(i_idx, k) + sth * s_mo(j_idx, k)
            s_mo(j_idx, k) = cth * s_mo(j_idx, k) - sth * s_mo(i_idx, k)
            s_mo(i_idx, k) = tmp
        end do
    else
        ! Rotate columns of mo_b: Givens rotation inline
        do k = 1, norb
            tmp = cth * mo_b(k, i_idx) + sth * mo_b(k, j_idx)
            mo_b(k, j_idx) = cth * mo_b(k, j_idx) - sth * mo_b(k, i_idx)
            mo_b(k, i_idx) = tmp
        end do
        ! Rotate columns of s_mo
        do k = 1, nbf
            tmp = cth * s_mo(k, i_idx) + sth * s_mo(k, j_idx)
            s_mo(k, j_idx) = cth * s_mo(k, j_idx) - sth * s_mo(k, i_idx)
            s_mo(k, i_idx) = tmp
        end do
    end if

  end subroutine jacobi_rotate_pair

!> @brief Swap sign of MO orbital
  subroutine jacobi_swap_sign(mo_a, mo_b, nbf, norb, swa, isegm)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(inout), dimension(:,:) :: mo_a, mo_b
    integer, intent(in) :: nbf, norb, swa, isegm

    integer :: i

    if (isegm == 0) then
        do i = 1, nbf
            mo_a(i, swa) = -mo_a(i, swa)
        end do
    else
        do i = 1, nbf
            mo_b(i, swa) = -mo_b(i, swa)
        end do
    end if

  end subroutine jacobi_swap_sign

!> @brief Check and fix signs of MO orbitals based on overlap
  subroutine jacobi_check_sign(mo_a, mo_b, smat, s_mo, nbf, norb)
    use precision, only: dp

    implicit none

    real(kind=dp), intent(inout), dimension(:,:) :: mo_a, mo_b
    real(kind=dp), intent(in), dimension(:,:) :: smat
    real(kind=dp), intent(inout), dimension(:,:) :: s_mo
    integer, intent(in) :: nbf, norb

    real(kind=dp), allocatable :: sq(:,:)
    integer :: i

    ! Fix signs where diagonal overlap is negative
    do i = 1, norb
        if (s_mo(i,i) < 0.0d0) then
            call jacobi_swap_sign(mo_a, mo_b, nbf, norb, i, 1)
        end if
    end do

    ! Recompute overlap matrix
    allocate(sq(nbf,nbf))
    call dgemm('t', 'n', nbf, nbf, nbf, 1.0d0, mo_a, nbf, smat, nbf, 0.0d0, sq, nbf)
    call dgemm('n', 'n', nbf, nbf, nbf, 1.0d0, sq, nbf, mo_b, nbf, 0.0d0, s_mo, nbf)
    deallocate(sq)

  end subroutine jacobi_check_sign

end module
