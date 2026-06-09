module tdhf_mrsf_z_vector_mod

  implicit none

  character(len=*), parameter :: module_name = "tdhf_mrsf_z_vector_mod"

  ! Module-level work arrays for GMRES to avoid repeated allocation
  ! UMRSF beta orbital energies for the gmres apply-operator wrapper (gmres_solve's
  ! apply_operator interface carries only mo_energy_a; umrsf_sfrolhs needs e_a AND e_b).
  real(kind=8), pointer :: umrsf_mo_energy_b(:) => null()
  real(kind=8), allocatable :: gmres_wrk1(:,:), gmres_wrk2(:,:), gmres_wrk3(:,:)
  real(kind=8), allocatable, target :: gmres_pa(:,:,:)
  real(kind=8), allocatable :: gmres_ab1_mo_a(:,:), gmres_ab1_mo_b(:,:)
  logical :: gmres_work_allocated = .false.
  integer :: gmres_nbf = 0
  integer :: gmres_nocca = 0
  integer :: gmres_noccb = 0

contains

  ! Initialize GMRES work arrays
  subroutine init_gmres_work(nbf, nocca, noccb)
    use messages, only: show_message, with_abort
    implicit none
    integer, intent(in) :: nbf, nocca, noccb
    integer :: nvira, nvirb, ok
    
    if (gmres_work_allocated) then
      ! Check if dimensions match
      if (gmres_nbf == nbf .and. gmres_nocca == nocca .and. gmres_noccb == noccb) then
        return  ! Arrays already allocated with correct size
      else
        ! Deallocate old arrays before reallocating
        call cleanup_gmres_work()
      end if
    end if
    
    nvira = nbf - nocca
    nvirb = nbf - noccb
    
    allocate(gmres_wrk1(nbf,nbf), &
             gmres_wrk2(nbf,nbf), &
             gmres_wrk3(nbf,nbf), &
             gmres_pa(nbf,nbf,2), &
             gmres_ab1_mo_a(nocca,nvira), &
             gmres_ab1_mo_b(noccb,nvirb), &
             stat=ok)
    
    if (ok /= 0) then
      call show_message('Cannot allocate GMRES work arrays', with_abort)
    end if
    
    gmres_work_allocated = .true.
    gmres_nbf = nbf
    gmres_nocca = nocca
    gmres_noccb = noccb
    
  end subroutine init_gmres_work
  
  ! Cleanup GMRES work arrays
  subroutine cleanup_gmres_work()
    implicit none
    
    if (allocated(gmres_wrk1)) deallocate(gmres_wrk1)
    if (allocated(gmres_wrk2)) deallocate(gmres_wrk2)
    if (allocated(gmres_wrk3)) deallocate(gmres_wrk3)
    if (allocated(gmres_pa)) deallocate(gmres_pa)
    if (allocated(gmres_ab1_mo_a)) deallocate(gmres_ab1_mo_a)
    if (allocated(gmres_ab1_mo_b)) deallocate(gmres_ab1_mo_b)
    
    gmres_work_allocated = .false.
    gmres_nbf = 0
    gmres_nocca = 0
    gmres_noccb = 0
    
  end subroutine cleanup_gmres_work

  ! GMRES solver for the z-vector equation
  subroutine gmres_solve(apply_operator, apply_precond, b, x, n, restart, max_iter, tol, &
                         infos, basis, molGrid, int2_driver, &
                         nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                         fa, fb, scale_exch, dft, error_out, iter_out, iw)
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use int2_compute, only: int2_compute_t
    use mod_dft_molgrid, only: dft_grid_t
    
    implicit none
    
    interface
      subroutine apply_operator(x_in, x_out, infos, basis, molGrid, int2_driver, &
                               nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                               fa, fb, scale_exch, dft)
        use precision, only: dp
        use types, only: information
        use basis_tools, only: basis_set
        use int2_compute, only: int2_compute_t
        use mod_dft_molgrid, only: dft_grid_t
        real(kind=dp), intent(in) :: x_in(:)
        real(kind=dp), intent(out) :: x_out(:)
        type(information), intent(inout) :: infos
        type(basis_set), pointer :: basis
        type(dft_grid_t), intent(inout) :: molGrid
        type(int2_compute_t), intent(inout) :: int2_driver
        integer, intent(in) :: nocca, noccb, nbf
        real(kind=dp), intent(in) :: mo_a(:,:), mo_b(:,:), mo_energy_a(:)
        real(kind=dp), intent(in) :: fa(:,:), fb(:,:), scale_exch
        logical, intent(in) :: dft
      end subroutine
      
      subroutine apply_precond(x_in, x_out)
        use precision, only: dp
        real(kind=dp), intent(in) :: x_in(:)
        real(kind=dp), intent(out) :: x_out(:)
      end subroutine
    end interface
    
    real(kind=dp), intent(in) :: b(:)
    real(kind=dp), intent(inout) :: x(:)
    integer, intent(in) :: n, restart, max_iter, iw
    real(kind=dp), intent(in) :: tol
    type(information), intent(inout) :: infos
    type(basis_set), pointer :: basis
    type(dft_grid_t), intent(inout) :: molGrid
    type(int2_compute_t), intent(inout) :: int2_driver
    integer, intent(in) :: nocca, noccb, nbf
    real(kind=dp), intent(in) :: mo_a(:,:), mo_b(:,:), mo_energy_a(:)
    real(kind=dp), intent(in) :: fa(:,:), fb(:,:), scale_exch
    logical, intent(in) :: dft
    real(kind=dp), intent(out) :: error_out
    integer, intent(out) :: iter_out
    
    ! Local variables
    real(kind=dp), allocatable :: V(:,:)     ! Krylov basis
    real(kind=dp), allocatable :: H(:,:)     ! Hessenberg matrix
    real(kind=dp), allocatable :: c(:), s(:) ! Givens rotation coefficients
    real(kind=dp), allocatable :: g(:)       ! RHS for least squares
    real(kind=dp), allocatable :: y(:)       ! Solution of least squares
    real(kind=dp), allocatable :: r(:)       ! Residual
    real(kind=dp), allocatable :: w(:)       ! Work vector
    real(kind=dp), allocatable :: Ax(:)      ! A*x
    
    real(kind=dp) :: beta, h_ij, temp, error, error_initial
    integer :: i, j, k, iter, m, restart_count, inner_iter
    logical :: converged
    
    ! Initialize GMRES work arrays ONCE at the beginning
    call init_gmres_work(nbf, nocca, noccb)
    
    ! Allocate workspace
    m = min(restart, n)
    allocate(V(n, m+1))
    allocate(H(m+1, m))
    allocate(c(m))
    allocate(s(m))
    allocate(g(m+1))
    allocate(y(m))
    allocate(r(n))
    allocate(w(n))
    allocate(Ax(n))
    
    iter_out = 0
    restart_count = 0
    converged = .false.
    
    write(iw,'(/," GMRES Solver Parameters:")')
    write(iw,'("   Problem size        : ", I8)') n
    write(iw,'("   Restart dimension   : ", I4)') m
    write(iw,'("   Max iterations      : ", I4)') max_iter
    write(iw,'("   Convergence tol     : ", 1p,e10.3)') tol
    write(iw,'(/," Iteration   Inner  Residual Norm   Reduction")')
    write(iw,'(" ---------   -----  -------------   ---------")')
    call flush(iw)
    
    ! Compute initial residual
    call apply_operator(x, Ax, infos, basis, molGrid, int2_driver, &
                       nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                       fa, fb, scale_exch, dft)
    r = b - Ax
    error_initial = sqrt(dot_product(r, r))
    
    ! Outer iteration loop (restarts)
    do iter = 1, max_iter
      
      restart_count = restart_count + 1
      
      ! Compute initial residual r = b - A*x
      call apply_operator(x, Ax, infos, basis, molGrid, int2_driver, &
                         nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                         fa, fb, scale_exch, dft)
      r = b - Ax
      
      ! Apply preconditioner to residual
      call apply_precond(r, V(:,1))
      
      beta = sqrt(dot_product(V(:,1), V(:,1)))
      
      ! Check for convergence
      error = beta
      if (iter == 1) then
        write(iw,'(I6,8x,"  0",2x,1p,F13.8,1x,F13.8)') &
              restart_count, error, error/error_initial
      end if
      
      if (error < tol) then
        error_out = error
        iter_out = iter
        converged = .true.
        exit
      end if
      
      V(:,1) = V(:,1) / beta
      g = 0.0_dp
      g(1) = beta
      
      ! Reset H matrix for this restart
      H = 0.0_dp
      
      ! Arnoldi process
      inner_iter = 0
      do j = 1, m
        inner_iter = j
        
        ! Apply operator to V_j
        call apply_operator(V(:,j), w, infos, basis, molGrid, int2_driver, &
                           nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                           fa, fb, scale_exch, dft)
        
        ! Apply preconditioner
        call apply_precond(w, V(:,j+1))
        
        ! Modified Gram-Schmidt orthogonalization
        do i = 1, j
          H(i,j) = dot_product(V(:,j+1), V(:,i))
          V(:,j+1) = V(:,j+1) - H(i,j) * V(:,i)
        end do
        
        ! Reorthogonalization for numerical stability
        do i = 1, j
          temp = dot_product(V(:,j+1), V(:,i))
          H(i,j) = H(i,j) + temp
          V(:,j+1) = V(:,j+1) - temp * V(:,i)
        end do
        
        H(j+1,j) = sqrt(dot_product(V(:,j+1), V(:,j+1)))
        
        ! Check for breakdown
        if (abs(H(j+1,j)) < 1.0d-12) then
          write(iw,'(" GMRES: Lucky breakdown at iteration ", I3)') j
          inner_iter = j
          exit
        end if
        
        V(:,j+1) = V(:,j+1) / H(j+1,j)
        
        ! Apply previous Givens rotations
        do i = 1, j-1
          temp = c(i) * H(i,j) + s(i) * H(i+1,j)
          H(i+1,j) = -s(i) * H(i,j) + c(i) * H(i+1,j)
          H(i,j) = temp
        end do
        
        ! Compute new Givens rotation
        call givens_rotation(H(j,j), H(j+1,j), c(j), s(j))
        
        ! Apply new Givens rotation
        H(j,j) = c(j) * H(j,j) + s(j) * H(j+1,j)
        H(j+1,j) = 0.0_dp
        
        temp = c(j) * g(j) + s(j) * g(j+1)
        g(j+1) = -s(j) * g(j) + c(j) * g(j+1)
        g(j) = temp
        
        ! Check convergence
        error = abs(g(j+1)) 
        iter_out = iter_out + 1
        
        ! Print progress every 5 inner iterations or at convergence
        if (mod(j, 5) == 0 .or. error < tol .or. j == m) then
          write(iw,'(I6,8x,I3,2x,1p,F13.8,1x,F13.8)') &
                restart_count, j, error, error/error_initial
          call flush(iw)
        end if
        
        if (error < tol) then
          converged = .true.
          inner_iter = j
          exit
        end if
        
        if (iter_out >= max_iter) then
          inner_iter = j
          exit
        end if
      end do
      
      ! Solve upper triangular system for y
      call back_substitution(H(1:inner_iter,1:inner_iter), g(1:inner_iter), &
                             y(1:inner_iter), inner_iter)
      
      ! Update solution: x = x + V*y
      do i = 1, inner_iter
        x = x + y(i) * V(:,i)
      end do
      
      error_out = error
      
      if (converged .or. iter_out >= max_iter) exit
      
      ! Print restart information
      if (.not. converged .and. inner_iter == m) then
        write(iw,'(" GMRES: Restarting (restart #", I3, ")")') restart_count
        call flush(iw)
      end if
      
    end do
    
    ! Final status
    write(iw,'(" ---------   -----  -------------   ---------")')
    if (converged) then
      write(iw,'(" GMRES converged in ", I4, " iterations (", I3, " restarts)")') &
            iter_out, restart_count-1
      write(iw,'(" Final residual norm: ", 1p,e13.6)') error_out
    else
      write(iw,'(" GMRES did not converge within ", I4, " iterations")') max_iter
      write(iw,'(" Final residual norm: ", 1p,e13.6)') error_out
    end if
    write(iw,'(" Relative reduction : ", 1p,e13.6)') error_out/error_initial
    call flush(iw)
    
    ! Clean up local arrays
    deallocate(V, H, c, s, g, y, r, w, Ax)
    
    ! NOTE: Do NOT clean up GMRES work arrays here - they will be cleaned in main routine
    
  contains
    
    subroutine givens_rotation(a, b, c, s)
      real(kind=dp), intent(in) :: a, b
      real(kind=dp), intent(out) :: c, s
      real(kind=dp) :: r
      
      if (abs(b) < 1.0d-14) then
        c = 1.0_dp
        s = 0.0_dp
      else
        r = sqrt(a*a + b*b)
        c = a / r
        s = b / r
      end if
    end subroutine givens_rotation
    
    subroutine back_substitution(A, b, x, n)
      integer, intent(in) :: n
      real(kind=dp), intent(in) :: A(n,n), b(n)
      real(kind=dp), intent(out) :: x(n)
      integer :: i, j
      
      x(n) = b(n) / A(n,n)
      do i = n-1, 1, -1
        x(i) = b(i)
        do j = i+1, n
          x(i) = x(i) - A(i,j) * x(j)
        end do
        x(i) = x(i) / A(i,i)
      end do
    end subroutine back_substitution
    
  end subroutine gmres_solve

  ! Apply the z-vector operator (A*x) - Modified to use module-level arrays
  subroutine apply_z_operator(x_in, x_out, infos, basis, molGrid, int2_driver, &
                              nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                              fa, fb, scale_exch, dft)
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use int2_compute, only: int2_compute_t
    use tdhf_lib, only: int2_tdgrd_data_t
    use tdhf_sf_lib, only: sfrogen, sfrolhs
    use mod_dft_gridint_fxc, only: utddft_fxc
    use mathlib, only: symmetrize_matrix, orthogonal_transform
    use mod_dft_molgrid, only: dft_grid_t
    use tdhf_lib, only: mntoia
    
    implicit none
    
    real(kind=dp), intent(in) :: x_in(:)
    real(kind=dp), intent(out) :: x_out(:)
    type(information), intent(inout) :: infos
    type(basis_set), pointer :: basis
    type(dft_grid_t), intent(inout) :: molGrid
    type(int2_compute_t), intent(inout) :: int2_driver
    integer, intent(in) :: nocca, noccb, nbf
    real(kind=dp), intent(in) :: mo_a(:,:), mo_b(:,:), mo_energy_a(:)
    real(kind=dp), intent(in) :: fa(:,:), fb(:,:), scale_exch
    logical, intent(in) :: dft
    
    ! Local variables
    real(kind=dp), pointer :: ab1(:,:,:)
    type(int2_tdgrd_data_t), allocatable, target :: int2_data
    integer :: nvira, nvirb
    
    nvira = nbf - nocca
    nvirb = nbf - noccb
    
    ! Ensure work arrays are initialized and have correct dimensions
    if (.not. gmres_work_allocated .or. &
        gmres_nbf /= nbf .or. &
        gmres_nocca /= nocca .or. &
        gmres_noccb /= noccb) then
      call init_gmres_work(nbf, nocca, noccb)
    end if
    
    ! Clear work arrays
    gmres_wrk1 = 0.0_dp
    gmres_wrk2 = 0.0_dp
    gmres_wrk3 = 0.0_dp
    gmres_pa = 0.0_dp
    gmres_ab1_mo_a = 0.0_dp
    gmres_ab1_mo_b = 0.0_dp
    
    ! Generate density matrices from x_in
    call sfrogen(gmres_wrk1, gmres_wrk2, x_in, nocca, noccb)
    
    ! Transform to AO basis
    call orthogonal_transform('t', nbf, mo_a, gmres_wrk1, gmres_pa(:,:,1), gmres_wrk3)
    call orthogonal_transform('t', nbf, mo_b, gmres_wrk2, gmres_pa(:,:,2), gmres_wrk3)
    
    ! Initialize ERI calculation with proper allocation
    allocate(int2_data)
    int2_data = int2_tdgrd_data_t( &
        d2 = gmres_pa, &
        int_apb = .true., &
        int_amb = .false., &
        tamm_dancoff = .false., &
        scale_exchange = scale_exch)
    
    call int2_driver%run(int2_data, &
          cam=dft.and.infos%dft%cam_flag, &
          alpha=infos%dft%cam_alpha, &
          beta=infos%dft%cam_beta,&
          mu=infos%dft%cam_mu)
    ab1 => int2_data%apb(:,:,:,1)
    
    call symmetrize_matrix(gmres_pa(:,:,1), nbf)
    call symmetrize_matrix(gmres_pa(:,:,2), nbf)
    
    if (dft) then
      call utddft_fxc( &
          basis = basis, &
          molGrid = molGrid, &
          isVecs = .true., &
          wfa = mo_a, &
          wfb = mo_b, &
          fxa = ab1(:,:,1:1), &
          fxb = ab1(:,:,2:2), &
          dxa = gmres_pa(:,:,1:1), &
          dxb = gmres_pa(:,:,2:2), &
          nmtx = 1, &
          threshold = 1.0d-15, &
          infos = infos)
    end if
    
    ! Transform to MO basis - Fixed to use correct mo_b for beta
    call mntoia(ab1(:,:,1), gmres_ab1_mo_a, mo_a, mo_a, nocca, nocca)
    call mntoia(ab1(:,:,2), gmres_ab1_mo_b, mo_b, mo_b, noccb, noccb)
    
    ! Apply the operator
    call sfrolhs(x_out, x_in, mo_energy_a, fa, fb, gmres_ab1_mo_a, gmres_ab1_mo_b, &
                 nocca, noccb)
    
    call int2_data%clean()
    deallocate(int2_data)
    
  end subroutine apply_z_operator

  ! UMRSF J operator: x_out = J . x_in, 4-block layout (OV_a, CV_a, CO_b, CV_b).
  ! Mirrors apply_z_operator but uses umrsf_sfrogen/umrsf_sfrolhs. The single shared
  ! int2_tdgrd(A+B) + utddft_fxc run on the independent alpha/beta perturbation
  ! densities supplies the cross-spin Coulomb + f^xc_alphabeta coupling.
  ! Соответствует §6 (eq:Zvec) umrsf_gradient_theory.tex.
  subroutine umrsf_apply_z_operator(x_in, x_out, infos, basis, molGrid, int2_driver, &
                              nocca, noccb, nbf, mo_a, mo_b, e_a, e_b, fa, fb, scale_exch, dft)
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use int2_compute, only: int2_compute_t
    use tdhf_lib, only: int2_tdgrd_data_t, mntoia
    use tdhf_mrsf_lib, only: umrsf_sfrogen, umrsf_sfrolhs
    use mod_dft_gridint_fxc, only: utddft_fxc
    use mathlib, only: symmetrize_matrix, orthogonal_transform
    use mod_dft_molgrid, only: dft_grid_t

    implicit none

    real(kind=dp), intent(in) :: x_in(:)
    real(kind=dp), intent(out) :: x_out(:)
    type(information), intent(inout) :: infos
    type(basis_set), pointer :: basis
    type(dft_grid_t), intent(inout) :: molGrid
    type(int2_compute_t), intent(inout) :: int2_driver
    integer, intent(in) :: nocca, noccb, nbf
    real(kind=dp), intent(in) :: mo_a(:,:), mo_b(:,:)
    real(kind=dp), intent(in) :: e_a(:), e_b(:)
    real(kind=dp), intent(in) :: fa(:,:), fb(:,:), scale_exch
    logical, intent(in) :: dft

    real(kind=dp), pointer :: ab1(:,:,:)
    type(int2_tdgrd_data_t), allocatable, target :: int2_data

    if (.not. gmres_work_allocated .or. gmres_nbf /= nbf .or. &
        gmres_nocca /= nocca .or. gmres_noccb /= noccb) then
      call init_gmres_work(nbf, nocca, noccb)
    end if

    gmres_wrk1 = 0.0_dp; gmres_wrk2 = 0.0_dp; gmres_wrk3 = 0.0_dp
    gmres_pa = 0.0_dp; gmres_ab1_mo_a = 0.0_dp; gmres_ab1_mo_b = 0.0_dp

  ! 4-block Z -> alpha/beta MO densities -> AO
    call umrsf_sfrogen(gmres_wrk1, gmres_wrk2, x_in, nocca, noccb)
    call orthogonal_transform('t', nbf, mo_a, gmres_wrk1, gmres_pa(:,:,1), gmres_wrk3)
    call orthogonal_transform('t', nbf, mo_b, gmres_wrk2, gmres_pa(:,:,2), gmres_wrk3)

  ! (A+B) ERI response (Coulomb couples alpha<->beta; exchange same-spin)
    allocate(int2_data)
    int2_data = int2_tdgrd_data_t(d2=gmres_pa, int_apb=.true., int_amb=.false., &
        tamm_dancoff=.false., scale_exchange=scale_exch)
    call int2_driver%run(int2_data, cam=dft.and.infos%dft%cam_flag, &
          alpha=infos%dft%cam_alpha, beta=infos%dft%cam_beta, mu=infos%dft%cam_mu)
    ab1 => int2_data%apb(:,:,:,1)

    call symmetrize_matrix(gmres_pa(:,:,1), nbf)
    call symmetrize_matrix(gmres_pa(:,:,2), nbf)

  ! cross-spin f^xc_alphabeta (UHF-aware kernel: fxa gets dxa AND dxb)
    if (dft) then
      call utddft_fxc(basis=basis, molGrid=molGrid, isVecs=.true., wfa=mo_a, wfb=mo_b, &
          fxa=ab1(:,:,1:1), fxb=ab1(:,:,2:2), dxa=gmres_pa(:,:,1:1), dxb=gmres_pa(:,:,2:2), &
          nmtx=1, threshold=1.0d-15, infos=infos)
    end if

  ! response AO -> MO (occ x vir, per spin)
    call mntoia(ab1(:,:,1), gmres_ab1_mo_a, mo_a, mo_a, nocca, nocca)
    call mntoia(ab1(:,:,2), gmres_ab1_mo_b, mo_b, mo_b, noccb, noccb)

  ! assemble J . z
    call umrsf_sfrolhs(x_out, x_in, e_a, e_b, fa, fb, gmres_ab1_mo_a, gmres_ab1_mo_b, nocca, noccb)

    call int2_data%clean()
    deallocate(int2_data)

  end subroutine umrsf_apply_z_operator

  ! Wrapper matching gmres_solve's apply_operator interface (which passes mo_energy_a
  ! only). Routes to umrsf_apply_z_operator with e_a=mo_energy_a and e_b from the module
  ! pointer umrsf_mo_energy_b (set by the caller before gmres_solve).
  subroutine umrsf_apply_z_op_gmres(x_in, x_out, infos, basis, molGrid, int2_driver, &
                              nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
                              fa, fb, scale_exch, dft)
    use precision, only: dp
    use types, only: information
    use basis_tools, only: basis_set
    use int2_compute, only: int2_compute_t
    use mod_dft_molgrid, only: dft_grid_t
    implicit none
    real(kind=dp), intent(in) :: x_in(:)
    real(kind=dp), intent(out) :: x_out(:)
    type(information), intent(inout) :: infos
    type(basis_set), pointer :: basis
    type(dft_grid_t), intent(inout) :: molGrid
    type(int2_compute_t), intent(inout) :: int2_driver
    integer, intent(in) :: nocca, noccb, nbf
    real(kind=dp), intent(in) :: mo_a(:,:), mo_b(:,:), mo_energy_a(:)
    real(kind=dp), intent(in) :: fa(:,:), fb(:,:), scale_exch
    logical, intent(in) :: dft
    call umrsf_apply_z_operator(x_in, x_out, infos, basis, molGrid, int2_driver, &
         nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, umrsf_mo_energy_b, &
         fa, fb, scale_exch, dft)
  end subroutine umrsf_apply_z_op_gmres

  ! Apply preconditioner (simple diagonal preconditioner)
  subroutine apply_z_precond(x_in, x_out, xminv)
    use precision, only: dp
    implicit none
    real(kind=dp), intent(in) :: x_in(:), xminv(:)
    real(kind=dp), intent(out) :: x_out(:)
    
    x_out = xminv * x_in
    
  end subroutine apply_z_precond

  subroutine tdhf_mrsf_z_vector_C(c_handle) bind(C, name="tdhf_mrsf_z_vector")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_mrsf_z_vector(inf)
  end subroutine tdhf_mrsf_z_vector_C

  subroutine tdhf_mrsf_z_vector(infos)
    use precision, only: dp
    use io_constants, only: iw
    use oqp_tagarray_driver

    use types, only: information
    use strings, only: Cstring, fstring
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use util, only: measure_time

    use int2_compute, only: int2_compute_t
    use tdhf_lib, only: int2_td_data_t
    use tdhf_lib, only: int2_tdgrd_data_t
    use tdhf_mrsf_lib, only: int2_mrsf_data_t
    use tdhf_lib, only: iatogen, mntoia
    use tdhf_sf_lib, only: sfrorhs, &
      sfromcal, sfrogen, sfrolhs, pcgrbpini, &
      pcgb, sfropcal, sfdmat
    use dft, only: dft_initialize, dftclean
    use mod_dft_gridint_fxc, only: utddft_fxc
    use mathlib, only: symmetrize_matrix, orthogonal_transform, &
            orthogonal_transform_sym
    use mod_dft_molgrid, only: dft_grid_t
    use mathlib, only: pack_matrix, unpack_matrix

    use tdhf_mrsf_lib, only: &
      mrinivec, mrsfcbc, mrsfxvec, mrsfsp, mrsfrowcal, &
      mrsfqrorhs, mrsfqropcal, mrsfqrowcal
    use oqp_linalg
    use printing, only: print_module_info


    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_mrsf_z_vector"

    type(basis_set), pointer :: basis
    type(information), target, intent(inout) :: infos

    integer :: ok

    real(kind=dp), allocatable :: ab1_mo_a(:,:)
    real(kind=dp), allocatable :: ab1_mo_b(:,:)
    real(kind=dp), allocatable :: xm(:)
    real(kind=dp), pointer :: ab2(:,:,:)
    real(kind=dp), pointer :: ab1(:,:,:)
    real(kind=dp), allocatable :: fa(:,:), fb(:,:)
    real(kind=dp), pointer :: bvec(:,:,:)
    real(kind=dp), pointer :: wmo(:,:)
    real(kind=dp), allocatable :: bvec_mo_d(:,:)
    real(kind=dp), allocatable, target :: &
      fmrst1(:,:,:,:)
    real(kind=dp), pointer :: fmrst2(:,:,:,:)

    integer :: nocca, nvira, noccb, nvirb
    integer :: nbf, nbf_tri
    integer :: iter, gmres_iter
    real(kind=dp) :: cnvtol, scale_exch, scale_exch2
    logical :: roref = .false.
    integer :: mrst

    type(int2_compute_t) :: int2_driver
    type(int2_mrsf_data_t), allocatable, target :: int2_data_st
    type(int2_td_data_t), allocatable, target :: int2_data_q
    class(int2_td_data_t), allocatable, target :: int2_data
    type(dft_grid_t) :: molGrid

  ! scr data
    real(kind=dp), allocatable, target :: wrk1(:,:), wrk2(:,:), wrk3(:,:)
    real(kind=dp), pointer :: wrk1t(:)

  ! SF-TD Gradient data
    real(kind=dp), allocatable :: &
      rhs(:), lhs(:), xminv(:), xk(:), pk(:), errv(:), &
      hxa(:,:), hxb(:,:), tij(:,:), ppija(:,:), ppijb(:,:), tab(:,:)
    real(kind=dp), allocatable, target :: pa(:,:,:)
    integer :: nsocc, lzdim, xvec_dim

  ! General data
    real(kind=dp) :: alpha, error
    character(len=10) :: solver_name

    logical :: dft
    integer :: scf_type, mol_mult, target_state

    ! tagarray
    real(kind=dp), contiguous, pointer :: &
      fock_a(:), mo_a(:,:), mo_energy_a(:), &
      fock_b(:), mo_b(:,:), &
      td_p(:,:), td_t(:,:), ta(:), tb(:), td_abxc(:,:), &
      td_mrsf_den(:,:,:), bvec_mo(:,:), wao(:), mrsf_energies(:)
    character(len=*), parameter :: tags_alloc(4) = (/ character(len=80) :: &
      OQP_WAO, OQP_td_mrsf_density, OQP_td_p, OQP_td_abxc /)
    character(len=*), parameter :: tags_required(8) = (/ character(len=80) :: &
      OQP_FOCK_A, OQP_E_MO_A, OQP_VEC_MO_A, OQP_FOCK_B, OQP_VEC_MO_B, OQP_td_bvec_mo, OQP_td_t, &
      OQP_td_energies /)

    mol_mult = infos%mol_prop%mult
    if (mol_mult/=3) call show_message(&
            'MRSF-TDDFT are available for ROHF/UHF ref.&
            &with ONLY triplet multiplicity(mult=3)', with_abort)

    scf_type = infos%control%scftype
    if (scf_type==3) roref = .true.

    dft = infos%control%hamilton == 20

  ! Files open
  ! 3. LOG: Write: Main output file
    open (unit=iw, file=infos%log_filename, position="append")
  !
    call print_module_info('MRSF_TDHF_Z_Vector','Solving Z-Vector for MRSF-TDDFT')

  ! Readings

  ! Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

    nbf = basis%nbf
    nbf_tri = nbf*(nbf+1)/2

    if (dft) call dft_initialize(infos, basis, molGrid)

  ! Parameter it should be inputed later
    mrst = infos%tddft%mult
    cnvtol = infos%tddft%zvconv

    nocca = infos%mol_prop%nelec_A
    nvira = nbf-noccA
    noccb = infos%mol_prop%nelec_B
    nvirb = nbf-noccb
    nsocc = nocca-noccb
    lzdim = noccb*(nsocc+nvira)+nsocc*nvira

    if(mrst==1 .or. mrst==3) then
      xvec_dim = nocca*nvirb
      allocate(&
    ! for Z-vector
        fmrst1(1,7,nbf,nbf), &
        bvec_mo_d(xvec_dim,1), &
        hxa(nbf,nocca), &
        hxb(nbf,nbf), &
    ! for gradient
        tij(nocca,nocca), &
        tab(nvirb,nvirb), &
        stat=ok, &
        source=0.0_dp)
    else if(mrst==5) then
      xvec_dim = noccb*nvira
      allocate(&
    ! for Z-vector
        hxa(nbf,nbf), &
        hxb(nbf,noccb), &
  ! for gradient
        tij(noccb,noccb), &
        tab(nvira,nvira), &
        stat=ok, &
        source=0.0_dp)
    endif
    if( ok/=0 ) call show_message('Cannot allocate memory', with_abort)

    allocate(&
  ! for Z-vector
      xminv(lzdim), &
      rhs(lzdim), &
      lhs(lzdim), &
      xm(lzdim), &
      xk(lzdim), &
      pk(lzdim), &
      errv(lzdim), &
   ! For gradient
      pa(nbf,nbf,2), &
      ppija(nocca,nocca), &
      ppijb(noccb,noccb), &
   ! Allocate TDDFT variables
      fa(nbf,nbf), &           ! Temporary matrix for diagonalization
      fb(nbf,nbf), &           ! Temporary matrix for diagonalization
      ab1_mo_a(nocca,nvira), &
      ab1_mo_b(noccb,nvirb), &
!   For scratch
      wrk1(nbf,nbf), &
      wrk2(nbf,nbf), &
      wrk3(nbf,nbf), &
      stat=ok, &
      source=0.0_dp)

    if( ok/=0 ) call show_message('Cannot allocate memory', with_abort)

    call infos%dat%remove_records(tags_alloc)

    call infos%dat%reserve_data(OQP_WAO, TA_TYPE_REAL64, nbf_tri, comment=OQP_WAO_comment)
    call infos%dat%reserve_data(OQP_td_mrsf_density, TA_TYPE_REAL64, nbf*nbf*7, (/7, nbf, nbf /), comment=OQP_td_mrsf_density)
    call infos%dat%reserve_data(OQP_td_p, TA_TYPE_REAL64, nbf_tri*2, (/ nbf_tri, 2 /), comment=OQP_td_p)
    call infos%dat%reserve_data(OQP_td_abxc, TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), comment=OQP_td_abxc)

    call data_has_tags(infos%dat, tags_alloc, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_WAO, wao)
    call tagarray_get_data(infos%dat, OQP_td_mrsf_density, td_mrsf_den)
    call tagarray_get_data(infos%dat, OQP_td_p, td_p)
    call tagarray_get_data(infos%dat, OQP_td_abxc, td_abxc)

    call data_has_tags(infos%dat, tags_required, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
    call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
    call tagarray_get_data(infos%dat, OQP_td_t, td_t)
    call tagarray_get_data(infos%dat, OQP_td_energies, mrsf_energies)

    ta => td_t(:,1)
    tb => td_t(:,2)

    target_state = min(infos%tddft%target_state, infos%tddft%nstate)
    if (target_state /=infos%tddft%target_state) then
      write(*,'(/1x,66("-")&
               &/1x,"WARNING: Target state has been changed to the max available nstates"/&
               &/1x,66("-")/)')
    end if

    ! Determine solver name for output
    if (infos%tddft%z_solver == 1) then
      solver_name = "GMRES"
    else
      solver_name = "CG"
    end if

    ! Save unrelaxed density matrices and the `b=A*x` vector for target state
    if (mrst==1 .or. mrst==3 ) then
      call mrsfxvec(infos, bvec_mo(:,target_state), bvec_mo_d(:,1))
      call sfdmat(bvec_mo_d(:,1), td_abxc, mo_a, ta, tb, nocca, noccb)
    else if (mrst==5 ) then
      call sfdmat(bvec_mo(:,target_state), td_abxc, mo_a, tb, ta, noccb, nocca)
    end if

    bvec(1:nbf,1:nbf,1:1) => td_abxc

  ! Initialize ERI calculations
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()

    write(*,'(/1x,71("-")&
             &/19x,"MRSF-DFT ENERGY GRADIENT CALCULATION"&
             &/1x,71("-")/)')

    write(iw,fmt='(5x,a/&
                  &5x,16("-")/&
                  &5x,a,x,i0,x,f17.10,x,"Hartree"/&
                  &5x,a,x,i0/&
                  &5x,a,x,e10.4/&
                  &5x,a,x,i0/&
                  &5x,a,x,a)') &
        'Z-vector options' &
      , 'Target state       is', target_state, infos%mol_energy%energy+mrsf_energies(target_state) &
      , 'Multiplicity       is', infos%tddft%mult &
      , 'Convergence        is', infos%tddft%zvconv &
      , 'Maximum iterations is', infos%control%maxit_zv &
      , 'Solver method      is', trim(solver_name)
    call flush(iw)

  ! Prepare for ROHF
    ! Fock matrices A and B
    if( roref )then
        wrk1t(1:nbf*nbf) => wrk1
  !   Alapha
      call orthogonal_transform_sym(nbf, nbf, fock_a, mo_a, nbf, wrk1)
      call unpack_matrix(wrk1t, fa)

  !   Beta
      call orthogonal_transform_sym(nbf, nbf, fock_b, mo_b, nbf, wrk1)
      call unpack_matrix(wrk1t, fb)
    end if

  ! Make density like part
    call unpack_matrix(ta, pa(:,:,1))
    call unpack_matrix(tb, pa(:,:,2))

  ! Initialize ERI calculations
    scale_exch = 1.0_dp
    scale_exch2 = 1.0_dp
    if (dft) then
       scale_exch = infos%dft%HFscale    !> Reference HF exchange
       scale_exch2 = infos%tddft%HFscale !> Response HF exchange
    end if

    if (mrst==1 .or. mrst==3 ) then

      int2_data_st = int2_mrsf_data_t( &
          d3 = fmrst1, &
          tamm_dancoff = .true., &
          scale_exchange = scale_exch2, &
          scale_coulomb = scale_exch2)

    else if( mrst==5  )then

      int2_data_q = int2_td_data_t( &
          d2=bvec, &
          int_apb = .false., &
          int_amb = .false., &
          tamm_dancoff = .true., &
          scale_exchange = scale_exch2)

    end if

    int2_data = int2_tdgrd_data_t( &
        d2 = pa, &
        int_apb = .true., &
        int_amb = .false., &
        tamm_dancoff = .false., &
        scale_exchange = scale_exch)

    call int2_driver%run(int2_data, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%dft%cam_alpha, &
            beta=infos%dft%cam_beta,&
            mu=infos%dft%cam_mu)
    ab1 => int2_data%apb(:,:,:,1)

    pa = pa*2
    call utddft_fxc( &
        basis = basis, &
        molGrid = molGrid, &
        isVecs = .true., &
        wfa = mo_a, &
        wfb = mo_b, &
        fxa = ab1(:,:,1:1), &
        fxb = ab1(:,:,2:2), &
        dxa = pa(:,:,1:1), &
        dxb = pa(:,:,2:2), &
        nmtx = 1, &
        threshold = 1.0d-15, &
        infos = infos)

!   ALPHA: AO(M,N) -> MO(IA+)
    call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)

    call mntoia(ab1(:,:,2), ab1_mo_b, mo_b, mo_b, noccb, noccb)

    if (mrst==1 .or. mrst==3) then

      call iatogen(bvec_mo(:,target_state), wrk1, nocca, noccb)
      call mrsfcbc(infos, mo_a, mo_a, wrk1, fmrst1(1,:,:,:))

      fmrst1(1,7,:,:) = td_abxc

      td_mrsf_den(1:7,:,:) = fmrst1(1,1:7,:,:)

    ! Initialize ERI calculations
      call int2_driver%run(int2_data_st, &
            cam = dft.and.infos%dft%cam_flag, &
            alpha = infos%tddft%cam_alpha, &
            alpha_coulomb = infos%tddft%cam_alpha, &
            beta = infos%tddft%cam_beta,&
            beta_coulomb = infos%tddft%cam_beta, &
            mu = infos%tddft%cam_mu)
      fmrst2 => int2_data_st%f3(:,:,:,:,1)! ado2v, ado1v, adco1, adco2, ao21v, aco12, agdlr

    ! Scaling factor if triplet
      if (mrst==3) fmrst2(:,1:6,:,:) = -1.0_dp*fmrst2(:,1:6,:,:)

      ! Spin pair coupling
      if (infos%tddft%spc_coco /= infos%tddft%hfscale) &
         fmrst2(:,6,:,:) = fmrst2(:,6,:,:) * infos%tddft%spc_coco / infos%tddft%hfscale
      if (infos%tddft%spc_ovov /= infos%tddft%hfscale) &
         fmrst2(:,5,:,:) = fmrst2(:,5,:,:) * infos%tddft%spc_ovov / infos%tddft%hfscale
      if (infos%tddft%spc_coov /= infos%tddft%hfscale) &
         fmrst2(:,1:4,:,:) = fmrst2(:,1:4,:,:) * infos%tddft%spc_coov / infos%tddft%hfscale

      call orthogonal_transform('n', nbf, mo_a, fmrst2(1,7,:,:), wrk2, wrk1)

      call mrsfxvec(infos, bvec_mo(:,target_state), bvec_mo_d(:,1))

      call iatogen(bvec_mo_d(:,1), wrk3, nocca, noccb)

      call dgemm('n', 't', nbf, nocca, nbf, &
                 2.0_dp, wrk2, nbf, &
                         wrk3, nbf, &
                 0.0_dp, hxa, nbf)
      call dgemm('t', 'n', nbf, nbf, nocca, &
                 2.0_dp, wrk2, nbf, &
                         wrk3, nbf, &
                 0.0_dp, hxb, nbf)

   ! spin pair ov-ov, co-co, co-ov coupling
      call mrsfsp(hxa, hxb, mo_a, mo_a, wrk3, fmrst2(1,:,:,:), nocca, noccb)

   !  Unrelaxed difference density matries T_ij and T_ab
   !  Ta(i+,j+):= -X(i+,a-)*X(j+,a-) for singlet and triplet
      call dgemm('n', 't', nocca, nocca, nvirb, &
                -1.0_dp, bvec_mo_d, nocca, &
                         bvec_mo_d, nocca, &
                 0.0_dp, tij, nocca)

   !  Tb(a-,b-):= X(i+,a-)*X(i+,b-) for singlet and triplet
      call dgemm('t', 'n', nvirb, nvirb, nocca, &
                 1.0_dp, bvec_mo_d, nocca, &
                         bvec_mo_d, nocca, &
                 0.0_dp, tab, nvirb)

      call sfrorhs(rhs, hxa, hxb, ab1_mo_a, ab1_mo_b, &
                   Tij, Tab, Fa, Fb, nocca, noccb)

    else if(mrst==5) then

   !  Initialize ERI calculations
      call int2_driver%run(int2_data_q, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%tddft%cam_alpha, &
            beta=infos%tddft%cam_beta,&
            mu=infos%tddft%cam_mu)

      call orthogonal_transform('n', nbf, mo_a, int2_data_q%amb(:,:,1,1), wrk2, wrk1)

      call iatogen(bvec_mo(:,target_state),wrk3,noccb,nocca)

      call dgemm('t', 'n', nbf, nbf, noccb, &
                 2.0_dp, wrk2, nbf, &
                         wrk3, nbf, &
                 0.0_dp, hxa, nbf)
      call dgemm('n', 't', nbf, noccb, nbf, &
                 2.0_dp, wrk2, nbf, &
                         wrk3, nbf, &
                 0.0_dp, hxb, nbf)

   !  Unrelaxed difference density matries T_ij and T_ab
   !  Ta(i+,j+):= -X(i+,a-)*X(j+,a-) for singlet and triplet
      call dgemm('n', 't', noccb, noccb, nvira, &
                -1.0_dp, bvec_mo(:,target_state), noccb, &
                         bvec_mo(:,target_state), noccb, &
                 0.0_dp, tij, noccb)

   !  Tb(a-,b-):= X(i+,a-)*X(i+,b-) for singlet and triplet
      call dgemm('t', 'n', nvira, nvira, noccb, &
                 1.0_dp, bvec_mo(:,target_state), noccb, &
                         bvec_mo(:,target_state), noccb, &
                 0.0_dp, tab, nvira)

      call mrsfqrorhs(rhs, hxa, hxb, ab1_mo_a, ab1_mo_b, &
                      tab, tij, fa, fb, nocca, noccb)
    end if

    write(*,'(/3x,25("-")&
             &/6x,"START Z-VECTOR LOOP (",A,")"&
             &/3x,25("-")/)') trim(solver_name)
    call flush(iw)

    call sfromcal(xm, xminv, mo_energy_a, fa, fb, nocca, noccb)

    ! Choose solver based on input option (0=CG, 1=GMRES)
    if (infos%tddft%z_solver == 1) then
      
      ! ============================================
      ! GMRES SOLVER
      ! ============================================
      
      ! Check preconditioner condition
      if (any(abs(xminv) < 1.0d-12) .or. any(abs(xminv) > 1.0d12)) then
        write(iw,'(" Warning: Preconditioner poorly conditioned, applying regularization")')
        where(abs(xminv) < 1.0d-12) xminv = sign(1.0d-12, xminv)
        where(abs(xminv) > 1.0d12) xminv = sign(1.0d12, xminv)
      end if
      
      ! Initial guess with same strategy as CG
      xk = 0.0_dp

      ! Call GMRES solver
      call gmres_solve( &
          apply_operator = apply_z_operator, &
          apply_precond = lambda_precond, &
          b = rhs, &
          x = xk, &
          n = lzdim, &
          restart = min(infos%tddft%gmres_dim, lzdim), &
          max_iter = infos%control%maxit_zv, &
          tol = cnvtol, &
          infos = infos, &
          basis = basis, &
          molGrid = molGrid, &
          int2_driver = int2_driver, &
          nocca = nocca, &
          noccb = noccb, &
          nbf = nbf, &
          mo_a = mo_a, &
          mo_b = mo_b, &
          mo_energy_a = mo_energy_a, &
          fa = fa, &
          fb = fb, &
          scale_exch = scale_exch, &
          dft = dft, &
          error_out = error, &
          iter_out = gmres_iter, &
          iw = iw)

      write(iw,'(/," Final Summary:")')
      write(iw,'(" GMRES total iterations: ", I4)') gmres_iter
      write(iw,'(" Final error norm      : ", 1p,e13.6)') error
      write(iw,'(" Convergence criterion : ", 1p,e13.6)') cnvtol
      call flush(iw)
      
      ! Clean up GMRES work arrays
      call cleanup_gmres_work()
      
    else
      
      ! ============================================
      ! ORIGINAL CONJUGATE GRADIENT SOLVER
      ! ============================================
      
      ! Initial guess with same strategy as CG
      xk = 0.0_dp

      call sfrogen(wrk1, wrk2, xk, nocca, noccb)
      ! Alpha
      call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)
      ! Beta
      call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)

      !****** INITIAL (A+B)*PK *************************************************
      ! Initialize ERI calculations
      call int2_data%clean()
      deallocate(int2_data)
      int2_data = int2_td_data_t( &
          d2 = pa, &
          int_apb = .true., &
          int_amb = .true., &
          tamm_dancoff = .false., &
          scale_exchange = scale_exch)

      call int2_driver%run(int2_data, &
              cam=dft.and.infos%dft%cam_flag, &
              alpha=infos%dft%cam_alpha, &
              beta=infos%dft%cam_beta,&
              mu=infos%dft%cam_mu)
      ab1 => int2_data%apb(:,:,:,1)
      ab2 => int2_data%amb(:,:,:,1)

      pa = pa*2
      call utddft_fxc( &
          basis = basis, &
          molGrid = molGrid, &
          isVecs = .true., &
          wfa = mo_a, &
          wfb = mo_b, &
          fxa = ab1(:,:,1:1), &
          fxb = ab1(:,:,2:2), &
          dxa = pa(:,:,1:1), &
          dxb = pa(:,:,2:2), &
          nmtx = 1, &
          threshold = 1.0d-15, &
          infos = infos)

      !   ALPHA: AO(M,N) -> MO(IA+) ... LPTMOA
      call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)

      wrk1 = 2*ab1(:,:,2) + ab2(:,:,2)
      call mntoia(wrk1, ab1_mo_b, mo_b, mo_b, noccb, noccb)

      call sfrolhs(lhs, xk, mo_energy_a, fa, fb, ab1_mo_a, ab1_mo_b, &
                   nocca, noccb)

      call pcgrbpini(errv, pk, error, rhs, xminv, lhs)

      write(iw,'(" Initial error =",3x,1p,e10.3,1x,"/",1p,e10.3)') error, cnvtol
      call flush(iw)

      ! -----------------------------------------------

      do iter = 1, infos%control%maxit_zv

        call sfrogen(wrk1, wrk2, pk, nocca, noccb)
        !     Alpha
        call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)
        !     Beta
        call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)

        !     (A+B)*PK
        call int2_data%clean()
        deallocate(int2_data)
        int2_data = int2_tdgrd_data_t( &
            d2 = pa, &
            int_apb = .true., &
            int_amb = .false., &
            tamm_dancoff = .false., &
            scale_exchange = scale_exch)

        call int2_driver%run(int2_data, &
              cam=dft.and.infos%dft%cam_flag, &
              alpha=infos%dft%cam_alpha, &
              beta=infos%dft%cam_beta,&
              mu=infos%dft%cam_mu)
        ab1 => int2_data%apb(:,:,:,1)

        !ab1 = ab1/2
        call symmetrize_matrix(pa(:,:,1), nbf)
        call symmetrize_matrix(pa(:,:,2), nbf)
        call utddft_fxc( &
            basis = basis, &
            molGrid = molGrid, &
            isVecs = .true., &
            wfa = mo_a, &
            wfb = mo_b, &
            fxa = ab1(:,:,1:1), &
            fxb = ab1(:,:,2:2), &
            dxa = pa(:,:,1:1), &
            dxb = pa(:,:,2:2), &
            nmtx = 1, &
            threshold = 1.0d-15, &
            infos = infos)

        !     ALPHA: AO(M,N) -> MO(IA+) ... LPTMOA
        call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)

        call mntoia(ab1(:,:,2), ab1_mo_b, mo_b, mo_b, noccb, noccb)

        call sfrolhs(lhs, pk, mo_energy_a, fa, fb, ab1_mo_a, ab1_mo_b, &
                     nocca, noccb)

        alpha = 1.0_dp/dot_product(pk, lhs)

        xk = xk + pk * alpha
        errv = errv - alpha*lhs

        error = dot_product(errv, errv)
        write(iw,'(" Iter#",I2," Error =",&
              &3x,1p,e10.3,1x,"/",1p,e10.3)') &
                iter, error, cnvtol
        call flush(iw)

        if (error<cnvtol) exit

        call pcgb(pk, errv, xminv)

      end do
      
    end if  ! End solver selection

! -----------------------------------------------
    if (error>cnvtol) then
       infos%mol_energy%Z_Vector_converged=.false.
       write(*,'(/3x,24("-")&
             &/6x,"Z-Vector not converged"&
             &/3x,24("-")/)')
    else
       infos%mol_energy%Z_Vector_converged=.true.
       write(*,'(/3x,24("-")&
             &/6x,"Z-Vector converged"&
             &/3x,24("-")/)')
    endif

    call flush(iw)

    if (mrst==1 .or. mrst==3) then

      call sfropcal(wrk1, wrk2, tij, tab, xk, nocca, noccb)

    else if (mrst==5) then

      call mrsfqropcal(wrk1, wrk2, tab, tij, xk, nocca, noccb)

    end if

 !  Update density for alpha
    call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)

 !  Update density for beta
    call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)
    call int2_data%clean()
    deallocate(int2_data)
    int2_data = int2_tdgrd_data_t( &
        d2 = pa, &
        int_apb = .true., &
        int_amb = .false., &
        tamm_dancoff = .false., &
        scale_exchange = scale_exch)

    call int2_driver%run(int2_data, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%dft%cam_alpha, &
            beta=infos%dft%cam_beta,&
            mu=infos%dft%cam_mu)
    ab1 => int2_data%apb(:,:,:,1)

    call symmetrize_matrix(pa(:,:,1), nbf)
    call symmetrize_matrix(pa(:,:,2), nbf)
    call pack_matrix(pa(:,:,1), td_p(:,1))
    call pack_matrix(pa(:,:,2), td_p(:,2))

    td_p = 0.5_dp*td_p

    call utddft_fxc( &
        basis = basis, &
        molGrid = molGrid, &
        isVecs = .true., &
        wfa = mo_a, &
        wfb = mo_b, &
        fxa = ab1(:,:,1:1), &
        fxb = ab1(:,:,2:2), &
        dxa = pa(:,:,1:1), &
        dxb = pa(:,:,2:2), &
        nmtx = 1, &
        threshold = 1.0d-15, &
        infos = infos)

!   ALPHA AO(M,N) -> MO(I-,J-) ... LPPIJA
    call dgemm('n', 'n', nbf, nocca, nbf,  &
               1.0_dp, ab1(:,:,1), nbf,  &
                       mo_a, nbf,  &
               0.0_dp, wrk2, nbf)
    call dgemm('t', 'n', nocca, nocca, nbf,  &
               1.0_dp, mo_a,  nbf,  &
                       wrk2,  nbf,  &
               0.0_dp, ppija, nocca)
!   BETA: AO(M,N) -> MO(I-,J-) ... LPPIJB
    call dgemm('n', 'n', nbf, noccb, nbf,  &
               1.0_dp, ab1(:,:,2), nbf,  &
                       mo_b, nbf,  &
               0.0_dp, wrk2, nbf)
    call dgemm('t', 'n', noccb, noccb, nbf,  &
               1.0_dp, mo_b,  nbf,  &
                       wrk2,  nbf,  &
               0.0_dp, ppijb, noccb)

!   Calculate W (in MO basis)
    wmo => wrk3
    wmo = 0
    if (mrst==1 .or. mrst==3) then

      call mrsfrowcal(wmo, mo_energy_a, fa, fb, xk, &
                      hxa, hxb, ppija, ppijb, &
                      nocca, noccb)

    else if (mrst==5) then

      call mrsfqrowcal(wmo, mo_energy_a, fa, fb, xk, &
                       hxa, hxb, ppija, ppijb, &
                       nocca, noccb)

    end if

    call orthogonal_transform('t', nbf, mo_a, wmo, wrk2, wrk1)
    call symmetrize_matrix(wrk2, nbf)
    call pack_matrix(wrk2, wao)
    wao = wao*0.5_dp
!   ROHF, half one more time:
    wao = wao*0.5_dp

    call int2_driver%clean()

    if (dft) call dftclean(infos)

    call measure_time(print_total=1, log_unit=iw)
    ! Clean up GMRES work arrays
    call cleanup_gmres_work()

    close(iw)

  contains

    ! Lambda wrapper for preconditioner
    subroutine lambda_precond(x_in, x_out)
      real(kind=dp), intent(in) :: x_in(:)
      real(kind=dp), intent(out) :: x_out(:)
      call apply_z_precond(x_in, x_out, xminv)
    end subroutine lambda_precond

  end subroutine tdhf_mrsf_z_vector

!###############################################################################
!  UMRSF Q/R test driver (Phase 2). NOT a production path.
!  Builds, for the target state, the UMRSF Z-vector RHS R (via umrsfqcal +
!  umrsfqrorhs) AND the golden RO reference R (via sfrorhs), and dumps both plus
!  Q_alpha/Q_beta to a text file `umrsf_qr_dump.txt`. In the RO limit (mo_a==mo_b,
!  ROHF reference) the two R's must agree per block (Check 1 / RO-limit test).
!  Theory: umrsf_gradient_theory.tex §8 (Q), §9 (R).
!###############################################################################
  subroutine tdhf_umrsf_qrtest_C(c_handle) bind(C, name="tdhf_umrsf_qrtest")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_umrsf_qrtest(inf)
  end subroutine tdhf_umrsf_qrtest_C

  subroutine tdhf_umrsf_qrtest(infos)
    use precision, only: dp
    use io_constants, only: iw
    use oqp_tagarray_driver
    use types, only: information
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use int2_compute, only: int2_compute_t
    use tdhf_lib, only: int2_td_data_t, int2_tdgrd_data_t, iatogen, mntoia
    use tdhf_mrsf_lib, only: int2_mrsf_data_t, int2_umrsf_data_t, &
      mrsfcbc, umrsfcbc, mrsfsp, umrsfsp, mrsfxvec, umrsfqcal, umrsfqrorhs, umrsfdmat, &
      umrsf_sfropcal, umrsf_sfrowcal, mrsfrowcal, umrsf_sfrogen, umrsf_dbg_jswap, umrsf_dbg_jswap2, &
      umrsf_dbg_wswap
    use tdhf_sf_lib, only: sfrorhs, sfdmat, sfropcal
    use tdhf_mrsf_gradient_mod, only: mrsf_2e_grad, umrsf_2e_grad, umrsf_dbg_nodt2, umrsf_dbg_zero
    use tdhf_sf_gradient_mod, only: sf_1e_grad
    use mod_dft_gridint_tdxc_grad, only: utddft_xc_gradient
    use mathlib, only: pack_matrix
    use dft, only: dft_initialize, dftclean
    use mod_dft_gridint_fxc, only: utddft_fxc
    use mathlib, only: orthogonal_transform, orthogonal_transform_sym, &
      unpack_matrix, symmetrize_matrix
    use mod_dft_molgrid, only: dft_grid_t
    use printing, only: print_module_info

    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_umrsf_qrtest"
    type(information), target, intent(inout) :: infos
    type(basis_set), pointer :: basis

    integer :: nbf, nbf_tri, nocca, noccb, nvira, nvirb, nsocc
    integer :: lzdim_u, lzdim_ro, mrst, target_state, ok, iu
    real(kind=dp) :: scale_exch, scale_exch2
    logical :: dft

    type(int2_compute_t) :: int2_driver
    type(int2_umrsf_data_t), allocatable, target :: int2_udata_st
    type(int2_mrsf_data_t),  allocatable, target :: int2_data_st
    class(int2_td_data_t), allocatable, target :: int2_data
    type(dft_grid_t) :: molGrid
    real(kind=dp), pointer :: ab1(:,:,:), fmrst2u(:,:,:,:), fmrst2r(:,:,:,:)

    real(kind=dp), allocatable :: fa(:,:), fb(:,:)
    real(kind=dp), allocatable :: bvec_mo_d(:), xmo(:,:)
    real(kind=dp), allocatable :: tij(:,:), tab(:,:)
    real(kind=dp), allocatable, target :: pa(:,:,:)
    real(kind=dp), allocatable :: hpta(:,:), hptb(:,:)
    real(kind=dp), allocatable :: ab1_mo_a(:,:), ab1_mo_b(:,:)
    real(kind=dp), allocatable :: hxa_u(:,:), hxb_u(:,:), hxa_r(:,:), hxb_r(:,:)
    real(kind=dp), allocatable :: spxa_u(:,:), spxb_u(:,:), spxa_r(:,:), spxb_r(:,:)
    real(kind=dp), allocatable :: h0a_u(:,:), h0b_u(:,:), h0a_r(:,:), h0b_r(:,:)
    real(kind=dp), allocatable :: qa(:,:), qb(:,:)
    real(kind=dp), allocatable :: rhs_u(:), rhs_ro(:)
    real(kind=dp), allocatable :: ju(:,:), jro(:,:), evec(:), xout(:)
    integer :: jc
    real(kind=dp), allocatable :: zu(:), zro(:), zuc(:)
    real(kind=dp), allocatable :: pau(:,:), pbu(:,:), pauc(:,:), pbuc(:,:), paro(:,:), pbro(:,:)
    real(kind=dp), allocatable :: ppija(:,:), ppijb(:,:), xhxa(:,:), xhxb(:,:)
    real(kind=dp), allocatable :: xhxa_u(:,:), xhxb_u(:,:), xhua(:,:), xhub(:,:)  ! UMRSF H[X,X]+Fock.T probe
    logical :: use_hxu
    real(kind=dp), allocatable :: hppmoa(:,:), hppmob(:,:)   ! full MO H+[P] (occ-virt term)
    real(kind=dp), allocatable :: wa_uc(:,:), wb_uc(:,:), wmo_ro(:,:), tbig(:,:)
    real(kind=dp), allocatable :: wa_u(:,:), wb_u(:,:)
    real(kind=dp), allocatable :: cab(:,:)   ! R1: beta cross-response to OV_alpha pert
    real(kind=dp) :: gmres_zdiff, gmres_resid
    integer :: gmres_iters
    real(kind=dp), allocatable :: de_u(:,:), de_r(:,:)
    real(kind=dp), allocatable, target :: d2g(:,:,:), pg(:,:,:), spca(:,:,:), spcb_(:,:,:), spcr(:,:,:)
    real(kind=dp), allocatable, target :: vdum(:,:)
    real(kind=dp), allocatable :: g1e_u(:,:), gxc_u(:,:), gtot_u(:,:)
    real(kind=dp), allocatable :: g1e_r(:,:), gxc_r(:,:), gtot_r(:,:)
    real(kind=dp), allocatable :: g1e_x0(:,:), gxc_x0(:,:), gtot_x0(:,:), gtot_t0(:,:), gtot_probe(:,:)
    real(kind=dp), pointer :: wao_p(:), tdp_p(:,:)
    integer :: natom
    real(kind=dp), allocatable, target :: densu(:,:,:,:), densr(:,:,:,:)
    real(kind=dp), allocatable :: wrk1(:,:), wrk2(:,:), wrk3(:,:), pmo(:,:)
    real(kind=dp), allocatable :: abxc(:,:), ta_p(:), tb_p(:)

    ! tagarray
    real(kind=dp), contiguous, pointer :: &
      fock_a(:), fock_b(:), mo_a(:,:), mo_b(:,:), mo_energy_a(:), mo_energy_b(:), bvec_mo(:,:)
    character(len=*), parameter :: tags_required(7) = (/ character(len=80) :: &
      OQP_FOCK_A, OQP_FOCK_B, OQP_E_MO_A, OQP_E_MO_B, OQP_VEC_MO_A, OQP_VEC_MO_B, OQP_td_bvec_mo /)

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf_tri = nbf*(nbf+1)/2
    dft = infos%control%hamilton == 20
    mrst = infos%tddft%mult

    nocca = infos%mol_prop%nelec_A
    noccb = infos%mol_prop%nelec_B
    nvira = nbf - nocca
    nvirb = nbf - noccb
    nsocc = nocca - noccb
    target_state = min(infos%tddft%target_state, infos%tddft%nstate)

    lzdim_ro = noccb*(nsocc+nvira) + nsocc*nvira          ! RO: CO + CV + OV
    lzdim_u  = nsocc*nvira + noccb*nvira + noccb*nsocc + noccb*nvira ! OVa+CVa+CObeta+CVbeta

    open (unit=iw, file=infos%log_filename, position="append")
    call print_module_info('UMRSF_QRtest','UMRSF Q/R builder test (Phase 2)')

    if (dft) call dft_initialize(infos, basis, molGrid)
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()

    scale_exch = 1.0_dp; scale_exch2 = 1.0_dp
    if (dft) then
      scale_exch  = infos%dft%HFscale
      scale_exch2 = infos%tddft%HFscale
    end if

    call data_has_tags(infos%dat, tags_required, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
    call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_E_MO_B, mo_energy_b)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)

    allocate(fa(nbf,nbf), fb(nbf,nbf), &
             bvec_mo_d(nocca*nvirb), xmo(nocca,nvirb), &
             tij(nocca,nocca), tab(nvirb,nvirb), pa(nbf,nbf,2), &
             hpta(nbf,nbf), hptb(nbf,nbf), &
             ab1_mo_a(nocca,nvira), ab1_mo_b(noccb,nvirb), &
             hxa_u(nbf,nbf), hxb_u(nbf,nbf), hxa_r(nbf,nbf), hxb_r(nbf,nbf), &
             spxa_u(nbf,nbf), spxb_u(nbf,nbf), spxa_r(nbf,nbf), spxb_r(nbf,nbf), &
             h0a_u(nbf,nbf), h0b_u(nbf,nbf), h0a_r(nbf,nbf), h0b_r(nbf,nbf), &
             qa(nbf,nbf), qb(nbf,nbf), rhs_u(lzdim_u), rhs_ro(lzdim_ro), &
             densu(1,11,nbf,nbf), densr(1,7,nbf,nbf), &
             wrk1(nbf,nbf), wrk2(nbf,nbf), wrk3(nbf,nbf), pmo(nbf,nbf), &
             abxc(nbf,nbf), ta_p(nbf_tri), tb_p(nbf_tri), &
             source=0.0_dp, stat=ok)
    if (ok/=0) call show_message('qrtest: cannot allocate', with_abort)

  ! UKS Fock in MO basis (alpha from mo_a, beta from mo_b; RO limit: mo_b==mo_a)
    call uks_fock_mo(fock_a, mo_a, fa, wrk1, nbf)
    call uks_fock_mo(fock_b, mo_b, fb, wrk1, nbf)

  ! Dimensional-transformed transition amplitude U.X  (nocca x nvirb)
    call mrsfxvec(infos, bvec_mo(:,target_state), bvec_mo_d)
    xmo = reshape(bvec_mo_d, (/ nocca, nvirb /))
    call iatogen(bvec_mo_d, wrk3, nocca, noccb)   ! X density (nbf x nbf), MO basis

  ! Unrelaxed difference density T (eq:Talpha/Tbeta)
  !   Ta(i,j) = -sum_a X(i,a) X(j,a);  Tb(a,b) = +sum_i X(i,a) X(i,b)
    call dgemm('n','t', nocca, nocca, nvirb, -1.0_dp, bvec_mo_d, nocca, &
               bvec_mo_d, nocca, 0.0_dp, tij, nocca)
    call dgemm('t','n', nvirb, nvirb, nocca,  1.0_dp, bvec_mo_d, nocca, &
               bvec_mo_d, nocca, 0.0_dp, tab, nvirb)

  ! AO T-density: pa(:,:,1) = mo_a * Ta * mo_a^T ; pa(:,:,2) = mo_b * Tb * mo_b^T
    pmo = 0.0_dp; pmo(1:nocca,1:nocca) = tij
    call orthogonal_transform('t', nbf, mo_a, pmo, pa(:,:,1), wrk1)
    pmo = 0.0_dp; pmo(noccb+1:nbf,noccb+1:nbf) = tab
    call orthogonal_transform('t', nbf, mo_b, pmo, pa(:,:,2), wrk1)

  ! H^+[T]: (A+B) operator on T-density + fxc, then AO->MO (full nbf x nbf)
    int2_data = int2_tdgrd_data_t(d2=pa, int_apb=.true., int_amb=.false., &
                                  tamm_dancoff=.false., scale_exchange=scale_exch)
    call int2_driver%run(int2_data, cam=dft.and.infos%dft%cam_flag, &
        alpha=infos%dft%cam_alpha, beta=infos%dft%cam_beta, mu=infos%dft%cam_mu)
    ab1 => int2_data%apb(:,:,:,1)
    pa = pa*2.0_dp
    if (dft) call utddft_fxc(basis=basis, molGrid=molGrid, isVecs=.true., &
        wfa=mo_a, wfb=mo_b, fxa=ab1(:,:,1:1), fxb=ab1(:,:,2:2), &
        dxa=pa(:,:,1:1), dxb=pa(:,:,2:2), nmtx=1, threshold=1.0d-15, infos=infos)
    call orthogonal_transform('n', nbf, mo_a, ab1(:,:,1), hpta, wrk1)
    call orthogonal_transform('n', nbf, mo_b, ab1(:,:,2), hptb, wrk1)
  ! RO-style occ x vir H^+[T] for sfrorhs reference
    call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)
    call mntoia(ab1(:,:,2), ab1_mo_b, mo_b, mo_b, noccb, noccb)

  ! ============ UMRSF H[X,X] (hxa_u, hxb_u) ============
  ! A0 SF transition density from the U-transformed amplitude (umrsfdmat), placed in
  ! the A0 slot (component 11) - mirrors the RO override densr(7)=sfdmat(bvec_mo_d).
    call umrsfdmat(bvec_mo_d, abxc, mo_a, mo_b, ta_p, tb_p, nocca, noccb)
    call umrsfcbc(infos, mo_a, mo_b, wrk3, densu(1,:,:,:))
    densu(1,11,:,:) = abxc
    int2_udata_st = int2_umrsf_data_t(d3=densu, tamm_dancoff=.true., &
        scale_exchange=scale_exch2, scale_coulomb=scale_exch2)
    call int2_driver%run(int2_udata_st, cam=dft.and.infos%dft%cam_flag, &
        alpha=infos%tddft%cam_alpha, alpha_coulomb=infos%tddft%cam_alpha, &
        beta=infos%tddft%cam_beta, beta_coulomb=infos%tddft%cam_beta, mu=infos%tddft%cam_mu)
    fmrst2u => int2_udata_st%f3(:,:,:,:,1)
    if (mrst==3) fmrst2u(:,1:10,:,:) = -fmrst2u(:,1:10,:,:)
  ! H0 part: A0 SF Fock (component 11) -> H0 in dedicated h0a_u (alpha), h0b_u (beta).
    call orthogonal_transform('n', nbf, mo_a, fmrst2u(1,11,:,:), wrk1, wrk2)
    call dgemm('n','t', nbf, nocca, nbf, 2.0_dp, wrk1, nbf, wrk3, nbf, 0.0_dp, h0a_u, nbf)
    call orthogonal_transform('n', nbf, mo_b, fmrst2u(1,11,:,:), wrk2, wrk1)
    call dgemm('t','n', nbf, nbf, nocca, 2.0_dp, wrk2, nbf, wrk3, nbf, 0.0_dp, h0b_u, nbf)
  ! spin-pairing part: compute into zeroed scratch (umrsfsp/mrsfsp declare xhxa/xhxb
  ! intent(out)) then assemble H[X,X] = H0 + sp by explicit sum into fresh hxa_u/hxb_u.
  ! Keeping H0 and sp in dedicated arrays avoids an in-place "hx = hx + sp" whose lhs is
  ! held live across the (intent(out)) sp call - observed to be miscompiled/corrupted.
    spxa_u = 0.0_dp; spxb_u = 0.0_dp
    call umrsfsp(spxa_u, spxb_u, mo_a, mo_b, wrk3, fmrst2u(1,:,:,:), nocca, noccb)
    hxa_u = h0a_u + spxa_u; hxb_u = h0b_u + spxb_u

  ! ============ RO reference H[X,X] (hxa_r, hxb_r) ============
    abxc = 0.0_dp
    call sfdmat(bvec_mo_d, abxc, mo_a, ta_p, tb_p, nocca, noccb)
    call mrsfcbc(infos, mo_a, mo_a, wrk3, densr(1,:,:,:))
    densr(1,7,:,:) = abxc
    int2_data_st = int2_mrsf_data_t(d3=densr, tamm_dancoff=.true., &
        scale_exchange=scale_exch2, scale_coulomb=scale_exch2)
    call int2_driver%run(int2_data_st, cam=dft.and.infos%dft%cam_flag, &
        alpha=infos%tddft%cam_alpha, alpha_coulomb=infos%tddft%cam_alpha, &
        beta=infos%tddft%cam_beta, beta_coulomb=infos%tddft%cam_beta, mu=infos%tddft%cam_mu)
    fmrst2r => int2_data_st%f3(:,:,:,:,1)
    if (mrst==3) fmrst2r(:,1:6,:,:) = -fmrst2r(:,1:6,:,:)
    call orthogonal_transform('n', nbf, mo_a, fmrst2r(1,7,:,:), wrk1, wrk2)
    call dgemm('n','t', nbf, nocca, nbf, 2.0_dp, wrk1, nbf, wrk3, nbf, 0.0_dp, h0a_r, nbf)
    call dgemm('t','n', nbf, nbf, nocca, 2.0_dp, wrk1, nbf, wrk3, nbf, 0.0_dp, h0b_r, nbf)
    spxa_r = 0.0_dp; spxb_r = 0.0_dp
    call mrsfsp(spxa_r, spxb_r, mo_a, mo_a, wrk3, fmrst2r(1,:,:,:), nocca, noccb)
    hxa_r = h0a_r + spxa_r; hxb_r = h0b_r + spxb_r

  ! ============ UMRSF Q and R ============
    call umrsfqcal(qa, qb, hpta, hptb, hxa_u, hxb_u, fa, fb, tij, tab, nocca, noccb)
    call umrsfqrorhs(rhs_u, qa, qb, nocca, noccb)

  ! ============ RO reference R via golden sfrorhs ============
  ! sfrorhs takes xhxa/xhxb intent(inout) and accumulates the Fock.T term in place
  ! (xhxa += 2*Fa*Tij): pass copies so the dumped hxa_r stays the bare H[X,X].
    spxa_r = hxa_r; spxb_r = hxb_r
    call sfrorhs(rhs_ro, spxa_r, spxb_r, ab1_mo_a, ab1_mo_b, tij, tab, fa, fb, nocca, noccb)

  ! ============ Phase 3: dense J operators (column-by-column J . e_i) ============
  ! umrsf 4-block J (lzdim_u) and RO 3-block J (lzdim_ro), for symmetry / alpha=beta
  ! perturbation tests and a dense LAPACK solve of J.Z = -R.
    block                                        ! probe hook: J same-spin Fock swap
      character(len=8) :: envv; integer :: envl, envs
      call get_environment_variable('UMRSF_JSWAP', envv, envl, envs)
      umrsf_dbg_jswap = 0
      if (envs == 0 .and. envl > 0) read(envv,*) umrsf_dbg_jswap
      call get_environment_variable('UMRSF_JSWAP2', envv, envl, envs)
      umrsf_dbg_jswap2 = 0
      if (envs == 0 .and. envl > 0) read(envv,*) umrsf_dbg_jswap2
      call get_environment_variable('UMRSF_WSWAP', envv, envl, envs)
      umrsf_dbg_wswap = 0
      if (envs == 0 .and. envl > 0) read(envv,*) umrsf_dbg_wswap
    end block
    allocate(ju(lzdim_u,lzdim_u), jro(lzdim_ro,lzdim_ro), &
             evec(max(lzdim_u,lzdim_ro)), xout(max(lzdim_u,lzdim_ro)), source=0.0_dp)
    do jc = 1, lzdim_u
      evec(1:lzdim_u) = 0.0_dp; evec(jc) = 1.0_dp
      call umrsf_apply_z_operator(evec(1:lzdim_u), xout(1:lzdim_u), infos, basis, &
           molGrid, int2_driver, nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, mo_energy_b, &
           fa, fb, scale_exch, dft)
      ju(:,jc) = xout(1:lzdim_u)
    end do
    do jc = 1, lzdim_ro
      evec(1:lzdim_ro) = 0.0_dp; evec(jc) = 1.0_dp
      call apply_z_operator(evec(1:lzdim_ro), xout(1:lzdim_ro), infos, basis, &
           molGrid, int2_driver, nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
           fa, fb, scale_exch, dft)
      jro(:,jc) = xout(1:lzdim_ro)
    end do

  ! ============ Phase 4: solve Z (dense LAPACK), build P = T + Z ============
    allocate(zu(lzdim_u), zro(lzdim_ro), zuc(lzdim_u), &
             pau(nbf,nbf), pbu(nbf,nbf), pauc(nbf,nbf), pbuc(nbf,nbf), &
             paro(nbf,nbf), pbro(nbf,nbf), source=0.0_dp)
    block
      real(kind=dp), allocatable :: jcp(:,:)
      integer, allocatable :: ipiv(:)
      integer :: info, o1, o2, o3
      allocate(ipiv(max(lzdim_u,lzdim_ro)))
      ! J.Z = -R ; rhs_u stores -R already, so solve J.Z = rhs_u.
      jcp = ju;  zu = rhs_u
      call dgesv(lzdim_u, 1, jcp, lzdim_u, ipiv, zu, lzdim_u, info)
      deallocate(jcp); jcp = jro; zro = rhs_ro
      call dgesv(lzdim_ro, 1, jcp, lzdim_ro, ipiv, zro, lzdim_ro, info)
      ! controlled zuc: map RO 3-block (CO,CV,OV) -> umrsf 4-block (OVa,CVa,COb,CVb)
      o1 = noccb*nsocc; o2 = noccb*nvira; o3 = nsocc*nvira     ! CO, CV, OV sizes
      zuc(1:o3)               = zro(o1+o2+1 : o1+o2+o3)        ! OVa <- OV
      zuc(o3+1 : o3+o2)       = zro(o1+1 : o1+o2)              ! CVa <- CV
      zuc(o3+o2+1 : o3+o2+o1) = zro(1 : o1)                    ! COb <- CO
      zuc(o3+o2+o1+1 : )      = zro(o1+1 : o1+o2)              ! CVb <- CV
    end block

  ! ---- GMRES production solve of J.Z = rhs_u (= -R) + equivalence vs dense zu ----
    umrsf_mo_energy_b => mo_energy_b    ! beta orbital energies (UHF: != mo_energy_a)
    block
      real(kind=dp), allocatable :: zg(:)
      real(kind=dp) :: gerr, zdiff(1)
      integer :: giter
      allocate(zg(lzdim_u), source=0.0_dp)
      call gmres_solve(umrsf_apply_z_op_gmres, ident_precond, rhs_u, zg, lzdim_u, &
           min(infos%tddft%gmres_dim, lzdim_u), infos%control%maxit_zv, infos%tddft%zvconv, &
           infos, basis, molGrid, int2_driver, nocca, noccb, nbf, mo_a, mo_b, mo_energy_a, &
           fa, fb, scale_exch, dft, gerr, giter, iw)
      gmres_zdiff = maxval(abs(zg - zu)); gmres_resid = gerr; gmres_iters = giter
    end block

    block   ! probe: SCALE a Z-block of zu by UMRSF_ZFAC (env UMRSF_ZBLOCK 1=OVa 2=CVa 3=COb
            ! 4=CVb; UMRSF_ZFAC default 0.0=zero) to localize/characterize a gradient component.
      character(len=16) :: envv; integer :: envl, envs, zblk, b1, b2, b3
      real(kind=dp) :: zfac
      call get_environment_variable('UMRSF_ZBLOCK', envv, envl, envs)
      zblk = 0; if (envs == 0 .and. envl > 0) read(envv,*) zblk
      zfac = 0.0_dp
      call get_environment_variable('UMRSF_ZFAC', envv, envl, envs)
      if (envs == 0 .and. envl > 0) read(envv,*) zfac
      b1 = nsocc*nvira; b2 = b1 + noccb*nvira; b3 = b2 + noccb*nsocc
      select case (zblk)
        case (1); zu(1:b1)      = zu(1:b1)      * zfac
        case (2); zu(b1+1:b2)   = zu(b1+1:b2)   * zfac
        case (3); zu(b2+1:b3)   = zu(b2+1:b3)   * zfac
        case (4); zu(b3+1:)     = zu(b3+1:)     * zfac
      end select
    end block
    call umrsf_sfropcal(pau,  pbu,  tij, tab, zu,  nocca, noccb)   ! real relaxed P
    call umrsf_sfropcal(pauc, pbuc, tij, tab, zuc, nocca, noccb)   ! controlled P
    call sfropcal(paro, pbro, tij, tab, zro, nocca, noccb)         ! RO golden P

  ! ============ Phase 4: H+[P], W builder (spin-resolved) + RO fold ============
    allocate(ppija(nocca,nocca), ppijb(noccb,noccb), xhxa(nbf,nbf), xhxb(nbf,nbf), &
             hppmoa(nbf,nbf), hppmob(nbf,nbf), &
             wa_uc(nbf,nbf), wb_uc(nbf,nbf), wmo_ro(nbf,nbf), tbig(nbf,nbf), source=0.0_dp)
  ! xhxa/xhxb = H[X,X] + Fock.T  (= Q without the masked H+[T])
    xhxa = hxa_r; xhxb = hxb_r
    tbig = 0.0_dp; tbig(1:nocca,1:nocca) = tij
    call dgemm('n','n', nbf, nbf, nbf, 2.0_dp, fa, nbf, tbig, nbf, 1.0_dp, xhxa, nbf)
    tbig = 0.0_dp; tbig(noccb+1:nbf,noccb+1:nbf) = tab
    call dgemm('n','n', nbf, nbf, nbf, 2.0_dp, fb, nbf, tbig, nbf, 1.0_dp, xhxb, nbf)
  ! PROBE: xhxa_u/xhxb_u use the UMRSF H[X,X] (hxa_u/hxb_u) instead of RO (hxa_r/hxb_r). The
  ! UMRSF W must be consistent with the R/Q which already uses hxa_u (line ~1611). RO-safe
  ! (hxa_u==hxa_r at alpha=beta). env UMRSF_HXU=1 routes the UMRSF W through xhxa_u.
    allocate(xhxa_u(nbf,nbf), xhxb_u(nbf,nbf), xhua(nbf,nbf), xhub(nbf,nbf))
    xhxa_u = hxa_u; xhxb_u = hxb_u
    tbig = 0.0_dp; tbig(1:nocca,1:nocca) = tij
    call dgemm('n','n', nbf, nbf, nbf, 2.0_dp, fa, nbf, tbig, nbf, 1.0_dp, xhxa_u, nbf)
    tbig = 0.0_dp; tbig(noccb+1:nbf,noccb+1:nbf) = tab
    call dgemm('n','n', nbf, nbf, nbf, 2.0_dp, fb, nbf, tbig, nbf, 1.0_dp, xhxb_u, nbf)
    block
      character(len=8) :: ev; integer :: el, es
      call get_environment_variable('UMRSF_HXU', ev, el, es)
      use_hxu = (es == 0 .and. el > 0)
    end block
    if (use_hxu) then
      xhua = xhxa_u; xhub = xhxb_u
    else
      xhua = xhxa;   xhub = xhxb
    end if
  ! H+[P] (occ-occ) on the RO relaxed density (== controlled P at alpha=beta)
    block
      type(int2_tdgrd_data_t), allocatable, target :: i2p
      real(kind=dp), pointer :: abp(:,:,:)
      call orthogonal_transform('t', nbf, mo_a, paro, pa(:,:,1), wrk1)
      call orthogonal_transform('t', nbf, mo_b, pbro, pa(:,:,2), wrk1)
      allocate(i2p)
      i2p = int2_tdgrd_data_t(d2=pa, int_apb=.true., int_amb=.false., &
            tamm_dancoff=.false., scale_exchange=scale_exch)
      call int2_driver%run(i2p, cam=dft.and.infos%dft%cam_flag, &
          alpha=infos%dft%cam_alpha, beta=infos%dft%cam_beta, mu=infos%dft%cam_mu)
      abp => i2p%apb(:,:,:,1)
      call symmetrize_matrix(pa(:,:,1), nbf); call symmetrize_matrix(pa(:,:,2), nbf)
      if (dft) call utddft_fxc(basis=basis, molGrid=molGrid, isVecs=.true., wfa=mo_a, wfb=mo_b, &
          fxa=abp(:,:,1:1), fxb=abp(:,:,2:2), dxa=pa(:,:,1:1), dxb=pa(:,:,2:2), &
          nmtx=1, threshold=1.0d-15, infos=infos)
      call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, abp(:,:,1), nbf, mo_a, nbf, 0.0_dp, wrk2, nbf)
      call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, mo_a, nbf, wrk2, nbf, 0.0_dp, hppmoa, nbf)
      call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, abp(:,:,2), nbf, mo_b, nbf, 0.0_dp, wrk2, nbf)
      call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, mo_b, nbf, wrk2, nbf, 0.0_dp, hppmob, nbf)
      ppija = hppmoa(1:nocca,1:nocca); ppijb = hppmob(1:noccb,1:noccb)
      call i2p%clean(); deallocate(i2p)
    end block
    call umrsf_sfrowcal(wa_uc, wb_uc, mo_energy_a, mo_energy_b, fa, fb, zuc, &
                        xhua, xhub, ppija, ppijb, hppmoa, hppmob, nocca, noccb)
    call mrsfrowcal(wmo_ro, mo_energy_a, fa, fb, zro, xhxa, xhxb, ppija, ppijb, nocca, noccb)

  ! Real W from the solved Zu (needs H+[P(Zu)] on the real relaxed density) -> Phase 5
    allocate(wa_u(nbf,nbf), wb_u(nbf,nbf), source=0.0_dp)
    block
      type(int2_tdgrd_data_t), allocatable, target :: i2p
      real(kind=dp), pointer :: abp(:,:,:)
      call orthogonal_transform('t', nbf, mo_a, pau, pa(:,:,1), wrk1)
      call orthogonal_transform('t', nbf, mo_b, pbu, pa(:,:,2), wrk1)
      allocate(i2p)
      i2p = int2_tdgrd_data_t(d2=pa, int_apb=.true., int_amb=.false., &
            tamm_dancoff=.false., scale_exchange=scale_exch)
      call int2_driver%run(i2p, cam=dft.and.infos%dft%cam_flag, &
          alpha=infos%dft%cam_alpha, beta=infos%dft%cam_beta, mu=infos%dft%cam_mu)
      abp => i2p%apb(:,:,:,1)
      call symmetrize_matrix(pa(:,:,1), nbf); call symmetrize_matrix(pa(:,:,2), nbf)
      if (dft) call utddft_fxc(basis=basis, molGrid=molGrid, isVecs=.true., wfa=mo_a, wfb=mo_b, &
          fxa=abp(:,:,1:1), fxb=abp(:,:,2:2), dxa=pa(:,:,1:1), dxb=pa(:,:,2:2), &
          nmtx=1, threshold=1.0d-15, infos=infos)
      call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, abp(:,:,1), nbf, mo_a, nbf, 0.0_dp, wrk2, nbf)
      call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, mo_a, nbf, wrk2, nbf, 0.0_dp, hppmoa, nbf)
      call dgemm('n','n', nbf, nbf, nbf, 1.0_dp, abp(:,:,2), nbf, mo_b, nbf, 0.0_dp, wrk2, nbf)
      call dgemm('t','n', nbf, nbf, nbf, 1.0_dp, mo_b, nbf, wrk2, nbf, 0.0_dp, hppmob, nbf)
      ppija = hppmoa(1:nocca,1:nocca); ppijb = hppmob(1:noccb,1:noccb)
      call i2p%clean(); deallocate(i2p)
    end block
    call umrsf_sfrowcal(wa_u, wb_u, mo_energy_a, mo_energy_b, fa, fb, zu, &
                        xhua, xhub, ppija, ppijb, hppmoa, hppmob, nocca, noccb)

  ! ============ R1 adjudication: beta cross-response H~^{beta,alpha}[Z_OV_alpha] =======
  ! Pure OV_alpha perturbation -> alpha density only; the beta response from the shared
  ! int2(A+B)+f^xc is the cross-spin H~^{beta,alpha}. Dump in full MO so the (i,a) vs
  ! (a,i) blocks give C_iaβ(OV) and C_aiβ(OV); their ratio adjudicates R1 (b) vs (c).
    allocate(cab(nbf,nbf), source=0.0_dp)
    block
      type(int2_tdgrd_data_t), allocatable, target :: i2c
      real(kind=dp), pointer :: abc(:,:,:)
      real(kind=dp), allocatable :: avt(:,:), bvt(:,:), ev(:)
      allocate(avt(nbf,nbf), bvt(nbf,nbf), ev(lzdim_u), source=0.0_dp)
      ev = 0.0_dp; ev(1:nsocc*nvira) = 1.0_dp           ! all OV_alpha components = 1
      call umrsf_sfrogen(avt, bvt, ev, nocca, noccb)    ! avt = OV_alpha density, bvt = 0
      call orthogonal_transform('t', nbf, mo_a, avt, pa(:,:,1), wrk1)
      call orthogonal_transform('t', nbf, mo_b, bvt, pa(:,:,2), wrk1)
      allocate(i2c)
      i2c = int2_tdgrd_data_t(d2=pa, int_apb=.true., int_amb=.false., &
            tamm_dancoff=.false., scale_exchange=scale_exch)
      call int2_driver%run(i2c, cam=dft.and.infos%dft%cam_flag, &
          alpha=infos%dft%cam_alpha, beta=infos%dft%cam_beta, mu=infos%dft%cam_mu)
      abc => i2c%apb(:,:,:,1)
      call symmetrize_matrix(pa(:,:,1), nbf); call symmetrize_matrix(pa(:,:,2), nbf)
      if (dft) call utddft_fxc(basis=basis, molGrid=molGrid, isVecs=.true., wfa=mo_a, wfb=mo_b, &
          fxa=abc(:,:,1:1), fxb=abc(:,:,2:2), dxa=pa(:,:,1:1), dxb=pa(:,:,2:2), &
          nmtx=1, threshold=1.0d-15, infos=infos)
      call orthogonal_transform('n', nbf, mo_b, abc(:,:,2), cab, wrk1)  ! beta resp -> full MO
      call i2c%clean(); deallocate(i2c)
    end block

  ! ============ Phase 6 Test A: FULL gradient, UMRSF (Zu) vs RO (Zro), alpha=beta ========
  ! Per-component capture (g1e = 1e+Pulay+nuc ; gxc = +xc ; gtot = +2e-Gamma) to localize.
  ! UMRSF uses its own solved Zu (P,W); RO uses Zro. At alpha=beta the physical gradient
  ! must match; the W_ia CV-collapse residual must wash out of gtot or it is a real bug.
    natom = size(infos%atoms%zn)
    allocate(d2g(nbf,nbf,2), pg(nbf,nbf,2), spca(7,nbf,nbf), spcb_(4,nbf,nbf), &
             spcr(7,nbf,nbf), vdum(nbf,nbf), de_u(3,natom), de_r(3,natom), &
             g1e_u(3,natom), gxc_u(3,natom), gtot_u(3,natom), &
             g1e_r(3,natom), gxc_r(3,natom), gtot_r(3,natom), &
             g1e_x0(3,natom), gxc_x0(3,natom), gtot_x0(3,natom), gtot_t0(3,natom), gtot_probe(3,natom), source=0.0_dp)
  ! block densities (orbital-based; same regardless of Z)
    spca(1,:,:)=densu(1,1,:,:); spca(2,:,:)=densu(1,3,:,:); spca(3,:,:)=densu(1,5,:,:)
    spca(4,:,:)=densu(1,7,:,:); spca(5,:,:)=densu(1,9,:,:); spca(6,:,:)=densu(1,10,:,:)
    spca(7,:,:)=densu(1,11,:,:)   ! intra swap probe reverted (NO-OP: ~2e-4 effect)
    spcb_(1,:,:)=densu(1,2,:,:); spcb_(2,:,:)=densu(1,4,:,:)
    spcb_(3,:,:)=densu(1,6,:,:); spcb_(4,:,:)=densu(1,8,:,:)
    spcr(1:7,:,:)=densr(1,1:7,:,:)
  ! reserve W/P tags consumed by sf_1e_grad
    call infos%dat%reserve_data(OQP_WAO, TA_TYPE_REAL64, nbf_tri, comment=OQP_WAO_comment)
    call infos%dat%reserve_data(OQP_td_p, TA_TYPE_REAL64, nbf_tri*2, [nbf_tri,2], comment=OQP_td_p_comment)
    call tagarray_get_data(infos%dat, OQP_WAO, wao_p)
    call tagarray_get_data(infos%dat, OQP_td_p, tdp_p)

  ! ----- UMRSF full gradient (from solved Zu) -----
    block                                                    ! probe hook: zero 2e-Gamma subterm
      character(len=8) :: envv; integer :: envl, envs
      call get_environment_variable('UMRSF_NODT2', envv, envl, envs)
      umrsf_dbg_nodt2 = (envs == 0 .and. envl > 0)
      call get_environment_variable('UMRSF_ZERO', envv, envl, envs)
      umrsf_dbg_zero = 0
      if (envs == 0 .and. envl > 0) read(envv,*) umrsf_dbg_zero
    end block
    call store_wp(wa_u, wb_u, pau, pbu, 0.25_dp)             ! wfac 0.25: spin-resolved fold
                                                             ! Wa+Wb==RO wmo so RO's 0.25 applies
    call sf_1e_grad(infos, basis);  g1e_u = infos%atoms%grad
    call build_d(d2g)
    if (dft) call utddft_xc_gradient(basis=basis, molGrid=molGrid, dedft=infos%atoms%grad, &
        da=d2g(:,:,1), db=d2g(:,:,2), pa=pg(:,:,1:1), pb=pg(:,:,2:2), nmtx=1, &
        threshold=0.0_dp, infos=infos)
    gxc_u = infos%atoms%grad
    call build_d(d2g)
    call umrsf_2e_grad(basis, infos, d2g, pg, spca, spcb_, vdum);  gtot_u = infos%atoms%grad

  ! PROBE: P-relaxation-only gradient = P=T+Z (full pau/pbu) but W=W(Z=0). Then
  !   gtot_u - gtot_probe = W-relaxation (Z in the Pulay W) ;
  !   gtot_probe - gtot_t0 = P-relaxation (Z in the density). Splits the Z-relaxation error
  !   into J/Z-vector (P) vs W-builder (W).
    block
      real(kind=dp), allocatable :: zz(:), wa_z0(:,:), wb_z0(:,:)
      allocate(zz(lzdim_u), wa_z0(nbf,nbf), wb_z0(nbf,nbf), source=0.0_dp)
      call umrsf_sfrowcal(wa_z0, wb_z0, mo_energy_a, mo_energy_b, fa, fb, zz, xhua, xhub, &
              hpta(1:nocca,1:nocca), hptb(1:noccb,1:noccb), hpta, hptb, nocca, noccb)
      call store_wp(wa_z0, wb_z0, pau, pbu, 0.25_dp)          ! W=W(Z=0), P=T+Z (full)
      call sf_1e_grad(infos, basis)
      call build_d(d2g)
      if (dft) call utddft_xc_gradient(basis=basis, molGrid=molGrid, dedft=infos%atoms%grad, &
          da=d2g(:,:,1), db=d2g(:,:,2), pa=pg(:,:,1:1), pb=pg(:,:,2:2), nmtx=1, &
          threshold=0.0_dp, infos=infos)
      call build_d(d2g)
      spca(1,:,:)=densu(1,1,:,:); spca(2,:,:)=densu(1,3,:,:); spca(3,:,:)=densu(1,5,:,:)
      spca(4,:,:)=densu(1,7,:,:); spca(5,:,:)=densu(1,9,:,:); spca(6,:,:)=densu(1,10,:,:)
      spca(7,:,:)=densu(1,11,:,:)
      spcb_(1,:,:)=densu(1,2,:,:); spcb_(2,:,:)=densu(1,4,:,:)
      spcb_(3,:,:)=densu(1,6,:,:); spcb_(4,:,:)=densu(1,8,:,:)
      call umrsf_2e_grad(basis, infos, d2g, pg, spca, spcb_, vdum);  gtot_probe = infos%atoms%grad
    end block

  ! ----- RO full gradient (from solved Zro) -----  (vdum=0 -> beta W transform vanishes)
    call store_wp(wmo_ro, vdum, paro, pbro, 0.25_dp)
    call sf_1e_grad(infos, basis);  g1e_r = infos%atoms%grad
    call build_d(d2g)
    if (dft) call utddft_xc_gradient(basis=basis, molGrid=molGrid, dedft=infos%atoms%grad, &
        da=d2g(:,:,1), db=d2g(:,:,2), pa=pg(:,:,1:1), pb=pg(:,:,2:2), nmtx=1, &
        threshold=0.0_dp, infos=infos)
    gxc_r = infos%atoms%grad
    call build_d(d2g)
    call mrsf_2e_grad(basis, infos, d2g, pg, spcr, vdum);  gtot_r = infos%atoms%grad
    de_u = gtot_u; de_r = gtot_r

  ! ----- BUG 1 probe: X->0 limit. Zero P, W, and the X-dependent block densities (~X^2);
  !       keep ground D. The total MUST reduce to the reference SCF gradient (Omega^xi->0).
    call store_wp(vdum, vdum, vdum, vdum, 0.5_dp)   ! W=0, P=0 -> pg=0
    call sf_1e_grad(infos, basis);  g1e_x0 = infos%atoms%grad
    call build_d(d2g)
    if (dft) call utddft_xc_gradient(basis=basis, molGrid=molGrid, dedft=infos%atoms%grad, &
        da=d2g(:,:,1), db=d2g(:,:,2), pa=pg(:,:,1:1), pb=pg(:,:,2:2), nmtx=1, &
        threshold=0.0_dp, infos=infos)
    gxc_x0 = infos%atoms%grad
    call build_d(d2g)
    spca = 0.0_dp; spcb_ = 0.0_dp            ! ball + inter + intra all ~X^2 -> 0
    call umrsf_2e_grad(basis, infos, d2g, pg, spca, spcb_, vdum);  gtot_x0 = infos%atoms%grad

  ! ----- T-only (unrelaxed) gradient: P=T (Z=0), W=W(Z=0). The unrelaxed gradient keeps the
  !       Q-part W (H[X,X]+Fock.T+H+[T]); W=0 was WRONG (apples-to-oranges). Splits the
  !       excitation gradient into the unrelaxed (T,W(Z=0)) part vs the relaxation Z part.
    block
      real(kind=dp), allocatable :: zt(:), pat(:,:), pbt(:,:), wa_t0(:,:), wb_t0(:,:)
      allocate(zt(lzdim_u), pat(nbf,nbf), pbt(nbf,nbf), wa_t0(nbf,nbf), wb_t0(nbf,nbf), source=0.0_dp)
      call umrsf_sfropcal(pat, pbt, tij, tab, zt, nocca, noccb)   ! P = T (Z=0)
    ! W(Z=0): umrsf_sfrowcal with z=0 -> Z-terms/couplings vanish, keeps xhx (2H[X,X]+Fock.T)
    ! and H+[T] occ-occ (hpta/hptb = H+[T] full MO). hppmoa/hppmob unused (occ-virt reverted).
      call umrsf_sfrowcal(wa_t0, wb_t0, mo_energy_a, mo_energy_b, fa, fb, zt, xhua, xhub, &
              hpta(1:nocca,1:nocca), hptb(1:noccb,1:noccb), hpta, hptb, nocca, noccb)
      call store_wp(wa_t0, wb_t0, pat, pbt, 0.25_dp)             ! W(Z=0), pg = T (AO), wfac 0.25
      call sf_1e_grad(infos, basis)
      call build_d(d2g)
      if (dft) call utddft_xc_gradient(basis=basis, molGrid=molGrid, dedft=infos%atoms%grad, &
          da=d2g(:,:,1), db=d2g(:,:,2), pa=pg(:,:,1:1), pb=pg(:,:,2:2), nmtx=1, &
          threshold=0.0_dp, infos=infos)
      call build_d(d2g)
      spca(1,:,:)=densu(1,1,:,:); spca(2,:,:)=densu(1,3,:,:); spca(3,:,:)=densu(1,5,:,:)
      spca(4,:,:)=densu(1,7,:,:); spca(5,:,:)=densu(1,9,:,:); spca(6,:,:)=densu(1,10,:,:)
      spca(7,:,:)=densu(1,11,:,:)
      spcb_(1,:,:)=densu(1,2,:,:); spcb_(2,:,:)=densu(1,4,:,:)
      spcb_(3,:,:)=densu(1,6,:,:); spcb_(4,:,:)=densu(1,8,:,:)
      call umrsf_2e_grad(basis, infos, d2g, pg, spca, spcb_, vdum);  gtot_t0 = infos%atoms%grad
    end block

  ! ============ Dump ============
    open(newunit=iu, file='umrsf_qr_dump.txt', status='replace', action='write')
    write(iu,'(a,5i6)') '# nbf nocca noccb lzdim_u lzdim_ro', nbf, nocca, noccb, lzdim_u, lzdim_ro
    call dump_mat(iu, 'QA', qa, nbf)
    call dump_mat(iu, 'QB', qb, nbf)
    call dump_mat(iu, 'HXA_U', hxa_u, nbf)    ! UMRSF H[X,X] alpha
    call dump_mat(iu, 'HXA_R', hxa_r, nbf)    ! RO golden H[X,X] alpha (clean, pre-sfrorhs)
    call dump_mat(iu, 'HXB_U', hxb_u, nbf)    ! UMRSF H[X,X] beta
    call dump_mat(iu, 'HXB_R', hxb_r, nbf)    ! RO golden H[X,X] beta
    call dump_mat(iu, 'HPTA', hpta, nbf)
    call dump_mat(iu, 'HPTB', hptb, nbf)
    call dump_mat(iu, 'FA', fa, nbf)
    call dump_mat(iu, 'FB', fb, nbf)
    call dump_vec(iu, 'EA', mo_energy_a, nbf)
    call dump_rect(iu, 'TIJ', tij, nocca, nocca)
    call dump_rect(iu, 'TAB', tab, nvirb, nvirb)
    call dump_rect(iu, 'AB1A', ab1_mo_a, nocca, nvira)
    call dump_rect(iu, 'AB1B', ab1_mo_b, noccb, nvirb)
    call dump_vec(iu, 'RHS_U',  rhs_u,  lzdim_u)
    call dump_vec(iu, 'RHS_RO', rhs_ro, lzdim_ro)
    call dump_rect(iu, 'JU',  ju,  lzdim_u,  lzdim_u)
    call dump_rect(iu, 'JRO', jro, lzdim_ro, lzdim_ro)
    call dump_vec(iu, 'ZU',  zu,  lzdim_u)
    call dump_vec(iu, 'ZRO', zro, lzdim_ro)
    call dump_mat(iu, 'PAU',  pau,  nbf)   ! real relaxed P alpha (from solved Zu)
    call dump_mat(iu, 'PBU',  pbu,  nbf)
    call dump_mat(iu, 'PAUC', pauc, nbf)   ! controlled-Z P alpha (for RO compare)
    call dump_mat(iu, 'PBUC', pbuc, nbf)
    call dump_mat(iu, 'PARO', paro, nbf)   ! RO golden P alpha
    call dump_mat(iu, 'PBRO', pbro, nbf)
    call dump_mat(iu, 'WA_UC', wa_uc, nbf) ! umrsf W alpha (controlled Z)
    call dump_mat(iu, 'WB_UC', wb_uc, nbf) ! umrsf W beta  (controlled Z)
    call dump_mat(iu, 'WMO_RO', wmo_ro, nbf) ! RO golden W (spin-summed)
    call dump_mat(iu, 'WA_U', wa_u, nbf)   ! real W alpha (from solved Zu) -> Phase 5
    call dump_mat(iu, 'WB_U', wb_u, nbf)   ! real W beta
    call dump_mat(iu, 'CAB', cab, nbf)     ! R1: beta cross-response (full MO) to OV_alpha
  ! GMRES-vs-dense: [max|z_gmres - z_dense|, gmres residual, gmres iters]
    call dump_vec(iu, 'GMRES', [gmres_zdiff, gmres_resid, real(gmres_iters,dp)], 3)
    call dump_rect(iu, 'G1E_U', g1e_u, 3, natom)   ! UMRSF: 1e+Pulay+nuc
    call dump_rect(iu, 'G1E_R', g1e_r, 3, natom)
    call dump_rect(iu, 'GXC_U', gxc_u, 3, natom)   ! +xc
    call dump_rect(iu, 'GXC_R', gxc_r, 3, natom)
    call dump_rect(iu, 'GTOT_U', gtot_u, 3, natom) ! +2e-Gamma (full gradient)
    call dump_rect(iu, 'GTOT_R', gtot_r, 3, natom)
    call dump_rect(iu, 'G1E_X0', g1e_x0, 3, natom)
    call dump_rect(iu, 'GXC_X0', gxc_x0, 3, natom)
    call dump_rect(iu, 'GTOT_X0', gtot_x0, 3, natom)
    call dump_rect(iu, 'GTOT_T0', gtot_t0, 3, natom)
    call dump_rect(iu, 'GTOT_PROBE', gtot_probe, 3, natom)
  ! Phase 5 inter-block fold: spin-resolved inter densities (umrsfcbc) vs RO (mrsfcbc)
  ! alpha inter: bco1a=densu5 bco2a=densu7 bo1va=densu3 bo2va=densu1
    call dump_mat(iu, 'BCO1A', densu(1,5,:,:), nbf); call dump_mat(iu, 'BCO2A', densu(1,7,:,:), nbf)
    call dump_mat(iu, 'BO1VA', densu(1,3,:,:), nbf); call dump_mat(iu, 'BO2VA', densu(1,1,:,:), nbf)
  ! beta inter: bco1b=densu6 bco2b=densu8 bo1vb=densu4 bo2vb=densu2
    call dump_mat(iu, 'BCO1B', densu(1,6,:,:), nbf); call dump_mat(iu, 'BCO2B', densu(1,8,:,:), nbf)
    call dump_mat(iu, 'BO1VB', densu(1,4,:,:), nbf); call dump_mat(iu, 'BO2VB', densu(1,2,:,:), nbf)
  ! RO inter: bco1=densr3 bco2=densr4 bo1v=densr2 bo2v=densr1
    call dump_mat(iu, 'BCO1R', densr(1,3,:,:), nbf); call dump_mat(iu, 'BCO2R', densr(1,4,:,:), nbf)
    call dump_mat(iu, 'BO1VR', densr(1,2,:,:), nbf); call dump_mat(iu, 'BO2VR', densr(1,1,:,:), nbf)
    close(iu)

    call int2_driver%clean()
    if (dft) call dftclean(infos)
    write(iw,'(1x,a)') 'UMRSF Q/R test: dumped umrsf_qr_dump.txt'
    close(iw)

  contains

    subroutine ident_precond(x_in, x_out)   ! identity preconditioner for the gmres test
      real(kind=dp), intent(in)  :: x_in(:)
      real(kind=dp), intent(out) :: x_out(:)
      x_out = x_in
    end subroutine ident_precond

    ! Store spin-resolved W (AO) -> wao_p and P (AO) -> tdp_p/pg, matching the RO
    ! z-vector convention (wao*0.25, td_p*0.5). RO: pass wb=0.
    subroutine store_wp(wa, wb, pamo, pbmo, wfac)
      use mathlib, only: pack_matrix, symmetrize_matrix, orthogonal_transform
      real(kind=dp), intent(in) :: wa(:,:), wb(:,:), pamo(:,:), pbmo(:,:), wfac
      call orthogonal_transform('t', nbf, mo_a, wa, wrk1, wrk3)
      call orthogonal_transform('t', nbf, mo_b, wb, wrk2, wrk3)
      wrk1 = wrk1 + wrk2
      call symmetrize_matrix(wrk1, nbf)
      call pack_matrix(wrk1, wao_p);  wao_p = wao_p*wfac
      call orthogonal_transform('t', nbf, mo_a, pamo, pg(:,:,1), wrk3)
      call orthogonal_transform('t', nbf, mo_b, pbmo, pg(:,:,2), wrk3)
    ! Mirror the golden RO density convention (z_vector line 1268-1293): AO P =
    ! 0.5*symmetrize(C P_MO C^T). symmetrize = M+M^T, so symmetric (T) blocks are
    ! preserved while the upper-only occ-virt (Z) block is properly symmetrized; the
    ! 0.5 must hit pg too, since the 2e Gamma consumes pg (golden unpacks td_p).
      call symmetrize_matrix(pg(:,:,1), nbf)
      call symmetrize_matrix(pg(:,:,2), nbf)
      pg = pg*0.5_dp
      call pack_matrix(pg(:,:,1), tdp_p(:,1));  call pack_matrix(pg(:,:,2), tdp_p(:,2))
    end subroutine store_wp

    subroutine build_d(d)   ! ground-state AO density D_sigma = C_occ C_occ^T
      real(kind=dp), intent(out) :: d(:,:,:)
      call dgemm('n','t', nbf, nbf, nocca, 1.0_dp, mo_a, nbf, mo_a, nbf, 0.0_dp, d(:,:,1), nbf)
      call dgemm('n','t', nbf, nbf, noccb, 1.0_dp, mo_b, nbf, mo_b, nbf, 0.0_dp, d(:,:,2), nbf)
    end subroutine build_d

    subroutine uks_fock_mo(focktri, mo, fmo, scr, n)
      use mathlib, only: orthogonal_transform_sym, unpack_matrix
      real(kind=dp), intent(in)  :: focktri(:), mo(:,:)
      real(kind=dp), intent(out) :: fmo(:,:)
      real(kind=dp), intent(inout) :: scr(:,:)
      integer, intent(in) :: n
      real(kind=dp), allocatable :: ftri(:)
      allocate(ftri(n*(n+1)/2))
      call orthogonal_transform_sym(n, n, focktri, mo, n, ftri)
      call unpack_matrix(ftri, fmo)
      deallocate(ftri)
    end subroutine uks_fock_mo

    subroutine dump_mat(u, name, a, n)
      integer, intent(in) :: u, n
      character(*), intent(in) :: name
      real(kind=dp), intent(in) :: a(:,:)
      integer :: i, j
      write(u,'(a,1x,a,2i6)') '@MAT', name, n, n
      do j = 1, n
        do i = 1, n
          write(u,'(es24.16)') a(i,j)
        end do
      end do
    end subroutine dump_mat

    subroutine dump_rect(u, name, a, n1, n2)
      integer, intent(in) :: u, n1, n2
      character(*), intent(in) :: name
      real(kind=dp), intent(in) :: a(:,:)
      integer :: i, j
      write(u,'(a,1x,a,2i6)') '@MAT', name, n1, n2
      do j = 1, n2
        do i = 1, n1
          write(u,'(es24.16)') a(i,j)
        end do
      end do
    end subroutine dump_rect

    subroutine dump_vec(u, name, v, m)
      integer, intent(in) :: u, m
      character(*), intent(in) :: name
      real(kind=dp), intent(in) :: v(:)
      integer :: i
      write(u,'(a,1x,a,i8)') '@VEC', name, m
      do i = 1, m
        write(u,'(es24.16)') v(i)
      end do
    end subroutine dump_vec

  end subroutine tdhf_umrsf_qrtest

end module tdhf_mrsf_z_vector_mod