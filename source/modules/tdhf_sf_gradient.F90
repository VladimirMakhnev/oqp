module tdhf_sf_gradient_mod

  use precision, only: dp
  use grd2, only: grd2_driver, grd2_compute_data_t
  use basis_tools, only: basis_set
  use types, only: information

  implicit none

  character(len=*), parameter :: module_name = "tdhf_sf_gradient_mod"

  public tdhf_sf_gradient

  type, extends(grd2_compute_data_t) :: grd2_sf_compute_data_t
    real(kind=dp), pointer :: d2(:,:,:) => null()
    real(kind=dp), pointer :: p2(:,:,:) => null()
    real(kind=dp), pointer :: v2(:,:) => null()
    integer :: nbf = 0
  contains
    procedure :: init => grd2_sf_compute_data_t_init
    procedure :: clean => grd2_sf_compute_data_t_clean
    procedure :: get_density => grd2_sf_compute_data_t_get_density
  end type

contains

  subroutine sf_gradient_C(c_handle) bind(C, name="tdhf_sf_gradient")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_sf_gradient(inf)
  end subroutine sf_gradient_c


!> @brief Compute SF-TDDFT analytical energy gradient
!>
!> The excited-state gradient consists of three contributions:
!>
!>   dE/dR = dE_1e/dR + dE_2e/dR + dE_xc/dR
!>
!> where:
!>   dE_1e/dR = Tr[(D+P) * dH/dR] - Tr[W * dS/dR] + nuclear terms
!>   dE_2e/dR = sum_ijkl Gamma_ijkl * d(ij|kl)/dR
!>   dE_xc/dR = integral f_xc[P] * drho/dR (DFT only)
!>
!> Input (from Z-vector calculation):
!>   - D: ground-state density matrix
!>   - P: relaxed difference density (T + Z)
!>   - W: energy-weighted density matrix
!>   - V: (A-B)*X response vector
!>
!> Reference: Furche & Ahlrichs, JCP 117, 7433 (2002)
!>
  subroutine tdhf_sf_gradient(infos)
    use io_constants, only: iw
    use oqp_tagarray_driver

    use types, only: information
    use strings, only: Cstring, fstring
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort

    use grd1, only: eijden, print_gradient
    use util, only: measure_time
    use tdhf_lib, only: &
      iatogen, mntoia
    use dft, only: dft_initialize, dftclean
    use mathlib, only: symmetrize_matrix
    use mod_dft_molgrid, only: dft_grid_t
    use mod_dft_gridint_tdxc_grad, only: utddft_xc_gradient
    use mathlib, only: unpack_matrix
    use printing, only: print_module_info

    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_sf_gradient"

    type(basis_set), pointer :: basis
    type(information), target, intent(inout) :: infos

    integer :: nbf, nbf_tri

    type(dft_grid_t) :: molGrid

  ! General data
    logical :: dft
    integer :: scf_type

    real(kind=dp), allocatable :: p(:,:,:), v(:,:,:), d(:,:,:)

    ! tagarray
    real(kind=dp), contiguous, pointer :: &
      dmat_a(:), dmat_b(:), td_abxc(:,:), td_p(:,:)
    character(len=*), parameter :: tags_general(*) = (/ character(len=80) :: &
      OQP_DM_A, OQP_DM_B, OQP_td_abxc, OQP_td_p /)

    scf_type = infos%control%scftype
    dft = infos%control%hamilton == 20

  ! Files open
    open (unit=IW, file=infos%log_filename, position="append")
  !
    call print_module_info('SF_Grad','Computing Gradient of SF-TDDFT')
!
    write(iw,'(/5X,"Gradient options"/&
                &5X,18("-")/&
                &5X,"Target State: ",I8/&
                &,5X,"*Note that the ground state of SF-TDDFT is 1.*"/)')&
                & infos%tddft%target_state

  ! Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
    call tagarray_get_data(infos%dat, OQP_td_abxc, td_abxc)
    call tagarray_get_data(infos%dat, OQP_td_p, td_p)

  ! Allocate H, S ,T and D matrices
    nbf = basis%nbf
    nbf_tri = nbf*(nbf+1)/2

! =============================================================================
!   STEP 1: One-electron gradient
!   Includes: nuclear repulsion, overlap (with W), kinetic, nuclear attraction
! =============================================================================
    call flush(iw)
    call sf_1e_grad(infos, basis)

    write(iw,"(' ..... End Of 1-Eelectron Gradient ......')")
    call measure_time(print_total=1, log_unit=iw)
    call flush(iw)

! =============================================================================
!   STEP 2: Prepare density matrices for 2e and XC gradients
!   d = ground-state density (alpha, beta)
!   p = relaxed difference density P = T + Z (from Z-vector)
! =============================================================================
    allocate(d(nbf,nbf,2), source=0.0d0)
    allocate(p(nbf,nbf,2), source=0.0d0)

    call unpack_matrix(td_p(:,1), p(:,:,1))
    call unpack_matrix(td_p(:,2), p(:,:,2))

    call unpack_matrix(dmat_a, d(:,:,1))
    call unpack_matrix(dmat_b, d(:,:,2))

! =============================================================================
!   STEP 3: DFT exchange-correlation gradient (if DFT)
!   dE_xc/dR = integral f_xc * dP/dR over grid
! =============================================================================
    if (dft) then
      call dft_initialize(infos, basis, molGrid, verbose=.true.)

      call utddft_xc_gradient(basis=basis, &
           molGrid=molGrid, &
           dedft=infos%atoms%grad, &
           da=d(:,:,1), &
           db=d(:,:,2), &
           pa=p(:,:,1:1), &
           pb=p(:,:,2:2), &
           nmtx=1, &
           threshold=0.0d0, &
           infos=infos)

      call dftclean(infos)
      call measure_time(print_total=1, log_unit=iw)
      call flush(iw)
    end if

! =============================================================================
!   STEP 4: Two-electron gradient
!   dE_2e/dR = sum_ijkl Gamma_ijkl * d(ij|kl)/dR
!   where Gamma is the two-particle density matrix
! =============================================================================
    allocate(v(nbf,nbf,2), source=0.0d0)
    v(:,:,1) = td_abxc  ! Response vector (A-B)*X
    call sf_2e_grad(basis, infos, d, p, v(:,:,1))

    call print_gradient(infos)

!   Print timings
    call measure_time(print_total=1, log_unit=iw)

    close(iw)

  end subroutine tdhf_sf_gradient

!###############################################################################

!> @brief Compute one-electron contribution to SF-TDDFT gradient
!>
!> One-electron gradient terms:
!>   dE_1e/dR = Tr[(D+P) * dH/dR] - Tr[(L+W) * dS/dR] + V_nn
!>
!> where:
!>   D = ground-state density
!>   P = relaxed difference density (T + Z)
!>   H = core Hamiltonian (kinetic + nuclear attraction)
!>   L = Lagrangian from SCF (orbital energy weighted density)
!>   W = excited-state energy-weighted density
!>   S = overlap matrix
!>   V_nn = nuclear repulsion gradient
!>
  subroutine sf_1e_grad(infos, basis)

    use oqp_tagarray_driver
    use types, only: information
    use basis_tools, only: basis_set
    use util, only: measure_time
    use messages, only: show_message, WITH_ABORT
    use precision, only: dp
    use constants, only: tol_int
    use grd1, only: eijden, print_gradient, &
            grad_nn, grad_ee_overlap, &
            grad_ee_kinetic, grad_en_hellman_feynman, grad_en_pulay, grad_1e_ecp

    use mathlib, only: symmetrize_matrix

    implicit none

    character(len=*), parameter :: subroutine_name = "sf_1e_grad"

    type(information), intent(inout) :: infos
    type(basis_set), intent(inout) :: basis

    real(kind=dp), allocatable :: dens(:)
    real(kind=dp) :: tol
    integer :: nbf, nbf_tri, ok

    ! tagarray
    real(kind=dp), pointer :: dmat_a(:), dmat_b(:), wao(:), td_p(:,:)
    character(len=*), parameter :: tags_general(4) = (/ character(len=80) :: &
      OQP_DM_A, OQP_DM_B, OQP_WAO, OQP_td_p /)

    nbf = basis%nbf
    nbf_tri = nbf*(nbf+1)/2

    tol = tol_int*log(10.0_dp)

!   initial memory allocation
    allocate(dens(nbf_tri), source=0.0_dp, stat=ok)
    if(ok/=0) call show_message('Cannot allocate memory', WITH_ABORT)

    call data_has_tags(infos%dat, tags_general, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)
    call tagarray_get_data(infos%dat, OQP_WAO, wao)
    call tagarray_get_data(infos%dat, OQP_td_p, td_p)

    associate( grad   => infos%atoms%grad &
             , xyz    => infos%atoms%xyz &
             , zn     => infos%atoms%zn &
             , p      => td_p &
             , w      => wao &
             , urohf  => infos%control%scftype>=2 &
        )

      grad = 0.0d0

      ! Nuclear repulsion gradient: dV_nn/dR
      call grad_nn(infos%atoms, infos%basis%ecp_zn_num)

      ! Build Lagrangian L from SCF orbital energies
      ! L_ij = sum_k n_k * epsilon_k * C_ik * C_jk
      call eijden(dens, nbf, infos)

      ! Add excited-state W matrix: total Lagrangian = L + W
      ! This gives -Tr[(L+W) * dS/dR] term in gradient
      dens = dens + w

      ! Overlap gradient: -Tr[(L+W) * dS/dR]
      call grad_ee_overlap(basis, dens, grad, logtol=tol)

      ! Build total density for core Hamiltonian gradient
      ! D_total = D_ground + P_excited (alpha + beta)
      dens = dmat_a + p(:,1)
      if (infos%control%scftype>=2) then
        dens = dens + dmat_b + p(:,2)
      end if

      ! Hellmann-Feynman: Tr[D_total * dV_ne/dR] (point charges)
      call grad_en_hellman_feynman(basis, xyz, zn, dens, grad, logtol=tol)

      ! Kinetic energy gradient: Tr[D_total * dT/dR]
      call grad_ee_kinetic(basis, dens, grad, logtol=tol)

      ! Pulay force: Tr[D_total * dV_ne/dR] (basis function derivatives)
      call grad_en_pulay(basis, xyz, zn, dens, grad, logtol=tol)

      ! ECP gradient (if present)
      call grad_1e_ecp(infos, basis, xyz, dens, grad, logtol=tol)

    end associate

   end subroutine

!###############################################################################

!> @brief Driver for two-electron contribution to SF-TDDFT gradient
!>
!> Computes: dE_2e/dR = sum_ijkl Gamma_ijkl * d(ij|kl)/dR
!>
!> The two-particle density matrix Gamma contains:
!>   - Coulomb terms: (D+P)_ij * D_kl + D_ij * P_kl
!>   - Exchange terms: (D+P)_ik * D_jl + D_ik * P_jl (scaled by HF exchange)
!>   - Response terms: V_ik * V_jl (from (A-B)*X, scaled by response HF exchange)
!>
!> @param[in] d  Ground-state density matrices (alpha, beta)
!> @param[in] p  Relaxed difference density P = T + Z (alpha, beta)
!> @param[in] v  Response vector (A-B)*X in AO basis
!>
  subroutine sf_2e_grad(basis, infos, d, p, v)

    use basis_tools, only: basis_set
    use precision, only: dp
    use messages, only: show_message, WITH_ABORT
    use types, only: information
    use io_constants, only: iw

    implicit none

    type(information), target, intent(inout) :: infos
    type(basis_set) :: basis
    real(kind=dp), contiguous, target :: p(:,:,:), d(:,:,:), v(:,:)

    logical :: urohf, dft
    real(kind=dp) :: scale_exch  !> HF scale in Reference
    real(kind=dp) :: scale_exch2 !> HF scale in Response

    integer :: ok
    real(kind=dp), allocatable :: de(:,:)
    class(grd2_compute_data_t), allocatable :: gcomp

    dft = infos%control%hamilton == 20 ! dft or hf
    urohf = infos%control%scftype >= 2

    scale_exch = 1.0_dp
    scale_exch2 = 1.0_dp
    if (dft) then
      scale_exch = infos%dft%HFscale
      scale_exch2 = infos%tddft%HFscale
    end if

    allocate(de(3,ubound(infos%atoms%zn,1)), &
            source=0.0d0, &
            stat=ok)

    if(ok/=0) call show_message('cannot allocate memory', WITH_ABORT)

    write(iw, '(/7x,"Fitting parameters")')
    if (.not.infos%dft%cam_flag) then
      write(iw, '(10x,"Exact HF exchange:")')
      write(iw, '(5x,"Reference: |", t20, f6.3, t29, "|")') scale_exch
      write(iw, '(5x,"Response:  |", t20, f6.3, t29, "|")') scale_exch2
    else
      write(iw, '(10x,"CAM parametres:")')
      write(iw, '(16x,"|   alpha   |    beta   |     mu    |")')
      write(iw, '(5x,"Reference: |", t20, f6.3, t29, "|", t32, f6.3, t41, "|", t44, f6.3, t53, "|")') &
         infos%dft%cam_alpha, infos%dft%cam_beta, infos%dft%cam_mu
      write(iw, '(5x,"Response:  |", t20, f6.3, t29, "|", t32, f6.3, t41, "|", t44, f6.3, t53, "|")') &
         infos%tddft%cam_alpha, infos%tddft%cam_beta, infos%tddft%cam_mu
    end if
    write(iw, '(10x,"Spin-pair coupling parametres:")')
    write(iw, '(16x,"|   CO-CO   |   OV-OV   |   CO-OV   |")')
    write(iw, '(16x,"|", t20, f6.3, t29, "|", t32, f6.3, t41, "|", t44, f6.3, t53, "|")') &
       infos%tddft%HFscale, infos%tddft%HFscale, infos%tddft%HFscale

    gcomp =  grd2_sf_compute_data_t( d2 = d &
                                   , p2 = p &
                                   , v2 = v &
                                   , nbf = basis%nbf )

    call gcomp%init()

    call grd2_driver(infos, basis, de, gcomp, &
                     cam = dft.and.infos%dft%cam_flag, &
                     alpha = infos%tddft%cam_alpha, &
                     beta = infos%tddft%cam_beta, &
                     mu = infos%tddft%cam_mu)

    infos%atoms%grad = infos%atoms%grad + de

    call gcomp%clean()

  end subroutine

!###############################################################################

!> @brief Initialize density matrices for 2e gradient calculation
!>
!> Transforms density matrices from (alpha, beta) to (total, spin) representation:
!>   D_total = D_alpha + D_beta  (Coulomb contribution)
!>   D_spin  = D_alpha - D_beta  (Exchange contribution)
!>
!> This representation simplifies the 2e gradient formulas:
!>   Coulomb: uses D_total only
!>   Exchange: uses both D_total and D_spin for proper spin coupling
!>
  subroutine grd2_sf_compute_data_t_init(this)
    implicit none
    class(grd2_sf_compute_data_t), target, intent(inout) :: this

    call this%clean()

    ! Convert ground-state density: (alpha, beta) -> (total, spin)
    ! d2(:,:,1) = D_alpha + D_beta (total density for Coulomb)
    ! d2(:,:,2) = D_alpha - D_beta (spin density for exchange)
    this%d2(:,:,1) = this%d2(:,:,1) +   this%d2(:,:,2)
    this%d2(:,:,2) = this%d2(:,:,1) - 2*this%d2(:,:,2)

    ! Same transformation for relaxed difference density P
    this%p2(:,:,1) = this%p2(:,:,1) +   this%p2(:,:,2)
    this%p2(:,:,2) = this%p2(:,:,1) - 2*this%p2(:,:,2)

  end subroutine

!###############################################################################

  subroutine grd2_sf_compute_data_t_clean(this)
    implicit none
    class(grd2_sf_compute_data_t), target, intent(inout) :: this
  end subroutine

!###############################################################################

!> @brief Compute two-particle density matrix elements for SF-TDDFT 2e gradient
!>
!> Forms the effective density Gamma_ijkl for 2e gradient:
!>   dE_2e/dR = sum_ijkl Gamma_ijkl * d(ij|kl)/dR
!>
!> Gamma consists of three parts:
!>
!> 1. Coulomb contribution (df1_coul):
!>    Gamma^J_ijkl = 4 * [(D+P)_ij * D_kl + D_ij * P_kl]
!>
!> 2. Exchange contribution (dq1, scaled by HF exchange):
!>    Gamma^K_ijkl = -hfscale * [(D+P)_ik * D_jl + D_ik * P_jl + (ik<->il)]
!>    Uses both total and spin densities for proper open-shell treatment
!>
!> 3. Response contribution (dt2, scaled by response HF exchange):
!>    Gamma^V_ijkl = -2 * hfscale2 * [V_ik * V_jl + V_il * V_jk + transposes]
!>    where V = (A-B)*X is the response vector
!>
!> @param[in]  id     Shell quartet indices (i,j,k,l)
!> @param[out] dab    Density matrix product for this shell quartet
!> @param[out] dabmax Maximum absolute value (for screening)
!>
  subroutine grd2_sf_compute_data_t_get_density(this, basis, id, dab, dabmax)

    implicit none

    class(grd2_sf_compute_data_t), target, intent(inout) :: this
    type(basis_set), intent(in) :: basis
    integer, intent(in) :: id(4)
    real(kind=dp), target, intent(out) :: dab(*)
    real(kind=dp), intent(out) :: dabmax

    real(kind=dp) :: df1, dq1, dt2, df1_coul
    real(kind=dp) :: coulfact, xcfact, xcfact2
    integer :: i, j, k, l
    integer :: loc(4)
    integer :: nbf(4)
    real(kind=dp), pointer :: ab(:,:,:,:)
    integer :: i1, j1, k1, l1

    coulfact = 4*this%coulscale    ! Coulomb prefactor
    xcfact = this%hfscale          ! Reference HF exchange
    xcfact2 = this%hfscale2        ! Response HF exchange
    dabmax = 0
    loc = basis%ao_offset(id)-1

    nbf = basis%naos(id)

    ab(1:nbf(4),1:nbf(3),1:nbf(2),1:nbf(1)) => dab(1:product(nbf))

    do i = 1, nbf(1)
      i1 = loc(1) + i

      do j = 1, nbf(2)
        j1 = loc(2) + j

        do k = 1, nbf(3)
          k1 = loc(3) + k

          do l = 1, nbf(4)
            l1 = loc(4) + l

            ! ----- Coulomb contribution -----
            ! Gamma^J = (D+P)_ij * D_kl + D_ij * P_kl
            ! Uses total density only (index 1)
            df1 = (this%d2(i1,j1,1)+this%p2(i1,j1,1))*this%d2(k1,l1,1) &
                +  this%d2(i1,j1,1)                  *this%p2(k1,l1,1)
            df1_coul = df1 * coulfact

            if (xcfact /= 0.0_dp .or. xcfact2 /= 0.0_dp) then
              ! ----- Exchange contribution -----
              ! Gamma^K = (D+P)_ik * D_jl + D_ik * P_jl + (ik<->il)
              ! Uses both total (1) and spin (2) densities
              dq1 = (this%d2(i1,k1,1)+this%p2(i1,k1,1))*this%d2(j1,l1,1) &
                  +  this%d2(i1,k1,1)                  *this%p2(j1,l1,1) &
                  + (this%d2(i1,l1,1)+this%p2(i1,l1,1))*this%d2(j1,k1,1) &
                  +  this%d2(i1,l1,1)                  *this%p2(j1,k1,1) &
                  + (this%d2(i1,k1,2)+this%p2(i1,k1,2))*this%d2(j1,l1,2) &
                  +  this%d2(i1,k1,2)                  *this%p2(j1,l1,2) &
                  + (this%d2(i1,l1,2)+this%p2(i1,l1,2))*this%d2(j1,k1,2) &
                  +  this%d2(i1,l1,2)                  *this%p2(j1,k1,2)

              ! ----- Response contribution -----
              ! Gamma^V = V_ik * V_jl + V_il * V_jk + transposes
              ! From (A-B)*X response vector
              dt2 = this%v2(i1,k1)*this%v2(j1,l1) &
                  + this%v2(k1,i1)*this%v2(l1,j1) &
                  + this%v2(i1,l1)*this%v2(j1,k1) &
                  + this%v2(l1,i1)*this%v2(k1,j1)

              ! Total: Gamma = Gamma^J - hfscale*Gamma^K - 2*hfscale2*Gamma^V
              df1 = df1_coul - xcfact*dq1 - xcfact2*2.0_dp*dt2
            else
              df1 = df1_coul
            end if

            dabmax = max(dabmax, abs(df1))
            ab(l,k,j,i) = df1*product(basis%bfnrm([i1,j1,k1,l1]))
          end do
        end do
      end do
    end do

  end subroutine grd2_sf_compute_data_t_get_density

!###############################################################################

end module tdhf_sf_gradient_mod
