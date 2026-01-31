module tdhf_sf_z_vector_mod

  implicit none

  character(len=*), parameter :: module_name = "tdhf_sf_z_vector_mod"

contains

  subroutine tdhf_sf_z_vector_C(c_handle) bind(C, name="tdhf_sf_z_vector")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call tdhf_sf_z_vector(inf)
  end subroutine tdhf_sf_z_vector_C


!> @brief Solve Z-vector equation for SF-TDDFT analytical gradients
!>
!> This subroutine solves the coupled-perturbed equation:
!>
!>   (A + B) * Z = -RHS
!>
!> for the orbital relaxation contribution to excited-state gradients.
!>
!> Algorithm (preconditioned conjugate gradient):
!>   1. Build RHS from H[T] and H[X]*X terms
!>   2. Initialize preconditioner M = diag(epsilon_a - epsilon_i)
!>   3. Iterate: pk -> (A+B)*pk -> update Z
!>   4. Converge when ||residual|| < tolerance
!>
!> After Z-vector converges:
!>   5. Build relaxed density P = T + Z
!>   6. Compute H[P] response
!>   7. Construct W matrix (energy-weighted density)
!>
!> Output:
!>   - OQP_td_p: relaxed density matrices (alpha, beta)
!>   - OQP_WAO: W matrix in AO basis for gradient
!>
!> Supports both UHF (separate alpha/beta) and ROHF (doc-socc-virt) references.
!>
!> Reference: Furche & Ahlrichs, JCP 117, 7433 (2002)
!>            Shao, Head-Gordon, Krylov, JCP 118, 4807 (2003)
!>
  subroutine tdhf_sf_z_vector(infos)

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
    use tdhf_lib, only: iatogen, mntoia
    use tdhf_sf_lib, only: sfrorhs, &
      sfromcal, sfrogen, sfrolhs, pcgrbpini, &
      pcgb, sfropcal, sfrowcal, sfdmat, &
      ! UHF-specific functions
      sfrcalc, xecalc, sfuesum, sfgen, sfpcal, sflhs, sfwcal
    use dft, only: dft_initialize, dftclean
    use mod_dft_gridint_fxc, only: utddft_fxc
    use mathlib, only: symmetrize_matrix, orthogonal_transform_sym, orthogonal_transform
    use mod_dft_molgrid, only: dft_grid_t
    use mathlib, only: pack_matrix, unpack_matrix
    use oqp_linalg
    use printing, only: print_module_info

    implicit none

    character(len=*), parameter :: subroutine_name = "tdhf_sf_z_vector"

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

    integer :: nocca, nvira, noccb, nvirb
    integer :: nbf, nbf_tri
    integer :: iter
    integer :: i, j, ij  ! loop indices for UHF block
    real(kind=dp) :: cnvtol, scale_exch, scale_exch2
    logical :: roref = .false.
    logical :: uhfref = .false.   ! True for pure UHF (not ROHF)

    type(int2_compute_t) :: int2_driver
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
    integer :: nsocc, lzdim

  ! General data
    real(kind=dp) :: alpha, error

    logical :: dft
    integer :: scf_type, mol_mult

    ! tagarray
    real(kind=dp), contiguous, pointer :: &
      fock_a(:), mo_a(:,:), mo_energy_a(:), mo_energy_b(:), td_abxc(:,:), &
      fock_b(:), mo_b(:,:), &
      wao(:), td_p(:,:), td_t(:,:), &
      ta(:), tb(:), bvec_mo(:,:), sf_energies(:)
    character(len=*), parameter :: tags_alloc(3) = (/ character(len=80) :: &
      OQP_WAO, OQP_td_p, OQP_td_abxc /)
    character(len=*), parameter :: tags_required(9) = (/ character(len=80) :: &
      OQP_FOCK_A, OQP_E_MO_A, OQP_E_MO_B, OQP_VEC_MO_A, OQP_FOCK_B, OQP_VEC_MO_B, OQP_td_bvec_mo, OQP_td_t, &
      OQP_td_energies /)

    mol_mult = infos%mol_prop%mult

    scf_type = infos%control%scftype
    if (scf_type==3) roref = .true.
    if (scf_type==2) uhfref = .true.   ! Pure UHF (not constrained)

    dft = infos%control%hamilton == 20

  ! Files open
    open (unit=IW, file=infos%log_filename, position="append")
    call print_module_info('SF_TDHF_Z_Vector','Solving Z-Vector for SF-TDDFT')

  ! Load basis set
    basis => infos%basis
    basis%atoms => infos%atoms

    nbf = basis%nbf
    nbf_tri = nbf*(nbf+1)/2

    if (dft) call dft_initialize(infos, basis, molGrid)

  ! Parameter it should be inputed later
  ! convergence tolerance in the iterative TD-DFT step.
    cnvtol = infos%tddft%zvconv

    nocca = infos%mol_prop%nelec_A
    nvira = nbf-noccA
    noccb = infos%mol_prop%nelec_B
    nvirb = nbf-noccB
    nsocc = nocca-noccb

  ! Z-vector dimension depends on SCF type
    if (uhfref) then
      ! UHF: separate alpha and beta blocks
      lzdim = nocca*nvira + noccb*nvirb
    else
      ! ROHF: doc-socc-virt structure
      lzdim = noccb*(nsocc+nvira)+nsocc*nvira
    end if

  ! Allocate common arrays
    allocate(&
  ! for Z-vector
      xminv(lzdim), &
      rhs(lzdim), &
      lhs(lzdim), &
      xm(lzdim), &
      xk(lzdim), &
      pk(lzdim), &
      errv(lzdim), &
  ! for gradient
      hxa(nbf,nocca), &
      hxb(nbf,nbf), &
      tij(nocca,nocca), &
      tab(nvirb,nvirb), &
      ppija(nocca,nocca), &
      ppijb(noccb,noccb), &
      pa(nbf,nbf,2), &
   ! Allocate TDDFT variables
      fa(nbf,nbf), &
      fb(nbf,nbf), &
!   For scratch
      wrk1(nbf,nbf), &
      wrk2(nbf,nbf), &
      wrk3(nbf,nbf), &
      stat=ok, &
      source=0.0_dp)

  ! Allocate MO transformation arrays (different sizes for UHF vs ROHF)
    if (uhfref) then
      ! UHF: alpha-occ->alpha-virt, beta-occ->beta-virt
      allocate(ab1_MO_a(nocca, nvira), ab1_MO_b(noccb, nvirb), source=0.0_dp)
    else
      ! ROHF: alpha-occ->beta-virt for SF excitations
      allocate(ab1_MO_a(nocca, nvirb), ab1_MO_b(noccb, nvirb), source=0.0_dp)
    end if

    if( ok/=0 ) call show_message('Cannot allocate memory', with_abort)

    call infos%dat%remove_records(tags_alloc)

    call infos%dat%reserve_data(OQP_WAO, TA_TYPE_REAL64, nbf_tri, comment=OQP_WAO_comment)
    call infos%dat%reserve_data(OQP_td_p, TA_TYPE_REAL64, nbf_tri*2, (/ nbf_tri, 2 /), comment=OQP_td_p_comment)
    call infos%dat%reserve_data(OQP_td_abxc, TA_TYPE_REAL64, nbf*nbf, (/ nbf, nbf /), comment=OQP_td_abxc)

    call data_has_tags(infos%dat, tags_alloc, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_WAO, wao)
    call tagarray_get_data(infos%dat, OQP_td_p, td_p)
    call tagarray_get_data(infos%dat, OQP_td_abxc, td_abxc)

    call data_has_tags(infos%dat, tags_required, module_name, subroutine_name, WITH_ABORT)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a)
    call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_energy_a)
    call tagarray_get_data(infos%dat, OQP_E_MO_B, mo_energy_b)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b)
    call tagarray_get_data(infos%dat, OQP_td_bvec_mo, bvec_mo)
    call tagarray_get_data(infos%dat, OQP_td_t, td_t)
    call tagarray_get_data(infos%dat, OQP_td_energies, sf_energies)

    ta          => td_t(:,1)
    tb          => td_t(:,2)

! =============================================================================
!   STEP 1: Save unrelaxed density T and response vector (A-B)*X for target state
!   T_ij = -X_ia*X_ja (occ-occ), T_ab = X_ia*X_ib (virt-virt)
!   These are needed for gradient: dE/dR contains Tr[T * dH/dR]
! =============================================================================
    if (uhfref) then
      call sfdmat(bvec_mo(:,infos%tddft%target_state), td_abxc, mo_a, ta, tb, nocca, noccb, mo_b)
    else
      call sfdmat(bvec_mo(:,infos%tddft%target_state), td_abxc, mo_a, ta, tb, nocca, noccb)
    end if

  ! Initialize ERI calculations
    call int2_driver%init(basis, infos)
    call int2_driver%set_screening()

    write(iw,'(/1x,71("-")&
             &/19x,"SF-DFT ENERGY GRADIENT CALCULATION"&
             &/1x,71("-")/)')
    write(iw,fmt='(5x,a/&
                  &5x,16("-")/&
                  &5x,a,x,i0,x,f17.10,x,"Hartree"/&
                  &5x,a,x,e10.4/&
                  &5x,a,x,i0)') &
        'Z-vector options' &
      , 'Target state       is', infos%tddft%target_state, infos%mol_energy%energy+sf_energies(infos%tddft%target_state) &
      , 'Convergence        is', infos%tddft%zvconv &
      , 'Maximum iterations is', infos%control%maxit_zv
    call flush(iw)

    bvec(1:nbf,1:nbf,1:1) => td_abxc

  ! Prepare for ROHF
    ! Fock matrices A and B
    if( roref )then
        wrk1t(1:nbf*nbf) => wrk1
  !   Alpha
      call orthogonal_transform_sym(nbf, nbf, fock_a, mo_a, nbf, wrk1)
      call unpack_matrix(wrk1t, fa)

  !   Beta
      call orthogonal_transform_sym(nbf, nbf, fock_b, mo_b, nbf, wrk1)
      call unpack_matrix(wrk1t, fb)
    end if

! =============================================================================
!   STEP 2: Compute H[T] contribution to RHS
!   H[T]_ia = sum_jb T_jb * <ij||ab> + XC contribution
!   This is the first term in RHS = H[T] + H[X]*X
! =============================================================================
    call unpack_matrix(ta, pa(:,:,1))
    call unpack_matrix(tb, pa(:,:,2))

    scale_exch = 1.0_dp
    scale_exch2 = 1.0_dp
    if (dft) then
       scale_exch = infos%dft%HFscale    ! Reference HF exchange scaling
       scale_exch2 = infos%tddft%HFscale ! Response HF exchange scaling
    end if

    ! Compute (A+B)[T] = 2*(ij|ab) - (ia|jb) contracted with T
    int2_data = int2_tdgrd_data_t(d2=pa, &
            int_apb=.true., &
            int_amb=.false., &
            tamm_dancoff=.false., &
            scale_exchange=scale_exch)

    call int2_driver%run(int2_data, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%dft%cam_alpha, &
            beta=infos%dft%cam_beta,&
            mu=infos%dft%cam_mu)
    ab1 => int2_data%apb(:,:,:,1)

    ! Add DFT XC kernel contribution: f_xc[T]
    pa = pa*2
    call utddft_fxc(basis=basis, &
           molGrid=molGrid, &
           isVecs=.true., &
           wfa=MO_A, &
           wfb=MO_B, &
           fxa=ab1(:,:,1:1), &
           fxb=ab1(:,:,2:2), &
           dxa=pa(:,:,1:1), &
           dxb=pa(:,:,2:2), &
           nmtx=1, &
           threshold=0.0d0, &
           infos=infos)

    ! Transform H[T] from AO to MO basis: ab1_mo = C^T * ab1 * C
    call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)
    call mntoia(ab1(:,:,2), ab1_mo_b, mo_b, mo_b, noccb, noccb)

! =============================================================================
!   STEP 3: Compute H[X]*X contribution to RHS
!   H[X] = (A-B)[X] is the response matrix applied to excitation amplitudes
!   hxa, hxb store the result for occ and virt blocks
! =============================================================================
    call int2_data%clean()
    deallocate(int2_data)

    ! Compute (A-B)[X] using TDA integrals (exchange only for SF)
    int2_data = int2_td_data_t(d2=bvec, &
            int_apb=.false., &
            int_amb=.true., &
            tamm_dancoff=.true., &
            scale_exchange=scale_exch2)

    call int2_driver%run(int2_data, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%tddft%cam_alpha, &
            beta=infos%tddft%cam_beta,&
            mu=infos%tddft%cam_mu)
    ab2 => int2_data%amb(:,:,:,1)

    ! Transform (A-B)[X] to MO basis
    if (uhfref) then
      ! UHF: mixed alpha-beta transformation (alpha_occ x beta_virt)
      call dgemm('n', 'n', nbf, nbf, nbf, 1.0_dp, ab2(:,:,1), nbf, mo_b, nbf, 0.0_dp, wrk1, nbf)
      call dgemm('t', 'n', nbf, nbf, nbf, 1.0_dp, mo_a, nbf, wrk1, nbf, 0.0_dp, wrk2, nbf)
    else
      ! ROHF: symmetric transformation with alpha MOs
      call orthogonal_transform('n', nbf, mo_a, ab2(:,:,1), wrk2, wrk1)
    end if

    ! Expand X from packed to full matrix form
    call iatogen(bvec_mo(:,infos%tddft%target_state), wrk3, nocca, noccb)

    ! hxa(p,i) = 2 * sum_q H[X]_pq * X_qi  (for occupied index i)
    call dgemm('n', 't', nbf, nocca, nbf,  &
               2.0_dp, wrk2, nbf,  &
                       wrk3, nbf,  &
               0.0_dp, hxa,  nbf)

    ! hxb(p,q) = 2 * X^T * H[X]  (full matrix for virtual block)
    call dgemm('t', 'n', nbf, nbf, nocca,  &
               2.0_dp, wrk2, nbf,  &
                       wrk3, nbf,  &
               0.0_dp, hxb,  nbf)

! =============================================================================
!   STEP 4: Build unrelaxed difference density matrices T_ij and T_ab
!   T_ij = -sum_a X_ia * X_ja  (occupied-occupied block, electron depletion)
!   T_ab = +sum_i X_ia * X_ib  (virtual-virtual block, electron population)
! =============================================================================
    call dgemm('n', 't', nocca, nocca, nvirb,  &
              -1.0_dp, bvec_mo(:,infos%tddft%target_state), nocca,  &
                       bvec_mo(:,infos%tddft%target_state), nocca,  &
               0.0_dp, tij,     nocca)

    call dgemm('t', 'n', nvirb, nvirb, nocca,  &
               1.0_dp, bvec_mo(:,infos%tddft%target_state), nocca,  &
                       bvec_mo(:,infos%tddft%target_state), nocca,  &
               0.0_dp, tab,     nvirb)

! =============================================================================
!   STEP 5: Build full RHS of Z-vector equation
!   RHS_ia = H[T]_ia + (H[X]*X)_ia
!   This is the right-hand side of (A+B)*Z = -RHS
! =============================================================================
    if (uhfref) then
      call sfrcalc(rhs, ab1_mo_a, ab1_mo_b, hxa, hxb, nocca, noccb)
    else
      call sfrorhs(rhs, hxa, hxb, ab1_mo_a, ab1_mo_b, &
                   Tij, Tab, Fa, Fb, nocca, noccb)
    end if

    write(iw,'(/3x,25("-")&
             &/6x,"START Z-VECTOR LOOP"&
             &/3x,25("-")/)')

! =============================================================================
!   STEP 6: Initialize preconditioned conjugate gradient (PCG) solver
!   Preconditioner M_ia = 1/(epsilon_a - epsilon_i) for fast convergence
!   Initial guess: Z = M * RHS
! =============================================================================
    if (uhfref) then
      ! UHF: orbital energy differences for alpha and beta blocks separately
      call xecalc(xm(1:nocca*nvira), mo_energy_a, nocca)
      call xecalc(xm(nocca*nvira+1:lzdim), mo_energy_b, noccb)
      xminv = 1.0_dp / xm
    else
      ! ROHF: uses Fock matrix elements for doc-socc-virt structure
      call sfromcal(xm, xminv, mo_energy_a, fa, fb, nocca, noccb)
    end if

    call pcgrbpini(errv, pk, error, rhs, xminv, lhs)

    write(iw,'(" INITIAL ERROR =",3X,1P,E10.3,1X,"/",1P,E10.3)') error, cnvtol

! =============================================================================
!   STEP 7: PCG iteration loop
!   Each iteration: pk -> (A+B)*pk (via 2e integrals) -> update Z
!   Converge when ||residual||^2 < tolerance
! =============================================================================
    do iter = 1, infos%control%maxit_zv
      if (uhfref) then
        ! UHF: convert Z-vector to AO for alpha and beta separately
        ! Alpha: pk(1:nocca*nvira) -> pa(:,:,1)
        call sfgen(wrk1, pk(1:nocca*nvira), nocca, nbf)
        call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)
        ! Beta: pk(nocca*nvira+1:lzdim) -> pa(:,:,2)
        call sfgen(wrk2, pk(nocca*nvira+1:lzdim), noccb, nbf)
        call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)
      else
        ! ROHF: doc-socc-virt structure
        call sfrogen(wrk1, wrk2, pk, nocca, noccb)
!       Alpha
        call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)
!       Beta
        call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)
      end if

!     (A+B)*PK
      call int2_data%clean()
      deallocate(int2_data)
      int2_data = int2_tdgrd_data_t(d2=pa, &
              int_apb=.true., &
              int_amb=.false., &
              tamm_dancoff=.false., &
              scale_exchange=scale_exch)

      call int2_driver%run(int2_data, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%dft%cam_alpha, &
            beta=infos%dft%cam_beta,&
            mu=infos%dft%cam_mu)
      ab1 => int2_data%apb(:,:,:,1)

      !ab1 = ab1/2
      call symmetrize_matrix(pa(:,:,1), nbf)
      call symmetrize_matrix(pa(:,:,2), nbf)
      call utddft_fxc(basis=basis, &
             molGrid=molGrid, &
             isVecs=.true., &
             wfa=MO_A, &
             wfb=MO_B, &
             fxa=ab1(:,:,1:1), &
             fxb=ab1(:,:,2:2), &
             dxa=pa(:,:,1:1), &
             dxb=pa(:,:,2:2), &
             nmtx=1, &
             !threshold=1.0d-15, &
             threshold=0.0d0, &
             infos=infos)

      if (uhfref) then
        ! UHF: AO -> MO transformation for alpha and beta
        ! Alpha: (nocca, nvira)
        call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)
        ! Beta: (noccb, nvirb)
        call mntoia(ab1(:,:,2), ab1_mo_b, mo_b, mo_b, noccb, noccb)

        ! UHF: LHS = (A+B)*pk + (E_a - E_i)*pk
        call sflhs(lhs, pk, mo_energy_a, mo_energy_b, ab1_mo_a, ab1_mo_b, &
                   nocca, noccb, nvira, nvirb)
      else
        ! ROHF: doc-socc-virt structure
!       ALPHA: AO(M,N) -> MO(IA+) ... LPTMOA
        call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)

        call mntoia(ab1(:,:,2), ab1_mo_b, mo_a, mo_a, noccb, noccb)

        call sfrolhs(lhs, pk, mo_energy_a, fa, fb, ab1_mo_a, ab1_mo_b, &
                     nocca, noccb)
      end if

      alpha = 1.0_dp/dot_product(pk, lhs)

      xk = xk + pk * alpha
      errv = errv - alpha*lhs

      error = dot_product(errv, errv)
      write(iw,'(" ITER#",I2," ERROR =",3X,1P,E10.3,1X,"/",1P,E10.3)') &
        iter, error, cnvtol

      if (error<cnvtol) exit

      call pcgb(pk, errv, xminv)

    end do


! -----------------------------------------------
    if (error>cnvtol) then
       infos%mol_energy%Z_Vector_converged=.false.
       write(iw,'(/3x,24("-")&
             &/6x,"Z-Vector not converged"&
             &/3x,24("-")/)')
    else
       infos%mol_energy%Z_Vector_converged=.true.
       write(iw,'(/3x,24("-")&
             &/6x,"Z-Vector converged"&
             &/3x,24("-")/)')
    endif

    call flush(iw)

! =============================================================================
!   STEP 8: Construct relaxed density matrix P = T + Z
!   P contains both unrelaxed (T) and orbital relaxation (Z) contributions
!   Used for: gradient = Tr[P * dH/dR] + Tr[W * dS/dR] + 2e terms
! =============================================================================
    if (uhfref) then
      call sfpcal(wrk1, wrk2, tij, tab, xk, nocca, noccb, nvira, nvirb)
    else
      call sfropcal(wrk1, wrk2, tij, tab, xk, nocca, noccb)
    end if

    ! Transform P from MO to AO basis: P_AO = C * P_MO * C^T
    call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)
    call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)

! =============================================================================
!   STEP 9: Compute H[P] response for W matrix construction
!   H[P]_ij = sum_kl P_kl * <ik||jl> + XC contribution
!   This enters the energy-weighted density W_ij
! =============================================================================
    call int2_data%clean()
    deallocate(int2_data)
    int2_data = int2_tdgrd_data_t(d2=pa, &
            int_apb=.true., int_amb=.false., tamm_dancoff=.false., &
            scale_exchange=scale_exch)

    call int2_driver%run(int2_data, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%dft%cam_alpha, &
            beta=infos%dft%cam_beta,&
            mu=infos%dft%cam_mu)
    ab1 => int2_data%apb(:,:,:,1)

    ! Save relaxed density to output array (symmetrized and packed)
    call symmetrize_matrix(pa(:,:,1), nbf)
    call symmetrize_matrix(pa(:,:,2), nbf)
    call pack_matrix(pa(:,:,1), td_p(:,1))
    call pack_matrix(pa(:,:,2), td_p(:,2))
    td_p = 0.5_dp*td_p

    ! Add DFT XC kernel contribution to H[P]
    call utddft_fxc(basis=basis, &
           molGrid=molGrid, &
           isVecs=.true., &
           wfa=MO_A, &
           wfb=MO_B, &
           fxa=ab1(:,:,1:1), &
           fxb=ab1(:,:,2:2), &
           dxa=pa(:,:,1:1), &
           dxb=pa(:,:,2:2), &
           nmtx=1, &
           threshold=0.0d0, &
           infos=infos)

    ! Transform H[P] to MO basis: ppij = C^T * H[P] * C (occ-occ block)
    ! ppija_ij = H[P]_ij for alpha occupied orbitals
    call dgemm('n', 'n', nbf, nocca, nbf,  &
               1.0_dp, ab1(:,:,1), nbf,  &
                       mo_a, nbf,  &
               0.0_dp, wrk2, nbf)
    call dgemm('t', 'n', nocca, nocca, nbf,  &
               1.0_dp, mo_a,  nbf,  &
                       wrk2,  nbf,  &
               0.0_dp, ppija, nocca)

    ! ppijb_ij = H[P]_ij for beta occupied orbitals
    call dgemm('n', 'n', nbf, noccb, nbf,  &
               1.0_dp, ab1(:,:,2), nbf,  &
                       mo_b, nbf,  &
               0.0_dp, wrk2, nbf)
    call dgemm('t', 'n', noccb, noccb, nbf,  &
               1.0_dp, mo_b,  nbf,  &
                       wrk2,  nbf,  &
               0.0_dp, ppijb, noccb)

! =============================================================================
!   STEP 10: Compute energy-weighted density matrix W
!   W enforces orbital orthonormality: dE/dR += sum_pq W_pq * dS_pq/dR
!   W_ij contains occupied orbital contributions
!   W_ab contains virtual orbital contributions
!   W_ia contains occ-virt coupling from Z-vector
! =============================================================================
    wmo => wrk3
    wmo = 0

    if (uhfref) then
      ! UHF: separate alpha and beta W matrices
      block
        real(kind=dp), allocatable :: wmo_a(:,:), wmo_b(:,:), wao_a(:), wao_b(:)
        allocate(wmo_a(nbf,nbf), wmo_b(nbf,nbf), source=0.0_dp)
        allocate(wao_a(nbf_tri), wao_b(nbf_tri), source=0.0_dp)

        ! Compute W in MO basis using sfwcal
        call sfwcal(wmo_a, wmo_b, sf_energies(infos%tddft%target_state), &
                    mo_energy_a, mo_energy_b, fa, fb, &
                    bvec_mo(:,infos%tddft%target_state), xk, &
                    hxb, ppija, ppijb, nocca, noccb)

        ! Transform W_alpha: MO -> AO, then symmetrize and pack
        call orthogonal_transform('t', nbf, mo_a, wmo_a, wrk2, wrk1)
        call symmetrize_matrix(wrk2, nbf)
        wrk2 = wrk2 * 0.5_dp
        call pack_matrix(wrk2, wao_a)

        ! Transform W_beta: MO -> AO, then symmetrize and pack
        call orthogonal_transform('t', nbf, mo_b, wmo_b, wrk2, wrk1)
        call symmetrize_matrix(wrk2, nbf)
        wrk2 = wrk2 * 0.5_dp
        call pack_matrix(wrk2, wao_b)

        ! Total W in AO = W_alpha + W_beta
        wao = wao_a + wao_b

        deallocate(wmo_a, wmo_b, wao_a, wao_b)
      end block
    else
      ! ROHF: single W matrix using shared MO basis
      call sfrowcal(wmo,sf_energies(infos%tddft%target_state), &
                    mo_energy_a, mo_energy_b, fa, fb, bvec_mo(:,infos%tddft%target_state), xk, &
                    hxa, hxb, ppija, ppijb, &
                    nocca, noccb)

      ! Transform W: MO -> AO, then symmetrize and pack
      call orthogonal_transform('t', nbf, mo_a, wmo, wrk2, wrk1)
      call symmetrize_matrix(wrk2, nbf)
      wrk2 = wrk2 * 0.5_dp
      call pack_matrix(wrk2, wao)
    end if

    call int2_driver%clean()

    if (dft) call dftclean(infos)

    call measure_time(print_total=1, log_unit=iw)
    close(iw)

  end subroutine tdhf_sf_z_vector

end module tdhf_sf_z_vector_mod
