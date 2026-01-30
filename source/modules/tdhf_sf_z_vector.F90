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
      sfrcalc, xecalc, sfuesum, sfugen, umrsfrowcal
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

    ! Save unrelaxed density matrices and the `b=A*x` vector for target state
    ! For UHF: pass mo_b for correct beta virtual transformation
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

!   ALPHA: AO(M,N) -> MO(IA+)
    call mntoia(ab1(:,:,1), ab1_mo_a, mo_a, mo_a, nocca, nocca)

    call mntoia(ab1(:,:,2), ab1_mo_b, mo_b, mo_b, noccb, noccb)

  ! Initialize ERI calculations
    call int2_data%clean()
    deallocate(int2_data)
    int2_data = int2_td_data_t(d2=bvec, &
            int_apb=.false., &
            int_amb=.true., &  ! MUST be true - ab2 uses amb!
            tamm_dancoff=.true., &
            scale_exchange=scale_exch2)

    call int2_driver%run(int2_data, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%tddft%cam_alpha, &
            beta=infos%tddft%cam_beta,&
            mu=infos%tddft%cam_mu)
    ab2 => int2_data%amb(:,:,:,1)

    if (uhfref) then
      ! UHF: H[X] connects alpha and beta spaces
      ! wrk2 = mo_a^T * ab2 * mo_b (not mo_a^T * ab2 * mo_a!)
      call dgemm('n', 'n', nbf, nbf, nbf, 1.0_dp, ab2(:,:,1), nbf, mo_b, nbf, 0.0_dp, wrk1, nbf)
      call dgemm('t', 'n', nbf, nbf, nbf, 1.0_dp, mo_a, nbf, wrk1, nbf, 0.0_dp, wrk2, nbf)
    else
      ! ROHF: standard symmetric transformation
      call orthogonal_transform('n', nbf, mo_a, ab2(:,:,1), wrk2, wrk1)
    end if

    call iatogen(bvec_mo(:,infos%tddft%target_state), wrk3, nocca, noccb)

    call dgemm('n', 't', nbf, nocca, nbf,  &
               2.0_dp, wrk2, nbf,  &
                       wrk3, nbf,  &
               0.0_dp, hxa,  nbf)

    call dgemm('t', 'n', nbf, nbf, nocca,  &
               2.0_dp, wrk2, nbf,  &
                       wrk3, nbf,  &
               0.0_dp, hxb,  nbf)

!   Unrelaxed difference density matries T_ij and T_ab
!     Ta(i+,j+):= -X(i+,a-)*X(j+,a-) for singlet and triplet
    call dgemm('n', 't', nocca, nocca, nvirb,  &
              -1.0_dp, bvec_mo(:,infos%tddft%target_state), nocca,  &
                       bvec_mo(:,infos%tddft%target_state), nocca,  &
               0.0_dp, tij,     nocca)

    ! Tb(a-,b-):= X(i+,a-)*X(i+,b-) for singlet and triplet
    call dgemm('t', 'n', nvirb, nvirb, nocca,  &
               1.0_dp, bvec_mo(:,infos%tddft%target_state), nocca,  &
                       bvec_mo(:,infos%tddft%target_state), nocca,  &
               0.0_dp, tab,     nvirb)

    if (uhfref) then
      ! UHF: simple alpha+beta blocks
      call sfrcalc(rhs, ab1_mo_a, ab1_mo_b, hxa, hxb, nocca, noccb)
    else
      ! ROHF: doc-socc-virt structure
      call sfrorhs(rhs, hxa, hxb, ab1_mo_a, ab1_mo_b, &
                   Tij, Tab, Fa, Fb, nocca, noccb)
    end if

    write(iw,'(/3x,25("-")&
             &/6x,"START Z-VECTOR LOOP"&
             &/3x,25("-")/)')

    if (uhfref) then
      ! UHF: XM = orbital energy differences for alpha and beta
      call xecalc(xm(1:nocca*nvira), mo_energy_a, nocca)
      call xecalc(xm(nocca*nvira+1:lzdim), mo_energy_b, noccb)
      xminv = 1.0_dp / xm
    else
      ! ROHF: uses Fock matrices
      call sfromcal(xm, xminv, mo_energy_a, fa, fb, nocca, noccb)
    end if

    call pcgrbpini(errv, pk, error, rhs, xminv, lhs)

    write(iw,'(" INITIAL ERROR =",3X,1P,E10.3,1X,"/",1P,E10.3)') error, cnvtol

! -----------------------------------------------

    do iter = 1, infos%control%maxit_zv
      if (uhfref) then
        ! UHF: convert Z-vector to AO for alpha and beta separately
        ! Alpha: pk(1:nocca*nvira) -> pa(:,:,1)
        call sfugen(wrk1, pk(1:nocca*nvira), nocca, nbf)
        call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)
        ! Beta: pk(nocca*nvira+1:lzdim) -> pa(:,:,2)
        call sfugen(wrk2, pk(nocca*nvira+1:lzdim), noccb, nbf)
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
        ! First copy (A+B)*pk part to lhs
        lhs = 0.0_dp
        ! Alpha part
        do j = 1, nvira
          do i = 1, nocca
            ij = (j-1)*nocca + i
            lhs(ij) = ab1_mo_a(i, j)
          end do
        end do
        ! Beta part
        do j = 1, nvirb
          do i = 1, noccb
            ij = nocca*nvira + (j-1)*noccb + i
            lhs(ij) = ab1_mo_b(i, j)
          end do
        end do
        ! Add diagonal (E_a - E_i)*pk
        call sfuesum(lhs(1:nocca*nvira), mo_energy_a, pk(1:nocca*nvira), nocca)
        call sfuesum(lhs(nocca*nvira+1:lzdim), mo_energy_b, pk(nocca*nvira+1:lzdim), noccb)
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

    ! ============ DEBUG: Z-VECTOR SOLUTION ============
    write(iw,'(/,"======== DEBUG: Z-VECTOR SOLUTION ========")')
    write(iw,'("[ZV_SOLN] xk: norm=",ES12.4," abssum=",ES12.4)') &
          sqrt(sum(xk**2)), sum(abs(xk))
    write(iw,'("[ZV_SOLN] xk(1:6)=",6ES11.3)') xk(1:6)
    write(iw,'("[ZV_SOLN] xk alpha part (1:",I0,"): norm=",ES12.4)') &
          nocca*nvira, sqrt(sum(xk(1:nocca*nvira)**2))
    write(iw,'("[ZV_SOLN] xk beta part (",I0,":",I0,"): norm=",ES12.4)') &
          nocca*nvira+1, lzdim, sqrt(sum(xk(nocca*nvira+1:lzdim)**2))
    call flush(iw)

    if (uhfref) then
      ! UHF: construct density from T + Z
      ! Alpha: wrk1(i,j) = Tij(i,j), wrk1(i,a) = xk(i,a)  [no factor of 0.5!]
      wrk1 = 0.0_dp
      wrk1(1:nocca, 1:nocca) = tij
      do j = 1, nvira
        do i = 1, nocca
          ij = (j-1)*nocca + i
          wrk1(i, nocca+j) = xk(ij)
        end do
      end do
      ! Beta: wrk2(a,b) = Tab(a,b), wrk2(i,a) = xk(nconfa+i,a)  [no factor of 0.5!]
      wrk2 = 0.0_dp
      do j = 1, nvirb
        do i = 1, nvirb
          wrk2(noccb+i, noccb+j) = tab(i, j)
        end do
      end do
      do j = 1, nvirb
        do i = 1, noccb
          ij = nocca*nvira + (j-1)*noccb + i
          wrk2(i, noccb+j) = xk(ij)
        end do
      end do

      ! ============ DEBUG: DENSITY MATRICES IN MO ============
      write(iw,'(/,"======== DEBUG: DENSITY IN MO ========")')
      write(iw,'("[DEN_MO_A] wrk1: norm=",ES12.4," abssum=",ES12.4)') &
            sqrt(sum(wrk1**2)), sum(abs(wrk1))
      write(iw,'("[DEN_MO_A] wrk1(1:3,1:3)=",3ES11.3)') wrk1(1,1:3), wrk1(2,1:3), wrk1(3,1:3)
      write(iw,'("[DEN_MO_A] Tij block norm=",ES12.4)') sqrt(sum(tij**2))
      write(iw,'("[DEN_MO_B] wrk2: norm=",ES12.4," abssum=",ES12.4)') &
            sqrt(sum(wrk2**2)), sum(abs(wrk2))
      write(iw,'("[DEN_MO_B] wrk2(1:3,1:3)=",3ES11.3)') wrk2(1,1:3), wrk2(2,1:3), wrk2(3,1:3)
      write(iw,'("[DEN_MO_B] Tab block norm=",ES12.4)') sqrt(sum(tab**2))

      ! ============ DEBUG: Separate T_AO and Z_AO ============
      ! Use wrk3 as temporary workspace for T_ao/Z_ao debug
      ! Alpha T: just occ-occ part (Tij)
      wrk3 = 0.0_dp
      wrk3(1:nocca, 1:nocca) = tij
      call orthogonal_transform('t', nbf, mo_a, wrk3, pa(:,:,1), wrk2)
      write(iw,'("[T_AO_A] norm=",ES12.4)') sqrt(sum(pa(:,:,1)**2))

      ! Alpha Z: just occ-vir part (xk, no factor of 0.5)
      wrk3 = 0.0_dp
      do j = 1, nvira
        do i = 1, nocca
          wrk3(i, nocca+j) = xk((j-1)*nocca + i)
        end do
      end do
      call orthogonal_transform('t', nbf, mo_a, wrk3, pa(:,:,2), wrk2)
      write(iw,'("[Z_AO_A] norm=",ES12.4)') sqrt(sum(pa(:,:,2)**2))

      ! Beta T: just vir-vir part (Tab)
      wrk3 = 0.0_dp
      do j = 1, nvirb
        do i = 1, nvirb
          wrk3(noccb+i, noccb+j) = tab(i, j)
        end do
      end do
      call orthogonal_transform('t', nbf, mo_b, wrk3, pa(:,:,1), wrk2)
      write(iw,'("[T_AO_B] norm=",ES12.4)') sqrt(sum(pa(:,:,1)**2))

      ! Beta Z: just occ-vir part (xk, no factor of 0.5)
      wrk3 = 0.0_dp
      do j = 1, nvirb
        do i = 1, noccb
          wrk3(i, noccb+j) = xk(nocca*nvira + (j-1)*noccb + i)
        end do
      end do
      call orthogonal_transform('t', nbf, mo_b, wrk3, pa(:,:,2), wrk2)
      write(iw,'("[Z_AO_B] norm=",ES12.4)') sqrt(sum(pa(:,:,2)**2))
      call flush(iw)
      ! ============ END DEBUG ============

      ! Need to restore wrk1, wrk2 for actual calculation (they got overwritten)
      ! Alpha: wrk1(i,j) = Tij(i,j), wrk1(i,a) = xk(i,a)
      wrk1 = 0.0_dp
      wrk1(1:nocca, 1:nocca) = tij
      do j = 1, nvira
        do i = 1, nocca
          wrk1(i, nocca+j) = xk((j-1)*nocca + i)
        end do
      end do
      ! Beta: wrk2(a,b) = Tab(a,b), wrk2(i,a) = xk(nconfa+i,a)
      wrk2 = 0.0_dp
      do j = 1, nvirb
        do i = 1, nvirb
          wrk2(noccb+i, noccb+j) = tab(i, j)
        end do
      end do
      do j = 1, nvirb
        do i = 1, noccb
          wrk2(i, noccb+j) = xk(nocca*nvira + (j-1)*noccb + i)
        end do
      end do
    else
      ! ROHF: doc-socc-virt structure
      call sfropcal(wrk1, wrk2, tij, tab, xk, nocca, noccb)
    end if

 !  Update density for alpha
    call orthogonal_transform('t', nbf, mo_a, wrk1, pa(:,:,1), wrk3)

 !  Update density for beta
    call orthogonal_transform('t', nbf, mo_b, wrk2, pa(:,:,2), wrk3)

    ! ============ DEBUG: DENSITY IN AO (before symmetrize) ============
    write(iw,'(/,"======== DEBUG: DENSITY IN AO ========")')
    write(iw,'("[DEN_AO_A] pa(:,:,1): norm=",ES12.4," abssum=",ES12.4," trace=",ES12.4)') &
          sqrt(sum(pa(:,:,1)**2)), sum(abs(pa(:,:,1))), sum([(pa(i,i,1), i=1,nbf)])
    write(iw,'("[DEN_AO_A] (1:3,1:3)=",3ES11.3)') pa(1,1:3,1), pa(2,1:3,1), pa(3,1:3,1)
    write(iw,'("[DEN_AO_B] pa(:,:,2): norm=",ES12.4," abssum=",ES12.4," trace=",ES12.4)') &
          sqrt(sum(pa(:,:,2)**2)), sum(abs(pa(:,:,2))), sum([(pa(i,i,2), i=1,nbf)])
    write(iw,'("[DEN_AO_B] (1:3,1:3)=",3ES11.3)') pa(1,1:3,2), pa(2,1:3,2), pa(3,1:3,2)
    ! Check symmetry of density: ||A - A^T||
    write(iw,'("[DEN_AO_A] asymmetry=",ES12.4)') &
          sqrt(sum((pa(:,:,1) - transpose(pa(:,:,1)))**2))
    write(iw,'("[DEN_AO_B] asymmetry=",ES12.4)') &
          sqrt(sum((pa(:,:,2) - transpose(pa(:,:,2)))**2))
    ! NOTE: Do NOT symmetrize! GAMESS uses asymmetric density for UTD2E
    call flush(iw)

    call int2_data%clean()
    deallocate(int2_data)
    int2_data = int2_tdgrd_data_t(d2=pa, &
            int_apb=.true., int_amb=.false., tamm_dancoff=.false., &
            scale_exchange=scale_exch)

    ! Print input to int2_driver
    write(iw,'("[int2_data%d2(1)] norm=",ES12.4)') sqrt(sum(int2_data%d2(:,:,1)**2))
    write(iw,'("[int2_data%d2(2)] norm=",ES12.4)') sqrt(sum(int2_data%d2(:,:,2)**2))
    call flush(iw)

    call int2_driver%run(int2_data, &
            cam=dft.and.infos%dft%cam_flag, &
            alpha=infos%dft%cam_alpha, &
            beta=infos%dft%cam_beta,&
            mu=infos%dft%cam_mu)
    ab1 => int2_data%apb(:,:,:,1)

    ! ============ DEBUG: (A+B)[P] 2e integrals ============
    write(iw,'(/,"======== DEBUG: (A+B)[P] ========")')
    write(iw,'("[APB_A] ab1(:,:,1): norm=",ES12.4," abssum=",ES12.4," trace=",ES12.4)') &
          sqrt(sum(ab1(:,:,1)**2)), sum(abs(ab1(:,:,1))), sum([(ab1(i,i,1), i=1,nbf)])
    write(iw,'("[APB_A] (1:3,1:3)=",3ES11.3)') ab1(1,1:3,1), ab1(2,1:3,1), ab1(3,1:3,1)
    write(iw,'("[APB_B] ab1(:,:,2): norm=",ES12.4," abssum=",ES12.4," trace=",ES12.4)') &
          sqrt(sum(ab1(:,:,2)**2)), sum(abs(ab1(:,:,2))), sum([(ab1(i,i,2), i=1,nbf)])
    write(iw,'("[APB_B] (1:3,1:3)=",3ES11.3)') ab1(1,1:3,2), ab1(2,1:3,2), ab1(3,1:3,2)
    call flush(iw)

    call symmetrize_matrix(pa(:,:,1), nbf)
    call symmetrize_matrix(pa(:,:,2), nbf)
    call pack_matrix(pa(:,:,1), td_p(:,1))
    call pack_matrix(pa(:,:,2), td_p(:,2))
    td_p = 0.5_dp*td_p

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

    ! ============ DEBUG: PPIJ matrices ============
    write(iw,'(/,"======== DEBUG: PPIJ ========")')
    write(iw,'("[PPIJA] ppija: norm=",ES12.4," abssum=",ES12.4," trace=",ES12.4)') &
          sqrt(sum(ppija**2)), sum(abs(ppija)), sum([(ppija(i,i), i=1,nocca)])
    write(iw,'("[PPIJA] (1:3,1:3)=",3ES11.3)') ppija(1,1:3), ppija(2,1:3), ppija(3,1:3)
    write(iw,'("[PPIJB] ppijb: norm=",ES12.4," abssum=",ES12.4," trace=",ES12.4)') &
          sqrt(sum(ppijb**2)), sum(abs(ppijb)), sum([(ppijb(i,i), i=1,noccb)])
    write(iw,'("[PPIJB] (1:3,1:3)=",3ES11.3)') ppijb(1,1:3), ppijb(2,1:3), ppijb(3,1:3)
    call flush(iw)

!   Calculate W (in MO basis)
    wmo => wrk3
    wmo = 0

    if (uhfref) then
      ! UHF: Use separate alpha/beta W matrices (GAMESS SFWCALC structure)
      block
        real(kind=dp), allocatable :: wmo_a(:,:), wmo_b(:,:), wao_a(:), wao_b(:)
        allocate(wmo_a(nbf,nbf), wmo_b(nbf,nbf), source=0.0_dp)
        allocate(wao_a(nbf_tri), wao_b(nbf_tri), source=0.0_dp)

        call umrsfrowcal(wmo_a, wmo_b, sf_energies(infos%tddft%target_state), &
                         mo_energy_a, mo_energy_b, fa, fb, &
                         bvec_mo(:,infos%tddft%target_state), xk, &
                         hxb, ppija, ppijb, nocca, noccb)

        ! ============ DEBUG: W matrices (MO basis) ============
        write(iw,'(/,"======== DEBUG: W MATRIX (UHF separate) ========")')
        write(iw,'("[WMOA] wmo_a: norm=",ES12.4)') sqrt(sum(wmo_a**2))
        write(iw,'("[WMOB] wmo_b: norm=",ES12.4)') sqrt(sum(wmo_b**2))
        ! DEBUG: Check MO coefficients
        write(iw,'("[MO_A] mo_a: norm=",ES12.4)') sqrt(sum(mo_a**2))
        write(iw,'("[MO_B] mo_b: norm=",ES12.4)') sqrt(sum(mo_b**2))
        call flush(iw)

        ! Transform alpha W: WMO_A -> WAO_A using mo_a
        ! DEBUG: Before symmetrization
        call orthogonal_transform('t', nbf, mo_a, wmo_a, wrk2, wrk1)
        write(iw,'("[WAOA_BEFORE_SYM] norm=",ES12.4)') sqrt(sum(wrk2**2))
        call symmetrize_matrix(wrk2, nbf)
        wrk2 = wrk2 * 0.5_dp  ! symmetrize_matrix does A+A^T, need 0.5 for (A+A^T)/2
        call pack_matrix(wrk2, wao_a)
        write(iw,'("[WAOA_FULL] norm=",ES12.4)') sqrt(sum(wrk2**2))
        write(iw,'("[WAOA_PACKED] norm=",ES12.4)') sqrt(sum(wao_a**2))

        ! Transform beta W: WMO_B -> WAO_B using mo_b
        call orthogonal_transform('t', nbf, mo_b, wmo_b, wrk2, wrk1)
        call symmetrize_matrix(wrk2, nbf)
        wrk2 = wrk2 * 0.5_dp  ! symmetrize_matrix does A+A^T, need 0.5 for (A+A^T)/2
        call pack_matrix(wrk2, wao_b)
        write(iw,'("[WAOB_FULL] norm=",ES12.4)') sqrt(sum(wrk2**2))
        write(iw,'("[WAOB_PACKED] norm=",ES12.4)') sqrt(sum(wao_b**2))

        ! Combine: wao = wao_a + wao_b (both contribute to gradient)
        ! No extra scaling - the 0.5 was already applied in symmetrization
        wao = wao_a + wao_b

        write(iw,'("[WAO_COMBINED] wao: norm=",ES12.4)') sqrt(sum(wao**2))
        call flush(iw)

        deallocate(wmo_a, wmo_b, wao_a, wao_b)
      end block
    else
      ! ROHF: Use original sfrowcal (single W matrix)
      call sfrowcal(wmo,sf_energies(infos%tddft%target_state), &
                    mo_energy_a, mo_energy_b, fa, fb, bvec_mo(:,infos%tddft%target_state), xk, &
                    hxa, hxb, ppija, ppijb, &
                    nocca, noccb)

      ! ============ DEBUG: W matrix ============
      write(iw,'(/,"======== DEBUG: W MATRIX ========")')
      write(iw,'("[WMO] wmo: norm=",ES12.4," abssum=",ES12.4," trace=",ES12.4)') &
            sqrt(sum(wmo**2)), sum(abs(wmo)), sum([(wmo(i,i), i=1,nbf)])
      write(iw,'("[WMO] (1:3,1:3)=",3ES11.3)') wmo(1,1:3), wmo(2,1:3), wmo(3,1:3)
      call flush(iw)

      call orthogonal_transform('t', nbf, mo_a, wmo, wrk2, wrk1)
      call symmetrize_matrix(wrk2, nbf)
      call pack_matrix(wrk2, wao)

      ! ============ DEBUG: W in AO (before scaling) ============
      write(iw,'(/,"======== DEBUG: W in AO ========")')
      write(iw,'("[WAO_BEFORE] wao: norm=",ES12.4," abssum=",ES12.4)') &
            sqrt(sum(wao**2)), sum(abs(wao))
      write(iw,'("[WAO_BEFORE] (1:6)=",6ES11.3)') wao(1:6)

      wao = wao*0.25_dp  ! ROHF factor

      write(iw,'("[WAO_AFTER] wao: norm=",ES12.4," abssum=",ES12.4)') &
            sqrt(sum(wao**2)), sum(abs(wao))
      write(iw,'("[WAO_AFTER] (1:6)=",6ES11.3)') wao(1:6)
      call flush(iw)
    end if

    call int2_driver%clean()

    if (dft) call dftclean(infos)

    call measure_time(print_total=1, log_unit=iw)
    close(iw)

  end subroutine tdhf_sf_z_vector

end module tdhf_sf_z_vector_mod
