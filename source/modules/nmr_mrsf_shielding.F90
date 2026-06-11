module nmr_mrsf_shielding_mod

  use precision, only: dp
  implicit none

  character(len=*), parameter :: module_name = "nmr_mrsf_shielding_mod"

  private
  public nmr_mrsf_shielding

contains

  subroutine nmr_mrsf_shielding_C(c_handle) bind(C, name="nmr_mrsf_shielding")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call nmr_mrsf_shielding(inf)
  end subroutine nmr_mrsf_shielding_C

!> @brief MRSF-TDDFT state NMR shielding -- Gate 2 prototype (state-density
!>  diamagnetic shielding + frozen-reference paramagnetic approximation).
!> @details Implements the first stage of the MRSF-NMR roadmap
!>  (NMR_MRSF_IMPLEMENTATION.md, Gate 2; theory: Overleaf THEORY.tex):
!>
!>  * DIAMAGNETIC term (EXACT for the state): by the master formula
!>    sigma_dia^K = Tr[D_rel^K d2h/dB dm], the GIAO diamagnetic kernels are
!>    density-independent AO objects, so contracting them with the gradient
!>    relaxed state density D_rel^K = DM_A + DM_B + td_p(:,1) + td_p(:,2)
!>    (the exact combination used by electric_moments_excited) gives the
!>    exact state diamagnetic tensor.
!>
!>  * PARAMAGNETIC term (APPROXIMATION, clearly labeled): the true state
!>    paramagnetic response requires the three magnetic response systems of
!>    THEORY.tex Sec. 7 (orbital + amplitude + multiplier), implemented in
!>    Gates 3-4.  This prototype reports the FROZEN-REFERENCE para: the
!>    ground-state GIAO CPHF response of the ROHF triplet reference (exactly
!>    the open-shell path of nmr_giao_shielding.F90).  In the zero-amplitude
!>    limit (X -> 0, z -> 0) this prototype therefore reduces term-by-term to
!>    the validated ground-state GIAO shielding -- the Gate 2 [CHECK].
!>
!>  Hard guards (no silent fallback to ground-state NMR): requires an ROHF
!>  triplet reference with the MRSF Z-vector relaxed density OQP::td_p
!>  present (i.e. a runtype=grad MRSF run); aborts for CAM/meta-GGA/ECP as
!>  the underlying magnetic kernels do.
  subroutine nmr_mrsf_shielding(infos)
    use io_constants, only: iw
    use oqp_tagarray_driver
    use basis_tools, only: basis_set
    use messages, only: show_message, with_abort
    use types, only: information
    use constants, only: tol_int
    use int1, only: giao_h10_core, giao_overlap_derivative, nmr_dia_shielding, &
                    giao_a11part_corr, giao_a01gp_contract
    use nmr_giao_debug_mod, only: giao_h10_twoe_matrix
    use nmr_giao_shielding_mod, only: giao_para_channel, semicanon_orbitals
    use dft, only: dft_initialize
    use mod_dft_molgrid, only: dft_grid_t
    use mod_dft_gridint_giao, only: giao_vxc

    implicit none

    character(len=*), parameter :: subroutine_name = "nmr_mrsf_shielding"
    real(kind=dp), parameter :: ALPHA = 1.0d0/137.035999084d0
    real(kind=dp), parameter :: a2ppm = ALPHA*ALPHA*1.0d6
    ! GIAO diamagnetic/London sign conventions; must stay identical to
    ! nmr_giao_shielding.F90 (single source of truth pending the Gate-0
    ! prefactor audit, THEORY.tex open item 5).
    real(kind=dp), parameter :: SG = -1.0d0, SC = 1.0d0, SA = -0.5d0
    real(kind=dp), parameter :: SX = -1.0d0

    type(information), target, intent(inout) :: infos
    type(basis_set), pointer :: basis

    integer :: nbf, nbf2, nat, nocc_a, nocc_b, nmo
    integer :: i, j, c, t, s, iat, target_state
    integer(4) :: status
    logical :: is_dft, iw_open
    real(kind=dp) :: tol, scale_exch, trg0, trc

    real(kind=dp), allocatable :: dens_state(:), dens_ref(:)   ! packed (nbf2)
    real(kind=dp), allocatable :: coords(:,:), zq(:)
    real(kind=dp), allocatable :: gdia0(:,:,:), corrpre(:,:,:), a01(:,:,:)
    real(kind=dp), allocatable :: sig_dia(:,:,:), sig_dia_ref(:,:,:)
    real(kind=dp), allocatable :: h10p(:,:), s10p(:,:)
    real(kind=dp), allocatable :: h1ao(:,:,:), h1ao_b(:,:,:), s1ao(:,:,:)
    real(kind=dp), allocatable :: twoe(:,:,:), twoe2(:,:,:), vj(:,:,:), vk(:,:,:), vkb(:,:,:)
    real(kind=dp), allocatable :: vxa(:,:,:), vxb(:,:,:)
    real(kind=dp), allocatable :: dm(:,:), dm_b(:,:)
    real(kind=dp), allocatable :: sig_u(:,:,:), sig_c(:,:,:)
    real(kind=dp), allocatable :: ca_sc(:,:), cb_sc(:,:), ea_sc(:), eb_sc(:)
    type(dft_grid_t) :: molGrid
    integer :: mxAngMom

    real(kind=dp), contiguous, pointer :: dmat_a(:), dmat_b(:)
    real(kind=dp), contiguous, pointer :: mo_a(:,:), mo_b(:,:), mo_e(:)
    real(kind=dp), contiguous, pointer :: fock_a(:), fock_b(:)
    real(kind=dp), contiguous, pointer :: td_p(:,:), td_energies(:)
    real(kind=dp), contiguous, pointer :: nmrout(:,:)

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf = basis%nbf
    nbf2 = nbf*(nbf+1)/2
    nat = ubound(basis%atoms%zn,1)
    tol = log(10.0d0)*tol_int

    ! Connect the log unit early so guard aborts land in the log.
    inquire(unit=iw, opened=iw_open)
    if (.not. iw_open) open(unit=iw, file=infos%log_filename, position="append")

    ! ------------------------------------------------------------------
    ! Hard guards: this is a state-specific property on an MRSF relaxed
    ! density; refuse anything else instead of silently degrading.
    ! ------------------------------------------------------------------
    if (infos%control%scftype /= 3) then
      call show_message('MRSF-NMR requires an ROHF reference (scf.type=rohf)', &
        with_abort)
    end if
    if (infos%mol_prop%mult /= 3) then
      call show_message('MRSF-NMR requires a triplet reference &
        &(scf.multiplicity=3)', with_abort)
    end if

    is_dft = infos%control%hamilton == 20
    scale_exch = 1.0d0
    if (is_dft) scale_exch = infos%dft%HFscale
    if (is_dft) then
      if (infos%dft%cam_flag) then
        call show_message('MRSF-NMR with range-separated (CAM) functionals &
          &is not implemented', with_abort)
      end if
      if (infos%functional%needtau) then
        call show_message('MRSF-NMR with meta-GGA (tau-dependent) functionals &
          &is not implemented', with_abort)
      end if
    end if
    if (allocated(infos%basis%ecp_zn_num)) then
      if (any(infos%basis%ecp_zn_num /= 0)) then
        call show_message('MRSF-NMR with ECP basis sets is not implemented', &
          with_abort)
      end if
    end if

    ! The MRSF relaxed difference density must exist (set by the MRSF
    ! Z-vector for the current target state).  Its absence means no MRSF
    ! response was solved -- abort, never fall back to ground-state NMR.
    call tagarray_get_data(infos%dat, OQP_td_p, td_p, status)
    if (status /= 0) then
      call show_message('MRSF-NMR: relaxed MRSF density (OQP::td_p) not &
        &found. Run runtype=grad with tdhf.type=mrsf so the MRSF Z-vector &
        &is solved for the target state; ground-state NMR fallback is &
        &deliberately not provided.', with_abort)
    end if

    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a, status)
    call check_status(status, module_name, subroutine_name, OQP_DM_A)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b, status)
    call check_status(status, module_name, subroutine_name, OQP_DM_B)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_A, mo_a, status)
    call check_status(status, module_name, subroutine_name, OQP_VEC_MO_A)
    call tagarray_get_data(infos%dat, OQP_VEC_MO_B, mo_b, status)
    call check_status(status, module_name, subroutine_name, OQP_VEC_MO_B)
    call tagarray_get_data(infos%dat, OQP_E_MO_A, mo_e, status)
    call check_status(status, module_name, subroutine_name, OQP_E_MO_A)
    call tagarray_get_data(infos%dat, OQP_FOCK_A, fock_a, status)
    call check_status(status, module_name, subroutine_name, OQP_FOCK_A)
    call tagarray_get_data(infos%dat, OQP_FOCK_B, fock_b, status)
    call check_status(status, module_name, subroutine_name, OQP_FOCK_B)

    target_state = infos%tddft%target_state
    nocc_a = int(infos%mol_prop%nelec_A)
    nocc_b = int(infos%mol_prop%nelec_B)
    nmo = size(mo_e)

    allocate(coords(3,nat), zq(nat))
    do iat = 1, nat
      coords(:,iat) = basis%atoms%xyz(:,iat)
    end do
    zq = infos%atoms%zn - infos%basis%ecp_zn_num

    ! ------------------------------------------------------------------
    ! Relaxed state density (packed): D_rel^K = DM_A + DM_B + td_p1 + td_p2
    ! (same combination as electric_moments_excited).  dens_ref is the
    ! reference (ROHF) density, kept for the Delta-dia diagnostic and the
    ! zero-amplitude check.
    ! ------------------------------------------------------------------
    allocate(dens_state(nbf2), dens_ref(nbf2))
    dens_ref = dmat_a + dmat_b
    dens_state = dens_ref + td_p(:,1) + td_p(:,2)

    ! ------------------------------------------------------------------
    ! Diamagnetic tensor (GIAO; EXACT for the state).  Assembly identical
    ! to nmr_giao_shielding.F90: a11part (trace-corrected) + a01gp.
    ! ------------------------------------------------------------------
    allocate(gdia0(3,3,nat), corrpre(3,3,nat), a01(3,3,nat), &
             sig_dia(3,3,nat), sig_dia_ref(3,3,nat), source=0.0d0)

    call dia_tensor(dens_state, sig_dia)
    call dia_tensor(dens_ref,   sig_dia_ref)

    ! ------------------------------------------------------------------
    ! Paramagnetic tensor: FROZEN-REFERENCE approximation (Gate 2).
    ! Recipe = open-shell branch of nmr_giao_shielding.F90 with the ROHF
    ! reference densities and semicanonical spin orbitals.
    ! ------------------------------------------------------------------
    allocate(h10p(nbf2,3), s10p(nbf2,3), source=0.0d0)
    call giao_h10_core(basis, coords, zq, h10p, debug=.false., logtol=tol)
    call giao_overlap_derivative(basis, s10p, debug=.false., logtol=tol)

    allocate(h1ao(nbf,nbf,3), h1ao_b(nbf,nbf,3), s1ao(nbf,nbf,3), source=0.0d0)
    do c = 1, 3
      call expand_antisym(h10p(:,c), h1ao(:,:,c), nbf)
      call expand_antisym(s10p(:,c), s1ao(:,:,c), nbf)
    end do

    allocate(dm(nbf,nbf), dm_b(nbf,nbf), source=0.0d0)
    call unpack_sym(dmat_a, dm, nbf)
    call unpack_sym(dmat_b, dm_b, nbf)

    allocate(vxa(3,nbf,nbf), vxb(3,nbf,nbf), source=0.0d0)
    if (is_dft) then
      mxAngMom = maxval(basis%am) + 2
      call dft_initialize(infos, basis, molGrid)
      block
        real(kind=dp), allocatable :: ca(:,:), cb(:,:)
        allocate(ca(nbf,nbf), cb(nbf,nbf))
        ca = mo_a(:,1:nbf)
        cb = mo_b(:,1:nbf)
        call giao_vxc(basis, molGrid, infos, ca, cb, .true., &
                      vxa, vxb, mxAngMom, nbf, infos%dft%grid_density_cutoff)
        deallocate(ca, cb)
      end block
    end if

    allocate(twoe(3,nbf,nbf), twoe2(3,nbf,nbf), vj(3,nbf,nbf), &
             vk(3,nbf,nbf), vkb(3,nbf,nbf), source=0.0d0)
    call giao_h10_twoe_matrix(basis, infos, dm+dm_b, vj,  twoe, twoe2) ! J[D_tot]
    call giao_h10_twoe_matrix(basis, infos, dm,      twoe, vk,  twoe2) ! K[D_a]
    call giao_h10_twoe_matrix(basis, infos, dm_b,    twoe, vkb, twoe2) ! K[D_b]
    do c = 1, 3
      do i = 1, nbf
        do j = 1, nbf
          h1ao_b(i,j,c) = h1ao(i,j,c) + vj(c,i,j) - scale_exch*vkb(c,i,j) + SX*vxb(c,i,j)
          h1ao(i,j,c)   = h1ao(i,j,c) + vj(c,i,j) - scale_exch*vk(c,i,j)  + SX*vxa(c,i,j)
        end do
      end do
    end do

    allocate(sig_u(3,3,nat), sig_c(3,3,nat), source=0.0d0)
    allocate(ca_sc(nbf,nmo), cb_sc(nbf,nmo), ea_sc(nmo), eb_sc(nmo), source=0.0d0)
    call semicanon_orbitals(fock_a, mo_a, nbf, nmo, nocc_a, ca_sc, ea_sc)
    call semicanon_orbitals(fock_b, mo_b, nbf, nmo, nocc_b, cb_sc, eb_sc)
    call giao_para_channel(infos, basis, ca_sc, ea_sc, nocc_a, nmo, nbf, nat, &
                           coords, h1ao,   s1ao, scale_exch, 1.0d0, sig_u, sig_c)
    call giao_para_channel(infos, basis, cb_sc, eb_sc, nocc_b, nmo, nbf, nat, &
                           coords, h1ao_b, s1ao, scale_exch, 1.0d0, sig_u, sig_c)
    sig_u = sig_u * a2ppm
    sig_c = sig_c * a2ppm

    ! ------------------------------------------------------------------
    ! Store: OQP::nmr_shielding_mrsf, shape (32, nat) per atom:
    !   1-9   dia tensor sigma(t,s), index (s-1)*3+t
    !   10-18 para uncoupled tensor (frozen-reference approx)
    !   19-27 para coupled tensor (frozen-reference approx)
    !   28-32 iso: dia, para_unc, para_cpl, total_unc, total_cpl
    ! ------------------------------------------------------------------
    call infos%dat%remove_records((/ character(len=80) :: OQP_nmr_shielding_mrsf /))
    call infos%dat%reserve_data(OQP_nmr_shielding_mrsf, TA_TYPE_REAL64, 32*nat, &
                                (/ 32, nat /), comment=OQP_nmr_shielding_mrsf_comment)
    call tagarray_get_data(infos%dat, OQP_nmr_shielding_mrsf, nmrout)
    do iat = 1, nat
      do s = 1, 3
        do t = 1, 3
          nmrout((s-1)*3+t,    iat) = sig_dia(t,s,iat)
          nmrout((s-1)*3+t+9,  iat) = sig_u(t,s,iat)
          nmrout((s-1)*3+t+18, iat) = sig_c(t,s,iat)
        end do
      end do
      nmrout(28,iat) = iso3(sig_dia(:,:,iat))
      nmrout(29,iat) = iso3(sig_u(:,:,iat))
      nmrout(30,iat) = iso3(sig_c(:,:,iat))
      nmrout(31,iat) = nmrout(28,iat) + nmrout(29,iat)
      nmrout(32,iat) = nmrout(28,iat) + nmrout(30,iat)
    end do

    ! ------------------------------------------------------------------
    ! Report (machine-parseable + human table)
    ! ------------------------------------------------------------------
    write(iw,'(/,A)') 'MRSF_NMR_BEGIN prototype state shielding (ppm)'
    write(iw,'(A,1X,I0)') 'MRSF_NMR_STATE', target_state
    write(iw,'(A,1X,I0)') 'MRSF_NMR_NATOM', nat
    write(iw,'(A,1X,F10.6)') 'MRSF_NMR_CX', scale_exch
    call tagarray_get_data(infos%dat, OQP_td_energies, td_energies, status)
    if (status == 0) then
      if (target_state >= 1 .and. target_state <= size(td_energies)) then
        write(iw,'(A,1X,ES24.16)') 'MRSF_NMR_OMEGA', td_energies(target_state)
      end if
    end if
    do iat = 1, nat
      do t = 1, 3
        do s = 1, 3
          write(iw,'(A,1X,I0,1X,I0,1X,I0,1X,ES24.16)') 'MRSF_NMR_DIA', &
            iat, t, s, sig_dia(t,s,iat)
          write(iw,'(A,1X,I0,1X,I0,1X,I0,1X,ES24.16)') 'MRSF_NMR_DIA_REF', &
            iat, t, s, sig_dia_ref(t,s,iat)
          write(iw,'(A,1X,I0,1X,I0,1X,I0,1X,ES24.16)') 'MRSF_NMR_PARA_UNC', &
            iat, t, s, sig_u(t,s,iat)
          write(iw,'(A,1X,I0,1X,I0,1X,I0,1X,ES24.16)') 'MRSF_NMR_PARA_CPL', &
            iat, t, s, sig_c(t,s,iat)
        end do
      end do
      write(iw,'(A,1X,I0,5(1X,ES24.16))') 'MRSF_NMR_ISO', iat, &
        iso3(sig_dia(:,:,iat)), iso3(sig_u(:,:,iat)), iso3(sig_c(:,:,iat)), &
        iso3(sig_dia(:,:,iat)) + iso3(sig_u(:,:,iat)), &
        iso3(sig_dia(:,:,iat)) + iso3(sig_c(:,:,iat))
    end do
    write(iw,'(A)') 'MRSF_NMR_END'

    write(iw,'(2/)')
    write(iw,'(4x,a)') '====================================================='
    write(iw,'(4x,a,i0,a)') 'MRSF-NMR shielding, state ', target_state, &
      '  (Gate 2 PROTOTYPE)'
    write(iw,'(4x,a)') '====================================================='
    write(iw,'(4x,a)') 'Diamagnetic term: EXACT state value (GIAO kernels x relaxed'
    write(iw,'(4x,a)') '  MRSF density DM+td_p; THEORY.tex master formula).'
    write(iw,'(4x,a)') 'Paramagnetic term: FROZEN-REFERENCE approximation (ROHF'
    write(iw,'(4x,a)') '  reference GIAO CPHF; true MRSF response lands in Gates 3-4).'
    write(iw,'(4x,a)') 'NOTE: dia/para split is origin-convention dependent; only the'
    write(iw,'(4x,a)') '  coupled total of the final method will be the observable.'
    write(iw,'(/4x,a,f8.4,a)') 'Isotropic shielding (ppm)   [exact-exchange c_x =', &
           scale_exch, ']'
    write(iw,'(4x,a)') '   Atom    Z   sigma_dia(state)  sigma_dia(ref)  Delta_dia'// &
           '   para_cpl(ref)   total(prototype)'
    do iat = 1, nat
      write(iw,'(4x,i6,f6.1,5f16.6)') iat, basis%atoms%zn(iat), &
        iso3(sig_dia(:,:,iat)), iso3(sig_dia_ref(:,:,iat)), &
        iso3(sig_dia(:,:,iat)) - iso3(sig_dia_ref(:,:,iat)), &
        iso3(sig_c(:,:,iat)), &
        iso3(sig_dia(:,:,iat)) + iso3(sig_c(:,:,iat))
    end do
    call flush(iw)
    close(iw)

    deallocate(dens_state, dens_ref, coords, zq, gdia0, corrpre, a01)
    deallocate(sig_dia, sig_dia_ref, h10p, s10p, h1ao, h1ao_b, s1ao)
    deallocate(dm, dm_b, vxa, vxb, twoe, twoe2, vj, vk, vkb)
    deallocate(sig_u, sig_c, ca_sc, cb_sc, ea_sc, eb_sc)

  contains

    !> GIAO diamagnetic tensor for a given (packed, symmetric) density:
    !> identical assembly to nmr_giao_shielding.F90 (a11part trace-corrected
    !> + a01gp), evaluated at the coordinate origin convention o0 = 0.
    subroutine dia_tensor(denp, sig)
      real(kind=dp), intent(in) :: denp(:)
      real(kind=dp), intent(out) :: sig(:,:,:)
      real(kind=dp) :: o0(3)
      integer :: ia, tt, ss
      o0 = 0.0d0
      gdia0 = 0.0d0; corrpre = 0.0d0; a01 = 0.0d0
      call nmr_dia_shielding(basis, denp, o0, coords, nat, gdia0, logtol=tol)
      call giao_a11part_corr(basis, denp, coords, nat, corrpre, logtol=tol)
      call giao_a01gp_contract(basis, denp, coords, nat, a01, logtol=tol)
      do ia = 1, nat
        trg0 = gdia0(1,1,ia)+gdia0(2,2,ia)+gdia0(3,3,ia)
        trc  = corrpre(1,1,ia)+corrpre(2,2,ia)+corrpre(3,3,ia)
        do tt = 1, 3
          do ss = 1, 3
            sig(tt,ss,ia) = ( SG*0.5d0*gdia0(ss,tt,ia) + SC*corrpre(tt,ss,ia) &
                   - merge(SG*0.5d0*trg0 + SC*trc, 0.0d0, tt==ss) &
                   + SA*a01(tt,ss,ia) ) * a2ppm
          end do
        end do
      end do
    end subroutine dia_tensor

  end subroutine nmr_mrsf_shielding

  pure function iso3(sig) result(v)
    real(kind=dp), intent(in) :: sig(:,:)
    real(kind=dp) :: v
    v = (sig(1,1) + sig(2,2) + sig(3,3))/3.0d0
  end function iso3

!> Expand packed lower-triangular antisymmetric matrix to full form.
  subroutine expand_antisym(packed, full, n)
    real(kind=dp), intent(in)  :: packed(:)
    real(kind=dp), intent(out) :: full(:,:)
    integer, intent(in) :: n
    integer :: p, q, idx
    full = 0.0d0
    do p = 1, n
      do q = 1, p
        idx = q + p*(p-1)/2
        full(p,q) =  packed(idx)
        full(q,p) = -packed(idx)
      end do
    end do
  end subroutine expand_antisym

  subroutine unpack_sym(packed, full, n)
    real(kind=dp), intent(in) :: packed(:)
    real(kind=dp), intent(out) :: full(:,:)
    integer, intent(in) :: n
    integer :: i, j, ij
    full = 0.0d0
    do i = 1, n
      do j = 1, i
        ij = j + i*(i-1)/2
        full(i,j) = packed(ij)
        full(j,i) = packed(ij)
      end do
    end do
  end subroutine unpack_sym

end module nmr_mrsf_shielding_mod
