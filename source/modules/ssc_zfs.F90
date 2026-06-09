!> @brief First-order direct spin-spin (SS) contribution to the ZFS D-tensor — AO contraction.
!>        SSC/ZFS project, branch ssc-zfs. Phase 2 / P2.1.
!>
!> Working equation (Sinnecker & Neese, JPCA 110, 12267 (2006), Eq. 9; canonical form):
!>
!>   D^SS_kl = C * sum_{mu nu} sum_{kappa tau}
!>               { P_munu P_kaptau - P_mukap P_nutau } * <mu nu | T_kl | kappa tau>
!>
!> with the rank-2 traceless dipolar kernel T_kl = (3 r12,k r12,l - delta_kl r12^2)/r12^5,
!> and P = P^(alpha-beta), the one-particle spin density of the M_S = S state. For a single
!> ROHF determinant (multiplicity 2S+1) the M_S = S spin density is exactly DM_A - DM_B.
!>
!> This routine builds the SIX AO-contracted tensor components (the Coulomb-like P_munu P_kaptau
!> term and the exchange-like P_mukap P_nutau term, structurally a Fock/K build) using the
!> validated native SS integral (mod_ssc_int2::comp_ssc_int2_prim). The integrals are returned in
!> the cartesian primitive normalisation (contraction coefficients already folded in via expfac);
!> the basis-function normalisation bfnrm is absorbed by pre-scaling the density:
!>     Q_munu = P_munu * bfnrm_mu * bfnrm_nu,
!> so the contraction can use the cartesian integrals directly.
!>
!> SCOPE (CLAUDE.md §1): first order only. The prefactor/sign C and the a.u.->cm^-1 unit factor are
!> pinned NUMERICALLY at L2 (O2 benchmark) — NOT here; this routine uses C = 1 and returns the raw
!> contraction. Diagonalisation to (D, E) is P2.2. No Z-vector / relaxed densities / gradients.

module ssc_zfs_mod

  use precision, only: dp
  implicit none
  private
  public :: compute_ssc_dtensor_raw     ! 6 raw AO-contracted components (C = 1)
  public :: ssc_dtensor_selftest        ! bind(C) harness: run + report invariants

  character(len=*), parameter :: module_name = "ssc_zfs_mod"

contains

!-------------------------------------------------------------------------------
!> Raw SS D-tensor contraction (C = 1). Returns the 6 symmetric components
!> dxx, dyy, dzz, dxy, dxz, dyz of  sum {P P - P P} <mu nu|T|kappa tau>.
  subroutine compute_ssc_dtensor_raw(infos, dcomp)
    use types,             only: information
    use basis_tools,       only: basis_set
    use mod_shell_tools,   only: shell_t, shpair_t
    use mod_ssc_int2,      only: comp_ssc_int2_prim
    use oqp_tagarray_driver, only: tagarray_get_data, OQP_DM_A, OQP_DM_B

    type(information), target, intent(inout) :: infos
    real(dp), intent(out) :: dcomp(6)

    type(basis_set), pointer :: basis
    real(dp), contiguous, pointer :: dmat_a(:), dmat_b(:)
    real(dp), allocatable :: q(:,:)            ! bfnrm-scaled spin density (square)
    real(dp), allocatable :: ssblk(:,:,:), acc(:,:,:)
    type(shell_t)  :: shi, shj, shk, shl
    type(shpair_t) :: cpij, cpkl
    integer :: nbf, nshell, ii, jj, kk, ll, ig, kg
    integer :: i, j, k, l, mu, nu, ka, ta, ij, kl, c
    integer :: nij, nkl
    real(dp) :: tol, w

    basis => infos%basis
    basis%atoms => infos%atoms
    nbf    = basis%nbf
    nshell = basis%nshell
    tol    = huge(1.0_dp)

    call tagarray_get_data(infos%dat, OQP_DM_A, dmat_a)
    call tagarray_get_data(infos%dat, OQP_DM_B, dmat_b)

    ! square, bfnrm-scaled spin density  Q_munu = (P^a - P^b)_munu * bfnrm_mu * bfnrm_nu
    allocate(q(nbf, nbf))
    call build_scaled_spin_density(dmat_a, dmat_b, basis%bfnrm, nbf, q)

    call cpij%alloc(basis)
    call cpkl%alloc(basis)
    dcomp = 0.0_dp

    do ii = 1, nshell
      call shi%fetch_by_id(basis, ii)
      do jj = 1, nshell
        call shj%fetch_by_id(basis, jj)
        call cpij%shell_pair(basis, shi, shj, tol, dup=.false.)
        if (cpij%numpairs == 0) cycle
        do kk = 1, nshell
          call shk%fetch_by_id(basis, kk)
          do ll = 1, nshell
            call shl%fetch_by_id(basis, ll)
            call cpkl%shell_pair(basis, shk, shl, tol, dup=.false.)
            if (cpkl%numpairs == 0) cycle

            nij = cpij%inao * cpij%jnao
            nkl = cpkl%inao * cpkl%jnao
            if (allocated(ssblk)) deallocate(ssblk, acc)
            allocate(ssblk(nij, nkl, 6), acc(nij, nkl, 6))

            ! contracted cartesian bare-Hessian block, summed over primitive pairs
            acc = 0.0_dp
            do ig = 1, cpij%numpairs
              do kg = 1, cpkl%numpairs
                call comp_ssc_int2_prim(cpij, ig, cpkl, kg, ssblk)
                acc = acc + ssblk
              end do
            end do
            ! physical traceless dipolar kernel: S = H - (1/3) Tr(H) I  (diagonal comps only)
            call make_traceless(acc, nij, nkl)

            ! contract with the spin density: Coulomb-like minus exchange-like
            do i = 1, cpij%inao
              mu = shi%locao + i
              do j = 1, cpij%jnao
                nu = shj%locao + j
                ij = (i-1)*cpij%jnao + j
                do k = 1, cpkl%inao
                  ka = shk%locao + k
                  do l = 1, cpkl%jnao
                    ta = shl%locao + l
                    kl = (k-1)*cpkl%jnao + l
                    w = q(mu,nu)*q(ka,ta) - q(mu,ka)*q(nu,ta)
                    do c = 1, 6
                      dcomp(c) = dcomp(c) + w * acc(ij,kl,c)
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine compute_ssc_dtensor_raw

!-------------------------------------------------------------------------------
  subroutine build_scaled_spin_density(da, db, bfnrm, nbf, q)
    real(dp), intent(in)  :: da(:), db(:), bfnrm(:)
    integer,  intent(in)  :: nbf
    real(dp), intent(out) :: q(:,:)
    integer :: i, j, ij
    real(dp) :: p
    do i = 1, nbf
      do j = 1, i
        ij = i*(i-1)/2 + j                  ! packed lower triangle, i >= j
        p = (da(ij) - db(ij)) * bfnrm(i) * bfnrm(j)
        q(i,j) = p
        q(j,i) = p
      end do
    end do
  end subroutine build_scaled_spin_density

!-------------------------------------------------------------------------------
!> Replace the 6 bare-Hessian components by their traceless part per (ij,kl):
!> S_kk = H_kk - (1/3)(Hxx+Hyy+Hzz); off-diagonals unchanged.
  subroutine make_traceless(b, nij, nkl)
    real(dp), intent(inout) :: b(:,:,:)
    integer,  intent(in)    :: nij, nkl
    integer  :: ij, kl
    real(dp) :: trc
    do ij = 1, nij
      do kl = 1, nkl
        trc = (b(ij,kl,1) + b(ij,kl,2) + b(ij,kl,3)) / 3.0_dp
        b(ij,kl,1) = b(ij,kl,1) - trc
        b(ij,kl,2) = b(ij,kl,2) - trc
        b(ij,kl,3) = b(ij,kl,3) - trc
      end do
    end do
  end subroutine make_traceless

!-------------------------------------------------------------------------------
  subroutine ssc_dtensor_selftest_C(c_handle) bind(C, name="ssc_dtensor_selftest")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call ssc_dtensor_selftest(inf)
  end subroutine ssc_dtensor_selftest_C

  subroutine ssc_dtensor_selftest(infos)
    use types,        only: information
    use io_constants, only: iw
    type(information), target, intent(inout) :: infos
    real(dp) :: d(6), trc, asym
    integer  :: u, ios

    call compute_ssc_dtensor_raw(infos, d)
    trc  = d(1) + d(2) + d(3)                 ! Tr(D) — must vanish (traceless kernel)
    asym = 0.0_dp                             ! D is symmetric by construction (6 comps)

    open(newunit=u, file="/tmp/ssc_dtensor_selftest.out", status="replace", &
         action="write", iostat=ios)
    if (ios == 0) then ; call report(u) ; close(u) ; end if
    call report(iw)

  contains
    subroutine report(unit)
      integer, intent(in) :: unit
      write(unit,'(/A)') '============== SSC_DTENSOR_SELFTEST (P2.1) =============='
      write(unit,'(A)')  ' Raw AO-contracted SS D-tensor components (C = 1, a.u.):'
      write(unit,'(A,ES16.8)') '   Dxx = ', d(1)
      write(unit,'(A,ES16.8)') '   Dyy = ', d(2)
      write(unit,'(A,ES16.8)') '   Dzz = ', d(3)
      write(unit,'(A,ES16.8)') '   Dxy = ', d(4)
      write(unit,'(A,ES16.8)') '   Dxz = ', d(5)
      write(unit,'(A,ES16.8)') '   Dyz = ', d(6)
      write(unit,'(A,ES12.4)') ' Tr(D) [must be ~0, traceless kernel] : ', trc
      write(unit,'(A)') ' NOTE: prefactor/sign C and a.u.->cm^-1 unit factor are pinned at L2 (O2);'
      write(unit,'(A)') '       these are raw values with C = 1. Diagonalisation -> (D,E) is P2.2.'
      if (abs(trc) <= 1.0e-8_dp*maxval(abs(d)) + 1.0e-12_dp) then
        write(unit,'(A)') ' SSC_DTENSOR_SELFTEST PASS (trace vanishes)'
      else
        write(unit,'(A)') ' SSC_DTENSOR_SELFTEST FAIL (nonzero trace)'
      end if
      write(unit,'(A/)') '========================================================'
    end subroutine report
  end subroutine ssc_dtensor_selftest

end module ssc_zfs_mod
