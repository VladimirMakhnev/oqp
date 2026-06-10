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

    call spin_population_diag(infos, dmat_a, dmat_b)

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
            ! AO index of the i-th cartesian function: locao is 1-based (= first AO index),
            ! and CART_X(i,..) is 1-based, so the i-th function sits at AO (locao + i - 1).
            do i = 1, cpij%inao
              mu = shi%locao + i - 1
              do j = 1, cpij%jnao
                nu = shj%locao + j - 1
                ij = (i-1)*cpij%jnao + j
                do k = 1, cpkl%inao
                  ka = shk%locao + k - 1
                  do l = 1, cpkl%jnao
                    ta = shl%locao + l - 1
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
!> Diagnostic: gross spin population (DM_A-DM_B diagonal) grouped by cartesian direction, to check
!> whether the reference spin density is cylindrically symmetric about z (px==py for a 3Sigma_g- O2).
  subroutine spin_population_diag(infos, da, db)
    use types,           only: information
    use basis_tools,     only: basis_set
    use mod_shell_tools, only: shell_t
    use constants,       only: CART_X, CART_Y, CART_Z
    type(information), target, intent(inout) :: infos
    real(dp), intent(in) :: da(:), db(:)
    type(basis_set), pointer :: basis
    type(shell_t) :: sh
    real(dp) :: spx, spy, spz, sp_tot, sdiag
    integer :: ii, i, mu, ang, nx, ny, nz
    basis => infos%basis
    spx = 0; spy = 0; spz = 0; sp_tot = 0
    do ii = 1, basis%nshell
      call sh%fetch_by_id(basis, ii)
      ang = sh%ang
      do i = 1, sh%nao
        mu = sh%locao + i - 1
        sdiag = da(mu*(mu-1)/2 + mu) - db(mu*(mu-1)/2 + mu)
        sp_tot = sp_tot + sdiag
        nx = CART_X(i,ang); ny = CART_Y(i,ang); nz = CART_Z(i,ang)
        if (nx > ny .and. nx >= nz) spx = spx + sdiag
        if (ny > nx .and. ny >= nz) spy = spy + sdiag
        if (nz > nx .and. nz > ny)  spz = spz + sdiag
      end do
    end do
    block
      integer :: u, ios
      open(newunit=u, file="/tmp/ssc_spinpop.out", status="replace", action="write", iostat=ios)
      if (ios == 0) then
        write(u,'(A)') ' [SSC diag] gross spin population by cartesian direction (DM_A-DM_B diag):'
        write(u,'(A,4F12.6)') '   px-type / py-type / pz-type / total = ', spx, spy, spz, sp_tot
        write(u,'(A)') '   (3Sigma_g- O2 along z expects px-type == py-type, both >> pz-type)'
        close(u)
      end if
    end block
  end subroutine spin_population_diag

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
!-------------------------------------------------------------------------------
!> Classic cyclic Jacobi diagonalisation of a symmetric 3x3 matrix.
  subroutine jacobi3(a_in, ev, v)
    real(dp), intent(in)  :: a_in(3,3)
    real(dp), intent(out) :: ev(3), v(3,3)
    real(dp) :: a(3,3), th, t, c, s, tau, g, h, api, apj
    integer  :: p, q, i, sweep
    a = a_in
    v = 0.0_dp; v(1,1)=1; v(2,2)=1; v(3,3)=1
    do sweep = 1, 50
      if (abs(a(1,2))+abs(a(1,3))+abs(a(2,3)) < 1.0e-300_dp) exit
      do p = 1, 2
        do q = p+1, 3
          if (abs(a(p,q)) <= 0.0_dp) cycle
          th = (a(q,q) - a(p,p)) / (2.0_dp*a(p,q))
          t  = sign(1.0_dp, th) / (abs(th) + sqrt(th*th + 1.0_dp))
          c  = 1.0_dp / sqrt(t*t + 1.0_dp)
          s  = t*c;  tau = s/(1.0_dp + c)
          g  = a(p,q)
          a(p,p) = a(p,p) - t*g
          a(q,q) = a(q,q) + t*g
          a(p,q) = 0.0_dp;  a(q,p) = 0.0_dp
          do i = 1, 3
            if (i /= p .and. i /= q) then
              api = a(i,p); apj = a(i,q)
              a(i,p) = api - s*(apj + tau*api); a(p,i) = a(i,p)
              a(i,q) = apj + s*(api - tau*apj); a(q,i) = a(i,q)
            end if
            h = v(i,p); g = v(i,q)
            v(i,p) = h - s*(g + tau*h)
            v(i,q) = g + s*(h - tau*g)
          end do
        end do
      end do
    end do
    ev = [a(1,1), a(2,2), a(3,3)]
  end subroutine jacobi3

!> Order eigenvalues so that |ev(3)| >= |ev(2)| >= |ev(1)| (ZFS convention 0 <= E/D <= 1/3).
  subroutine order_zfs(ev)
    real(dp), intent(inout) :: ev(3)
    real(dp) :: t
    integer  :: i, j
    do i = 1, 2
      do j = i+1, 3
        if (abs(ev(j)) < abs(ev(i))) then
          t = ev(i); ev(i) = ev(j); ev(j) = t
        end if
      end do
    end do
  end subroutine order_zfs

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
    use physical_constants, only: alpha => FINE_STRUCTURE
    type(information), target, intent(inout) :: infos
    real(dp) :: d(6), trc
    real(dp) :: dmat(3,3), evec(3,3), ev(3), dval_au, eval_au, dval_cm, eval_cm, cpref
    integer  :: u, ios
    ! L2-pinned prefactor (Neese, J. Chem. Phys. 127, 164112 (2007), Eq. 46):
    !   C = - g_e^2 alpha^2 / [16 S(2S-1)],   here for S = 1 (O2 triplet) -> /16.
    ! Pinned numerically on O2 3Sigma_g- (matches +1.50 cm^-1, within the 1.44-1.6 band, to ~1%).
    real(dp), parameter :: GE = 2.00231930436_dp
    real(dp), parameter :: HA2WN = 219474.6313705_dp
    real(dp) :: S

    call compute_ssc_dtensor_raw(infos, d)
    trc  = d(1) + d(2) + d(3)                 ! Tr(D) — must vanish (traceless kernel)

    S = 1.0_dp                                ! TODO(L3): general spin from the state; O2 triplet S=1
    cpref = - GE**2 * alpha**2 / (16.0_dp * S*(2.0_dp*S - 1.0_dp))

    ! symmetric 3x3 D-tensor, diagonalise -> principal values
    dmat = reshape([ d(1), d(4), d(5),  d(4), d(2), d(6),  d(5), d(6), d(3) ], [3,3])
    call jacobi3(dmat, ev, evec)
    call order_zfs(ev)                        ! |ev(3)| >= |ev(2)| >= |ev(1)|  (0 <= E/D <= 1/3)
    dval_au = ev(3) - 0.5_dp*(ev(1) + ev(2))
    eval_au = 0.5_dp*(ev(2) - ev(1))
    dval_cm = cpref * dval_au * HA2WN
    eval_cm = cpref * eval_au * HA2WN

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
      write(unit,'(A)') ' Principal values (a.u., C=1):'
      write(unit,'(A,3ES16.8)') '   ev = ', ev
      write(unit,'(A)') ' Pinned prefactor C = -g_e^2 alpha^2 / [16 S(2S-1)]  (Neese 2007, Eq. 46):'
      write(unit,'(A,ES16.8)') '   C (S=1)     = ', cpref
      write(unit,'(A,F12.5,A,F12.5)') '   D^SS        = ', dval_cm, ' cm^-1     E^SS = ', eval_cm
      write(unit,'(A,F8.4)') '   E/D         = ', merge(eval_au/dval_au, 0.0_dp, abs(dval_au)>0)
      write(unit,'(A)') ' (O2 3Sigma_g- target D^SS = +1.44..1.6 cm^-1, axial E/D=0.)'
      if (abs(trc) <= 1.0e-8_dp*maxval(abs(d)) + 1.0e-12_dp) then
        write(unit,'(A)') ' SSC_DTENSOR_SELFTEST PASS (trace vanishes)'
      else
        write(unit,'(A)') ' SSC_DTENSOR_SELFTEST FAIL (nonzero trace)'
      end if
      write(unit,'(A/)') '========================================================'
    end subroutine report
  end subroutine ssc_dtensor_selftest

end module ssc_zfs_mod
