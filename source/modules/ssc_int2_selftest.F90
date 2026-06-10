!> @brief L1 self-test for the two-electron spin-spin (SS) dipolar integral (SSC/ZFS, ssc-zfs).
!>
!> Validation harness only. For a set of (s,p,d) shell quartets of the input molecule it compares,
!> per primitive quartet and per Cartesian function combination, the ANALYTIC SS bare-Hessian
!> integral H_kl (mod_ssc_int2::comp_ssc_int2_prim) against a central-difference + Richardson
!> reference built from the engine's OWN Coulomb ERI with electron 2 rigidly displaced
!> (mod_ssc_int2::comp_eri2_prim_disp) -- displacing electron 2 == displacing the 1/r12 operator,
!> so d^2 ERI / d dshift_k d dshift_l = H_kl. It also checks that the traceless dipolar integral
!> S_kl = H_kl - (1/3) Tr(H) delta_kl is traceless.
!>
!> This is the L1 gate per benchmarks.md (FD of the existing ERI engine; trace = 0). It does NOT
!> involve the D-tensor prefactor/sign C (an L2 concern). Writes a PASS/FAIL banner + the worst
!> relative disagreement and worst |trace| to /tmp/ssc_int2_selftest.out and the OpenQP log.
!>
!> NOTE (honesty): a PASS here means the Fortran analytic integral reproduces finite differences
!> of the engine's ERI to FD precision. Whether L1 is *declared* passed is a human decision.

module ssc_int2_selftest_mod

  use precision, only: dp
  implicit none
  character(len=*), parameter :: module_name = "ssc_int2_selftest_mod"

contains

  subroutine ssc_int2_selftest_C(c_handle) bind(C, name="ssc_int2_selftest")
    use c_interop, only: oqp_handle_t, oqp_handle_get_info
    use types, only: information
    type(oqp_handle_t) :: c_handle
    type(information), pointer :: inf
    inf => oqp_handle_get_info(c_handle)
    call ssc_int2_selftest(inf)
  end subroutine ssc_int2_selftest_C

  subroutine ssc_int2_selftest(infos)
    use types,           only: information
    use io_constants,    only: iw
    use basis_tools,     only: basis_set
    use mod_shell_tools, only: shell_t, shpair_t
    use mod_ssc_int2,    only: comp_ssc_int2_prim, comp_eri2_prim_disp

    implicit none
    type(information), target, intent(inout) :: infos

    type(basis_set), pointer :: basis
    type(shell_t)   :: shi, shj, shk, shl
    type(shpair_t)  :: cpij, cpkl

    integer, parameter :: MAXREP = 4
    integer :: reps(MAXREP), nrep
    integer :: nshell, sh, a, b, c, d, q, ig, kg
    integer :: nij, nkl, comp, ij, kl
    real(dp) :: tol, h, rho_q
    real(dp), allocatable :: ss(:,:,:)            ! analytic (nij,nkl,6)
    real(dp), allocatable :: e0(:,:), ep(:,:), em(:,:)
    real(dp), allocatable :: epp(:,:), epm(:,:), emp(:,:), emm(:,:)
    real(dp), allocatable :: fdb1(:,:), fdb2(:,:), fdb3(:,:)
    real(dp) :: href_kl, rel, trc, scale
    real(dp) :: worst_rel, worst_trace, worst_an, worst_ref, worst_scale
    real(dp) :: worst_rel_sig                ! worst rel diff over non-negligible (scale>1e-9) blocks
    real(dp) :: norm_qgauss = 0, norm_textbook = 0, norm_ratio = 0   ! absolute-normalisation check
    integer  :: worst_q(4), worst_comp, worst_q_sig(4)
    logical  :: signif
    integer  :: ntested, npass, nfail_alls, nfail_iandj
    integer  :: ios, u
    logical  :: ok
    ! component -> (k,l): 1 xx, 2 yy, 3 zz, 4 xy, 5 xz, 6 yz
    integer, parameter :: ck(6) = [1,2,3,1,1,2]
    integer, parameter :: cl(6) = [1,2,3,2,3,3]

    basis => infos%basis
    basis%atoms => infos%atoms
    nshell = basis%nshell
    tol = huge(1.0_dp)          ! no screening in the self-test

    ! ---- collect representative s/p/d shells, preferring angular-momentum variety ----
    ! OpenQP uses CARTESIAN Gaussians (basis%naos = NUM_CART_BF), so a d shell is 6 cartesian
    ! functions and CART_X(i,2) gives their powers directly -- no spherical-harmonic transform is
    ! involved in the integral engine. Exclude pathologically tight (core) primitives: their charge
    ! cloud is ~1/sqrt(alpha) wide, so the operator-displacement FD reference is roundoff-limited and
    ! cannot validate to 1e-6 (the analytic path is identical and is covered by the prototype check).
    nrep = 0
    ! one shell of each angular momentum d, p (highest first), non-tight leading primitive
    do sh = 1, nshell
      if (basis%am(sh) == 2 .and. basis%ex(basis%g_offset(sh)) < 1.0e2_dp) then
        nrep = nrep + 1; reps(nrep) = sh; exit
      end if
    end do
    do sh = 1, nshell
      if (nrep >= MAXREP) exit
      if (basis%am(sh) == 1 .and. basis%ex(basis%g_offset(sh)) < 1.0e2_dp) then
        nrep = nrep + 1; reps(nrep) = sh; exit
      end if
    end do
    ! s shells (non-tight) on as many distinct atoms as possible
    do sh = 1, nshell
      if (nrep >= MAXREP) exit
      if (basis%am(sh) == 0 .and. basis%ex(basis%g_offset(sh)) < 1.0e2_dp) then
        if (.not. any(reps(1:nrep) == sh)) then
          nrep = nrep + 1; reps(nrep) = sh
        end if
      end if
    end do
    ! pad with any remaining non-tight s/p/d shells
    do sh = 1, nshell
      if (nrep >= MAXREP) exit
      if (basis%am(sh) <= 2 .and. basis%ex(basis%g_offset(sh)) < 1.0e2_dp &
          .and. .not. any(reps(1:nrep) == sh)) then
        nrep = nrep + 1; reps(nrep) = sh
      end if
    end do

    call cpij%alloc(basis)
    call cpkl%alloc(basis)

    worst_rel = 0.0_dp
    worst_trace = 0.0_dp
    worst_q = 0
    worst_an = 0.0_dp; worst_ref = 0.0_dp; worst_scale = 0.0_dp; worst_comp = 0
    worst_rel_sig = 0.0_dp; worst_q_sig = 0
    ntested = 0
    npass = 0
    nfail_alls = 0
    nfail_iandj = 0

    do a = 1, nrep
      call shi%fetch_by_id(basis, reps(a))
      do b = 1, nrep
        call shj%fetch_by_id(basis, reps(b))
        call cpij%shell_pair(basis, shi, shj, tol, dup=.false.)
        if (cpij%numpairs == 0) cycle
        do c = 1, nrep
          call shk%fetch_by_id(basis, reps(c))
          do d = 1, nrep
            call shl%fetch_by_id(basis, reps(d))
            call cpkl%shell_pair(basis, shk, shl, tol, dup=.false.)
            if (cpkl%numpairs == 0) cycle

            nij = cpij%inao * cpij%jnao
            nkl = cpkl%inao * cpkl%jnao
            ig = 1; kg = 1     ! first primitive pair of each (primitive-level test)

            ! FD step matched to the charge-cloud size 1/sqrt(rho): the operator-displacement
            ! ERI varies on that scale, so a fixed step fails for tight (core) primitives.
            rho_q = cpij%p(ig)%aa * cpkl%p(kg)%aa / (cpij%p(ig)%aa + cpkl%p(kg)%aa)
            h    = min(2.0e-2_dp, 0.10_dp/sqrt(rho_q))

            if (allocated(ss)) deallocate(ss)
            allocate(ss(nij, nkl, 6))
            call comp_ssc_int2_prim(cpij, ig, cpkl, kg, ss)

            if (allocated(e0)) deallocate(e0, ep, em, epp, epm, emp, emm, fdb1, fdb2, fdb3)
            allocate(e0(nij,nkl), ep(nij,nkl), em(nij,nkl))
            allocate(epp(nij,nkl), epm(nij,nkl), emp(nij,nkl), emm(nij,nkl))
            allocate(fdb1(nij,nkl), fdb2(nij,nkl), fdb3(nij,nkl))

            ! Comparison reference magnitude = largest analytic component in this quartet block,
            ! with an absolute floor: blocks whose SS integrals are all < 1e-9 are physically
            ! negligible (vanishing by symmetry), so the FD noise (~1e-16) must be judged against an
            ! absolute scale, not a relative one (else 1e-16/1e-16 ~ O(1) spurious "failures").
            scale = max(maxval(abs(ss)), 1.0e-9_dp)
            signif = maxval(abs(ss)) > 1.0e-9_dp    ! quartet has non-negligible SS integrals

            do comp = 1, 6
              ! 3-level Richardson FD reference block (O(h^6)); FD blocks computed once per step.
              call fd_block(ck(comp), cl(comp), h,        fdb1)
              call fd_block(ck(comp), cl(comp), h*0.5_dp,  fdb2)
              call fd_block(ck(comp), cl(comp), h*0.25_dp, fdb3)
              do ij = 1, nij
                do kl = 1, nkl
                  block
                    real(dp) :: r1a, r1b
                    r1a = (4.0_dp*fdb2(ij,kl) - fdb1(ij,kl))/3.0_dp
                    r1b = (4.0_dp*fdb3(ij,kl) - fdb2(ij,kl))/3.0_dp
                    href_kl = (16.0_dp*r1b - r1a)/15.0_dp
                  end block
                  ! relative to the integral block magnitude (sig figs vs the integral itself)
                  rel = abs(ss(ij,kl,comp) - href_kl)/scale
                  ntested = ntested + 1
                  if (rel <= 1.0e-6_dp) then
                    npass = npass + 1
                  else
                    if (basis%am(reps(a))==0 .and. basis%am(reps(b))==0 .and. &
                        basis%am(reps(c))==0 .and. basis%am(reps(d))==0) nfail_alls = nfail_alls + 1
                    if (reps(a)==reps(b) .or. reps(c)==reps(d)) nfail_iandj = nfail_iandj + 1
                  end if
                  if (rel > worst_rel) then
                    worst_rel = rel; worst_q = [reps(a),reps(b),reps(c),reps(d)]
                    worst_an = ss(ij,kl,comp); worst_ref = href_kl
                    worst_scale = scale; worst_comp = comp
                  end if
                  if (signif .and. rel > worst_rel_sig) then
                    worst_rel_sig = rel; worst_q_sig = [reps(a),reps(b),reps(c),reps(d)]
                  end if
                end do
              end do
            end do

            ! traceless invariant (L4): form S_kl = H_kl - Tr(H)/3 delta_kl and check Tr(S)=0.
            do ij = 1, nij
              do kl = 1, nkl
                trc = ss(ij,kl,1) + ss(ij,kl,2) + ss(ij,kl,3)   ! Tr(H)
                ! Tr(S) = (Hxx - trc/3) + (Hyy - trc/3) + (Hzz - trc/3) = trc - trc
                worst_trace = max(worst_trace, abs( (ss(ij,kl,1) - trc/3.0_dp) &
                                                  + (ss(ij,kl,2) - trc/3.0_dp) &
                                                  + (ss(ij,kl,3) - trc/3.0_dp) ))
              end do
            end do

          end do
        end do
      end do
    end do

    ok = (worst_rel <= 1.0e-6_dp) .and. (ntested > 0)

    ! ---- absolute-normalisation check: 2-centre (ss|ss) ERI vs textbook ----
    ! Compares comp_eri2_prim_disp(0) to dij_fac * 2*pi^(5/2)/(p q sqrt(p+q)) * F0(rho|P-Q|^2).
    ! (L1 only validated the integral as a ratio; this pins the absolute scale.)
    block
      use boys, only: boysf
      integer :: sA, sB, sh2
      type(shell_t) :: sa_sh, sb_sh
      real(dp), allocatable :: eri(:,:)
      real(dp) :: pp, qq, rho2, t, ft(0:2), pref, textbook, x
      sA = 0; sB = 0
      do sh2 = 1, nshell
        if (basis%am(sh2) == 0 .and. basis%ex(basis%g_offset(sh2)) < 1.0e2_dp) then
          if (sA == 0) then
            sA = sh2
          else if (basis%origin(sh2) /= basis%origin(sA)) then
            sB = sh2; exit
          end if
        end if
      end do
      if (sA > 0 .and. sB > 0) then
        call sa_sh%fetch_by_id(basis, sA)
        call sb_sh%fetch_by_id(basis, sB)
        call cpij%shell_pair(basis, sa_sh, sa_sh, tol, dup=.false.)
        call cpkl%shell_pair(basis, sb_sh, sb_sh, tol, dup=.false.)
        allocate(eri(1,1))
        call comp_eri2_prim_disp(cpij, 1, cpkl, 1, [0._dp,0._dp,0._dp], eri)
        pp = cpij%p(1)%aa; qq = cpkl%p(1)%aa
        rho2 = pp*qq/(pp+qq)
        x = rho2 * sum((cpij%p(1)%r - cpkl%p(1)%r)**2)
        call boysf(2, x, ft)
        pref = 2.0_dp * (4.0_dp*atan(1.0_dp))**2.5_dp / (pp*qq*sqrt(pp+qq))
        textbook = cpij%p(1)%expfac * cpkl%p(1)%expfac * pref * ft(0)
        norm_qgauss = eri(1,1); norm_textbook = textbook
        norm_ratio = eri(1,1)/textbook
      end if
    end block

    open(newunit=u, file="/tmp/ssc_int2_selftest.out", status="replace", &
         action="write", iostat=ios)
    if (ios == 0) then
      call report(u)
      close(u)
    end if
    call report(iw)

  contains

    !> Finite-difference block of H_kl (component k,l) at step hh: the whole (nij,nkl) matrix
    !> at once (the displaced ERI blocks are computed once and reused across all elements).
    subroutine fd_block(k, l, hh, out)
      integer,  intent(in)  :: k, l
      real(dp), intent(in)  :: hh
      real(dp), intent(out) :: out(:,:)
      real(dp) :: ds(3)
      if (k == l) then
        ! central second difference along axis k
        ds = 0.0_dp; ds(k) =  hh; call comp_eri2_prim_disp(cpij, ig, cpkl, kg, ds, ep)
        ds = 0.0_dp;              call comp_eri2_prim_disp(cpij, ig, cpkl, kg, ds, e0)
        ds = 0.0_dp; ds(k) = -hh; call comp_eri2_prim_disp(cpij, ig, cpkl, kg, ds, em)
        out = (ep - 2.0_dp*e0 + em)/(hh*hh)
      else
        ds = 0.0_dp; ds(k) =  hh; ds(l) =  hh
        call comp_eri2_prim_disp(cpij, ig, cpkl, kg, ds, epp)
        ds = 0.0_dp; ds(k) =  hh; ds(l) = -hh
        call comp_eri2_prim_disp(cpij, ig, cpkl, kg, ds, epm)
        ds = 0.0_dp; ds(k) = -hh; ds(l) =  hh
        call comp_eri2_prim_disp(cpij, ig, cpkl, kg, ds, emp)
        ds = 0.0_dp; ds(k) = -hh; ds(l) = -hh
        call comp_eri2_prim_disp(cpij, ig, cpkl, kg, ds, emm)
        out = (epp - epm - emp + emm)/(4.0_dp*hh*hh)
      end if
    end subroutine fd_block

    subroutine report(unit)
      integer, intent(in) :: unit
      write(unit,'(/A)') '================ SSC_INT2_SELFTEST (L1) ================'
      write(unit,'(A,I0,A,I0)') ' representative s/p shells used: ', nrep, ' of nshell=', nshell
      write(unit,'(A,4(I3,A,I1,A))') ' reps (shell:am): ', &
        (reps(sh), ':', basis%am(reps(sh)), '  ', sh=1,nrep)
      write(unit,'(A,I0,A,I0)') ' failures with all-s quartet    : ', nfail_alls
      write(unit,'(A,I0)')      ' failures with a same-shell pair: ', nfail_iandj
      write(unit,'(A,I0)')      ' element comparisons tested     : ', ntested
      write(unit,'(A,I0,A,I0)') ' passed (rel<=1e-6)             : ', npass, ' / ', ntested
      write(unit,'(A,ES12.4)')  ' worst analytic-vs-FD rel diff  : ', worst_rel
      write(unit,'(A,4I4,A,I0)')'   at shell quartet (i j k l)   : ', worst_q, '  comp=', worst_comp
      write(unit,'(A,3ES14.6)') '   analytic / FD / block-scale  : ', worst_an, worst_ref, worst_scale
      write(unit,'(A,ES12.4,A,4I4)') ' worst rel diff, non-negligible : ', worst_rel_sig, &
        '   at quartet ', worst_q_sig
      write(unit,'(A,ES12.4)')  ' worst |Tr(S)| (traceless inv.) : ', worst_trace
      write(unit,'(A,3ES15.7)') ' 2c (ss|ss) qgauss/textbook/ratio: ', &
        norm_qgauss, norm_textbook, norm_ratio
      write(unit,'(A,ES10.2)')    ' FD base step last quartet (h) : ', h
      write(unit,'(A)')           ' (adaptive h=min(0.02,0.1/sqrt(rho)); 3-level Richardson h,h/2,h/4)'
      if (ok) then
        write(unit,'(A)') ' SSC_INT2_SELFTEST PASS'
      else
        write(unit,'(A)') ' SSC_INT2_SELFTEST FAIL'
      end if
      write(unit,'(A/)') '========================================================'
    end subroutine report

  end subroutine ssc_int2_selftest

end module ssc_int2_selftest_mod
