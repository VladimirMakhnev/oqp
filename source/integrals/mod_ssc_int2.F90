!> @brief Two-electron spin-spin (SS) dipolar integral primitive for the first-order
!>        direct dipolar contribution to the ZFS D-tensor (SSC/ZFS project, branch ssc-zfs).
!>
!> The SS dipolar operator is the rank-2 traceless kernel
!>
!>     T_kl(r12) = (3 r12,k r12,l - delta_kl r12^2) / r12^5 ,   r12 = r1 - r2,
!>
!> which is the *traceless part* of the Hessian of 1/r12. Since d/dr12,k = d/dr1,k and
!> = -d/dr2,k, two integrations by parts give the working identity used here:
!>
!>     H_kl = integral of  d_k d_l (1/r12)  contracted with the two charge clouds
!>          = -<  d_k(mu nu) | 1/r12 | d_l(kappa tau) > ,
!>
!> i.e. one first derivative of the electron-1 charge distribution (direction k) and one of
!> the electron-2 distribution (direction l). This needs only +1 angular-momentum padding on
!> EACH electron, computed with the standard Rys/VRR/HRR ERI machinery (Path A). The physical
!> traceless dipolar integral is then S_kl = H_kl - (1/3) Tr(H) delta_kl.
!>
!> This module is deliberately self-contained (its own padded Rys builder) so it does not
!> perturb the tested SOC 2e path (comp_soc_int2_prim / QGaussRys2e in mod_1e_primitives).
!> The closed form for the (ss|ss) case was validated independently in
!> tests/ssc_prototype_ssss.py before this Fortran was written.
!>
!> NB: this is L1-stage integral code. The prefactor/sign C of the D-tensor is an L2 concern
!> (pinned numerically on O2); it is NOT applied here. These routines return the bare 2e
!> integral H (and its traceless part), in the engine's primitive normalisation.

module mod_ssc_int2

  use ISO_FORTRAN_ENV, only: real64
  use mod_shell_tools, only: shpair_t
  use rys,             only: rys_root_t
  use constants,       only: CART_X, CART_Y, CART_Z, MAX_ANG => BAS_MXANG

  implicit none

  private
  public :: comp_ssc_int2_prim     ! analytic 6-component SS bare-Hessian integral block
  public :: comp_eri2_prim_disp    ! plain ERI block with electron-2 rigidly displaced (FD ref)
  public :: ssc_nroots             ! Rys roots needed for an SS quartet

  real(real64), parameter :: PI252 = 34.986836655250_real64   ! 2*pi^(5/2)

contains

!-------------------------------------------------------------------------------
!> Number of Rys roots for the SS integral over a shell quartet.
!> The two derivatives raise the effective total angular momentum by up to 2.
  pure integer function ssc_nroots(iang, jang, kang, lang) result(nr)
    integer, intent(in) :: iang, jang, kang, lang
    nr = (iang + jang + kang + lang + 2)/2 + 1
  end function ssc_nroots

!-------------------------------------------------------------------------------
!> Padded Rys 2e 1D-integral table for one primitive quartet, with electron 2 rigidly
!> displaced by `dshift` (used both for the analytic build, dshift=0, and the FD reference).
!> Builds gfull(nj,ni,nl,nk,xyz,t) for nj=0..jang+1, ni=0..iang+1, nl=0..lang+1, nk=0..kang+1.
!> Faithful copy of QGaussRys2e (mod_1e_primitives) with: (i) electron-2 padded +1 as well as
!> electron 1; (ii) no GAMESS VRR truncation (full rectangle built); (iii) operator displacement
!> via pq -> (P-Q) - dshift. Dynamically dimensioned to avoid any index-range surprises.
  subroutine qgauss_ss(ryscomp, cpij, idij, cpkl, idkl, dshift, gfull)
    type(rys_root_t), intent(inout) :: ryscomp
    type(shpair_t),   intent(in)    :: cpij, cpkl
    integer,          intent(in)    :: idij, idkl
    real(real64),     intent(in)    :: dshift(3)
    real(real64),     intent(out)   :: gfull(0:, 0:, 0:, 0:, :, :)   ! (nj,ni,nl,nk,xyz,t)

    integer      :: t, n, m, ni, nj, nk, nl, nmax
    integer      :: maxij, maxkl, maxp1, maxp, nrow, ncol
    real(real64) :: f00, aandb, rho, expe
    real(real64) :: b00, b10, bp01, c10, cp01, cp10, c01
    real(real64) :: c00(3), cp00(3), dij(3), dkl(3), pq(3), pq0(3)
    real(real64), allocatable :: g(:,:,:)

    associate ( ppij => cpij%p(idij), ppkl => cpkl%p(idkl), &
                iang => cpij%iang,    jang => cpij%jang,    &
                kang => cpkl%iang,    lang => cpkl%jang )

    aandb = ppij%aa + ppkl%aa
    rho   = ppij%aa * ppkl%aa / aandb
    ! Operator displacement (e2 shifted by dshift) must enter ONLY the Boys argument (interelectronic
    ! separation). The `expe` prefactor carries an exp(-x/rho) Gaussian in |P-Q| that belongs to the
    ! fixed orbital geometry; letting dshift leak into it corrupts the displaced ERI (adds a spurious
    ! second-derivative term). So evaluate the Boys roots at the displaced separation, but build the
    ! prefactor from the undisplaced one.
    pq0   = ppij%r - ppkl%r                      ! undisplaced (fixed geometry)
    pq    = pq0 - dshift                         ! displaced separation (operator)
    ryscomp%x = rho * SUM(pq**2)
    call ryscomp%evaluate()

    maxij = iang + jang + 1        ! electron-1 padded by +1
    maxkl = kang + lang + 1        ! electron-2 padded by +1
    maxp1 = maxij + 1
    maxp  = maxkl + 1

    nrow = (maxij + 1)*(jang + 2)  ! safe upper bound for ni + maxp1*nj
    ncol = (maxkl + 1)*(lang + 2)  ! safe upper bound for nk + maxp *nl
    allocate(g(0:nrow, 0:ncol, 3))

    dij = cpij%ri - cpij%rj        ! A - B
    dkl = cpkl%ri - cpkl%rj        ! C - D
    expe = PI252 / (ppij%aa * ppkl%aa * SQRT(aandb)) * EXP(-rho * SUM(pq0**2) / rho)

    do t = 1, ryscomp%nroots
      g = 0.0_real64
      f00 = expe * ryscomp%w(t)
      associate (uu => ryscomp%u(t))
        b00  = uu*rho                 / (2*(ppij%aa*ppkl%aa + uu*rho*aandb))
        b10  = (ppkl%aa + uu*rho)     / (2*(ppij%aa*ppkl%aa + uu*rho*aandb))
        bp01 = (ppij%aa + uu*rho)     / (2*(ppij%aa*ppkl%aa + uu*rho*aandb))
      end associate
      ! VRR centres use the DISPLACED separation (they encode the operator-electron coupling that
      ! genuinely shifts with the operator). Only the `expe` Gaussian prefactor above is frozen at
      ! the undisplaced geometry -- it is a fixed-orbital normalisation, not part of the kernel.
      c00  = (ppij%r - cpij%ri) + 2*b00*ppkl%aa * pq
      cp00 = (ppkl%r - cpkl%ri) - 2*b00*ppij%aa * pq

      ! seeds
      g(0,0,1) = 1.0_real64 ; g(0,0,2) = 1.0_real64 ; g(0,0,3) = f00
      g(1,0,1) = c00(1)     ; g(1,0,2) = c00(2)     ; g(1,0,3) = c00(3)*f00
      g(0,1,1) = cp00(1)    ; g(0,1,2) = cp00(2)    ; g(0,1,3) = cp00(3)*f00
      g(1,1,1) = c00(1)*cp00(1) + b00
      g(1,1,2) = c00(2)*cp00(2) + b00
      g(1,1,3) = (c00(3)*cp00(3) + b00)*f00

      ! VRR electron 1 (N up, M=0,1)
      c10 = 0.0_real64 ; cp10 = b00
      do n = 2, maxij
        c10  = c10  + b10
        cp10 = cp10 + b00
        g(n,0,:) = c10 *g(n-2,0,:) + c00*g(n-1,0,:)
        g(n,1,:) = cp10*g(n-1,0,:) + cp00*g(n,0,:)
      end do

      ! VRR electron 2 (M up), full rectangle (no GAMESS truncation)
      cp01 = 0.0_real64 ; c01 = b00
      do m = 2, maxkl
        cp01 = cp01 + bp01
        c01  = c01  + b00
        g(0,m,:) = cp01*g(0,m-2,:) + cp00*g(0,m-1,:)
        g(1,m,:) = c01 *g(0,m-1,:) + c00 *g(0,m,:)
        cp10 = b00
        nmax = maxij
        do n = 2, nmax
          cp10 = cp10 + b00
          g(n,m,:) = cp01*g(n,m-2,:) + cp10*g(n-1,m-1,:) + cp00*g(n,m-1,:)
        end do
      end do

      ! HRR electron 1 (build NJ rows)
      do nj = 1, jang + 1
        do ni = maxij - nj, 0, -1
          g(ni + maxp1*nj, 0:maxkl, :) = &
              g(ni + maxp1*(nj-1) + 1, 0:maxkl, :) &
            + spread(dij, 1, maxkl+1) * g(ni + maxp1*(nj-1), 0:maxkl, :)
        end do
      end do

      ! HRR electron 2 (build NL cols) for every electron-1 packed row
      do nl = 1, lang + 1
        do nk = maxkl - nl, 0, -1
          do nj = 0, jang + 1
            ni = min(iang + 1, maxij - nj)
            if (ni < 0) cycle
            g(maxp1*nj : ni + maxp1*nj, nk + maxp*nl, :) = &
                g(maxp1*nj : ni + maxp1*nj, nk + maxp*(nl-1) + 1, :) &
              + spread(dkl, 1, ni + 1) * g(maxp1*nj : ni + maxp1*nj, nk + maxp*(nl-1), :)
          end do
        end do
      end do

      ! unpack into gfull(nj,ni,nl,nk,xyz,t)
      do nl = 0, lang + 1
        do nk = 0, kang + 1
          do nj = 0, jang + 1
            do ni = 0, iang + 1
              gfull(nj,ni,nl,nk,:,t) = g(ni + maxp1*nj, nk + maxp*nl, :)
            end do
          end do
        end do
      end do
    end do   ! roots

    end associate
  end subroutine qgauss_ss

!-------------------------------------------------------------------------------
!> Electron-1 first derivative of the 1D table for fixed direction d (acts on ni,nj),
!> returned as a per-root vector. di(m,n) = n*g(m,n-1) - 2ai*g(m,n+1) [bra]
!>                                         + m*g(m-1,n) - 2aj*g(m+1,n) [ket].
  pure function e1d(gfull, d, pj, pi, pl, pk, ai, aj, nr) result(v)
    real(real64), intent(in) :: gfull(0:,0:,0:,0:,:,:)
    integer,      intent(in) :: d, pj, pi, pl, pk, nr
    real(real64), intent(in) :: ai, aj
    real(real64) :: v(nr)
    v = -2*ai*gfull(pj,pi+1,pl,pk,d,1:nr) - 2*aj*gfull(pj+1,pi,pl,pk,d,1:nr)
    if (pi > 0) v = v + pi*gfull(pj,pi-1,pl,pk,d,1:nr)
    if (pj > 0) v = v + pj*gfull(pj-1,pi,pl,pk,d,1:nr)
  end function e1d

!> Electron-2 first derivative (acts on nk,nl).
  pure function e2d(gfull, d, pj, pi, pl, pk, ak, al, nr) result(v)
    real(real64), intent(in) :: gfull(0:,0:,0:,0:,:,:)
    integer,      intent(in) :: d, pj, pi, pl, pk, nr
    real(real64), intent(in) :: ak, al
    real(real64) :: v(nr)
    v = -2*ak*gfull(pj,pi,pl,pk+1,d,1:nr) - 2*al*gfull(pj,pi,pl+1,pk,d,1:nr)
    if (pk > 0) v = v + pk*gfull(pj,pi,pl,pk-1,d,1:nr)
    if (pl > 0) v = v + pl*gfull(pj,pi,pl-1,pk,d,1:nr)
  end function e2d

!> Mixed electron-1 x electron-2 second derivative in the SAME direction d (diagonal kl).
  pure function e12d(gfull, d, pj, pi, pl, pk, ai, aj, ak, al, nr) result(v)
    real(real64), intent(in) :: gfull(0:,0:,0:,0:,:,:)
    integer,      intent(in) :: d, pj, pi, pl, pk, nr
    real(real64), intent(in) :: ai, aj, ak, al
    real(real64) :: v(nr)
    v = -2*ai*e2d(gfull,d,pj,pi+1,pl,pk,ak,al,nr) &
        -2*aj*e2d(gfull,d,pj+1,pi,pl,pk,ak,al,nr)
    if (pi > 0) v = v + pi*e2d(gfull,d,pj,pi-1,pl,pk,ak,al,nr)
    if (pj > 0) v = v + pj*e2d(gfull,d,pj-1,pi,pl,pk,ak,al,nr)
  end function e12d

!-------------------------------------------------------------------------------
!> Analytic SS bare-Hessian integral block for one primitive quartet.
!> ssblk(ij, kl, c) for c = 1..6 = (xx, yy, zz, xy, xz, yz);
!> ij = (i-1)*jnao + j over electron-1 cart functions, kl = (k-1)*lnao + l over electron 2.
!> H_kl = -dij_fac * sum_t [factors], with d_k on the electron-1 cloud and d_l on electron 2.
  subroutine comp_ssc_int2_prim(cpij, idij, cpkl, idkl, ssblk)
    type(shpair_t), intent(in)    :: cpij, cpkl
    integer,        intent(in)    :: idij, idkl
    real(real64),   intent(out)   :: ssblk(:, :, :)   ! (inao*jnao, knao*lnao, 6)

    type(rys_root_t) :: ryscomp
    real(real64), allocatable :: gfull(:,:,:,:,:,:)
    integer :: nr, ij, kl, i, j, k, l
    integer :: nxi,nyi,nzi, nxj,nyj,nzj, nxk,nyk,nzk, nxl,nyl,nzl
    real(real64) :: dij_fac
    real(real64) :: vx(64), vy(64), vz(64)   ! per-root work (64 >> any nroots)

    associate ( iang => cpij%iang, jang => cpij%jang, &
                kang => cpkl%iang, lang => cpkl%jang, &
                inao => cpij%inao, jnao => cpij%jnao, &
                knao => cpkl%inao, lnao => cpkl%jnao, &
                ai => cpij%p(idij)%ai, aj => cpij%p(idij)%aj, &
                ak => cpkl%p(idkl)%ai, al => cpkl%p(idkl)%aj )

    nr = ssc_nroots(iang, jang, kang, lang)
    ryscomp%nroots = nr
    allocate(gfull(0:jang+1, 0:iang+1, 0:lang+1, 0:kang+1, 3, nr))
    call qgauss_ss(ryscomp, cpij, idij, cpkl, idkl, [0._real64,0._real64,0._real64], gfull)

    dij_fac = cpij%p(idij)%expfac * cpkl%p(idkl)%expfac
    ssblk = 0.0_real64

    do i = 1, inao
      nxi = CART_X(i,iang); nyi = CART_Y(i,iang); nzi = CART_Z(i,iang)
      do j = 1, jnao
        nxj = CART_X(j,jang); nyj = CART_Y(j,jang); nzj = CART_Z(j,jang)
        ij = (i-1)*jnao + j
        do k = 1, knao
          nxk = CART_X(k,kang); nyk = CART_Y(k,kang); nzk = CART_Z(k,kang)
          do l = 1, lnao
            nxl = CART_X(l,lang); nyl = CART_Y(l,lang); nzl = CART_Z(l,lang)
            kl = (k-1)*lnao + l

            ! xx : d_x(e1) d_x(e2) on x; plain y,z
            vx(1:nr) = e12d(gfull,1,nxj,nxi,nxl,nxk,ai,aj,ak,al,nr)
            vy(1:nr) = gfull(nyj,nyi,nyl,nyk,2,1:nr)
            vz(1:nr) = gfull(nzj,nzi,nzl,nzk,3,1:nr)
            ssblk(ij,kl,1) = -dij_fac*sum(vx(1:nr)*vy(1:nr)*vz(1:nr))

            ! yy
            vx(1:nr) = gfull(nxj,nxi,nxl,nxk,1,1:nr)
            vy(1:nr) = e12d(gfull,2,nyj,nyi,nyl,nyk,ai,aj,ak,al,nr)
            vz(1:nr) = gfull(nzj,nzi,nzl,nzk,3,1:nr)
            ssblk(ij,kl,2) = -dij_fac*sum(vx(1:nr)*vy(1:nr)*vz(1:nr))

            ! zz
            vx(1:nr) = gfull(nxj,nxi,nxl,nxk,1,1:nr)
            vy(1:nr) = gfull(nyj,nyi,nyl,nyk,2,1:nr)
            vz(1:nr) = e12d(gfull,3,nzj,nzi,nzl,nzk,ai,aj,ak,al,nr)
            ssblk(ij,kl,3) = -dij_fac*sum(vx(1:nr)*vy(1:nr)*vz(1:nr))

            ! xy : d_x(e1) on x, d_y(e2) on y
            vx(1:nr) = e1d(gfull,1,nxj,nxi,nxl,nxk,ai,aj,nr)
            vy(1:nr) = e2d(gfull,2,nyj,nyi,nyl,nyk,ak,al,nr)
            vz(1:nr) = gfull(nzj,nzi,nzl,nzk,3,1:nr)
            ssblk(ij,kl,4) = -dij_fac*sum(vx(1:nr)*vy(1:nr)*vz(1:nr))

            ! xz : d_x(e1) on x, d_z(e2) on z
            vx(1:nr) = e1d(gfull,1,nxj,nxi,nxl,nxk,ai,aj,nr)
            vy(1:nr) = gfull(nyj,nyi,nyl,nyk,2,1:nr)
            vz(1:nr) = e2d(gfull,3,nzj,nzi,nzl,nzk,ak,al,nr)
            ssblk(ij,kl,5) = -dij_fac*sum(vx(1:nr)*vy(1:nr)*vz(1:nr))

            ! yz : d_y(e1) on y, d_z(e2) on z
            vx(1:nr) = gfull(nxj,nxi,nxl,nxk,1,1:nr)
            vy(1:nr) = e1d(gfull,2,nyj,nyi,nyl,nyk,ai,aj,nr)
            vz(1:nr) = e2d(gfull,3,nzj,nzi,nzl,nzk,ak,al,nr)
            ssblk(ij,kl,6) = -dij_fac*sum(vx(1:nr)*vy(1:nr)*vz(1:nr))
          end do
        end do
      end do
    end do

    end associate
  end subroutine comp_ssc_int2_prim

!-------------------------------------------------------------------------------
!> Plain Coulomb ERI block for one primitive quartet with electron 2 rigidly displaced by
!> `dshift` (= displacing the 1/r12 operator). eriblk(ij, kl). Used as the FD reference:
!> H_kl = d^2/d dshift_k d dshift_l [ ERI(dshift) ] at dshift=0.
  subroutine comp_eri2_prim_disp(cpij, idij, cpkl, idkl, dshift, eriblk)
    type(shpair_t), intent(in)  :: cpij, cpkl
    integer,        intent(in)  :: idij, idkl
    real(real64),   intent(in)  :: dshift(3)
    real(real64),   intent(out) :: eriblk(:, :)      ! (inao*jnao, knao*lnao)

    type(rys_root_t) :: ryscomp
    real(real64), allocatable :: gfull(:,:,:,:,:,:)
    integer :: nr, ij, kl, i, j, k, l
    integer :: nxi,nyi,nzi, nxj,nyj,nzj, nxk,nyk,nzk, nxl,nyl,nzl
    real(real64) :: dij_fac

    associate ( iang => cpij%iang, jang => cpij%jang, &
                kang => cpkl%iang, lang => cpkl%jang, &
                inao => cpij%inao, jnao => cpij%jnao, &
                knao => cpkl%inao, lnao => cpkl%jnao )

    nr = ssc_nroots(iang, jang, kang, lang)
    ryscomp%nroots = nr
    allocate(gfull(0:jang+1, 0:iang+1, 0:lang+1, 0:kang+1, 3, nr))
    call qgauss_ss(ryscomp, cpij, idij, cpkl, idkl, dshift, gfull)

    dij_fac = cpij%p(idij)%expfac * cpkl%p(idkl)%expfac
    eriblk = 0.0_real64

    do i = 1, inao
      nxi = CART_X(i,iang); nyi = CART_Y(i,iang); nzi = CART_Z(i,iang)
      do j = 1, jnao
        nxj = CART_X(j,jang); nyj = CART_Y(j,jang); nzj = CART_Z(j,jang)
        ij = (i-1)*jnao + j
        do k = 1, knao
          nxk = CART_X(k,kang); nyk = CART_Y(k,kang); nzk = CART_Z(k,kang)
          do l = 1, lnao
            nxl = CART_X(l,lang); nyl = CART_Y(l,lang); nzl = CART_Z(l,lang)
            kl = (k-1)*lnao + l
            eriblk(ij,kl) = dij_fac * sum( gfull(nxj,nxi,nxl,nxk,1,1:nr) &
                                         * gfull(nyj,nyi,nyl,nyk,2,1:nr) &
                                         * gfull(nzj,nzi,nzl,nzk,3,1:nr) )
          end do
        end do
      end do
    end do

    end associate
  end subroutine comp_eri2_prim_disp

end module mod_ssc_int2
