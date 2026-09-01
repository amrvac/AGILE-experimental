! Minimal reproducer for AGILE issue #154
!   https://github.com/amrvac/AGILE/issues/154
!   "calling to_primitive or to_conservative from specialbound"
! and the nvfortran OpenACC miscompile fixed in commit bb31e7ad.
!
! WHAT IT SHOWS
!   specialbound_usr fills a block's boundary ghost cells from an analytic
!   state on the GPU, in the same OpenACC call chain the *_pole test cases use:
!
!     !$acc parallel loop gang                  (over blocks)
!       fill_bnd          !$acc routine vector   <- fill_boundary_before_gc
!        bc_phys          !$acc routine vector
!         specialbound_usr  !$acc routine vector
!            !$acc loop collapse(3) vector       (over the nghost-wide ghost slab)
!              call analytic_state !$acc routine seq   <- to_primitive
!
!   The driver runs FOUR variants of that innermost loop and checks each
!   against the analytic state computed on the host:
!
!     variant 1  BUGGY       !$acc loop collapse(3) + call analytic_state,
!                            ghost slab sized by the RUNTIME nw
!     variant 2  WORKAROUND  !$acc loop vector on the innermost loop, no
!                            collapse, outer loops sequential   (call kept)
!     variant 3  WORKAROUND  !$acc loop collapse(3), analytic_state inlined
!                            so the collapsed body has no call
!     variant 4  WORKAROUND  identical to variant 1 (collapse(3) + call) except
!                            the slab/wpt are sized by a compile-time parameter
!                            nw_phys instead of the runtime nw
!
!   gfortran -fopenacc and nvfortran without -acc are exact for all four
!   variants - which is what makes the failures a compiler bug, not a bad
!   program.  On nvfortran + OpenACC, as shipped (radius_dependent_flow =
!   .true.): variants 1 and 2 are wrong, variants 3 and 4 exact.
!
!     - variant 1  (collapse + call) is off by O(1): the ghost layer is
!       index-rotated.
!     - variant 2  (drop collapse, keep call) swaps the two iterations of the
!       trip-2 radial `!$acc loop vector`.  With radius_dependent_flow =
!       .false. - the *_pole cases' situation, where the boundary state does
!       not vary with radius - i1=1 and i1=2 hold the same value, the swap is
!       invisible, and variant 2 passes.  That is why dropping the collapse
!       is the fix those cases ship; it is not robust in general.
!     - variant 3  inlines the callee: no `call` in the vector loop, exact
!       either way.  Removing the call is the robust fix.
!     - variant 4  is variant 1 with the slab/wpt sized by a compile-time
!       parameter (nw_phys) instead of the runtime nw.  Not a practical fix
!       for AGILE (nw is genuinely runtime), but it isolates the trigger: a
!       runtime extent inside the collapsed nest is part of what nvfortran
!       gets wrong.
!
!   So the `call` from inside the vectorised loop is the real trigger, not
!   the `collapse` on its own:
!
!                          | collapse(3)          | vector on innermost loop
!       ------------------- | -------------------- | ------------------------
!       call analytic_state | FAIL (index rotate)  | FAIL (adjacent-cell swap)
!       callee inlined      | pass                 | pass
!
!   Variants 2 and 3 are the two fixes noted in CLAUDE.md ("Bug-hunting
!   notes"), variant 2 being what the *_pole cases now use.
!
! BUILD & RUN  (see run.sh)
!   gfortran  -O2 -fopenacc                                repro154.f90 -o r && ./r
!   nvfortran -O3 -fast                                    repro154.f90 -o r && ./r
!   nvfortran -O3 -fast -acc=gpu -Mvect=levels:5 -Minline  repro154.f90 -o r && ./r

module repro154_mod
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: nblocks = 8
  integer, parameter :: nghost  = 2                 ! ghost-slab width (degenerate index)
  integer, parameter :: nx      = 4                 ! interior cells per dim
  integer, parameter :: nt      = nx + 2*nghost     ! full block width incl. ghosts

  integer :: nw = 3                                 ! runtime, as in AGILE
  integer, parameter :: nw_phys = 3                 ! same value, compile-time
  integer :: variant = 1                            ! 1 buggy, 2/3/4 workarounds

  ! .false.: boundary state depends on the angular coordinates only, as in the
  ! *_pole cases - variant 2's residual radial-loop swap is then invisible.
  ! .true.:  state also varies with the radius x(1) - variant 2 then FAILs
  ! too, variants 3 and 4 stay exact.
  logical :: radius_dependent_flow = .true.
  !$acc declare create(nw, variant, radius_dependent_flow)
contains

  !> Stands in for analytic_state / to_primitive: the second `!$acc routine`,
  !> sized/shaped like the real one so -Minline leaves the call in place.
  pure subroutine analytic_state(x, w)
    !$acc routine seq
    real(dp), intent(in)  :: x(3)
    real(dp), intent(out) :: w(nw)
    real(dp) :: sinp, cosp, v(3), rho
    ! A uniform Cartesian vector written in curvilinear components: with
    ! radius_dependent_flow = .false. it is a function of the angular
    ! coordinates x(2), x(3) only, exactly as in the *_pole cases.
    sinp = sin(x(3)); cosp = cos(x(3))
    v(1) =  0.7_dp*cosp + 0.2_dp*sinp
    v(2) = -0.7_dp*sinp + 0.2_dp*cosp
    v(3) =  0.3_dp*x(2)
    rho = 1.0_dp
    if (radius_dependent_flow) rho = 1.0_dp + 0.5_dp*x(1)
    w(1) = rho
    w(2) = v(1) + 2.0_dp*v(3)                       ! varies with x(2) and x(3)
    w(3) = 1.0_dp + 0.5_dp*(v(1)**2 + v(2)**2 + v(3)**2)
  end subroutine analytic_state

  ! Three stand-ins for specialbound_usr, each an `!$acc routine vector` with
  ! a SINGLE loop nest (as the real specialbound_usr has).  Runtime ixI (the
  ! array) and ixO (the loop) bounds; ixO is a corner-shaped strict sub-box.

  !> variant 1  BUGGY: collapse(3) with a call to analytic_state in the body.
  subroutine sb_collapse_call(ilo1,ihi1,ilo2,ihi2,ilo3,ihi3, &
                              olo1,ohi1,olo2,ohi2,olo3,ohi3, w, x)
    !$acc routine vector
    integer,  intent(in)    :: ilo1,ihi1,ilo2,ihi2,ilo3,ihi3
    integer,  intent(in)    :: olo1,ohi1,olo2,ohi2,olo3,ohi3
    real(dp), intent(inout) :: w(ilo1:ihi1,ilo2:ihi2,ilo3:ihi3,1:nw)
    real(dp), intent(in)    :: x(ilo1:ihi1,ilo2:ihi2,ilo3:ihi3,1:3)
    integer  :: ix1, ix2, ix3
    real(dp) :: wpt(nw), x_loc(3)
    !$acc loop collapse(3) vector private(wpt, x_loc)
    do ix3 = olo3, ohi3
       do ix2 = olo2, ohi2
          do ix1 = olo1, ohi1
             x_loc(1:3)          = x(ix1,ix2,ix3,1:3)
             call analytic_state(x_loc, wpt)
             w(ix1,ix2,ix3,1:nw) = wpt(1:nw)
          end do
       end do
    end do
  end subroutine sb_collapse_call

  !> variant 2  WORKAROUND: no collapse, `!$acc loop vector` on the innermost
  !> loop, outer loops sequential (call kept).  This is what the *_pole cases
  !> now use.
  subroutine sb_vector_call(ilo1,ihi1,ilo2,ihi2,ilo3,ihi3, &
                            olo1,ohi1,olo2,ohi2,olo3,ohi3, w, x)
    !$acc routine vector
    integer,  intent(in)    :: ilo1,ihi1,ilo2,ihi2,ilo3,ihi3
    integer,  intent(in)    :: olo1,ohi1,olo2,ohi2,olo3,ohi3
    real(dp), intent(inout) :: w(ilo1:ihi1,ilo2:ihi2,ilo3:ihi3,1:nw)
    real(dp), intent(in)    :: x(ilo1:ihi1,ilo2:ihi2,ilo3:ihi3,1:3)
    integer  :: ix1, ix2, ix3
    real(dp) :: wpt(nw), x_loc(3)
    do ix3 = olo3, ohi3
       do ix2 = olo2, ohi2
          !$acc loop vector private(wpt, x_loc)
          do ix1 = olo1, ohi1
             x_loc(1:3)          = x(ix1,ix2,ix3,1:3)
             call analytic_state(x_loc, wpt)
             w(ix1,ix2,ix3,1:nw) = wpt(1:nw)
          end do
       end do
    end do
  end subroutine sb_vector_call

  !> variant 3  WORKAROUND: collapse(3) kept, analytic_state inlined so the
  !> collapsed body has no call.
  subroutine sb_collapse_inline(ilo1,ihi1,ilo2,ihi2,ilo3,ihi3, &
                                olo1,ohi1,olo2,ohi2,olo3,ohi3, w, x)
    !$acc routine vector
    integer,  intent(in)    :: ilo1,ihi1,ilo2,ihi2,ilo3,ihi3
    integer,  intent(in)    :: olo1,ohi1,olo2,ohi2,olo3,ohi3
    real(dp), intent(inout) :: w(ilo1:ihi1,ilo2:ihi2,ilo3:ihi3,1:nw)
    real(dp), intent(in)    :: x(ilo1:ihi1,ilo2:ihi2,ilo3:ihi3,1:3)
    integer  :: ix1, ix2, ix3
    real(dp) :: x_loc(3), sinp, cosp, v(3)
    !$acc loop collapse(3) vector private(x_loc, sinp, cosp, v)
    do ix3 = olo3, ohi3
       do ix2 = olo2, ohi2
          do ix1 = olo1, ohi1
             x_loc(1:3) = x(ix1,ix2,ix3,1:3)
             sinp = sin(x_loc(3)); cosp = cos(x_loc(3))
             v(1) =  0.7_dp*cosp + 0.2_dp*sinp
             v(2) = -0.7_dp*sinp + 0.2_dp*cosp
             v(3) =  0.3_dp*x_loc(2)
             w(ix1,ix2,ix3,1) = merge(1.0_dp + 0.5_dp*x_loc(1), 1.0_dp, &
                                      radius_dependent_flow)
             w(ix1,ix2,ix3,2) = v(1) + 2.0_dp*v(3)
             w(ix1,ix2,ix3,3) = 1.0_dp + 0.5_dp*(v(1)**2 + v(2)**2 + v(3)**2)
          end do
       end do
    end do
  end subroutine sb_collapse_inline

  !> variant 4  WORKAROUND: byte-for-byte variant 1 (collapse(3) + call), the
  !> ONLY change being that the ghost-slab array and the wpt private buffer
  !> are sized by the compile-time parameter nw_phys instead of the runtime
  !> nw.  That is enough to make nvfortran generate the collapsed index
  !> arithmetic correctly - a runtime trip/extent in the collapsed nest is
  !> part of what trips it.
  subroutine sb_collapse_param(ilo1,ihi1,ilo2,ihi2,ilo3,ihi3, &
                               olo1,ohi1,olo2,ohi2,olo3,ohi3, w, x)
    !$acc routine vector
    integer,  intent(in)    :: ilo1,ihi1,ilo2,ihi2,ilo3,ihi3
    integer,  intent(in)    :: olo1,ohi1,olo2,ohi2,olo3,ohi3
    real(dp), intent(inout) :: w(ilo1:ihi1,ilo2:ihi2,ilo3:ihi3,1:nw_phys)
    real(dp), intent(in)    :: x(ilo1:ihi1,ilo2:ihi2,ilo3:ihi3,1:3)
    integer  :: ix1, ix2, ix3
    real(dp) :: wpt(nw_phys), x_loc(3)
    !$acc loop collapse(3) vector private(wpt, x_loc)
    do ix3 = olo3, ohi3
       do ix2 = olo2, ohi2
          do ix1 = olo1, ohi1
             x_loc(1:3)               = x(ix1,ix2,ix3,1:3)
             call analytic_state(x_loc, wpt)
             w(ix1,ix2,ix3,1:nw_phys) = wpt(1:nw_phys)
          end do
       end do
    end do
  end subroutine sb_collapse_param

  !> Stands in for bc_phys: computes the ghost-slab bounds for face iB at
  !> runtime, then dispatches to the selected variant - another
  !> `!$acc routine vector` level.
  subroutine bc_phys(iB, w, x)
    !$acc routine vector
    integer,  intent(in)    :: iB
    real(dp), intent(inout) :: w(nt,nt,nt,nw)
    real(dp), intent(in)    :: x(nt,nt,nt,3)
    integer :: o1lo, o1hi
    o1lo = 1; o1hi = nt
    if (iB == 1) o1hi = nghost               ! r-min ghost slab, nghost wide
    if (iB == 2) o1lo = nt - nghost + 1      ! r-max ghost slab
    select case (variant)
    case (1); call sb_collapse_call  (1,nt,1,nt,1,nt, o1lo,o1hi,1,nt,1,nt, w, x)
    case (2); call sb_vector_call    (1,nt,1,nt,1,nt, o1lo,o1hi,1,nt,1,nt, w, x)
    case (3); call sb_collapse_inline(1,nt,1,nt,1,nt, o1lo,o1hi,1,nt,1,nt, w, x)
    case (4); call sb_collapse_param (1,nt,1,nt,1,nt, o1lo,o1hi,1,nt,1,nt, w, x)
    end select
  end subroutine bc_phys

  !> Stands in for fill_boundary_before_gc: outermost `!$acc routine vector`,
  !> called straight from the gang loop, with the idims/iside seq loops.
  subroutine fill_bnd(w, x)
    !$acc routine vector
    real(dp), intent(inout) :: w(nt,nt,nt,nw)
    real(dp), intent(in)    :: x(nt,nt,nt,3)
    integer :: idims, iside
    do idims = 1, 3
       do iside = 1, 2
          if (idims == 1) call bc_phys(iside, w, x)
       end do
    end do
  end subroutine fill_bnd

  subroutine run_variant(v, ok)
    integer, intent(in)  :: v
    logical, intent(out) :: ok
    real(dp) :: w(nt,nt,nt,nw,nblocks), x(nt,nt,nt,3,nblocks), ref(nw)
    integer  :: ib, i1, i2, i3, nbad
    real(dp) :: maxerr

    variant = v
    !$acc update device(nw, variant, radius_dependent_flow)

    do ib = 1, nblocks
       do i3 = 1, nt
          do i2 = 1, nt
             do i1 = 1, nt
                x(i1,i2,i3,1,ib) = 0.5_dp + 0.1_dp*ib + 0.01_dp*i1   ! radius
                x(i1,i2,i3,2,ib) = 0.15_dp*i2 + 0.03_dp*ib           ! theta / z
                x(i1,i2,i3,3,ib) = 0.20_dp*i3                        ! phi
             end do
          end do
       end do
    end do
    w = -1.0_dp

    !$acc parallel loop gang copyin(x) copyout(w)
    do ib = 1, nblocks
       call fill_bnd(w(:,:,:,:,ib), x(:,:,:,:,ib))
    end do

    maxerr = 0.0_dp; nbad = 0
    do ib = 1, nblocks
       do i3 = 1, nt
          do i2 = 1, nt
             do i1 = 1, nghost                     ! the r-min ghost slab we filled
                call analytic_state(x(i1,i2,i3,1:3,ib), ref)
                if (maxval(abs(w(i1,i2,i3,1:nw,ib) - ref)) > 1.0e-10_dp) then
                   nbad = nbad + 1
                   maxerr = max(maxerr, maxval(abs(w(i1,i2,i3,1:nw,ib) - ref)))
                end if
             end do
          end do
       end do
    end do

    ok = (nbad == 0)
    print '(a,i0,a,i0,a,i0,a,es9.2,a)', '  variant ', v, ':  wrong cells ', &
         nbad, ' / ', nghost*nt*nt*nblocks, ',  max error ', maxerr, &
         merge('   PASS', '   FAIL', ok)
  end subroutine run_variant

  subroutine driver()
    logical :: ok(4)
    integer :: v
    print '(a)', 'issue #154 reproducer  (1 = buggy, 2/3/4 = workarounds)'
    do v = 1, 4
       call run_variant(v, ok(v))
    end do
    if (all(ok)) then
       print '(a)', 'no miscompile seen on this toolchain'
    else
       if (all(ok(2:4))) then
          print '(a)', 'issue #154 reproduced: variant 1 wrong, workarounds 2-4 fix it'
       else
          print '(a)', 'issue #154 reproduced: variant 1 wrong; a workaround also &
               &misbehaves on this toolchain'
       end if
       call exit(1)
    end if
  end subroutine driver
end module repro154_mod

program repro154
  use repro154_mod
  implicit none
  call driver()
end program repro154
