! Dummy routines which can be overwritten by a physics-dependent implementation

#:def get_cs2()
!> obtain the squared sound speed
pure real(dp) function get_cs2(u) result(cs2)
  !$acc routine seq
  real(dp), intent(in)  :: u(nw_phys)

  ! provide nonsensical soundspeed so the user notices:
  cs2 = -1.0d0
  
end function get_cs2
#:enddef


#:def estimate_speeds_minmax()
subroutine estimate_speeds_minmax(uL, uR, xC, flux_dim, wL, wR)
  !$acc routine seq
  real(dp), intent(in)  :: uL(nw_phys), uR(nw_phys)
  real(dp), intent(in)  :: xC(ndim)
  integer, intent(in)   :: flux_dim
  real(dp), intent(out) :: wL, wR

  wL = wR = -1._dp

end subroutine estimate_speeds_minmax
#:enddef


#:def estimate_speeds_toro_pvrs()
subroutine estimate_speeds_toro_pvrs(uL, uR, xC, flux_dim, sL, sR)
  !$acc routine seq
  real(dp), intent(in)  :: uL(nw_phys), uR(nw_phys)
  real(dp), intent(in)  :: xC(ndim)
  integer,  intent(in)  :: flux_dim
  real(dp), intent(out) :: sL, sR

  sL = sR = -1._dp

end subroutine estimate_speeds_toro_pvrs
#:enddef


#:def addsource_nonlocal()
subroutine addsource_nonlocal(qdt, dtfactor, qtC, wCTprim, qt, wnew, x, dx, idir, &
     qsourcesplit)
  !$acc routine seq

  real(dp), intent(in)     :: qdt, dtfactor, qtC, qt
  real(dp), intent(in)     :: wCTprim(nw_phys,5)
  real(dp), intent(in)     :: x(1:ndim), dx(1:ndim)
  real(dp), intent(inout)  :: wnew(nw_phys)
  integer, intent(in)      :: idir
  logical, intent(in)      :: qsourcesplit


end subroutine addsource_nonlocal
#:enddef
