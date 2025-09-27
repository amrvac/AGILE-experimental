! setup.pl -d=2
! test thermal conduction in a ring in ffhd module
module mod_usr
  use mod_amrvac
  use mod_physics
  implicit none

contains

  subroutine usr_init()
    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 !cm-3

    usr_init_one_grid   => initonegrid_usr

    call set_coordinate_system("Cartesian")
    call phys_activate()

  end subroutine usr_init

  subroutine initonegrid_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)

    double precision :: r(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        theta(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    ! set all velocity to zero
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(1)) = zero
    ! uniform pressure
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_) =1.d0
    r(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=dsqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)**2+x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)**2)
    where(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)>0.d0)
      theta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=atan(x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1))
    elsewhere(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)<0.d0)
      theta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=dpi-atan(x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)/abs(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)))
    elsewhere
      theta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=0.d0
    endwhere
    ! hot central circular spot with uniform pressure
    where(r(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)>0.5d0 .and. r(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)<0.7d0 .and. theta(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)>11.d0/12.d0*dpi .and. theta(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)<13.d0/12.d0*dpi)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)=4.d0
    elsewhere
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)=1.d0
    endwhere

    call phys_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)

  end subroutine initonegrid_usr

  !> analytical formula for the unit vectors along B
  pure real(dp) function bfield(x, idim) result(field)
    !$acc routine seq
    real(dp), intent(in)    :: x(1:ndim)
    integer, value, intent(in)     :: idim
    real(dp) :: r, theta, wB0(2)

    r=dsqrt(x(1)**2+x(2)**2)
    if(x(1)>0.d0) then
      theta=atan(x(2)/x(1))
    else if(x(1)<0.d0) then
      theta=dpi-atan(x(2)/abs(x(1)))
    else
      theta=0.d0
    end if
    ! straight line current
    wB0(1)=dcos(theta+0.5*dpi)/r
    wB0(2)=dsin(theta+0.5*dpi)/r
    if(idim==1) field=wB0(1)/dsqrt(wB0(1)**2+wB0(2)**2)
    if(idim==2) field=wB0(2)/dsqrt(wB0(1)**2+wB0(2)**2)
    if(idim==3) field=0.d0

  end function bfield

end module mod_usr
