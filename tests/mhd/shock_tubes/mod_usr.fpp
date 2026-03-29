module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

contains

  subroutine usr_init()

    usr_init_one_grid => initonegrid_usr

    call set_coordinate_system("Cartesian_3D")

    call phys_activate()

  end subroutine usr_init


  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)

    double precision :: rhoL, rhoR, pL, pR
    double precision :: vxL, vxR, vyL, vyR, vzL, vzR
    double precision :: BxL, BxR, ByL, ByR, BzL, BzR
    double precision :: x0

    ! Discontinuity location
    x0 = 0.5_dp

    select case(iprob)

    case(1)
      ! Brio-Wu (1988), gamma = 2.0
      rhoL = 1.0_dp;       rhoR = 0.125_dp
      pL   = 1.0_dp;       pR   = 0.1_dp
      vxL  = zero;          vxR  = zero
      vyL  = zero;          vyR  = zero
      vzL  = zero;          vzR  = zero
      BxL  = 0.75_dp;       BxR  = 0.75_dp
      ByL  = 1.0_dp;        ByR  = -1.0_dp
      BzL  = zero;           BzR  = zero

    case(2)
      ! Ryu-Jones 1a (1995), gamma = 5/3
      rhoL = 1.0_dp;       rhoR = 0.1_dp
      pL   = 20.0_dp;      pR   = 1.0_dp
      vxL  = 10.0_dp;      vxR  = zero
      vyL  = zero;          vyR  = zero
      vzL  = zero;          vzR  = zero
      BxL  = 5.0_dp / sqrt(4.0_dp * dpi)
      BxR  = BxL
      ByL  = 5.0_dp / sqrt(4.0_dp * dpi)
      ByR  = 5.0_dp / sqrt(4.0_dp * dpi)
      BzL  = zero;           BzR  = zero

    case default
      error stop "Unknown iprob for Riemann test"

    end select

    ! Set left and right states
    where (x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1) < x0)
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)   = rhoL
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(1))  = vxL
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(2))  = vyL
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(3))  = vzL
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,p_)      = pL
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(1))  = BxL
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(2))  = ByL
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(3))  = BzL
    elsewhere
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)   = rhoR
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(1))  = vxR
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(2))  = vyR
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(3))  = vzR
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,p_)      = pR
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(1))  = BxR
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(2))  = ByR
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(3))  = BzR
    end where

    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

end module mod_usr
