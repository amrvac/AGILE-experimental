module mod_usr
  use mod_rho
  use mod_particles

  implicit none

  integer :: i_sol,i_err ! indices for extra variables
  integer :: isol   ! index for extra payload in gridvars for particles

contains

  subroutine usr_init()
    use mod_variables

    usr_init_one_grid => initonegrid_usr
    usr_process_grid => store_sol_err
    usr_print_log => print_min_max

    ! Particle stuff
    usr_create_particles => generate_particles
    usr_particle_position => move_particle
    particles_define_additional_gridvars => define_additional_gridvars_usr
    particles_fill_additional_gridvars => fill_additional_gridvars_usr
    usr_update_payload => update_payload_usr

    call set_coordinate_system("Cartesian_2D")
    call rho_activate()

    ! extra variables to store solution and error (for I/O purposes)
    i_sol = var_set_extravar("solution", "solution")
    i_err = var_set_extravar("error", "error")
  end subroutine usr_init

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid 
    use mod_physics
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: rhoprofile(ixI^S)
    logical, save:: first=.true.

    if (first) then
       if (mype==0) then
          if(ndim==2)then
             print *,'advection of VAC logo in 2D periodic box'
          else
             call mpistop("VAC logo advection is 2D, setup.pl -d=2")
          endif
       end if
       first=.false.
    end if
    call set_density_profile(ixI^L,ixO^L,0.0d0,x,rhoprofile)
    w(ixO^S,rho_)=rhoprofile(ixO^S)

  end subroutine initonegrid_usr

  subroutine set_density_profile(ixI^L,ixO^L,qt,x,rhoprofile)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: qt,x(ixI^S,1:ndim)
    double precision, intent(out) :: rhoprofile(ixI^S)

    logical::          maskv(ixI^S),maska(ixI^S),maskc(ixI^S)
    double precision:: rhoflat,rhosquare, xshift(ixI^S,1:ndim)
    double precision:: xc1,yc1,xa1,xa2,ya1,ya2,xb1,xb2,yb1,yb2,xc2,yc2, &
         rad,rad2,alp,nsig
    !----------------------------------------------------------------------------

    rhoflat  = 0.5d0 
    rhosquare= 2.0d0 
    {^IFONED   call mpistop("problem is 2D!") }
    {^IFTHREED call mpistop("problem is 2D!") }
    {^IFTWOD
     ! account for periodicity, at least during one cycle....
     xshift(ixO^S,1)=x(ixO^S,1)-rho_v(1)*qt
     xshift(ixO^S,2)=x(ixO^S,2)-rho_v(2)*qt
     ! when v_x,v_y positive
     maskv(ixO^S)=(xshift(ixO^S,1)<xprobmin1)
     where(maskv(ixO^S)) xshift(ixO^S,1)=x(ixO^S,1)-rho_v(1)*qt+(xprobmax1-xprobmin1)
     maskv(ixO^S)=(xshift(ixO^S,2)<xprobmin2)
     where(maskv(ixO^S)) xshift(ixO^S,2)=x(ixO^S,2)-rho_v(2)*qt+(xprobmax2-xprobmin2)
     ! when v_x,v_y negative
     maskv(ixO^S)=(xshift(ixO^S,1)>xprobmax1)
     where(maskv(ixO^S)) xshift(ixO^S,1)=x(ixO^S,1)-rho_v(1)*qt-(xprobmax1-xprobmin1)
     maskv(ixO^S)=(xshift(ixO^S,2)>xprobmax2)
     where(maskv(ixO^S)) xshift(ixO^S,2)=x(ixO^S,2)-rho_v(2)*qt-(xprobmax2-xprobmin2)
     xc1=0.25d0
     yc1=0.50d0
     rad=0.23d0
     rad2=0.13d0
     alp=dpi/3.0d0
     xa1=xc1
     ya1=yc1-rad
     xa2=xc1-rad*cos(alp)
     ya2=yc1+rad*sin(alp)
     xb1=xa1
     yb1=ya1
     xb2=xc1+rad*cos(alp)
     yb2=yc1+rad*sin(alp)
     xc2=xc1
     yc2=ya2+sqrt(rad2**2-(xa2-xc2)**2)
     maskv(ixO^S)= ((xshift(ixO^S,1)-xc1)**2+(xshift(ixO^S,2)-yc1)**2 <= rad**2) &
          .and.(xshift(ixO^S,2)>= (ya2-ya1)*(xshift(ixO^S,1)-xa1)/(xa2-xa1)+ya1) & 
          .and.(xshift(ixO^S,2)>= (yb2-yb1)*(xshift(ixO^S,1)-xb1)/(xb2-xb1)+yb1) & 
          .and.((xshift(ixO^S,1)-xc2)**2+(xshift(ixO^S,2)-yc2)**2 > rad2**2) 
     xc1=0.45d0
     yc1=0.475d0
     xa1=xc1
     ya1=yc1+rad
     xa2=xc1-rad*cos(alp)
     ya2=yc1-rad*sin(alp)
     xb1=xa1
     yb1=ya1
     xb2=xc1+rad*cos(alp)
     yb2=yc1-rad*sin(alp)
     xc2=xc1
     yc2=ya2-sqrt(rad2**2-(xa2-xc2)**2)
     maska(ixO^S)= ((xshift(ixO^S,1)-xc1)**2+(xshift(ixO^S,2)-yc1)**2 <= rad**2) &
            .and.(xshift(ixO^S,2)<= (ya2-ya1)*(xshift(ixO^S,1)-xa1)/(xa2-xa1)+ya1) & 
            .and.(xshift(ixO^S,2)<= (yb2-yb1)*(xshift(ixO^S,1)-xb1)/(xb2-xb1)+yb1) & 
            .and.((xshift(ixO^S,1)-xc2)**2+(xshift(ixO^S,2)-yc2)**2 > rad2**2) 
     xc1=0.75d0
     yc1=0.50d0
     alp=half*dpi-alp
     xa1=xc1-rad
     ya1=yc1
     xa2=xc1+rad*cos(alp)
     ya2=yc1+rad*sin(alp)
     xb1=xa1
     yb1=ya1
     xb2=xc1+rad*cos(alp)
     yb2=yc1-rad*sin(alp)
     yc2=yc1
     xc2=xa2+sqrt(rad2**2-(ya2-yc2)**2)
     maskc(ixO^S)= ((xshift(ixO^S,1)-xc1)**2+(xshift(ixO^S,2)-yc1)**2 <= rad**2) &
            .and.(xshift(ixO^S,2)<= (ya2-ya1)*(xshift(ixO^S,1)-xa1)/(xa2-xa1)+ya1) & 
            .and.(xshift(ixO^S,2)>= (yb2-yb1)*(xshift(ixO^S,1)-xb1)/(xb2-xb1)+yb1) & 
            .and.((xshift(ixO^S,1)-xc2)**2+(xshift(ixO^S,2)-yc2)**2 > rad2**2) 
     where(maskv(ixO^S).or.maska(ixO^S).or.maskc(ixO^S))
        rhoprofile(ixO^S)     = rhosquare
     elsewhere
        rhoprofile(ixO^S)     = rhoflat
     endwhere
     }

  end subroutine set_density_profile

  subroutine store_sol_err(igrid,level,ixI^L,ixO^L,qt,w,x)
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: rhoprofile(ixI^S)

    call set_density_profile(ixI^L,ixO^L,qt,x,rhoprofile)
    w(ixO^S,i_sol) = rhoprofile(ixO^S)
    w(ixO^S,i_err) = dabs(w(ixO^S,rho_) - w(ixO^S,i_sol))
  end subroutine store_sol_err

  subroutine print_min_max
    use mod_input_output, only: get_global_minima, get_global_maxima, &
                                get_volume_average
    double precision   :: minvals(nw),maxvals(nw)

    integer :: iw
    double precision :: modes(nw,2), volume

    call get_global_minima(minvals)
    call get_global_maxima(maxvals)
    call get_volume_average(1,modes(:,1),volume)
    call get_volume_average(2,modes(:,2),volume)

    if (mype == 0) then
       write(*, "(A,4E16.8,A,3E16.8)") " time + rho min-max-tot:", &
              global_time, minvals(rho_), maxvals(rho_), modes(rho_,1), &
             " Error: Linf-L1-L2:", &
              maxvals(i_err),modes(i_err,1),dsqrt(modes(i_err,2))
    end if
  end subroutine print_min_max

  subroutine generate_particles(n_particles, x, v, q, m, follow)
    use mod_particles
    integer, intent(in)           :: n_particles
    double precision, intent(out) :: x(3, n_particles)
    double precision, intent(out) :: v(3, n_particles)
    double precision, intent(out) :: q(n_particles)
    double precision, intent(out) :: m(n_particles)
    logical, intent(out)          :: follow(n_particles)
    integer                       :: n

    do n = 1, n_particles
      call get_particle(x(:, n), v(:, n), q(n), m(n), n, n_particles)
    end do

    follow(:) = .true.

  end subroutine generate_particles

  ! Set particle properties (in code units)
  subroutine get_particle(x, v, q, m, ipart, n_particles)
    double precision, intent(out) :: x(3), v(3), q, m
    integer, intent(in)           :: ipart, n_particles

    q = 0.d0
    m = 0.d0
    v(:) = 0.d0

    ! 3 particles in a line
    x(1) = xprobmin1 + (xprobmax1 - xprobmin1)/4.d0*DBLE(ipart)
    x(2) = xprobmin2 + (xprobmax2 - xprobmin2)/2.d0
    print *, "Particle", ipart, x

  end subroutine get_particle

  ! User-defined particle movement
  subroutine move_particle(x, n, told, tnew)
    double precision, intent(inout) :: x(3)
    double precision, intent(in)  :: told,tnew
    integer, intent(in)           :: n
    double precision              :: v_move
    double precision              :: rad, th, xc,yc, freq

    v_move = 0.25d0
    rad = 0.0398d0
    xc = (xprobmax1 - xprobmin1)*.75d0-rad
    yc = (xprobmax2 - xprobmin2)/2.d0
    freq = v_move/rad

    ! particle 1 does not move;
    ! particle 2 moves on a straight line with vy=v_move
    ! particle 3 moves in a circle with vtang=v_move
    if (n .eq. 2) x(2) = x(2) + v_move*(tnew-told)
    if (n .eq. 3) then
      th = atan2(x(2)-yc,x(1)-xc)
      th = th + freq*(tnew-told)
      x(1) = xc+rad*cos(th)
      x(2) = yc+rad*sin(th)
    end if

  end subroutine move_particle

  ! Custom particle gridvars definition
  subroutine define_additional_gridvars_usr(ngridvars)
    use mod_global_parameters
    integer, intent(inout) :: ngridvars
 
    ! extra variable as payload: add the actual solution
    isol = ngridvars+1
    ngridvars = ngridvars+1

  end subroutine define_additional_gridvars_usr

  ! Custom particle gridvars filler
  subroutine fill_additional_gridvars_usr
    use mod_global_parameters
    use mod_usr_methods, only: usr_particle_fields

    !integer :: igrid, iigrid
    !double precision :: rhoprofile(ixG^T)

    !do iigrid=1,igridstail; igrid=igrids(iigrid);
    !  call set_density_profile(ixG^LL,ixG^LL,global_time,ps(igrid)%x(ixG^T,1:ndim),rhoprofile)
    !  gridvars(igrid)%w(ixG^T,isol) = rhoprofile(ixG^T)
    !end do

  end subroutine fill_additional_gridvars_usr

  subroutine update_payload_usr(igrid,w,wold,xgrid,xpart,upart,qpart,mpart,payload,npayload,particle_time)
    use mod_global_parameters
    integer, intent(in)           :: igrid,npayload
    double precision, intent(in)  :: w(ixG^T,1:nw),wold(ixG^T,1:nw)
    double precision, intent(in)  :: xgrid(ixG^T,1:ndim),xpart(1:ndir),upart(1:ndir),qpart,mpart,particle_time
    double precision, intent(out) :: payload(npayload)
 
    double precision :: rhoprofile(ixG^T)

    ! put the solution at particle_time for comparison
    if (npayload > 0) then
      call set_density_profile(ixG^LL,ixG^LL,particle_time,xgrid,rhoprofile)
      gridvars(igrid)%w(ixG^T,isol) = rhoprofile(ixG^T)
      call interpolate_var(igrid,ixG^LL,ixG^LL,gridvars(igrid)%w(ixG^T,isol),xgrid,xpart,payload(1))
    end if

  end subroutine update_payload_usr

end module mod_usr
