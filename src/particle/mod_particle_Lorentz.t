module mod_particle_Lorentz
  use mod_particle_base

  private

  !> Variable index for magnetic field
  integer, allocatable, protected :: bp(:)
  !> Variable index for electric field
  integer, allocatable, protected :: ep(:)

  public :: Lorentz_init
  public :: Lorentz_create_particles
  public :: bp, ep

contains

  subroutine Lorentz_init()
    use mod_global_parameters
    integer :: idir, nwx

    if(physics_type/='mhd') call mpistop("Lorentz particles need magnetic field!")
    if(ndir/=3) call mpistop("Lorentz particles need ndir=3!")
    dtsave_particles=dtsave_particles*unit_time
    ngridvars=ndir*2
    nwx = 0

    allocate(bp(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      bp(idir) = nwx
    end do

    allocate(ep(ndir))
    do idir = 1, ndir
      nwx = nwx + 1
      ep(idir) = nwx
    end do

    particles_fill_gridvars => Lorentz_fill_gridvars
    particles_integrate     => Lorentz_integrate_particles
  end subroutine Lorentz_init

  subroutine Lorentz_create_particles()

    use mod_global_parameters
    use mod_usr_methods, only: usr_create_particles

    integer          :: n, igrid_particle, ipe_particle
    double precision :: lfac
    double precision :: x(3, num_particles)
    double precision :: v(3, num_particles)
    double precision :: q(num_particles)
    double precision :: m(num_particles)
    logical          :: follow(num_particles)

    if (.not. associated(usr_create_particles)) then
      call mpistop("Error: no usr_create_particles method specified")
    else if (mype == 0) then
      call usr_create_particles(num_particles, x, v, q, m, follow)
    end if

    call MPI_BCAST(x,3*num_particles,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(v,3*num_particles,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(q,num_particles,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(m,num_particles,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    call MPI_BCAST(follow,num_particles,MPI_LOGICAL,0,icomm,ierrmpi)

    nparticles = num_particles

    ! Find ipe and igrid responsible for particle
    do n = 1, num_particles
      call find_particle_ipe(x(:, n),igrid_particle,ipe_particle)

      particle(n)%igrid = igrid_particle
      particle(n)%ipe   = ipe_particle

      if (ipe_particle == mype) then
        call push_particle_into_particles_on_mype(n)

        call get_lfac_from_velocity(v(:, n), lfac)

        allocate(particle(n)%self)
        particle(n)%self%x      = x(:, n)
        particle(n)%self%u      = v(:, n) * lfac
        particle(n)%self%q      = q(n)
        particle(n)%self%m      = m(n)
        particle(n)%self%follow = follow(n)
        particle(n)%self%index  = n
        particle(n)%self%t      = 0.0d0
        particle(n)%self%dt     = 0.0d0

        ! initialise payloads for Lorentz module
        allocate(particle(n)%payload(npayload))
        particle(n)%payload(:) = 0.0d0
      end if
    end do

  end subroutine Lorentz_create_particles

  subroutine Lorentz_fill_gridvars()
    use mod_global_parameters
    use mod_usr_methods, only: usr_particle_fields

    integer                                   :: igrid, iigrid, idir
    double precision, dimension(ixG^T,1:ndir) :: beta
    double precision, dimension(ixG^T,1:nw)   :: w,wold
    double precision                          :: current(ixG^T,7-2*ndir:3)
    integer                                   :: idirmin
    double precision                          :: E_field(ixG^T, 1:ndir)
    double precision                          :: B_field(ixG^T, 1:ndir)


    do iigrid=1,igridstail; igrid=igrids(iigrid);
      gridvars(igrid)%w(ixG^T,1:ngridvars) = 0.0d0

      if (associated(usr_particle_fields)) then
        call usr_particle_fields(pw(igrid)%x, E_field, B_field)
        gridvars(igrid)%w(ixG^T,ep(:)) = E_field
        gridvars(igrid)%w(ixG^T,bp(:)) = B_field
      else
        ! Determine fields from MHD variables
        ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
        w(ixG^T,1:nw) = pw(igrid)%w(ixG^T,1:nw)
        call phys_to_primitive(ixG^LL,ixG^LL,w,pw(igrid)%x)

        ! fill with magnetic field:
        gridvars(igrid)%w(ixG^T,bp(:)) = w(ixG^T,iw_mag(:))

        ! fill with electric field
        current = zero
        call get_current(w,ixG^LL,ixG^LLIM^D^LSUB1,idirmin,current)
        gridvars(igrid)%w(ixG^T,ep(1)) = gridvars(igrid)%w(ixG^T,bp(2)) * w(ixG^T,iw_mom(3)) &
             - gridvars(igrid)%w(ixG^T,bp(3)) * w(ixG^T,iw_mom(2)) + particles_eta * current(ixG^T,1)
        gridvars(igrid)%w(ixG^T,ep(2)) = gridvars(igrid)%w(ixG^T,bp(3)) * w(ixG^T,iw_mom(1)) &
             - gridvars(igrid)%w(ixG^T,bp(1)) * w(ixG^T,iw_mom(3)) + particles_eta * current(ixG^T,2)
        gridvars(igrid)%w(ixG^T,ep(3)) = gridvars(igrid)%w(ixG^T,bp(1)) * w(ixG^T,iw_mom(2)) &
             - gridvars(igrid)%w(ixG^T,bp(2)) * w(ixG^T,iw_mom(1)) + particles_eta * current(ixG^T,3)

        ! scale to cgs units:
        gridvars(igrid)%w(ixG^T,bp(:)) = &
             gridvars(igrid)%w(ixG^T,bp(:)) * sqrt(4.0d0*dpi*unit_velocity**2 * unit_density)
        gridvars(igrid)%w(ixG^T,ep(:)) = &
             gridvars(igrid)%w(ixG^T,ep(:)) * sqrt(4.0d0*dpi*unit_velocity**2 * unit_density) * unit_velocity / const_c
      end if

      if (time_advance) then
        gridvars(igrid)%wold(ixG^T,1:ngridvars) = 0.0d0

        if (associated(usr_particle_fields)) then
          !> @todo Compute these fields at a different time?
          call usr_particle_fields(pw(igrid)%x, E_field, B_field)
          gridvars(igrid)%wold(ixG^T,ep(:)) = E_field
          gridvars(igrid)%wold(ixG^T,bp(:)) = B_field
        else
          wold(ixG^T,1:nw) = pw(igrid)%wold(ixG^T,1:nw)
          call phys_to_primitive(ixG^LL,ixG^LL,wold,pw(igrid)%x)
          ! fill with magnetic field:
          gridvars(igrid)%wold(ixG^T,bp(:)) = wold(ixG^T,iw_mag(:))
          ! fill with electric field
          current = zero
          call get_current(wold,ixG^LL,ixG^LLIM^D^LSUB1,idirmin,current)
          gridvars(igrid)%wold(ixG^T,ep(1)) = gridvars(igrid)%wold(ixG^T,bp(2)) * wold(ixG^T,iw_mom(3)) &
               - gridvars(igrid)%wold(ixG^T,bp(3)) * wold(ixG^T,iw_mom(2)) + particles_eta * current(ixG^T,1)
          gridvars(igrid)%wold(ixG^T,ep(2)) = gridvars(igrid)%wold(ixG^T,bp(3)) * wold(ixG^T,iw_mom(1)) &
               - gridvars(igrid)%wold(ixG^T,bp(1)) * wold(ixG^T,iw_mom(3)) + particles_eta * current(ixG^T,2)
          gridvars(igrid)%wold(ixG^T,ep(3)) = gridvars(igrid)%wold(ixG^T,bp(1)) * wold(ixG^T,iw_mom(2)) &
               - gridvars(igrid)%wold(ixG^T,bp(2)) * wold(ixG^T,iw_mom(1)) + particles_eta * current(ixG^T,3)

          ! scale to cgs units:
          gridvars(igrid)%wold(ixG^T,bp(:)) = &
               gridvars(igrid)%wold(ixG^T,bp(:)) * sqrt(4.0d0*dpi*unit_velocity**2 * unit_density)
          gridvars(igrid)%wold(ixG^T,ep(:)) = &
               gridvars(igrid)%wold(ixG^T,ep(:)) * sqrt(4.0d0*dpi*unit_velocity**2 * unit_density) * unit_velocity / const_c
        end if
      end if

    end do

  end subroutine Lorentz_fill_gridvars

  !> Relativistic Boris scheme
  subroutine Lorentz_integrate_particles(end_time)
    use mod_global_parameters
    double precision, intent(in)      :: end_time
    integer                           :: ipart, iipart
    double precision                  :: lfac, q, m, dt_p, cosphi, sinphi, phi1, phi2, r, re
    double precision, dimension(ndir) :: b, e, emom, uminus, t_geom, s, udash, tmp, uplus, xcart1, xcart2, ucart2, radmom

    do iipart=1,nparticles_active_on_mype
      ipart = particles_active_on_mype(iipart);
      q     = particle(ipart)%self%q
      m     = particle(ipart)%self%m

      dt_p = Lorentz_get_particle_dt(particle(ipart), end_time)
      particle(ipart)%self%dt = dt_p

      ! Push particle over half time step
      call get_lfac(particle(ipart)%self%u,lfac)
      particle(ipart)%self%x(1:ndir) = particle(ipart)%self%x(1:ndir) &
           + 0.5d0 * dt_p * particle(ipart)%self%u(1:ndir)/lfac &
           * const_c / unit_length

      ! Get E, B at new position
      call get_vec(bp, particle(ipart)%igrid, &
           particle(ipart)%self%x,particle(ipart)%self%t,b)
      call get_vec(ep, particle(ipart)%igrid, &
           particle(ipart)%self%x,particle(ipart)%self%t,e)

      ! 'Kick' particle (update velocity)
      select case(typeaxial)

        ! CARTESIAN COORDINATES
      case('slab')
        ! Momentum update
        emom = q * e * dt_p /(2.0d0 * m * const_c)
        if(losses) then
          call get_lfac(particle(ipart)%self%u,lfac)
          re = abs(q)**2 / (m * const_c**2)
          call cross(particle(ipart)%self%u,b,tmp)
          radmom = - third * re**2 * lfac &
               * ( sum((e(:)+tmp(:)/lfac)**2)  &
               -  (sum(e(:)*particle(ipart)%self%u(:))/lfac)**2 ) &
               * particle(ipart)%self%u / m / const_c * dt_p
        else
          radmom = 0.0d0
        end if

        uminus = particle(ipart)%self%u + emom + radmom

        call get_lfac(uminus,lfac)
        call get_t(b,lfac,dt_p,q,m,t_geom)
        call get_s(t_geom,s)

        call cross(uminus,t_geom,tmp)
        udash = uminus + tmp

        call cross(udash,s,tmp)
        uplus = uminus + tmp

        if(losses) then
          call cross(uplus,b,tmp)
          radmom = - third * re**2 * lfac &
               * ( sum((e(:)+tmp(:)/lfac)**2)  &
               -  (sum(e(:)*uplus(:))/lfac)**2 ) &
               * uplus / m / const_c * dt_p
        else
          radmom = 0.0d0
        end if

        ! Update the velocity
        particle(ipart)%self%u = uplus + emom + radmom

        ! Payload update
        ! current gyroradius
        call cross(particle(ipart)%self%u,b,tmp)
        tmp = tmp / sqrt(sum(b(:)**2))
        particle(ipart)%payload(1) = sqrt(sum(tmp(:)**2)) / sqrt(sum(b(:)**2)) * &
             m / abs(q) * const_c**2

        ! e.b
        if(npayload>1) particle(ipart)%payload(2) = sum(e(:)*b(:))/sqrt(sum(b(:)**2)*sum(e(:)**2))

      case ('cylindrical')

        !  Momentum update
        emom = q * e * dt_p /(2.0d0 * m * const_c)

        if(losses) then
          call get_lfac(particle(ipart)%self%u,lfac)
          re = abs(q)**2 / (m * const_c**2)
          call cross(particle(ipart)%self%u,b,tmp)
          radmom = - third * re**2 * lfac &
               * ( sum((e(:)+tmp(:)/lfac)**2)  &
               -  (sum(e(:)*particle(ipart)%self%u(:))/lfac)**2 ) &
               * particle(ipart)%self%u / m / const_c * dt_p
        else
          radmom = 0.0d0
        end if

        uminus = particle(ipart)%self%u + emom + radmom

        call get_lfac(uminus,lfac)
        call get_t(b,lfac,dt_p,q,m,t_geom)
        call get_s(t_geom,s)

        call cross(uminus,t_geom,tmp)
        udash = uminus + tmp

        call cross(udash,s,tmp)
        uplus = uminus + tmp

        if(losses) then
          call cross(uplus,b,tmp)
          radmom = - third * re**2 * lfac &
               * ( sum((e(:)+tmp(:)/lfac)**2)  &
               -  (sum(e(:)*uplus(:))/lfac)**2 ) &
               * uplus / m / const_c * dt_p
        else
          radmom = 0.0d0
        end if

        particle(ipart)%self%u = uplus + emom + radmom
        ! Position update
        ! Get cartesian coordinates:
        phi1       = particle(ipart)%self%x(phi_)
        cosphi     = cos(phi1)
        sinphi     = sin(phi1)

        xcart1(1)  = particle(ipart)%self%x(r_) * cosphi
        xcart1(2)  = particle(ipart)%self%x(r_) * sinphi
        xcart1(3)  = particle(ipart)%self%x(z_)

        ucart2(1)   = cosphi * particle(ipart)%self%u(r_) - sinphi * particle(ipart)%self%u(phi_)
        ucart2(2)   = cosphi * particle(ipart)%self%u(phi_) + sinphi * particle(ipart)%self%u(r_)
        ucart2(3)   = particle(ipart)%self%u(z_)

        ! update position
        xcart2(1:ndir) = xcart1(1:ndir) &
             + dt_p * ucart2(1:ndir)/lfac &
             * const_c/unit_length

        ! back to cylindrical coordinates
        phi2 = atan2(xcart2(2),xcart2(1))
        if(phi2 .lt. 0.0d0) phi2 = 2.0d0*dpi + phi2
        r = sqrt(xcart2(1)**2 + xcart2(2)**2)
        particle(ipart)%self%x(r_)   = r
        particle(ipart)%self%x(phi_) = phi2
        particle(ipart)%self%x(z_)   = xcart2(3)

        ! Rotate the momentum to the new cooridnates
        ! rotate velocities
        cosphi     = cos(phi2-phi1)
        sinphi     = sin(phi2-phi1)

        tmp = particle(ipart)%self%u
        particle(ipart)%self%u(r_)   = cosphi * tmp(r_)   + sinphi * tmp(phi_)
        particle(ipart)%self%u(phi_) = cosphi * tmp(phi_) - sinphi * tmp(r_)
        particle(ipart)%self%u(z_)   = tmp(z_)

        ! Payload update
        ! current gyroradius
        call cross(particle(ipart)%self%u,b,tmp)
        tmp = tmp / sqrt(sum(b(:)**2))
        particle(ipart)%payload(1) = sqrt(sum(tmp(:)**2)) / sqrt(sum(b(:)**2)) * m / abs(q) * 8.9875d+20

      end select

      call get_lfac(particle(ipart)%self%u,lfac)

      ! Push particle over half time step at end
      particle(ipart)%self%x(1:ndir) = particle(ipart)%self%x(1:ndir) &
           + 0.5d0 * dt_p * particle(ipart)%self%u(1:ndir)/lfac &
           * const_c / unit_length

      ! Time update
      particle(ipart)%self%t = particle(ipart)%self%t + dt_p

    end do ! ipart loop

  contains

    subroutine get_t(b,lfac,dt,q,m,t)
      implicit none
      double precision, dimension(ndir), intent(in)      :: b
      double precision, intent(in)                      :: lfac, dt, q, m
      double precision, dimension(ndir), intent(out)     :: t

      t = q * b * dt / (2.0d0 * lfac * m * const_c)

    end subroutine get_t

    subroutine get_s(t,s)
      implicit none
      double precision, dimension(ndir), intent(in)   :: t
      double precision, dimension(ndir), intent(out)  :: s

      s = 2.0d0 * t / (1.0d0+sum(t(:)**2))

    end subroutine get_s

  end subroutine Lorentz_integrate_particles

  function Lorentz_get_particle_dt(partp, end_time) result(dt_p)
    use mod_global_parameters
    type(particle_ptr), intent(in)   :: partp
    double precision, intent(in)     :: end_time
    double precision                 :: dt_p
    integer                          :: ipart, iipart, nout
    double precision,dimension(ndir) :: b,v
    double precision                 :: lfac,absb,dt_cfl
    double precision                 :: tout
    double precision, parameter      :: cfl = 0.5d0

    if (const_dt_particles > 0) then
      dt_p = const_dt_particles
      return
    end if

    call get_vec(bp, partp%igrid,partp%self%x,partp%self%t,b)
    absb = sqrt(sum(b(:)**2))
    call get_lfac(partp%self%u,lfac)

    ! CFL timestep
    ! make sure we step only one cell at a time:
    v(:) = abs(const_c * partp%self%u(:) / lfac)

    ! convert to angular velocity:
    if(typeaxial =='cylindrical'.and.phi_>0) then
      v(phi_) = abs(v(phi_)/partp%self%x(r_))
    end if

    dt_cfl = min({rnode(rpdx^D_,partp%igrid)/max(v(^D), smalldouble)})

    if(typeaxial =='cylindrical'.and.phi_>0) then
      ! phi-momentum leads to radial velocity:
      if(phi_ .gt. ndim) dt_cfl = min(dt_cfl, &
           sqrt(rnode(rpdx1_,partp%igrid)/partp%self%x(r_)) &
           / v(phi_))
      ! limit the delta phi of the orbit (just for aesthetic reasons):
      dt_cfl = min(dt_cfl,0.1d0/max(v(phi_), smalldouble))
      ! take some care at the axis:
      dt_cfl = min(dt_cfl,(partp%self%x(r_)+smalldouble)/max(v(r_), smalldouble))
    end if

    dt_cfl = dt_cfl * cfl

    ! bound by gyro-rotation:
    dt_p = abs( dtheta * const_c/unit_length * partp%self%m * lfac &
         / (partp%self%q * absb) )

    dt_p = min(dt_p, dt_cfl)*unit_length

    ! Make sure we don't advance beyond end_time
    call limit_dt_endtime(end_time - partp%self%t, dt_p)

  end function Lorentz_get_particle_dt

end module mod_particle_Lorentz
