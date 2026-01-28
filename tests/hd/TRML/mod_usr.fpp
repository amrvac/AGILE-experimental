module mod_usr

  use mod_amrvac
  use mod_physics
! TODO: Be more explicit about what you're using because I don't know what is
!   coming from where.
  use mod_global_parameters
  use mod_constants
  implicit none

  double precision, parameter :: pi = 4d0*atan(1d0)

! User parameters.
! NOTE: These are in CGS units.
  integer :: seed
  double precision :: xi
  double precision :: chi
  double precision :: mach_H
  double precision :: rho_C
  double precision :: L_x
  double precision :: v_sh
  double precision :: v_off
  double precision, dimension(2:3) :: v_per
  double precision :: z0
  double precision :: d_z
  double precision :: sig_z
  double precision :: s
  double precision, dimension(2:3) :: f_til
  integer :: mode_root
  integer :: n_modes

! Precomputed derived quantities.
! NOTE: These are in CGS units.
  double precision, dimension(3) :: L
  double precision :: rho_H
  double precision :: P
  double precision :: v_C
  double precision :: v_H
  double precision :: t_sh
  double precision :: c_H
  double precision :: t_cool
  double precision :: Lam_peak
  double precision :: T_H
  double precision :: T_C
  double precision :: T_mix
  double precision :: rho_mix
  double precision :: b

! Precomputed random numbers.
  double precision, dimension(:,:,:), allocatable :: rand_1

contains

  subroutine params_read_usr(files)
    character(len=*), dimension(:), intent(in) :: files
    integer :: n
    namelist /usr_list/&
      seed, xi, chi, mach_H, rho_C, L_x, v_sh, v_off, v_per, z0, d_z, sig_z, s,&
      f_til, mode_root, n_modes
    do n = 1, size(files)
      open(unitpar, file=trim(files(n)), status='old')
      read(unitpar, usr_list, end=111)
111   close(unitpar)
    end do
  end subroutine params_read_usr

  subroutine usr_init()
    call set_coordinate_system('Cartesian_3D')
    call params_read_usr(par_files)
! mp_cgs = 1.67e-24
! kb_cgs = 1.38e-16
    ! unit_temperature = mp_cgs/kb_cgs/b ! 6.05e-9 for b = 2
    ! unit_numberdensity = 1d0/mp_cgs ! 5.99e23
    unit_length = 7.17d5 ! L_x
    unit_density = 1d-10 ! rho_C
    unit_temperature = 2.00d4 ! T_C
    usr_set_parameters => set_parameters_usr
    usr_init_one_grid  => init_one_grid_usr
    usr_special_bc     => specialbound_usr
    call phys_activate()
    if (mype == 0) call print_units()
    ! call mpistop('usr_init successful')
  end subroutine usr_init

! TODO: Nicely print the user parameters, good ordering and formatting.
  subroutine print_params()
    print *, 'xi', xi
    print *, 'chi', chi
    print *, 'mach_H', mach_H
    print *, 'rho_C', rho_C
    print *, 'L_x', L_x
    print *, 't_sh', t_sh
    print *, 'v_off', v_off
    print *, 'v_per(2)', v_per(2)
    print *, 'v_per(3)', v_per(3)
    print *, 'z0', z0
    print *, 'd_z', d_z
    print *, 'sig_z', sig_z
    print *, 's', s
    print *, 'f_til(2)', f_til(2)
    print *, 'f_til(3)', f_til(3)
    print *, 'mode_root', mode_root
    print *, 'n_modes', n_modes
    print *, 'L(1)', L(1)
    print *, 'L(2)', L(2)
    print *, 'L(3)', L(3)
    print *, 'rho_H', rho_H
    print *, 'P', P
    print *, 'v_C', v_C
    print *, 'v_H', v_H
    print *, 'v_sh', v_sh
    print *, 'c_H', c_H
    print *, 't_cool', t_cool
    print *, 'Lam_peak', Lam_peak
    print *, 'T_H', T_H
    print *, 'T_C', T_C
    print *, 'T_mix', T_mix
    print *, 'rho_mix', rho_mix
  end subroutine

! TODO: Remove this, it's for debugging purposes only.
  subroutine print_units()
    print *, 'print_units'
    print *, 'unit_length', unit_length
    print *, 'unit_time', unit_time
    print *, 'unit_velocity', unit_velocity
    print *, 'unit_density', unit_density
    print *, 'unit_numberdensity', unit_numberdensity
    print *, 'unit_pressure', unit_pressure
    print *, 'unit_temperature', unit_temperature
  end subroutine

! TODO: Stuff is not OpenACC updated, may not work on GPU, see other examples.
!   Every variable needs to be copied to the GPU, they are not passed through
!   the kernels since Fortran promotes the usage of globals.
  subroutine set_parameters_usr()
    integer :: i
    integer :: s_r_state
    integer, dimension(:), allocatable :: r_state

! `rand_1` is random uniform in [-1, 1).
! TODO: I'm not sure if this even works.
    call random_seed(size=s_r_state)
    allocate(r_state(s_r_state))
    r_state = seed
    call random_seed(put=r_state)
    allocate(rand_1(n_modes,3,2:3))
    do i = 1, 10
      call random_number(rand_1(1,1,1))
    end do
    call random_number(rand_1)
    rand_1 = rand_1*2-1

    L = [xprobmax1-xprobmin1,&
         xprobmax2-xprobmin2,&
         xprobmax3-xprobmin3]
    L = L*L_x

! TODO: I have been assuming b=1, need to redo some of this. But it's not a big
!   deal, shouldn't make or break it.
    t_sh = L(1)/v_sh
    v_C = -v_sh/2+v_off
    v_H = +v_sh/2+v_off
    c_H = v_sh/mach_H
    t_cool = t_sh/xi
    rho_H = rho_C/chi
    rho_mix = sqrt(rho_C*rho_H)
    P = c_H**2*rho_H/hd_gamma
    T_H = mp_cgs/(hd_gamma*kb_cgs)*c_H**2
    T_C = T_H/chi
    T_mix = sqrt(T_C*T_H)
    Lam_peak = hd_gamma/(hd_gamma-1)*mp_cgs**2*P/(rho_C*rho_H*t_cool)

    if (mype == 0) call print_params()
  end subroutine set_parameters_usr

  subroutine init_one_grid_usr(&
      ixGmin1, ixGmin2, ixGmin3, ixGmax1, ixGmax2, ixGmax3,&
      ixmin1,  ixmin2,  ixmin3,  ixmax1,  ixmax2,  ixmax3,&
      w, x)
    integer, intent(in) :: ixGmin1, ixGmin2, ixGmin3, ixGmax1, ixGmax2, ixGmax3
    integer, intent(in) :: ixmin1,  ixmin2,  ixmin3,  ixmax1,  ixmax2,  ixmax3
    double precision, intent(inout),&
      dimension(ixGmin1:ixGmax1, ixGmin2:ixGmax2, ixGmin3:ixGmax3, 1:nw) :: w
    double precision, intent(in),&
      dimension(ixGmin1:ixGmax1, ixGmin2:ixGmax2, ixGmin3:ixGmax3, 1:ndim) :: x
    double precision, dimension(3) :: c
    integer :: n

    associate(&
      w_ => w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,:),&
      x_ => x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,:))

! STRATEGY:
!   We are given empty w_, and x_ in code units, and w_ should be returned as
!   conserved.
!   1. Compute w_ in CGS units. Take care, x_ is in code units.
!   2. Convert w_ to code units.
!   3. Convert w_ to conserved.

      w_(:,:,:,iw_rho) = rho_H+(rho_C-rho_H)/2*(1-tanh((x_(:,:,:,3)*unit_length-z0)/d_z))
      w_(:,:,:,iw_e) = P
      w_(:,:,:,iw_mom(1)) = v_H-v_sh/2*(1-tanh((x_(:,:,:,3)*unit_length-z0)/d_z))
! We first compute the the noise terms in the y- and z-velocity buffers, then
! the actual velocities from those.
      w_(:,:,:,iw_mom(2)) = 0
      w_(:,:,:,iw_mom(3)) = 0
      do n = mode_root+1, mode_root-1+n_modes
        c(:) = rand_1(n-mode_root,:,2)
        w_(:,:,:,iw_mom(2)) = w_(:,:,:,iw_mom(2))+f_til(2)*c(1)*(&
          sin(2*pi*x_(:,:,:,1)*unit_length/L(1)*n+f_til(2)*c(2))/2+&
          sin(2*pi*x_(:,:,:,2)*unit_length/L(2)*n+f_til(2)*c(3))/2)
      end do
      do n = mode_root+1, mode_root-1+n_modes
        c(:) = rand_1(n-mode_root,:,3)
        w_(:,:,:,iw_mom(3)) = w_(:,:,:,iw_mom(3))+f_til(3)*c(1)*(&
          sin(2*pi*x_(:,:,:,1)*unit_length/L(1)*n+f_til(3)*c(2))/2+&
          sin(2*pi*x_(:,:,:,2)*unit_length/L(2)*n+f_til(3)*c(3))/2)
      end do
      w_(:,:,:,iw_mom(2)) = v_per(2)*(&
        sin(2*pi*x_(:,:,:,1)*unit_length/L(1)*mode_root)/2+&
        sin(2*pi*x_(:,:,:,2)*unit_length/L(2)*mode_root)/2+&
        w_(:,:,:,iw_mom(2)))*&
        exp(-((x_(:,:,:,3)*unit_length-z0)/sig_z)**2)
      w_(:,:,:,iw_mom(3)) = v_per(3)*(&
        sin(2*pi*x_(:,:,:,1)*unit_length/L(1)*mode_root/2)+&
        sin(2*pi*x_(:,:,:,2)*unit_length/L(2)*mode_root/2)+&
        w_(:,:,:,iw_mom(3)))*&
        exp(-((x_(:,:,:,3)*unit_length-z0)/sig_z)**2)

      w_(:,:,:,iw_rho)    = w_(:,:,:,iw_rho)/unit_density
      w_(:,:,:,iw_e)      = w_(:,:,:,iw_e)/unit_pressure
      w_(:,:,:,iw_mom(1)) = w_(:,:,:,iw_mom(1))/unit_velocity
      w_(:,:,:,iw_mom(2)) = w_(:,:,:,iw_mom(2))/unit_velocity
      w_(:,:,:,iw_mom(3)) = w_(:,:,:,iw_mom(3))/unit_velocity

    end associate

    call phys_to_conserved(&
      ixGmin1, ixGmin2, ixGmin3, ixGmax1, ixGmax2, ixGmax3,&
      ixmin1,  ixmin2,  ixmin3,  ixmax1,  ixmax2,  ixmax3,&
      w, x)
  end subroutine init_one_grid_usr

! NOTE: This must be named exactly `specialbound_usr` in AGILE since during
!   compilation it deals with this function.
  subroutine specialbound_usr(&
      qt,&
      ixImin1, ixImin2, ixImin3, ixImax1, ixImax2, ixImax3,&
      ixOmin1, ixOmin2, ixOmin3, ixOmax1, ixOmax2, ixOmax3,&
      iB, w, x)
!$acc routine vector
    integer, intent(in) :: ixImin1, ixImin2, ixImin3, ixImax1, ixImax2, ixImax3
    integer, intent(in) :: ixOmin1, ixOmin2, ixOmin3, ixOmax1, ixOmax2, ixOmax3
    integer, intent(in) :: iB
    double precision, intent(in) :: qt
    double precision, intent(in),&
      dimension(ixImin1:ixImax1, ixImin2:ixImax2, ixImin3:ixImax3, 1:ndim) :: x
    double precision, intent(inout),&
      dimension(ixImin1:ixImax1, ixImin2:ixImax2, ixImin3:ixImax3, 1:nw) :: w
    integer :: ix1, ix2, ix3

! NOTE:
!   Boundary conditions work on conserved variables. For simplicity I convert
!   to primitive, but this is less performant, Hao Wu writes this in primitive.
!   So this should be updated.

    call phys_to_primitive(&
      ixImin1, ixImin2, ixImin3, ixImax1, ixImax2, ixImax3,&
      ixOmin1, ixOmin2, ixOmin3, ixOmax1, ixOmax2, ixOmax3,&
      w, x)

! NOTE: Zero-gradient for vy and vz, i.e. copy the closest interior cell into
!   the halo. Corners don't matter.
    select case(iB)
    case(5)
!$acc loop collapse(3) vector
      do ix3 = ixOmin3, ixOmax3
        do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1
            w(ix1,ix2,ix3,iw_rho) = rho_C/unit_density
            w(ix1,ix2,ix3,iw_e) = P/unit_pressure
            w(ix1,ix2,ix3,iw_mom(1)) = v_C/unit_velocity
            w(ix1,ix2,ix3,iw_mom(2)) = w(ix1,ix2,ixOmax3+1,iw_mom(2))
            w(ix1,ix2,ix3,iw_mom(3)) = w(ix1,ix2,ixOmax3+1,iw_mom(3))
          end do
        end do
      end do
    case(6)
!$acc loop collapse(3) vector
      do ix3 = ixOmin3, ixOmax3
        do ix2 = ixOmin2, ixOmax2
          do ix1 = ixOmin1, ixOmax1
            w(ix1,ix2,ix3,iw_rho) = rho_H/unit_density
            w(ix1,ix2,ix3,iw_e) = P/unit_pressure
            w(ix1,ix2,ix3,iw_mom(1)) = v_H/unit_velocity
            w(ix1,ix2,ix3,iw_mom(2)) = w(ix1,ix2,ixOmin3-1,iw_mom(2))
            w(ix1,ix2,ix3,iw_mom(3)) = w(ix1,ix2,ixOmin3-1,iw_mom(3))
          end do
        end do
      end do
    end select

    call phys_to_conserved(&
      ixImin1, ixImin2, ixImin3, ixImax1, ixImax2, ixImax3,&
      ixOmin1, ixOmin2, ixOmin3, ixOmax1, ixOmax2, ixOmax3,&
      w, x)
  end subroutine specialbound_usr

end module mod_usr
