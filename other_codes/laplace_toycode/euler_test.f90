#include "macros.h"
! Author: Jannis Teunissen (jannis.teunissen@cwi.nl)
program euler_test
  use m_blockgrid
  use m_euler
  use m_config

  implicit none
  integer, parameter :: dp = kind(0.0d0)

  type(block_grid_t) :: bg
  type(CFG_t)        :: cfg

  integer            :: n_iter   = 100
  logical            :: periodic(2)    = .true.
  integer            :: nx(2)   = 512
  integer            :: bx(2)  = 16
  integer            :: n_gc           = 2
  logical            :: write_output   = .false.
  logical            :: write_all_output   = .false.
  real(dp)           :: grid_length(2) = [1.0_dp, 1.0_dp]
  real(dp)           :: dt = 1e-3_dp
  real(dp)           :: cfl_number     = 0.5_dp

  integer            :: n, n_steps_time_integrator
  integer            :: t_start, t_end, count_rate
  real(dp)           :: t_total, unknowns_per_ns

  call CFG_update_from_arguments(cfg)
  call CFG_add_get(cfg, 'n_iter', n_iter, 'Number of iterations')
  call CFG_add_get(cfg, 'periodic', periodic, 'Whether the domain is periodic')
  call CFG_add_get(cfg, 'nx', nx, 'Size of uniform grid')
  call CFG_add_get(cfg, 'bx', bx, 'Size of grid blocks')
  call CFG_add_get(cfg, 'n_gc', n_gc, 'Number of ghost cells')
  call CFG_add_get(cfg, 'write_output', write_output, 'Whether to write output')
  call CFG_add_get(cfg, 'write_all_output', write_all_output, &
       'Whether to write all time states to output')
  call CFG_add_get(cfg, 'cfl_number', cfl_number, 'CFL number to use')
  call CFG_add_get(cfg, 'dt', dt, 'Time step')

  if (n_gc < 2) error stop "n_gc < 2"

  call initialize_grid(nx, grid_length, bx, n_gc, n_variables, &
       periodic, bg)
  call set_initial_conditions()

  write(*, "(A,E12.4)") "Using dt = ", dt

  if (write_output .or. write_all_output) then
     !$acc update self(bg%uu)
     call write_brick_of_values(bg, "test_rho", "rho", 0, 0.0_dp, i_rho)
     call write_brick_of_values(bg, "test_ux", "ux", 0, 0.0_dp, i_momx)
     call write_brick_of_values(bg, "test_uy", "uy", 0, 0.0_dp, i_momy)
     call write_brick_of_values(bg, "test_e", "e", 0, 0.0_dp, i_e)
  end if

  call system_clock(t_start, count_rate)
  do n = 1, n_iter
     call advance_heuns_method(bg, dt)

     if (write_all_output) then
        !$acc update self(bg%uu)
        call write_brick_of_values(bg, "test_rho", "rho", n, n*dt, i_rho)
        call write_brick_of_values(bg, "test_ux", "ux", n, n*dt, i_momx)
        call write_brick_of_values(bg, "test_uy", "uy", n, n*dt, i_momy)
        call write_brick_of_values(bg, "test_e", "e", n, n*dt, i_e)
     end if
  end do
  call system_clock(t_end, count_rate)

  if (write_output .and. .not. write_all_output) then
     !$acc update self(bg%uu)
     call write_brick_of_values(bg, "test_rho", "rho", 1, n_iter*dt, i_rho)
     call write_brick_of_values(bg, "test_ux", "ux", 1, n_iter*dt, i_momx)
     call write_brick_of_values(bg, "test_uy", "uy", 1, n_iter*dt, i_momy)
     call write_brick_of_values(bg, "test_e", "e", 1, n_iter*dt, i_e)
  end if

  t_total = (t_end - t_start)/real(count_rate, dp)
  n_steps_time_integrator = 2
  unknowns_per_ns = 1e-9_dp * n_iter/t_total * product(bg%nx) * &
       n_steps_time_integrator

  write(*, "(5(A6,' '),2(A8,' '),2A12)") "nx", "ny", "box_nx", "box_ny", &
       "n_gc", "n_iter", "n_boxes", "t_total", "cells/ns"
  write(*, "(5(I6,' '),2(I8,' '),2F12.6)") bg%nx, bg%bx, bg%n_gc, &
       n_iter, bg%n_blocks, t_total, unknowns_per_ns

contains

  pure subroutine to_primitive(u)
    !$acc routine seq
    real(dp), intent(inout) :: u(n_vars_euler)
    integer :: idim

    u(ix_momx) = u(ix_momx)/u(ix_rho)
    u(ix_momy) = u(ix_momy)/u(ix_rho)

    u(ix_e) = (euler_gamma-1.0_dp) * (u(ix_e) - &
         0.5_dp * u(ix_rho)* (u(ix_momx)**2 + u(ix_momy)**2))

  end subroutine to_primitive

  pure subroutine to_conservative(u)
    !$acc routine seq
    real(dp), intent(inout) :: u(n_vars_euler)

    ! Compute energy from pressure and kinetic energy
    u(ix_e) = u(ix_e) * inv_gamma_m1 + &
         0.5_dp * u(ix_rho) * (u(ix_momx)**2 + u(ix_momy)**2)

    ! Compute momentum from density and velocity components
    u(ix_momx) = u(ix_rho) * u(ix_momx)
    u(ix_momy) = u(ix_rho) * u(ix_momy)
  end subroutine to_conservative

  subroutine muscl_flux_euler_prim(u, flux_dim, flux)
    !$acc routine seq
    real(dp), intent(in)  :: u(5, n_vars_euler)
    integer, intent(in)   :: flux_dim
    real(dp), intent(out) :: flux(n_vars_euler, 2)
    real(dp)              :: uL(n_vars_euler), uR(n_vars_euler), wL, wR, wmax
    real(dp)              :: flux_l(n_vars_euler), flux_r(n_vars_euler)

    ! Construct uL, uR for first cell face
    uL = u(2, :) + 0.5_dp * vanleer(u(2, :) - u(1, :), u(3, :) - u(2, :))
    uR = u(3, :) - 0.5_dp * vanleer(u(3, :) - u(2, :), u(4, :) - u(3, :))

    call euler_flux(uL, flux_dim, flux_l, wL)
    call euler_flux(uR, flux_dim, flux_r, wR)
    wmax = max(wL, wR)

    call to_conservative(uL)
    call to_conservative(uR)
    flux(:, 1) = 0.5_dp * (flux_l + flux_r - wmax * (uR - uL))

    ! Construct uL, uR for second cell face
    uL = u(3, :) + 0.5_dp * vanleer(u(3, :) - u(2, :), u(4, :) - u(3, :))
    uR = u(4, :) - 0.5_dp * vanleer(u(4, :) - u(3, :), u(5, :) - u(4, :))

    call euler_flux(uL, flux_dim, flux_l, wL)
    call euler_flux(uR, flux_dim, flux_r, wR)
    wmax = max(wL, wR)

    call to_conservative(uL)
    call to_conservative(uR)
    flux(:, 2) = 0.5_dp * (flux_l + flux_r - wmax * (uR - uL))

  end subroutine muscl_flux_euler_prim

  subroutine euler_flux(u, flux_dim, flux, w)
    !$acc routine seq
    real(dp), intent(in)  :: u(n_vars_euler)
    integer, intent(in)   :: flux_dim
    real(dp), intent(out) :: flux(n_vars_euler)
    real(dp), intent(out) :: w

    ! Density flux
    flux(ix_rho) = u(ix_rho) * u(ix_mom(flux_dim))

    ! Momentum flux with pressure term
    flux(ix_momx) = u(ix_rho) * u(ix_momx) * u(ix_mom(flux_dim))
    flux(ix_momy) = u(ix_rho) * u(ix_momy) * u(ix_mom(flux_dim))
    flux(ix_mom(flux_dim)) = flux(ix_mom(flux_dim)) + u(ix_e)

    ! Energy flux
    flux(ix_e) = u(ix_mom(flux_dim)) * (u(ix_e) * inv_gamma_m1 + &
         0.5_dp * u(ix_rho) * (u(ix_momx)**2 + u(ix_momy)**2) + u(ix_e))

    w = sqrt(euler_gamma * u(ix_e) / u(ix_rho)) + abs(u(ix_mom(flux_dim)))
  end subroutine euler_flux

  subroutine advance_heuns_method(bg, dt)
    type(block_grid_t), intent(inout) :: bg
    real(dp), intent(in)              :: dt

    call forward_euler(bg, bg%bx, bg%ilo, bg%ihi, bg%n_vars, bg%n_blocks, &
         dt, bg%uu, 0, 1, [0], [1.0_dp], 1)
    call forward_euler(bg, bg%bx, bg%ilo, bg%ihi, bg%n_vars, bg%n_blocks, &
         0.5_dp*dt, bg%uu, 1, 2, [0, 1], [0.5_dp, 0.5_dp], 0)
  end subroutine advance_heuns_method

  subroutine set_initial_conditions()
    integer  :: n, i, j, ix(2)
    real(dp) :: r(2)
    real(dp) :: u0(n_vars_euler, 4)

    select case ("sod")
    case ("first")
       u0(ix_e, :)      = [1.0_dp, 0.4_dp, 0.0439_dp, 0.15_dp]
       u0(ix_rho, :)    = [1.0_dp, 0.5197_dp, 0.1072_dp, 0.2579_dp]
       u0(ix_momx, :) = [0.0_dp, -0.7259_dp, -0.7259_dp, 0.0_dp]
       u0(ix_momy, :) = [0.0_dp, 0.0_dp, -1.4045_dp, -1.4045_dp]
    case ("sixth")
       u0(ix_e, :)      = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
       u0(ix_rho, :)    = [1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp]
       u0(ix_momx, :) = [0.75_dp, 0.75_dp, -0.75_dp, -0.75_dp]
       u0(ix_momy, :) = [-0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp]
    case ("sod")
       ! 1D Sod shock test case
       u0(ix_rho, :)    = [0.125_dp, 1.0_dp, 1.0_dp, 0.125_dp]
       u0(ix_e, :)      = [0.1_dp, 1.0_dp, 1.0_dp, 0.1_dp]
       u0(ix_momx, :) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
       u0(ix_momy, :) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    case default
       error stop "Unknown initial condition"
    end select

    do n = 1, 4
       call to_conservative(u0(:, n))
    end do

    !$acc parallel loop private(r, ix)
    do n = 1, bg%n_blocks
       !$acc loop collapse(2)
       do j = 1, bg%bx(2)
          do i = 1, bg%bx(1)
             ix = [i, j] + (bg%id_to_index(:, n) - 1) * bg%bx
             r = bg%dr * (ix - 0.5_dp)

             if (r(1) > 0.5_dp .and. r(2) > 0.5_dp) then
                bg%uu(IX(i, j, n, i_vars_grid)) = u0(:, 1)
             elseif (r(1) <= 0.5_dp .and. r(2) >= 0.5_dp) then
                bg%uu(IX(i, j, n, i_vars_grid)) = u0(:, 2)
             elseif (r(1) <= 0.5_dp .and. r(2) <= 0.5_dp) then
                bg%uu(IX(i, j, n, i_vars_grid)) = u0(:, 3)
             else
                bg%uu(IX(i, j, n, i_vars_grid)) = u0(:, 4)
             end if
          end do
       end do
    end do
  end subroutine set_initial_conditions

  subroutine forward_euler(bg, bx, lo, hi, n_vars, n_blocks, dt, uu, &
       s_deriv, n_prev, s_prev, w_prev, s_out)
    type(block_grid_t), intent(inout) :: bg
    integer, intent(in)               :: n_blocks, bx(2), lo(2), hi(2), n_vars
    real(dp), intent(in)              :: dt
    real(dp), intent(inout)           :: uu(IX(lo(1):hi(1), lo(2):hi(2), n_blocks, n_vars))
    integer, intent(in)               :: s_deriv        !< State to compute derivatives from
    integer, intent(in)               :: n_prev         !< Number of previous states
    integer, intent(in)               :: s_prev(n_prev) !< Previous states
    real(dp), intent(in)              :: w_prev(n_prev) !< Weights of previous states
    integer, intent(in)               :: s_out          !< Output state
    integer                           :: n, i, j, m, iv
    real(dp)                          :: inv_dr(2)
    real(dp)                          :: fx(n_vars_euler, 2), fy(n_vars_euler, 2)
    real(dp)                          :: tmp(5, n_vars_euler), u(n_vars_euler)
    real(dp)                          :: uprim(lo(1):hi(1), lo(2):hi(2), n_vars_euler)
    real(dp)                          :: dvar(bx(1), bx(2), n_vars_euler)

    do n = 1, n_vars_euler
       call update_ghostcells(bg, i_vars_grid(n)+s_deriv)
    end do

    inv_dr = 1/bg%dr

    !$acc parallel loop private(fx, fy, tmp, u, uprim, dvar, m, iv)
    do n = 1, n_blocks

       !$acc loop collapse(2)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             ! Convert to primitive
             u = uu(IX(i, j, n, i_vars_grid+s_deriv))
             call to_primitive(u)
             uprim(i, j, :) = u
          end do
       end do

       !$acc loop collapse(2)
       do j = 1, bx(2)
          do i = 1, bx(1)
             ! Compute x and y fluxes
             tmp = uprim(i-2:i+2, j, :)
             call muscl_flux_euler_prim(tmp, 1, fx)

             tmp = uprim(i, j-2:j+2, :)
             call muscl_flux_euler_prim(tmp, 2, fy)

             ! Keep track of changes in variables
             dvar(i, j, :) = dt * &
                  ((fx(:, 1) - fx(:, 2)) * inv_dr(1) + &
                  (fy(:, 1) - fy(:, 2)) * inv_dr(2))
          end do
       end do

       ! Set output state after computations are done, since s_out can be
       ! equal to s_deriv and s_prev
       !$acc loop collapse(3)
       do iv = 1, n_vars_euler
          do j = 1, bx(2)
             do i = 1, bx(1)
                do m = 1, n_prev
                   ! Add weighted previous states
                   dvar(i, j, iv) = dvar(i, j, iv) + &
                        uu(IX(i, j, n, i_vars_grid(iv)+s_prev(m))) * w_prev(m)
                end do
                uu(IX(i, j, n, i_vars_grid(iv)+s_out)) = dvar(i, j, iv)
             end do
          end do
       end do
    end do
  end subroutine forward_euler

  elemental pure real(dp) function minmod(a, b)
    real(dp), intent(in) :: a, b

    if (a * b <= 0) then
       minmod = 0.0_dp
    else if (abs(a) < abs(b)) then
       minmod = a
    else
       minmod = b
    end if
  end function minmod

  elemental pure real(dp) function vanleer(a, b) result(phi)
    real(dp), intent(in) :: a, b
    real(dp)             :: ab

    ab = a * b
    if (ab > 0) then
       phi = 2 * ab / (a + b)
    else
       phi = 0
    end if
  end function vanleer

end program
