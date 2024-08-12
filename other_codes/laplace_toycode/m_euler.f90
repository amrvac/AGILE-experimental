#include "macros.h"
! Module for Euler gas dynamics functionality
module m_euler

  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

  real(dp), parameter :: euler_gamma  = 1.4_dp
  real(dp), parameter :: inv_gamma_m1 = 1/(euler_gamma-1.0_dp)
  integer, parameter  :: n_vars_euler = 4
  integer, parameter  :: ix_rho       = 1
  integer, parameter  :: ix_momx      = 2
  integer, parameter  :: ix_momy      = 3
  integer, parameter  :: ix_mom(2)    = [2, 3]
  integer, parameter  :: ix_e         = 4

  integer, parameter :: n_variables               = 2 * n_vars_euler
  integer, parameter :: i_rho                     = 1
  integer, parameter :: i_momx                    = 3
  integer, parameter :: i_momy                    = 5
  integer, parameter :: i_e                       = 7
  integer, parameter :: i_vars_grid(n_vars_euler) = [i_rho, i_momx, i_momy, i_e]

end module m_euler
