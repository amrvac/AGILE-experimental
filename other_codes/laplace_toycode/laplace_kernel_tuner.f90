#include "macros.h"
! Toy code for solving a 2D Laplace equation on a uniform grid with OpenACC. The
! grid is divided into blocks with a layer of ghost cells around them.
!
! Author: Jannis Teunissen (jannis.teunissen@cwi.nl)
program ghostcell
  use m_blockgrid

  implicit none
  !$tuner initialize
  integer, parameter :: dp = kind(0.0d0)
  !$tuner stop
  integer, parameter :: NDIM = 2

  integer :: grid_size(NDIM)
  integer :: block_size(NDIM)
  integer :: n_variables
  integer :: n_iterations
  integer :: n_gc
  logical :: periodic(2) = .false.

  integer :: n_args
  character(len=100) :: tmp

  type(block_grid_t) :: bg

  integer :: n
  integer :: i_new, i_old
  integer :: t_start, t_end, count_rate
  real(dp) :: t_total, unknowns_per_ns

  ! Set default values
  grid_size = 512
  block_size = 16
  n_gc = 1                      ! number of ghost cells
  n_iterations = 100
  n_variables = 2

  call initialize_grid(grid_size, block_size, n_gc, n_variables, periodic, bg)
  call set_initial_conditions()

  ! Copy data structure and allocatable components to device
  !$acc enter data copyin(bg, bg%uu, bg%bnd_x, bg%bnd_y, bg%bc_xlo, &
  !$acc bg%bc_xhi, bg%bc_ylo, bg%bc_yhi)

  call system_clock(t_start, count_rate)
  do n = 1, n_iterations
     i_new = mod(n, 2) + 1
     i_old = mod(n+1, 2) + 1
     call update_ghostcells(bg, i_old)
     call update_solution(bg%bx, bg%ilo, bg%ihi, bg%n_vars, bg%n_blocks, &
          i_new, i_old, bg%uu)
  end do
  call system_clock(t_end, count_rate)

  !$acc exit data copyout(bg%uu)

  t_total = (t_end - t_start)/real(count_rate, dp)
  unknowns_per_ns = 1e-9_dp * n_iterations/t_total * product(grid_size)

  write(*, "(5(A6,' '),2(A8,' '),2A12)") "nx", "ny", "box_nx", "box_ny", &
       "n_gc", "n_iter", "n_boxes", "t_total", "unknowns/ns"
  write(*, "(5(I6,' '),2(I8,' '),2F12.6)") bg%nx, bg%bx, bg%n_gc, &
       n_iterations, bg%n_blocks, t_total, unknowns_per_ns

contains

  subroutine set_initial_conditions()
    ! Zero everywhere
    bg%uu(IX(:, :, :, :)) = 0.0_dp

    ! Fix boundary value
    bg%bc(:) = 1.0_dp
  end subroutine set_initial_conditions

  subroutine update_solution(bx, lo, hi, n_vars, n_blocks, i_new, i_old, uu)
    integer, intent(in)     :: n_blocks, bx(2), lo(2), hi(2), i_new, i_old, n_vars
    real(dp), intent(inout) :: uu(IX(lo(1):hi(1), lo(2):hi(2), n_blocks, n_vars))
    !$tuner initialize
    integer                 :: n, i, j
    !$tuner stop

    !$tuner start update_solution uu(double*:block_size,block_size,n_blocks,n_vars) i_old(int:1) i_new(int:2)
    !$acc parallel loop vector_length(vlength) collapse(collapse_factor)
    do n = 1, n_blocks
       do j = 1, block_size
          do i = 1, block_size
             uu(IX(i, j, n, i_new)) = 0.25_dp * ( &
                  uu(IX(i+1, j, n, i_old)) + uu(IX(i-1, j, n, i_old)) + &
                  uu(IX(i, j-1, n, i_old)) + uu(IX(i, j+1, n, i_old)))
          end do
       end do
    end do
    !$tuner stop
  end subroutine update_solution

end program ghostcell
