program openacc_test
  implicit none
  integer, parameter :: n = 1000000
  double precision, dimension(n) :: a, b, c
  integer :: i
  double precision :: sum_cpu, sum_gpu

  ! Initialize arrays on CPU
  do i = 1, n
     a(i) = real(i)
     b(i) = 2.0 * real(i)
  end do

  !---------------------------------------------------------
  ! Part 1: Run on CPU only
  !---------------------------------------------------------
  sum_cpu = 0.0
  do i = 1, n
     sum_cpu = sum_cpu + a(i) * b(i)
  end do
  print *, "Dot product on CPU: ", sum_cpu

  !---------------------------------------------------------
  ! Part 2: Run on GPU using OpenACC
  !---------------------------------------------------------
  sum_gpu = 0.0
  !$acc data copyin(a,b) copy(sum_gpu)
  !$acc parallel loop reduction(+:sum_gpu)
  do i = 1, n
     sum_gpu = sum_gpu + a(i) * b(i)
  end do
  !$acc end parallel loop
  !$acc end data

  print *, "Dot product on GPU: ", sum_gpu
  print *, "Relative difference: ", abs(sum_cpu - sum_gpu) / abs(sum_cpu)

end program openacc_test

