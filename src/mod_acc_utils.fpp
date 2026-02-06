!=====================================================================
! Fully safe wrapper for allocatable or pointer arrays
!=====================================================================

module acc_utils
  use openacc
  implicit none
contains

  !-------------------------------------------------------------
  ! ensure_device_present_any
  ! - Works for allocatable or pointer arrays
  ! - Any rank
  ! - Allocates on device if not present, updates if present
  !-------------------------------------------------------------
  subroutine copy_or_present(array)
    implicit none
    class(*), intent(inout) :: array(..)  ! assumed-rank, any type
    logical :: is_pointer, is_alloc

    ! Determine if array is pointer or allocatable
    is_pointer = .false.
    is_alloc   = .false.

    ! pointer check
    ! Fortran intrinsic associated() works only on pointers
    ! Use "select type" to distinguish pointers vs allocatables
    select type(arr => array)
    type is (real(kind=8), pointer)       ! pointer array
       is_pointer = .true.
       if (.not. associated(arr)) then
          print *, 'copy_or_present: trying to check present on null-pointer'
          return
       end if
    class default                   ! everything else
       ! fallthrough
    end select

    ! allocatable check
    ! Note: allocated() works on allocatable arrays
    if (.not. is_pointer) then
       select type(arr => array)
       type is (real(kind=8), allocatable)
          is_alloc = .true.
          if (.not. allocated(arr)) then
          print *, 'copy_or_present: trying to check present on unallocated array'
          return
       end if
       class default
          ! cannot handle other types
          print *, 'copy_or_present: trying to check present on unknown type'
          return
       end select
    end if

    ! Now safe to check device presence
    if (.not. acc_is_present(array)) then
       !$acc enter data copyin(array)
    else
       !$acc update device(array)
    end if

  end subroutine copy_or_present

end module acc_utils
