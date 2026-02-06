!=====================================================================
! Fully safe wrapper for allocatable or pointer arrays
!=====================================================================

module acc_utils
  use openacc
  implicit none


  ! Adjust this as needed
  #:set MAXRANK = 6
  #:set TYPES = [("double", "double precision"), ("logical", "logical"), ("integer", "integer")]
                 
  ! Generic interface for arrays
  interface copy_or_update
  #:for tname, tdecl in TYPES
  #:for rank in range(1, MAXRANK+1)
     module procedure copy_or_update_${tname}$_nonpointer_r${rank}$
  #:endfor
  #:endfor
     
 ! Generic interface for various scalar types
  #:for tname, tdecl in TYPES
    module procedure copy_or_update_${tname}$_nonpointer
  #:endfor
  end interface copy_or_update


  
  ! Generic interface for arrays
  interface copy_or_update_pointer
  #:for tname, tdecl in TYPES
  #:for rank in range(1, MAXRANK+1)
     module procedure copy_or_update_${tname}$_pointer_r${rank}$
  #:endfor
  #:endfor
     
 ! Generic interface for various scalar types
  #:for tname, tdecl in TYPES
    module procedure copy_or_update_${tname}$_pointer
  #:endfor
 end interface copy_or_update_pointer


 
  interface copy_or_update_alloc
  #:for tname, tdecl in TYPES
  #:for rank in range(1, MAXRANK+1)
     module procedure copy_or_update_${tname}$_alloc_r${rank}$
  #:endfor
  #:endfor
  end interface copy_or_update_alloc

     
contains

  ! Versions for scalars
  #:for tname, tdecl in TYPES

  subroutine copy_or_update_${tname}$_pointer(scalar)
    implicit none
    ${tdecl}$, pointer, intent(inout) :: scalar
    
    if (.not. associated(scalar)) then
       print *, 'copy_or_update: null-pointer encountered'
       return
    end if

    if (.not. acc_is_present(scalar, sizeof(scalar))) then
       !$acc enter data copyin(scalar)
    else
       !$acc update device(scalar)
    end if
    
  end subroutine copy_or_update_${tname}$_pointer

  subroutine copy_or_update_${tname}$_nonpointer(scalar)
    implicit none
    ${tdecl}$, intent(inout) :: scalar

    if (.not. acc_is_present(scalar, sizeof(scalar))) then
       !$acc enter data copyin(scalar)
    else
       !$acc update device(scalar)
    end if
    
  end subroutine copy_or_update_${tname}$_nonpointer

  #:endfor

  ! Versions for various array ranks
  
  #:for tname, tdecl in TYPES
  ! Generate non-pointer and pointer versions for ranks 1..MAXRANK
  #:for rank in range(1, MAXRANK+1)
  #:set shp = '(' + ','.join([':']*rank) + ')'
  
  !------------------------------
  ! Non-pointer, non-allocatable
  !------------------------------
  subroutine copy_or_update_${tname}$_nonpointer_r${rank}$(array)
    implicit none
    double precision, intent(inout) :: array${shp}$

    if (.not. acc_is_present(array)) then
       !$acc enter data copyin(array)
    else
       !$acc update device(array)
    end if
    
  end subroutine copy_or_update_${tname}$_nonpointer_r${rank}$

  !-------------
  ! Pointer case
  !-------------
  subroutine copy_or_update_${tname}$_pointer_r${rank}$(array)
    implicit none
    double precision, pointer, intent(inout) :: array${shp}$

    if (.not. associated(array)) then
       print *, 'copy_or_update: null-pointer encountered (rank ${rank}$)'
       return
    end if

    if (.not. acc_is_present(array)) then
       !$acc enter data copyin(array)
       !$acc enter data attach(array)
    else
       ! Assume the TARGET has been updated, so skip update:
!       !$acc update device(array)
    end if
    
  end subroutine copy_or_update_${tname}$_pointer_r${rank}$

  !----------------
  ! Allocatable case
  !----------------
  subroutine copy_or_update_${tname}$_alloc_r${rank}$(array)
    implicit none
    double precision, allocatable, intent(inout) :: array${shp}$

    if (.not. allocated(array)) then
       print *, 'copy_or_update: unallocated allocatable (rank ${rank}$)'
       return
    end if
    
    if (.not. acc_is_present(array)) then
       !$acc enter data copyin(array)
    else
       !$acc update device(array)
    end if
    
  end subroutine copy_or_update_${tname}$_alloc_r${rank}$

  #:endfor
  #:endfor


end module acc_utils
