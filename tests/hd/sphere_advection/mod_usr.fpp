module mod_usr
  use mod_amrvac
  use mod_physics

  implicit none

  double precision :: rhodens  = 10.0d0
  double precision :: rholight = 1.0d0
  !$acc declare copyin(rhodens, rholight)
  double precision :: p        = 10.0d0
  double precision :: rsphere  = 0.25d0
  double precision :: v1 = 1.0d0, v2 = 1.0d0, v3 = 1.0d0

contains

  subroutine usr_init()

    call set_coordinate_system("Cartesian_3D")

    usr_init_one_grid => initonegrid_usr

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
    ! .. local ..
    double precision                :: r2
    integer                         :: ix1, ix2, ix3

    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1

             r2 = x(ix1,ix2,ix3,1)**2 + x(ix1,ix2,ix3,2)**2 + x(ix1,ix2,ix3,3)**2

             if ( r2 < rsphere**2 ) then
                w(ix1,ix2,ix3,rho_) = rhodens
             else
                w(ix1,ix2,ix3,rho_) = rholight

             end if
             
             w(ix1,ix2,ix3,p_)     = p
             w(ix1,ix2,ix3,mom(1)) = v1
             w(ix1,ix2,ix3,mom(2)) = v2
             w(ix1,ix2,ix3,mom(3)) = v3
                
          end do
       end do
    end do
    
    call phys_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

  end subroutine initonegrid_usr

  subroutine usr_refine_grid(igrid,level,ixGmin1,ixGmin2,ixGmin3,&
    ixGmax1,ixGmax2,ixGmax3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,&
    qt,w,x,refine,coarsen)

    use mod_global_parameters
    !$acc routine vector
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

    ! you must set consistent values for integers refine/coarsen:

    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement

    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen

    integer, intent(in)             :: igrid, level, ixGmin1,ixGmin2,&
        ixGmin3,ixGmax1,ixGmax2,ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,&
        ixmax2,ixmax3
    double precision, intent(in)    :: qt, x(ixGmin1:ixGmax1,&
         ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim)
    double precision, intent(in)    :: w(ixGmin1:ixGmax1,&
         ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:nw)
    integer, intent(inout) :: refine, coarsen
    ! .. local ..
    integer                         :: ix1, ix2, ix3
    logical                         :: has_sphere

    has_sphere = .false.
    !$acc loop collapse(3) reduction(.or.:has_sphere)
    do ix3 = ixmin3, ixmax3
       do ix2 = ixmin2, ixmax2
          do ix1 = ixmin1, ixmax1
             if ( w(ix1, ix2, ix3, rho_) > (rhodens-rholight)/2.0d0 ) then
                has_sphere = .true.
             end if
          end do
       end do
    end do

    if (has_sphere) then
       coarsen = -1
       refine  = 1
    else 
       coarsen = 0
       refine  = 0
    end if

  end subroutine usr_refine_grid

end module mod_usr
