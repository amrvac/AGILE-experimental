! TODO: make sure methods do nothing by default
! TODO: document methods
module mod_usr_methods

  implicit none
  public

  ! Initialize grid data
  procedure(set_parameters), pointer  :: usr_set_parameters   => null()
  procedure(init_one_grid), pointer   :: usr_init_one_grid    => null()

  ! Boundary condition related
  procedure(special_bc), pointer      :: usr_special_bc       => null()
  procedure(internal_bc), pointer     :: usr_internal_bc      => null()

  ! Output related
  procedure(p_no_args), pointer       :: usr_print_log        => null()
  procedure(transform_w), pointer     :: usr_transform_w      => null()
  procedure(aux_output), pointer      :: usr_aux_output       => null()
  procedure(add_aux_names), pointer   :: usr_add_aux_names    => null()
  procedure(special_convert), pointer :: usr_special_convert  => null()

  ! Called every time step
  procedure(process_grid), pointer    :: usr_process_grid     => null()
  procedure(process_global), pointer  :: usr_process_global   => null()

  ! Called before the start of the simulation
  procedure(p_no_args), pointer       :: usr_before_main_loop => null()

  ! Source terms
  procedure(source), pointer          :: usr_source           => null()
  procedure(get_dt), pointer          :: usr_get_dt           => null()
  procedure(phys_gravity), pointer    :: usr_gravity          => null()

  ! Refinement related procedures
  procedure(refine_grid), pointer     :: usr_refine_grid      => null()
  procedure(var_for_errest), pointer  :: usr_var_for_errest   => null()
  procedure(a_refine_threshold), pointer :: usr_refine_threshold => null()
  procedure(flag_grid), pointer       :: usr_flag_grid        => null()

  ! Set time-independent magnetic field for B0 splitting
  procedure(set_B0), pointer          :: usr_set_B0           => null()
  ! Set time-independent current density for B0 splitting
  procedure(set_J0), pointer          :: usr_set_J0           => null()
  procedure(special_resistivity), pointer :: usr_special_resistivity => null()

  ! Particles module related
  procedure(update_payload), pointer    :: usr_update_payload    => null()
  procedure(create_particles), pointer  :: usr_create_particles  => null()
  procedure(particle_fields), pointer   :: usr_particle_fields   => null()
  procedure(particle_analytic), pointer :: usr_particle_analytic => null()

  abstract interface

     subroutine p_no_args()
     end subroutine p_no_args

     !> Initialize the user's settings (after initializing amrvac)
     subroutine set_parameters()
       use mod_global_parameters
     end subroutine set_parameters

     !> Initialize one grid
     subroutine init_one_grid(ixI^L,ixO^L,w,x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: x(ixI^S,1:ndim)
       double precision, intent(inout) :: w(ixI^S,1:nw)
     end subroutine init_one_grid

     !> special boundary types, user defined user must assign conservative
     !> variables in bounderies
     subroutine special_bc(qt,ixI^L,ixO^L,iB,w,x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, iB
       double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
       double precision, intent(inout) :: w(ixI^S,1:nw)
     end subroutine special_bc

     !> internal boundary, user defined
     !
     !> This subroutine can be used to artificially overwrite ALL conservative 
     !> variables in a user-selected region of the mesh, and thereby act as
     !> an internal boundary region. It is called just before external (ghost cell)
     !> boundary regions will be set by the BC selection. Here, you could e.g. 
     !> want to introduce an extra variable (nwextra, to be distinguished from nwaux)
     !> which can be used to identify the internal boundary region location.
     !> Its effect should always be local as it acts on the mesh.
     subroutine internal_bc(level,qt,ixI^L,ixO^L,w,x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L,ixO^L,level
       double precision, intent(in)    :: qt
       double precision, intent(inout) :: w(ixI^S,1:nw)
       double precision, intent(in)    :: x(ixI^S,1:ndim)
     end subroutine internal_bc

     !> this subroutine is ONLY to be used for computing auxiliary variables
     !> which happen to be non-local (like div v), and are in no way used for
     !> flux computations. As auxiliaries, they are also not advanced
     subroutine process_grid(igrid,level,ixI^L,ixO^L,qt,w,x)
       use mod_global_parameters
       integer, intent(in)             :: igrid,level,ixI^L,ixO^L
       double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
       double precision, intent(inout) :: w(ixI^S,1:nw)
     end subroutine process_grid

     !> This subroutine is called at the beginning of each time step 
     !> by each processor. No communication is specified, so the user
     !> has to implement MPI routines if information has to be shared
     subroutine process_global(iit,qt)
       use mod_global_parameters
       integer, intent(in)          :: iit
       double precision, intent(in) :: qt
     end subroutine process_global

     !> this subroutine can be used in convert, to add auxiliary variables to the
     !> converted output file, for further analysis using tecplot, paraview, ....
     !> these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
     !
     !> the array normconv can be filled in the (nw+1:nw+nwauxio) range with 
     !> corresponding normalization values (default value 1)
     subroutine aux_output(ixI^L,ixO^L,w,x,normconv)
       use mod_global_parameters
       integer, intent(in)          :: ixI^L,ixO^L
       double precision, intent(in) :: x(ixI^S,1:ndim)
       double precision             :: w(ixI^S,nw+nwauxio)
       double precision             :: normconv(0:nw+nwauxio)
     end subroutine aux_output

     !> Add names for the auxiliary variables
     subroutine add_aux_names(varnames)
       use mod_global_parameters
       character(len=*) :: varnames
     end subroutine add_aux_names

     !> Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
     !> iw=iwmin...iwmax.  wCT is at time qCT
     subroutine source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
       double precision, intent(in)    :: qdt, qtC, qt
       double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
       double precision, intent(inout) :: w(ixI^S,1:nw)
     end subroutine source

     !> Limit "dt" further if necessary, e.g. due to the special source terms.
     !> The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
     !> module have already been called.
     subroutine get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim)
       double precision, intent(inout) :: w(ixI^S,1:nw), dtnew
     end subroutine get_dt

     !> Calculate gravitational acceleration in each dimension
     subroutine phys_gravity(ixI^L,ixO^L,wCT,x,gravity_field)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       double precision, intent(in)    :: x(ixI^S,1:ndim)
       double precision, intent(in)    :: wCT(ixI^S,1:nw)
       double precision, intent(out)   :: gravity_field(ixI^S,ndim)
     end subroutine phys_gravity

     !> Set the "eta" array for resistive MHD based on w or the
     !> "current" variable which has components between idirmin and 3.
     subroutine special_resistivity(w,ixI^L,ixO^L,idirmin,x,current,eta)
       use mod_global_parameters
       integer, intent(in)          :: ixI^L, ixO^L, idirmin
       double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)
       double precision             :: current(ixI^S,7-2*ndir:3), eta(ixI^S)
     end subroutine special_resistivity

     !> Enforce additional refinement or coarsening
     !> One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
     !> you must set consistent values for integers refine/coarsen:
     !> refine = -1 enforce to not refine
     !> refine =  0 doesn't enforce anything
     !> refine =  1 enforce refinement
     !> coarsen = -1 enforce to not coarsen
     !> coarsen =  0 doesn't enforce anything
     !> coarsen =  1 enforce coarsen
     !> e.g. refine for negative first coordinate x < 0 as
     !> if (any(x(ix^S,1) < zero)) refine=1
     subroutine refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
       use mod_global_parameters
       integer, intent(in)          :: igrid, level, ixI^L, ixO^L
       double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
       integer, intent(inout)       :: refine, coarsen
     end subroutine refine_grid

     !> this is the place to compute a local auxiliary variable to be used
     !> as refinement criterion for the Lohner error estimator only
     !>  -->it is then requiring and iflag>nw
     !> note that ixO=ixI=ixG, hence the term local (gradients need special attention!)
     subroutine var_for_errest(ixI^L,ixO^L,iflag,w,var)
       use mod_global_parameters
       integer, intent(in)           :: ixI^L,ixO^L,iflag
       double precision, intent(in)  :: w(ixI^S,1:nw)
       double precision, intent(out) :: var(ixI^S)
     end subroutine var_for_errest

     !> Here one can add a steady (time-independent) potential background field
     subroutine set_B0(ixI^L,ixO^L,x,wB0)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L,ixO^L
       double precision, intent(in)    :: x(ixI^S,1:ndim)
       double precision, intent(inout) :: wB0(ixI^S,1:ndir)
     end subroutine set_B0

     !> Here one can add a time-independent background current density
     subroutine set_J0(ixI^L,ixO^L,x,wJ0)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L,ixO^L
       double precision, intent(in)    :: x(ixI^S,1:ndim)
       double precision, intent(inout) :: wJ0(ixI^S,7-2*ndir:ndir)
     end subroutine set_J0

     !> regenerate w and eqpartf arrays to output into *tf.dat
     subroutine transform_w(ixI^L,ixO^L,nw_in,w_in,x,w_out)
       use mod_global_parameters
       integer, intent(in)           :: ixI^L, ixO^L, nw_in
       double precision, intent(in)  :: w_in(ixI^S,1:nw_in)
       double precision, intent(in)  :: x(ixI^S, 1:ndim)
       double precision, intent(out) :: w_out(ixI^S,1:nw)
     end subroutine transform_w

     !> use different threshold in special regions for AMR to
     !> reduce/increase resolution there where nothing/something interesting happens.
     subroutine a_refine_threshold(wlocal,xlocal,threshold,qt)
       use mod_global_parameters
       double precision, intent(in)    :: wlocal(1:nw),xlocal(1:ndim),qt
       double precision, intent(inout) :: threshold
     end subroutine a_refine_threshold

     !> Allow user to use their own data-postprocessing procedures
     subroutine special_convert(qunitconvert)
       use mod_global_parameters
       integer, intent(in) :: qunitconvert
       character(len=20)   :: userconvert_type
     end subroutine special_convert

     subroutine fix_p(ixI^L,ixO^L,w,x)
       use mod_global_parameters
       integer, intent(in)                :: ixI^L,ixO^L
       double precision, intent(inout)    :: w(ixI^S,1:nw)
       double precision, intent(in)       :: x(ixI^S,1:ndim)
     end subroutine fix_p

     !> flag=-1 : Treat all cells active, omit deactivation (onentry, default)
     !> flag=0  : Treat as normal domain
     !> flag=1  : Treat as passive, but reduce by safety belt
     !> flag=2  : Always treat as passive
     subroutine flag_grid(qt,ixI^L,ixO^L,w,x,flag)
       use mod_global_parameters
       integer, intent(in)             :: ixI^L, ixO^L
       integer, intent(inout)          :: flag
       double precision, intent(in)    :: qt
       double precision, intent(inout) :: w(ixI^S,1:nw)
       double precision, intent(in)    :: x(ixI^S,1:ndim)
     end subroutine flag_grid

     !> Update payload of particles
     subroutine update_payload(igrid,w,wold,xgrid,xpart,payload,npayload,particle_time)
       use mod_global_parameters
       integer, intent(in)           :: igrid,npayload
       double precision, intent(in)  :: w(ixG^T,1:nw),wold(ixG^T,1:nw)
       double precision, intent(in)  :: xgrid(ixG^T,1:ndim),xpart(1:ndir),particle_time
       double precision, intent(out) :: payload(npayload)
     end subroutine update_payload

     subroutine create_particles(n_particles, x, v, q, m, follow)
       integer, intent(in)           :: n_particles
       double precision, intent(out) :: x(3, n_particles)
       double precision, intent(out) :: v(3, n_particles)
       double precision, intent(out) :: q(n_particles)
       double precision, intent(out) :: m(n_particles)
       logical, intent(out)          :: follow(n_particles)
     end subroutine create_particles

     subroutine particle_fields(w, x, E, B)
       use mod_global_parameters
       double precision, intent(in)  :: w(ixG^T,1:nw)
       double precision, intent(in)  :: x(ixG^T,1:ndim)
       double precision, intent(out) :: E(ixG^T, ndir)
       double precision, intent(out) :: B(ixG^T, ndir)
     end subroutine particle_fields

     subroutine particle_analytic(ix, x, tloc, vec)
       use mod_global_parameters
       integer, intent(in)           :: ix(ndir) !< Indices in gridvars
       double precision, intent(in)  :: x(ndir)
       double precision, intent(in)  :: tloc
       double precision, intent(out) :: vec(ndir)
     end subroutine particle_analytic

  end interface

end module mod_usr_methods
