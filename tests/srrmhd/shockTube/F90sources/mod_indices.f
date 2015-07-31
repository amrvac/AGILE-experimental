module mod_indices
   implicit none

   !AMR specific parameters
   !
   ! ngridshi:  maximum number of grids
   ! nlevelshi: maximum number of levels in grid refinement

   integer, parameter :: ngridshi  = 4000
   integer, parameter :: nlevelshi = 13

   ! number of interleaving sending buffers for ghostcells
   integer, parameter :: npwbuf=2

   integer, save :: ixMlo1,ixMhi1

   integer, dimension(nlevelshi), save :: ng1
   double precision, dimension(nlevelshi), save :: dg1

   logical, save :: slab

   integer, save :: npe, mype, icomm
   integer, save :: log_fh
   integer, save :: type_block, type_coarse_block, type_sub_block(2)
   integer, save :: type_block_io, size_block_io

   integer, save :: type_subblock_io, type_subblock_x_io
   integer, save :: type_block_xc_io,type_block_xcc_io
   integer, save :: type_block_wc_io,type_block_wcc_io
   integer, save :: itag, ierrmpi
   integer, save :: irecv, isend
   integer, dimension(:), allocatable, save :: recvrequest, sendrequest
   integer, dimension(:,:), allocatable, save :: recvstatus, sendstatus

   integer, save :: snapshot, snapshotnext, slice, slicenext, collapseNext,&
       icollapse

   logical, allocatable, dimension(:), save :: patchfalse

   logical, save :: B0field
   double precision, save :: Bdip, Bquad, Boct, Busr

end module mod_indices