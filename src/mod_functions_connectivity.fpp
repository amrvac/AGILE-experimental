module mod_functions_connectivity


  implicit none
  private

  public :: get_level_range
  public :: getigrids
  public :: build_connectivity
  contains

  subroutine get_level_range
    use mod_forest
    use mod_global_parameters

    integer :: level

    ! determine new finest level
    do level=refine_max_level,1,-1
       if (associated(level_tail(level)%node)) then
          levmax=level
          exit
       end if
    end do

    ! determine coarsest level
    do level=1,levmax
       if (associated(level_tail(level)%node)) then
          levmin=level
          exit
       end if
    end do

  end subroutine get_level_range

  subroutine getigrids
    use mod_forest
    use mod_global_parameters
    implicit none

    integer :: iigrid, igrid

    iigrid=0
    do igrid=1,max_blocks
       if (igrid_inuse(igrid,mype)) then
          iigrid=iigrid+1
          igrids(iigrid)=igrid
       end if
    end do

    igridstail=iigrid
    !$acc update device(igrids, igridstail)

  end subroutine getigrids

  subroutine build_connectivity
    use mod_forest
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_amr_neighbors, only: find_neighbor

    integer :: iigrid, igrid, i1,i2,i3, my_neighbor_type
    integer :: iside, idim, ic1,ic2,ic3, inc1,inc2,inc3, ih1,ih2,ih3, icdim
    type(tree_node_ptr) :: tree, my_neighbor, child
    logical, dimension(3) :: pole
    logical :: nopole
    ! Variables to detect special corners for stagger grid
    integer :: idir,pi1,pi2,pi3, mi1,mi2,mi3, ph1,ph2,ph3, mh1,mh2,mh3,&
        ipe_neighbor
    integer :: nrecvs,nsends

    ! total size of buffer arrays
    integer :: nbuff_bc_recv_srl, nbuff_bc_send_srl, nbuff_bc_recv_r,&
         nbuff_bc_send_r, nbuff_bc_recv_p, nbuff_bc_send_p

    nrecv_bc_srl=0; nsend_bc_srl=0
    nrecv_bc_r=0; nsend_bc_r=0
    nrecv_bc_p=0; nsend_bc_p=0
    nrecv_fc=0; nsend_fc=0
    nbuff_bc_recv_srl=0; nbuff_bc_send_srl=0
    nbuff_bc_recv_r=0; nbuff_bc_send_r=0
    nbuff_bc_recv_p=0; nbuff_bc_send_p=0
    if(stagger_grid) nrecv_cc=0; nsend_cc=0

    call nbprocs_info%reset

    do iigrid=1,igridstail; igrid=igrids(iigrid);
       tree%node => igrid_to_node(igrid,mype)%node

       do i3=-1,1
       do i2=-1,1
       do i1=-1,1
          ! skip the grid itself
          if (i1==0.and.i2==0.and.i3==0) then
             neighbor_type(0,0,0,igrid)=0
             neighbor(1,0,0,0,igrid)=igrid
             neighbor(2,0,0,0,igrid)=mype
          else
             call find_neighbor(my_neighbor,my_neighbor_type,tree,i1,i2,i3,&
                pole)
             nopole=.not.any(pole)

             select case (my_neighbor_type)
             ! adjacent to physical boundary
             case (neighbor_boundary)
                neighbor(1,i1,i2,i3,igrid)=0
                neighbor(2,i1,i2,i3,igrid)=-1
             ! fine-coarse transition
             case (neighbor_coarse)
                neighbor(1,i1,i2,i3,igrid)=my_neighbor%node%igrid
                neighbor(2,i1,i2,i3,igrid)=my_neighbor%node%ipe
                if (my_neighbor%node%ipe/=mype) then
                   ic1=1+modulo(tree%node%ig1-1,2)
                   ic2=1+modulo(tree%node%ig2-1,2)
                   ic3=1+modulo(tree%node%ig3-1,2);
                   ! avoid double counting of coarse neighbors
                   if ((i1==0.or.i1==2*ic1-3).and.(i2==0.or.i2==&
                        2*ic2-3).and.(i3==0.or.i3==2*ic3-3)) then
                     ! This is the local index of the prolonged ghost zone
                     inc1=ic1+i1;inc2=ic2+i2;inc3=ic3+i3;
                      call nbprocs_info%add_to_c(my_neighbor%node%ipe, igrid, i1, i2, i3, inc1, inc2, inc3)
                     nrecv_bc_p=nrecv_bc_p+1
                     nsend_bc_r=nsend_bc_r+1
                     nbuff_bc_send_r=nbuff_bc_send_r+sizes_r_send_total(i1,i2,&
                        i3)
                     nbuff_bc_recv_p=nbuff_bc_recv_p+sizes_p_recv_total(inc1,&
                        inc2,inc3)
                   end if
                end if
             ! same refinement level
             case (neighbor_sibling)
                neighbor(1,i1,i2,i3,igrid)=my_neighbor%node%igrid
                neighbor(2,i1,i2,i3,igrid)=my_neighbor%node%ipe
                if (my_neighbor%node%ipe/=mype) then
                   call nbprocs_info%add_to_srl(my_neighbor%node%ipe, igrid, i1, i2, i3)
                   nrecv_bc_srl=nrecv_bc_srl+1
                   nsend_bc_srl=nsend_bc_srl+1
                   nbuff_bc_send_srl=nbuff_bc_send_srl+sizes_srl_send_total(i1,&
                        i2,i3)
                   nbuff_bc_recv_srl=nbuff_bc_recv_srl+sizes_srl_recv_total(i1,&
                        i2,i3)
                end if
                ! coarse-fine transition
             case (neighbor_fine)
                neighbor(1,i1,i2,i3,igrid)=0
                neighbor(2,i1,i2,i3,igrid)=-1
                ! Loop over the local indices of children ic^D
                ! and calculate local indices of ghost zone inc^D.
                do ic3=1+int((1-i3)/2),2-int((1+i3)/2)
                   inc3=2*i3+ic3
                   if (pole(3)) then
                      ih3=3-ic3
                   else
                      ih3=ic3
                   end if
                do ic2=1+int((1-i2)/2),2-int((1+i2)/2)
                   inc2=2*i2+ic2
                   if (pole(2)) then
                      ih2=3-ic2
                   else
                      ih2=ic2
                   end if
                do ic1=1+int((1-i1)/2),2-int((1+i1)/2)
                   inc1=2*i1+ic1
                   if (pole(1)) then
                      ih1=3-ic1
                   else
                      ih1=ic1
                   end if
                   child%node => my_neighbor%node%child(ih1,ih2,ih3)%node
                   neighbor_child(1,inc1,inc2,inc3,igrid)=child%node%igrid
                   neighbor_child(2,inc1,inc2,inc3,igrid)=child%node%ipe
                   if (child%node%ipe/=mype) then
                      call nbprocs_info%add_to_f(child%node%ipe, igrid, inc1, inc2, inc3)
                      nrecv_bc_r=nrecv_bc_r+1
                      nsend_bc_p=nsend_bc_p+1
                      nbuff_bc_send_p=nbuff_bc_send_p+sizes_p_send_total(inc1,&
                           inc2,inc3)
                      nbuff_bc_recv_r=nbuff_bc_recv_r+sizes_r_recv_total(inc1,&
                           inc2,inc3)
                   end if
                end do
                end do
                end do
             end select

             ! flux fix for conservation only for pure directional shifts
             if (abs(i1)+abs(i2)+abs(i3)==1) then
                if (i1/=0) then
                   idim=1
                   iside=int((i1+3)/2)
                end if
                if (i2/=0) then
                   idim=2
                   iside=int((i2+3)/2)
                end if
                if (i3/=0) then
                   idim=3
                   iside=int((i3+3)/2)
                end if
                select case (my_neighbor_type)
                ! only across fine-coarse or coarse-fine boundaries
                case (neighbor_coarse)
                   if (my_neighbor%node%ipe/=mype) then
                      if (.not.pole(idim)) nsend_fc(idim)=nsend_fc(idim)+1
                   end if
                case (neighbor_fine)
                   if (pole(idim)) then
                      icdim=iside
                   else
                      icdim=3-iside
                   end if
                   select case (idim)
                   case (1)
                      do ic1=icdim,icdim
                   do ic2=1,2
                   do ic3=1,2
                         child%node => my_neighbor%node%child(ic1,ic2,&
                            ic3)%node
                         if (child%node%ipe/=mype) then
                            if (.not.pole(1)) nrecv_fc(1)=nrecv_fc(1)+1
                         end if
                      end do
                   end do
                   end do
                   case (2)
                      do ic1=1,2
                   do ic2=icdim,icdim
                   do ic3=1,2
                         child%node => my_neighbor%node%child(ic1,ic2,&
                            ic3)%node
                         if (child%node%ipe/=mype) then
                            if (.not.pole(2)) nrecv_fc(2)=nrecv_fc(2)+1
                         end if
                      end do
                   end do
                   end do
                   case (3)
                      do ic1=1,2
                   do ic2=1,2
                   do ic3=icdim,icdim
                         child%node => my_neighbor%node%child(ic1,ic2,&
                            ic3)%node
                         if (child%node%ipe/=mype) then
                            if (.not.pole(3)) nrecv_fc(3)=nrecv_fc(3)+1
                         end if
                      end do
                   end do
                   end do
                   end select
                end select
             end if

             if (phi_ > 0) then
               neighbor_pole(i1,i2,i3,igrid)=0
               if (my_neighbor_type>1) then
                 do idim=1,3
                   if (pole(idim)) then
                     neighbor_pole(i1,i2,i3,igrid)=idim
                     exit ! there can only be one pole between two meshes
                   end if
                 end do
               end if
             end if
             neighbor_type(i1,i2,i3,igrid)=my_neighbor_type

          end if
       end do
       end do
    end do

    ! Now all the neighbour information is known.
    ! already fill the idphyb structure
    call identifyphysbound_connectivity(igrid)


    if(stagger_grid) then
       !Check if there are special corners that need to be communicated
       !To determine whether to send/receive, we must check three neighbours
        do i3=-1,1
        do i2=-1,1
        do i1=-1,1
           if (abs(i1)+abs(i2)+abs(i3)==1) then
             if (neighbor_pole(i1,i2,i3,igrid)/=0) cycle
            ! Assign value to idim and iside
            if (i1/=0) then
                idim=1
                iside=int((i1+3)/2)
             end if
            if (i2/=0) then
                idim=2
                iside=int((i2+3)/2)
             end if
            if (i3/=0) then
                idim=3
                iside=int((i3+3)/2)
             end if
             ! Fine block surrounded by coarse blocks
             if (neighbor_type(i1,i2,i3,igrid)==2) then
               do idir=idim+1,ndim
                 pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
                 mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
                 ph1=pi1-kr(idim,1)*(2*iside-3)
                 ph2=pi2-kr(idim,2)*(2*iside-3)
                 ph3=pi3-kr(idim,3)*(2*iside-3);
                 mh1=mi1-kr(idim,1)*(2*iside-3)
                 mh2=mi2-kr(idim,2)*(2*iside-3)
                 mh3=mi3-kr(idim,3)*(2*iside-3);

                 if (neighbor_type(pi1,pi2,pi3,igrid)==2.and.neighbor_type(ph1,&
                    ph2,ph3,igrid)==2.and.mype/=neighbor(2,pi1,pi2,pi3,&
                    igrid).and.neighbor_pole(pi1,pi2,pi3,igrid)==0) then
                    nsend_cc(idim) = nsend_cc(idim) + 1
                 end if

                 if (neighbor_type(mi1,mi2,mi3,igrid)==2.and.neighbor_type(mh1,&
                    mh2,mh3,igrid)==2.and.mype/=neighbor(2,mi1,mi2,mi3,&
                    igrid).and.neighbor_pole(mi1,mi2,mi3,igrid)==0) then
                    nsend_cc(idim) = nsend_cc(idim) + 1
                 end if
               end do
             end if
             ! Coarse block diagonal to fine block(s)
             if (neighbor_type(i1,i2,i3,igrid)==3) then
               do idir=idim+1,ndim
                 pi1=i1+kr(idir,1);pi2=i2+kr(idir,2);pi3=i3+kr(idir,3);
                 mi1=i1-kr(idir,1);mi2=i2-kr(idir,2);mi3=i3-kr(idir,3);
                 ph1=pi1-kr(idim,1)*(2*iside-3)
                 ph2=pi2-kr(idim,2)*(2*iside-3)
                 ph3=pi3-kr(idim,3)*(2*iside-3);
                 mh1=mi1-kr(idim,1)*(2*iside-3)
                 mh2=mi2-kr(idim,2)*(2*iside-3)
                 mh3=mi3-kr(idim,3)*(2*iside-3);

                 if (neighbor_type(pi1,pi2,pi3,igrid)==4.and.neighbor_type(ph1,&
                    ph2,ph3,igrid)==3.and.neighbor_pole(pi1,pi2,pi3,&
                    igrid)==0) then
                  ! Loop on children (several in 3D)
                  do ic3=1+int((1-pi3)/2),2-int((1+pi3)/2)
                      inc3=2*pi3+ic3
                  do ic2=1+int((1-pi2)/2),2-int((1+pi2)/2)
                      inc2=2*pi2+ic2
                  do ic1=1+int((1-pi1)/2),2-int((1+pi1)/2)
                      inc1=2*pi1+ic1
                     if (mype.ne.neighbor_child(2,inc1,inc2,inc3,igrid)) then
                       nrecv_cc(idim) = nrecv_cc(idim) + 1
                     end if
                  end do
                  end do
                  end do
                 end if

                 if (neighbor_type(mi1,mi2,mi3,igrid)==4.and.neighbor_type(mh1,&
                    mh2,mh3,igrid)==3.and.neighbor_pole(mi1,mi2,mi3,&
                    igrid)==0) then
                  ! Loop on children (several in 3D)
                  do ic3=1+int((1-mi3)/2),2-int((1+mi3)/2)
                      inc3=2*mi3+ic3
                  do ic2=1+int((1-mi2)/2),2-int((1+mi2)/2)
                      inc2=2*mi2+ic2
                  do ic1=1+int((1-mi1)/2),2-int((1+mi1)/2)
                      inc1=2*mi1+ic1
                     if (mype.ne.neighbor_child(2,inc1,inc2,inc3,igrid)) then
                       nrecv_cc(idim) = nrecv_cc(idim) + 1
                     end if
                  end do
                  end do
                  end do
                 end if
               end do
             end if
           end if
        end do
        end do
        end do
       end if

    end do

    ! allocate with new nbstructure
    call nbprocs_info%alloc_buffers_srl(nwgc, &
         ixS_srl_min1, ixS_srl_max1, &
         ixS_srl_min2, ixS_srl_max2, &
         ixS_srl_min3, ixS_srl_max3, &
         ixR_srl_min1, ixR_srl_max1, &
         ixR_srl_min2, ixR_srl_max2, &
         ixR_srl_min3, ixR_srl_max3 &
         )
    call nbprocs_info%alloc_buffers_f(nwgc, &
         ixS_p_min1, ixS_p_max1, &
         ixS_p_min2, ixS_p_max2, &
         ixS_p_min3, ixS_p_max3, &
         ixR_r_min1, ixR_r_max1, &
         ixR_r_min2, ixR_r_max2, &
         ixR_r_min3, ixR_r_max3 &
         )
    call nbprocs_info%alloc_buffers_c(nwgc, &
         ixS_r_min1, ixS_r_max1, &
         ixS_r_min2, ixS_r_max2, &
         ixS_r_min3, ixS_r_max3, &
         ixR_p_min1, ixR_p_max1, &
         ixR_p_min2, ixR_p_max2, &
         ixR_p_min3, ixR_p_max3 &
         )

    ! allocate nbstructure srl requests and status
    ! SRL:
    if (allocated(recvstatus_srl_nb)) then
       deallocate(recv_srl_nb, recvstatus_srl_nb)
    end if
    allocate(recv_srl_nb(nbprocs_info%nbprocs_srl*2), recvstatus_srl_nb(MPI_STATUS_SIZE,nbprocs_info%nbprocs_srl*2))
    if (allocated(sendstatus_srl_nb)) then
       deallocate(send_srl_nb, sendstatus_srl_nb)
    end if
    allocate(send_srl_nb(nbprocs_info%nbprocs_srl*2), sendstatus_srl_nb(MPI_STATUS_SIZE,nbprocs_info%nbprocs_srl*2))

    ! F:
    if (allocated(recvstatus_f_nb)) then
       deallocate(recv_f_nb, recvstatus_f_nb)
    end if
    allocate(recv_f_nb(nbprocs_info%nbprocs_f*2), recvstatus_f_nb(MPI_STATUS_SIZE,nbprocs_info%nbprocs_f*2))
    if (allocated(sendstatus_f_nb)) then
       deallocate(send_f_nb, sendstatus_f_nb)
    end if
    allocate(send_f_nb(nbprocs_info%nbprocs_f*2), sendstatus_f_nb(MPI_STATUS_SIZE,nbprocs_info%nbprocs_f*2))

    ! C:
    if (allocated(recvstatus_c_nb)) then
       deallocate(recv_c_nb, recvstatus_c_nb)
    end if
    allocate(recv_c_nb(nbprocs_info%nbprocs_c*2), recvstatus_c_nb(MPI_STATUS_SIZE,nbprocs_info%nbprocs_c*2))
    if (allocated(sendstatus_c_nb)) then
       deallocate(send_c_nb, sendstatus_c_nb)
    end if
    allocate(send_c_nb(nbprocs_info%nbprocs_c*2), sendstatus_c_nb(MPI_STATUS_SIZE,nbprocs_info%nbprocs_c*2))

    
    ! allocate space for mpi recieve for siblings and restrict ghost cell filling
    nrecvs=nrecv_bc_srl+nrecv_bc_r
    if (allocated(recvstatus_c_sr)) then
      deallocate(recvstatus_c_sr,recvrequest_c_sr)
      allocate(recvstatus_c_sr(MPI_STATUS_SIZE,nrecvs),&
         recvrequest_c_sr(nrecvs))
    else
      allocate(recvstatus_c_sr(MPI_STATUS_SIZE,nrecvs),&
         recvrequest_c_sr(nrecvs))
    end if
    recvrequest_c_sr=MPI_REQUEST_NULL

    ! allocate space for mpi send for siblings and restrict ghost cell filling
    nsends=nsend_bc_srl+nsend_bc_r
    if (allocated(sendstatus_c_sr)) then
      deallocate(sendstatus_c_sr,sendrequest_c_sr)
      allocate(sendstatus_c_sr(MPI_STATUS_SIZE,nsends),&
         sendrequest_c_sr(nsends))
    else
      allocate(sendstatus_c_sr(MPI_STATUS_SIZE,nsends),&
         sendrequest_c_sr(nsends))
    end if
    sendrequest_c_sr=MPI_REQUEST_NULL

    ! allocate space for mpi recieve for prolongation ghost cell filling
    if (allocated(recvstatus_c_p)) then
      deallocate(recvstatus_c_p,recvrequest_c_p)
      allocate(recvstatus_c_p(MPI_STATUS_SIZE,nrecv_bc_p),&
         recvrequest_c_p(nrecv_bc_p))
    else
      allocate(recvstatus_c_p(MPI_STATUS_SIZE,nrecv_bc_p),&
         recvrequest_c_p(nrecv_bc_p))
    end if
    recvrequest_c_p=MPI_REQUEST_NULL

    ! allocate space for mpi send for prolongation ghost cell filling
    if (allocated(sendstatus_c_p)) then
      deallocate(sendstatus_c_p,sendrequest_c_p)
      allocate(sendstatus_c_p(MPI_STATUS_SIZE,nsend_bc_p),&
         sendrequest_c_p(nsend_bc_p))
    else
      allocate(sendstatus_c_p(MPI_STATUS_SIZE,nsend_bc_p),&
         sendrequest_c_p(nsend_bc_p))
    end if
    sendrequest_c_p=MPI_REQUEST_NULL

    if(stagger_grid) then
      ! allocate space for recieve buffer for siblings ghost cell filling
      if (allocated(recvbuffer_srl)) then
        if (nbuff_bc_recv_srl /= size(recvbuffer_srl)) then
          deallocate(recvbuffer_srl)
          allocate(recvbuffer_srl(nbuff_bc_recv_srl))
        end if
      else
        allocate(recvbuffer_srl(nbuff_bc_recv_srl))
      end if
      if (allocated(recvstatus_srl)) then
        deallocate(recvstatus_srl,recvrequest_srl)
        allocate(recvstatus_srl(MPI_STATUS_SIZE,nrecv_bc_srl),&
           recvrequest_srl(nrecv_bc_srl))
      else
        allocate(recvstatus_srl(MPI_STATUS_SIZE,nrecv_bc_srl),&
           recvrequest_srl(nrecv_bc_srl))
      end if
      recvrequest_srl=MPI_REQUEST_NULL

      ! allocate space for send buffer for siblings ghost cell filling
      if (allocated(sendbuffer_srl)) then
        if (nbuff_bc_send_srl /= size(sendbuffer_srl)) then
          deallocate(sendbuffer_srl)
          allocate(sendbuffer_srl(nbuff_bc_send_srl))
        end if
      else
        allocate(sendbuffer_srl(nbuff_bc_send_srl))
      end if
      if (allocated(sendstatus_srl)) then
        deallocate(sendstatus_srl,sendrequest_srl)
        allocate(sendstatus_srl(MPI_STATUS_SIZE,nsend_bc_srl),&
           sendrequest_srl(nsend_bc_srl))
      else
        allocate(sendstatus_srl(MPI_STATUS_SIZE,nsend_bc_srl),&
           sendrequest_srl(nsend_bc_srl))
      end if
      sendrequest_srl=MPI_REQUEST_NULL

      ! allocate space for recieve buffer for restrict ghost cell filling
      if (allocated(recvbuffer_r)) then
        if (nbuff_bc_recv_r /= size(recvbuffer_r)) then
          deallocate(recvbuffer_r)
          allocate(recvbuffer_r(nbuff_bc_recv_r))
        end if
      else
        allocate(recvbuffer_r(nbuff_bc_recv_r))
      end if
      if (allocated(recvstatus_r)) then
        deallocate(recvstatus_r,recvrequest_r)
        allocate(recvstatus_r(MPI_STATUS_SIZE,nrecv_bc_r),&
           recvrequest_r(nrecv_bc_r))
      else
        allocate(recvstatus_r(MPI_STATUS_SIZE,nrecv_bc_r),&
           recvrequest_r(nrecv_bc_r))
      end if
      recvrequest_r=MPI_REQUEST_NULL

      ! allocate space for send buffer for restrict ghost cell filling
      if (allocated(sendbuffer_r)) then
        if (nbuff_bc_send_r /= size(sendbuffer_r)) then
          deallocate(sendbuffer_r)
          allocate(sendbuffer_r(nbuff_bc_send_r))
        end if
      else
        allocate(sendbuffer_r(nbuff_bc_send_r))
      end if
      if (allocated(sendstatus_r)) then
        deallocate(sendstatus_r,sendrequest_r)
        allocate(sendstatus_r(MPI_STATUS_SIZE,nsend_bc_r),&
           sendrequest_r(nsend_bc_r))
      else
        allocate(sendstatus_r(MPI_STATUS_SIZE,nsend_bc_r),&
           sendrequest_r(nsend_bc_r))
      end if
      sendrequest_r=MPI_REQUEST_NULL

      ! allocate space for recieve buffer for prolong ghost cell filling
      if (allocated(recvbuffer_p)) then
        if (nbuff_bc_recv_p /= size(recvbuffer_p)) then
          deallocate(recvbuffer_p)
          allocate(recvbuffer_p(nbuff_bc_recv_p))
        end if
      else
        allocate(recvbuffer_p(nbuff_bc_recv_p))
      end if
      if (allocated(recvstatus_p)) then
        deallocate(recvstatus_p,recvrequest_p)
        allocate(recvstatus_p(MPI_STATUS_SIZE,nrecv_bc_p),&
           recvrequest_p(nrecv_bc_p))
      else
        allocate(recvstatus_p(MPI_STATUS_SIZE,nrecv_bc_p),&
           recvrequest_p(nrecv_bc_p))
      end if
      recvrequest_p=MPI_REQUEST_NULL

      ! allocate space for send buffer for restrict ghost cell filling
      if (allocated(sendbuffer_p)) then
        if (nbuff_bc_send_p /= size(sendbuffer_p)) then
          deallocate(sendbuffer_p)
          allocate(sendbuffer_p(nbuff_bc_send_p))
        end if
      else
        allocate(sendbuffer_p(nbuff_bc_send_p))
      end if
      if (allocated(sendstatus_p)) then
        deallocate(sendstatus_p,sendrequest_p)
        allocate(sendstatus_p(MPI_STATUS_SIZE,nsend_bc_p),&
           sendrequest_p(nsend_bc_p))
      else
        allocate(sendstatus_p(MPI_STATUS_SIZE,nsend_bc_p),&
           sendrequest_p(nsend_bc_p))
      end if
      sendrequest_p=MPI_REQUEST_NULL
    end if

    
    !update the neighbor information on the device
    !$acc update device(neighbor, neighbor_type, neighbor_pole, neighbor_child, idphyb)
    call nbprocs_update_device


  end subroutine build_connectivity

  subroutine identifyphysbound_connectivity(igrid)
    use mod_global_parameters
    integer, intent(in) :: igrid

    if (neighbor_type(-1,0,0,igrid) == neighbor_boundary .and. neighbor_type(+1,0,0,igrid) == neighbor_boundary ) then
       idphyb(1,igrid) = 2
    else if (neighbor_type(-1,0,0,igrid) == neighbor_boundary) then
       idphyb(1,igrid) = -1
    else if (neighbor_type(+1,0,0,igrid) == neighbor_boundary) then
       idphyb(1,igrid) = +1
    else
       idphyb(1,igrid) = 0
    end if

    if (neighbor_type(0,-1,0,igrid) == neighbor_boundary .and. neighbor_type(0,+1,0,igrid) == neighbor_boundary ) then
       idphyb(2,igrid) = 2
    else if (neighbor_type(0,-1,0,igrid) == neighbor_boundary) then
       idphyb(2,igrid) = -1
    else if (neighbor_type(0,+1,0,igrid) == neighbor_boundary) then
       idphyb(2,igrid) = +1
    else
       idphyb(2,igrid) = 0
    end if

    if (neighbor_type(0,0,-1,igrid) == neighbor_boundary .and. neighbor_type(0,0,+1,igrid) == neighbor_boundary ) then
       idphyb(3,igrid) = 2
    else if (neighbor_type(0,0,-1,igrid) == neighbor_boundary) then
       idphyb(3,igrid) = -1
    else if (neighbor_type(0,0,+1,igrid) == neighbor_boundary) then
       idphyb(3,igrid) = +1
    else
       idphyb(3,igrid) = 0
    end if

  end subroutine identifyphysbound_connectivity

end module mod_functions_connectivity
