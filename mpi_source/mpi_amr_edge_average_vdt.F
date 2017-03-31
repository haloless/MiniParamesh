!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_edge_average_vdt
!! NAME
!!
!!   amr_edge_average_vdt
!!
!! SYNOPSIS
!!
!!   call amr_edge_average_vdt(mype, nsub)
!!
!!   call amr_edge_average_vdt(integer, integer)
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: mype
!!     The calling processors id.
!!
!!   integer, intent(in) :: nsub
!!     ???
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   paramesh_interfaces
!!   paramesh_mpi_interfaces
!!
!! CALLS
!!
!!   amr_restrict_edge_data_vdt
!!   mpi_amr_comm_setup
!!   mpi_put_edge_buffer_1blk
!!
!!
!! RETURNS
!!
!!   Nothing returned.  Upon return the edge data at refinement jumps is
!!   properly averaged.
!!
!! DESCRIPTION
!! 
!!   This routine gets cell edge-based data at block boundaries from 
!!   neighbors who are parents of leaf blocks. 
!!
!!   The data structure used to store and pass this data is defined
!!   in the module 'physicaldata'.
!!
!!   This version is used when variable timesteps are allowed across the
!!   blocks in the computation. It is called by the wrapper routine
!!   'amr_edge_average'.
!! 
!! AUTHORS
!!
!!  Written by Peter MacNeice (July 1997).
!!
!!***

#include "paramesh_preprocessor.fh"


      subroutine amr_edge_average_vdt(mype,nsub)


      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_interfaces, only : amr_restrict_edge_data_vdt
#ifdef MPI_USED
      use paramesh_mpi_interfaces, only : mpi_amr_comm_setup,
     .                                    mpi_put_edge_buffer_1blk
#endif

      implicit none

#include "amr_shmem.fh"

      integer, intent(in)  ::  mype,nsub

!------------------------------------
! local variables

      integer :: remote_pe,remote_block
      integer :: remote_pe2,remote_block2
      integer,save :: anodetype(1)
      integer,save :: cnodetype

      integer :: tag_offset,nprocs
      integer, external :: shmem_n_pes
      logical :: lfound

      logical :: lcc,lfc,lec,lnc
      logical :: lfulltree, lrestrict, lprolong, ledge, lflux
      logical :: lguard
      integer :: iopt, lb, lcycle, phase0, phase1
      integer :: iface

      integer :: nguard0
      integer :: ng_off
      integer :: kup,klo,kup1
      integer :: jf, iblk

      integer :: ierrorcode

!------------------------------------

      if (var_dt) then

      nguard0 = nguard*npgs
      ng_off = nguard0+iface_off

      klo  = 1+k3d*nguard0
      kup  = 1+k3d*(nzb+nguard0-1)
      kup1 = k3d+nzb+k3d*nguard0

      if (ndim >= 2) then

      lcc = .false.
      lfc = .false.
      lec = .false.
      lnc = .false.
      iopt = 1

#ifdef MPI_USED
      nprocs = shmem_n_pes()
      tag_offset = 100
!      call mpi_morton_bnd_restrict
!     .             (mype,nprocs,tag_offset)
      tag_offset = 100
!      call mpi_amr_edgeaver_setup
!     .             (mype,nprocs,tag_offset)


! Note, both ledge and lrestrict are true so that the fluxes
! are acquired which are needed in the restriction operation.
      lguard    = .false.
      lprolong  = .false.
      lflux     = .false.
      ledge     = .true.
      lrestrict = .true.
      lfulltree = .false.
      call mpi_amr_comm_setup(mype,nprocs,lguard,lprolong,
     .                        lflux,ledge,lrestrict,lfulltree,
     .                        iopt,lcc,lfc,lec,lnc,tag_offset)

#endif /* MPI_USED */



        if(lnblocks.gt.0) then
        do lb = 1,lnblocks

! Is this a parent of at least one leaf block ?
      if(nodetype(lb).eq.2) then


! Set timestep phases for the current block, and for the next finer level.
        lcycle = loc_cycle(lrefine(lb))
        phase0 = phase_dt(lrefine(lb))
        phase1 = phase_dt(lrefine(lb)+1)

! At start of the current blocks timestep zero out the arrays used to
! accumulate boundary edge data from its children.
        if(lcycle.eq.1) then
           ttbedge_facex_y(:,:,:,:,lb) = 0.
           ttbedge_facey_x(:,:,:,:,lb) = 0.
           if(ndim.eq.3) then
             ttbedge_facex_z(:,:,:,:,lb) = 0.
             ttbedge_facey_z(:,:,:,:,lb) = 0.
             ttbedge_facez_x(:,:,:,:,lb) = 0.
             ttbedge_facez_y(:,:,:,:,lb) = 0.
           endif
        endif

      endif
      enddo
      endif
!------------------------------------

! Leaf blocks which have completed their timestep provide reduced
! boundary edge data to their parents.
! Edge values are accumulated in the ttbedge_face arrays.
      call amr_restrict_edge_data_vdt(mype)

#ifdef MPI_USED
      tag_offset = 100
!      call mpi_morton_bnd_fluxcon
!     .             (mype,nprocs,tag_offset)
      tag_offset = 100
!      call mpi_amr_edgeaver_setup
!     .             (mype,nprocs,tag_offset)


      lguard    = .false.
      lprolong  = .false.
      lflux     = .false.
      ledge     = .true.
      lrestrict = .false.
      lfulltree = .false.
      call mpi_amr_comm_setup(mype,nprocs,lguard,lprolong,
     .                        lflux,ledge,lrestrict,lfulltree,
     .                        iopt,lcc,lfc,lec,lnc,tag_offset)

#endif /* MPI_USED */
!------------------------------------

! Parents who have completed their timestep and border a leaf block
! update their edges.
      do lb = 1,lnblocks


! Is this a parent block of at least one leaf node?
      if((nodetype(lb).eq.2).and.ldtcomplete(lb)) then

! If yes then cycle through its neighbors.
        do iface=1,nfaces

! If this neighbor is a leaf block or an external boundary then
! replace fluxes with restricted fluxes.
          cnodetype = 1
          if(neigh(1,iface,lb).ge.1) then
            remote_pe    = neigh(2,iface,lb)
            remote_block = neigh(1,iface,lb)

            call shmem_integer_get(anodetype,nodetype(remote_block),
     .                             1,remote_pe)
            cnodetype = anodetype(1)
          endif
          if(cnodetype.eq.1) then

            if(iface.eq.1) then

              bedge_facex_y(:,1,:,:,lb)=ttbedge_facex_y(:,1,:,:,lb)

              if((ndim.eq.3).or.(l2p5d.eq.1)) 
     .        bedge_facex_z(:,1,:,:,lb)=ttbedge_facex_z(:,1,:,:,lb)

! make common variables on an edge consistent
              bedge_facey_z(:,1+nguard0,1,klo:kup,lb) =
     .               bedge_facex_z(:,1,1+nguard0,klo:kup,lb)

              bedge_facey_z(:,1+nguard0,2,klo:kup,lb) =
     .               bedge_facex_z(:,1,1+nyb+nguard0,klo:kup,lb)

              bedge_facez_y(:,1+nguard0,1+nguard0:nyb+nguard0,1,lb)
     .             = bedge_facex_y(:,1,1+nguard0:nyb+nguard0,klo,lb)

              if(ndim.eq.3) then
              bedge_facez_y(:,1+nguard0,1+nguard0:nyb+nguard0,
     .                                                      1+k3d,lb)
     .             = bedge_facex_y(:,1,1+nguard0:nyb+nguard0,kup1,lb)
              endif

            elseif(iface.eq.2) then

              bedge_facex_y(:,2,:,:,lb)=ttbedge_facex_y(:,2,:,:,lb)

              if((ndim.eq.3).or.(l2p5d.eq.1)) 
     .        bedge_facex_z(:,2,:,:,lb)=ttbedge_facex_z(:,2,:,:,lb)

! make common variables on an edge consistent
              bedge_facey_z(:,1+nxb+nguard0,1,klo:kup,lb) =
     .            bedge_facex_z(:,2,1+nguard0,klo:kup,lb)

              bedge_facey_z(:,1+nxb+nguard0,2,klo:kup,lb) =
     .            bedge_facex_z(:,2,1+nyb+nguard0,klo:kup,lb)

              bedge_facez_y(:,1+nxb+nguard0,1+nguard0:nyb+nguard0,
     .                                                       1,lb)=
     .            bedge_facex_y(:,2,1+nguard0:nyb+nguard0,klo,lb)

              if(ndim.eq.3) then
              bedge_facez_y(:,1+nxb+nguard0,1+nguard0:nyb+nguard0,
     .                                                   1+k3d,lb)=
     .            bedge_facex_y(:,2,1+nguard0:nyb+nguard0,kup1,lb)
              endif

            elseif(iface.eq.3) then

              bedge_facey_x(:,:,1,:,lb)=ttbedge_facey_x(:,:,1,:,lb)

              if((ndim.eq.3).or.(l2p5d.eq.1)) 
     .        bedge_facey_z(:,:,1,:,lb)=ttbedge_facey_z(:,:,1,:,lb)

! make common variables on an edge consistent
              bedge_facex_z(:,1,1+nguard0,klo:kup,lb) =
     .          bedge_facey_z(:,1+nguard0,1,klo:kup,lb)

              bedge_facex_z(:,2,1+nguard0,klo:kup,lb) =
     .          bedge_facey_z(:,1+nxb+nguard0,1,klo:kup,lb)

              bedge_facez_x(:,1+nguard0:nxb+nguard0,1+nguard0,1,lb)=
     .          bedge_facey_x(:,1+nguard0:nxb+nguard0,1,klo,lb)

              if(ndim.eq.3) then
              bedge_facez_x(:,1+nguard0:nxb+nguard0,1+nguard0,
     .                                                    1+k3d,lb)=
     .          bedge_facey_x(:,1+nguard0:nxb+nguard0,1,kup1,lb)
              endif

            elseif(iface.eq.4) then

              bedge_facey_x(:,:,2,:,lb)=ttbedge_facey_x(:,:,2,:,lb)

              if((ndim.eq.3).or.(l2p5d.eq.1)) 
     .        bedge_facey_z(:,:,2,:,lb)=ttbedge_facey_z(:,:,2,:,lb)

! make common variables on an edge consistent
               bedge_facex_z(:,1,1+nyb+nguard0,klo:kup,lb) =
     .              bedge_facey_z(:,1+nguard0,2,klo:kup,lb)

              bedge_facex_z(:,2,1+nyb+nguard0,klo:kup,lb) =
     .              bedge_facey_z(:,1+nxb+nguard0,2,klo:kup,lb)

              bedge_facez_x(:,1+nguard0:nxb+nguard0,1+nyb+nguard0,
     .                                                       1,lb)=
     .              bedge_facey_x(:,1+nguard0:nxb+nguard0,2,klo,lb)

              if(ndim.eq.3) then
              bedge_facez_x(:,1+nguard0:nxb+nguard0,1+nyb+nguard0,
     .                                                   1+k3d,lb)=
     .              bedge_facey_x(:,1+nguard0:nxb+nguard0,2,kup1,lb)
              endif

            elseif(iface.eq.5) then

              bedge_facez_x(:,:,:,1,lb)=ttbedge_facez_x(:,:,:,1,lb)
              bedge_facez_y(:,:,:,1,lb)=ttbedge_facez_y(:,:,:,1,lb)

! make common variables on an edge consistent
              bedge_facey_x(:,1+nguard0:nxb+nguard0,1,klo,lb)=
     .          bedge_facez_x(:,1+nguard0:nxb+nguard0,1+nguard0,
     .                                                     1,lb)

              bedge_facey_x(:,1+nguard0:nxb+nguard0,2,klo,lb)=
     .          bedge_facez_x(:,1+nguard0:nxb+nguard0,1+nyb+nguard0
     .                                                     ,1,lb)

              bedge_facex_y(:,1,1+nguard0:nyb+nguard0,klo,lb)=
     .          bedge_facez_y(:,1+nguard0,1+nguard0:nyb+nguard0,
     .                                                     1,lb)

              bedge_facex_y(:,2,1+nguard0:nyb+nguard0,klo,lb)=
     .          bedge_facez_y(:,1+nxb+nguard0,1+nguard0:nyb+nguard0
     .                                                    ,1,lb)

            elseif(iface.eq.6) then

              bedge_facez_x(:,:,:,2,lb)=ttbedge_facez_x(:,:,:,2,lb)
              bedge_facez_y(:,:,:,2,lb)=ttbedge_facez_y(:,:,:,2,lb)

! make common variables on an edge consistent
              bedge_facey_x(:,1+nguard0:nxb+nguard0,1,kup1,lb)=
     .          bedge_facez_x(:,1+nguard0:nxb+nguard0,1+nguard0,2,lb)

              bedge_facey_x(:,1+nguard0:nxb+nguard0,2,kup1,lb)=
     .          bedge_facez_x(:,1+nguard0:nxb+nguard0,1+nyb+nguard0,
     .                                                        2,lb)

              bedge_facex_y(:,1,1+nguard0:nyb+nguard0,kup1,lb)=
     .          bedge_facez_y(:,1+nguard0,1+nguard0:nyb+nguard0,2,lb)

              bedge_facex_y(:,2,1+nguard0:nyb+nguard0,kup1,lb)=
     .          bedge_facez_y(:,1+nxb+nguard0,1+nguard0:nyb+nguard0,
     .                                                          2,lb)

            endif
          endif
        enddo
      endif
      enddo

!------------------------------------

! cycle through the grid blocks on this processor
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks

! Is this a leaf block and not at the original refinement level ?
!      if(nodetype(lb).eq.1.and.lrefine(lb).gt.1) then
      if(nodetype(lb).eq.1) then

! Has this block completed its timestep?
      if(ldtcomplete(lb)) then

! Cycle over the blocks faces
       do jf = 1,nfaces

          remote_pe = neigh(2,jf,lb)
          remote_block  = neigh(1,jf,lb)
          remote_pe2 = neigh(2,jf,lb)
          remote_block2  = neigh(1,jf,lb)
          cnodetype = 0
          lfound = .false.

          if(remote_block.gt.0) then

#ifdef MPI_USED
! (remote_block,remote_pe) may be a local block, a remote block,
! or it may not exist.
! If it is a local block then check its nodetype.
! If it is found in the list of remote blocks stored in buffer space
! then check its nodetype.
! If it is not found in either of these places, then set its nodetype
! to 0.
         if(remote_pe2.ne.mype) then

            lfound = .false.
            do iblk = strt_buffer,last_buffer
               if(remote_block2.eq.laddress(1,iblk).and.
     .              remote_pe2 .eq.laddress(2,iblk) ) then
                  remote_block2 = iblk
                  remote_pe2    = mype
                  lfound = .true.
               endif
            enddo

         elseif(remote_pe2.eq.mype) then

            lfound = .true.

         endif
#else  /* MPI_USED */
         lfound = .true.
#endif /* MPI_USED */

! Is the neighbor to this face a parent of a leaf block?
          if(lfound) then
             call shmem_integer_get(anodetype,
     .                     nodetype(remote_block2),1,remote_pe2)
             cnodetype = anodetype(1)
          endif

          endif ! end if (remote_block

       if(cnodetype.eq.2) then

          if (remote_pe == mype .and. remote_block <= lnblocks) then

          if(jf.eq.1) then

             call shmem_real_get(recvarx1e(1,1,1,1),
     .            bedge_facex_y(1,1,1,1,remote_block),
     .            len_block_ex*nedges,remote_pe)
             bedge_facex_y(:,1,:,:,lb) = recvarx1e(:,2,:,:)

             if((ndim.eq.3).or.(l2p5d.eq.1)) then
                call shmem_real_get(recvarx1e(1,1,1,1),
     .               bedge_facex_z(1,1,1,1,remote_block),
     .               len_block_ex*nedges,remote_pe)
                bedge_facex_z(:,1,:,:,lb) = recvarx1e(:,2,:,:)
             endif

          elseif(jf.eq.2) then

             call shmem_real_get(recvarx1e(1,1,1,1),
     .            bedge_facex_y(1,1,1,1,remote_block),
     .            len_block_ex*nedges,remote_pe)
             bedge_facex_y(:,2,:,:,lb) = recvarx1e(:,1,:,:)
             
             if((ndim.eq.3).or.(l2p5d.eq.1)) then
                call shmem_real_get(recvarx1e(1,1,1,1),
     .               bedge_facex_z(1,1,1,1,remote_block),
     .               len_block_ex*nedges,remote_pe)
                bedge_facex_z(:,2,:,:,lb) = recvarx1e(:,1,:,:)
             endif

          elseif(jf.eq.3) then

             if((ndim.eq.3).or.(l2p5d.eq.1)) then
                call shmem_real_get(recvary1e(1,1,1,1),
     .               bedge_facey_z(1,1,1,1,remote_block),
     .               len_block_ey*nedges,remote_pe)
                bedge_facey_z(:,:,1,:,lb) = recvary1e(:,:,2,:)
             endif

             call shmem_real_get(recvary1e(1,1,1,1),
     .            bedge_facey_x(1,1,1,1,remote_block),
     .            len_block_ey*nedges,remote_pe)
             bedge_facey_x(:,:,1,:,lb) = recvary1e(:,:,2,:)
             
          elseif(jf.eq.4) then

            if((ndim.eq.3).or.(l2p5d.eq.1)) then
              call shmem_real_get(recvary1e(1,1,1,1),
     .           bedge_facey_z(1,1,1,1,remote_block),
     .           len_block_ey*nedges,remote_pe)
              bedge_facey_z(:,:,2,:,lb) = recvary1e(:,:,1,:)
            endif

            call shmem_real_get(recvary1e(1,1,1,1),
     .           bedge_facey_x(1,1,1,1,remote_block),
     .           len_block_ey*nedges,remote_pe)
            bedge_facey_x(:,:,2,:,lb) = recvary1e(:,:,1,:)

          elseif(jf.eq.5) then

             call shmem_real_get(recvarz1e(1,1,1,1),
     .            bedge_facez_x(1,1,1,1,remote_block),
     .            len_block_ez*nedges,remote_pe)
             bedge_facez_x(:,:,:,1,lb) = recvarz1e(:,:,:,2)
             call shmem_real_get(recvarz1e(1,1,1,1),
     .            bedge_facez_y(1,1,1,1,remote_block),
     .            len_block_ez*nedges,remote_pe)
             bedge_facez_y(:,:,:,1,lb) = recvarz1e(:,:,:,2) 
             
          elseif(jf.eq.6) then 

             call shmem_real_get(recvarz1e(1,1,1,1),
     .            bedge_facez_x(1,1,1,1,remote_block),
     .            len_block_ez*nedges,remote_pe)
             bedge_facez_x(:,:,:,2,lb) = recvarz1e(:,:,:,1)
            call shmem_real_get(recvarz1e(1,1,1,1),
     .            bedge_facez_y(1,1,1,1,remote_block),
     .            len_block_ez*nedges,remote_pe)
            bedge_facez_y(:,:,:,2,lb) = recvarz1e(:,:,:,1)

          endif


          else                     ! if (remote_pe


            call mpi_put_edge_buffer_1blk(lb,remote_block,remote_pe)

            if(jf == 1) then
               bedge_facex_y(:,1,:,:,lb) = recvarx1e(:,2,:,:)
            elseif(jf == 2) then
               bedge_facex_y(:,2,:,:,lb) = recvarx1e(:,1,:,:)
            end if

            if(jf == 1) then
               if((ndim.eq.3).or.(l2p5d.eq.1)) then
                  bedge_facex_z(:,1,:,:,lb) = recvarx2e(:,2,:,:)
               endif
            elseif(jf == 2) then
               if((ndim.eq.3).or.(l2p5d.eq.1)) then
                  bedge_facex_z(:,2,:,:,lb) = recvarx2e(:,1,:,:)
               endif
            end if

            if(jf.eq.3) then
               bedge_facey_x(:,:,1,:,lb) = recvary1e(:,:,2,:)
            elseif(jf.eq.4) then
               bedge_facey_x(:,:,2,:,lb) = recvary1e(:,:,1,:)
            end if

            if(jf.eq.3) then
               if((ndim.eq.3).or.(l2p5d.eq.1)) then
                  bedge_facey_z(:,:,1,:,lb) = recvary2e(:,:,2,:)
               endif
            elseif(jf.eq.4) then
               if((ndim.eq.3).or.(l2p5d.eq.1)) then
                  bedge_facey_z(:,:,2,:,lb) = recvary2e(:,:,1,:)
               endif
            end if

            if(jf.eq.5) then
               bedge_facez_x(:,:,:,1,lb) = recvarz1e(:,:,:,2)
            elseif(jf.eq.6) then
               bedge_facez_x(:,:,:,2,lb) = recvarz1e(:,:,:,1)
            endif

            if(jf.eq.5) then
               bedge_facez_y(:,:,:,1,lb) = recvarz2e(:,:,:,2)
            elseif(jf.eq.6) then
               bedge_facez_y(:,:,:,2,lb) = recvarz2e(:,:,:,1)
            endif

         endif                  ! end if (remote_pe

! make common variables on an edge consistent

         if(jf.eq.1) then

            bedge_facey_z(:,1+nguard0,1,klo:kup,lb) =
     .           bedge_facex_z(:,1,1+nguard0*k2d,klo:kup,lb)

            bedge_facey_z(:,1+nguard0,2,klo:kup,lb) =
     .           bedge_facex_z(:,1,k2d+nyb+nguard0*k2d,klo:kup,lb)

            if((ndim.eq.3).or.(l2p5d.eq.1)) then
               bedge_facez_y(:,1+nguard0,
     .              1+nguard0*k2d:nyb+nguard0*k2d,1,lb)
     .              = bedge_facex_y(:,1,
     .              1+nguard0*k2d:nyb+nguard0*k2d,klo,lb)

               if(ndim.eq.3)
     .              bedge_facez_y(:,1+nguard0,
     .              1+nguard0*k2d:nyb+nguard0*k2d,2,lb)
     .              = bedge_facex_y(:,1,
     .              1+nguard0*k2d:nyb+nguard0*k2d,kup1,lb)
            endif

         elseif(jf.eq.2) then

            bedge_facey_z(:,1+nxb+nguard0,1,klo:kup,lb) =
     .           bedge_facex_z(:,2,1+nguard0*k2d,klo:kup,lb)

            bedge_facey_z(:,1+nxb+nguard0,2,klo:kup,lb) =
     .           bedge_facex_z(:,2,k2d+nyb+nguard0*k2d,
     .                          klo:kup,lb)

            if((ndim.eq.3).or.(l2p5d.eq.1)) then
               bedge_facez_y(:,1+nxb+nguard0,
     .              1+nguard0*k2d:nyb+nguard0*k2d,
     .              1,lb)=
     .              bedge_facex_y(:,2,1+nguard0*k2d:nyb+nguard0*k2d,
     .              klo,lb)

               if(ndim.eq.3)
     .              bedge_facez_y(:,1+nxb+nguard0,
     .              1+nguard0*k2d:nyb+nguard0*k2d,
     .              2,lb)=
     .              bedge_facex_y(:,2,1+nguard0*k2d:nyb+nguard0*k2d,
     .              kup1,lb)

            endif

         elseif(jf.eq.3) then

            bedge_facex_z(:,1,1+nguard0*k2d,klo:kup,lb) =
     .           bedge_facey_z(:,1+nguard0,1,klo:kup,lb)

            bedge_facex_z(:,2,1+nguard0*k2d,klo:kup,lb) =
     .           bedge_facey_z(:,1+nxb+nguard0,1,klo:kup,lb)

            if((ndim.eq.3).or.(l2p5d.eq.1)) then
               bedge_facez_x(:,1+nguard0:nxb+nguard0,
     .              1+nguard0*k2d,1,lb)=
     .              bedge_facey_x(:,1+nguard0:nxb+nguard0,1,klo,lb)

               if(ndim.eq.3)
     .              bedge_facez_x(:,1+nguard0:nxb+nguard0,
     .              1+nguard0*k2d,2,lb)=
     .              bedge_facey_x(:,1+nguard0:nxb+nguard0,1,kup1,lb)
            endif

         elseif(jf.eq.4) then

            bedge_facex_z(:,1,k2d+nyb+nguard0*k2d,klo:kup,lb) =
     .           bedge_facey_z(:,1+nguard0,2,klo:kup,lb)

            bedge_facex_z(:,2,k2d+nyb+nguard0*k2d,klo:kup,lb) =
     .           bedge_facey_z(:,1+nxb+nguard0,2,klo:kup,lb)

            if((ndim.eq.3).or.(l2p5d.eq.1)) then
               bedge_facez_x(:,1+nguard0:nxb+nguard0,
     .              k2d+nyb+nguard0*k2d,
     .              1,lb)=
     .              bedge_facey_x(:,1+nguard0:nxb+nguard0,2,klo,lb)

               if(ndim.eq.3)
     .              bedge_facez_x(:,1+nguard0:nxb+nguard0,
     .              k2d+nyb+nguard0*k2d,
     .              2,lb)=
     .              bedge_facey_x(:,1+nguard0:nxb+nguard0,2,kup1,lb)
            endif

         elseif(jf.eq.5) then

            bedge_facey_x(:,1+nguard0:nxb+nguard0,1,klo,lb)=
     .           bedge_facez_x(:,1+nguard0:nxb+nguard0,
     .           1+nguard0*k2d,
     .           1,lb)

            bedge_facey_x(:,1+nguard0:nxb+nguard0,2,klo,lb)=
     .           bedge_facez_x(:,1+nguard0:nxb+nguard0,
     .           k2d+nyb+nguard0*k2d,1,lb)

            bedge_facex_y(:,1,1+nguard0*k2d:nyb+nguard0*k2d,
     .           klo,lb)=
     .           bedge_facez_y(:,1+nguard0,
     .           1+nguard0*k2d:nyb+nguard0*k2d,
     .           1,lb)

            bedge_facex_y(:,2,1+nguard0*k2d:nyb+nguard0*k2d,
     .           klo,lb)=
     .           bedge_facez_y(:,1+nxb+nguard0,
     .           1+nguard0*k2d:nyb+nguard0*k2d
     .           ,1,lb)

         elseif(jf.eq.6) then

            bedge_facey_x(:,1+nguard0:nxb+nguard0,1,kup1,lb)=
     .           bedge_facez_x(:,1+nguard0:nxb+nguard0,1+nguard0*k2d,
     .           2,lb)

            bedge_facey_x(:,1+nguard0:nxb+nguard0,2,kup1,lb)=
     .           bedge_facez_x(:,1+nguard0:nxb+nguard0,
     .           k2d+nyb+nguard0*k2d,
     .           2,lb)

            bedge_facex_y(:,1,1+nguard0*k2d:nyb+nguard0*k2d,
     .           kup1,lb)=
     .           bedge_facez_y(:,1+nguard0,
     .           1+nguard0*k2d:nyb+nguard0*k2d,2,lb)

            bedge_facex_y(:,2,1+nguard0*k2d:nyb+nguard0*k2d,
     .           kup1,lb)=
     .           bedge_facez_y(:,1+nxb+nguard0,
     .           1+nguard0*k2d:nyb+nguard0*k2d,
     .           2,lb)

        endif

        endif                     ! end if (cnodetype

        enddo                     ! end loop do jf =

      endif                       ! end of ldtcomplete if test

      endif ! end if (nodetype
      enddo ! end loop over blocks
      endif

!------------------------------------

      end if ! if (ndim >= 2)

      endif

      return
      end subroutine amr_edge_average_vdt
