!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_edge_diagonal_check
!! NAME
!!
!!   amr_edge_diagonal_check (mype)
!!
!! SYNOPSIS
!!
!!   call amr_edge_diagonal_check (mype)
!!
!!   call amr_edge_diagonal_check (integer)
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: mype
!!     The calling processors id.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   mpi_morton
!!   paramesh_interfaces
!!   paramesh_mpi_interfaces
!!
!! CALLS
!!
!!   amr_mpi_find_blk_in_buffer
!!   mpi_set_message_limits
!!   mpi_put_edge_buffer_1blk
!!
!! RETURNS
!!
!!   Nothing returned. 
!!
!! DESCRIPTION
!! 
!!  This routine checks to see if the diagonal block between two
!!  leaf-neighbors at the same refinement level as the current block,
!!  is refined. If it is then the edge-based variables along the edge
!!  shared with that diagonal block is given the edge values
!!  form the refined diagonal block, to insure conservation properties.
!! 
!! AUTHORS
!!
!!  Written by Peter MacNeice (October 1997).
!!
!!***

#include "paramesh_preprocessor.fh"

!#define DEBUG

      subroutine amr_edge_diagonal_check(mype)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

      use paramesh_interfaces, only : amr_mpi_find_blk_in_buffer
      use paramesh_mpi_interfaces, only : mpi_set_message_limits,
     .                                    mpi_put_edge_buffer_1blk

      implicit none 

#include "amr_shmem.fh"
      include 'mpif.h'

      integer, intent(in)  ::  mype

!------------------------------------
! local variables

      integer :: nguard0
      integer :: klo,kup
      integer :: jlo,jup
      integer :: ilo,iup

      integer :: remote_pe,remote_block
      integer :: remote_pe2,remote_block2
      integer :: mark_edge(12,maxblocks)
      integer :: i, ie, iblk, lb, k, j, ierrcode
      integer :: ia, ib, ja, jb, ka, kb
      integer :: ia0, ib0, ja0, jb0, ka0, kb0
      integer :: dtype, vtype
      integer :: index, index0, n

      logical :: lfound

      real,allocatable :: receive(:,:)

!------------------------------------

      nguard0 = nguard*npgs
      klo=1+nguard0*k3d
      kup=klo*k3d+nzb
      jlo=1+nguard0*k2d
      jup=jlo*k2d+nyb
      ilo=1+nguard0
      iup=ilo+nxb

      allocate(receive(nedges,maxdim+2*nguard0))

#ifdef DEBUG
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, ierrcode)
       write(*,*) 'starting amr_edge_diagonal_check : pe ',mype
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, ierrcode)
#endif /* DEBUG */

       if (ndim >= 2) then


! Initialize array marking edges for diagonal patching.
      mark_edge(:,:) = 0

      do i = 1,no_of_diagonal_edges
        ie   = edge_mark(6,1,i)
        iblk = edge_mark(6,2,i)
        mark_edge(ie,iblk) = i
      enddo


! Loop over the blocks on this processor.
      if(lnblocks.gt.0) then
      do lb=1,lnblocks

! Is this a leaf block which has finished its current timestep?
      if((nodetype(lb).eq.1.and..not.var_dt) .or.
     .   (nodetype(lb).eq.1.and.ldtcomplete(lb))) then

#ifdef DEBUG
         write(*,*) 'amr_edge_diagonal_check : ',
     .              'checking edges on blk ',lb
         write(*,*) 'nbedges ',nbedges
#endif /* DEBUG */


! Any edges on this block which are still marked need a diagonal patch.
! Note that in the shmem_gets below, we can always assume that a
! neighbor block exists, since the edge would not have been marked
! earlier if that was not so.

! Loop over the edges on this block.
       do ie=1,nbedges


       if(mark_edge(ie,lb).ge.1) then

        lfound = .false.
        remote_block  = edge_mark(6,3,mark_edge(ie,lb))
        remote_pe     = edge_mark(6,4,mark_edge(ie,lb))
        remote_block2 = edge_mark(6,3,mark_edge(ie,lb))
        remote_pe2    = edge_mark(6,4,mark_edge(ie,lb))

#ifdef DEBUG
         write(*,*) 'amr_edge_diagonal_check : ',
     .              'data source edge for ',lb, 
     .              ' is on blk ',remote_block,remote_pe
         write(*,*) 'nbedges ',nbedges
#endif /* DEBUG */

! (remote_block,remote_pe) may be a local block or a remote block.


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


! The edge data on the neighboring faces can be assumed to have been averaged
! correctly from the refined diagonal blocks.

      if (remote_pe == mype .and. remote_block <= lnblocks) then

! Now copy over the edge data from one of the neighbors.
       if(ie.eq.1) then                    ! edge: x low edge, y low edge
         do k=klo,kup-k3d
           call shmem_real_get(
     .         receive(1:nedgevar,k),
     .         bedge_facex_z(1:nedgevar,2,jup,k,remote_block),
     .         nedgevar,remote_pe)
           bedge_facex_z(:,1,1+nguard0*k2d,k,lb)= receive(:,k)
         enddo
         bedge_facey_z(:,1+nguard0,1,klo:kup-k3d,lb)=
     .                bedge_facex_z(:,1,1+nguard0*k2d,klo:kup-k3d,lb)

       elseif(ie.eq.2) then               ! edge: x low edge, y high edge
         do k=klo,kup-k3d
           call shmem_real_get(
     .         receive(1:nedgevar,k),
     .         bedge_facex_z(1:nedgevar,2,jlo,k,remote_block),
     .         nedgevar,remote_pe)
           bedge_facex_z(:,1,k2d+nguard0*k2d+nyb,k,lb)= 
     .         receive(:,k)
         enddo
         bedge_facey_z(:,1+nguard0,2,klo:kup-k3d,lb)=
     .     bedge_facex_z(:,1,k2d+nguard0*k2d+nyb,klo:kup-k3d,lb)

       elseif(ie.eq.3) then               ! edge: x high edge, y low edge
         do k=klo,kup-k3d
           call shmem_real_get(
     .         receive(1:nedgevar,k),
     .         bedge_facex_z(1:nedgevar,1,jup,k,remote_block),
     .         nedgevar,remote_pe)
           bedge_facex_z(:,2,1+nguard0*k2d,k,lb)= receive(:,k)
         enddo
         bedge_facey_z(:,1+nguard0+nxb,1,klo:kup-k3d,lb)=
     .     bedge_facex_z(:,2,1+nguard0*k2d,klo:kup-k3d,lb)


       elseif(ie.eq.4) then               ! edge: x high edge, y high edge
         do k=klo,kup-k3d
           call shmem_real_get(
     .         receive(1:nedgevar,k),
     .         bedge_facex_z(1:nedgevar,1,jlo,k,remote_block),
     .         nedgevar,remote_pe)
           bedge_facex_z(:,2,k2d+nguard0*k2d+nyb,k,lb)= 
     .         receive(:,k)
         enddo
         bedge_facey_z(:,1+nguard0+nxb,2,klo:kup-k3d,lb)=
     .     bedge_facex_z(:,2,k2d+nguard0*k2d+nyb,klo:kup-k3d,lb)


       elseif(ie.eq.5) then                ! edge: y low edge, z low edge
         do i=ilo,iup-1
           call shmem_real_get(
     .         receive(1:nedgevar,i),
     .         bedge_facey_x(1:nedgevar,i,2,kup,remote_block),
     .         nedgevar,remote_pe)
           bedge_facey_x(:,i,1,klo,lb)= receive(:,i)
         enddo
         bedge_facez_x(:,ilo:iup-1,1+nguard0*k3d,1,lb)=
     .                bedge_facey_x(:,ilo:iup-1,1,klo,lb)


       elseif(ie.eq.6) then                ! edge: y high edge, z low edge
         do i=ilo,iup-1
           call shmem_real_get(
     .         receive(1:nedgevar,i),
     .         bedge_facey_x(1:nedgevar,i,1,kup,remote_block),
     .         nedgevar,remote_pe)
           bedge_facey_x(:,i,2,klo,lb)= receive(:,i)
         enddo
         bedge_facez_x(:,ilo:iup-1,k2d+nguard0*k2d+nyb,1,lb)=
     .                bedge_facey_x(:,ilo:iup-1,2,klo,lb)


       elseif(ie.eq.7) then                ! edge: y low edge, z high edge
         do i=ilo,iup-1
           call shmem_real_get(
     .         receive(1:nedgevar,i),
     .         bedge_facey_x(1:nedgevar,i,2,klo,remote_block),
     .         nedgevar,remote_pe)
           bedge_facey_x(:,i,1,kup,lb)= receive(:,i)
         enddo
         bedge_facez_x(:,ilo:iup-1,1+nguard0*k2d,2,lb)=
     .                bedge_facey_x(:,ilo:iup-1,1,kup,lb)


       elseif(ie.eq.8) then                ! edge: y high edge, z high edge
         do i=ilo,iup-1
           call shmem_real_get(
     .         receive(1:nedgevar,i),
     .         bedge_facey_x(1:nedgevar,i,1,klo,remote_block),
     .         nedgevar,remote_pe)
           bedge_facey_x(:,i,2,kup,lb)= receive(:,i)
         enddo
         bedge_facez_x(:,ilo:iup-1,k2d+nguard0*k2d+nyb,2,lb)=
     .                bedge_facey_x(:,ilo:iup-1,2,kup,lb)

       elseif(ie.eq.9) then                ! edge: x low edge, z low edge
         do j=jlo,jup-k2d
           call shmem_real_get(
     .         receive(1:nedgevar,j),
     .         bedge_facex_y(1:nedgevar,2,j,kup,remote_block),
     .         nedgevar,remote_pe)
           bedge_facex_y(:,1,j,klo,lb)= receive(:,j)
         enddo
         bedge_facez_y(:,1+nguard0,jlo:jup-k2d,1,lb)=
     .                bedge_facex_y(:,1,jlo:jup-k2d,klo,lb)


       elseif(ie.eq.10) then                ! edge: x low edge, z high edge
         do j=jlo,jup-k2d
           call shmem_real_get(
     .         receive(1:nedgevar,j),
     .         bedge_facex_y(1:nedgevar,2,j,klo,remote_block),
     .         nedgevar,remote_pe)
           bedge_facex_y(:,1,j,kup,lb)= receive(:,j)
         enddo
         bedge_facez_y(:,1+nguard0,jlo:jup-k2d,2,lb)=
     .                bedge_facex_y(:,1,jlo:jup-k2d,kup,lb)


       elseif(ie.eq.11) then                ! edge: x high edge, z low edge
         do j=jlo,jup-k2d
           call shmem_real_get(
     .         receive(1:nedgevar,j),
     .         bedge_facex_y(1:nedgevar,1,j,kup,remote_block),
     .         nedgevar,remote_pe)
           bedge_facex_y(:,2,j,klo,lb)= receive(:,j)
         enddo
         bedge_facez_y(:,1+nguard0+nxb,jlo:jup-k2d,1,lb)=
     .                bedge_facex_y(:,2,jlo:jup-k2d,klo,lb)


       elseif(ie.eq.12) then                ! edge: x high edge, z high edge
         do j=jlo,jup-k2d
           call shmem_real_get(
     .         receive(1:nedgevar,j),
     .         bedge_facex_y(1:nedgevar,1,j,klo,remote_block),
     .         nedgevar,remote_pe)
           bedge_facex_y(:,2,j,kup,lb)= receive(:,j)
         enddo
         bedge_facez_y(:,1+nguard0+nxb,jlo:jup-k2d,2,lb)=
     .                bedge_facex_y(:,2,jlo:jup-k2d,kup,lb)
       endif


      else                      ! if (remote_pe

!#ifdef NEVER

         call mpi_put_edge_buffer_1blk(lb,remote_block,remote_pe)

         if(ie.eq.1) then       ! edge: x low edge, y low edge
            do k=klo,kup-k3d
               bedge_facex_z(:,1,jlo,k,lb)= 
     .              recvarx2e(:,2,jup,k)
            enddo
            bedge_facey_z(:,ilo,1,klo:kup-k3d,lb)=
     .           bedge_facex_z(:,1,jlo,klo:kup-k3d,lb)
            
         elseif(ie.eq.2) then   ! edge: x low edge, y high edge
            do k=klo,kup-k3d
               bedge_facex_z(:,1,jup,k,lb)= 
     .              recvarx2e(:,2,jlo,k)
            enddo
            bedge_facey_z(:,ilo,2,klo:kup-k3d,lb)=
     .           bedge_facex_z(:,1,jup,klo:kup-k3d,lb)
            
         elseif(ie.eq.3) then   ! edge: x high edge, y low edge
            do k=klo,kup-k3d
               bedge_facex_z(:,2,jlo,k,lb)= 
     .              recvarx2e(:,1,jup,k)
            enddo
            bedge_facey_z(:,iup,1,klo:kup-k3d,lb)=
     .           bedge_facex_z(:,2,jlo,klo:kup-k3d,lb)
            
         elseif(ie.eq.4) then   ! edge: x high edge, y high edge
            do k=klo,kup-k3d
               bedge_facex_z(:,2,jup,k,lb)= 
     .              recvarx2e(:,1,jlo,k)
            enddo
            bedge_facey_z(:,iup,2,klo:kup-k3d,lb)=
     .           bedge_facex_z(:,2,jup,klo:kup-k3d,lb)
         elseif(ie.eq.5) then   ! edge: y low edge, z low edge
            do i=ilo,iup-1
               bedge_facey_x(:,i,1,klo,lb)= recvary1e(:,i,2,kup)
            enddo
            bedge_facez_x(:,ilo:iup-1,jlo,1,lb)=
     .           bedge_facey_x(:,ilo:iup-1,1,klo,lb)
            
         elseif(ie.eq.6) then   ! edge: y high edge, z low edge
            do i=ilo,iup-1
               bedge_facey_x(:,i,2,klo,lb)= recvary1e(:,i,1,kup)
            enddo
            bedge_facez_x(:,ilo:iup-1,jup,1,lb)=
     .           bedge_facey_x(:,ilo:iup-1,2,klo,lb)
            
         elseif(ie.eq.7) then   ! edge: y low edge, z high edge
            do i=ilo,iup-1
               bedge_facey_x(:,i,1,kup,lb)= recvary1e(:,i,2,klo)
            enddo
            bedge_facez_x(:,ilo:iup-1,jlo,2,lb)=
     .           bedge_facey_x(:,ilo:iup-1,1,kup,lb)
            
         elseif(ie.eq.8) then   ! edge: y high edge, z high edge
            do i=ilo,iup-1
               bedge_facey_x(:,i,2,kup,lb)= recvary1e(:,i,1,klo)
            enddo
            bedge_facez_x(:,ilo:iup-1,jup,2,lb)=
     .           bedge_facey_x(:,ilo:iup-1,2,kup,lb)
            
         elseif(ie.eq.9) then   ! edge: x low edge, z low edge
            do j=jlo,jup-k2d
               bedge_facex_y(:,1,j,klo,lb)= recvarx1e(:,2,j,kup)
            enddo
            bedge_facez_y(:,1+nguard0,jlo:jup-k2d,1,lb)=
     .           bedge_facex_y(:,1,jlo:jup-k2d,klo,lb)
            
         elseif(ie.eq.10) then  ! edge: x low edge, z high edge
            do j=jlo,jup-k2d
               bedge_facex_y(:,1,j,kup,lb)= recvarx1e(:,2,j,klo)
            enddo
            bedge_facez_y(:,1+nguard0,jlo:jup-k2d,2,lb)=
     .           bedge_facex_y(:,1,jlo:jup-k2d,kup,lb)
            
         elseif(ie.eq.11) then  ! edge: x high edge, z low edge
            do j=jlo,jup-k2d
               bedge_facex_y(:,2,j,klo,lb)= recvarx1e(:,1,j,kup)
            enddo
            bedge_facez_y(:,1+nguard0+nxb,jlo:jup-k2d,1,lb)=
     .           bedge_facex_y(:,2,jlo:jup-k2d,klo,lb)
            
         elseif(ie.eq.12) then  ! edge: x high edge, z high edge
            do j=jlo,jup-k2d
               bedge_facex_y(:,2,j,kup,lb)= recvarx1e(:,1,j,klo)
            enddo
            bedge_facez_y(:,1+nguard0+nxb,jlo:jup-k2d,2,lb)=
     .           bedge_facex_y(:,2,jlo:jup-k2d,kup,lb)
            
         end if
!#endif

#ifdef NEVER

         lfound = .false.
         call amr_mpi_find_blk_in_buffer(mype,remote_block,
     .        remote_pe,1,dtype,index0,lfound)
         vtype = 8
         call mpi_set_message_limits(dtype,
     .        ia0,ib0,ja0,jb0,ka0,kb0,vtype)

         index = index0 + 1

         if(dtype.eq.13.or.dtype.eq.15.or.dtype.eq.14) then

            ia = ia0
            ib = ib0
            ja = ja0
            jb = jb0
            ka = ka0
            kb = kb0
            
            if(dtype.eq.13) then
               ia = 1
               ib = 1
            elseif(dtype.eq.15) then
               ia = 2
               ib = 2
            elseif(dtype.eq.14) then
               ia = 1
               ib = 2
            endif
               
                  
            do k = ka,kb
               do j = ja,jb
                  do i = ia,ib
                     do n=1,nedges
                        recvarx1e(n,i,j,k) =
     .                          temprecv_buf(index)
                        index  = index + 1
                     enddo
                  enddo
               enddo
            enddo
            
            if(ie.eq.9) then    ! edge: x low edge, z low edge
               do j=jlo,jup-k2d
                  bedge_facex_y(:,1,j,klo,lb)= recvarx1e(:,2,j,kup)
               enddo
               bedge_facez_y(:,1+nguard0,jlo:jup-k2d,1,lb)=
     .           bedge_facex_y(:,1,jlo:jup-k2d,klo,lb)

            elseif(ie.eq.10) then ! edge: x low edge, z high edge
               do j=jlo,jup-k2d
                  bedge_facex_y(:,1,j,kup,lb)= recvarx1e(:,2,j,klo)
               enddo
               bedge_facez_y(:,1+nguard0,jlo:jup-k2d,2,lb)=
     .              bedge_facex_y(:,1,jlo:jup-k2d,kup,lb)

            elseif(ie.eq.11) then ! edge: x high edge, z low edge
               do j=jlo,jup-k2d
                  bedge_facex_y(:,2,j,klo,lb)= recvarx1e(:,1,j,kup)
               enddo
               bedge_facez_y(:,1+nguard0+nxb,jlo:jup-k2d,1,lb)=
     .              bedge_facex_y(:,2,jlo:jup-k2d,klo,lb)
               
            elseif(ie.eq.12) then ! edge: x high edge, z high edge
               do j=jlo,jup-k2d
                  bedge_facex_y(:,2,j,kup,lb)= recvarx1e(:,1,j,klo)
               enddo
               bedge_facez_y(:,1+nguard0+nxb,jlo:jup-k2d,2,lb)=
     .              bedge_facex_y(:,2,jlo:jup-k2d,kup,lb)

            end if

            if(ndim.eq.3.or.l2p5d.eq.1) then
               do k = ka , kb
                  do j = ja , jb
                     do i = ia , ib
                        do n=1,nedges
                           recvarx1e(n,i,j,k) =
     .                          temprecv_buf(index)
                           index  = index + 1
                        enddo
                     enddo
                  enddo
               enddo
               
               if(ie.eq.1) then ! edge: x low edge, y low edge
                  do k=klo,kup-k3d
                     bedge_facex_z(:,1,jlo,k,lb)= 
     .                    recvarx1e(:,2,jup,k)
                  enddo
                  bedge_facey_z(:,ilo,1,klo:kup-k3d,lb)=
     .              bedge_facex_z(:,1,jlo,klo:kup-k3d,lb)

               elseif(ie.eq.2) then ! edge: x low edge, y high edge
                  do k=klo,kup-k3d
                     bedge_facex_z(:,1,jup,k,lb)= 
     .                    recvarx1e(:,2,jlo,k)
                  enddo
                  bedge_facey_z(:,ilo,2,klo:kup-k3d,lb)=
     .             bedge_facex_z(:,1,jup,klo:kup-k3d,lb)
                  
               elseif(ie.eq.3) then ! edge: x high edge, y low edge
                  do k=klo,kup-k3d
                     bedge_facex_z(:,2,jlo,k,lb)= 
     .                    recvarx1e(:,1,jup,k)
                  enddo
                  bedge_facey_z(:,iup,1,klo:kup-k3d,lb)=
     .              bedge_facex_z(:,2,jlo,klo:kup-k3d,lb)

               elseif(ie.eq.4) then ! edge: x high edge, y high edge
                  do k=klo,kup-k3d
                     bedge_facex_z(:,2,jup,k,lb)= 
     .                    recvarx1e(:,1,jlo,k)
                  enddo
                  bedge_facey_z(:,iup,2,klo:kup-k3d,lb)=
     .             bedge_facex_z(:,2,jup,klo:kup-k3d,lb)
               end if

            endif               ! end if (ndim

         end if

         if(ndim.ge.2) then
         if(dtype.eq.11.or.dtype.eq.17.or.dtype.eq.14) then

            ia = ia0
            ib = ib0
            ja = ja0
            jb = jb0
            ka = ka0
            kb = kb0

            if(dtype.eq.11) then
               ja = 1
               jb = 1
            elseif(dtype.eq.17) then
               ja = 2
               jb = 2
            elseif(dtype.eq.14) then
               ja = 1
               jb = 2
            endif

            do k = ka,kb
               do j = ja,jb
                  do i = ia,ib
                     do n=1,nedges
                        recvary1e(n,i,j,k) = temprecv_buf(index)
                        index  = index + 1
                     enddo
                  enddo
               enddo
            enddo
               
            if(ie.eq.5) then    ! edge: y low edge, z low edge
               do i=ilo,iup-1
                  bedge_facey_x(:,i,1,klo,lb)= recvary1e(:,i,2,kup)
               enddo
               bedge_facez_x(:,ilo:iup-1,jlo,1,lb)=
     .              bedge_facey_x(:,ilo:iup-1,1,klo,lb)
               
            elseif(ie.eq.6) then ! edge: y high edge, z low edge
               do i=ilo,iup-1
                  bedge_facey_x(:,i,2,klo,lb)= recvary1e(:,i,1,kup)
               enddo
               bedge_facez_x(:,ilo:iup-1,jup,1,lb)=
     .              bedge_facey_x(:,ilo:iup-1,2,klo,lb)
               
            elseif(ie.eq.7) then ! edge: y low edge, z high edge
               do i=ilo,iup-1
                  bedge_facey_x(:,i,1,kup,lb)= recvary1e(:,i,2,klo)
               enddo
               bedge_facez_x(:,ilo:iup-1,jlo,2,lb)=
     .              bedge_facey_x(:,ilo:iup-1,1,kup,lb)
               
            elseif(ie.eq.8) then ! edge: y high edge, z high edge
               do i=ilo,iup-1
                  bedge_facey_x(:,i,2,kup,lb)= recvary1e(:,i,1,klo)
               enddo
               bedge_facez_x(:,ilo:iup-1,jup,2,lb)=
     .              bedge_facey_x(:,ilo:iup-1,2,kup,lb)
               
            end if
            
            if(ndim.eq.3.or.l2p5d.eq.1) then
               do k = ka , kb
                  do j = ja , jb
                     do i = ia , ib
                        do n=1,nedges
                           index  = index + 1
                        enddo
                     enddo
                  enddo
               enddo
            endif
            
         end if
         end if

#endif

      end if                    ! end if (remote_pe

      endif


      enddo                     ! loop over edges

      endif

      enddo
      endif

      end if ! if (ndim >= 2)

      deallocate(receive)

      return
      end subroutine amr_edge_diagonal_check
