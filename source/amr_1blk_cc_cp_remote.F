!#define DEBUG

#define GCELL_ON
!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#define TESTING_RESTRICT
#include "paramesh_preprocessor.fh"


      subroutine amr_1blk_cc_cp_remote(mype,remote_pe,remote_block,
     .    idest,iopt,id,jd,kd,is,js,ks,ilays,jlays,klays,
     .    nblk_ind,ipolar)



!------------------------------------------------------------------------
!
! This routine copies guard cell information to face iface in layer
! idest of the working block, from the appropriate face of the neighboring 
! block, assuming that the neighboring block is on a different processor.
!
! Arguments:
!      mype             local processor
!      remote_pe        remote processor
!      remote_block     local block id of the block to be copied from
!                        the remote processor
!      idest            selects the storage space in data_1blk.fh which is to
!                        be used in this call. If the leaf node is having its
!                        guardcells filled then set this to 1, if its parent
!                        is being filled set it to 2.
!      iopt             a switch to control which data source is to be used
!                        iopt=1 will use 'unk'
!                        iopt>=2 will use 'work'
!      id               lower limit of index range of points in x direction
!                        on destination block
!      jd               lower limit of index range of points in y direction
!                        on destination block
!      kd               lower limit of index range of points in z direction
!                        on destination block
!      is               lower limit of index range of points in x direction
!                        on source block
!      js               lower limit of index range of points in y direction
!                        on source block
!      ks               lower limit of index range of points in z direction
!                        on source block
!      ilay             no. of mesh points in x direction to be copied
!      jlay             no. of mesh points in y direction to be copied
!      klay             no. of mesh points in z direction to be copied
!      nblk_ind        index, running from 1-27 denoting location of neighbor block
!
!
!
! Written :     Peter MacNeice          July 1998
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
#ifdef MPI_USED
      use mpi_morton

      use paramesh_interfaces, only : amr_mpi_find_blk_in_buffer
      use paramesh_mpi_interfaces, only : mpi_set_message_limits
      use timings

#endif /* MPI_USED */

      implicit none

#include "amr_shmem.fh"



!-------------------------
      integer, intent(in) :: mype,remote_pe,remote_block
      integer, intent(in) :: idest,iopt,id,jd,kd,is,js,ks
      integer, intent(in) :: ilays,jlays,klays,nblk_ind
      integer, intent(in) :: ipolar(:)


! local variables and arrays

#ifndef MPI_USED
      real recv(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .       kl_bnd:ku_bnd)
      save recv
#endif

      integer :: il,jl,kl,iopt0
      integer :: ill,jll,kll
      integer :: ia, ib, ja, jb, ka, kb, jstride, js0, jsl
      integer :: js1,js2
      integer :: i, j, k, ii, jj, kk, ivar, ivar_next
      integer :: indx

#ifdef MPI_USED
      include 'mpif.h'
      double precision :: time1
      integer :: ierr,dtype
      integer :: vtype
      logical :: lfound
#endif /* MPI_USED */

!-------------------------

#ifdef DEBUG
       if(l_f_to_c) then
       write(*,*) 'amr_1blk_cc_cp_remote args : ',
     .    mype,remote_pe,remote_block,
     .    idest,iopt,id,jd,kd,is,js,ks,ilays,jlays,klays
       endif
#endif /* DEBUG */

!
! Adjust index ranges
      il = ilays-1
      jl = (jlays-1)*k2d
      kl = (klays-1)*k3d

#ifdef MPI_USED
      if(remote_block.le.lnblocks.and.remote_pe.eq.mype) then

         if (timing_mpix) then
            time1 = mpi_wtime()
         endif

#endif /* MPI_USED */

         jstride = 1


!-------------------------
      if(iopt.eq.1) then
!-------------------------

#ifdef MPI_USED
         ia = is
         ib = is+il
         ja = js
         jb = js+jl
         ka = ks
         kb = ks+kl
#endif

!#ifdef SPHERICAL
      if(spherical_pm) then
      if(lsingular_line) then
      if(ipolar(1).eq.-1.and.jd.le.nguard) then
        jstride = -1
        if (no_permanent_guardcells) then
          ja = nguard-jd+1
          jb = ja - jl
        else
          ja = 2*nguard     
          jb = ja -nguard     
        endif
      elseif(ipolar(2).eq.+1.and.jd.gt.nyb+nguard) then
        jstride = -1
        if (no_permanent_guardcells) then
          ja = (nyb+1)-(jd-(nyb+nguard))
          jb = ja - jl
        else
          ja = (nyb+nguard)-(jd-(nyb+nguard+1))
          jb = (nyb+nguard)-((jd+jl)-(nyb+nguard+1))
        endif
      endif
      endif
      endif
!#endif /* SPHERICAL */

! Copy complete remote block into a buffer block called recv.
        if (no_permanent_guardcells) then      

#ifndef MPI_USED
        call shmem_real_get( recv(1,1,1,1),
     .                       gt_unk(1,1,1,1,remote_block),
     .                       len_block,remote_pe) 
#else
        do ivar=1,nvar
          if(int_gcell_on_cc(ivar)) then

        if(.not.l_f_to_c) then
          unk1(ivar,id:id+il,jd:jd+jl,kd:kd+kl,idest) =
     .       gt_unk(ivar,ia:ib,ja:jb:jstride,ka:kb,remote_block)
        else
          unk1_fl(ivar,id:id+il,jd:jd+jl,kd:kd+kl) =
     .       gt_unk(ivar,ia:ib,ja:jb:jstride,ka:kb,remote_block)
        endif

          endif
        enddo
#endif

        else ! no_permanent_guardcells

#ifndef MPI_USED
        call shmem_real_get( recv(1,1,1,1),
     .                       unk(1,1,1,1,remote_block),
     .                       len_block,remote_pe) 
#else
        do ivar=1,nvar
          if(int_gcell_on_cc(ivar)) then

        if(.not.l_f_to_c) then
          unk1(ivar,id:id+il,jd:jd+jl,kd:kd+kl,idest) =
     .       unk(ivar,ia:ib,ja:jb:jstride,ka:kb,remote_block)
        else
          unk1_fl(ivar,id:id+il,jd:jd+jl,kd:kd+kl) =
     .       unk(ivar,ia:ib,ja:jb:jstride,ka:kb,remote_block)
        endif

          endif
        enddo
#endif
        endif ! no_permanent_guardcells

#ifndef MPI_USED
        do ivar=1,nvar
          if(int_gcell_on_cc(ivar)) then

        if(.not.l_f_to_c) then
          unk1(ivar,id:id+il,jd:jd+jl,kd:kd+kl,idest) = 
     .             recv(ivar,is:is+il,js:js+jl,ks:ks+kl)
        else
          unk1_fl(ivar,id:id+il,jd:jd+jl,kd:kd+kl) = 
     .             recv(ivar,is:is+il,js:js+jl,ks:ks+kl)
        endif

          endif
        enddo
#endif

!-------------------------
      elseif(iopt.ge.2) then
!-------------------------
      jstride = 1
      js0 = js
      jsl = js0 + jl
!#ifdef SPHERICAL
      if(spherical_pm) then
      if(lsingular_line) then
      if(ipolar(1).eq.-1.and.jd.le.nguard_work) then
        jstride = -1
        if (no_permanent_guardcells) then
          js0 = nguard_work-jd+1
          jsl = js0 - jl
        else
          js0 = 2*nguard_work
          jsl = js0 -nguard_work
        endif
      elseif(ipolar(2).eq.+1.and.jd.gt.nyb+nguard) then
        jstride = -1
        if (no_permanent_guardcells) then
          js0 = (nyb+1)-(jd-(nyb+nguard_work))
          jsl = js0 - jl
        else
          js0 = (nyb+nguard_work)-(jd-(nyb+nguard_work+1))
          jsl = (nyb+nguard_work)-((jd+jl)-(nyb+nguard_work+1))
        endif
      endif
      endif
      endif
!#endif /* SPHERICAL */

        iopt0 = iopt-1
! Copy complete remote block into a buffer block called recvw.
#ifndef MPI_USED
        call shmem_real_get( recvw(1,1,1),
     .                       work(1,1,1,remote_block,iopt0),
     .                       len_wblock,remote_pe) 
        if(.not.l_f_to_c) then
          work1(id:id+il,jd:jd+jl,kd:kd+kl,idest) = 
     .               recvw(is:is+il,js0:js0+jsl:jstride,ks:ks+kl)
        else
          work1_fl(id:id+il,jd:jd+jl,kd:kd+kl) = 
     .               recvw(is:is+il,js0:js0+jsl:jstride,ks:ks+kl)
        endif
#else
        if(.not.l_f_to_c) then
          work1(id:id+il,jd:jd+jl,kd:kd+kl,idest) =
     .     work(is:is+il,js0:jsl:jstride,ks:ks+kl,remote_block,iopt0)
        else
! probably doesnt work for spherical polar blocks
          work1_fl(id:id+il,jd:jd+jl,kd:kd+kl) =
     .     work(is:is+il,js0:js0+jl,ks:ks+kl,remote_block,iopt0)
        endif
#endif

!-------------------------
      endif
!-------------------------

      if (timing_mpix) then
         timer_amr_1blk_cc_cp_remote(1) = 
     .     timer_amr_1blk_cc_cp_remote(1) + mpi_wtime() - time1
      else
         timer_amr_1blk_cc_cp_remote(1) = -999.
      endif

#ifdef MPI_USED
      else                        ! parallel section

         if (timing_mpix) then
            time1 = mpi_wtime()
         endif

        call amr_mpi_find_blk_in_buffer(mype,remote_block,
     .                        remote_pe,idest,dtype,indx,lfound)

        if (timing_mpix) then
         timer_amr_1blk_cc_cp_remote(2) = 
     .     timer_amr_1blk_cc_cp_remote(2) + mpi_wtime() - time1
         time1 = mpi_wtime()
        else
         timer_amr_1blk_cc_cp_remote(2) = -999.
        endif

! If this routine is executing a copy to fill guardcells of a
! leaf blocks^s parent, and the remote block is not found, then
! it is assumed that it is not in the list of buffered remote blocks
! because it is not really needed. Therefore in this case we
! return without copying anything.
        if(idest.eq.2.and.(.not.lfound)) then
           return
        end if



        if(iopt.eq.1) then
          vtype = 1
          if (l_f_to_c) then
             if (ilays == 2*nguard) then
                ill = nguard
             else
                ill = (ilays-1)/2
             end if
             if (jlays == 2*nguard) then
                jll = nguard
             else
                jll = (jlays-1)/2
             end if
             if (klays == 2*nguard) then
                kll = nguard
             else
                kll = (klays-1)/2
             end if
          else
             ill = ilays
             jll = jlays
             kll = klays
          end if

          call mpi_set_message_limits(
     .                 dtype,ia,ib,ja,jb,ka,kb,vtype,
     .                 ill,jll,kll)

          kk = kd
          do k = ka,kb
          jj = jd
          jstride = 1


      js2 = js
      js1 = js+jl


!#ifdef SPHERICAL
      if(spherical_pm) then
      if(lsingular_line) then
      if(ipolar(1).eq.-1.and.jd.le.nguard) then
        jstride = -1
        if (no_permanent_guardcells) then
          js2 = jd + jl + (nguard - 2*(jd+jl)) +1
          js1 = jd      + (nguard - 2* jd    ) +1
          jj  = jd + jl
        else
          js2 = jd + jl + 2*(nguard - (jd+jl)) +1
          js1 = jd      + 2*(nguard -  jd    ) +1
          jj  = jd + jl
        endif
      elseif(ipolar(2).eq.+1.and.jd.gt.nyb+nguard) then
        jstride = -1
        if (no_permanent_guardcells) then
          js1 = nyb - ( jd    -(nyb+nguard+1))
          js2 = nyb - ((jd+jl)-(nyb+nguard+1))
          jj  = jd + jl
        else
          js1 = (nyb+nguard)-( jd    -(nyb+nguard+1))
          js2 = (nyb+nguard)-((jd+jl)-(nyb+nguard+1))
          jj  = jd + jl
        endif
      endif
      endif
      endif
!#endif /* SPHERICAL */
          do j = ja,jb
          ii = id
          do i = ia,ib
            if (k >= ks .and. k <= ks + kl) then
!            if (j >= js .and. j <= js + jl) then
            if (j >= js2.and. j <= js1) then
            if (i >= is .and. i <= is + il) then

        do ivar=1,ngcell_on_cc
              ivar_next = gcell_on_cc_pointer(ivar)
              if(.not.l_f_to_c) then
                unk1(ivar_next,ii,jj,kk,idest) =
     .                 temprecv_buf(indx+ivar)
              else
                unk1_fl(ivar_next,ii,jj,kk) =
     .                 temprecv_buf(indx+ivar)
              endif
        enddo
            end if
            end if
            end if
            if (i >= is .and. i <= is + il) ii = ii + 1
#ifdef GCELL_ON
            indx = indx+ngcell_on_cc
#else
            indx = indx+nvar
#endif /* GCELL_ON */
          enddo
!          if (j >= js .and. j <= js + jl) jj = jj 
          if (j >= js2 .and. j <= js1) jj = jj + jstride
          enddo
          if (k >= ks .and. k <= ks + kl) kk = kk + 1
          enddo

        elseif(iopt.gt.1) then
          vtype = 0
          if (l_f_to_c) then
             if (ilays == 2*nguard_work) then
                ill = nguard_work
             else
                ill = (ilays-1)/2
             end if
             if (jlays == 2*nguard_work) then
                jll = nguard_work
             else
                jll = (jlays-1)/2
             end if
             if (klays == 2*nguard_work) then
                kll = nguard_work
             else
                kll = (klays-1)/2
             end if
          else
             ill = ilays
             jll = jlays
             kll = klays
          end if
          call mpi_set_message_limits(
     .                 dtype,ia,ib,ja,jb,ka,kb,vtype,
     .                 ill,jll,kll)
          kk = kd
          do k = ka,kb
          jj = jd
          jstride = 1


      js2 = js
      js1 = js+jl

!#ifdef SPHERICAL
      if(spherical_pm) then
      if(lsingular_line) then
      if(ipolar(1).eq.-1.and.jd.le.nguard_work) then
        jstride = -1
        if (no_permanent_guardcells) then
          js2 = jd + jl + (nguard_work - 2*(jd+jl)) +1
          js1 = jd      + (nguard_work - 2* jd    ) +1
          jj  = jd + jl
        else
          js2 = jd + jl + 2*(nguard_work - (jd+jl)) +1
          js1 = jd      + 2*(nguard_work -  jd    ) +1
          jj  = jd + jl
        endif
      elseif(ipolar(2).eq.+1.and.jd.gt.nyb+nguard_work) then
        jstride = -1
        if (no_permanent_guardcells) then
          js1 = nyb - ( jd    -(nyb+nguard_work+1))
          js2 = nyb - ((jd+jl)-(nyb+nguard_work+1))
          jj  = jd + jl
        else
          js1 = (nyb+nguard_work)-( jd    -(nyb+nguard_work+1))
          js2 = (nyb+nguard_work)-((jd+jl)-(nyb+nguard_work+1))
          jj  = jd + jl
        endif
      endif
      endif
      endif
!#endif /* SPHERICAL */

#ifdef NOTNOW
!#ifdef SPHERICAL
      if(spherical_pm) then
      if(lsingular_line) then
      if(ipolar(1).ne.0.and.ipolar(2).ne.0) then
        jj = jd + jl
        jstride = -1
      endif
      endif
      endif
!#endif /* SPHERICAL */
#endif /* NOTNOW */
          do j = ja,jb
          ii = id
          do i = ia,ib
            if (k >= ks .and. k <= ks + kl) then
!            if (j >= js .and. j <= js + jl) then
            if (j >= js2.and. j <= js1) then
            if (i >= is .and. i <= is + il) then
              if(.not.l_f_to_c) then
                work1(ii,jj,kk,idest) =
     .                 temprecv_buf(indx+1)
              else
                work1_fl(ii,jj,kk) =
     .                 temprecv_buf(indx+1)
              end if
            end if
            end if
            end if
            if (i >= is .and. i <= is + il) ii = ii + 1
            indx = indx+1
          enddo
!          if (j >= js .and. j <= js + jl) jj = jj + jstride
          if (j >= js2 .and. j <= js1) jj = jj + jstride
          enddo
          if (k >= ks .and. k <= ks + kl) kk = kk + 1
          enddo

        endif
      
        if (timing_mpix) then
         timer_amr_1blk_cc_cp_remote(3) = 
     .     timer_amr_1blk_cc_cp_remote(3) + mpi_wtime() - time1
        else
         timer_amr_1blk_cc_cp_remote(3) = -999.
        endif

      endif
#endif /* MPI_USED */


#ifdef TESTING_RESTRICT
!
! Record index ranges for controlling restriction of finelayer 
! data to fill guardcells where finer neighbors are found.
        if(l_f_to_c) then
!-------------------------
        if(iopt.eq.1) then
!-------------------------



        f2c_ind_unk(1,1,nblk_ind) = min( id,
     .                                   f2c_ind_unk(1,1,nblk_ind))
        f2c_ind_unk(2,1,nblk_ind) = max( id+il,
     .                                   f2c_ind_unk(2,1,nblk_ind))
        f2c_ind_unk(1,2,nblk_ind) = min( jd,
     .                                   f2c_ind_unk(1,2,nblk_ind))
        f2c_ind_unk(2,2,nblk_ind) = max( jd+jl,
     .                                   f2c_ind_unk(2,2,nblk_ind))
        f2c_ind_unk(1,3,nblk_ind) = min( kd,
     .                                   f2c_ind_unk(1,3,nblk_ind))
        f2c_ind_unk(2,3,nblk_ind) = max( kd+kl,
     .                                   f2c_ind_unk(2,3,nblk_ind))


!-------------------------
        elseif(iopt.gt.1) then
!-------------------------

        f2c_ind_work(1,1,nblk_ind) = min( id,
     .                                   f2c_ind_work(1,1,nblk_ind))
        f2c_ind_work(2,1,nblk_ind) = max( id+il,
     .                                   f2c_ind_work(2,1,nblk_ind))
        f2c_ind_work(1,2,nblk_ind) = min( jd,
     .                                   f2c_ind_work(1,2,nblk_ind))
        f2c_ind_work(2,2,nblk_ind) = max( jd+jl,
     .                                   f2c_ind_work(2,2,nblk_ind))
        f2c_ind_work(1,3,nblk_ind) = min( kd,
     .                                   f2c_ind_work(1,3,nblk_ind))
        f2c_ind_work(2,3,nblk_ind) = max( kd+kl,
     .                                   f2c_ind_work(2,3,nblk_ind))


!-------------------------
       endif                   ! end of iopt if test
!-------------------------
       endif                   ! end of l_f_to_c if test
#endif

      return
      end subroutine amr_1blk_cc_cp_remote





