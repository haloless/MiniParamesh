
#include "paramesh_preprocessor.fh"



      subroutine amr_initial_soln





!
! This file is a template describing how the solution can be
! initialized on the initial grid. Modify it for your own use.
!
!--------------------------------------------------------------
! include files for amr
      use paramesh_dimensions
      use physicaldata
      use tree

#ifdef SHMEM_USED
#include "amr_shmem.fh"
      include 'shmem_reduce.fh'
#endif

      integer :: nguard0

!--------------------------------------------------------------


      nguard0 = nguard*npgs

! loop over leaf grid blocks
      if(lnblocks.gt.0) then
      do lb=1,lnblocks

      if(nodetype(lb).eq.1 .or. advance_all_levels) then


        if(nvar.gt.0) then

! set values for unk
      dx = bsize(1,lb)/real(nxb)
      dy = bsize(2,lb)/real(nyb)
        do k=kl_bnd+nguard0*k3d,ku_bnd-nguard0*k3d
          do j=jl_bnd+nguard0*k2d,ju_bnd-nguard0*k2d
            do i=il_bnd+nguard0,iu_bnd-nguard0
              unk(1,i,j,k,lb) = 1.0
              xi =  bnd_box(1,1,lb) + dx*(real(i-nguard0)-.5)
              yi =  bnd_box(1,2,lb) + dy*(real(j-nguard0)-.5)
              if( abs(xi).lt.1.0 .and. abs(yi).lt.1.0) then
                          unk(1,i,j,k,lb) = 10.0
              endif
            enddo
          enddo
        enddo

        endif

        if(nvarcorn.gt.0) then

! set values for unk_n
        do k=kl_bnd+nguard0*k3d,ku_bnd+(nguard0+1)*k3d
          do j=jl_bnd+nguard0*k2d,ju_bnd+(nguard0+1)*k2d
            do i=il_bnd+nguard0,iu_bnd+(nguard0+1)
!              unk_n(1,i,j,k,lb) = ???
!              unk_n(2,i,j,k,lb) = ???
!              .
!              .
!              .
            enddo
          enddo
        enddo

        endif


        if(nfacevar.gt.0) then

! set values for facevarx
        do k=kl_bnd+nguard0*k3d,ku_bnd+nguard0*k3d
          do j=jl_bnd+nguard0*k2d,ju_bnd+nguard0*k2d
            do i=il_bnd+nguard0,iu_bnd+nguard0+1
!              facevarx(1,i,j,k,lb) = ???
!              facevarx(2,i,j,k,lb) = ???
!              .
!              .
!              .
            enddo
          enddo
        enddo

! set values for facevary
        do k=kl_bnd+nguard0*k3d,ku_bnd+nguard0*k3d
          do j=jl_bnd+nguard0*k2d,ju_bnd+(nguard0+1)*k2d
            do i=il_bnd+nguard0,iu_bnd+nguard0
!              facevary(1,i,j,k,lb) = ???
!              facevary(2,i,j,k,lb) = ???
!              .
!              .
!              .
            enddo
          enddo
        enddo

! set values for facevarz
        do k=kl_bnd+nguard0*k3d,ku_bnd+(nguard0+1)*k3d
          do j=jl_bnd+nguard0*k2d,ju_bnd+nguard0*k2d
            do i=il_bnd+nguard0,iu_bnd+nguard0
!              facevarz(1,i,j,k,lb) = ???
!              facevarz(2,i,j,k,lb) = ???
!              .
!              .
!              .
            enddo
          enddo
        enddo

        endif



        if(nvaredge.gt.0) then

! set values for unk_e_x
        do k=kl_bnd+nguard0*k3d,ku_bnd+(nguard0+1)*k3d
          do j=jl_bnd+nguard0*k2d,ju_bnd+(nguard0+1)*k2d
            do i=il_bnd+nguard0,iu_bnd+nguard0
!              unk_e_x(1,i,j,k,lb) = ???
!              unk_e_x(2,i,j,k,lb) = ???
!              .
!              .
!              .
            enddo
          enddo
        enddo

! set values for unk_e_y
        do k=kl_bnd+nguard0*k3d,ku_bnd+(nguard0+1)*k3d
          do j=jl_bnd+nguard0*k2d,ju_bnd+nguard0*k2d
            do i=il_bnd+nguard0,iu_bnd+(nguard0+1)
!              unk_e_y(1,i,j,k,lb) = ???
!              unk_e_y(2,i,j,k,lb) = ???
!              .
!              .
!              .
            enddo
          enddo
        enddo

! set values for unk_e_z
        do k=kl_bnd+nguard0*k3d,ku_bnd+nguard0*k3d
          do j=jl_bnd+nguard0*k2d,ju_bnd+(nguard0+1)*k2d
            do i=il_bnd+nguard0,iu_bnd+(nguard0+1)
!              unk_e_z(1,i,j,k,lb) = ???
!              unk_e_z(2,i,j,k,lb) = ???
!              .
!              .
!              .
            enddo
          enddo
        enddo

        endif

      endif

      enddo ! end loop over grid blocks
      endif

      return
      end subroutine amr_initial_soln
