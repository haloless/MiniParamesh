!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_refine_derefine
!! NAME
!!
!!   amr_refine_derefine
!!
!! SYNOPSIS
!!
!!   call amr_refine_derefine()
!!
!!
!! ARGUMENTS
!!
!!   NO ARGUMENTS
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
!!   timings
!!   io
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_check_refine
!!   amr_check_derefine
!!   amr_refine_blocks
!!   amr_derefine_blocks
!!   amr_morton_order
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit blocks marked for refinement have been 
!!   created, blocks marked for derefinement have been eliminated, and the blocks 
!!   have been reordered to acheive load balance.
!!
!! DESCRIPTION
!! 
!!  Subroutine to refine or derefine blocks.  This routine is called by the user
!!  who sets the refine or derefine flags to be true or false.  These flags are
!!  logical variables called 'refine' and 'derefine' and are stored for each block.
!!  If a block is marked for refinement, amr_refine_derefine will create its new 
!!  child blocks.  Also, tests will be executed to see if any other blocks need to 
!!  also refine (by calling amr_check_refine) to ensure that the a jump in refinement of
!!  more than one level is not created.
!!
!!  If a block is marked for derefinement, amr_refine_derefine first checks to make 
!!  that the block can derefine and not create a jump in refinement of more than one
!!  level.  A check is also run to check that all the siblings of the derefining
!!  block's siblings are also marked for derefinement.  If these tests succeed, the
!!  block is removed from the list of blocks.
!!
!!  Once these operations are completed, the routine 'amr_morton_order' is called
!!  and the tree data structure is reorganized to acheive load balance using a
!!  morton space filling curve.  After this routine is called, the routine 
!!  'amr_redist_blk' is called, which actually moves the block data into the correct
!!  positions in the morton order list of blocks.
!!
!!  Finally, the routine 'amr_morton_process' is called.  This routine computes the
!!  communications patterns needed for guardcell filling, restriction, and 
!!  prologation and stores them for later use. 
!!
!! AUTHORS
!!
!!  Kevin Olson (1997)
!!
!!***

#include "paramesh_preprocessor.fh"
!#define DEBUG_AMR
!#define DEBUG
      subroutine amr_refine_derefine


      use paramesh_dimensions
      use physicaldata
      use tree
#ifdef SAVE_MORTS
      use mpi_morton
#endif /* SAVE_MORTS */
      use timings
      use io

      use paramesh_interfaces, only : amr_check_refine,
     &                                amr_check_derefine,
     &                                amr_refine_blocks,
     &                                amr_derefine_blocks,
     &                                amr_morton_order

      use paramesh_mpi_interfaces, only : mpi_amr_singular_line


      implicit none

      include 'mpif.h'


! local variables and arrays

      integer :: lnblocks2,tot_blocks,tot_blocksa,icontinue
      integer :: icontinue_ref,icontinue_deref
      integer :: icontinue2,max_blocks
      integer :: min_blocks
      integer :: nprocs,mype
      integer :: i,l
      integer :: lnblocks_old

      integer :: istrategy

      logical :: l_move_solution
      logical,save :: first_call = .true.
      integer :: ierrorcode, ierr

      logical refinet(maxblocks_tr)

#ifdef SAVE_MORTS
      integer :: lb
      integer :: mort_neigh(6,3,3,3)
      real    :: xmin,ymin,zmin,xmax,ymax,zmax
#endif /* SAVE_MORTS */

      double precision :: time1

! ---------------------------------------------------------------------------

      call MPI_COMM_SIZE (MPI_COMM_WORLD,nprocs,ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,mype,ierr)

! error trap for lrefine_max and lrefine_min
      if(lrefine_max.lt.1.or.lrefine_max.gt.100) then
        write(*,*) 'PARAMESH error : lrefine_max has a bad value'
     .             ,lrefine_max
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      if(lrefine_min.lt.1.or.lrefine_min.gt.100.or.
     .   lrefine_min.gt.lrefine_max) then
        write(*,*) 'PARAMESH error : lrefine_min or lrefine_max ',
     .   'has a bad value : lrefine_min= ',lrefine_min,
     .   ' lrefine_max= ',lrefine_max
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif

#ifdef SAVE_MORTS
! Find the coordinate ranges

      if (timing_mpi) then
         time1 = mpi_wtime()
      endif

      call mpi_amr_global_domain_limits
      if (timing_mpi) then
      timer_amr_global_domain_limits = 
     .     timer_amr_global_domain_limits + mpi_wtime() - time1
      endif

! Compute xmin,ymin,zmin,xmax,ymax,zmax or get them from storage
      xmin = grid_xmin
      ymin = grid_ymin
      zmin = grid_zmin
      xmax = grid_xmax
      ymax = grid_ymax
      zmax = grid_zmax

      if(first_call) then
        do lb = 1,lnblocks
          call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
     .                          lperiodicx,lperiodicy,lperiodicz,
     .                          coord(:,lb),bsize(:,lb),ndim,
     .                          lrefine(lb),lrefine_max,mort_neigh)
          surr_morts(:,:,:,:,lb) = 
     .                      mort_neigh(:,:,2-k2d:2+k2d,2-k3d:2+k3d)
        enddo
      endif
#endif /* SAVE_MORTS */

! enforce refinement level limits
! first upper limit
      where(lrefine==lrefine_max) refine = .false.
! then lower limit
      where(lrefine==lrefine_min) derefine = .false.
! finally force grid to refine toward base level if too coarse
      where( (lrefine<lrefine_min) .and.
     .       (nodetype==1) ) refine = .true.

      refinet(1:lnblocks) = refine(1:lnblocks)
      newchild(:) = .FALSE.

! CHECK derefinements and refinements

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test to see if any refinements have been requested.
        icontinue=0
        icontinue_deref = 0
        if(lnblocks.gt.0) then
           do l = 1,lnblocks
              if(nodetype(l).eq.1.and.refine(l)) then
                 icontinue=1
                 goto 10
              endif
           enddo
        endif
10      continue
        call MPI_ALLREDUCE (icontinue,icontinue2,1,MPI_INTEGER,
     &                      MPI_MAX,MPI_COMM_WORLD,ierr)

        icontinue_ref = icontinue2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#ifdef SPHERICAL
      if(spherical_pm) then
      istrategy = 0
      if(ndim.eq.3.and.lsingular_line) istrategy = 1
      call mpi_amr_singular_line(istrategy,nprocs)
      endif
!#endif /* SPHERICAL */

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file=amr_log_file,status='unknown',
     &        position='append')
         write (30,*) ' starting CHECK_REFINE '
         close(30)
         print *, ' starting CHECK_REFINE '
      end if
#endif

      if (timing_mpi) then
         time1 = mpi_wtime()
      endif
      call amr_check_refine (nprocs,mype,icontinue_ref)
      if (timing_mpi) then
      timer_amr_check_refine =  timer_amr_check_refine
     .                          + mpi_wtime() - time1
      endif

 20   continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test to see if any derefinements have been requested
      icontinue=0
      icontinue_deref = 0
      if(lnblocks.gt.0) then
         do l=1,lnblocks
            if(nodetype(l).eq.1.and.derefine(l)) then
               icontinue=1
               goto 101
            endif
         enddo
      endif
 101  continue
      call MPI_ALLREDUCE (icontinue,icontinue2,1,MPI_INTEGER,
     &     MPI_MAX,MPI_COMM_WORLD,ierr)
      icontinue_deref = icontinue2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file=amr_log_file,status='unknown',
     &        position='append')
         write (30,*)' starting CHECK_DEREFINE '
         close(30)
         print *,' starting CHECK_DEREFINE '
      end if
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (timing_mpi) then
         time1 = mpi_wtime()
      endif
      if(icontinue_deref.gt.0) call amr_check_derefine (mype)
      if (timing_mpi) then
      timer_amr_check_derefine(0) =  timer_amr_check_derefine(0)
     .                          + mpi_wtime() - time1
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test to see if any derefinements have been requested which passed the tests
! in amr_check_derefine
      if(icontinue_deref.gt.0) then

      icontinue=0
      icontinue_deref = 0
      if(lnblocks.gt.0) then
         do l=1,lnblocks
            if(nodetype(l).eq.1.and.derefine(l)) then
               icontinue=1
               goto 102
            endif
         enddo
      endif
 102  continue
      call MPI_ALLREDUCE (icontinue,icontinue2,1,MPI_INTEGER,
     &     MPI_MAX,MPI_COMM_WORLD,ierr)
      icontinue_deref = icontinue2

      endif

      if(icontinue_ref.eq.0   .and. 
     &   icontinue_deref.eq.0 .and. 
     &   .not.first_call) return
      first_call = .false.
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! NOW Actually refine and derefine the mesh

#ifdef DEBUG_AMR
      if (mype.eq.0) then         
         open (unit=30,file=amr_log_file,status='unknown',
     &        position='append')
         write (30,*)' starting REFINE_BLOCKS'
         close(30)
         print *,' starting REFINE_BLOCKS'
      end if
#endif

      lnblocks_old = lnblocks
      if (timing_mpi) then
         time1 = mpi_wtime()
      endif
      if(icontinue_ref.gt.0) call amr_refine_blocks (nprocs,mype)
      if (timing_mpi) then
      timer_amr_refine_blocks =  timer_amr_refine_blocks
     .                          + mpi_wtime() - time1
      endif

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file=amr_log_file,status='unknown',
     &        position='append')
         write (30,*) ' starting DEREFINE_BLOCKS'
         close(30)
         print *, ' starting DEREFINE_BLOCKS'
      end if
#endif     

      if (timing_mpi) then
         time1 = mpi_wtime()
      endif
      if(icontinue_deref.gt.0.or.icontinue_ref.gt.0) 
     .               call amr_derefine_blocks(lnblocks_old,mype)

      if (timing_mpi) then
      timer_amr_derefine_blocks =  timer_amr_derefine_blocks
     .                          + mpi_wtime() - time1
      endif


#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file=amr_log_file,status='unknown',
     &        position='append')
         close(30)
      end if
#endif

      call MPI_ALLREDUCE (lnblocks,tot_blocks,1,MPI_INTEGER,
     &                    MPI_SUM,MPI_COMM_WORLD,ierr)


#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file=amr_log_file,status='unknown',
     &        position='append')
         write(30,*) ' tot_blocks before ',tot_blocks
         close(30)
         print *,' tot_blocks before ',tot_blocks
      end if

! I copy lnblocks to lnblocks2 since lnblocks2 can be put in a save statement.
      lnblocks2 = lnblocks 
      call MPI_ALLREDUCE (lnblocks2,max_blocks,1,MPI_INTEGER,
     &                    MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE (lnblocks2,min_blocks,1,MPI_INTEGER,
     &                    MPI_MIN,MPI_COMM_WORLD,ierr)

      if (mype.eq.0) then
         open (unit=30,file=amr_log_file,status='unknown',
     &        position='append')
         write (30,*) ' max_blocks 1',max_blocks
         write (30,*) ' min_blocks 1',min_blocks
         print *, ' max_blocks 1',max_blocks
         print *, ' min_blocks 1',min_blocks
         close(30)
      end if
#endif

! set work values

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file=amr_log_file,status='unknown',
     &        position='append')
         write (30,*) ' starting MORTON ORDERING'
         print *, ' starting MORTON ORDERING '
         close(30)
      end if
#endif
      
      work_block(:) = 0.
      do i = 1,lnblocks
         if (nodetype(i).eq.1) work_block(i) = 2.        !<<< USER EDIT
         if (nodetype(i).ge.2) work_block(i) = 1.        !<<< USER EDIT
      end do
      l_move_solution = .true.
      if (timing_mpi) then
         time1 = mpi_wtime()
      endif
      call amr_morton_order (lnblocks_old,nprocs,mype,
     .                       l_move_solution)

      if (timing_mpi) then
      timer_amr_morton_order =  timer_amr_morton_order
     .                          + mpi_wtime() - time1
      endif
#ifdef DEBUG_AMR
      if (mype == 0) print *,' exited amr_morton_order'
#endif

      
! I copy lnblocks to lnblocks2 since lnblocks2 can be put in a save statement.
      lnblocks2 = lnblocks 
      call MPI_ALLREDUCE (lnblocks2,tot_blocksa,1,MPI_INTEGER,
     &                    MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE (lnblocks2,max_blocks,1,MPI_INTEGER,
     &                    MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE (lnblocks2,min_blocks,1,MPI_INTEGER,
     &                    MPI_MIN,MPI_COMM_WORLD,ierr)

#ifdef DEBUG_AMR
      if (mype == 0) print *,' exited MPI_ALLREDUCE'
#endif

      if (mype.eq.0) then
#ifdef DEBUG_AMR
         open (unit=30,file=amr_log_file,status='unknown',
     &        position='append')
         write (30,*) ' tot_blocks after ',tot_blocksa
         write (30,*) ' max_blocks 2',max_blocks
         write (30,*) ' min_blocks 2',min_blocks
         print *, ' tot_blocks after ',tot_blocksa
         print *, ' max_blocks 2',max_blocks
         print *, ' min_blocks 2',min_blocks
         close (30)
#endif
      end if

#ifdef DEBUG_AMR
      if (tot_blocksa.ne.tot_blocks) then
         print *,' ERROR: tot_blocksa.ne.tot_blocks ',
     &        tot_blocksa,tot_blocks
         call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
      end if
#endif

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file=amr_log_file,status='unknown',
     &        position='append')
         write (30,*) ' done REFINE_DEREFINE '
         write (30,*) ' '
         print *, ' done REFINE_DEREFINE '
         print *,' '
         close(30)
      end if
#endif

      if (timing_mpi) then
         time1 = mpi_wtime()
      endif
      call amr_morton_process()
      if (timing_mpi) then
      timer_amr_morton_process =  timer_amr_morton_process
     .                          + mpi_wtime() - time1
      endif

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file=amr_log_file,status='unknown',
     &        position='append')
         write (30,*) ' done MORTON_PROCESS '
         write (30,*) ' '
         print *, ' done MORTON_PROCESS '
         print *,' '
         close(30)
      end if
#endif


!---------------------------------------------
!
! St up an array of cell sizes for each grid refinement level.
! These can be used to minimize variation due to roundoff, but
! should ONLY be used with a uniformly spaced grid.
      level_cell_sizes = 0.
      level_cell_sizes(1,1) = (grid_xmax-grid_xmin)/real(nxb)
      if(ndim.gt.1)
     .  level_cell_sizes(2,1) = (grid_ymax-grid_ymin)/real(nyb)
      if(ndim.eq.3)
     .  level_cell_sizes(3,1) = (grid_zmax-grid_zmin)/real(nzb)
      do i=2,lrefine_max
        level_cell_sizes(1:ndim,i) = .5*level_cell_sizes(1:ndim,i-1)
      enddo
!---------------------------------------------

!
! set grid modification flag
      grid_changed = 1
      grid_analysed_mpi = 1

#ifdef DEBUG
      write(*,*) 'exiting amr_refine_derefine : pe ',mype
#endif /* DEBUG */




#ifdef MPI_USED
      if (timing_mpi) then
         time1 = mpi_wtime()
      endif
      call mpi_amr_boundary_block_info(mype,nprocs)
      if (timing_mpi) then
      timer_amr_boundary_block_info =  timer_amr_boundary_block_info
     .                          + mpi_wtime() - time1
      endif
#endif /* MPI_USED */


#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file=amr_log_file,status='unknown',
     &        position='append')
         write (30,*) ' done mpi_amr_boundary_block_info '
         write (30,*) ' '
         print *, ' done mpi_amr_boundary_block_info '
         print *,' '
         close(30)
      end if
#endif


      return
      end subroutine amr_refine_derefine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_morton_process()

      use paramesh_dimensions
      use physicaldata
      use tree
      use timings
      use io


      use paramesh_mpi_interfaces, only : mpi_morton_bnd,
     &                                    mpi_morton_bnd_prolong1,
     &                                    mpi_morton_bnd_fluxcon,
     &                                    mpi_morton_bnd_restrict,
     &                                    mpi_mort_comm_for_surrblks,
     &                                    mpi_amr_morton_limits,
     &                                    mpi_amr_gsurr_blks,
     &                                    mpi_setup


      implicit none

      include 'mpif.h'

      integer :: nprocs,mype,tag_offset,ierr
      logical :: lec, lnc, lfulltree

      double precision :: time1

#ifdef DEBUG_AMR
      write(*,*) 'entered mort_process '
#endif /* DEBUG */

      call MPI_COMM_SIZE (MPI_COMM_WORLD,nprocs,ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,mype,ierr)


#ifdef DEBUG_AMR
      write(*,*) 'before mort_comm_for_surrblks : pe ',mype
#endif /* DEBUG */

!---------------------------
! call setup routines in preparation for calling 
! mpi_mort_comm_for_surrblks all the mpi_morton_bnd_XXX
! routines.


! Find the coordinate ranges
      if (timing_mpi) then
         time1 = mpi_wtime()
      endif
      call mpi_amr_global_domain_limits
      if (timing_mpi) then
      timer_amr_global_domain_limits = 
     .     timer_amr_global_domain_limits + mpi_wtime() - time1
      endif

      if (timing_mpi) then
         time1 = mpi_wtime()
      endif
      call mpi_setup(mype,nprocs)
      if (timing_mpi) then
      timer_mpi_setup = 
     .     timer_mpi_setup + mpi_wtime() - time1
      endif
!
! Compute and save morton number range for each processor
      if (timing_mpi) then
         time1 = mpi_wtime()
      endif
      call mpi_amr_morton_limits(mype)
      if (timing_mpi) then
      timer_amr_morton_limits = 
     .     timer_amr_morton_limits + mpi_wtime() - time1
      endif
!---------------------------
#ifdef DEBUG_AMR
      write(*,*) 'before mort_comm_for_surrblks : pe ',mype
#endif /* DEBUG */
!
! Set up surrounding blocks of all local blocks (must not precede
! setting of grid_xmin,... etc)
      tag_offset = 100
      if (timing_mpi) then
         time1 = mpi_wtime()
      endif
      call mpi_mort_comm_for_surrblks(mype,nprocs,tag_offset)
      if (timing_mpi) then
      timer_mort_comm_for_surrblks = 
     .     timer_mort_comm_for_surrblks + mpi_wtime() - time1
      endif
#ifdef DEBUG_AMR
      write(*,*) 'after mort_comm_for_surrblks : pe ',mype
#endif /* DEBUG */

! call mpi_amr_gsurr_blks immediately after call to mpi_morton_bnd
! because it uses pe_source and r_mortonbnd which are reset in the
! other morton_bnd_?? routines. Note this is temporary. Will
! redesign so mpi_amr_gsurr_blks is less context sensitive.

#ifdef DEBUG_AMR
      write(*,*) 'before gsurr : pe ',mype
#endif /* DEBUG */
      call mpi_amr_gsurr_blks(mype,nprocs)
#ifdef DEBUG_AMR
      write(*,*) 'after gsurr : pe ',mype
#endif /* DEBUG */

      tag_offset = 100
#ifdef DEBUG_AMR
      write(*,*) 'before morton_bnd : pe ',mype
#endif /* DEBUG */
      call mpi_morton_bnd(mype,nprocs,tag_offset)
#ifdef DEBUG_AMR
      write(*,*) 'after morton_bnd : pe ',mype
#endif /* DEBUG */

      tag_offset = 100
#ifdef DEBUG_AMR
      write(*,*) 'before morton_bnd_prol : pe ',mype
#endif /* DEBUG */
      call mpi_morton_bnd_prolong1(mype,nprocs,tag_offset)
#ifdef DEBUG_AMR
      write(*,*) 'after morton_bnd_prol : pe ',mype
      call amr_flush(6)
      call shmem_barrier_all()
#endif /* DEBUG */


      tag_offset = 100
#ifdef DEBUG_AMR
      write(*,*) 'before morton_bnd_flux : pe ',mype
      call amr_flush(6)
      call shmem_barrier_all()
#endif /* DEBUG */
      call mpi_morton_bnd_fluxcon(mype,nprocs,tag_offset)
#ifdef DEBUG_AMR
      write(*,*) 'after morton_bnd_flux : pe ',mype
      call amr_flush(6)
      call shmem_barrier_all()
#endif /* DEBUG */


      lec = .false.
      lnc = .false.
      if(nvaredge.gt.0) lec = .true.
      if(nvarcorn.gt.0) lnc = .true.
      tag_offset = 100
      lfulltree = .false.
#ifdef DEBUG_AMR
      write(*,*) 'before morton_bnd_restr : pe ',mype
#endif /* DEBUG */
      call mpi_morton_bnd_restrict(mype,nprocs,
     .                             lfulltree,lec,lnc,tag_offset)


!        lfulltree = .true.
!        call mpi_morton_bnd_restrict(mype,nprocs,
!     .                               lfulltree,lec,lnc,tag_offset)
#ifdef DEBUG_AMR
      write(*,*) 'after morton_bnd_restr : pe ',mype
#endif /* DEBUG */

#ifdef DEBUG
      write(*,*) 'exiting amr_morton_process : pe ',mype
#endif /* DEBUG */

      return
      end subroutine amr_morton_process
