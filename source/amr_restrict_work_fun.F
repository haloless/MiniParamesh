!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_restrict_work_fun(datain,dataout,iopt)




!------------------------------------------------------------------------
!
! This routine performs restriction on the array datain and
! returns the result in dataout. Note that this does not updata
! guard cell elements of dataout.
!
! This particular version is only appropriate for 2nd order schemes 
! using linear interpolation with even number of mesh points along 
! each block axis.
!
! Written :     Peter MacNeice          January 1997
!------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata
      use workspace

      use paramesh_interfaces, only : amr_restrict_work_genorder,
     .                                amr_restrict_work_user

      implicit none

      real, intent(in)    :: datain(:,:,:)
      real, intent(inout) :: dataout(:,:,:)
      integer, intent(in) :: iopt

      integer :: order

!------------------------------------


      if (interp_mask_work_res(iopt-1) < 20) then
         
! default interpolation routine for restriction of 'work' array

         order = interp_mask_work_res(iopt-1)
         if (order <= 0 .or. order > 5) order = 1
         call amr_restrict_work_genorder (datain,dataout,iopt,order)

      else

! user defined routine for restriction of 'work' array
         call amr_restrict_work_user()

      end if

      return
      end subroutine amr_restrict_work_fun


