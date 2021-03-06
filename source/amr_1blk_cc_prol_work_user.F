!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_1blk_cc_prol_work_user()


!
!------------------------------------------------------------------------
!
! This is a stub routine and is meant to be a place holder to allow
! a user to write their own prolocation routine for cell centered 
! data.  
!
! NOTE: To use this feature you must define interp_mask_unk to be >= 20.
!
! NOTE2: Use one of the other routines which are provided for doing this
! operation as an example.
!
!------------------------------------------------------------------------
!

      use paramesh_dimensions
      use physicaldata
      use tree
      use prolong_arrays

      implicit none

      include 'mpif.h'

      return
      end subroutine amr_1blk_cc_prol_work_user
