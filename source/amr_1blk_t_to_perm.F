!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

      subroutine amr_1blk_t_to_perm( lcc,lfc,lec,lnc,lb,idest)




!------------------------------------------------------------------------
!
! This routine copies data from the 1-block working arrays with guardcells
! to the permanent data arrays, which may or may not have permanent
! guardcells, depending on whether NO_PERMANENT_GUARDCELLS is defined 
! in physicaldata.fh.
!
!
! Arguments :
!      lcc          logical       copies cell centered data if true
!      lfc          logical       copies cell face-centered data if true
!      lec          logical       copies cell edge-centered data if true
!      lnc          logical       copies cell corner data if true
!      lb           integer       block into which data is to be copied
!      idest        integer       sets value for last dimension index
!                                  in the 1-blk data arrays
!
!
! Written :     Peter MacNeice          February 1999
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      implicit none

#include "amr_shmem.fh"


!------------------------------------

      integer, intent(in) :: lb,idest
      logical, intent(in) :: lcc,lfc,lec,lnc
      integer :: nguard0


!------------------------------------

      if (var_dt .or. pred_corr) then

      nguard0 = nguard*(1-npgs)

! cell-centered data
       if(lcc) then

           if (no_permanent_guardcells) then
           t_unk(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,lb)
     .       = unk1(:,il_bnd+nguard0:iu_bnd+nguard0,
     .               jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d,
     .               kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,idest)
           else
           t_unk(:,:,:,:,lb) = unk1(:,:,:,:,idest)
           end if
       endif                             ! end of lcc if test



! cell face-centered data
       if(lfc) then

! x-face
         if (no_permanent_guardcells) then
         tfacevarx(1:nfacevar,il_bnd:iu_bnd+1,
     .                       jl_bnd:ju_bnd,kl_bnd:ku_bnd,lb)
     .       = facevarx1(1:nfacevar,il_bnd+nguard0:iu_bnd+nguard0+1,
     .               jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d,
     .               kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,idest)
         else
         tfacevarx(1:nfacevar,:,:,:,lb) =
     .                        facevarx1(1:nfacevar,:,:,:,idest)
         endif

         if (ndim > 1) then
! y-face
         if (no_permanent_guardcells) then
         tfacevary(1:nfacevar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,
     .                         kl_bnd:ku_bnd,lb)
     .       = facevary1(1:nfacevar,il_bnd+nguard0:iu_bnd+nguard0,
     .               jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d+k2d,
     .               kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,idest)
         else
         tfacevary(1:nfacevar,:,:,:,lb) =
     .                        facevary1(1:nfacevar,:,:,:,idest)
         endif
         end if

         if (ndim == 3) then
! z-face
         if (no_permanent_guardcells) then
         tfacevarz(1:nfacevar,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .                         kl_bnd:ku_bnd+k3d,lb)
     .       = facevarz1(1:nfacevar,il_bnd+nguard0:iu_bnd+nguard0,
     .               jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d,
     .               kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d+k3d,idest)
         else
         tfacevarz(1:nfacevar,:,:,:,lb) =
     .                        facevarz1(1:nfacevar,:,:,:,idest)
         endif
         end if

        endif                     ! end of lfc if test


! cell edge-centered data
       if(lec) then
         if (ndim > 1) then
! x-edge
         if (no_permanent_guardcells) then
         t_unk_e_x(1:nvaredge,il_bnd:iu_bnd,
     .               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,lb)
     .       = unk_e_x1(1:nvaredge,il_bnd+nguard0:iu_bnd+nguard0,
     .               jl_bnd+nguard0*k2d:ju_bnd+(nguard0+1)*k2d,
     .               kl_bnd+nguard0*k3d:ku_bnd+(nguard0+1)*k3d,idest)
         else
         t_unk_e_x(1:nvaredge,:,:,:,lb) =
     .                        unk_e_x1(1:nvaredge,:,:,:,idest)
         endif
! y-edge
         if (no_permanent_guardcells) then
         t_unk_e_y(1:nvaredge,il_bnd:iu_bnd+1,
     .               jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d,lb)
     .       = unk_e_y1(1:nvaredge,il_bnd+nguard0:iu_bnd+1+nguard0,
     .               jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d,
     .               kl_bnd+nguard0*k3d:ku_bnd+(nguard0+1)*k3d,idest)
         else
         t_unk_e_y(1:nvaredge,:,:,:,lb) =
     .                        unk_e_y1(1:nvaredge,:,:,:,idest)
         endif

         if (ndim == 3) then
! z-edge
         if (no_permanent_guardcells) then
         t_unk_e_z(1:nvaredge,il_bnd:iu_bnd+1,
     .               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd,lb)
     .       = unk_e_z1(1:nvaredge,il_bnd+nguard0:iu_bnd+1+nguard0,
     .               jl_bnd+nguard0*k2d:ju_bnd+(nguard0+1)*k2d,
     .               kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,idest)
         else
         t_unk_e_z(1:nvaredge,:,:,:,lb) =
     .                        unk_e_z1(1:nvaredge,:,:,:,idest)
         endif
         end if

         end if

        endif                     ! end of lec if test

! cell corner data
       if(lnc) then
         if (no_permanent_guardcells) then
         t_unk_n(1:nvarcorn,il_bnd:iu_bnd+1,
     .               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,lb)
     .       = unk_n1(1:nvarcorn,il_bnd+nguard0:iu_bnd+1+nguard0,
     .               jl_bnd+nguard0*k2d:ju_bnd+(nguard0+1)*k2d,
     .               kl_bnd+nguard0*k3d:ku_bnd+(nguard0+1)*k3d,idest)
         else
         t_unk_n(1:nvarcorn,:,:,:,lb) =
     .                        unk_n1(1:nvarcorn,:,:,:,idest)
        endif
        endif                     ! end of lnc if test

        endif ! var_dt


      return
      end subroutine amr_1blk_t_to_perm
