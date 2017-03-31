!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"
        subroutine amr_set_runtime_parameters()

        use paramesh_dimensions
        use physicaldata
        use tree
        use timings
        use io

        implicit none

        integer :: ierr

#ifdef LIBRARY
        open (unit=35,
     &        file='amr_runtime_parameters', 
     &        status='old',
     &        action='READ',
     &        form='formatted')
! integers
        read (35,*) maxblocks
        read (35,*) ndim
        read (35,*) l2p5d
        read (35,*) nxb
        read (35,*) nyb
        read (35,*) nzb
        read (35,*) nvar
        read (35,*) nfacevar
        read (35,*) nvaredge
        read (35,*) nvarcorn
        read (35,*) nvar_work
        read (35,*) nguard
        read (35,*) nguard_work
        read (35,*) nfluxvar
        read (35,*) nedgevar1
        read (35,*) iface_off
        read (35,*) mflags
        read (35,*) nfield_divf
        read (35,*) nboundaries
! logicals
        read (35,*) diagonals
        read (35,*) amr_error_checking
        read (35,*) no_permanent_guardcells
        read (35,*) advance_all_levels
        read (35,*) force_consistency
        read (35,*) consv_fluxes
        read (35,*) consv_flux_densities
        read (35,*) edge_value
        read (35,*) edge_value_integ
        read (35,*) var_dt
        read (35,*) pred_corr
        read (35,*) empty_cells
        read (35,*) conserve
        read (35,*) divergence_free
        read (35,*) curvilinear
        read (35,*) curvilinear_conserve
        read (35,*) cartesian_pm
        read (35,*) cylindrical_pm
        read (35,*) spherical_pm
        read (35,*) polar_pm
        read (35,*) lsingular_line
        read (35,*) timing_mpi
        read (35,*) timing_mpix
! characters
        read (35,*) output_dir

        close(35)
#else
        output_dir = OUTPUT_DIR
#endif

        amr_log_file = trim(output_dir) // 'amr.log'

        return
        end subroutine amr_set_runtime_parameters
