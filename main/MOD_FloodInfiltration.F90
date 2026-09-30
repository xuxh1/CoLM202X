#include <define.h>

MODULE MOD_FloodInfiltration
!-----------------------------------------------------------------------
! Re-infiltration of floodplain water into the soil of
! flooded soil patches when, in a single point, a site flood series floods
! them. The soil side follows the CaMa-Flood re-infiltration (LWINFILT)
! block of WATER_VSF. Gridded runs use DEF_GridRiverLake_FloodFeedback.
!-----------------------------------------------------------------------
   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   real(r8), allocatable :: fld_frc_p    (:)  ! flooded fraction of the patch [-]
   real(r8), allocatable :: fld_dph_p    (:)  ! flood water depth over the flooded part [mm]
   real(r8), allocatable :: fld_qinfl_p  (:)  ! patch-mean re-infiltration of this step [mm/s]
   ! Largest flooded fraction of the patch the land may see,
   ! filled by the methane module from the GLWD floodplain classes of its
   ! cell; it bounds fld_frc_p and, in gridded runs, the flood fraction
   ! the grid-based routing publishes to the patches.
   ! Not allocated: no bound.
   real(r8), allocatable :: fld_cap_p    (:)  ! upper bound of fld_frc_p [-]

CONTAINS

   SUBROUTINE flood_infil_alloc (npatch)
      integer, intent(in) :: npatch
      IF (.not. allocated(fld_frc_p)) THEN
         allocate (fld_frc_p(max(npatch,0)), fld_dph_p(max(npatch,0)), fld_qinfl_p(max(npatch,0)))
         fld_frc_p   = 0._r8
         fld_dph_p   = 0._r8
         fld_qinfl_p = 0._r8
      ENDIF
   END SUBROUTINE flood_infil_alloc

END MODULE MOD_FloodInfiltration
