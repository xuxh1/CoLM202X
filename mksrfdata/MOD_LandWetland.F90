#include <define.h>

#ifdef CH4
MODULE MOD_LandWetland

!------------------------------------------------------------------------------------
! DESCRIPTION:
!
!    Build wetland patches.
!
! Created by Xionghui Xu, Dec 2024
!    porting codes from Hua Yuan's OpenMP version to MPI parallel version.
!------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Grid
   IMPLICIT NONE

   ! ---- Instance ----
   type(grid_type) :: grid_wetland
   ! type(pixelset_type) :: landwetland
   integer,  allocatable :: wetlandclass (:)
   real(r8), allocatable :: wetlandfrac (:)

CONTAINS

   ! -------------------------------
   SUBROUTINE landwetland_build (lc_year)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Block
   USE MOD_DataType
   USE MOD_LandElm
#ifdef CATCHMENT
   USE MOD_LandHRU
#endif
   USE MOD_LandPatch
   USE MOD_NetCDFBlock
   USE MOD_PixelsetShared
   USE MOD_5x5DataReadin
#ifdef SinglePoint
   USE MOD_SingleSrfdata
#endif

   IMPLICIT NONE

   integer, intent(in) :: lc_year

   ! Local Variables
   character(len=255) :: cyear, file_patch, dir_5x5, suffix
   integer :: npatch_glb
   type(block_data_real8_2d) :: pctwetland_xy
   type(block_data_real8_3d) :: pctshared_xy
   type(block_data_real8_3d) :: wetlanddata
   integer :: sharedfilter(1), wetlandfilter(1)
   integer :: iblkme, ib, jb
   real(r8), allocatable :: pctshared  (:)
   integer , allocatable :: classshared(:)

      write(cyear,'(i4.4)') lc_year
      IF (p_is_master) THEN
         write(*,'(A)') 'Making patches (wetland shared) :'
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_io) THEN

         dir_5x5 = trim(DEF_dir_rawdata) // '/plant_15s'
         suffix  = 'MOD'//trim(cyear)

         CALL allocate_block_data (grid_patch, pctwetland_xy)
         CALL read_5x5_data (dir_5x5, suffix, grid_patch, 'PCT_WETLAND', pctwetland_xy)
         
         CALL allocate_block_data (grid_patch, pctshared_xy, 2)
         DO iblkme = 1, gblock%nblkme
            ib = gblock%xblkme(iblkme)
            jb = gblock%yblkme(iblkme)
            pctshared_xy%blk(ib,jb)%val(1,:,:) = 1. - pctwetland_xy%blk(ib,jb)%val/100.
            pctshared_xy%blk(ib,jb)%val(2,:,:) = pctwetland_xy%blk(ib,jb)%val/100.
         ENDDO
      ENDIF

      allocate (wetlandclass(numpatch))
      allocate (wetlandfrac(numpatch))

      wetlandclass = 1
      wetlandfrac = 1.
   
   END SUBROUTINE landwetland_build

END MODULE MOD_LandWetland
#endif
