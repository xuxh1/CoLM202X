#include <define.h>

#ifdef CH4
MODULE MOD_ch4_readin

	!-----------------------------------------------------------------------
	USE MOD_Precision
	IMPLICIT NONE
	SAVE

	! PUBLIC MEMBER FUNCTIONS:
	PUBLIC :: ch4_readin

CONTAINS

	SUBROUTINE ch4_readin ()
		! ===========================================================
		! ! DESCRIPTION:
		! Read in ch4 data.
		!
		! Original: Xionghui Xu, 2024
		! ===========================================================

		USE MOD_Precision
		USE MOD_Namelist
		USE MOD_SPMD_Task
		USE MOD_LandPatch
		USE MOD_NetCDFSerial
		USE MOD_NetCDFBlock
		USE MOD_SpatialMapping
		USE MOD_Vars_TimeInvariants
		USE MOD_Vars_TimeVariables

		USE MOD_Vars_Global
		! USE MOD_LandPFT
		! USE MOD_Vars_PFTimeVariables
		USE MOD_RangeCheck
		USE MOD_Block

		IMPLICIT NONE

		character(len=256) :: file_ch4
		type(grid_type)    :: grid_ch4
		type(block_data_real8_2d)  :: f_xy_ch4
		type(spatial_mapping_type) :: mg2patch_ch4

		! Local variables
		real(r8), allocatable :: lat(:), lon(:)
		real(r8) :: missing_value
		integer  :: npatch
		character(len=2) :: cx
		integer  :: iblkme, iblk, jblk
		integer  :: maxvalue, minvalue

      	! READ in ch4

		file_ch4 = trim(DEF_dir_runtime) // '/ch4/conc_ch4.nc'

		CALL ncio_read_bcast_serial (file_ch4, 'lat', lat)
		CALL ncio_read_bcast_serial (file_ch4, 'lon', lon)
  
		CALL grid_ch4%define_by_center (lat, lon)
  
		IF (p_is_io) THEN
		   CALL allocate_block_data  (grid_ch4, f_xy_ch4)
		ENDIF
  
		! missing value
		IF (p_is_master) THEN
		   CALL ncio_get_attr (file_ch4, 'conc_ch4', 'missing_value', missing_value)
		ENDIF

		IF (p_is_io) THEN
		   CALL ncio_read_block (file_ch4, 'conc_ch4', grid_ch4, f_xy_ch4)
		ENDIF
  
		CALL mg2patch_ch4%build_arealweighted (grid_ch4, landpatch)
		CALL mg2patch_ch4%set_missing_value   (f_xy_ch4, missing_value)
  
		IF (allocated(lon)) deallocate(lon)
		IF (allocated(lat)) deallocate(lat)
  
		IF (p_is_worker) THEN
		   IF (numpatch > 0)  allocate(conc_ch4_tmp   (numpatch))
		ENDIF
  
		! (1) Read in plant date for rice2.
		file_ch4 = trim(DEF_dir_runtime) // '/ch4/conc_ch4.nc'
		IF (p_is_io) THEN
			CALL ncio_read_block (file_ch4, 'conc_ch4', grid_ch4, f_xy_ch4)
		ENDIF
  
		CALL mg2patch_ch4%grid2pset (f_xy_ch4, conc_ch4_tmp)
  
		IF (p_is_worker) THEN
		   DO npatch = 1, numpatch
			  IF (conc_ch4_tmp(npatch) /= spval) THEN
				 conc_ch4 (npatch) = int(conc_ch4_tmp (npatch))
			  ELSE
				 conc_ch4 (npatch) = 0
			  ENDIF
		   ENDDO
		ENDIF
  
#ifdef RangeCheck
	CALL check_vector_data ('ch4 concentration in soil', conc_ch4)
#endif

	END SUBROUTINE ch4_readin

END MODULE MOD_ch4_readin
#endif

