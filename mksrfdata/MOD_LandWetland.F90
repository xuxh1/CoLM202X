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
   type(grid_type) :: gwetland

   ! type(pixelset_type) :: landwetland
   integer,  allocatable :: wetlandclass (:)
   real(r8), allocatable :: pctshrpwh (:)

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

#if (defined SinglePoint && defined CH4)
      IF ((SITE_landtype == WETLAND) .and. (USE_SITE_pctwetland)) THEN

         numpatch = count(SITE_pctwetland > 0.)

         allocate (pctshrpwh(numpatch))
         allocate (wetlandclass(numpatch))
         wetlandclass = pack(SITE_wetlandtyp, SITE_pctwetland > 0.)
         pctshrpwh = pack(SITE_pctwetland, SITE_pctwetland > 0.)

         pctshrpwh = pctshrpwh / sum(pctshrpwh)

         IF (allocated(landpatch%eindex))  deallocate(landpatch%eindex)
         IF (allocated(landpatch%ipxstt))  deallocate(landpatch%ipxstt)
         IF (allocated(landpatch%ipxend))  deallocate(landpatch%ipxend)
         IF (allocated(landpatch%settyp))  deallocate(landpatch%settyp)
         IF (allocated(landpatch%ielm  ))  deallocate(landpatch%ielm  )

         allocate (landpatch%eindex (numpatch))
         allocate (landpatch%ipxstt (numpatch))
         allocate (landpatch%ipxend (numpatch))
         allocate (landpatch%settyp (numpatch))
         allocate (landpatch%ielm   (numpatch))

         landpatch%eindex(:) = 1
         landpatch%ielm  (:) = 1
         landpatch%ipxstt(:) = 1
         landpatch%ipxend(:) = 1
         landpatch%settyp(:) = WETLAND
         
         landpatch%has_shared = .true.
         allocate (landpatch%pctshared(numpatch))
         landpatch%pctshared = pctshrpwh

         landpatch%nset = numpatch
         CALL landpatch%set_vecgs

         RETURN
      ENDIF
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_io) THEN

         dir_5x5 = trim(DEF_dir_rawdata) // '/plant_15s'
         suffix  = 'MOD'//trim(cyear)

         CALL allocate_block_data (gpatch, pctwetland_xy)
         CALL read_5x5_data (dir_5x5, suffix, gpatch, 'PCT_WETLAND', pctwetland_xy)
         
         CALL allocate_block_data (gpatch, pctshared_xy, 2)
         DO iblkme = 1, gblock%nblkme
            ib = gblock%xblkme(iblkme)
            jb = gblock%yblkme(iblkme)
            pctshared_xy%blk(ib,jb)%val(1,:,:) = 1. - pctwetland_xy%blk(ib,jb)%val/100.
            pctshared_xy%blk(ib,jb)%val(2,:,:) = pctwetland_xy%blk(ib,jb)%val/100.
         ENDDO
      ENDIF
      
      sharedfilter = (/ 1 /)

      IF (landpatch%has_shared) then
         CALL pixelsetshared_build (landpatch, gpatch, pctshared_xy, 2, sharedfilter, &
            pctshared, classshared, fracin = landpatch%pctshared)
      ELSE
         CALL pixelsetshared_build (landpatch, gpatch, pctshared_xy, 2, sharedfilter, &
            pctshared, classshared)
      ENDIF

      IF (p_is_worker) THEN
         IF (landpatch%nset > 0) THEN
            WHERE (classshared == 2) landpatch%settyp = WETLAND
         ENDIF
      ENDIF

      IF (p_is_io) THEN
         file_patch = trim(DEF_dir_rawdata) // '/global_WFT_surface_data.nc'
         CALL allocate_block_data (gwetland, wetlanddata, N_WFT)
         CALL ncio_read_block (file_patch, 'PCT_WFT', gwetland, N_WFT, wetlanddata)
      ENDIF

      wetlandfilter = (/ WETLAND /)
      
      CALL pixelsetshared_build (landpatch, gwetland, wetlanddata, N_WFT, wetlandfilter, &
         pctshrpwh, wetlandclass, fracin = pctshared)

      wetlandclass = wetlandclass + N_PFT
      
      numpatch = landpatch%nset

      landpatch%has_shared = .true.
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            IF (allocated(landpatch%pctshared)) THEN
               deallocate(landpatch%pctshared)
            ENDIF

            allocate(landpatch%pctshared(numpatch))
            landpatch%pctshared = pctshrpwh
         ENDIF
      ENDIF

      IF (allocated(pctshared  )) deallocate(pctshared  )
      IF (allocated(classshared)) deallocate(classshared)

#ifdef USEMPI
      IF (p_is_worker) THEN
         CALL mpi_reduce (numpatch, npatch_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', npatch_glb, ' patches (with wetland).'
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numpatch, ' patches.'
#endif

      CALL elm_patch%build (landelm, landpatch, use_frac = .true.)
#ifdef CATCHMENT
      CALL hru_patch%build (landhru, landpatch, use_frac = .true.)
#endif

      CALL write_patchfrac (DEF_dir_landdata, lc_year)
   
   END SUBROUTINE landwetland_build

END MODULE MOD_LandWetland
#endif
