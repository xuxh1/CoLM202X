#include <define.h>

#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Reactive_Methane_WetlandVeg

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Vegetation of the permanent-wetland tile from its GLWD v2 make-up,
!   replacing the five-zone climate proxy when
!   DEF_METHANE%wetland_veg_glwd is true.
!
!   Per patch, from the grid cell nearest to the patch centre:
!   - wetveg_forest: forested share of the cell's wetland-tile classes,
!     (16, 18, 22, 24, 26) / (16-19, 22-27);
!   - wetveg_laicap: peak LAI of the non-forested share, the area mean of
!     wetland_lai_open_peat over the open-peatland classes (23, 25) and
!     wetland_lai_marsh over the marsh classes (17, 19, 27).
!   A cell without any tile class gives forest share 0 and no LAI cap.
!   With DEF_METHANE%wetland_lai_shape, wetveg_laipeak holds the annual peak
!   of the patch's remote-sensing LAI for the LAI year in use, and the
!   non-forested share is scaled by cap / peak instead of clipped at the cap.
!   With DEF_METHANE%wetland_dim_class also the process-dimension
!   shares of the tile classes (16-19, 22-27): wetveg_emerg, emergent marsh
!   (17); wetveg_dome, tropical peat dome (26, 27); wetveg_moss, moss-carpeted
!   peatland (22-25). A cell without any tile class gives 0 for all three.
!   The water-source dimension D2 is the rain-fed share
!   s_r = s_p + s_b, fed by precipitation alone, held as its permafrost
!   peat plateau part s_p (wetveg_ombro, the share without lateral inflow)
!   and its open-bog part s_b (wetveg_bog); globally rainfed_share and
!   bog_share of the file, at a tower wetland_share_rain_site and
!   wetland_share_bog_site.
!
! INPUT FILE: DEF_METHANE%wetland_veg_file, a regular lat/lon grid with
!   variables lat, lon, forested_share and area_class_NN (km2), as written
!   by v2/scripts/glwd_class_2deg.py; optionally rainfed_share and
!   bog_share, as written by v2/scripts/mk_rainfed_share.py
!   and v2/scripts/mk_bog_share.py, 0 when the file has none.
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: PI
   USE, INTRINSIC :: IEEE_ARITHMETIC, only: ieee_is_finite

   IMPLICIT NONE
   SAVE
   PRIVATE

   real(r8), allocatable, public :: wetveg_forest(:)   ! forested share [-]
   real(r8), allocatable, public :: wetveg_ombro(:)    ! share without lateral inflow (permafrost bog, PEB) [-]
   real(r8), allocatable, public :: wetveg_bog(:)      ! open-bog share [-]
   real(r8), allocatable, public :: wetveg_emerg(:)    ! emergent marsh share (GLWD 17) [-]
   real(r8), allocatable, public :: wetveg_dome(:)     ! tropical peat dome share (GLWD 26-27) [-]
   real(r8), allocatable, public :: wetveg_moss(:)     ! moss-carpeted peatland share (GLWD 22-25) [-]
   real(r8), allocatable, public :: wetveg_laicap(:)   ! LAI cap of the non-forested share [m2/m2]
   logical,  public :: wetveg_active = .false.
   real(r8), allocatable, public :: wetveg_laipeak(:)  ! annual peak of the remote-sensing LAI [m2/m2]
   integer :: wetveg_peak_year = -huge(1)              ! LAI year wetveg_laipeak belongs to

   PUBLIC :: read_methane_wetveg
   PUBLIC :: deallocate_methane_wetveg
   PUBLIC :: wetveg_cap_lai
   PUBLIC :: wetveg_herb_lai_ratio
   PUBLIC :: wetveg_forest_input_ratio
   PUBLIC :: wetveg_lai_peak_needed
   PUBLIC :: wetveg_track_lai_peak
   PUBLIC :: wetveg_bg_frac
   PUBLIC :: wetland_root_profile
   PUBLIC :: set_wetland_veg_glwd

CONTAINS

   SUBROUTINE deallocate_methane_wetveg ()
      IF (allocated(wetveg_forest)) deallocate(wetveg_forest)
      IF (allocated(wetveg_ombro))  deallocate(wetveg_ombro)
      IF (allocated(wetveg_bog))    deallocate(wetveg_bog)
      IF (allocated(wetveg_emerg))  deallocate(wetveg_emerg)
      IF (allocated(wetveg_dome))   deallocate(wetveg_dome)
      IF (allocated(wetveg_moss))   deallocate(wetveg_moss)
      IF (allocated(wetveg_laicap)) deallocate(wetveg_laicap)
      IF (allocated(wetveg_laipeak)) deallocate(wetveg_laipeak)
      wetveg_peak_year = -huge(1)
      wetveg_active = .false.
   END SUBROUTINE deallocate_methane_wetveg

   SUBROUTINE read_methane_wetveg (file_veg, lai_open_peat, lai_marsh, patchlatr_in, patchlonr_in, numpatch)
#ifdef USEMPI
      USE MPI
#endif
      USE netcdf
      USE MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE

      character(len=*), intent(in) :: file_veg
      real(r8), intent(in) :: lai_open_peat, lai_marsh
      real(r8), intent(in) :: patchlatr_in(:), patchlonr_in(:)   ! radians
      integer,  intent(in) :: numpatch

      integer, parameter :: nopen = 2, nmarsh = 3
      integer, parameter :: open_classes(nopen) = (/23, 25/), marsh_classes(nmarsh) = (/17, 19, 27/)
      integer, parameter :: ntile = 10
      integer, parameter :: tile_classes(ntile) = (/16, 17, 18, 19, 22, 23, 24, 25, 26, 27/)
      integer :: ncid, vid, ierr, nlat, nlon, ilat, ilon, ip, k, dims(2), bad, bad_dim
      real(r8), allocatable :: lat_g(:), lon_g(:), forest_g(:,:), cap_g(:,:), a(:,:), aopen(:,:), amarsh(:,:)
      real(r8), allocatable :: ombro_g(:,:), bog_g(:,:)
      real(r8), allocatable :: emerg_g(:,:), dome_g(:,:), moss_g(:,:), atile(:,:), acls(:,:)
      real(r8) :: lat_deg, lon_deg, d, dmin, dlon, srain
      character(len=16) :: vname

      CALL deallocate_methane_wetveg ()
      bad = 0
      bad_dim = 0
      dims = 0

      IF (p_is_master) THEN
         ierr = nf90_open(trim(file_veg), NF90_NOWRITE, ncid)
         IF (ierr /= NF90_NOERR) THEN
            write(*,'(A,A,A,A)') ' ERROR: wetland vegetation file ', trim(file_veg), ': ', trim(nf90_strerror(ierr))
            bad = 1
         ELSE
            ierr = nf90_inq_dimid(ncid, 'lat', vid)
            IF (ierr == NF90_NOERR) ierr = nf90_inquire_dimension(ncid, vid, len = dims(1))
            IF (ierr == NF90_NOERR) ierr = nf90_inq_dimid(ncid, 'lon', vid)
            IF (ierr == NF90_NOERR) ierr = nf90_inquire_dimension(ncid, vid, len = dims(2))
            IF (ierr /= NF90_NOERR .or. dims(1) <= 0 .or. dims(2) <= 0) bad = 1
         ENDIF
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (bad,  1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (dims, 2, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
      IF (bad /= 0) CALL CoLM_stop (' ***** ERROR: cannot read DEF_METHANE%wetland_veg_file.')
      nlat = dims(1); nlon = dims(2)
      allocate (lat_g(nlat), lon_g(nlon), forest_g(nlat,nlon), cap_g(nlat,nlon), ombro_g(nlat,nlon))
      allocate (bog_g(nlat,nlon))
      ombro_g = 0._r8
      bog_g   = 0._r8
      IF (DEF_METHANE%wetland_dim_class) THEN
         allocate (emerg_g(nlat,nlon), dome_g(nlat,nlon), moss_g(nlat,nlon))
         emerg_g = 0._r8
         dome_g  = 0._r8
         moss_g  = 0._r8
      ENDIF

      IF (p_is_master) THEN
         allocate (a(nlon,nlat), aopen(nlat,nlon), amarsh(nlat,nlon))
         ! Every variable must be read; a failed read would otherwise leave the
         ! previous class's values in a and change LAI and CH4 transport silently.
         vname = 'lat'
         ierr = nf90_inq_varid(ncid, trim(vname), vid)
         IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, lat_g)
         IF (ierr == NF90_NOERR) THEN
            vname = 'lon'
            ierr = nf90_inq_varid(ncid, trim(vname), vid)
            IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, lon_g)
         ENDIF
         IF (ierr == NF90_NOERR) THEN
            vname = 'forested_share'
            ierr = nf90_inq_varid(ncid, trim(vname), vid)
            IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, a)
            IF (ierr == NF90_NOERR) forest_g = transpose(a)
         ENDIF
         ! rain-fed share of the wetland classes, when the file has it
         IF (ierr == NF90_NOERR) THEN
            IF (nf90_inq_varid(ncid, 'rainfed_share', vid) == NF90_NOERR) THEN
               vname = 'rainfed_share'
               ierr = nf90_get_var(ncid, vid, a)
               IF (ierr == NF90_NOERR) ombro_g = transpose(a)
            ENDIF
         ENDIF
         ! open-bog share of the wetland classes, when the file has it
         IF (ierr == NF90_NOERR .and. nf90_inq_varid(ncid, 'bog_share', vid) == NF90_NOERR) THEN
            vname = 'bog_share'
            ierr = nf90_get_var(ncid, vid, a)
            IF (ierr == NF90_NOERR) bog_g = transpose(a)
         ENDIF
         aopen = 0._r8; amarsh = 0._r8
         DO k = 1, nopen
            IF (ierr /= NF90_NOERR) EXIT
            write(vname,'(A,I2.2)') 'area_class_', open_classes(k)
            ierr = nf90_inq_varid(ncid, trim(vname), vid)
            IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, a)
            IF (ierr == NF90_NOERR) aopen = aopen + max(transpose(a), 0._r8)
         ENDDO
         DO k = 1, nmarsh
            IF (ierr /= NF90_NOERR) EXIT
            write(vname,'(A,I2.2)') 'area_class_', marsh_classes(k)
            ierr = nf90_inq_varid(ncid, trim(vname), vid)
            IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, a)
            IF (ierr == NF90_NOERR) amarsh = amarsh + max(transpose(a), 0._r8)
         ENDDO
         IF (DEF_METHANE%wetland_dim_class) THEN
            ! area of each dimension over the area of the tile classes
            allocate (atile(nlat,nlon), acls(nlat,nlon))
            atile = 0._r8
            DO k = 1, ntile
               write(vname,'(A,I2.2)') 'area_class_', tile_classes(k)
               ierr = nf90_inq_varid(ncid, trim(vname), vid)
               IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, a)
               IF (ierr /= NF90_NOERR) THEN
                  bad_dim = 1
                  EXIT
               ENDIF
               acls = transpose(a)
               WHERE (.not. ieee_is_finite(acls) .or. acls < 0._r8 .or. acls > 1.e30_r8) acls = 0._r8
               atile = atile + acls
               SELECT CASE (tile_classes(k))
               CASE (17)
                  emerg_g = emerg_g + acls
               CASE (22:25)
                  moss_g  = moss_g + acls
               CASE (26:27)
                  dome_g  = dome_g + acls
               END SELECT
            ENDDO
            WHERE (atile > 0._r8)
               emerg_g = emerg_g / atile
               dome_g  = dome_g  / atile
               moss_g  = moss_g  / atile
            ELSEWHERE
               emerg_g = 0._r8
               dome_g  = 0._r8
               moss_g  = 0._r8
            END WHERE
            deallocate (atile, acls)
         ENDIF
         IF (ierr /= NF90_NOERR) THEN
            write(*,'(A,A,A,A,A,A)') ' ERROR: wetland vegetation file ', trim(file_veg), &
               ', variable ', trim(vname), ': ', trim(nf90_strerror(ierr))
            bad = 1
         ENDIF
         ierr = nf90_close(ncid)
         ! netCDF fill for cells without tile classes: no forest, no cap
         WHERE (.not. ieee_is_finite(forest_g) .or. forest_g < 0._r8 .or. forest_g > 1._r8) forest_g = 0._r8
         WHERE (.not. ieee_is_finite(ombro_g) .or. ombro_g < 0._r8 .or. ombro_g > 1._r8) ombro_g = 0._r8
         WHERE (.not. ieee_is_finite(bog_g) .or. bog_g < 0._r8 .or. bog_g > 1._r8) bog_g = 0._r8
         WHERE (aopen + amarsh > 0._r8)
            cap_g = (aopen * lai_open_peat + amarsh * lai_marsh) / (aopen + amarsh)
         ELSEWHERE
            cap_g = huge(1._r8)
         END WHERE
         deallocate (a, aopen, amarsh)
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (bad,      1,         MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
      IF (bad /= 0) CALL CoLM_stop (' ***** ERROR: cannot read DEF_METHANE%wetland_veg_file.')
#ifdef USEMPI
      CALL mpi_bcast (lat_g,    nlat,      MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (lon_g,    nlon,      MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (forest_g, nlat*nlon, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (cap_g,    nlat*nlon, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (ombro_g,  nlat*nlon, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (bog_g,    nlat*nlon, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      IF (DEF_METHANE%wetland_dim_class) THEN
         CALL mpi_bcast (bad_dim, 1,         MPI_INTEGER, p_address_master, p_comm_glb, p_err)
         CALL mpi_bcast (emerg_g, nlat*nlon, MPI_REAL8,   p_address_master, p_comm_glb, p_err)
         CALL mpi_bcast (dome_g,  nlat*nlon, MPI_REAL8,   p_address_master, p_comm_glb, p_err)
         CALL mpi_bcast (moss_g,  nlat*nlon, MPI_REAL8,   p_address_master, p_comm_glb, p_err)
      ENDIF
#endif
      IF (bad_dim /= 0) CALL CoLM_stop (' ***** ERROR: wetland_dim_class needs area_class_16 to 19 '// &
         'and 22 to 27 in DEF_METHANE%wetland_veg_file.')

      allocate (wetveg_forest(max(numpatch,0)), wetveg_laicap(max(numpatch,0)), wetveg_ombro(max(numpatch,0)))
      allocate (wetveg_bog(max(numpatch,0)))
      IF (DEF_METHANE%wetland_dim_class) &
         allocate (wetveg_emerg(max(numpatch,0)), wetveg_dome(max(numpatch,0)), wetveg_moss(max(numpatch,0)))
      DO ip = 1, numpatch
         lat_deg = patchlatr_in(ip) * 180._r8 / PI
         lon_deg = patchlonr_in(ip) * 180._r8 / PI
         dmin = huge(1._r8); ilat = 1
         DO k = 1, nlat
            d = abs(lat_g(k) - lat_deg)
            IF (d < dmin) THEN; dmin = d; ilat = k; ENDIF
         ENDDO
         dmin = huge(1._r8); ilon = 1
         DO k = 1, nlon
            dlon = abs(modulo(lon_g(k) - lon_deg + 180._r8, 360._r8) - 180._r8)
            IF (dlon < dmin) THEN; dmin = dlon; ilon = k; ENDIF
         ENDDO
         wetveg_forest(ip) = forest_g(ilat, ilon)
         wetveg_laicap(ip) = cap_g(ilat, ilon)
         wetveg_ombro(ip)  = ombro_g(ilat, ilon)
         wetveg_bog(ip)    = bog_g(ilat, ilon)
         IF (DEF_METHANE%wetland_dim_class) THEN
            wetveg_emerg(ip) = emerg_g(ilat, ilon)
            wetveg_dome(ip)  = dome_g(ilat, ilon)
            wetveg_moss(ip)  = moss_g(ilat, ilon)
         ENDIF
      ENDDO
      deallocate (lat_g, lon_g, forest_g, cap_g, ombro_g, bog_g)
      IF (allocated(emerg_g)) deallocate (emerg_g, dome_g, moss_g)
      ! A tower sits in one wetland, not in the cell's mix: single-point runs
      ! give its own make-up through these two keys (< 0: keep the file).
      IF (DEF_METHANE%wetland_forest_share_site >= 0._r8) &
         wetveg_forest(:) = min(DEF_METHANE%wetland_forest_share_site, 1._r8)
      IF (DEF_METHANE%wetland_lai_cap_site > 0._r8) wetveg_laicap(:) = DEF_METHANE%wetland_lai_cap_site
      IF (DEF_METHANE%wetland_ombro_share_site >= 0._r8) &
         wetveg_ombro(:) = min(DEF_METHANE%wetland_ombro_share_site, 1._r8)
      IF (DEF_METHANE%wetland_bog_share_site >= 0._r8) &
         wetveg_bog(:) = min(DEF_METHANE%wetland_bog_share_site, 1._r8)
      IF (DEF_METHANE%wetland_dim_class) THEN
         IF (DEF_METHANE%wetland_share_emerg_site >= 0._r8) &
            wetveg_emerg(:) = min(DEF_METHANE%wetland_share_emerg_site, 1._r8)
         IF (DEF_METHANE%wetland_share_dome_site >= 0._r8) &
            wetveg_dome(:) = min(DEF_METHANE%wetland_share_dome_site, 1._r8)
         IF (DEF_METHANE%wetland_share_moss_site >= 0._r8) &
            wetveg_moss(:) = min(DEF_METHANE%wetland_share_moss_site, 1._r8)
         ! water source D2: the tower's open-bog share s_b, and its
         ! rain-fed share s_r, of which the permafrost plateau part s_r - s_b
         ! takes no lateral inflow; s_b is bounded by s_r.
         IF (DEF_METHANE%wetland_share_bog_site >= 0._r8) &
            wetveg_bog(:) = min(DEF_METHANE%wetland_share_bog_site, 1._r8)
         IF (DEF_METHANE%wetland_share_rain_site >= 0._r8) THEN
            srain = min(DEF_METHANE%wetland_share_rain_site, 1._r8)
            wetveg_bog(:)   = min(wetveg_bog(:), srain)
            wetveg_ombro(:) = srain - wetveg_bog(:)
         ENDIF
      ENDIF
      ! canopy height of the forested share of tropical tiles
      IF (DEF_METHANE%wetland_forest_htop_trop > 0._r8) THEN
         CALL wetveg_set_forest_htop (patchlatr_in, numpatch)
         IF (p_is_master) write(*,'(A,F7.2)') &
            ' I-16 forested share of tropical wetland tiles at canopy height [m]: ', &
            DEF_METHANE%wetland_forest_htop_trop
      ENDIF
      wetveg_active = .true.
      IF (p_is_master) write(*,'(A,A)') ' C-13 wetland vegetation read from ', trim(file_veg)
      IF (p_is_master .and. DEF_METHANE%wetland_dim_class) write(*,'(A,5F7.3)') &
         ' C-88 wetland dimension shares on; site emergent, dome, moss, rain-fed, bog (< 0: from the file): ', &
         DEF_METHANE%wetland_share_emerg_site, DEF_METHANE%wetland_share_dome_site, &
         DEF_METHANE%wetland_share_moss_site, DEF_METHANE%wetland_share_rain_site, &
         DEF_METHANE%wetland_share_bog_site

   END SUBROUTINE read_methane_wetveg

   SUBROUTINE wetveg_set_forest_htop (patchlatr_in, numpatch)
      ! Canopy top height of the tropical wetland tiles: the forested
      ! share f stands at wetland_forest_htop_trop and the rest at the land
      ! class height htop0, averaged by share as the PFT heights of a soil
      ! patch (HTOP_readin). Tropical is within 23.5 degrees of the equator,
      ! as in get_biome_f_methane. The height is rebuilt from the land-class
      ! table, so a second call (ch4_reactive_init re-entered) gives the same
      ! value. hbot is not read on the canopy path of a wetland tile
      ! (LeafTemperature) and is left as it is.
      USE MOD_Vars_TimeInvariants, only: htop, patchtype, patchclass
      USE MOD_Const_LC, only: htop0
      USE MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE
      real(r8), intent(in) :: patchlatr_in(:)   ! radians
      integer,  intent(in) :: numpatch
      integer  :: ip
      real(r8) :: f

      IF (.not. p_is_worker) RETURN
      IF (.not. allocated(htop) .or. .not. allocated(patchtype) .or. .not. allocated(patchclass)) RETURN
      DO ip = 1, min(numpatch, size(htop), size(wetveg_forest), size(patchlatr_in))
         IF (patchtype(ip) /= 2) CYCLE
         IF (abs(patchlatr_in(ip)) * 180._r8 / PI > 23.5_r8) CYCLE
         f = min(max(wetveg_forest(ip), 0._r8), 1._r8)
         htop(ip) = f * DEF_METHANE%wetland_forest_htop_trop + (1._r8 - f) * htop0(patchclass(ip))
      ENDDO
   END SUBROUTINE wetveg_set_forest_htop

   SUBROUTINE wetveg_cap_lai ()
      ! Called after each LAI read: the non-forested share of a wetland patch
      ! keeps the remote-sensing LAI only up to its measured peak. With
      ! wetland_lai_shape the whole year is scaled by cap / annual peak
      ! instead, so the peak is the measured one and green-up and senescence
      ! keep the remote-sensing timing; a patch whose peak is below the cap
      ! is left as read, as under the clip.
      USE MOD_Vars_TimeVariables, only: tlai
      USE MOD_Vars_TimeInvariants, only: patchtype
      USE MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE
      integer :: i
      logical :: use_shape
      IF (.not. (DEF_METHANE%wetland_veg_glwd .and. wetveg_active)) RETURN
      IF (.not. p_is_worker) RETURN
      IF (.not. allocated(tlai) .or. .not. allocated(patchtype)) RETURN
      use_shape = DEF_METHANE%wetland_lai_shape .and. allocated(wetveg_laipeak)
      DO i = 1, min(size(tlai), size(wetveg_forest))
         IF (patchtype(i) /= 2) CYCLE
         IF (use_shape) THEN
            IF (i <= size(wetveg_laipeak)) THEN
               IF (wetveg_laipeak(i) > wetveg_laicap(i)) THEN
                  tlai(i) = wetveg_forest(i) * tlai(i) + (1._r8 - wetveg_forest(i)) &
                     * tlai(i) * wetveg_laicap(i) / wetveg_laipeak(i)
                  CYCLE
               ENDIF
            ENDIF
         ENDIF
         tlai(i) = wetveg_forest(i) * tlai(i) + (1._r8 - wetveg_forest(i)) * min(tlai(i), wetveg_laicap(i))
      ENDDO
   END SUBROUTINE wetveg_cap_lai

   real(r8) FUNCTION wetveg_herb_lai_ratio (ipatch)
      ! LAI of the non-forested share over the tile LAI left by
      ! wetveg_cap_lai (<= 1), for wetland_forest_input_herb. With the tile
      ! assimilation split between the two shares by LAI, rescaling the
      ! forested share to the non-forested LAI leaves the tile's assimilation
      ! times this ratio. The clip gives min(1, cap / LAI) for any forested
      ! share f; the shape scaling, with c = cap / peak, gives c / (f + (1 - f) c).
      ! Until the first LAI read of a run the peak is unknown and the clip
      ! relation stands in. 1 without wetland_veg_glwd or without a cap in the cell.
      USE MOD_Vars_TimeVariables, only: tlai
      USE MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE
      integer, intent(in) :: ipatch
      real(r8) :: f, c

      wetveg_herb_lai_ratio = 1._r8
      IF (.not. (DEF_METHANE%wetland_veg_glwd .and. wetveg_active)) RETURN
      IF (.not. allocated(tlai) .or. .not. allocated(wetveg_forest)) RETURN
      IF (ipatch < 1 .or. ipatch > min(size(tlai), size(wetveg_forest))) RETURN
      IF (.not. ieee_is_finite(tlai(ipatch)) .or. tlai(ipatch) <= 0._r8) RETURN
      f = min(max(wetveg_forest(ipatch), 0._r8), 1._r8)
      IF (DEF_METHANE%wetland_lai_shape .and. allocated(wetveg_laipeak)) THEN
         IF (ipatch <= size(wetveg_laipeak)) THEN
            IF (wetveg_laipeak(ipatch) > wetveg_laicap(ipatch)) THEN
               c = wetveg_laicap(ipatch) / wetveg_laipeak(ipatch)
               wetveg_herb_lai_ratio = c / (f + (1._r8 - f) * c)
               RETURN
            ENDIF
         ENDIF
      ENDIF
      IF (tlai(ipatch) > wetveg_laicap(ipatch)) &
         wetveg_herb_lai_ratio = wetveg_laicap(ipatch) / tlai(ipatch)
   END FUNCTION wetveg_herb_lai_ratio

   real(r8) FUNCTION wetveg_forest_input_ratio (ipatch)
      ! Plant input of a wetland tile under wetland_forest_input_herb over its
      ! assimilation-based input: the herb layer ratio r for both shares.
      ! With wetland_dim_class the forested share f keeps at least its
      ! moss layer where the ground is moss-carpeted peat,
      !   (1 - f) r + f max(r, wetland_moss_input_frac s_m),
      ! with the tile's moss share s_m standing for the moss share inside the
      ! forested share.
      USE MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE
      integer, intent(in) :: ipatch
      real(r8) :: r, f, m

      r = wetveg_herb_lai_ratio(ipatch)
      wetveg_forest_input_ratio = r
      IF (.not. DEF_METHANE%wetland_dim_class) RETURN
      IF (.not. allocated(wetveg_moss) .or. .not. allocated(wetveg_forest)) RETURN
      IF (ipatch < 1 .or. ipatch > min(size(wetveg_moss), size(wetveg_forest))) RETURN
      f = min(max(wetveg_forest(ipatch), 0._r8), 1._r8)
      m = DEF_METHANE%wetland_moss_input_frac * min(max(wetveg_moss(ipatch), 0._r8), 1._r8)
      wetveg_forest_input_ratio = r + f * max(m - r, 0._r8)
   END FUNCTION wetveg_forest_input_ratio

   logical FUNCTION wetveg_lai_peak_needed (year)
      ! True when wetland_lai_shape needs the annual peak of the LAI year
      ! about to be read: the first LAI read of a run, or a new LAI year.
      ! The value is the same on every process, so the extra LAI reads of
      ! the peak pass stay collective.
      USE MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE
      integer, intent(in) :: year
      wetveg_lai_peak_needed = DEF_METHANE%wetland_veg_glwd .and. DEF_METHANE%wetland_lai_shape &
         .and. wetveg_active .and. year /= wetveg_peak_year
   END FUNCTION wetveg_lai_peak_needed

   SUBROUTINE wetveg_track_lai_peak (first, year)
      ! Called after each LAI read of the peak pass (every month, or every
      ! 8-day composite, of one LAI year): running maximum of the wetland
      ! patches' remote-sensing LAI before any capping.
      USE MOD_Vars_TimeVariables, only: tlai
      USE MOD_Vars_TimeInvariants, only: patchtype
      logical, intent(in) :: first
      integer, intent(in) :: year
      integer :: i
      wetveg_peak_year = year
      IF (.not. p_is_worker) RETURN
      IF (.not. allocated(tlai) .or. .not. allocated(patchtype) .or. .not. allocated(wetveg_forest)) RETURN
      IF (.not. allocated(wetveg_laipeak)) allocate (wetveg_laipeak(size(wetveg_forest)))
      IF (first) wetveg_laipeak(:) = 0._r8
      DO i = 1, min(size(tlai), size(wetveg_laipeak))
         IF (patchtype(i) /= 2) CYCLE
         IF (ieee_is_finite(tlai(i))) wetveg_laipeak(i) = max(wetveg_laipeak(i), tlai(i))
      ENDDO
   END SUBROUTINE wetveg_track_lai_peak

   real(r8) FUNCTION wetveg_bg_frac (ipatch)
      ! Belowground share of the wetland plant input, weighted by the
      ! forested share; wetland_bg_frac alone without wetland_veg_glwd.
      USE MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE
      integer, intent(in) :: ipatch
      wetveg_bg_frac = DEF_METHANE%wetland_bg_frac
      IF (.not. (DEF_METHANE%wetland_veg_glwd .and. wetveg_active)) RETURN
      IF (ipatch < 1 .or. ipatch > size(wetveg_forest)) RETURN
      wetveg_bg_frac = wetveg_forest(ipatch) * DEF_METHANE%wetland_bg_frac_forest &
         + (1._r8 - wetveg_forest(ipatch)) * DEF_METHANE%wetland_bg_frac
   END FUNCTION wetveg_bg_frac

   SUBROUTINE wetland_root_profile (ipatch, prof)
      ! Root profile of a wetland tile, as layer shares summing to one: the
      ! land-class profile, or with wetland_root_efold > 0
      ! exp(-z/wetland_root_efold) integrated over each layer.
      USE MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE
      USE MOD_Vars_Global, only: nl_soil, zi_soi
      USE MOD_Vars_TimeInvariants, only: patchclass
      USE MOD_Const_LC, only: rootfr
      integer,  intent(in)  :: ipatch
      real(r8), intent(out) :: prof(1:nl_soil)
      integer  :: j
      real(r8) :: ztop, s

      IF (DEF_METHANE%wetland_root_efold > 0._r8) THEN
         ztop = 0._r8
         DO j = 1, nl_soil
            prof(j) = exp(-ztop / DEF_METHANE%wetland_root_efold) &
                    - exp(-zi_soi(j) / DEF_METHANE%wetland_root_efold)
            ztop = zi_soi(j)
         ENDDO
      ELSE
         prof(:) = max(rootfr(1:nl_soil, patchclass(ipatch)), 0._r8)
      ENDIF
      s = sum(prof)
      IF (s > 0._r8) THEN
         prof(:) = prof(:) / s
      ELSE
         prof(:) = 0._r8
         prof(1) = 1._r8
      ENDIF
   END SUBROUTINE wetland_root_profile

   SUBROUTINE set_wetland_veg_glwd (ipatch, lai_in, lai_out, annsum_npp_out, agnpp_out, bgnpp_out, rootfr_out)
      ! Replaces get_wetland_veg_proxy: NPP from the tile's own
      ! assimilation times wetland_npp_frac (as the wetland_plant_input), split by the
      ! forest-weighted belowground share; the land class's root profile (as
      ! the wetland_plant_input); CLM4Me grass aerenchyma defaults, the forested share
      ! at nongrassporosratio of the grass porosity (Riley et al. 2011), or
      ! with woody_conduit_area on its own conduit and the grass
      ! porosity for the tillers.
      ! annsum_npp_out is 0 so the aerenchyma uses the annual means of agnpp
      ! and bgnpp accumulated by methane_annualupdate.
      USE MOD_Vars_Global, only: nl_soil
      USE MOD_Vars_1DFluxes, only: assim
      USE MOD_Vars_TimeInvariants, only: patchclass
      USE MOD_Const_LC, only: rootfr
      USE MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE
      USE MOD_Tracer_Reactive_Methane_VegOverride, only: wetland_aere_poros, wetland_aere_radius, &
         wetland_aere_tillerC, wetland_aere_scale, wetland_aere_active
      integer,  intent(in)  :: ipatch
      real(r8), intent(in)  :: lai_in
      real(r8), intent(out) :: lai_out, annsum_npp_out, agnpp_out, bgnpp_out
      real(r8), intent(out) :: rootfr_out(1:nl_soil)
      real(r8) :: npp, bg, f, s

      f = 0._r8
      IF (ipatch >= 1 .and. ipatch <= size(wetveg_forest)) f = wetveg_forest(ipatch)
      lai_out = max(lai_in, 0._r8)
      npp = 0._r8
      IF (allocated(assim)) THEN
         IF (ieee_is_finite(assim(ipatch)) .and. assim(ipatch) > 0._r8 .and. assim(ipatch) < 1.e-2_r8) &
            npp = assim(ipatch) * 12.011_r8 * DEF_METHANE%wetland_npp_frac       ! [gC m-2 s-1]
      ENDIF
      bg = wetveg_bg_frac(ipatch)
      agnpp_out = (1._r8 - bg) * npp
      bgnpp_out = bg * npp
      annsum_npp_out = 0._r8

      CALL wetland_root_profile (ipatch, rootfr_out)

      IF (allocated(wetland_aere_active) .and. ipatch >= 1 .and. ipatch <= size(wetland_aere_active)) THEN
         IF (DEF_METHANE%woody_conduit_area >= 0._r8) THEN
            ! the forested share has a conduit of its own
            ! (woody_conduit_area, methane_aere), so the tillers of the
            ! non-forested share keep the grass porosity
            wetland_aere_poros  (ipatch) = DEF_METHANE%poros_tiller
         ELSE
         wetland_aere_poros  (ipatch) = DEF_METHANE%poros_tiller * &
            (1._r8 - f + f * DEF_METHANE%nongrassporosratio)
         ENDIF
         wetland_aere_radius (ipatch) = DEF_METHANE%aere_radius
         wetland_aere_tillerC(ipatch) = DEF_METHANE%tiller_C
         wetland_aere_scale  (ipatch) = DEF_METHANE%scale_factor_aere
         wetland_aere_active (ipatch) = .true.
      ENDIF
   END SUBROUTINE set_wetland_veg_glwd

END MODULE MOD_Tracer_Reactive_Methane_WetlandVeg
#endif
