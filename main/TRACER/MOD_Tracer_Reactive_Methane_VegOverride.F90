#include <define.h>

#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Reactive_Methane_VegOverride
!=======================================================================
! Per-patch aerenchyma parameter overrides for wetland (patchtype==2).
!
! Written by get_wetland_veg_proxy in MOD_Tracer_Reactive_Methane_BgcLink based
! on 5-zone climate classification (tropical reed/papyrus, tropical
! swamp, temperate marsh, boreal fen, Sphagnum bog).
!
! Read by methane_aere / SiteOxAere in MOD_Tracer_Reactive_Methane_Physics via
! getters that fall through to DEF_METHANE%* defaults when the
! per-patch override is inactive.
!
! Module exists separately from BgcLink so Physics can read these
! arrays without a Physics->BgcLink dependency (BgcLink compiles
! after Physics in the Makefile).
!
!=======================================================================

   USE MOD_Precision
   IMPLICIT NONE
   SAVE
   PRIVATE

   real(r8), allocatable, public :: wetland_aere_poros  (:)   ! tiller porosity (-)
   real(r8), allocatable, public :: wetland_aere_radius (:)   ! tiller radius (m)
   real(r8), allocatable, public :: wetland_aere_tillerC(:)   ! gC per tiller
   real(r8), allocatable, public :: wetland_aere_scale  (:)   ! scale_factor_aere multiplier (-)
   logical,  allocatable, public :: wetland_aere_active (:)   ! .true. = use override for patch

   ! Wetland functional type of each patch, the axis for per-class parameter
   ! sets calibrated at the FLUXNET-CH4 towers:
   !   0 = unset (every get_wft_param falls through to the global default)
   !   1 = permafrost peatland   2 = peatland   3 = mineral wetland
   ! Single point: derived from SITE_wetland_class.  Global: to be filled by a
   ! static class map (permafrost cut: MAAT/permafrost map; peat cut:
   ! Peat-ML / GLWD peat classes / OM_density); until that reader exists the
   ! array stays 0 and behaviour is bit-identical to the global-scalar model.
   integer,  allocatable, public :: wetland_wft_class   (:)
   ! Topsoil pH sampled from the environment map for wetland patches
   ! (SoilGrids layer in the wft class file); <=0 = unset.  Priority at the
   ! consumer: SITE_ph > spatial pH vector > this map > ph_fallback.
   real(r8), allocatable, public :: wetland_ph_map      (:)
   ! Salinity assigned from the environment map's tidal-marsh fraction
   ! (Worthington layer); <0 = unset.  Priority: SITE_salinity > this map.
   real(r8), allocatable, public :: wetland_salinity_map(:)

   PUBLIC :: allocate_wetland_aere_overrides
   PUBLIC :: deallocate_wetland_aere_overrides
   PUBLIC :: get_aere_poros, get_aere_radius, get_aere_tillerC, get_aere_scale
   PUBLIC :: get_wft_param
   PUBLIC :: load_wft_class_map

CONTAINS

   SUBROUTINE allocate_wetland_aere_overrides(numpatch)
      integer, intent(in) :: numpatch
      IF (allocated(wetland_aere_poros)) RETURN
      ! Keep zero-length arrays allocated on ranks with numpatch==0; downstream
      ! getters guard by bounds, and collective setup expects allocation state
      ! to be consistent across MPI ranks.
      allocate(wetland_aere_poros  (numpatch))
      allocate(wetland_aere_radius (numpatch))
      allocate(wetland_aere_tillerC(numpatch))
      allocate(wetland_aere_scale  (numpatch))
      allocate(wetland_aere_active (numpatch))
      allocate(wetland_wft_class   (numpatch))
      allocate(wetland_ph_map      (numpatch))
      allocate(wetland_salinity_map(numpatch))
      wetland_aere_poros  (:) = 0._r8
      wetland_aere_radius (:) = 0._r8
      wetland_aere_tillerC(:) = 0._r8
      wetland_aere_scale  (:) = 0._r8
      wetland_aere_active (:) = .false.
      wetland_wft_class   (:) = 0
      wetland_ph_map      (:) = -1._r8
      wetland_salinity_map(:) = -1._r8
   END SUBROUTINE allocate_wetland_aere_overrides

   SUBROUTINE deallocate_wetland_aere_overrides()
      IF (allocated(wetland_aere_poros))   deallocate(wetland_aere_poros)
      IF (allocated(wetland_aere_radius))  deallocate(wetland_aere_radius)
      IF (allocated(wetland_aere_tillerC)) deallocate(wetland_aere_tillerC)
      IF (allocated(wetland_aere_scale))   deallocate(wetland_aere_scale)
      IF (allocated(wetland_aere_active))  deallocate(wetland_aere_active)
      IF (allocated(wetland_wft_class))    deallocate(wetland_wft_class)
      IF (allocated(wetland_ph_map))       deallocate(wetland_ph_map)
      IF (allocated(wetland_salinity_map)) deallocate(wetland_salinity_map)
   END SUBROUTINE deallocate_wetland_aere_overrides

   real(r8) FUNCTION get_aere_poros(ipatch, default_val)
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: default_val
      get_aere_poros = default_val
      IF (.not. allocated(wetland_aere_active)) RETURN
      IF (ipatch < 1 .or. ipatch > size(wetland_aere_active)) RETURN
      IF (wetland_aere_active(ipatch)) get_aere_poros = wetland_aere_poros(ipatch)
   END FUNCTION get_aere_poros

   real(r8) FUNCTION get_aere_radius(ipatch, default_val)
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: default_val
      get_aere_radius = default_val
      IF (.not. allocated(wetland_aere_active)) RETURN
      IF (ipatch < 1 .or. ipatch > size(wetland_aere_active)) RETURN
      IF (wetland_aere_active(ipatch)) get_aere_radius = wetland_aere_radius(ipatch)
   END FUNCTION get_aere_radius

   real(r8) FUNCTION get_aere_tillerC(ipatch, default_val)
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: default_val
      get_aere_tillerC = default_val
      IF (.not. allocated(wetland_aere_active)) RETURN
      IF (ipatch < 1 .or. ipatch > size(wetland_aere_active)) RETURN
      IF (wetland_aere_active(ipatch)) get_aere_tillerC = wetland_aere_tillerC(ipatch)
   END FUNCTION get_aere_tillerC

   real(r8) FUNCTION get_aere_scale(ipatch, default_val)
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: default_val
      get_aere_scale = default_val
      IF (.not. allocated(wetland_aere_active)) RETURN
      IF (ipatch < 1 .or. ipatch > size(wetland_aere_active)) RETURN
      IF (wetland_aere_active(ipatch)) get_aere_scale = wetland_aere_scale(ipatch)
   END FUNCTION get_aere_scale

   SUBROUTINE load_wft_class_map(file_class, patchlatr_in, patchlonr_in, patchtype_in, numpatch, &
         tidal_frac_min, tidal_salinity_psu)
      ! Populate wetland_wft_class from the static global class map for grid
      ! runs.  The map is small (5-arcmin byte grid, ~8 MB) so every rank
      ! reads its own copy from /share -- no MPI choreography needed.
      !
      ! Precedence: a class already set (single-point path writes it from
      ! SITE_wetland_class before this runs) always wins; the map only fills
      ! patches still at 0, and only wetland patches (patchtype 2).  A missing
      ! or unreadable file logs a warning and leaves everything unset, which
      ! downstream means "fall through to the global scalars" -- the safe
      ! failure mode.
      USE netcdf
      character(len=*), intent(in) :: file_class
      real(r8), intent(in) :: patchlatr_in(:), patchlonr_in(:)   ! radians
      integer,  intent(in) :: patchtype_in(:)
      integer,  intent(in) :: numpatch
      real(r8), intent(in) :: tidal_frac_min      ! tidal fraction above which the cell counts as tidal
      real(r8), intent(in) :: tidal_salinity_psu  ! salinity assigned to tidal wetland patches

      real(r8), parameter :: PI = 3.14159265358979323846_r8
      integer :: ncid, vid, dimid, nlat, nlon, ierr
      integer :: ipatch, ilat, ilon, nset
      real(r8), allocatable :: mlat(:), mlon(:)
      integer(1), allocatable :: mclass(:,:)
      real(r8), allocatable :: mph(:,:), mtd(:,:)
      logical :: has_ph, has_td
      real(r8) :: lat_deg, lon_deg

      IF (len_trim(file_class) == 0 .or. trim(file_class) == 'null') RETURN
      IF (.not. allocated(wetland_wft_class)) RETURN
      IF (numpatch < 1) RETURN

      ierr = nf90_open(trim(file_class), NF90_NOWRITE, ncid)
      IF (ierr /= NF90_NOERR) THEN
         write(6,*) 'WARNING load_wft_class_map: cannot open ', trim(file_class), &
            '; wetland WFT classes stay unset (global scalars apply).'
         RETURN
      ENDIF
      ierr = nf90_inq_dimid(ncid, 'lat', dimid)
      IF (ierr == NF90_NOERR) ierr = nf90_inquire_dimension(ncid, dimid, len=nlat)
      IF (ierr == NF90_NOERR) ierr = nf90_inq_dimid(ncid, 'lon', dimid)
      IF (ierr == NF90_NOERR) ierr = nf90_inquire_dimension(ncid, dimid, len=nlon)
      IF (ierr /= NF90_NOERR) THEN
         write(6,*) 'WARNING load_wft_class_map: bad lat/lon dims in ', trim(file_class)
         ierr = nf90_close(ncid); RETURN
      ENDIF
      allocate(mlat(nlat), mlon(nlon), mclass(nlon, nlat))
      ierr = nf90_inq_varid(ncid, 'lat', vid)
      IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, mlat)
      IF (ierr == NF90_NOERR) ierr = nf90_inq_varid(ncid, 'lon', vid)
      IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, mlon)
      IF (ierr == NF90_NOERR) ierr = nf90_inq_varid(ncid, 'wft_class', vid)
      ! netCDF stores (lat, lon); Fortran reads it transposed as (lon, lat)
      IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, mclass)
      has_ph = .false.
      IF (ierr == NF90_NOERR) THEN
         IF (nf90_inq_varid(ncid, 'ph', vid) == NF90_NOERR) THEN
            allocate(mph(nlon, nlat))
            IF (nf90_get_var(ncid, vid, mph) == NF90_NOERR) has_ph = .true.
         ENDIF
      ENDIF
      has_td = .false.
      IF (ierr == NF90_NOERR) THEN
         IF (nf90_inq_varid(ncid, 'tidal_frac', vid) == NF90_NOERR) THEN
            allocate(mtd(nlon, nlat))
            IF (nf90_get_var(ncid, vid, mtd) == NF90_NOERR) has_td = .true.
         ENDIF
      ENDIF
      IF (ierr /= NF90_NOERR) THEN
         write(6,*) 'WARNING load_wft_class_map: read failed on ', trim(file_class)
         deallocate(mlat, mlon, mclass)
         ierr = nf90_close(ncid); RETURN
      ENDIF
      ierr = nf90_close(ncid)

      nset = 0
      DO ipatch = 1, min(numpatch, size(wetland_wft_class))
         IF (patchtype_in(ipatch) /= 2) CYCLE
         IF (wetland_wft_class(ipatch) /= 0) CYCLE   ! site-assigned class wins
         lat_deg = patchlatr_in(ipatch) * 180._r8 / PI
         lon_deg = patchlonr_in(ipatch) * 180._r8 / PI
         IF (lon_deg > 180._r8) lon_deg = lon_deg - 360._r8
         ilat = minloc(abs(mlat - lat_deg), dim=1)
         ilon = minloc(abs(mlon - lon_deg), dim=1)
         IF (mclass(ilon, ilat) >= 1 .and. mclass(ilon, ilat) <= 3) THEN
            wetland_wft_class(ipatch) = int(mclass(ilon, ilat))
            nset = nset + 1
         ENDIF
         IF (has_ph .and. allocated(wetland_ph_map)) THEN
            IF (mph(ilon, ilat) > 0._r8 .and. mph(ilon, ilat) < 14._r8) &
               wetland_ph_map(ipatch) = mph(ilon, ilat)
         ENDIF
         IF (has_td .and. allocated(wetland_salinity_map)) THEN
            IF (mtd(ilon, ilat) >= tidal_frac_min) &
               wetland_salinity_map(ipatch) = tidal_salinity_psu
         ENDIF
      ENDDO
      write(6,*) 'load_wft_class_map: assigned WFT class to ', nset, ' wetland patches from ', &
         trim(file_class)
      deallocate(mlat, mlon, mclass)
      IF (allocated(mph)) deallocate(mph)
      IF (allocated(mtd)) deallocate(mtd)
   END SUBROUTINE load_wft_class_map

   real(r8) FUNCTION get_wft_param(ipatch, default_val, per_class)
      ! Per-WFT parameter lookup: return per_class(class of ipatch) when that
      ! class is set (1..3) and the value is not the -999 "unset" sentinel;
      ! otherwise the global default.  Zero is a legitimate calibrated value
      ! (e.g. aereoxid), so the sentinel test is `< 0`, not `<= 0`.
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: default_val
      real(r8), intent(in) :: per_class(3)
      integer :: c
      get_wft_param = default_val
      IF (.not. allocated(wetland_wft_class)) RETURN
      IF (ipatch < 1 .or. ipatch > size(wetland_wft_class)) RETURN
      c = wetland_wft_class(ipatch)
      IF (c < 1 .or. c > 3) RETURN
      IF (per_class(c) >= 0._r8) get_wft_param = per_class(c)
   END FUNCTION get_wft_param

END MODULE MOD_Tracer_Reactive_Methane_VegOverride
#endif
