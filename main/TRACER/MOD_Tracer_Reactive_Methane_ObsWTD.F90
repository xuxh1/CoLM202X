#include <define.h>

#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Reactive_Methane_ObsWTD
!=======================================================================
! Tower-observed water table series driving the methane column
! (PLUMBER-style observed-driver mode for site evaluation).
!
! The class-median prescribed table is a year-round constant, so the
! seasonal flux shape carries only the temperature signal; at towers
! that measured WTD (38 of the FLUXNET-CH4 sites) the observed series
! restores the water-side half of the seasonality and the site's own
! magnitude.  Input is the pre-digested daily text written by the
! harness (data/FLUXNET-CH4/WTDdaily/<site>_WTD.txt: year doy metres-
! below-surface, obs positive-above already converted and clamped);
! date arithmetic stays in python where it is safe.
!
! Priority at the consumer: observed series (when a record covers the
! current day) > per-WFT table > site-class table > scalar.
!=======================================================================
   USE MOD_Precision
   IMPLICIT NONE
   SAVE
   PRIVATE

   logical, public :: obs_wtd_loaded = .false.
   integer,  allocatable :: owtd_key(:)   ! year*1000 + doy, ascending
   real(r8), allocatable :: owtd_val(:)   ! m below surface, >= 0

   PUBLIC :: read_obs_wtd, get_obs_wtd, free_obs_wtd

CONTAINS

   SUBROUTINE read_obs_wtd(fname)
      character(len=*), intent(in) :: fname
      integer :: u, ios, n, y, doy
      real(r8) :: w
      character(len=256) :: line

      CALL free_obs_wtd()
      OPEN(newunit=u, file=trim(fname), status='old', action='read', iostat=ios)
      IF (ios /= 0) THEN
         write(6,*) 'ObsWTD: no series at ', trim(fname), &
            ' -- prescribed/class water table applies.'
         RETURN
      ENDIF
      n = 0
      DO
         READ(u, '(A)', iostat=ios) line
         IF (ios /= 0) EXIT
         IF (line(1:1) == '#') CYCLE
         n = n + 1
      ENDDO
      IF (n < 30) THEN
         ! under a month of records cannot drive a run
         write(6,*) 'ObsWTD: only ', n, ' records in ', trim(fname), '; ignored.'
         CLOSE(u)
         RETURN
      ENDIF
      allocate(owtd_key(n), owtd_val(n))
      REWIND(u)
      n = 0
      DO
         READ(u, '(A)', iostat=ios) line
         IF (ios /= 0) EXIT
         IF (line(1:1) == '#') CYCLE
         READ(line, *, iostat=ios) y, doy, w
         IF (ios /= 0) CYCLE
         n = n + 1
         owtd_key(n) = y * 1000 + doy
         owtd_val(n) = max(0._r8, w)
      ENDDO
      CLOSE(u)
      obs_wtd_loaded = .true.
      write(6,*) 'ObsWTD: loaded ', n, ' daily records from ', trim(fname)
   END SUBROUTINE read_obs_wtd

   SUBROUTINE get_obs_wtd(year, doy, wtd_out, ok)
      ! Nearest record within 3 days of (year, doy); outside coverage the
      ! caller falls back to the prescribed chain.
      integer,  intent(in)  :: year, doy
      real(r8), intent(out) :: wtd_out
      logical,  intent(out) :: ok
      integer :: key, lo, hi, mid, best, d1, d2

      ok = .false.
      wtd_out = -1._r8
      IF (.not. obs_wtd_loaded) RETURN
      key = year * 1000 + doy
      lo = 1
      hi = size(owtd_key)
      DO WHILE (lo < hi)
         mid = (lo + hi) / 2
         IF (owtd_key(mid) < key) THEN
            lo = mid + 1
         ELSE
            hi = mid
         ENDIF
      ENDDO
      best = lo
      IF (best > 1) THEN
         d1 = abs(owtd_key(best) - key)
         d2 = abs(owtd_key(best - 1) - key)
         IF (d2 < d1) best = best - 1
      ENDIF
      ! key distance within a year maps 1:1 to day distance; the year*1000
      ! stride keeps cross-year neighbours far apart, which is fine -- a gap
      ! spanning New Year loses at most a few edge days to the tolerance.
      IF (abs(owtd_key(best) - key) <= 3) THEN
         wtd_out = owtd_val(best)
         ok = .true.
      ENDIF
   END SUBROUTINE get_obs_wtd

   SUBROUTINE free_obs_wtd()
      IF (allocated(owtd_key)) deallocate(owtd_key)
      IF (allocated(owtd_val)) deallocate(owtd_val)
      obs_wtd_loaded = .false.
   END SUBROUTINE free_obs_wtd

END MODULE MOD_Tracer_Reactive_Methane_ObsWTD
#endif
