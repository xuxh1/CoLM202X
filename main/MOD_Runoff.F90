#include <define.h>

MODULE MOD_Runoff

!-----------------------------------------------------------------------
   USE MOD_Precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: SurfaceRunoff_TOPMOD
   PUBLIC :: SubsurfaceRunoff_TOPMOD
   PUBLIC :: Runoff_XinAnJiang
   PUBLIC :: Runoff_SimpleVIC


!-----------------------------------------------------------------------

CONTAINS

   SUBROUTINE SurfaceRunoff_TOPMOD (nl_soil,wimp,porsl,psi0,hksati,&
                                    fsatmax,fsatdcf,&
                                    z_soisno,dz_soisno,zi_soisno,&
                                    eff_porosity,icefrac,zwt,gwat,&
                                    rsur,rsur_se,rsur_ie,&
                                    topoweti,alp_twi,chi_twi,mu_twi,frcsat,eta_out)

!=======================================================================
!  the original code was provide by Robert E. Dickinson based on
!  following clues: a water table level determination level added
!  including highland and lowland levels and fractional area of wetland
!  (water table above the surface.  Runoff is parametrized from the
!  lowlands in terms of precip incident on wet areas and a base flow,
!  where these are estimated using ideas from TOPMODEL.
!
!  Author : Yongjiu Dai, 07/29/2002, Guoyue Niu, 06/2012
!=======================================================================

   USE MOD_Namelist,        only: DEF_TOPMOD_method
   USE MOD_IncompleteGamma, only: GRATIO
   USE MOD_SPMD_Task
   IMPLICIT NONE

!-------------------------- Dummy Arguments ----------------------------

   integer, intent(in) :: nl_soil   ! number of soil layers
   real(r8), intent(in) :: &
        ! wtfact,                 &! (updated to gridded 'fsatmax' data)
                                   ! fraction of model area with high water table
        wimp,                     &! water impermeable if porosity less than wimp
        porsl(1:nl_soil),         &! saturated volumetric soil water content(porosity)
        psi0(1:nl_soil),          &! saturated soil suction (mm) (NEGATIVE)
        hksati(1:nl_soil),        &! hydraulic conductivity at saturation (mm h2o/s)
        fsatmax,                  &! maximum fraction of saturation area [-]
        fsatdcf,                  &! decay factor in calc of fraction of saturation area [1/m]
        z_soisno(1:nl_soil),      &! layer depth (m)
        dz_soisno(1:nl_soil),     &! layer thickness (m)
        zi_soisno(0:nl_soil),     &! interface level below a "z" level (m)
        eff_porosity(1:nl_soil),  &! effective porosity = porosity - vol_ice
        icefrac(1:nl_soil),       &! ice fraction (-)
        gwat,                     &! net water input from top
        zwt                        ! the depth from ground (soil) surface to water table [m]

   real(r8), intent(out) :: rsur   ! surface runoff (mm h2o/s)
   real(r8), intent(out), optional :: rsur_se! saturation excess surface runoff (mm h2o/s)
   real(r8), intent(out), optional :: rsur_ie! infiltration excess surface runoff (mm h2o/s)

   real(r8), intent(in),  optional :: topoweti
   real(r8), intent(in),  optional :: alp_twi, chi_twi, mu_twi
   real(r8), intent(out), optional :: frcsat
   real(r8), intent(out), optional :: eta_out

!-------------------------- Local Variables ----------------------------

   real(r8), parameter :: vdcf = 2.0

   real(r8) qinmax       ! maximum infiltration capability
   real(r8) fsat         ! fractional area with water table at surface

   real(r8) eta, pgr0, pgr1, qgr, gfun
   integer  niter

   ! updated to gridded 'fsatdcf' (by Shupeng Zhang)
   ! real(r8), parameter :: fff = 0.5   ! runoff decay factor (m-1)

!-----------------------------------------------------------------------

!  fraction of saturated area (updated to gridded 'fsatmax' and 'fsatdcf')
      !fsat = wtfact*min(1.0,exp(-0.5*fff*zwt))
      IF ((DEF_TOPMOD_method == 0) .or. (DEF_TOPMOD_method == 1)) THEN

         fsat = fsatmax * exp(- fsatdcf * vdcf * zwt)

      ELSE

         IF (zwt <= 0.) THEN

            fsat = 1.
            eta  = mu_twi

         ELSE

            eta = topoweti
            niter = 0
            DO WHILE (niter < 20)
               niter = niter + 1
               CALL GRATIO (alp_twi+1, (eta-mu_twi)/chi_twi, pgr1, qgr, 0)
               CALL GRATIO (alp_twi,   (eta-mu_twi)/chi_twi, pgr0, qgr, 0)
               gfun = ((eta-mu_twi)*pgr0 - chi_twi*alp_twi*pgr1)/vdcf - zwt

               IF (abs(gfun) > 1.e-6) THEN
                  eta = mu_twi + (chi_twi * alp_twi * pgr1 + vdcf*zwt) / pgr0
               ELSE
                  EXIT
               ENDIF
            ENDDO

            IF (abs(gfun) > 1.e-6) THEN
               write(*,*) 'Fail to converge in TOPModel: (alp,chi,mu,twi,zwt,gfun) = ', &
                  alp_twi, chi_twi, mu_twi, topoweti, zwt, gfun
            ENDIF

            CALL GRATIO (alp_twi, (eta-mu_twi)/chi_twi, pgr0, qgr, 0)

            fsat = qgr

         ENDIF

         IF (present(eta_out)) THEN
            eta_out = eta
         ENDIF

      ENDIF

      IF (present(frcsat)) THEN
         frcsat = fsat
      ENDIF

! Maximum infiltration capacity
      qinmax = minval(10.**(-6.0*icefrac(1:min(3,nl_soil)))*hksati(1:min(3,nl_soil)))
      IF(eff_porosity(1)<wimp) qinmax = 0.

! Surface runoff
      rsur = fsat*max(0.0,gwat) + (1.-fsat)*max(0.,gwat-qinmax)

      IF (present(rsur_se)) THEN
         rsur_se = fsat*max(0.0,gwat)
      ENDIF

      IF (present(rsur_ie)) THEN
         rsur_ie = (1.-fsat)*max(0.,gwat-qinmax)
      ENDIF

   END SUBROUTINE SurfaceRunoff_TOPMOD

   SUBROUTINE SurfaceRunoff_TOPMOD_CLM (nl_soil,wimp,porsl,psi0,hksati,&
                                 fsatmax,fsatdcf,&
                                 z_soisno,dz_soisno,zi_soisno,&
                                 eff_porosity,icefrac,gwat,zwt,slpratio,deltim,pondmx,&
                                 f_h2osfc,wdsrf,rsur,q_soil,q_excess,q_h2osfc,q_drain_h2osfc,q_h2osfc_surf,rsur_se,rsur_ie,&
                                 topoweti,alp_twi,chi_twi,mu_twi,frcsat,eta_out)

!=======================================================================
!  the original code was provide by Robert E. Dickinson based on
!  following clues: a water table level determination level added
!  including highland and lowland levels and fractional area of wetland
!  (water table above the surface.  Runoff is parametrized from the
!  lowlands in terms of precip incident on wet areas and a base flow,
!  where these are estimated using ideas from TOPMODEL.
!
!  Author : Yongjiu Dai, 07/29/2002, Guoyue Niu, 06/2012
!=======================================================================

   USE MOD_Namelist,        only: DEF_TOPMOD_method
   USE MOD_IncompleteGamma, only: GRATIO
   USE MOD_SPMD_Task
   USE MOD_Const_ch4, only: DEF_CH4_hydrology
   USE MOD_Vars_Global, only : PI
   IMPLICIT NONE

!-------------------------- Dummy Arguments ----------------------------

   integer, intent(in) :: nl_soil   ! number of soil layers
   real(r8), intent(in) :: &
        ! wtfact,                 &! (updated to gridded 'fsatmax' data)
                                   ! fraction of model area with high water table
        wimp,                     &! water impermeable if porosity less than wimp
        porsl(1:nl_soil),         &! saturated volumetric soil water content(porosity)
        psi0(1:nl_soil),          &! saturated soil suction (mm) (NEGATIVE)
        hksati(1:nl_soil),        &! hydraulic conductivity at saturation (mm h2o/s)
        fsatmax,                  &! maximum fraction of saturation area [-]
        fsatdcf,                  &! decay factor in calc of fraction of saturation area [1/m]
        z_soisno(1:nl_soil),      &! layer depth (m)
        dz_soisno(1:nl_soil),     &! layer thickness (m)
        zi_soisno(0:nl_soil),     &! interface level below a "z" level (m)
        eff_porosity(1:nl_soil),  &! effective porosity = porosity - vol_ice
        icefrac(1:nl_soil),       &! ice fraction (-)
        gwat,                     &! net water input from top
        zwt                      ,&! the depth from ground (soil) surface to water table [m]
        slpratio,                 &! the slope ratio
        deltim,                   &
        pondmx

   real(r8), intent(inout) :: &
         f_h2osfc ,& ! fractional inundated area
         wdsrf       ! surface water (mm)

   real(r8), intent(out) :: rsur   ! surface runoff (mm h2o/s)
   real(r8), intent(out) :: q_soil, &
         q_excess, &
         q_h2osfc, &
         q_drain_h2osfc, &
         q_h2osfc_surf

   real(r8), intent(out), optional :: rsur_se! saturation excess surface runoff (mm h2o/s)
   real(r8), intent(out), optional :: rsur_ie! infiltration excess surface runoff (mm h2o/s)

   real(r8), intent(in),  optional :: topoweti
   real(r8), intent(in),  optional :: alp_twi, chi_twi, mu_twi
   real(r8), intent(out), optional :: frcsat
   real(r8), intent(out), optional :: eta_out

!-------------------------- Local Variables ----------------------------

   real(r8), parameter :: vdcf = 2.0

   real(r8) qinmax       ! maximum infiltration capability
   real(r8) fsat         ! fractional area with water table at surface

   real(r8) eta, pgr0, pgr1, qgr, gfun
   integer  niter

   real(r8) :: micro_sigma, sigma, d
   real(r8) :: fd, dfdd
   real(r8) :: pondmin
   real(r8) :: f_connected
   real(r8) :: k_h2osfc
   integer :: l, p

   real(r8) :: q_liq0, q_over, q_in_surface
   real(r8) :: q_soil_pot, q_soil_final
   real(r8) :: sigma_mm, fd_tol, wdsrf_thresh
   ! updated to gridded 'fsatdcf' (by Shupeng Zhang)
   ! real(r8), parameter :: fff = 0.5   ! runoff decay factor (m-1)
   real(r8) :: wdsrf_old, f_h2osfc_old
   real(r8) :: chk_bal, q_in, q_out, q_store
   real(r8) :: rsur_old
!-----------------------------------------------------------------------

!  fraction of saturated area (updated to gridded 'fsatmax' and 'fsatdcf')
      !fsat = wtfact*min(1.0,exp(-0.5*fff*zwt))
      IF ((DEF_TOPMOD_method == 0) .or. (DEF_TOPMOD_method == 1)) THEN

         fsat = fsatmax * exp(- fsatdcf * vdcf * zwt)

      ELSE

         IF (zwt <= 0.) THEN

            fsat = 1.
            eta  = mu_twi

         ELSE

            eta = topoweti
            niter = 0
            DO WHILE (niter < 20)
               niter = niter + 1
               CALL GRATIO (alp_twi+1, (eta-mu_twi)/chi_twi, pgr1, qgr, 0)
               CALL GRATIO (alp_twi,   (eta-mu_twi)/chi_twi, pgr0, qgr, 0)
               gfun = ((eta-mu_twi)*pgr0 - chi_twi*alp_twi*pgr1)/vdcf - zwt

               IF (abs(gfun) > 1.e-6) THEN
                  eta = mu_twi + (chi_twi * alp_twi * pgr1 + vdcf*zwt) / pgr0
               ELSE
                  EXIT
               ENDIF
            ENDDO

            IF (abs(gfun) > 1.e-6) THEN
               write(*,*) 'Fail to converge in TOPModel: (alp,chi,mu,twi,zwt,gfun) = ', &
                  alp_twi, chi_twi, mu_twi, topoweti, zwt, gfun
            ENDIF

            CALL GRATIO (alp_twi, (eta-mu_twi)/chi_twi, pgr0, qgr, 0)

            fsat = qgr

         ENDIF

         IF (present(eta_out)) THEN
            eta_out = eta
         ENDIF

      ENDIF

      IF (present(frcsat)) THEN
         frcsat = fsat
      ENDIF

! Maximum infiltration capacity
      qinmax = minval(10.**(-6.0*icefrac(1:min(3,nl_soil)))*hksati(1:min(3,nl_soil)))
      IF(eff_porosity(1)<wimp) qinmax = 0.
      wdsrf_old   = wdsrf
      f_h2osfc_old = f_h2osfc
      rsur_old    = rsur

      ! Surface Water

      fd_tol = 1.e-10_r8

      !--------------------------------------------------
      ! 0) Define available liquid input (no negative rain input into runoff)
      !--------------------------------------------------
      q_liq0 = max(0._r8, gwat)

      !--------------------------------------------------
      ! 1) Saturation-excess runoff (TOPMODEL)
      !--------------------------------------------------
      q_over = fsat * q_liq0
      q_in_surface = (1._r8 - fsat) * q_liq0

      !--------------------------------------------------
      ! 2) Diagnose f_h2osfc from current wdsrf (OLD state) before flux partition
      !    (Newton solve with guards)
      !--------------------------------------------------
      pondmin = 1.e-8_r8

      ! compute micro_sigma (unit should be meters); keep your formula but clamp
      micro_sigma = (atan(slpratio) + DEF_CH4_hydrology%slopemax**(1._r8/DEF_CH4_hydrology%slopebeta))**DEF_CH4_hydrology%slopebeta
      micro_sigma = max(0._r8, min(DEF_CH4_hydrology%slopemax, micro_sigma))
      sigma_mm = 1.0e3_r8 * micro_sigma   ! m -> mm

      if (sigma_mm > 1.e-3_r8) then   ! 1e-6 m == 1e-3 mm

         d = 0.0_r8   ! d is now in mm

         do p = 1, 4
            fd   = 0.5_r8 * (1.0_r8 + erf(d / (sigma_mm * sqrt(2.0_r8)))) &
               - DEF_CH4_hydrology%pc

            dfdd = exp( -d**2 / (2.0_r8 * sigma_mm**2) ) &
               / ( sigma_mm * sqrt(2.0_r8 * PI) )

            d = d - fd / dfdd
         enddo

         ! W(d) in mm (CLM fill-and-spill formula)
         wdsrf_thresh = 0.5_r8 * d * (1.0_r8 + erf(d / (sigma_mm * sqrt(2.0_r8)))) + &
                        sigma_mm / sqrt(2.0_r8 * PI) * &
                        exp( -d**2 / (2.0_r8 * sigma_mm**2) )

      else
         wdsrf_thresh = 0._r8
      endif

      ! if (wdsrf > pondmin .and. sigma_mm > 1.e-9_r8) then
      !    ! initial guess: d ~ wdsrf (both in mm)
      !    d = max(0._r8, wdsrf)

      !    do l = 1, 20
      !       fd = 0.5_r8*d*(1.0_r8+erf(d/(sigma_mm*sqrt(2.0_r8)))) &
      !          + sigma_mm/sqrt(2.0_r8*PI)*exp(-d**2/(2.0_r8*sigma_mm**2)) &
      !          - wdsrf

      !       dfdd = 0.5_r8*(1.0_r8+erf(d/(sigma_mm*sqrt(2.0_r8))))

      !       if (abs(fd) < fd_tol) exit
      !       if (dfdd < 1.e-12_r8) exit

      !       d = d - fd/dfdd

      !       ! guard against runaway
      !       d = max(-10._r8*sigma_mm, min(10._r8*sigma_mm, d))
      !    enddo

      !    f_h2osfc = 0.5_r8*(1.0_r8+erf(d/(sigma_mm*sqrt(2.0_r8))))
      !    f_h2osfc = min(1._r8, max(0._r8, f_h2osfc))
      ! else
      !    f_h2osfc = 0._r8
      ! endif

      !--------------------------------------------------
      ! 3) Partition incoming water between soil and h2osfc (CLM-consistent)
      !    - q_soil: FINAL soil infiltration flux (excluding h2osfc drainage)
      !    - q_excess: Hortonian excess moved into h2osfc
      !--------------------------------------------------
      q_soil_pot = (1._r8 - f_h2osfc) * q_in_surface

      ! IMPORTANT: threshold must include non-inundated fraction
      q_excess   = max( q_soil_pot - (1._r8 - f_h2osfc)*qinmax, 0._r8 )

      q_h2osfc   = f_h2osfc * q_in_surface + q_excess
      q_soil_final = q_soil_pot - q_excess

      ! return q_soil as FINAL soil-side infiltration (excluding drainage)
      q_soil = max(0._r8, q_soil_final)

      !--------------------------------------------------
      ! 4) Connectivity function (CLM-style)
      !--------------------------------------------------
      if (f_h2osfc <= DEF_CH4_hydrology%pc) then
         f_connected = 0._r8
      else
         f_connected = (f_h2osfc - DEF_CH4_hydrology%pc)**0.14
      endif

      !--------------------------------------------------
      ! 5) h2osfc surface outflow (only when ponding above pondmx)
      !    NOTE: q_h2osfc_surf is a RUNOFF component and must be counted in rsur.
      !--------------------------------------------------
      if (wdsrf > wdsrf_thresh) then
         k_h2osfc = 1.0e-4_r8 * sin(atan(slpratio))
         q_h2osfc_surf = k_h2osfc * f_connected * (wdsrf - wdsrf_thresh)
         q_h2osfc_surf = min(q_h2osfc_surf, (wdsrf - wdsrf_thresh)/deltim)
      else
         q_h2osfc_surf = 0._r8
      endif
      if (q_h2osfc_surf < 1.e-12_r8) q_h2osfc_surf = 0._r8

      !--------------------------------------------------
      ! 6) Update wdsrf by adding h2osfc inflow and subtracting surface outflow
      !--------------------------------------------------
      wdsrf = wdsrf + (q_h2osfc - q_h2osfc_surf) * deltim
      wdsrf = max(0._r8, wdsrf)

      !--------------------------------------------------
      ! 7) h2osfc drainage to soil (adds to soil infiltration later)
      !--------------------------------------------------
      q_drain_h2osfc = min( f_h2osfc*qinmax, wdsrf/deltim )
      q_drain_h2osfc = max(0._r8, q_drain_h2osfc)

      wdsrf = wdsrf - q_drain_h2osfc * deltim
      wdsrf = max(0._r8, wdsrf)

      !--------------------------------------------------
      ! 8) Re-diagnose f_h2osfc from UPDATED wdsrf (NEW state)
      !--------------------------------------------------
      if (wdsrf > pondmin .and. micro_sigma > 1.e-12_r8) then
         sigma_mm = 1.0e3_r8 * micro_sigma
         d = max(0._r8, wdsrf)

         do l = 1, 20
            fd = 0.5_r8*d*(1.0_r8+erf(d/(sigma_mm*sqrt(2.0_r8)))) &
               + sigma_mm/sqrt(2.0_r8*PI)*exp(-d**2/(2.0_r8*sigma_mm**2)) &
               - wdsrf
            dfdd = 0.5_r8*(1.0_r8+erf(d/(sigma_mm*sqrt(2.0_r8))))

            if (abs(fd) < fd_tol) exit
            if (dfdd < 1.e-12_r8) exit

            d = d - fd/dfdd
            d = max(-10._r8*sigma_mm, min(10._r8*sigma_mm, d))
         enddo

         f_h2osfc = 0.5_r8*(1.0_r8+erf(d/(sigma_mm*sqrt(2.0_r8))))
         f_h2osfc = min(1._r8, max(0._r8, f_h2osfc))
      else
         f_h2osfc = 0._r8
      endif

      !--------------------------------------------------
      ! 9) FINAL surface runoff returned to WATER_2014
      !--------------------------------------------------
      rsur = q_over + q_h2osfc_surf

      if (present(rsur_se)) rsur_se = q_over
      if (present(rsur_ie)) rsur_ie = q_h2osfc_surf


      ! q_in    = q_liq0
      ! q_out   = rsur + q_soil + q_drain_h2osfc
      ! q_store = (wdsrf - wdsrf_old) / deltim
      ! chk_bal = q_in - q_out - q_store   ! mm/s, should be ~0
      ! write(6,'(A)') '---------------- SurfaceRunoff_TOPMOD_CLM DEBUG ----------------'
      ! write(6,'(A,1PE12.4)') 'deltim [s]              = ', deltim
      ! write(6,'(A,1PE12.4)') 'slpratio [-]            = ', slpratio
      ! write(6,'(A,1PE12.4)') 'zwt [m]                 = ', zwt

      ! write(6,'(A,1PE12.4)') 'fsat [-]                = ', fsat
      ! write(6,'(A,1PE12.4)') 'qinmax [mm/s]           = ', qinmax
      ! write(6,'(A,1PE12.4)') 'q_liq0 [mm/s]           = ', q_liq0
      ! write(6,'(A,1PE12.4)') 'q_over [mm/s]           = ', q_over
      ! write(6,'(A,1PE12.4)') 'q_in_surface [mm/s]     = ', q_in_surface

      ! write(6,'(A,1PE12.4)') 'micro_sigma [m]         = ', micro_sigma
      ! write(6,'(A,1PE12.4)') 'sigma_mm [mm]           = ', sigma_mm
      ! write(6,'(A,1PE12.4)') 'wdsrf_thresh [mm]       = ', wdsrf_thresh
      ! write(6,'(A,1PE12.4)') 'pondmx [mm]             = ', pondmx
      ! write(6,'(A,1PE12.4)') 'pondmin [mm]            = ', pondmin

      ! write(6,'(A,1PE12.4)') 'f_h2osfc_old [-]        = ', f_h2osfc_old
      ! write(6,'(A,1PE12.4)') 'f_h2osfc_new [-]        = ', f_h2osfc
      ! write(6,'(A,1PE12.4)') 'f_connected [-]         = ', f_connected

      ! write(6,'(A,1PE12.4)') 'q_soil_pot [mm/s]       = ', q_soil_pot
      ! write(6,'(A,1PE12.4)') 'q_excess [mm/s]         = ', q_excess
      ! write(6,'(A,1PE12.4)') 'q_h2osfc [mm/s]         = ', q_h2osfc
      ! write(6,'(A,1PE12.4)') 'q_h2osfc_surf [mm/s]    = ', q_h2osfc_surf
      ! write(6,'(A,1PE12.4)') 'q_drain_h2osfc [mm/s]   = ', q_drain_h2osfc
      ! write(6,'(A,1PE12.4)') 'q_soil(final) [mm/s]    = ', q_soil

      ! write(6,'(A,1PE12.4)') 'wdsrf_old [mm]          = ', wdsrf_old
      ! write(6,'(A,1PE12.4)') 'wdsrf_new [mm]          = ', wdsrf
      ! write(6,'(A,1PE12.4)') 'dwdsrf/dt [mm/s]        = ', q_store

      ! write(6,'(A,1PE12.4)') 'rsur(final) [mm/s]      = ', rsur
      ! write(6,'(A,1PE12.4)') 'CHECK: q_in-q_out-q_sto = ', chk_bal
   END SUBROUTINE SurfaceRunoff_TOPMOD_CLM

! -------------------------------------------------------------------------
   SUBROUTINE SubsurfaceRunoff_TOPMOD (nl_soil, icefrac, dz_soisno, zi_soisno, zwt, rsubst, &
         hksati, topoweti, eta)

   USE MOD_Namelist, only: DEF_TOPMOD_method
   IMPLICIT NONE

!-------------------------- Dummy Arguments ----------------------------
   integer,  intent(in) :: nl_soil                 !
   real(r8), intent(in) :: icefrac(1:nl_soil)      ! ice fraction (-)

   real(r8), intent(in) :: dz_soisno  (1:nl_soil)  ! layer depth (m)
   real(r8), intent(in) :: zi_soisno  (0:nl_soil)  ! interface level below a "z" level (m)

   real(r8), intent(in)  :: zwt    ! the depth from ground (soil) surface to water table [m]
   real(r8), intent(out) :: rsubst ! subsurface runoff (positive = out of soil column) (mm H2O /s)

   real(r8), intent(in), optional :: hksati (1:nl_soil)
   real(r8), intent(in), optional :: topoweti
   real(r8), intent(in), optional :: eta

!-------------------------- Local Variables ----------------------------

   real(r8), parameter :: vdcf = 2.0

   integer  :: j                ! indices
   integer  :: jwt              ! index of the soil layer right above the water table (-)
   real(r8) :: dzmm(1:nl_soil)  ! layer thickness (mm)

   real(r8) :: dzsum
   real(r8) :: icefracsum
   real(r8) :: fracice_rsub
   real(r8) :: imped
!-----------------------------------------------------------------------

      DO j = 1,nl_soil
         dzmm(j) = dz_soisno(j)*1000.
      ENDDO

      jwt = nl_soil
      ! allow jwt to equal zero when zwt is in top layer
      DO j = 1, nl_soil
         IF(zwt <= zi_soisno(j)) THEN
            jwt = j-1
            EXIT
         ENDIF
      ENDDO

      !-- Topographic runoff  --
      dzsum = 0.
      icefracsum = 0.
      DO j = max(jwt,1), nl_soil
         dzsum = dzsum + dzmm(j)
         icefracsum = icefracsum + icefrac(j) * dzmm(j)
      ENDDO
      ! add ice impedance factor to baseflow
      fracice_rsub = max(0.,exp(-3.*(1.-(icefracsum/dzsum)))-exp(-3.))/(1.0-exp(-3.))
      imped = max(0.,1.-fracice_rsub)

      IF ((DEF_TOPMOD_method == 1) .and. present(hksati) .and. present(topoweti)) THEN
         rsubst = imped * 3.e4 * sum(hksati(1:nl_soil))/nl_soil / vdcf * exp(-topoweti) * exp(-vdcf*zwt)
      ELSEIF ((DEF_TOPMOD_method == 2) .and. present(hksati) .and. present(eta)) THEN
         rsubst = imped * 3.e3 * sum(hksati(1:nl_soil))/nl_soil / vdcf * exp(-eta)
      ELSE
         rsubst = imped * 5.5e-3 * exp(-2.5*zwt)
      ENDIF

   END SUBROUTINE SubsurfaceRunoff_TOPMOD

! -------------------------------------------------------------------------
   SUBROUTINE Runoff_XinAnJiang ( &
         nl_soil, dz_soisno, eff_porosity, vol_liq, elvstd, gwat, deltim, &
         rsur, rsubst, frcsat)

   USE MOD_Precision
   IMPLICIT NONE

   integer,  intent(in) :: nl_soil ! number of soil layers

   real(r8), intent(in) :: &
        dz_soisno   (1:nl_soil),  &! layer thickness (m)
        eff_porosity(1:nl_soil),  &! effective porosity = porosity - vol_ice
        vol_liq     (1:nl_soil),  &! partial volume of liquid water in layer
        elvstd,                   &! standard deviation of elevation (m)
        gwat,                     &! net water input from top
        deltim                     ! time step (s)

   real(r8), intent(out) :: rsur   ! surface runoff (mm h2o/s)
   real(r8), intent(out) :: rsubst ! subsurface runoff (mm h2o/s)

   real(r8), intent(out), optional :: frcsat

   ! Local Variables
   real(r8) :: btopo, watin, w_int, wsat_int, wtmp, infil
   real(r8), parameter :: sigmin = 100.
   real(r8), parameter :: sigmax = 1000.

      watin = gwat * deltim / 1000.

      btopo = (elvstd - sigmin) / (elvstd + sigmax)
      btopo = min(max(btopo, 0.01), 0.5)

      w_int    = sum(vol_liq     (1:6) * dz_soisno(1:6))
      wsat_int = sum(eff_porosity(1:6) * dz_soisno(1:6))

      w_int = max(min(w_int, wsat_int), 0.)

      IF (present(frcsat)) THEN
         frcsat = 1.-(1.-w_int/wsat_int)**(btopo/(1.+btopo))
      ENDIF

      IF (watin <= 0.) THEN

         rsur   = 0.
         rsubst = 0.

      ELSE

         wtmp  = (1-w_int/wsat_int)**(1/(btopo+1)) - watin/((btopo+1)*wsat_int)
         infil = wsat_int - w_int - wsat_int * (max(0., wtmp))**(btopo+1)

         infil = min(infil, watin)

         rsur   = (watin - infil) * 1000. / deltim
         rsubst = 0.

      ENDIF

   END SUBROUTINE Runoff_XinAnJiang


   ! -------------------------------------------------------------------------
   SUBROUTINE Runoff_SimpleVIC ( &
      nl_soil, dz_soisno, eff_porosity, vol_liq, BVIC, gwat, deltim, &
      rsur, rsubst, frcsat)

   USE MOD_Precision
   IMPLICIT NONE

   integer,  intent(in) :: nl_soil ! number of soil layers

   real(r8), intent(in) :: &
      dz_soisno   (1:nl_soil),  &  ! layer thickness (m)
      eff_porosity(1:nl_soil),  &  ! effective porosity = porosity - vol_ice
      vol_liq     (1:nl_soil),  &  ! partial volume of liquid water in layer
      BVIC,                     &  ! VIC infiltration parameter
      gwat,                     &  ! net water input from top
      deltim                       ! time step (s)

   real(r8), intent(out) :: rsur   ! surface runoff (mm h2o/s)
   real(r8), intent(out) :: rsubst ! subsurface runoff (mm h2o/s)

   real(r8), intent(out), optional :: frcsat

   ! Local Variables
   real(r8) :: btopo, watin, w_int, wsat_int, wtmp, infil
   real(r8) :: InfilExpFac, WaterDepthMax, WaterDepthInit, RunoffSurface, InfilVarTmp
   real(r8) :: SoilSaturateFrac

      watin = gwat * deltim / 1000. ! convert mm/s to m

      w_int    = sum(vol_liq     (1:6) * dz_soisno(1:6))
      wsat_int = sum(eff_porosity(1:6) * dz_soisno(1:6))

      w_int = max(min(w_int, wsat_int), 0.)

      InfilExpFac = BVIC / ( 1.0 + BVIC )

      IF (present(frcsat)) THEN
         frcsat = 1.-(1.-w_int/wsat_int)**InfilExpFac
      ENDIF

      IF (watin <= 0.) THEN
         rsur   = 0.
         rsubst = 0.
      ELSE
         ! fractional saturated area from soil moisture
         SoilSaturateFrac = 1.0 - (max(0.0, (1.0-(w_int/wsat_int))))**InfilExpFac
         SoilSaturateFrac = max(0.0, SoilSaturateFrac)
         SoilSaturateFrac = min(1.0, SoilSaturateFrac)

         ! Infiltration for the previous time-step soil moisture based on SoilSaturateFrac
         WaterDepthMax  = (1.0 + BVIC) * wsat_int
         WaterDepthInit = WaterDepthMax * (1.0 - (1.0 - SoilSaturateFrac)**(1.0/BVIC))

         ! Solve for surface runoff
         if ( WaterDepthMax <= 0.0 ) then
            RunoffSurface = watin
         ELSEIF   ( (WaterDepthInit + watin) > WaterDepthMax ) then
            !RunoffSurface = (WaterDepthInit + w_int) - WaterDepthMax
            RunoffSurface = watin - wsat_int + w_int
         ELSE
            InfilVarTmp  = 1.0 - ((WaterDepthInit +watin ) / WaterDepthMax)
            RunoffSurface =watin - wsat_int + w_int + wsat_int * (InfilVarTmp**(1.0+BVIC))
         ENDIF

         IF ( RunoffSurface < 0.0 ) RunoffSurface = 0.0
         IF ( RunoffSurface > watin) RunoffSurface = watin

         infil = watin - RunoffSurface
         rsur= RunoffSurface * 1000. / deltim
         rsubst = 0.
      ENDIF

   END SUBROUTINE Runoff_SimpleVIC

END MODULE MOD_Runoff
! ---------- EOP ------------
