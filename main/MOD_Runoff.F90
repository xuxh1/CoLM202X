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
      qinmax = (1-fsat)*minval(10.**(-6.0*icefrac(1:min(3,nl_soil)))*hksati(1:min(3,nl_soil)))
      IF(eff_porosity(1)<wimp) qinmax = 0.

! Surface runoff
      rsur = fsat*max(0.0,gwat)

      IF (present(rsur_se)) THEN
         rsur_se = fsat*max(0.0,gwat)
      ENDIF

      IF (present(rsur_ie)) THEN
         rsur_ie = 0
      ENDIF

   END SUBROUTINE SurfaceRunoff_TOPMOD_CLM

   SUBROUTINE split_h2osfc_from_surface_soil(nl_soil, wimp, hksati, z_soisno, dz_soisno,   &
                                          zi_soisno, eff_porosity, icefrac,  &
                                          gwat, rsur, frcsat, slpratio, deltim, pondmx, f_h2osfc, wdsrf, &                                                                                 
                                          q_soil, q_excess, q_h2osfc, q_drain_h2osfc, q_h2osfc_surf)
   IMPLICIT NONE

   integer, intent(in) :: nl_soil   ! number of soil layers
   real(r8), intent(in) :: &
         ! wtfact,                 &! (updated to gridded 'fsatmax' data)
                                    ! fraction of model area with high water table
         wimp,                     &! water impermeable if porosity less than wimp
         hksati(1:nl_soil),        &! hydraulic conductivity at saturation (mm h2o/s)
         z_soisno(1:nl_soil),      &! layer depth (m)
         dz_soisno(1:nl_soil),     &! layer thickness (m)
         zi_soisno(0:nl_soil),     &! interface level below a "z" level (m)
         eff_porosity(1:nl_soil),  &! effective porosity = porosity - vol_ice
         icefrac(1:nl_soil),       &! ice fraction (-)
         gwat,                     &! net water input from top
         rsur,                     &! surface runoff (mm h2o/s)
         frcsat,                   &! fractional area with water table at surface
         slpratio,                 &! the slope ratio
         deltim,                   &
         pondmx

   real(r8), intent(inout) :: &
         f_h2osfc ,& ! fractional inundated area
         wdsrf       ! surface water (mm)

   real(r8), intent(out) :: q_soil, &
         q_excess, &
         q_h2osfc, &
         q_drain_h2osfc, &
         q_h2osfc_surf

   real(r8) :: qinmax
   real(r8) :: micro_sigma, sigma, d
   real(r8) :: fd, dfdd
   real(r8) :: pondmin
   real(r8) :: f_connected
   real(r8) :: k_h2osfc
   real(r8) :: 
   integer :: l
   
      !--------------------------------------------------
      ! 1. Maximum infiltration capacity
      !--------------------------------------------------
      qinmax = (1-frcsat)*minval(10.0_r8**(-6.0_r8*icefrac(1:min(3,nl_soil))) * &
                                 hksati(1:min(3,nl_soil)))

      IF(eff_porosity(1) < wimp) qinmax = 0.0_r8
      
      !--------------------------------------------------
      ! 2. Partition incoming water
      !--------------------------------------------------
      q_soil = (1.0_r8 - f_h2osfc) * (gwat - rsur)
      q_excess = max(q_soil - qinmax, 0.0_r8)
      q_h2osfc = f_h2osfc*(gwat - rsur) + q_excess

      !--------------------------------------------------
      ! 3. Connectivity function (CLM-style)
      !--------------------------------------------------
      IF (f_h2osfc <= 0.4_r8) THEN
         f_connected = 0.0_r8
      ELSE
         f_connected = (f_h2osfc - 0.4_r8)**0.14_r8
      END IF

      !--------------------------------------------------
      ! 4. h2osfc runoff related to slope
      !--------------------------------------------------
      if (wdsrf > pondmx) then
         k_h2osfc = 1.0e-4_r8 * sin(atan(slpratio))
         q_h2osfc_surf = k_h2osfc*f_connected*(wdsrf-pondmx)
         q_h2osfc_surf = min(q_h2osfc_surf, (wdsrf-pondmx)/deltim)
      else
         q_h2osfc_surf = 0._r8
      endif

      if (q_h2osfc_surf<1.e-8) q_h2osfc_surf = 0._r8

      wdsrf = wdsrf + (q_h2osfc - q_h2osfc_surf)*deltim

      !--------------------------------------------------
      ! 5. h2osfc drainage
      !--------------------------------------------------
      q_drain_h2osfc = min(f_h2osfc*qinmax,wdsrf/deltim)

      wdsrf = wdsrf - q_drain_h2osfc*deltim

      !--------------------------------------------------
      ! 6. Update inundation fraction (using CLM h2osfc scheme)
      !--------------------------------------------------
      micro_sigma = (atan(slpratio) + DEF_CH4_hydrology%slopemax**(1._r8/DEF_CH4_hydrology%slopebeta))**DEF_CH4_hydrology%slopebeta

      pondmin = 1.e-8_r8
      if (wdsrf > pondmin) then
         ! a cutoff is needed for numerical reasons...(nonconvergence after 5 iterations)
         d=0.0_r8

         sigma=1.0e3 * micro_sigma ! convert to mm
         do l=1,10
            fd = 0.5_r8*d*(1.0_r8+erf(d/(sigma*sqrt(2.0_r8)))) &
               +sigma/sqrt(2.0_r8*PI)*exp(-d**2/(2.0_r8*sigma**2)) &
               -wdsrf
            dfdd = 0.5_r8*(1.0_r8+erf(d/(sigma*sqrt(2.0_r8))))

            d = d - fd/dfdd
         enddo
         !--  update the submerged areal fraction using the new d value
         f_h2osfc = 0.5_r8*(1.0_r8+erf(d/(sigma*sqrt(2.0_r8))))

      else
         f_h2osfc = 0._r8
         ! The update of h2osfc is deferred to later, keeping with our standard
         ! separation of flux calculations from state updates, and because the state
         ! update needs to happen for tracers as well as bulk. However, it's important
         ! that this flux be applied soon after this routine, so that h2osfc remains in
         ! sync with frac_h2osfc.
      endif
   END SUBROUTINE
   
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
