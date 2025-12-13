   SUBROUTINE SurfaceRunoff_TOPMOD_h2osfc (nl_soil,wimp,porsl,psi0,hksati,&
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

   END SUBROUTINE SurfaceRunoff_TOPMOD_h2osfc

   SUBROUTINE split_h2osfc_from_surface_soil(nl_soil, wimp, hksati, z_soisno, dz_soisno,   &
                                          zi_soisno, eff_porosity, icefrac,  &
                                          gwat, rsur, frcsat, slpratio, deltim, pondmx, f_h2osfc, wdsrf, &                                                                                 
                                          q_soil, q_excess, q_h2osfc, q_drain_h2osfc, q_h2osfc_surf)
   
   USE MOD_Const_ch4
   USE MOD_Vars_Global, only : PI
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
         slpratio,                 &! the slope ratioc
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