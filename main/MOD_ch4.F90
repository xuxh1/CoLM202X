#include "define.h"

module MOD_ch4
    !=======================================================================
	! !DESCRIPTION:
	! Module holding routines to calculate methane fluxes
	! The driver averages up to gridcell, weighting by finundated, and checks for balance errors.
	! Sources, sinks, "competition" for CH4 & O2, & transport are resolved in ch4_tran.

	! !ORIGINAL:
	! The Community Land Model version 5.0 (CLM5.0)

	! !REFERENCES:
	! Lawrence, D.M., Fisher, R.A., Koven, C.D., Oleson, K.W., Swenson, S.C., Bonan, G., Collier, N., 
	! Ghimire, B., van Kampenhout, L., Kennedy, D. and Kluzek, E., 2019. 
	! The Community Land Model version 5: Description of new features, benchmarking,
	! and impact of forcing uncertainty. Journal of Advances in Modeling Earth Systems, 11(12), 4245-4287.

	! !REVISION:
	! Xionghui Xu, 2025, 1) Modify original CLM5 to be compatible with CoLM code structure. 
    !                    2) Fix some bugs based on the original code. 
	!=======================================================================
	use MOD_Precision
	use MOD_SPMD_Task
	use MOD_TimeManager
   	use MOD_Vars_TimeInvariants, only: wetwatmax
	use MOD_Vars_Global, only : maxsnl,nl_soil,nl_lake,spval,PI,deg2rad
	use MOD_Const_Physical, only: rgas, denh2o, denice, tfrz, grav
	use MOD_Const_ch4
	! use MOD_ch4varcon
	!-----------------------------------------------------------------------
	implicit none
	save

	INTERFACE print_var
		MODULE procedure print_var_int32
		MODULE procedure print_var_real8
   	END INTERFACE print_var

	public  :: ch4
	
	private :: ch4_annualupdate
	private :: ch4_prod
	private :: ch4_oxid
	private :: ch4_aere
	private :: ch4_ebul
	private :: ch4_tran

contains

	!-----------------------------------------------------------------------
	subroutine ch4 (idate,patchtype,&!input
		lb,snl,&
		dlon,dlat,&
		deltim,&
		z_soisno,dz_soisno,zi_soisno,t_soisno,t_grnd,wliq_soisno,wice_soisno,&
		forc_t,forc_pbot,forc_po2m,forc_pco2m,&
		zwt,rootfr,snowdp,wat,rsur,etr,wdsrf,bsw,&
		smp,porsl,lai,rootr,&
		annsum_npp,rr,&
		fsatmax,fsatdcf,frcsat,&
		agnpp,bgnpp,somhr,&
		crootfr,lithr,hr_vr,o_scalar,fphr,pot_f_nit_vr,pH,&
		cellorg,t_h2osfc,organic_max,&
		ch4_first_time,&
		!!!! --------------------------------------------------------------------------------------------------------
		!!!!                                         sum data   
		!!!! --------------------------------------------------------------------------------------------------------
		net_methane, &
		ch4_prod_depth, o2_decomp_depth, ch4_oxid_depth, o2_oxid_depth, &
		ch4_aere_depth, ch4_tran_depth, o2_aere_depth, ch4_ebul_depth, &
		o2stress, ch4stress, &
		ch4_surf_flux_tot, ch4_surf_aere, ch4_surf_ebul, ch4_surf_diff, &
		ch4_ebul_tot, ch4_prod_tot, ch4_oxid_tot, &
		totcolch4, grnd_ch4_cond, conc_o2, conc_ch4, &
		!!!! --------------------------------------------------------------------------------------------------------
		!!!! --------------------------------------------------------------------------------------------------------
		!!!!                                         sum data (unsaturated / saturated)
		!!!! --------------------------------------------------------------------------------------------------------
		net_methane_unsat, net_methane_sat, &
		ch4_prod_depth_unsat, ch4_prod_depth_sat, o2_decomp_depth_unsat, o2_decomp_depth_sat, &
		ch4_oxid_depth_unsat, ch4_oxid_depth_sat, o2_oxid_depth_unsat, o2_oxid_depth_sat, &
		ch4_aere_depth_unsat, ch4_aere_depth_sat, ch4_tran_depth_unsat, ch4_tran_depth_sat, &
		o2_aere_depth_unsat, o2_aere_depth_sat, ch4_ebul_depth_unsat, ch4_ebul_depth_sat, &
		o2stress_unsat, o2stress_sat, ch4stress_unsat, ch4stress_sat, &
		ch4_surf_flux_tot_unsat, ch4_surf_flux_tot_sat, ch4_surf_aere_unsat, ch4_surf_aere_sat, &
		ch4_surf_ebul_unsat, ch4_surf_ebul_sat, ch4_surf_diff_unsat, ch4_surf_diff_sat, &
		ch4_ebul_tot_unsat, ch4_ebul_tot_sat, ch4_prod_tot_unsat, ch4_prod_tot_sat, &
		ch4_oxid_tot_unsat, ch4_oxid_tot_sat, &
		totcolch4_unsat, totcolch4_sat, grnd_ch4_cond_unsat, grnd_ch4_cond_sat, &
		conc_o2_unsat, conc_o2_sat, conc_ch4_unsat, conc_ch4_sat, &
		!!!! --------------------------------------------------------------------------------------------------------
		c_atm, forc_pch4m, layer_sat_lag, lake_soilc, &
		annavg_agnpp, annavg_bgnpp, annavg_somhr, annavg_finrw, &
		tempavg_agnpp, tempavg_bgnpp, annsum_counter, tempavg_somhr, tempavg_finrw)

		!=======================================================================
		! !DESCRIPTION:
		! Driver for the methane emissions model
		!=======================================================================

		!===================== input ===========================================
		integer, intent(in) :: &
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			patchtype        , &! land patch type (0=soil, 1=urban or built-up, 2=wetland, 3=land ice, 4=land water bodies, 99=ocean)
			
			lb               , &! lower bound of array   (snl+1)
			snl				     ! number of snow layers (-5~-1)

		real(r8), intent(in) :: &
			dlon                    , &! latitude (degrees)
			dlat                    , &! longitude (degrees)

			deltim                  , &! land model time step [sec]
			z_soisno (maxsnl+1:nl_soil)    , &! layer depth [m]
			dz_soisno(maxsnl+1:nl_soil)    , &! layer thickness [m]
			zi_soisno(maxsnl:nl_soil)      , &! interface level below a "z" level [m]

			t_soisno (maxsnl+1:nl_soil)    , &! soil temperature [K]
			t_grnd                 		    , &! ground surface temperature [K]
			wliq_soisno(maxsnl+1:nl_soil)  , &! liquid water in layers [kg/m2]
			wice_soisno(maxsnl+1:nl_soil)  , &! ice lens in layers [kg/m2]

			forc_t                  , &! temperature at reference height [K]
			forc_pbot               , &! atm bottom level pressure (or reference height) [Pa]
			forc_po2m               , &! O2 concentration in atmos. [Pa]
			forc_pco2m              , &! CO2 concentration in atmos. [Pa]

			zwt                     , &! the depth from ground (soil) surface to water table [m]

			rootfr   (1:nl_soil)    , &! fraction of roots in each soil layer (the sum of all layer is 1)

			snowdp                  , &! snow depth [m]
			wat                     , &! total water storage [mm]
			rsur                    , &! surface runoff [mm/s]
			etr                     , &! transpiration rate [mm/s]
			wdsrf                   , &! depth of surface water [mm]
			bsw      (1:nl_soil)   	, &! Clapp and Hornberger "b" (nlevgrnd)             

			smp      (1:nl_soil)    , &! soil matrix potential [mm]
			porsl    (1:nl_soil)    , &! volumetric soil water at saturation (porosity)
			lai                     , &! leaf area index [m2/m2]

			annsum_npp              , &! annual sum NPP (g C/m2/yr)
			rr                      , &! root respiration (fine root MR + total root GR) (gC/m2/s)

			fsatmax                 , &! maximum saturated area fraction [-]
			fsatdcf                 , &! decay factor in calculation of saturated area fraction [1/m]
			frcsat                     ! fraction of saturation area
		!------------------- ch4_annualupdate ------------------------------------
		real(r8), intent(in) :: &			
			agnpp                   , &! aboveground NPP (gC/m2/s)
			bgnpp                   , &! belowground NPP (gC/m2/s)
			somhr                      ! (gC/m2/s) soil organic matter heterotrophic respiration

		!------------------- ch4_prod --------------------------------------------
		real(r8), intent(in) :: &
			crootfr  (1:nl_soil)    , &! fraction of roots for carbon in each soil layer (the sum of all layer is 1)
			lithr                   , &! (gC/m2/s) litter heterotrophic respiration        
			hr_vr    (1:nl_soil)    , &! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
			o_scalar (1:nl_soil)    , &! fraction by which decomposition is limited by DEF_CH4%anoxia
			fphr     (1:nl_soil)    , &! fraction of potential heterotrophic respiration 

			pot_f_nit_vr(1:nl_soil) , &! (gN/m3/s) potential soil nitrification flux 
			pH                         ! soil water pH                                     

		!------------------- ch4_aere --------------------------------------------
		real(r8), intent(in) :: &
			rootr    (1:nl_soil)       ! effective fraction of roots in each soil layer (the sum of all layer is 1)
		
		!------------------- ch4_tran --------------------------------------------
		real(r8), intent(in) :: &
			cellorg  (1:nl_soil)   		, &! column 3D org (kg/m^3 organic matter)
			t_h2osfc               		, &! surface water temperature               
			organic_max               		! organic matter content (kg m-3) where soil is assumed to act like peat

		logical, intent(inout) ::&
			ch4_first_time

		!!!! --------------------------------------------------------------------------------------------------------
		!!!!                                         sum data   
		!!!! --------------------------------------------------------------------------------------------------------
		!------------------- ch4_flux, balance, depth variables ------------------------------
		real(r8), intent(out) :: &
			net_methane                     , & ! average net methane correction to CO2 flux (mol/m2/s)
			ch4_prod_depth    (1:nl_soil)   , & ! production of CH4 in each soil layer (mol/m3/s)
			o2_decomp_depth   (1:nl_soil)   , & ! O2 consumption during decomposition in each soil layer (mol/m3/s)
			ch4_oxid_depth    (1:nl_soil)   , & ! CH4 consumption rate via oxidation in each soil layer (mol/m3/s)
			o2_oxid_depth     (1:nl_soil)   , & ! O2 consumption rate via oxidation in each soil layer (mol/m3/s)
			ch4_aere_depth    (1:nl_soil)   , & ! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s)
			ch4_tran_depth    (1:nl_soil)   , & ! CH4 loss rate via transpiration in each soil layer (mol/m3/s)
			o2_aere_depth     (1:nl_soil)   , & ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s)
			ch4_ebul_depth    (1:nl_soil)   , & ! CH4 loss rate via ebullition in each soil layer (mol/m3/s)
			o2stress          (1:nl_soil)   , & ! Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs
			ch4stress         (1:nl_soil)   , & ! Ratio of methane available to total per-timestep methane sinks
			ch4_surf_flux_tot               , & ! CH4 flux to atmosphere (mol/m2/s)
			ch4_surf_aere                   , & ! CH4 surface flux via aerenchyma (mol/m2/s)
			ch4_surf_ebul                   , & ! CH4 ebullition flux (mol/m2/s)
			ch4_surf_diff                   , & ! CH4 diffusion flux (mol/m2/s)
			ch4_ebul_tot                    , & ! Total CH4 ebullition (mol/m2/s)
			ch4_prod_tot                    , & ! Total CH4 production (mol/m2/s)
			ch4_oxid_tot                        ! Total CH4 oxidation (mol/m2/s)

		!------------------- total and concentration variables ------------------------------
		real(r8), intent(inout) :: &
			totcolch4               , & ! total methane in soil column (mol/m2)
			grnd_ch4_cond           , & ! tracer conductance for boundary layer [m/s]
			conc_o2  (1:nl_soil)    , & ! O2 conc in each soil layer (mol/m3)
			conc_ch4 (1:nl_soil)        ! CH4 conc in each soil layer (mol/m3)
		!!!! --------------------------------------------------------------------------------------------------------

		!!!! --------------------------------------------------------------------------------------------------------
		!!!!                                         sum data (unsaturated / saturated)
		!!!! --------------------------------------------------------------------------------------------------------
		!------------------- ch4_flux, balance, depth variables ------------------------------
		real(r8), intent(out) :: &
			net_methane_unsat               , & ! average unsaturated net methane correction to CO2 flux (mol/m2/s)
			net_methane_sat                 , & ! average saturated net methane correction to CO2 flux (mol/m2/s)
			ch4_prod_depth_unsat (1:nl_soil), & ! CH4 production rate in unsaturated soil layer (mol/m3/s)
			ch4_prod_depth_sat   (1:nl_soil), & ! CH4 production rate in saturated soil layer (mol/m3/s)
			o2_decomp_depth_unsat(1:nl_soil), & ! O2 consumption during decomposition (unsaturated) (mol/m3/s)
			o2_decomp_depth_sat  (1:nl_soil), & ! O2 consumption during decomposition (saturated) (mol/m3/s)
			ch4_oxid_depth_unsat (1:nl_soil), & ! CH4 oxidation rate in unsaturated soil layer (mol/m3/s)
			ch4_oxid_depth_sat   (1:nl_soil), & ! CH4 oxidation rate in saturated soil layer (mol/m3/s)
			o2_oxid_depth_unsat  (1:nl_soil), & ! O2 oxidation rate in unsaturated soil layer (mol/m3/s)
			o2_oxid_depth_sat    (1:nl_soil), & ! O2 oxidation rate in saturated soil layer (mol/m3/s)
			ch4_aere_depth_unsat (1:nl_soil), & ! CH4 loss rate via aerenchyma (unsaturated) (mol/m3/s)
			ch4_aere_depth_sat   (1:nl_soil), & ! CH4 loss rate via aerenchyma (saturated) (mol/m3/s)
			ch4_tran_depth_unsat (1:nl_soil), & ! CH4 loss rate via transpiration (unsaturated) (mol/m3/s)
			ch4_tran_depth_sat   (1:nl_soil), & ! CH4 loss rate via transpiration (saturated) (mol/m3/s)
			o2_aere_depth_unsat  (1:nl_soil), & ! O2 gain via aerenchyma (unsaturated) (mol/m3/s)
			o2_aere_depth_sat    (1:nl_soil), & ! O2 gain via aerenchyma (saturated) (mol/m3/s)
			ch4_ebul_depth_unsat (1:nl_soil), & ! CH4 ebullition loss (unsaturated) (mol/m3/s)
			ch4_ebul_depth_sat   (1:nl_soil), & ! CH4 ebullition loss (saturated) (mol/m3/s)
			o2stress_unsat       (1:nl_soil), & ! O2 stress ratio (unsaturated)
			o2stress_sat         (1:nl_soil), & ! O2 stress ratio (saturated)
			ch4stress_unsat      (1:nl_soil), & ! CH4 stress ratio (unsaturated)
			ch4stress_sat        (1:nl_soil), & ! CH4 stress ratio (saturated)
			ch4_surf_flux_tot_unsat         , & ! CH4 surface flux to atmosphere (unsaturated) (mol/m2/s)
			ch4_surf_flux_tot_sat           , & ! CH4 surface flux to atmosphere (saturated) (mol/m2/s)
			ch4_surf_aere_unsat             , & ! CH4 surface flux via aerenchyma (unsaturated) (mol/m2/s)
			ch4_surf_aere_sat               , & ! CH4 surface flux via aerenchyma (saturated) (mol/m2/s)
			ch4_surf_ebul_unsat             , & ! CH4 ebullition flux (unsaturated) (mol/m2/s)
			ch4_surf_ebul_sat               , & ! CH4 ebullition flux (saturated) (mol/m2/s)
			ch4_surf_diff_unsat             , & ! CH4 diffusion flux (unsaturated) (mol/m2/s)
			ch4_surf_diff_sat               , & ! CH4 diffusion flux (saturated) (mol/m2/s)
			ch4_ebul_tot_unsat              , & ! Total CH4 ebullition (unsaturated) (mol/m2/s)
			ch4_ebul_tot_sat                , & ! Total CH4 ebullition (saturated) (mol/m2/s)
			ch4_prod_tot_unsat              , & ! Total CH4 production (unsaturated) (mol/m2/s)
			ch4_prod_tot_sat                , & ! Total CH4 production (saturated) (mol/m2/s)
			ch4_oxid_tot_unsat              , & ! Total CH4 oxidation (unsaturated) (mol/m2/s)
			ch4_oxid_tot_sat                    ! Total CH4 oxidation (saturated) (mol/m2/s)

		!------------------- total and concentration variables ------------------------------
		real(r8), intent(inout) :: &
			totcolch4_unsat         , & ! total methane in soil column (unsaturated) (mol/m2)
			totcolch4_sat           , & ! total methane in soil column (saturated) (mol/m2)
			grnd_ch4_cond_unsat     , & ! tracer conductance for boundary layer (unsaturated) [m/s]
			grnd_ch4_cond_sat       , & ! tracer conductance for boundary layer (saturated) [m/s]
			conc_o2_unsat (1:nl_soil), & ! O2 conc in unsaturated soil layer (mol/m3)
			conc_o2_sat   (1:nl_soil), & ! O2 conc in saturated soil layer (mol/m3)
			conc_ch4_unsat(1:nl_soil), & ! CH4 conc in unsaturated soil layer (mol/m3)
			conc_ch4_sat  (1:nl_soil)   ! CH4 conc in saturated soil layer (mol/m3)
		!!!! --------------------------------------------------------------------------------------------------------

		!------------------- atmospheric and structural variables ------------------------------
		real(r8), intent(out) :: &
			c_atm      (1:3)             , & ! CH4, O2, CO2 atmospheric conc (mol/m3)
			forc_pch4m                    , & ! CH4 concentration in atmosphere (Pa)
			layer_sat_lag(1:nl_soil)      , & ! lagged saturation ratio per layer
			lake_soilc  (1:nl_soil)         ! total soil organic matter per layer (gC/m3)

		!------------------- annual accumulators ------------------------------
		real(r8), intent(out) :: &
			annavg_agnpp            , & ! annual average above-ground NPP (gC/m2/s)
			annavg_bgnpp            , & ! annual average below-ground NPP (gC/m2/s)
			annavg_somhr            , & ! annual average SOM heterotrophic respiration (gC/m2/s)
			annavg_finrw              ! respiration-weighted annual average of inundated zones (gC/m2/s)

		!------------------- temporary accumulators ------------------------------
		real(r8), intent(inout) :: &
			tempavg_agnpp           , & ! temporary average above-ground NPP (gC/m2/s)
			tempavg_bgnpp           , & ! temporary average below-ground NPP (gC/m2/s)
			annsum_counter          , & ! seconds since last annual accumulator turnover
			tempavg_somhr           , & ! temporary average SOM heterotrophic respiration (gC/m2/s)
			tempavg_finrw              ! respiration-weighted temporary average of inundated zones (gC/m2/s)


		!=================== Local Variables ============================================
		integer  :: i,j,s                     ! indices
		integer  :: sat                     ! 0 = unsatured, 1 = saturated
		integer  :: finundated              ! fractional inundated area, =sat(0 or 1)
		integer  :: jwt                     ! index of the soil layer right above the water table (-)
		real(r8) :: lon,lat                 ! lon,lat
		real(r8) :: total                   ! diff + aere + ebul
		real(r8) :: total_sat                   ! diff + aere + ebul
		real(r8) :: total_unsat                 ! diff + aere + ebul

		real(r8) :: dfsat
		real(r8) :: fsat_bef                ! finundated from previous timestep
		real(r8) :: errch4                  ! g C / m^2
		real(r8) :: redoxlags_vertical      ! Vertical redox lag time in s
		integer  :: dummyfilter(1)          ! empty filter
		real(r8) :: totcolch4_bef           ! total methane in soil column, start of timestep (mol/m2)

		real(r8) :: k_h_cc(0:nl_soil,ngases)! ratio of mol/m3 in liquid to mol/m3 in gas [-]

		real(r8) :: vol_aqu  (1:nl_soil)    ! liquid volumetric water content ---- water volume/all volume [m3/m3]
		real(r8) :: vol_gas  (1:nl_soil)    ! air volumetric water content ---- air volume/all volume [m3/m3]
		! real(r8) :: vol_sol  (1:nl_soil)    ! ice volumetric water content ---- ice volume/all volume [m3/m3]
		real(r8) :: f_aqu    (1:nl_soil)    ! water-filled proportion
		real(r8) :: f_gas    (1:nl_soil)    ! air-filled proportion
		! real(r8) :: f_sol    (1:nl_soil)    ! ice-filled proportion

		real(r8) :: conc_ch4_gas     (1:nl_soil) ! gas phase CH4 conc in each soil layer (mol/m3)
		real(r8) :: conc_ch4_aqu     (1:nl_soil) ! aqueous phase CH4 conc in each soil layer (mol/m3)
		! real(r8) :: conc_ch4_sol     (1:nl_soil) ! solid phase CH4 conc in each soil layer (mol/m3)
		real(r8) :: conc_ch4_porsl     (1:nl_soil) ! CH4 conc in each porosity (mol/m3)
		real(r8) :: conc_ch4_gas_porsl (1:nl_soil) ! gas phase CH4 conc in each porosity (mol/m3)
		real(r8) :: conc_ch4_aqu_porsl (1:nl_soil) ! aqueous phase CH4 conc in each porosity (mol/m3)
		! real(r8) :: conc_ch4_sol_porsl (1:nl_soil) ! solid phase CH4 conc in each porosity (mol/m3)

		real(r8) :: conc_o2_gas     (1:nl_soil) ! gas phase O2 conc in each soil layer (mol/m3)
		real(r8) :: conc_o2_aqu     (1:nl_soil) ! aqueous phase O2 conc in each soil layer (mol/m3)
		! real(r8) :: conc_o2_sol     (1:nl_soil) ! solid phase O2 conc in each soil layer (mol/m3)
		real(r8) :: conc_o2_porsl     (1:nl_soil) ! O2 conc in each porosity (mol/m3)
		real(r8) :: conc_o2_gas_porsl (1:nl_soil) ! gas phase O2 conc in each porosity (mol/m3)
		real(r8) :: conc_o2_aqu_porsl (1:nl_soil) ! aqueous phase O2 conc in each porosity (mol/m3)
		! real(r8) :: conc_o2_sol_porsl (1:nl_soil) ! solid phase O2 conc in each porosity (mol/m3)
		real(r8) :: err
	   	real(r8) :: fsat_wetland         ! fractional area with water table at surface

		real(r8) :: zwt_sat, wice_soisno_sat(1:nl_soil), wliq_soisno_sat(1:nl_soil), wdsrf_sat
		real(r8) :: zwt_unsat, wice_soisno_unsat(1:nl_soil), wliq_soisno_unsat(1:nl_soil), wdsrf_unsat
		!-----------------------------------------------------------------------

		! Set parameters
		redoxlags_vertical = DEF_CH4%redoxlag_vertical*secspday ! days --> s
		
		totcolch4_bef = totcolch4
		totcolch4 = 0
		totcolch4_sat   = 0
		totcolch4_unsat = 0

		finundated = frcsat
		! Call print_var(totcolch4_bef,'ch4 totcolch4_bef',idate)

		! ! fsat_wetland = fsatmax * exp(- fsatdcf * DEF_CH4_hydrology%vdcf * zwt)
		! ! fsat_wetland = 0.38 * exp(- 0.5 * DEF_CH4_hydrology%vdcf * zwt)
		! Call print_var(fsat_wetland,'ch4 fsat_wetland',idate)
		! Call print_var(fsatdcf,'ch4 fsatdcf',idate)
		! Call print_var(fsatmax,'ch4 fsatmax',idate)

		! Initialize fluxes to zero
		ch4_surf_flux_tot       = 0._r8
		ch4_surf_flux_tot_sat   = 0._r8
		ch4_surf_flux_tot_unsat = 0._r8

		ch4_prod_tot            = 0._r8
		ch4_prod_tot_sat        = 0._r8
		ch4_prod_tot_unsat      = 0._r8

		ch4_oxid_tot            = 0._r8
		ch4_oxid_tot_sat        = 0._r8
		ch4_oxid_tot_unsat      = 0._r8

		! Adjustment to NEE for methane production - oxidation
		net_methane             = 0._r8
		net_methane_sat         = 0._r8
		net_methane_unsat       = 0._r8

		! Check if offline. If offline, the default atmospheric methane concentration will be adopted globally (1700ppb)
		if (DEF_CH4%ch4offline) then
			forc_pch4m = DEF_CH4%atmch4*forc_pbot
			! [Pa]     =[mol/mol]*[Pa]
		else
			if (forc_pch4m == 0._r8) then
				write(6,*) 'not using DEF_CH4%ch4offline, but methane concentration not passed from the atmosphere', &
				'to land model! CoLM Model is stopping.'
				CALL CoLM_stop ()
			end if
		end if
		Call print_var(forc_pbot,'ch4 forc_pbot',idate)
		c_atm(1) =  forc_pch4m / rgasm / forc_t
		c_atm(2) =  forc_po2m  / rgasm / forc_t
		c_atm(3) =  forc_pco2m / rgasm / forc_t
		! n/V = P/RT
		![mol/m3]=[Pa]         /[J/K/mol]/[K]
		![mol/m3]=[J/m3]       *[K*mol/J]*[1/K]
		do i = 1,3
			Call print_var(c_atm(i),'ch4 c_atm',idate,i)
		enddo
		!!!! Begin biochemistry
		! First for soil
		! Do CH4 Annual Averages
		call ch4_annualupdate(idate, finundated, deltim,  agnpp, bgnpp, somhr, &
			annavg_agnpp, annavg_bgnpp, annavg_somhr,  annavg_finrw, &
			tempavg_agnpp,tempavg_bgnpp,annsum_counter,tempavg_somhr, tempavg_finrw)

		call henry_law(t_grnd,t_soisno,k_h_cc)

		layer_sat_lag = 1.
		!-------------------------------------------------
		! Loop
		!-------------------------------------------------
		do sat= 0, 1
			if (sat==0) then ! unsaturated
				zwt_unsat = zwt
				wliq_soisno_unsat = wliq_soisno
				wice_soisno_unsat = wice_soisno
				wdsrf_unsat = wdsrf
				jwt_unsat = nl_soil
				! allow jwt to equal zero when zwt is in top layer
				do j = 1, nl_soil
					if(zwt_unsat <= zi_soisno(j)) then
						jwt_unsat = j-1
						exit
					end if
				end do

				do j=1,nl_soil
					if (DEF_CH4%use_vertical_redoxlag .and. j > jwt .and. redoxlags_vertical > 0._r8) then ! saturated currently
						layer_sat_lag(j) = layer_sat_lag(j) * exp(-deltim/redoxlags_vertical) &
							+ (1._r8 - exp(-deltim/redoxlags_vertical))
					else if (DEF_CH4%use_vertical_redoxlag .and. redoxlags_vertical > 0._r8) then
						layer_sat_lag(j) = layer_sat_lag(j) * exp(-deltim/redoxlags_vertical)
					else if (j > jwt) then  ! redoxlags_vertical = 0
						layer_sat_lag(j) = 1._r8
					else
						layer_sat_lag(j) = 0._r8
					end if
				end do

				call split_ch4_o2_phases( dz_soisno, wliq_soisno_unsat, porsl, &
					conc_ch4_unsat, conc_o2_unsat, k_h_cc, idate, &
					vol_aqu_unsat, vol_gas_unsat, f_aqu_unsat, f_gas_unsat, &
					conc_ch4_gas_unsat, conc_ch4_aqu_unsat, conc_ch4_porsl_unsat, conc_ch4_gas_porsl_unsat, conc_ch4_aqu_porsl_unsat, &
					conc_o2_gas_unsat, conc_o2_aqu_unsat, conc_o2_porsl_unsat, conc_o2_gas_porsl_unsat, conc_o2_aqu_porsl_unsat )

				! Calculate CH4 production in each soil layer
				call ch4_prod ( idate, patchtype, sat, finundated, jwt_unsat, rr, deltim, &
					z_soisno, dz_soisno, zi_soisno, t_soisno, &
					lai, conc_o2_unsat, rootfr, annavg_finrw, &
					crootfr, somhr, lithr, hr_vr, o_scalar, fphr, pot_f_nit_vr, pH, layer_sat_lag, &
					ch4_prod_depth_unsat, o2_decomp_depth_unsat )

				! Calculate CH4 oxidation in each soil layer
				call ch4_oxid ( idate, jwt_unsat, sat, t_soisno, smp, vol_aqu_unsat, &
					conc_o2_aqu_porsl_unsat, conc_ch4_aqu_porsl_unsat, &
					ch4_oxid_depth_unsat, o2_oxid_depth_unsat )

				! Calculate CH4 aerenchyma losses in each soil layer
				call ch4_aere ( idate, jwt_unsat, sat, lai, deltim, &
					z_soisno, dz_soisno, zi_soisno, t_soisno, &
					rootfr, rootr, etr, grnd_ch4_cond_unsat, c_atm, annsum_npp, &
					annavg_agnpp, annavg_bgnpp, conc_ch4_unsat, ch4_prod_depth_unsat, conc_ch4_aqu_porsl_unsat, conc_ch4_gas_porsl_unsat, conc_o2_aqu_porsl_unsat, conc_o2_gas_porsl_unsat, &
					ch4_aere_depth_unsat, ch4_tran_depth_unsat, o2_aere_depth_unsat )

				! Calculate CH4 ebullition losses in each soil layer
				call ch4_ebul ( idate, jwt_unsat, sat, finundated, deltim, &
					z_soisno, dz_soisno, zi_soisno, forc_pbot, &
					t_soisno, wdsrf_unsat, conc_ch4_unsat, conc_ch4_gas_porsl_unsat, &
					ch4_ebul_depth_unsat )

				! Solve CH4 reaction/diffusion equation 
				! Competition for oxygen will occur here.
				call ch4_tran ( idate, patchtype, &
					lb, snl, jwt_unsat, sat, finundated, &
					lon, lat, deltim, z_soisno, dz_soisno, zi_soisno, t_soisno, t_grnd, &
					porsl, wliq_soisno_unsat, wice_soisno_unsat, wdsrf_unsat, bsw, c_atm, ch4_prod_depth_unsat, o2_aere_depth_unsat, &
					cellorg, t_h2osfc, organic_max, k_h_cc, conc_ch4_gas_porsl_unsat, conc_ch4_aqu_porsl_unsat, conc_o2_gas_porsl_unsat, conc_o2_aqu_porsl_unsat, vol_aqu_unsat, vol_gas_unsat, &
					o2stress_unsat, ch4stress_unsat, ch4_surf_aere_unsat, ch4_surf_ebul_unsat, ch4_surf_diff_unsat, ch4_ebul_tot_unsat, &
					ch4_oxid_depth_unsat, ch4_aere_depth_unsat, ch4_ebul_depth_unsat, &
					grnd_ch4_cond_unsat, o2_oxid_depth_unsat, o2_decomp_depth_unsat, conc_o2_unsat, conc_ch4_unsat )

			elseif (sat==1) then ! saturated
				zwt_sat=0

				DO j = 1, nl_soil
					IF(t_soisno(j)>tfrz)THEN
						wliq_soisno_sat(j) = porsl(j)*dz_soisno(j)*denh2o
						wice_soisno_sat(j) = 0.
					ELSE
						wliq_soisno_sat(j) = 0.
						wice_soisno_sat(j) = porsl(j)*dz_soisno(j)*denice
					ENDIF
				ENDDO
			
				wdsrf_sat = 0.
				jwt_sat = 0

				call split_ch4_o2_phases( dz_soisno, wliq_soisno_sat, porsl, &
					conc_ch4_sat, conc_o2_sat, k_h_cc, idate, &
					vol_aqu_sat, vol_gas_sat, f_aqu_sat, f_gas_sat, &
					conc_ch4_gas_sat, conc_ch4_aqu_sat, conc_ch4_porsl_sat, conc_ch4_gas_porsl_sat, conc_ch4_aqu_porsl_sat, &
					conc_o2_gas_sat, conc_o2_aqu_sat, conc_o2_porsl_sat, conc_o2_gas_porsl_sat, conc_o2_aqu_porsl_sat )

				! Calculate CH4 production in each soil layer
				call ch4_prod ( idate, patchtype, sat, finundated, jwt_sat, rr, deltim, &
					z_soisno, dz_soisno, zi_soisno, t_soisno, &
					lai, conc_o2_sat, rootfr, annavg_finrw, &
					crootfr, somhr, lithr, hr_vr, o_scalar, fphr, pot_f_nit_vr, pH, layer_sat_lag, &
					ch4_prod_depth_sat, o2_decomp_depth_sat )

				! Calculate CH4 oxidation in each soil layer
				call ch4_oxid ( idate, jwt_sat, sat, t_soisno, smp, vol_aqu_sat, &
					conc_o2_aqu_porsl_sat, conc_ch4_aqu_porsl_sat, &
					ch4_oxid_depth_sat, o2_oxid_depth_sat )

				! Calculate CH4 aerenchyma losses in each soil layer
				call ch4_aere ( idate, jwt_sat, sat, lai, deltim, &
					z_soisno, dz_soisno, zi_soisno, t_soisno, &
					rootfr, rootr, etr, grnd_ch4_cond_sat, c_atm, annsum_npp, &
					annavg_agnpp, annavg_bgnpp, conc_ch4_sat, ch4_prod_depth_sat, conc_ch4_aqu_porsl_sat, conc_ch4_gas_porsl_sat, conc_o2_aqu_porsl_sat, conc_o2_gas_porsl_sat, &
					ch4_aere_depth_sat, ch4_tran_depth_sat, o2_aere_depth_sat )

				! Calculate CH4 ebullition losses in each soil layer
				call ch4_ebul ( idate, jwt_sat, sat, finundated, deltim, &
					z_soisno, dz_soisno, zi_soisno, forc_pbot, &
					t_soisno, wdsrf_sat, conc_ch4_sat, conc_ch4_gas_porsl_sat, &
					ch4_ebul_depth_sat )

				! Solve CH4 reaction/diffusion equation 
				! Competition for oxygen will occur here.
				call ch4_tran ( idate, patchtype, &
					lb, snl, jwt_sat, sat, finundated, &
					lon, lat, deltim, z_soisno, dz_soisno, zi_soisno, t_soisno, t_grnd, &
					porsl, wliq_soisno_sat, wice_soisno_sat, wdsrf_sat, bsw, c_atm, ch4_prod_depth_sat, o2_aere_depth_sat, &
					cellorg, t_h2osfc, organic_max, k_h_cc, conc_ch4_gas_porsl_sat, conc_ch4_aqu_porsl_sat, conc_o2_gas_porsl_sat, conc_o2_aqu_porsl_sat, vol_aqu_sat, vol_gas_sat, &
					o2stress_sat, ch4stress_sat, ch4_surf_aere_sat, ch4_surf_ebul_sat, ch4_surf_diff_sat, ch4_ebul_tot_sat, &
					ch4_oxid_depth_sat, ch4_aere_depth_sat, ch4_ebul_depth_sat, &
					grnd_ch4_cond_sat, o2_oxid_depth_sat, o2_decomp_depth_sat, conc_o2_sat, conc_ch4_sat )

			endif

		enddo


		do j=1,nl_soil
			if (j == 1) then
				ch4_surf_flux_tot_sat = ch4_surf_diff_sat + ch4_surf_aere_sat + ch4_surf_ebul_sat
				ch4_surf_flux_tot_unsat = ch4_surf_diff_unsat + ch4_surf_aere_unsat + ch4_surf_ebul_unsat
				ch4_surf_flux_tot = finundated*ch4_surf_flux_tot_sat + (1-finundated)*ch4_surf_flux_tot_unsat
				! Call print_var(total,'ch4 total',idate)
				! Call print_var(ch4_surf_diff,'ch4 ch4_surf_diff',idate)
				! Call print_var(ch4_surf_aere,'ch4 ch4_surf_aere',idate)
				! Call print_var(ch4_surf_ebul,'ch4 ch4_surf_ebul',idate)
			end if
			ch4_oxid_tot_unsat = ch4_oxid_tot_unsat + ch4_oxid_depth_unsat(j) * dz_soisno(j)
			ch4_prod_tot_unsat = ch4_prod_tot_unsat + ch4_prod_depth_unsat(j) * dz_soisno(j)

			ch4_oxid_tot_sat   = ch4_oxid_tot_sat   + ch4_oxid_depth_sat(j)   * dz_soisno(j)
			ch4_prod_tot_sat   = ch4_prod_tot_sat   + ch4_prod_depth_sat(j)   * dz_soisno(j)

			ch4_oxid_tot       = ch4_oxid_tot_sat * finundated + ch4_oxid_tot_unsat * (1.0_r8 - finundated)
			ch4_prod_tot       = ch4_prod_tot_sat * finundated + ch4_prod_tot_unsat * (1.0_r8 - finundated)
			! [mol/m2/s]    = [mol/m2/s] + [mol/m3/s]       * [m]
			! Call print_var(ch4_oxid_tot,'ch4 ch4_oxid_tot',idate,j)
			! Call print_var(ch4_prod_tot,'ch4 ch4_prod_tot',idate,j)
			! Call print_var(ch4_oxid_depth(j),'ch4 ch4_oxid_depth(j)',idate,j)
			! Call print_var(ch4_prod_depth(j),'ch4 ch4_prod_depth(j)',idate,j)

			if (j == nl_soil) then
				! Adjustment to NEE flux to atm. for methane production and oxidation
				net_methane_unsat = net_methane_unsat - ch4_prod_tot_unsat
				net_methane_unsat = net_methane_unsat + ch4_oxid_tot_unsat

				net_methane_sat   = net_methane_sat   - ch4_prod_tot_sat
				net_methane_sat   = net_methane_sat   + ch4_oxid_tot_sat

				! Combine unsaturated and saturated contributions
				net_methane = net_methane_sat * finundated + net_methane_unsat * (1.0_r8 - finundated)
			end if
		end do

		do j = 1, nl_soil
			! Accumulate total column CH4 for unsaturated and saturated zones
			totcolch4_unsat = totcolch4_unsat + conc_ch4_unsat(j) * dz_soisno(j)
			totcolch4_sat   = totcolch4_sat   + conc_ch4_sat(j)   * dz_soisno(j)
			! [mol/m2]      = [mol/m2]        + [mol/m3]          * [m]

			totcolch4 = totcolch4_sat * finundated + totcolch4_unsat * (1.0_r8 - finundated)
			! Call print_var(totcolch4,'ch4 totcolch4',idate,j)
			! Call print_var(conc_ch4(j),'ch4 2 conc_ch4(j)',idate,j)
		end do

		! Column level balance
		if (.not. ch4_first_time) then
			! Check balance
			errch4 = totcolch4 - totcolch4_bef - deltim*(ch4_prod_tot - ch4_oxid_tot - ch4_surf_flux_tot/fsat_wetland) 
			! [g CH4/m2]    = [g CH4/m2] - [g CH4/m2] + [s]*[g CH4/m2/s]
			if (abs(errch4) > 1.e-7_r8) then
				write(6,*)'Lat,Lon,Patchtype        = ', dlat,dlon, patchtype
				write(6,*)'totcolch4                = ', totcolch4
				write(6,*)'totcolch4_bef            = ', totcolch4_bef
				write(6,*)'deltim*ch4_prod_tot      = ', deltim*ch4_prod_tot
				write(6,*)'deltim*ch4_oxid_tot      = ', deltim*ch4_oxid_tot
				write(6,*)'deltim*ch4_surf_flux_tot = ', deltim*ch4_surf_flux_tot
				CALL CoLM_stop ()
			end if
		end if

		ch4_first_time = .false.
	end subroutine ch4

	!-----------------------------------------------------------------------
	subroutine ch4_annualupdate(idate, finundated, deltim,  agnpp, bgnpp, somhr, &
		annavg_agnpp, annavg_bgnpp, annavg_somhr,  annavg_finrw, &
		tempavg_agnpp,tempavg_bgnpp,annsum_counter,tempavg_somhr, tempavg_finrw)
		!-----------------------------------------------------------------------
		! !DESCRIPTION:
		! Annual mean fields. 
		! Only the annavg is useful, the tempavg is temp.
		!-----------------------------------------------------------------------
		use MOD_Precision
		implicit none
		!-----------------------Argument----------------------------------------
		integer, intent(in) :: &
			idate(3)             , &! model calendar for next time step (year, days of the year, seconds of the day)
			finundated              ! fractional inundated area, =sat(0 or 1)

		real(r8), intent(in) :: &
			deltim                  , &! land model time step (sec)

			agnpp                   , &! aboveground NPP (gC/m2/s)
			bgnpp                   , &! belowground NPP (gC/m2/s)
			somhr                      ! soil organic matter heterotrophic respiration (gC/m2/s)

		real(r8), intent(out) :: &
			annavg_agnpp            , &! annual average aboveground NPP (gC/m2/s)         
			annavg_bgnpp            , &! annual average belowground NPP (gC/m2/s)         
			annavg_somhr            , &! annual average SOM heterotrophic resp. (gC/m2/s)  
			annavg_finrw               ! respiration-weighted annual average of finundated (1e-2*%) 
    		! definition different with tempavg_finrw 
		
		real(r8), intent(inout) :: &
			! Cumulative data (from year start to now)   
			tempavg_agnpp           , &! temporary average aboveground NPP (gC/m2/s)   
			tempavg_bgnpp           , &! temporary average belowground NPP (gC/m2/s)      
			annsum_counter          , &! seconds since last annual accumulator turnover    
			tempavg_somhr           , &! temporary average SOM heterotrophic resp. (gC/m2/s)
			tempavg_finrw              ! respiration-weighted annual average of finundated (gC/m2/s)
		!-----------------------Local Variables------------------------------         
		real(r8):: secsperyear       ! total number of seconds this year
		!-----------------------------------------------------------------------
		! set time steps
		if ( isleapyear(idate(1)) ) then
			secsperyear = 366*secspday
		else
			secsperyear = 365*secspday
		endif

		annsum_counter = annsum_counter + deltim

		if (annsum_counter >= secsperyear) then

			annsum_counter = 0._r8

			! update annual average somhr
			annavg_somhr      =  tempavg_somhr
			tempavg_somhr     = 0._r8

			! update annual average finrw
			if (annavg_somhr > 0._r8) then
				annavg_finrw      =  tempavg_finrw / annavg_somhr
			else
				annavg_finrw      = 0._r8
			end if
			tempavg_finrw  = 0._r8

			annavg_agnpp = tempavg_agnpp
			tempavg_agnpp = 0._r8
			
			annavg_bgnpp = tempavg_bgnpp
			tempavg_bgnpp = 0._r8

		else

			tempavg_somhr  = tempavg_somhr + deltim/secsperyear * somhr
			tempavg_finrw  = tempavg_finrw + deltim/secsperyear * finundated * somhr

			tempavg_agnpp = tempavg_agnpp + deltim/secsperyear * agnpp
			tempavg_bgnpp = tempavg_bgnpp + deltim/secsperyear * bgnpp
		
		end if

	end subroutine ch4_annualupdate

	!-----------------------------------------------------------------------
	subroutine ch4_prod (idate,patchtype,sat,finundated,jwt,rr,deltim,& !input
		z_soisno,dz_soisno,zi_soisno,t_soisno,&
		lai,conc_o2,rootfr,annavg_finrw,&
		crootfr,somhr,lithr,hr_vr,o_scalar,fphr,pot_f_nit_vr,pH,layer_sat_lag,&
		ch4_prod_depth,o2_decomp_depth)!output
		!-----------------------------------------------------------------------
		! !DESCRIPTION:
		! Production is done below the water table, based on CN heterotrophic respiration.
		! O2 is consumed by roots & by heterotrophic aerobes.
		! Production is done separately for sat & unsat, and is adjusted for temperature, seasonal inundation,
		! pH (optional), & redox lag factor.
		!-----------------------------------------------------------------------

		!-----------------------Argument----------------------------------------  
		integer , intent(in) :: &
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			patchtype                   , &! land patch type (0=soil, 1=urban or built-up, 
			! 2=wetland, 3=land ice, 4=land water bodies, 99=ocean)
			sat                         , &! 0 = unsaturated; 1 = saturated 
			finundated                  , &! fractional inundated area in soil column 
			jwt                            ! index of the soil layer right above the water table (-) 
  
		real(r8), intent(in) :: &
			rr                          , &! root respiration (fine root MR + total root GR) (gC/m2/s)

			deltim                      , &! land model time step (sec)
			z_soisno (maxsnl+1:nl_soil) , &! layer depth (m)
			dz_soisno(maxsnl+1:nl_soil) , &! layer thickness (m)
			zi_soisno(maxsnl:nl_soil)   , &! interface level below a "z" level (m)
	
			t_soisno (maxsnl+1:nl_soil) , &! soil temperature (K)
	
			lai                         , &! leaf area index [m2/m2]
			conc_o2  (1:nl_soil)        , &! O2 conc in each soil layer (mol/m3) (nl_soil)   
			rootfr   (1:nl_soil)        , &! fraction of roots in each soil layer (the sum of all layer is 1)
	
			annavg_finrw                , &! respiration-weighted annual average of finundated 
	
			crootfr  (1:nl_soil)        , &! fraction of roots for carbon in each soil layer (the sum of all layer is 1)
	
			somhr                       , &! soil organic matter heterotrophic respiration (gC/m2/s)
			lithr                       , &! litter heterotrophic respiration (gC/m2/s)        
			hr_vr    (1:nl_soil)        , &! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
			o_scalar (1:nl_soil)        , &! fraction by which decomposition is limited by DEF_CH4%anoxia
			fphr     (1:nl_soil)        , &! fraction of potential heterotrophic respiration 
	
			pot_f_nit_vr  (1:nl_soil)   , &! potential soil nitrification flux (gN/m3/s)
	
			pH                          , &! soil water pH                                     
  			layer_sat_lag   (1:nl_soil)    ! Lagged saturation status of soil layer in the unsaturated zone (1 = sat)

		real(r8), intent(out) :: &            
			ch4_prod_depth (1:nl_soil)  , &! production of CH4 in each soil layer (nl_soil) (mol/m3/s)
			o2_decomp_depth(1:nl_soil)     ! O2 consumption during decomposition in each soil layer (nl_soil) (mol/m3/s)

		!-----------------------Local Variables---------------------------------
		integer  :: j,s              ! indices
		real(r8) :: base_decomp      ! base heterotrophic respiration rate [mol C/m2/s]
		real(r8) :: partition_z
		real(r8) :: pH_fact_ch4      ! pH factor in methane production
  
		! Factors for methanogen temperature dependence being greater than soil aerobes
		real(r8) :: f_ch4_adj        ! Adjusted DEF_CH4%f_ch4
		real(r8) :: t_fact_ch4       ! Temperature factor calculated using additional Q10
  
		! O2 limitation on decomposition and methanogenesis
		real(r8) :: seasonalfin      ! finundated in excess of respiration-weighted annual average
  
		! For calculating column average (rootfrac(j)*rr(j))
		real(r8) :: rr_vr(1:nl_soil) ! vertically resolved column-mean root respiration (g C/m2/s)
  
		real(r8) :: sif              ! (unitless) ratio applied to sat. prod. to account for seasonal inundation
  		!-----------------------------------------------------------------------  
		! PATCH loop to calculate vertically resolved column-averaged root respiration
		rr_vr(:) = 0.0_r8

		if (lai > 0._r8) then
			do j=1,nl_soil
				rr_vr(j) = rr_vr(j) + rr*crootfr(j)
				! [g C/m2/s]=[g C/m2/s] + [g C/m2/s]*[-]
			end do
		end if
		  
		partition_z = 1._r8
		base_decomp = 0.0_r8
  
		! column loop to partition decomposition_rate into each soil layer
		do j=1,nl_soil
			! Use soil heterotrophic respiration (based on Wania)
			base_decomp = (somhr+lithr) / catomw
			! [mol C/m2/s]=[g C/m2/s]   / [g C/mol C]

			! Multiply base_decomp by factor accounting for lower carbon stock in seasonally inundated areas than
			! if it were inundated all year.
			! This is to reduce emissions in seasonally inundated zones, because the eq.
			! C-flux will be less than predicted by a non-O2-lim model
			if (sat == 1) then
				sif = 1._r8
				if (.not. DEF_CH4%anoxia) then
					if (annavg_finrw /= spval) then
						seasonalfin = max(finundated-annavg_finrw, 0._r8)
						if (seasonalfin > 0._r8) then
							sif = (annavg_finrw + DEF_CH4%mino2lim*seasonalfin) / finundated
							base_decomp = base_decomp * sif
						end if
					end if
				end if ! DEF_CH4%anoxia
			end if

			! For sensitivity studies
			base_decomp = base_decomp * DEF_CH4%cnscalefactor
  
			! For all landunits, prevent production or oxygen consumption when soil is at or below freezing.
			! If using VERTSOILC, it is OK to use base_decomp as given because liquid water stress will limit decomp.
			! if (t_soisno(j) <= tfrz .and. (patchtype == 4)) base_decomp = 0._r8
			if (t_soisno(j) <= tfrz .and. (patchtype == 4)) base_decomp = 0._r8
  
			! depth dependence of production either from rootfr or decomp model
			if ( (somhr + lithr) > 0._r8) then
				partition_z = hr_vr(j) * dz_soisno(j) / (somhr + lithr)
				!   [-]     = [g C/m3/s]*[m]/[g C/m2/s]
			else
				partition_z = 1._r8
			end if


			! Adjust DEF_CH4%f_ch4 to account for the fact that methanogens may have a higher Q10 than aerobic decomposers.
			! Note this is crude and should ideally be applied to all anaerobic decomposition rather than just the
			! DEF_CH4%f_ch4.
			f_ch4_adj = 1.0_r8

			t_fact_ch4 = DEF_CH4%q10ch4**((t_soisno(j) - DEF_CH4%q10ch4base)/10._r8)
			! Adjust DEF_CH4%f_ch4 by the ratio
			f_ch4_adj = DEF_CH4%f_ch4 * t_fact_ch4

			! Remove CN nitrogen limitation, as methanogenesis is not N limited.
			! Also remove (low) moisture limitation
			if (DEF_CH4%ch4rmcnlim) then
				if (fphr(j) > 0._r8) then
					f_ch4_adj = f_ch4_adj / fphr(j)
				end if
			end if
	
	
			! If switched on, use pH factor for production based on spatial pH data defined in surface data.
			if (patchtype /= 4 .and. DEF_CH4%usephfact)then 
				if (  pH >  DEF_CH4%pHmin .and.pH <  DEF_CH4%pHmax) then
					pH_fact_ch4 = 10._r8**(-0.2235_r8*pH*pH + 2.7727_r8*pH - 8.6_r8)
					! fitted function using data from Dunfield et al. 1993  
					! Strictly less than one, with optimum at 6.5
					! From Lei Meng
					f_ch4_adj = f_ch4_adj * pH_fact_ch4
				end if
			else
				! if no data, then no pH effects
			end if
	
			! Redox factor     
			if ((patchtype /= 4) .and. sat==1 .and. finundated_lag < finundated)  then
				f_ch4_adj = f_ch4_adj * finundated_lag / finundated
			elseif (sat == 0 .and. j > jwt) then ! Assume lag in decay of alternative electron acceptors vertically
				f_ch4_adj = f_ch4_adj * layer_sat_lag(j)
			end if
			! Alternative electron acceptors will be consumed first after soil is inundated.
	
			f_ch4_adj = min(f_ch4_adj, 0.5_r8)
			! Must be less than 0.5 because otherwise the actual implied aerobic respiration would be negative.
			! The total of aer. respiration + methanogenesis must remain equal to the SOMHR calculated in CN,
			! so that the NEE is sensible. Even perfectly anaerobic conditions with no alternative
			! electron acceptors would predict no more than 0.5 b/c some oxygen is present in organic matter.
			! e.g. 2CH2O --> CH4 + CO2.	
	
			! Decomposition uses 1 mol O2 per mol CO2 produced (happens below WT also, to deplete O2 below WT)
			! o2_decomp_depth is the demand in the absense of O2 supply limitation, in addition to autotrophic respiration.
			! Competition will be done in ch4_prod
	
			o2_decomp_depth(j) = base_decomp * partition_z / dz_soisno(j)
			if (DEF_CH4%anoxia) then
				! Divide off o_scalar to use potential O2-unlimited HR to represent aerobe demand for oxygen competition
				if (o_scalar(j) > 0._r8) then
					o2_decomp_depth(j) = o2_decomp_depth(j) / o_scalar(j)
				end if
			end if ! DEF_CH4%anoxia
	    
			! Add root respiration
			if (patchtype /= 4) then
				o2_decomp_depth(j) = o2_decomp_depth(j) + rr_vr(j)/catomw/dz_soisno(j)
				! [mol/m3/s]       = [mol/m3/s]         + [g C/m2/s]/[g C/mol C]/[m]
			end if

			! Add oxygen demand for nitrification
			if (DEF_CH4%use_nitrif_denitrif) then
				o2_decomp_depth(j) = o2_decomp_depth(j) + pot_f_nit_vr(j) * 2.0_r8/14.0_r8
				! [mol/m3/s]       = [mol/m3/s]         + [g N/m3/s]/[g N/mol N]
			end if
	
			if (j  >  jwt) then ! Below the water table so anaerobic CH4 production can occur
				! partition decomposition to layer
				! turn into per volume-total by dz
				ch4_prod_depth(j) = f_ch4_adj * base_decomp * partition_z / dz_soisno(j)
				! [mol C/m3/s]    = [-]       * [mol C/m2/s]* [-]         / [m]
			else ! Above the WT
				if (DEF_CH4%anoxicmicrosites) then
					ch4_prod_depth(j) = f_ch4_adj * base_decomp * partition_z / dz_soisno(j) &
						/ (1._r8 + DEF_CH4%oxinhib*conc_o2(j))
				else
					ch4_prod_depth(j) = 0._r8 ! [mol/m3 total/s]
				endif ! DEF_CH4%anoxicmicrosites
			endif ! WT
			Call print_var(f_ch4_adj,'ch4_prod f_ch4_adj',idate,j)
			Call print_var(base_decomp,'ch4_prod base_decomp',idate,j)
			Call print_var(partition_z,'ch4_prod partition_z',idate,j)
			Call print_var(ch4_prod_depth(j),'ch4_prod ch4_prod_depth(j)',idate,j)
	
		end do ! nl_soil
	
	end subroutine ch4_prod

	!---------------------------------------------------------------------------
	subroutine ch4_oxid (idate,jwt,  sat, t_soisno, smp, vol_aqu, &
		conc_o2_aqu_porsl, conc_ch4_aqu_porsl, &
		ch4_oxid_depth, o2_oxid_depth) 
		!-----------------------------------------------------------------------
		! !DESCRIPTION:
		! Oxidation is based on double Michaelis-Mentin kinetics, and is adjusted for low soil moisture.
		! Oxidation will be limited by available oxygen in ch4_tran.
		!-----------------------------------------------------------------------

		!-----------------------Argument---------- -----------------------------
		integer , intent(in) :: &
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			jwt                    , &! index of the soil layer right above the water table (-) 
			sat                       ! 0 = unsaturated; 1 = saturated 

		real(r8), intent(in) :: &
			t_soisno (maxsnl+1:nl_soil)    , &! soil temperature (Kelvin)
			smp      (1:nl_soil)   , &! soil matrix potential [mm]
			vol_aqu  (1:nl_soil)   , &! liquid volumetric water content

			conc_o2_aqu_porsl  (1:nl_soil)   , &! aqueous phase O2 conc in each porosity (mol/m3)
			conc_ch4_aqu_porsl (1:nl_soil)      ! aqueous phase CH4 conc in each porosity (mol/m3)

		real(r8), intent(out) :: &
			ch4_oxid_depth (1:nl_soil)   , &! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) 
			o2_oxid_depth  (1:nl_soil)      ! O2 consumption rate via oxidation in each soil layer (mol/m3/s) 

		!-----------------------Local Variables---------------------------------
		integer :: j                              ! indices
		real(r8):: t0                             ! Base temperature for Q10
		real(r8):: oxid_a                         ! Oxidation predicted by method A (temperature & enzyme limited) (mol CH4/m3/s)
		real(r8):: smp_fact                       ! factor for reduction based on soil moisture (unitless)
		real(r8):: k_m_eff                        ! effective DEF_CH4%k_m
		real(r8):: vmax_eff                       ! effective vmax 
		!-----------------------------------------------------------------------
		t0 = tfrz + 12._r8 ! Walter, for Michigan site where the 45 M/h comes from

		! Loop to determine oxidation in each layer
		do j=1,nl_soil
			if (sat == 1 .or. j > jwt) then ! Below the water table
				! Literature (e.g. Bender & Conrad, 1992) suggests lower DEF_CH4%k_m and vmax for high-CH4-affinity methanotrophs in
				! upland soils consuming ambient methane.
				k_m_eff = DEF_CH4%k_m
				vmax_eff = DEF_CH4%vmax_ch4_oxid
			else
				k_m_eff = DEF_CH4%k_m_unsat
				vmax_eff = DEF_CH4%vmax_oxid_unsat
			end if

			if (j <= jwt .and. smp(j) < 0._r8) then
				smp_fact = exp(-smp(j)/DEF_CH4%smp_crit)
				! Schnell & King, 1996, Figure 3
			else
				smp_fact = 1._r8
			end if

			oxid_a              = vmax_eff     * vol_aqu(j)* conc_ch4_aqu_porsl(j) / (k_m_eff + conc_ch4_aqu_porsl(j)) &
			! [mol/m3/s]        = [mol/m3/s]   * [-]     [mol/m3-w]    [mol/m3-w]  [mol/m3-w]
				* conc_o2_aqu_porsl(j) / (DEF_CH4%k_m_o2 + conc_o2_aqu_porsl(j)) &
				* DEF_CH4%q10_ch4_oxid ** ((t_soisno(j) - t0) / 10._r8) * smp_fact

			! For all landunits / levels, prevent oxidation if at or below freezing
			if (t_soisno(j) <= tfrz) oxid_a = 0._r8

			ch4_oxid_depth(j) = oxid_a
			o2_oxid_depth(j) = ch4_oxid_depth(j) * 2._r8
			Call print_var(smp_fact,'ch4_oxid smp_fact',idate,j)
			Call print_var(conc_ch4_aqu_porsl(j),'ch4_oxid conc_ch4_aqu_porsl(j)',idate,j)
			Call print_var(conc_o2_aqu_porsl(j),'ch4_oxid conc_o2_aqu_porsl(j)',idate,j)
			Call print_var(vol_aqu(j),'ch4_oxid vol_aqu(j)',idate,j)
			Call print_var(ch4_oxid_depth(j),'ch4_oxid ch4_oxid_depth(j)',idate,j)
			Call print_var(o2_oxid_depth(j),'ch4_oxid o2_oxid_depth(j)',idate,j)
		end do  
	end subroutine ch4_oxid

		!---------------------------------------------------------------------------
	subroutine ch4_aere (idate,jwt, sat, lai, deltim,&
		z_soisno, dz_soisno, zi_soisno, t_soisno,&
		rootfr, rootr, etr, grnd_ch4_cond, c_atm, annsum_npp,&
		annavg_agnpp, annavg_bgnpp, conc_ch4, ch4_prod_depth,conc_ch4_aqu_porsl,conc_ch4_gas_porsl,conc_o2_aqu_porsl,conc_o2_gas_porsl,&
		ch4_aere_depth, ch4_tran_depth, o2_aere_depth)
		!-----------------------------------------------------------------------
		! !DESCRIPTION:
		! Arctic c3 grass (which is often present in fens) and all vegetation in inundated areas is assumed to have
		! some root porosity. Currently, root porosity is allowed to be different for grasses & non-grasses.
		! CH4 diffuses out and O2 diffuses into the soil.  CH4 is also lossed via transpiration, which is both
		! included in the "aere" variables and output separately.  In practice this value is small.
		! By default upland veg. has small 5% porosity but this can be switched to be equal to inundated porosity.
		!-----------------------------------------------------------------------

		!-----------------------Argument----------------------------------------
		integer , intent(in) :: &
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			jwt                    , &! index of the soil layer right above the water table (-) 
			sat                       ! 0 = unsatured, 1 = saturated 

		real(r8), intent(in) :: &
			lai                    , &! adjusted leaf area index for seasonal variation [-]

			deltim                 , &! land model time step (sec)
			z_soisno (maxsnl+1:nl_soil)    , &! layer depth (m)
			dz_soisno(maxsnl+1:nl_soil)    , &! layer thickness (m)
			zi_soisno(maxsnl:nl_soil)      , &! interface level below a "z" level (m)

			t_soisno (maxsnl+1:nl_soil)    , &! soil temperature (Kelvin)
			rootfr   (1:nl_soil)   , &! fraction of roots in each soil layer
			rootr    (1:nl_soil)   , &! effective fraction of roots in each soil layer (SMS method only)
			! rootr here for effective per-layer transpiration, which may not be the same as rootfr
			etr                    , &! transpiration rate [mm/s]
			grnd_ch4_cond          , &! tracer conductance for boundary layer [m/s]

			c_atm(3)               , &! CH4, O2, CO2 atmospheric conc  (mol/m3)

			! These variables help us swap between big-leaf and fates boundary conditions
			annsum_npp             , &! annual sum NPP (g C/m2/yr)
			annavg_agnpp           , &! annual avg aboveground NPP (gC/m2/s)
			annavg_bgnpp           , &! annual avg belowground NPP (gC/m2/s)

			! These variables help us swap between saturated and unsaturated boundary conditions
			conc_ch4 (1:nl_soil)        , &! CH4 conc in each soil layer (mol/m3) 
			ch4_prod_depth (1:nl_soil)  , &! production of CH4 in each soil layer (mol/m3/s) 

			conc_ch4_aqu_porsl (1:nl_soil) , &! aqueous phase CH4 conc in each porosity [mol/m3]
			conc_ch4_gas_porsl (1:nl_soil) , &! gas phase CH4 conc in each porosity [mol/m3]
			conc_o2_aqu_porsl  (1:nl_soil) , &! aqueous phase O2 conc in each porosity [mol/m3]
			conc_o2_gas_porsl  (1:nl_soil)    ! gas phase O2 conc in each porosity [mol/m3]


		real(r8), intent(out) :: &
			ch4_aere_depth  (1:nl_soil)  , &! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) 
			ch4_tran_depth  (1:nl_soil)  , &! CH4 loss rate via transpiration in each soil layer (mol/m3/s) 
			o2_aere_depth   (1:nl_soil)     ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s) 

		!-----------------------Local Variables---------------------------------
		integer  :: j              ! indices

		! ch4 aerenchyma parameters
		real(r8) :: tranloss(1:nl_soil)    ! loss due to transpiration (mol / m3 /s)
		real(r8) :: aere    (1:nl_soil) 
		real(r8) :: oxaere  (1:nl_soil)    ! (mol / m3 /s)

		real(r8) :: aeretran
		!-----------------------------------------------------------------------
		! Initialize ch4_aere_depth
		do j=1,nl_soil
			ch4_aere_depth(j) = 0._r8
			ch4_tran_depth(j) = 0._r8
			o2_aere_depth(j) = 0._r8
		end do


		call SiteOxAere(idate,jwt,  sat, lai, z_soisno, dz_soisno,  zi_soisno,  t_soisno,  &
			rootfr, rootr, grnd_ch4_cond, etr, &
			annsum_npp, annavg_agnpp, annavg_bgnpp, c_atm, conc_ch4_aqu_porsl, conc_ch4_gas_porsl, conc_o2_aqu_porsl,conc_o2_gas_porsl,&
			tranloss, aere, oxaere)

		do j = 1,nl_soil
			! Impose limitation based on available methane during timestep
			! By imposing the limitation here, don't allow aerenchyma access to methane from other Patches.
			aeretran = min(aere(j)+tranloss(j), conc_ch4(j)/deltim + ch4_prod_depth(j))
			ch4_aere_depth (j) = ch4_aere_depth(j) + aeretran
			ch4_tran_depth (j) = ch4_tran_depth(j) + min(tranloss(j), aeretran)
			o2_aere_depth  (j) = o2_aere_depth (j) + oxaere(j)
			Call print_var(aeretran,'ch4_aere aeretran',idate,j)
			Call print_var(aere(j),'ch4_aere aere(j)',idate,j)
			Call print_var(tranloss(j),'ch4_aere tranloss(j)',idate,j)
			Call print_var(conc_ch4(j),'ch4_aere conc_ch4(j)',idate,j)
			Call print_var(ch4_prod_depth(j),'ch4_aere ch4_prod_depth(j)',idate,j)
			Call print_var(ch4_aere_depth(j),'ch4_aere ch4_aere_depth(j)',idate,j)
			Call print_var(ch4_tran_depth(j),'ch4_aere ch4_tran_depth(j)',idate,j)
			Call print_var(o2_aere_depth(j),'ch4_aere o2_aere_depth(j)',idate,j)
		end do ! over levels

	end subroutine ch4_aere


	!--------------------------------------------------------------------------- 
	subroutine SiteOxAere(idate,jwt,  sat, lai, z_soisno, dz_soisno,  zi_soisno,  t_soisno,  &
	 	rootfr, rootr, grnd_ch4_cond, etr, &
		annsum_npp, annavg_agnpp, annavg_bgnpp, c_atm, conc_ch4_aqu_porsl, conc_ch4_gas_porsl, conc_o2_aqu_porsl,conc_o2_gas_porsl,&
		tranloss, aere, oxaere)
		!-----------------------------------------------------------------------
		! !DESCRIPTION:
		! Site(column) level fluxes for O2 gain rate via
		! aerenchyma and ch4 losss rates from transpiration
		!-----------------------------------------------------------------------

		!-----------------------Argument----------------------------------------
		integer , intent(in) :: &
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			jwt                    , &! index of the soil layer right above the water table (-) 
			sat                       ! 0 = unsatured, 1 = saturated 

		real(r8), intent(in) :: &
			lai                    , &! leaf area index [m2/m2]

			z_soisno (maxsnl+1:nl_soil)    , &! layer depth (m)
			dz_soisno(maxsnl+1:nl_soil)    , &! layer thickness (m)
			zi_soisno(maxsnl:nl_soil)      , &! interface level below a "z" level (m)

			t_soisno (maxsnl+1:nl_soil)    , &! soil temperature (Kelvin)
			rootfr   (1:nl_soil)   , &! fraction of roots in each soil layer
			rootr    (1:nl_soil)   , &! root resistance of a layer, all layers sum to 1
			grnd_ch4_cond          , &! tracer conductance for boundary layer [m/s] 
			etr                    , &! transpiration rate [mm/s]

			annsum_npp             , &! annual sum NPP (g C/m2/yr)
			annavg_agnpp           , &! annual average aboveground NPP (g C/m2/s)
			annavg_bgnpp           , &! annual average belowground NPP (g C/m2/s)

			c_atm(3)               , &! CH4, O2, CO2 atmospheric conc  (mol/m3)

			conc_ch4_aqu_porsl (1:nl_soil) , &! aqueous phase CH4 conc in each porosity [mol/m3]
			conc_ch4_gas_porsl (1:nl_soil) , &! gas phase CH4 conc in each porosity [mol/m3]
			conc_o2_aqu_porsl  (1:nl_soil) , &! aqueous phase O2 conc in each porosity [mol/m3]
			conc_o2_gas_porsl  (1:nl_soil)    ! gas phase O2 conc in each porosity [mol/m3]

		real(r8), intent(out) :: &
			tranloss        (1:nl_soil)  , &! CH4 in soil water tran rate via plant transpiration in each soil layer (mol/m3/s) 
			aere            (1:nl_soil)  , &! CH4 tran rate via aerenchyma in each soil layer (mol/m3/s) 
			oxaere          (1:nl_soil)     ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s) 

		!-----------------------Local Variables---------------------------------         
		integer  :: j
		real(r8) :: area_tiller ! cross-sectional area of tillers (m2/m2)
		real(r8) :: m_tiller 
		real(r8) :: n_tiller 
		real(r8) :: anpp        ! annual sum NPP (gC/m2/yr) 
		real(r8) :: nppratio    ! bg/sum NPP
		real(r8) :: aere_ch4_resis    ! aerenchyma resistance [s/m]
		real(r8) :: grnd_ch4_resis    ! boundary layer resistance [s/m]
		real(r8) :: aere_o2_resis    ! aerenchyma resistance [s/m]
		real(r8) :: grnd_o2_resis    ! boundary layer resistance [s/m]
		real(r8) :: aerecond    ! aerenchyma conductance [m/s]
		real(r8) :: poros_tiller_real
		real(r8), parameter :: smallnumber = 1.e-12_r8
		!-----------------------------------------------------------------------
		! This parameter is poorly constrained and should be done on a patch-specific basis...

		! Attn EK: This calculation of aerenchyma properties is very uncertain. Let's check in once all
		! the new components are in; if there is any tuning to be done to get a realistic global flux,
		! this would probably be the place.  We will have to document clearly in the Tech Note
		! any major changes from the Riley et al. 2011 version. (There are a few other minor ones.)

		anpp = annsum_npp 
		anpp = max(anpp, 0._r8) ! NPP can be negative b/c of consumption of storage pools

		if (annavg_agnpp /= spval .and. annavg_bgnpp /= spval .and. &
			annavg_agnpp > 0._r8 .and. annavg_bgnpp > 0._r8) then
			nppratio = annavg_bgnpp / (annavg_agnpp + annavg_bgnpp)
		else
			nppratio = 0.5_r8
		end if

		do j=1,nl_soil
			! Calculate transpiration loss
			if (DEF_CH4%transpirationloss .and. lai > 0) then
				tranloss(j) = conc_ch4_aqu_porsl(j) * rootr(j)*etr / dz_soisno(j) / 1000._r8
				! [mol/m3/s]= [mol/m3]           * [-]     *[mm/s]/ [m]        /   [mm/m]
				! Use rootr here for effective per-layer transpiration, which may not be the same as rootfr
				tranloss(j) = max(tranloss(j), 0._r8) ! in case transpiration is pathological
			else
				tranloss(j) = 0._r8
			end if

			! Calculate aerenchyma diffusion   
			if (j > jwt .and. t_soisno(j) > tfrz .and. lai > 0) then ! Below water table
				! Estimate area of tillers (see Wania thesis)
				!m_tiller = anpp * r_leaf_root * lai ! (4.17 Wania)
				!m_tiller = 600._r8 * 0.5_r8 * 2._r8  ! used to be 300
				! Note: this calculation is based on Arctic graminoids, and should be refined for woody plants, if not
				! done on a PFT-specific basis.

				m_tiller = anpp * nppratio * DEF_CH4%wet_lai  !replace the elai(p) by constant 4 (by Xiyan Xu, 05/2016)
				! anpp as the sum carbon storage    [g C/m2/yr] as [g C/m2]
				! [g C/m2] = [g C/m2] * [-] * [m2/m2]

				n_tiller = m_tiller / DEF_CH4%tiller_C
				! [tiller/m2] = [g C/m2]/ [g C/tiller]
				if (sat == 0) then ! unsaturate
					poros_tiller_real = DEF_CH4%poros_tiller_unsat
				else
					poros_tiller_real = DEF_CH4%poros_tiller
				end if

				area_tiller = DEF_CH4%scale_factor_aere * n_tiller * poros_tiller_real * PI * DEF_CH4%aere_radius**2._r8
				! [m2/m2]   = [-]               * [tiller/m2] * [-]       * [-]* [m2/tiller]
				Call print_var(nppratio,'SiteOxAere nppratio',idate,j)
				Call print_var(anpp,'SiteOxAere anpp',idate,j)
				Call print_var(m_tiller,'SiteOxAere m_tiller',idate,j)
				Call print_var(n_tiller,'SiteOxAere n_tiller',idate,j)
				Call print_var(area_tiller,'SiteOxAere area_tiller',idate,j)

				aere_ch4_resis = 1._r8/((area_tiller * rootfr(j) * d_con_g(1,1) * 1e-4_r8 / (z_soisno(j)*DEF_CH4%rob))+smallnumber)
				! [s/m]         = /([m2/m2]          * [-]       * [m2/s]                 / [m]         /[-])
				! Add in boundary layer resistance
				! grnd_ch4_resis = 1._r8/(grnd_ch4_cond+smallnumber)
				grnd_ch4_resis = 1._r8/(grnd_ch4_cond+smallnumber)
				aerecond = 1._r8/(aere_ch4_resis + grnd_ch4_resis)
				! aerecond = max(aerecond,1.e-8_r8)
				aere(j) = aerecond*(conc_ch4_gas_porsl(j) - c_atm(1)) / dz_soisno(j) 
				! [mol/m3/s] = [mol/m3]                      / [m]          / [s/m]
				!ZS: Added porsl & Henry's const.
				aere(j) = max(aere(j), 0._r8) ! prevent backwards diffusion
				Call print_var(aere_ch4_resis,'SiteOxAere aere_ch4_resis',idate,j)
				Call print_var(grnd_ch4_resis,'SiteOxAere grnd_ch4_resis',idate,j)
				Call print_var(grnd_ch4_cond,'SiteOxAere grnd_ch4_cond',idate)
				Call print_var(c_atm(1),'SiteOxAere c_atm(1)',idate,j)
				Call print_var(conc_ch4_gas_porsl(j),'SiteOxAere conc_ch4_gas_porsl(j)',idate,j)
				Call print_var(aere(j),'SiteOxAere aere(j)',idate,j)

				! Do oxygen diffusion into layer
				aere_o2_resis = 1._r8/((area_tiller * rootfr(j) * d_con_g(2,1) * 1e-4_r8 / (z_soisno(j)*DEF_CH4%rob)) + smallnumber)
				grnd_o2_resis = 1._r8/(grnd_ch4_cond+smallnumber)
		
				oxaere(j) = -(conc_o2_gas_porsl(j) - c_atm(2)) / (dz_soisno(j)*(aere_o2_resis + grnd_o2_resis)) ![mol/m3-total/s]
				oxaere(j) = max(oxaere(j), 0._r8)
				! Diffusion in is positive; prevent backwards diffusion
				if ( .not. DEF_CH4%use_aereoxid_prog ) then ! fixed aere oxid proportion; will be done in ch4_tran
					oxaere(j) = 0._r8
				end if
			else
				aere(j) = 0._r8
				oxaere(j) = 0._r8
			end if ! veg type, below water table, & above freezing
		end do

  	end subroutine SiteOxAere

	!---------------------------------------------------------------------------
	subroutine ch4_ebul (idate,jwt, sat,finundated, deltim, &
		z_soisno, dz_soisno, zi_soisno, forc_pbot, &
		t_soisno, wdsrf, conc_ch4, conc_ch4_gas_porsl,&
		ch4_ebul_depth)
		!-----------------------------------------------------------------------
		! !DESCRIPTION:
		! Bubbling is based on temperature & pressure dependent solubility (k_h_cc), 
		! with assumed proportion of bubbles
		! which are CH4, and assumed early nucleation at DEF_CH4%vgc_max sat (Wania).
		! Bubbles are released to the water table surface in ch4_tran.  
		!-----------------------------------------------------------------------

		!-----------------------Argument---------- -----------------------------
		integer , intent(in) :: &
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			jwt                        , &! index of the soil layer right above the water table (-) 
			sat                        , &! 0 = unsaturated; 1 = saturated 
			finundated

		real(r8), intent(in) :: &
			deltim                     , &! land model time step (sec)
			z_soisno (maxsnl+1:nl_soil), &! layer depth (m)
			dz_soisno(maxsnl+1:nl_soil), &! layer thickness (m)
			zi_soisno(maxsnl:nl_soil)  , &! interface level below a "z" level (m)

			forc_pbot                  , &! atm bottom level pressure (or reference height) (Pa)
			t_soisno (maxsnl+1:nl_soil), &! soil temperature (Kelvin)
			wdsrf                      , &! depth of surface water [mm]
			conc_ch4       (1:nl_soil) , &! CH4 conc in each soil layer (mol/m3)

		   	conc_ch4_gas_porsl (1:nl_soil)! gas phase CH4 conc in each porosity (mol/m3)

		real(r8), intent(out) :: &
			ch4_ebul_depth (1:nl_soil)    ! CH4 loss rate via ebullition in each soil layer (mol/m3/s)

		!-----------------------Local Variables---------------------------------
		integer :: j      ! indices

		real(r8) :: vgc     ! gas phase volumetric CH4 content (m3 CH4/m3 pore air)
		real(r8) :: vgc_min ! minimum aqueous CH4 content when ebullition ceases
		real(r8) :: k_h     ! 
		real(r8) :: pressure! sum atmospheric and hydrostatic pressure
		real(r8) :: ebul_timescale
		!-----------------------------------------------------------------------
		vgc_min = DEF_CH4%vgc_max
		ebul_timescale = deltim ! Allow fast bubbling

		! column loop to estimate ebullition CH4 flux from each soil layer
		do j=1,nl_soil
			if (j  >  jwt .and. t_soisno(j) > tfrz) then ! Ebullition occurs only below the water table
				pressure = forc_pbot + denh2o * grav * (z_soisno(j)-zi_soisno(jwt)) 
				! [Pa]   = [Pa]      + [kg/m3]* [m/s2]* [m]
				! [Pa]   = [N/m2] = [kg]*[m/s2]/[m2] = [kg/m/s2]
				if (sat == 1 .and. finundated>0._r8) then ! Add ponding pressure head
					pressure = pressure + denh2o * grav * wdsrf/1000._r8/finundated
					! [Pa]   = [Pa]     + [kg/m3]* [m/s2]* [mm]/[mm/m]
				end if

				! Compare partial pressure to ambient pressure.
				vgc = conc_ch4_gas_porsl(j) * rgasm * t_soisno(j) / pressure
				! [-]= [mol/m3]           * [Pa*m3/K/mol]*[K]   / [Pa]

				if (vgc > DEF_CH4%vgc_max * DEF_CH4%bubble_f) then ! If greater than max value, remove amount down to vgc_min
					ch4_ebul_depth(j) = (vgc - vgc_min * DEF_CH4%bubble_f) * conc_ch4(j) / ebul_timescale
					! [mol/m3/s]      = [-]                        * [mol/m3]    / [s]
				else
					ch4_ebul_depth(j) = 0._r8
				endif
				Call print_var(ch4_ebul_depth(j),'ch4_ebul ch4_ebul_depth(j)',idate,j)
				Call print_var(vgc,'ch4_ebul vgc',idate,j)
				Call print_var(conc_ch4(j),'ch4_ebul conc_ch4(j)',idate,j)
				Call print_var(conc_ch4_gas_porsl(j),'ch4_ebul conc_ch4_gas_porsl(j)',idate,j)
				Call print_var(t_soisno(j),'ch4_ebul t_soisno(j)',idate,j)
				Call print_var(pressure,'ch4_ebul pressure',idate,j)
				Call print_var(wdsrf,'ch4_ebul wdsrf',idate,j)
				Call print_var(forc_pbot,'ch4_ebul forc_pbot',idate,j)

			else ! above the water table or freezing
				ch4_ebul_depth(j) = 0._r8
			endif ! below the water table and not freezing
		end do ! j

	end subroutine ch4_ebul

		!---------------------------------------------------------------------------
	subroutine ch4_tran (idate,patchtype, &
		lb, snl, jwt, sat, finundated,&
		lon, lat, deltim, z_soisno, dz_soisno, zi_soisno,  t_soisno, t_grnd, &
	 	porsl, wliq_soisno, wice_soisno, wdsrf, bsw, c_atm, ch4_prod_depth, o2_aere_depth,&
		cellorg,t_h2osfc, organic_max, k_h_cc, conc_ch4_gas_porsl,conc_ch4_aqu_porsl,conc_o2_gas_porsl,conc_o2_aqu_porsl,vol_aqu,vol_gas,&
		o2stress, ch4stress, ch4_surf_aere, ch4_surf_ebul, ch4_surf_diff, ch4_ebul_tot, &
		ch4_oxid_depth, ch4_aere_depth, ch4_ebul_depth, &
		grnd_ch4_cond, o2_oxid_depth, o2_decomp_depth, conc_o2, conc_ch4 )
		!-----------------------------------------------------------------------
		! !DESCRIPTION:
		! Solves the reaction & diffusion equation for the timestep.  
		! 1  "Competition" between processes for CH4 & O2 demand is done.  
		! 2  Concentrations are apportioned into gas & liquid fractions; 
		!    only the gas fraction is considered for diffusion in unsat.  
		! 3  Snow and lake water resistance to diffusion is added as a bulk term in the ground conductance 
		!    (which is really a surface layer conductance), but concentrations are not tracked and oxidation 
		!    is not allowed inside snow and lake water.
		! 4  Diffusivity is set based on soil texture and organic matter fraction. 
		!    A Crank-Nicholson solution is used.
		! 5  CH4 diffusive flux is calculated and consistency is checked.
		!-----------------------------------------------------------------------

		!-----------------------Argument----------------------------------------
		integer, intent(in) :: &
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			patchtype        	! land patch type (0=soil, 1=urban or built-up, 2=wetland,
										! 3=land ice, 4=land water bodies, 99=ocean
			! istep             , &! the i time step


		integer , intent(in) :: &
			lb                , &! lower bound of array (snl+1)
			snl				  , &!  number of snow layers     (-5~-1)
			jwt               , &! index of the soil layer right above the water table (-) 
			sat               , &! 0 = unsaturated; 1 = saturated 
			finundated

		real(r8), intent(in) :: &
			lon   	   				        , &! logitude 
			lat     	   			        , &! latitude 
 
			deltim                  	    , &! land model time step (sec)
			z_soisno (maxsnl+1:nl_soil)   	, &! layer depth (m)
			dz_soisno(maxsnl+1:nl_soil)   	, &! layer thickness (m)
			zi_soisno(maxsnl:nl_soil)   	, &! interface level below a "z" level (m)


			t_soisno (maxsnl+1:nl_soil)    	, &! soil temperature (Kelvin)
			t_grnd                 		    , &! ground surface temperature [k]

			porsl             (1:nl_soil)   , &! volumetric soil water at saturation (porosity)
			wliq_soisno(maxsnl+1:nl_soil)	, &! liquid water in layers [kg/m2]
			wice_soisno(maxsnl+1:nl_soil) 	, &! ice lens in layers [kg/m2]
			wdsrf                  		    , &! depth of surface water [mm]
			bsw               (1:nl_soil)   , &! Clapp and Hornberger "b" (nlevgrnd)             

			c_atm(3)               		    , &! CH4, O2, CO2 atmospheric conc  (mol/m3)

			ch4_prod_depth    (1:nl_soil)   , &! production of CH4 in each soil layer (mol/m3/s) 
			o2_aere_depth     (1:nl_soil)   , &! O2 gain rate via aerenchyma in each soil layer (mol/m3/s) 


			cellorg           (1:nl_soil)   , &! column 3D org (kg/m^3 organic matter)
			t_h2osfc               		    , &! surface water temperature               
			organic_max               	    , &! organic matter content (kg m-3) where soil is assumed to act like peat

			k_h_cc(0:nl_soil,ngases)        , &! ratio of mol/m3 in liquid to mol/m3 in gas
			conc_ch4_gas_porsl(1:nl_soil)   , &! gas phase CH4 conc in each porosity (mol/m3)
			conc_ch4_aqu_porsl(1:nl_soil)   , &! aqueous phase CH4 conc in each porosity (mol/m3)
			conc_o2_gas_porsl (1:nl_soil)   , &! gas phase O2 conc in each porosity (mol/m3)
			conc_o2_aqu_porsl (1:nl_soil)   , &! aqueous phase O2 conc in each porosity (mol/m3)
			vol_aqu           (1:nl_soil)   , &
			vol_gas           (1:nl_soil)

		real(r8), intent(out) :: &
			o2stress          (1:nl_soil)   , &! Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs
			ch4stress         (1:nl_soil)   , &! Ratio of methane available to the total per-timestep methane sinks 
			ch4_surf_aere                   , &! Total column CH4 aerenchyma (mol/m2/s)
			ch4_surf_ebul                   , &! CH4 ebullition to atmosphere (mol/m2/s)
			ch4_surf_diff                   , &! CH4 surface flux (mol/m2/s)
			ch4_ebul_tot                    ! Total column CH4 ebullition (mol/m2/s)

		real(r8), intent(inout) :: &
			ch4_oxid_depth    (1:nl_soil)   , &! InOut: CH4 consumption rate via oxidation in each soil layer (mol/m3/s) 
			ch4_aere_depth    (1:nl_soil)   , &! InOut: CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) 
			ch4_ebul_depth    (1:nl_soil)   , &! InOut: CH4 loss rate via ebullition in each soil layer (mol/m3/s)
			o2_oxid_depth     (1:nl_soil)   , &! InOut: O2 loss rate via ebullition in each soil layer (mol/m3/s) 
			o2_decomp_depth   (1:nl_soil)   , &! InOut: O2 consumption during decomposition in each soil layer (mol/m3/s)

			grnd_ch4_cond                   , &! InOut: tracer conductance for boundary layer [m/s]  
			conc_o2           (1:nl_soil)   , &! InOut: O2 conc in each soil layer (mol/m3) 
			conc_ch4          (1:nl_soil)      ! InOut: CH4 conc in each soil layer (mol/m3) 

  		!-----------------------Local Variables------------------------------
		integer :: j,s,i			                                               ! indices
		integer :: jtop                                                        ! top level at each column
		real(r8) :: at (0:nl_soil)                     ! "a" vector for tridiagonal matrix
		real(r8) :: bt (0:nl_soil)                     ! "b" vector for tridiagonal matrix
		real(r8) :: ct (0:nl_soil)                     ! "c" vector for tridiagonal matrix
		real(r8) :: rt (0:nl_soil)                     ! "r" vector for tridiagonal solution
		real(r8) :: f_a                                                        ! air-filled fraction of available pore space
		real(r8) :: diffus (0:nl_soil)                 ! diffusivity (m2/s)
		real(r8) :: dzj                                                        ! 
		real(r8) :: dp1_zp1 (0:nl_soil)                ! diffusivity/delta_z for next j
		real(r8) :: dm1_zm1 (0:nl_soil)                ! diffusivity/delta_z for previous j
		real(r8) :: t_soisno_c                                                 ! soil temperature   (maxsnl+1:nl_soil)
		real(r8) :: deficit                                                    ! mol CH4 /m^2 that must be subtracted from diffusive flux to atm. to make up
		! for keeping concentrations always above zero
		real(r8) :: conc_ch4_bef(1:nl_soil)            ! concentration at the beginning of the timestep
		real(r8) :: errch4                            ! Error (Mol CH4 /m^2) [+ = too much CH4]
		real(r8) :: conc_ch4_rel(0:nl_soil)            ! Concentration per volume of air or water
		real(r8) :: conc_o2_rel(0:nl_soil)             ! Concentration per volume of air or water
		real(r8) :: conc_ch4_rel_old(0:nl_soil)        ! Concentration during last Crank-Nich. loop
		real(r8), parameter :: smallnumber = 1.e-12_r8
		real(r8) :: snowdiff                                                   ! snow diffusivity (m^2/s)
		real(r8) :: snow_resis                           ! Cumulative Snow resistance (s/m). Also includes
		real(r8) :: pond_resis                                                    ! Additional resistance from ponding, up to pondmx water on top of top soil layer (s/m)
		real(r8) :: pondz                                                      ! Depth of ponding (m)
		real(r8) :: ponddiff                                                   ! Pondwater diffusivity (m^2/s)
		real(r8) :: spec_grnd_cond(1:ngases)            ! species grnd conductance (s/m)
		real(r8) :: airfrac                                                    ! air fraction in snow
		real(r8) :: waterfrac                                                  ! water fraction in snow
		real(r8) :: icefrac                                                    ! ice fraction in snow
		real(r8) :: epsilon_t (1:nl_soil,1:ngases)     !
		real(r8) :: epsilon_t_old (1:nl_soil,1:ngases) ! epsilon_t from last time step !Currently deprecated
		real(r8) :: source (1:nl_soil,1:ngases)        ! source
		real(r8) :: source_old (1:nl_soil,1:ngases)    ! source from last time step !Currently deprecated
		real(r8) :: om_frac                                                    ! organic matter fraction
		real(r8) :: o2demand, ch4demand                                        ! mol/m^3/s

		real(r8) :: err1,err2,err3,err4,err5,err6,err7,err8,err9,err10
  		!-----------------------------------------------------------------------
		! Perform competition for oxygen and methane in each soil layer if demands over the course of the timestep
		! exceed that available. Assign to each process in proportion to the quantity demanded in the absense of
		! the limitation.
		do j = 1,nl_soil
			o2demand = o2_decomp_depth(j) + o2_oxid_depth(j) ! o2_decomp_depth includes autotrophic root respiration
			if (o2demand > 0._r8) then
				o2stress(j) = min((conc_o2(j) / deltim + o2_aere_depth(j)) / o2demand, 1._r8)
			else
				o2stress(j) = 1._r8
			end if
	
			ch4demand = ch4_oxid_depth(j) + ch4_aere_depth(j) + ch4_ebul_depth(j)
			if (ch4demand > 0._r8) then
				ch4stress(j) = min((conc_ch4(j) / deltim + ch4_prod_depth(j)) / ch4demand, 1._r8)
			else
				ch4stress(j) = 1._r8
			end if

			Call print_var(o2stress(j),'ch4_tran 1 o2stress(j)',idate,j)
			Call print_var(ch4stress(j),'ch4_tran 1 ch4stress(j)',idate,j)
			Call print_var(o2demand,'ch4_tran o2demand',idate,j)
			Call print_var(ch4demand,'ch4_tran ch4demand',idate,j)
			Call print_var(ch4_oxid_depth(j),'ch4_tran 1 ch4_oxid_depth(j)',idate,j)
			Call print_var(o2_oxid_depth(j),'ch4_tran 1 o2_oxid_depth(j)',idate,j)
			Call print_var(ch4_aere_depth(j),'ch4_tran 1 ch4_aere_depth(j)',idate,j)
			Call print_var(ch4_ebul_depth(j),'ch4_tran 1 ch4_ebul_depth(j)',idate,j)
			Call print_var(o2_decomp_depth(j),'ch4_tran 1 o2_decomp_depth(j)',idate,j)
			! Resolve methane oxidation
			if (o2stress(j) < 1._r8 .or. ch4stress(j) < 1._r8) then
				if (ch4stress(j) <= o2stress(j)) then 
					! methane limited
					if (o2stress(j) < 1._r8) then
						! Recalculate oxygen limitation
						o2demand = o2_decomp_depth(j)
						if (o2demand > 0._r8) then
							o2stress(j) = min((conc_o2(j)/deltim + o2_aere_depth(j) - ch4stress(j)*o2_oxid_depth(j))/o2demand, 1._r8)
						else
							o2stress(j) = 1._r8
						end if
					end if
					! Reset oxidation
					ch4_oxid_depth(j) = ch4_oxid_depth(j) * ch4stress(j)
					o2_oxid_depth(j)  = o2_oxid_depth(j) * ch4stress(j)
				else                                      
					! oxygen limited
					if (ch4stress(j) < 1._r8) then
						! Recalculate methane limitation
						ch4demand = ch4_aere_depth(j) + ch4_ebul_depth(j)
						if (ch4demand > 0._r8) then
							ch4stress(j) = min( (conc_ch4(j) / deltim + ch4_prod_depth(j) - &
									o2stress(j)*ch4_oxid_depth(j)) / ch4demand, 1._r8)
						else
							ch4stress(j) = 1._r8
						end if
					end if
					! Reset oxidation
					ch4_oxid_depth(j) = ch4_oxid_depth(j) * o2stress(j)
					o2_oxid_depth(j) = o2_oxid_depth(j) * o2stress(j)
				end if
			end if
	
			! Reset non-methanotroph demands
			ch4_aere_depth(j) = ch4_aere_depth(j) * ch4stress(j)
			ch4_ebul_depth(j) = ch4_ebul_depth(j) * ch4stress(j)
			o2_decomp_depth(j) = o2_decomp_depth(j) * o2stress(j)

			Call print_var(o2stress(j),'ch4_tran 2 o2stress(j)',idate,j)
			Call print_var(ch4stress(j),'ch4_tran 2 ch4stress(j)',idate,j)

			Call print_var(ch4_oxid_depth(j),'ch4_tran 2 ch4_oxid_depth(j)',idate,j)
			Call print_var(o2_oxid_depth(j),'ch4_tran 2 o2_oxid_depth(j)',idate,j)
			Call print_var(ch4_aere_depth(j),'ch4_tran 2 ch4_aere_depth(j)',idate,j)
			Call print_var(ch4_ebul_depth(j),'ch4_tran 2 ch4_ebul_depth(j)',idate,j)
			Call print_var(o2_decomp_depth(j),'ch4_tran 2 o2_decomp_depth(j)',idate,j)
		end do !j
  
  
		! Accumulate ebullition to place in first layer above water table, or directly to atmosphere
		do j = 1,nl_soil
			if (j == 1) ch4_ebul_tot = 0._r8
			ch4_ebul_tot = ch4_ebul_tot + ch4_ebul_depth(j) * dz_soisno(j)
			Call print_var(ch4_ebul_depth(j),'ch4_tran 3 ch4_ebul_depth(j)',idate,j)
			Call print_var(ch4_ebul_tot,'ch4_tran 3 ch4_ebul_tot',idate,j)
		end do
  
		! Set the source term for each species (no need to do j=0, since epsilon_t and source not used there)
		! Note that because of the semi-implicit diffusion and the 30 min timestep combined with explicit
		! sources, occasionally negative concentration will result. In this case it is brought to zero and the
		! surface flux is adjusted to conserve. This results in some inaccuracy as compared to a shorter timestep
		! or iterative solution.
		do j = 1,nl_soil
	
			if ( .not. DEF_CH4%use_aereoxid_prog ) then
				! First remove the CH4 oxidation that occurs at the base of root tissues (aere), and add to oxidation
				ch4_oxid_depth(j) = ch4_oxid_depth(j) + DEF_CH4%aereoxid * ch4_aere_depth(j)
				ch4_aere_depth(j) = ch4_aere_depth(j) - DEF_CH4%aereoxid * ch4_aere_depth(j)
			end if ! else oxygen is allowed to diffuse in via aerenchyma
	
			source(j,1) = ch4_prod_depth(j) - ch4_oxid_depth(j) - &
					ch4_aere_depth(j) - ch4_ebul_depth(j) ! [mol/m3-total/s]
			Call print_var(source(j,1),'ch4_tran source(j,1)',idate,j)

			! aerenchyma added to surface flux below
			! ebul added to soil depth just above WT
			if (source(j,1) + conc_ch4(j) / deltim < -1.e-12_r8)then 
	
				write(6,*) 'Methane demands exceed methane available. Error in methane competition (mol/m^3/s), j:', &
						source(j,1) + conc_ch4(j) / deltim, j
				write(6,*)'Lat,Lon=',lat,lon
				CALL CoLM_stop ()
	
			else if (ch4stress(j) < 1._r8 .and. source(j,1) + conc_ch4(j) / deltim > 1.e-12_r8) then  
	
				write(6,*) 'Methane limited, yet some left over. Error in methane competition (mol/m^3/s), j:', &
						source(j,1) + conc_ch4(j) / deltim, j
				write(6,*)'Lat,Lon=',lat,lon
				CALL CoLM_stop ()
	
			end if
	
			source(j,2) = -o2_oxid_depth(j) - o2_decomp_depth(j) + o2_aere_depth(j) ! O2 [mol/m3/s]
			Call print_var(source(j,2),'ch4_tran source(j,2)',idate,j)

			if (source(j,2) + conc_o2(j) / deltim < -1.e-12_r8) then
	
				write(6,*) 'Oxygen demands exceed oxygen available. Error in oxygen competition (mol/m^3/s), j:', &
						source(j,2) + conc_o2(j) / deltim, j
				write(6,*)'Lat,Lon=',lat,lon
				CALL CoLM_stop ()
	
			else if (o2stress(j) < 1._r8 .and. source(j,2) + conc_o2(j) / deltim > 1.e-12_r8) then
	
				write(6,*) 'Oxygen limited, yet some left over. Error in oxygen competition (mol/m^3/s), j:', &
						source(j,2) + conc_o2(j) / deltim, j
				write(6,*)'Lat,Lon=',lat,lon
				CALL CoLM_stop ()
	
			end if
			Call print_var(conc_ch4_bef(j),'ch4_tran 1 conc_ch4_bef(j)',idate,j)
			Call print_var(conc_ch4(j),'ch4_tran 1 conc_ch4(j)',idate,j,s,.True.)
			conc_ch4_bef(j) = conc_ch4(j) !For Balance Check

		enddo ! j
  
  
		! Accumulate aerenchyma to add directly to atmospheric flux
		do j = 1,nl_soil
			if (j==1) ch4_surf_aere = 0._r8
			ch4_surf_aere = ch4_surf_aere + ch4_aere_depth(j) * dz_soisno(j)
			Call print_var(ch4_aere_depth(j),'ch4_tran 3 ch4_aere_depth(j)',idate,j)
			Call print_var(ch4_surf_aere,'ch4_tran 3 ch4_surf_aere',idate,j)
		enddo
  
		! Add in ebullition to source at depth just above WT
		if (jwt /= 0) then
			source(jwt,1) = source(jwt,1) + ch4_ebul_tot/dz_soisno(jwt)
			Call print_var(source(jwt,1),'ch4_tran 3 source(jwt,1)',idate,j)
		endif

		! Calculate concentration relative to m^3 of air or water: needed for the diffusion
		do j = 0,nl_soil
			if (j == 0) then
				conc_ch4_rel(j) = c_atm(1)
				conc_o2_rel(j)  = c_atm(2)
			else
				if (j <= jwt) then  ! Above the WT
					do s =1,2
						epsilon_t(j,s) = porsl(j)- (1._r8-k_h_cc(j,s))*vol_aqu(j)
					end do
					conc_ch4_rel(j) = conc_ch4_gas_porsl(j)
					conc_o2_rel(j)  = conc_o2_gas_porsl(j)
					! Partition between the liquid and gas phases. The gas phase will drive the diffusion.
				else ! Below the WT
					do s =1,2
						epsilon_t(j,s) = porsl(j)
						Call print_var(epsilon_t(j,s),'ch4_tran epsilon_t(j,s)',idate,j,s,.True.)
					end do
					conc_ch4_rel(j) = conc_ch4(j)/epsilon_t(j,1)
					conc_o2_rel(j)  = conc_o2(j) /epsilon_t(j,2)
					Call print_var(conc_ch4_aqu_porsl(j),'ch4_tran conc_ch4_aqu_porsl(j)',idate,j,s,.True.)
					Call print_var(conc_o2_aqu_porsl(j),'ch4_tran conc_o2_aqu_porsl(j)',idate,j,s,.True.)
				end if
			end if
			Call print_var(conc_ch4_rel(j),'ch4_tran 1 conc_ch4_rel(j)',idate,j,s,.True.)
			Call print_var(conc_o2_rel(j),'ch4_tran 1 conc_o2_rel(j)',idate,j,s,.True.)
		end do

  
		! Loop over species
		do s = 1, 2 ! 1=CH4; 2=O2; 3=CO2
			! Adjust the grnd_ch4_cond to keep it positive, and add the snow resistance & pond resistance
			do j = maxsnl + 1,0
				if (j == maxsnl + 1) then
					if (grnd_ch4_cond < smallnumber .and. s==1) grnd_ch4_cond = smallnumber
					! Needed to prevent overflow when ground is frozen, e.g. for lakes
					snow_resis = 0._r8
				end if
				Call print_var(snl,'ch4_tran snl',idate,j)
	
				! Add snow resistance
				if (j >= snl + 1) then
					! For the ice layer, all = ice + water + air, no soil
					t_soisno_c = t_soisno(j) - tfrz
					icefrac = wice_soisno(j)/denice/dz_soisno(j)
					! [-]   = [kg/m2]       / [kg/m3] / [m]
					waterfrac = wliq_soisno(j)/denh2o/dz_soisno(j)
					airfrac = max(1._r8 - icefrac - waterfrac, 0._r8)
					Call print_var(t_soisno_c,'ch4_tran t_soisno_c',idate,j)
					Call print_var(icefrac,'ch4_tran icefrac',idate,j)
					Call print_var(waterfrac,'ch4_tran waterfrac',idate,j)
					Call print_var(airfrac,'ch4_tran airfrac',idate,j)

					! Calculate snow diffusivity
					if (airfrac > 0.05_r8) then
						! Use Millington-Quirk Expression, as hydraulic properties (bsw) not available
						snowdiff = (d_con_g(s,1) + d_con_g(s,2)*t_soisno_c) * 1.e-4_r8 * &
								airfrac**(10._r8/3._r8) / (airfrac+waterfrac)**2 &
								* DEF_CH4%scale_factor_gasdiff
					else !solute diffusion in water only, airfrac -> 0
						snowdiff = (airfrac+waterfrac)**DEF_CH4%satpow * (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
								* DEF_CH4%scale_factor_liqdiff
					end if
					Call print_var(snowdiff,'ch4_tran snowdiff',idate,j)
					Call print_var(snow_resis,'ch4_tran 1 snow_resis',idate,j)

					snowdiff = max(snowdiff, smallnumber)
					snow_resis = snow_resis + dz_soisno(j)/snowdiff
					Call print_var(snow_resis,'ch4_tran 2 snow_resis',idate,j)
				end if
	
				if (j == 0) then ! final loop
					! Add pond resistance
					pond_resis = 0._r8

					! First old pond formulation up to pondmx
					if (patchtype /= 4 .and. snl == 0 .and. (wliq_soisno(1)/(dz_soisno(1)*denh2o)) > porsl(1)) then
						t_soisno_c = t_soisno(1) - tfrz
						if (t_soisno(1) <= tfrz) then
							ponddiff = (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
									* (wliq_soisno(1)/denh2o+smallnumber)/ &
									(wliq_soisno(1)/denh2o+wice_soisno(1)/denice+smallnumber) &
									* DEF_CH4%scale_factor_liqdiff
						else ! Unfrozen
							ponddiff = (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
									* DEF_CH4%scale_factor_liqdiff
						end if
						pondz = dz_soisno(1) * ((wliq_soisno(1)/(dz_soisno(1)*denh2o)) - porsl(1))
						pond_resis = pondz / ponddiff
						Call print_var(pondz,'ch4_tran 1 pondz',idate,j)
						Call print_var(pond_resis,'ch4_tran 1 pond_resis',idate,j)
						Call print_var(ponddiff,'ch4_tran 1 ponddiff',idate,j)
						Call print_var(t_soisno_c,'ch4_tran 1 t_soisno_c',idate,j)
					end if

					! Now add new wdsrf form
					if (patchtype /= 4 .and. sat == 1 .and. finundated>0._r8) then
						if (t_h2osfc >= tfrz) then
							t_soisno_c = t_h2osfc - tfrz
							ponddiff = (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
									* DEF_CH4%scale_factor_liqdiff
							pondz = wdsrf / 1000._r8/finundated ! Assume all wdsrf corresponds to sat area
							! [m] = [mm]  /  [mm/m]
							pond_resis = pond_resis + pondz / ponddiff
						else if (wdsrf/finundated > DEF_CH4%capthick) then
							! assume surface ice is impermeable
							pond_resis = 1/smallnumber
						end if
						Call print_var(wdsrf,'ch4_tran 2 wdsrf',idate,j)
						Call print_var(pondz,'ch4_tran 2 pondz',idate,j)
						Call print_var(pond_resis,'ch4_tran 2 pond_resis',idate,j)
						Call print_var(ponddiff,'ch4_tran 2 ponddiff',idate,j)
						Call print_var(t_soisno_c,'ch4_tran 2 t_soisno_c',idate,j)
					end if

					spec_grnd_cond(s) = 1._r8/(1._r8/grnd_ch4_cond + snow_resis + pond_resis)
					Call print_var(grnd_ch4_cond,'ch4_tran grnd_ch4_cond',idate,j)
					Call print_var(snow_resis,'ch4_tran snow_resis',idate,j)
					Call print_var(pond_resis,'ch4_tran pond_resis',idate,j)
					Call print_var(spec_grnd_cond(s),'ch4_tran spec_grnd_cond(s)',idate,s)
				end if
			end do ! j
	
			! Determine gas diffusion and fraction of open pore (f_a)
			do j = 1,nl_soil
				t_soisno_c = t_soisno(j) - tfrz
	
				if (j <= jwt) then  ! Above the WT
					if (organic_max > 0._r8) then
						om_frac = min(DEF_CH4%om_frac_sf*cellorg(j)/organic_max, 1._r8)
						! Use first power, not square as in iniTimeConst
					else
						om_frac = 1._r8
					end if
					diffus (j) = (d_con_g(s,1) + d_con_g(s,2)*t_soisno_c) * 1.e-4_r8 * &
								(om_frac * vol_gas(j)**(10._r8/3._r8) / porsl(j)**2._r8 + &
								(1._r8-om_frac) * vol_gas(j)**2._r8 * (vol_gas(j)/porsl(j))**(3._r8 / bsw(j)) ) &
								* DEF_CH4%scale_factor_gasdiff
				else ! Below the WT use saturated diffusivity and only water in epsilon_t
					! Note the following is not currently corrected for the effect on diffusivity of excess ice in soil under
					! lakes (which is currently experimental only).
					diffus (j) = porsl(j)**DEF_CH4%satpow * (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
						* DEF_CH4%scale_factor_liqdiff
					if (t_soisno(j)<=tfrz) then
						diffus(j) = diffus(j)*(wliq_soisno(j)/denh2o+smallnumber)/ &
							(wliq_soisno(j)/denh2o+wice_soisno(j)/denice+smallnumber)
					end if
				end if ! Above/below the WT
				diffus(j) = max(diffus(j), smallnumber) ! Prevent overflow
				Call print_var(diffus(j) ,'ch4_tran diffus(j)',idate,j)
				Call print_var(t_soisno_c ,'ch4_tran t_soisno_c',idate,j)
				Call print_var(wliq_soisno(j) ,'ch4_tran wliq_soisno(j)',idate,j)
				Call print_var(wice_soisno(j) ,'ch4_tran wice_soisno(j)',idate,j)
			enddo ! j
	
			do j = 1,nl_soil
	
				! Set up coefficients for tridiagonal solver.
				if (j == 1 .and. j /= jwt .and. j /= jwt+1) then
					dm1_zm1(j) = 1._r8/(1._r8/spec_grnd_cond(s)+dz_soisno(j)/(diffus(j)*2._r8))
					! replace Diffusivity / Delta_z by conductance (grnd_ch4_cond) for top layer
					dp1_zp1(j) = 2._r8/(dz_soisno(j)/diffus(j)+dz_soisno(j+1)/diffus(j+1))
				else if (j == 1 .and. j == jwt) then
					dm1_zm1(j) = 1._r8/(1._r8/spec_grnd_cond(s)+dz_soisno(j)/(diffus(j)*2._r8))
					! layer resistance mult. by k_h_cc for dp1_zp1 term
					dp1_zp1(j) = 2._r8/(dz_soisno(j)*k_h_cc(j,s)/diffus(j)+dz_soisno(j+1)/diffus(j+1))
				else if (j == 1) then ! water table at surface: multiply ground resistance by k_h_cc
					dm1_zm1(j) = 1._r8/(k_h_cc(j-1,s)/spec_grnd_cond(s)+dz_soisno(j)/(diffus(j)*2._r8))
					! air concentration will be mult. by k_h_cc below
					dp1_zp1(j) = 2._r8/(dz_soisno(j)/diffus(j)+dz_soisno(j+1)/diffus(j+1))
				else if (j <= nl_soil-1 .and. j /= jwt .and. j /= jwt+1) then
					dm1_zm1(j) = 2._r8/(dz_soisno(j)/diffus(j)+dz_soisno(j-1)/diffus(j-1))
					dp1_zp1(j) = 2._r8/(dz_soisno(j)/diffus(j)+dz_soisno(j+1)/diffus(j+1))
				else if (j <= nl_soil-1 .and. j == jwt) then ! layer resistance mult. by k_h_cc for dp1_zp1 term
					dm1_zm1(j) = 2._r8/(dz_soisno(j)/diffus(j)+dz_soisno(j-1)/diffus(j-1))
					dp1_zp1(j) = 2._r8/(dz_soisno(j)*k_h_cc(j,s)/diffus(j)+dz_soisno(j+1)/diffus(j+1))
					! Concentration in layer will be mult. by k_h_cc below
				else if (j <= nl_soil-1) then ! j==jwt+1: layer above resistance mult. by k_h_cc for dm1_zm1 term
					dm1_zm1(j) = 2._r8/(dz_soisno(j)/diffus(j)+dz_soisno(j-1)*k_h_cc(j-1,s)/diffus(j-1))
					! Concentration in layer above will be mult. by k_h_cc below
					dp1_zp1(j) = 2._r8/(dz_soisno(j)/diffus(j)+dz_soisno(j+1)/diffus(j+1))
				else if (j /= jwt+1) then ! j ==nl_soil
					dm1_zm1(j) = 2._r8/(dz_soisno(j)/diffus(j)+dz_soisno(j-1)/diffus(j-1))
				else                    ! jwt == nl_soil-1: layer above resistance mult. by k_h_cc for dm1_zm1 term
					dm1_zm1(j) = 2._r8/(dz_soisno(j)/diffus(j)+dz_soisno(j-1)*k_h_cc(j-1,s)/diffus(j-1))
				end if
				Call print_var(dm1_zm1(j) ,'ch4_tran dm1_zm1(j)',idate,j)
				Call print_var(dp1_zp1(j) ,'ch4_tran dp1_zp1(j)',idate,j)
			end do ! j; nl_soil
	
			! Perform a second loop for the tridiagonal coefficients since need dp1_zp1 and dm1_z1 at each depth
			do j = 0,nl_soil
				conc_ch4_rel_old(j) = conc_ch4_rel(j)
	
				if (j > 0) dzj = dz_soisno(j)
				if (j == 0) then ! top layer (atmosphere) doesn't change regardless of where WT is
					at(j) = 0._r8
					bt(j) = 1._r8
					ct(j) = 0._r8
					rt(j) = c_atm(s) ! 0th level stays at constant atmospheric conc
				elseif (j < nl_soil .and. j == jwt) then ! concentration inside needs to be mult. by k_h_cc for dp1_zp1 term
					at(j) = -0.5_r8 / dzj * dm1_zm1(j)
					bt(j) = epsilon_t(j,s) / deltim + 0.5_r8 / dzj * (dp1_zp1(j)*k_h_cc(j,s) + dm1_zm1(j))
					ct(j) = -0.5_r8 / dzj * dp1_zp1(j)
				elseif (j < nl_soil .and. j == jwt+1) then
					! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
					at(j) = -0.5_r8 / dzj * dm1_zm1(j) * k_h_cc(j-1,s)
					bt(j) = epsilon_t(j,s) / deltim + 0.5_r8 / dzj * (dp1_zp1(j) + dm1_zm1(j))
					ct(j) = -0.5_r8 / dzj * dp1_zp1(j)
				elseif (j < nl_soil) then
					at(j) = -0.5_r8 / dzj * dm1_zm1(j)
					bt(j) = epsilon_t(j,s) / deltim + 0.5_r8 / dzj * (dp1_zp1(j) + dm1_zm1(j))
					ct(j) = -0.5_r8 / dzj * dp1_zp1(j)
				else if (j == nl_soil .and. j== jwt+1) then
					! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
					at(j) = -0.5_r8 / dzj * dm1_zm1(j) * k_h_cc(j-1,s)
					bt(j) = epsilon_t(j,s) / deltim + 0.5_r8 / dzj * dm1_zm1(j)
					ct(j) = 0._r8
				else ! j==nl_soil and jwt<nl_soil-1 or jwt==nl_soil: 0 flux at bottom
					at(j) = -0.5_r8 / dzj * dm1_zm1(j)
					bt(j) = epsilon_t(j,s) / deltim + 0.5_r8 / dzj * dm1_zm1(j)
					ct(j) = 0._r8
				endif
			enddo ! j; nl_soil
	
	
			jtop = 0
	
	
			if (s == 1) then  ! CH4
	
				! Set rt, since it depends on conc
				do j = 1,nl_soil
	
					! For correct balance, deprecate source_old.
					! The source terms are effectively constant over the timestep.
					source_old(j,s) = source(j,s)
					! source_old could be removed later
					epsilon_t_old(j,s) = epsilon_t(j,s)
					! epsilon_t acts like source also
					dzj = dz_soisno(j)
					if (j < nl_soil .and. j == jwt) then ! concentration inside needs to be mult. by k_h_cc for dp1_zp1 term
						rt(j) = epsilon_t_old(j,s) / deltim * conc_ch4_rel(j) +           &
							0.5_r8 / dzj * (dp1_zp1(j) * (conc_ch4_rel(j+1)-conc_ch4_rel(j)*k_h_cc(j,s)) - &
							dm1_zm1(j) * (conc_ch4_rel(j)  -conc_ch4_rel(j-1))) + &
							0.5_r8 * (source(j,s) + source_old(j,s))
					elseif (j < nl_soil .and. j == jwt+1) then
						! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
						rt(j) = epsilon_t_old(j,s) / deltim * conc_ch4_rel(j) +           &
							0.5_r8 / dzj * (dp1_zp1(j) * (conc_ch4_rel(j+1)-conc_ch4_rel(j)) - &
							dm1_zm1(j) * (conc_ch4_rel(j) -conc_ch4_rel(j-1)*k_h_cc(j-1,s))) + &
							0.5_r8 * (source(j,s) + source_old(j,s))
					elseif (j < nl_soil) then
						rt(j) = epsilon_t_old(j,s) / deltim * conc_ch4_rel(j) +           &
							0.5_r8 / dzj * (dp1_zp1(j) * (conc_ch4_rel(j+1)-conc_ch4_rel(j)) - &
							dm1_zm1(j) * (conc_ch4_rel(j)  -conc_ch4_rel(j-1))) + &
							0.5_r8 * (source(j,s) + source_old(j,s))
					else if (j == nl_soil .and. j== jwt+1) then
						! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
						rt(j) = epsilon_t_old(j,s) / deltim * conc_ch4_rel(j) +           &
							0.5_r8 / dzj * ( - dm1_zm1(j) * (conc_ch4_rel(j) -conc_ch4_rel(j-1)*k_h_cc(j-1,s))) + &
							0.5_r8 * (source(j,s) + source_old(j,s))
					else  !j==nl_soil
						rt(j) = epsilon_t_old(j,s) / deltim * conc_ch4_rel(j) +           &
							0.5_r8 / dzj * ( - dm1_zm1(j) * (conc_ch4_rel(j)  -conc_ch4_rel(j-1))) + &
							0.5_r8 * (source(j,s) + source_old(j,s))
					endif
					epsilon_t_old(j,s) = epsilon_t(j,s)
					source_old(j,s) = source(j,s)
					Call print_var(at(j) ,'ch4_tran at(j)',idate,j,s)
					Call print_var(bt(j) ,'ch4_tran bt(j)',idate,j,s)
					Call print_var(ct(j) ,'ch4_tran ct(j)',idate,j,s)
					Call print_var(rt(j) ,'ch4_tran rt(j)',idate,j,s)
					Call print_var(epsilon_t_old(j,s) ,'ch4_tran epsilon_t_old(j,s)',idate,j,s)
					Call print_var(source_old(j,s) ,'ch4_tran source_old(j,s)',idate,j,s)
	
				enddo ! j; nl_soil
	
				call Tridiagonal(0, nl_soil, &
					jtop, &
					at(:), &
					bt(:), &
					ct(:), &
					rt(:), &
					conc_ch4_rel(0:nl_soil))
	
					
	
				! Calculate net ch4 flux to the atmosphere from the surface (+ to atm)
				if (jwt /= 0) then ! WT not at the surface
					ch4_surf_diff = dm1_zm1(1) * ( (conc_ch4_rel(1)+conc_ch4_rel_old(1))/2._r8 &
						- c_atm(s)) ! [mol/m2/s]
					ch4_surf_ebul = 0._r8 ! all the ebullition has already come out in the soil column (added to source)
					! Try adding directly to atm. to prevent destabilization of diffusion
					!ch4_surf_ebul = ch4_ebul_tot ! [mol/m2/s]
				else ! WT at the surface; i.e., jwt==0
					ch4_surf_diff = dm1_zm1(1) * ( (conc_ch4_rel(1)+conc_ch4_rel_old(1))/2._r8 &
						- c_atm(s)*k_h_cc(0,s)) ! [mol/m2/s]
					! atmospheric concentration gets mult. by k_h_cc as above
					ch4_surf_ebul = ch4_ebul_tot ! [mol/m2/s]

					err1 = (conc_ch4_rel(1)+conc_ch4_rel_old(1))/2._r8-c_atm(s)*k_h_cc(0,s)
					err2 = (conc_ch4_rel(1)+conc_ch4_rel_old(1))/2._r8
					err3 = c_atm(s)*k_h_cc(0,s)
					Call print_var(err1,'ch4_tran err1',idate,1)
					Call print_var(err2,'ch4_tran err2',idate,1)
					Call print_var(err3,'ch4_tran err3',idate,1)
					Call print_var(conc_ch4_rel(1),'ch4_tran conc_ch4_rel(1)',idate,1)
					Call print_var(conc_ch4_rel_old(1),'ch4_tran conc_ch4_rel_old(1)',idate,1)
					Call print_var(c_atm(s),'ch4_tran c_atm(s)',idate,1,s)
					Call print_var(k_h_cc(0,s),'ch4_tran k_h_cc(0,s)',idate,1,s)
					Call print_var(ch4_surf_diff,'ch4_tran ch4_surf_diff',idate,1)
					Call print_var(ch4_surf_ebul,'ch4_tran ch4_surf_ebul',idate,1)
				endif
	
				! Ensure that concentrations stay above 0
				! This should be done after the flux, so that the flux calculation is consistent.
				do j = 1,nl_soil
					if (conc_ch4_rel(j) < 0._r8) then
						deficit = - conc_ch4_rel(j)*epsilon_t(j,1)*dz_soisno(j)  ! Mol/m^2 added
						if (deficit > 1.e-3_r8 * DEF_CH4%scale_factor_gasdiff) then
							if (deficit > 1.e-2_r8) then
								write(6,*)'Note: sink > source in ch4_tran, sources are changing '// &
										' quickly relative to diffusion timestep, and/or diffusion is rapid.'
										write(6,*)'Lat,Lon=',lat,lon
								write(6,*)'This typically occurs when there is a larger than normal '// &
										' diffusive flux.'
								write(6,*)'If this occurs frequently, consider reducing land model (or '// &
										' methane model) timestep, or reducing the max. sink per timestep in the methane model.'
							end if
							write(6,*) 'Negative conc. in ch4tran. j,deficit (mol):',j,deficit
						end if
						conc_ch4_rel(j) = 0._r8
						! Subtract deficit
						ch4_surf_diff = ch4_surf_diff - deficit/deltim
					end if
				enddo
  
  
		   	elseif (s == 2) then  ! O2
  
				! Set rt, since it depends on conc
				do j = 1,nl_soil
	
					! For correct balance, deprecate source_old.
					source_old(j,s) = source(j,s)
					! source_old could be removed later
					epsilon_t_old(j,s) = epsilon_t(j,s)
					! epsilon_t acts like source also
					dzj     = dz_soisno(j)
					if (j < nl_soil .and. j == jwt) then ! concentration inside needs to be mult. by k_h_cc for dp1_zp1 term
						rt(j) = epsilon_t_old(j,s) / deltim * conc_o2_rel(j) +           &
							0.5_r8 / dzj * (dp1_zp1(j) * (conc_o2_rel(j+1)-conc_o2_rel(j)*k_h_cc(j,s)) - &
							dm1_zm1(j) * (conc_o2_rel(j)  -conc_o2_rel(j-1))) + &
							0.5_r8 * (source(j,s) + source_old(j,s))
					elseif (j < nl_soil .and. j == jwt+1) then
						! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
						rt(j) = epsilon_t_old(j,s) / deltim * conc_o2_rel(j) +           &
							0.5_r8 / dzj * (dp1_zp1(j) * (conc_o2_rel(j+1)-conc_o2_rel(j)) - &
							dm1_zm1(j) * (conc_o2_rel(j) -conc_o2_rel(j-1)*k_h_cc(j-1,s))) + &
							0.5_r8 * (source(j,s) + source_old(j,s))
					elseif (j < nl_soil) then
						rt(j) = epsilon_t_old(j,s) / deltim * conc_o2_rel(j) +           &
							0.5_r8 / dzj * (dp1_zp1(j) * (conc_o2_rel(j+1)-conc_o2_rel(j)) - &
							dm1_zm1(j) * (conc_o2_rel(j)  -conc_o2_rel(j-1))) + &
							0.5_r8 * (source(j,s) + source_old(j,s))
					else if (j == nl_soil .and. j== jwt+1) then
						! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
						rt(j) = epsilon_t_old(j,s) / deltim * conc_o2_rel(j) +           &
							0.5_r8 / dzj * ( - dm1_zm1(j) * (conc_o2_rel(j) -conc_o2_rel(j-1)*k_h_cc(j-1,s))) + &
							0.5_r8 * (source(j,s) + source_old(j,s))
					else  !j==nl_soil
						rt(j) = epsilon_t_old(j,s) / deltim * conc_o2_rel(j) +           &
							0.5_r8 / dzj * ( - dm1_zm1(j) * (conc_o2_rel(j)  -conc_o2_rel(j-1))) + &
							0.5_r8 * (source(j,s) + source_old(j,s))
					endif
					epsilon_t_old(j,s) = epsilon_t(j,s)
					source_old(j,s) = source(j,s)
					Call print_var(rt(j) ,'ch4_tran rt(j)',idate,j)
					Call print_var(epsilon_t_old(j,s) ,'ch4_tran epsilon_t_old(j,s)',idate,j,s)
					Call print_var(source_old(j,s) ,'ch4_tran source_old(j,s)',idate,j,s)
				enddo ! j; nl_soil
  
				call Tridiagonal(0, nl_soil, jtop, &
					at(:), &
					bt(:), &
					ct(:), &
					rt(:), &
					conc_o2_rel(0:nl_soil))
  
				! Ensure that concentrations stay above 0
				do j = 1,nl_soil
					Call print_var(conc_o2_rel(j) ,'ch4_tran 1 conc_o2_rel(j)',idate,j)

					conc_o2_rel(j) = max (conc_o2_rel(j), 1.e-12_r8)
					Call print_var(conc_o2_rel(j) ,'ch4_tran 2 conc_o2_rel(j)',idate,j)

					! In case of pathologically large aerenchyma conductance. Should be OK in general but
					! this will maintain stability even if a PATCH with very small weight somehow has an absurd NPP or LAI.
					! Also, oxygen above ambient will probably bubble.
					conc_o2_rel(j) = min (conc_o2_rel(j), c_atm(2)/epsilon_t(j,2))
					Call print_var(conc_o2_rel(j) ,'ch4_tran 3 conc_o2_rel(j)',idate,j)
				enddo
  
		   endif  ! species
  
		enddo  ! species
  
		! Update absolute concentrations per unit volume
		do j = 1,nl_soil ! No need to update the atm. level concentrations
			conc_ch4(j) = conc_ch4_rel(j)*epsilon_t(j,1)
			conc_o2(j)  = conc_o2_rel(j) *epsilon_t(j,2)
			Call print_var(conc_ch4(j) ,'ch4_tran 2 conc_ch4(j)',idate,j)
			Call print_var(conc_o2(j) ,'ch4_tran 2 conc_o2(j)',idate,j)
			Call print_var(conc_ch4_rel(j) ,'ch4_tran conc_ch4_rel(j)',idate,j)
			Call print_var(conc_o2_rel(j) ,'ch4_tran conc_o2_rel(j)',idate,j)
		end do
  
		! Do Balance Check and absorb small
		!    discrepancy into surface flux.
		err5 = 0
		err6 = 0
		err7 = 0
		do j = 1,nl_soil
		   	if (j == 1) errch4 = 0._r8
			errch4 = errch4 + (conc_ch4(j) - conc_ch4_bef(j))*dz_soisno(j)
			errch4 = errch4 - ch4_prod_depth(j)*dz_soisno(j)*deltim
			errch4 = errch4 + ch4_oxid_depth(j)*dz_soisno(j)*deltim

			err1 = (conc_ch4(j) - conc_ch4_bef(j))*dz_soisno(j)
			err2 = ch4_prod_depth(j)*dz_soisno(j)*deltim
			err3 = ch4_oxid_depth(j)*dz_soisno(j)*deltim
			err4 = err1 - err2 + err3

			err5 = err1 + err5
			err6 = err2 + err6
			err7 = err3 + err7
			Call print_var(err1,'ch4_tran err1',idate,j,s,.True.)
			Call print_var(err2,'ch4_tran err2',idate,j,s,.True.)
			Call print_var(err3,'ch4_tran err3',idate,j,s,.True.)
			Call print_var(err4,'ch4_tran err4',idate,j,s,.True.)

			Call print_var(err5,'ch4_tran err5',idate,j,s,.True.)
			Call print_var(err6,'ch4_tran err6',idate,j,s,.True.)
			Call print_var(err7,'ch4_tran err7',idate,j,s,.True.)

			Call print_var(conc_ch4(j),'ch4_tran conc_ch4(j)',idate,j,s,.True.)
			Call print_var(conc_ch4_bef(j),'ch4_tran conc_ch4_bef(j)',idate,j,s,.True.)
			Call print_var(ch4_prod_depth(j),'ch4_tran ch4_prod_depth(j)',idate,j,s,.True.)
			Call print_var(ch4_oxid_depth(j),'ch4_tran ch4_oxid_depth(j)',idate,j,s,.True.)
			Call print_var(errch4,'ch4_tran errch4',idate,j,s,.True.)
		end do
		
		! For history make sure that grnd_ch4_cond includes snow, for methane diffusivity
		grnd_ch4_cond = spec_grnd_cond(1)
  
		errch4 = errch4 + (ch4_surf_aere + ch4_surf_ebul + ch4_surf_diff)*deltim
		err8 = (ch4_surf_aere + ch4_surf_ebul + ch4_surf_diff)*deltim
		Call print_var(err8,'ch4_tran err8',idate)

		Call print_var(ch4_surf_aere,'ch4_tran ch4_surf_aere',idate)
		Call print_var(ch4_surf_ebul,'ch4_tran ch4_surf_ebul',idate)
		Call print_var(ch4_surf_diff,'ch4_tran ch4_surf_diff',idate)
		Call print_var(errch4,'ch4_tran errch4',idate)
		Call print_var(lat,'ch4_tran lat',idate)

		if (abs(errch4) < 1.e-6_r8) then 
		   ch4_surf_diff = ch4_surf_diff - errch4/deltim
		else ! errch4 > 1e-8 mol / m^2 / timestep
		   	write(6,*)'CH4 Conservation Error in CH4Mod during diffusion, istep, errch4 (mol /m^2.timestep)', &
			idate,errch4
			write(6,*)'Lat,Lon=',lat,lon
			CALL CoLM_stop ()
		end if
  
	end subroutine ch4_tran

	!---------------------------------------------------------------------------
	subroutine Tridiagonal (lbj, ubj, jtop, a, b, c, r, u)
		!-----------------------------------------------------------------------
		! !DESCRIPTION:
		! Tridiagonal matrix solution
		!-----------------------------------------------------------------------

		!-----------------------Argument---------- -----------------------------
		implicit none         
		integer , intent(in)    :: lbj, ubj   ! lbinning and ubing level indices
		integer , intent(in)    :: jtop       ! top level for each column [col]
		real(r8), intent(in)    :: a(lbj:ubj) ! "a" left off diagonal of tridiagonal matrix [col, j]
		real(r8), intent(in)    :: b(lbj:ubj) ! "b" diagonal column for tridiagonal matrix [col, j]
		real(r8), intent(in)    :: c(lbj:ubj) ! "c" right off diagonal tridiagonal matrix [col, j]
		real(r8), intent(in)    :: r(lbj:ubj) ! "r" forcing term of tridiagonal matrix [col, j]
		real(r8), intent(inout) :: u(lbj:ubj) ! solution [col, j]

		!-----------------------Local Variables---------------------------------
		integer  :: j                 !indices
		real(r8) :: gam(lbj:ubj)      !temporary
		real(r8) :: bet               !temporary

		!-----------------------------------------------------------------------
		! Solve the matrix

		bet = b(jtop)

		do j = lbj, ubj
			if (j >= jtop) then
				if (j == jtop) then
					u(j) = r(j) / bet
				else
					gam(j) = c(j-1) / bet
					bet = b(j) - a(j) * gam(j)
					u(j) = (r(j) - a(j)*u(j-1)) / bet
				end if
			end if
		end do

		do j = ubj-1,lbj,-1
			if (j >= jtop) then
				u(j) = u(j) - gam(j+1) * u(j+1)
			end if
		end do
		
	end subroutine Tridiagonal

	subroutine henry_law (t_grnd,t_soisno,k_h_cc)
		real(r8), intent(in) :: &
			t_grnd                 		   , &! ground surface temperature [k]
			t_soisno (maxsnl+1:nl_soil)       ! soil temperature [K]

		real(r8), intent(out) :: &
			k_h_cc(0:nl_soil,ngases)          ! Dimensionless Henry's coefficient [-]

		integer :: j,s  			          ! Indices
		real(r8):: k_h						  ! Henry's coefficient [mol/L/atm]

		do j = 0,nl_soil
			do s=1,2
				if (j==0) then
					k_h = kh_theta(s)*exp(c_h(s) * (1._r8 / t_grnd - 1._r8 / kh_tbase))
					! [mol/L/atm] = [mol/L/atm]*e**([K]*[1/K])
					k_h_cc(j,s) = k_h * rgasLatm * t_grnd
					! [-]  = [mol/L/atm]*[L*atm/mol/K]*[K]
					call print_var(k_h_cc(j,s),'ch4 k_h_cc(j,s)',idate,j,s)
				else
					k_h = kh_theta(s)*exp(c_h(s) * (1._r8 / t_soisno(j) - 1._r8 / kh_tbase))
					! [mol/L/atm] = [mol/L/atm]*e**([K]*[1/K])
					k_h_cc(j,s) = k_h * rgasLatm * t_soisno(j)
					! [-]  = [mol/L/atm]*[L*atm/mol/K]*[K]
					call print_var(k_h_cc(j,s),'ch4 k_h_cc(j,s)',idate,j,s)
				end if
			end do
		end do
	end subroutine henry_law

	subroutine split_ch4_o2_phases( dz_soisno, wliq_soisno, porsl, &
									conc_ch4, conc_o2, k_h_cc, idate,&
									vol_aqu,vol_gas,f_aqu,f_gas,&
									conc_ch4_gas,conc_ch4_aqu,conc_ch4_porsl,conc_ch4_gas_porsl,conc_ch4_aqu_porsl,&
									conc_o2_gas,conc_o2_aqu,conc_o2_porsl,conc_o2_gas_porsl,conc_o2_aqu_porsl)

		implicit none
		real(r8), intent(in) :: dz_soisno(1:nl_soil)     ! layer thickness [m]
		real(r8), intent(in) :: wliq_soisno(1:nl_soil)   ! liquid water in layers [kg/m2]
		real(r8), intent(in) :: porsl(1:nl_soil)                ! volumetric soil water at saturation (porosity)
		real(r8), intent(in) :: conc_ch4(1:nl_soil)          ! CH4 concentration in each soil layer [mol/m3]
		real(r8), intent(in) :: conc_o2(1:nl_soil)           ! O2 concentration in each soil layer [mol/m3]
		real(r8), intent(in) :: k_h_cc(0:nl_soil,ngases)        ! ratio of mol/m3 in liquid to mol/m3 in gas [-]
		integer, intent(in) :: idate(3)

		integer :: j,s
		real(r8), intent(out) :: vol_aqu(1:nl_soil)          ! liquid volumetric water content [m3/m3]
		real(r8), intent(out) :: vol_gas(1:nl_soil)          ! air volumetric water content [m3/m3]
		real(r8), intent(out) :: f_aqu(1:nl_soil)            ! water-filled proportion [-]
		real(r8), intent(out) :: f_gas(1:nl_soil)            ! air-filled proportion [-]
		real(r8), intent(out) :: conc_ch4_gas(1:nl_soil)     ! gas phase CH4 conc [mol/m3]
		real(r8), intent(out) :: conc_ch4_aqu(1:nl_soil)     ! aqueous phase CH4 conc [mol/m3]
		real(r8), intent(out) :: conc_ch4_porsl(1:nl_soil)   ! CH4 conc per porosity [mol/m3]
		real(r8), intent(out) :: conc_ch4_gas_porsl(1:nl_soil) ! gas phase CH4 conc per porosity [mol/m3]
		real(r8), intent(out) :: conc_ch4_aqu_porsl(1:nl_soil) ! aqueous phase CH4 conc per porosity [mol/m3]
		real(r8), intent(out) :: conc_o2_gas(1:nl_soil)      ! gas phase O2 conc [mol/m3]
		real(r8), intent(out) :: conc_o2_aqu(1:nl_soil)      ! aqueous phase O2 conc [mol/m3]
		real(r8), intent(out) :: conc_o2_porsl(1:nl_soil)    ! O2 conc per porosity [mol/m3]
		real(r8), intent(out) :: conc_o2_gas_porsl(1:nl_soil)! gas phase O2 conc per porosity [mol/m3]
		real(r8), intent(out) :: conc_o2_aqu_porsl(1:nl_soil)! aqueous phase O2 conc per porosity [mol/m3]

		!-----------------------------------------------------------
		! Main computation: calculate phase partitioning for CH4 and O2
		!-----------------------------------------------------------
		do j = 1, nl_soil 

			! ---- Compute volumetric water content ----
			vol_aqu(j) = min(wliq_soisno(j)/(dz_soisno(j)*denh2o), porsl(j))
			! [m3/m3] = [kg/m2] / ([m] * [kg/m3])

			! ---- Compute volumetric gas content ----
			vol_gas(j) = max(porsl(j) - vol_aqu(j), 0._r8)

			! ---- Compute filled proportions ----
			f_aqu(j) = vol_aqu(j)/porsl(j)
			f_gas(j) = vol_gas(j)/porsl(j)

			! ---- CH4 partitioning between gas and aqueous phases ----
			conc_ch4_aqu(j) = conc_ch4(j)/(f_aqu(j)+f_gas(j)/k_h_cc(j,1)) 
			conc_ch4_gas(j) = conc_ch4(j)/(k_h_cc(j,1)*f_aqu(j)+f_gas(j)) 

			! ---- O2 partitioning between gas and aqueous phases ----
			conc_o2_aqu(j) = conc_o2(j)/(f_aqu(j)+f_gas(j)/k_h_cc(j,1)) 
			conc_o2_gas(j) = conc_o2(j)/(k_h_cc(j,1)*f_aqu(j)+f_gas(j)) 

			! ---- Concentrations normalized by porosity ----
			conc_ch4_porsl(j)      = conc_ch4(j)/porsl(j)
			conc_ch4_aqu_porsl(j)  = conc_ch4_aqu(j)/porsl(j)
			conc_ch4_gas_porsl(j)  = conc_ch4_gas(j)/porsl(j)
			conc_o2_porsl(j)       = conc_o2(j)/porsl(j)
			conc_o2_aqu_porsl(j)   = conc_o2_aqu(j)/porsl(j)
			conc_o2_gas_porsl(j)   = conc_o2_gas(j)/porsl(j)

			! ---- Print results for debugging ----
			call print_var(vol_aqu(j),'ch4 vol_aqu(j)',idate,j)
			call print_var(vol_gas(j),'ch4 vol_gas(j)',idate,j)
			call print_var(f_aqu(j),'ch4 f_aqu(j)',idate,j)
			call print_var(f_gas(j),'ch4 f_gas(j)',idate,j)
			call print_var(conc_ch4(j),'ch4 1 conc_ch4(j)',idate,j)
			call print_var(conc_ch4_aqu(j),'ch4 conc_ch4_aqu(j)',idate,j)
			call print_var(conc_ch4_gas(j),'ch4 conc_ch4_gas(j)',idate,j)
			call print_var(conc_ch4_aqu_porsl(j),'ch4 conc_ch4_aqu_porsl(j)',idate,j)
			call print_var(conc_ch4_gas_porsl(j),'ch4 conc_ch4_gas_porsl(j)',idate,j)

		end do
	end subroutine split_ch4_o2_phases

#ifdef SinglePoint
	subroutine print_var_real8(var1,varname1,idate,j,s,layer_ok)
		real(r8), intent(in) :: var1
		character(len=*), intent(in) :: varname1
		integer, intent(in) :: idate(3)
		integer, intent(in), optional :: j,s
		logical, intent(in), optional :: layer_ok
		logical :: print_all_layers
		
		! Determine if we should print all layers
		if (present(layer_ok)) then
			print_all_layers = layer_ok
		else
			print_all_layers = .false.
		endif
		
		if (idate(2)==3) then
			! Print if: j not present, OR j==0 or j==1, OR layer_ok is true
			if (.not. present(j) .or. (j==0 .or. j==1 .or. j==10)) then
				print*, "===================================================================================="
				print*, "year, day, second",idate(1),idate(2),idate(3)
				if (present(j)) then
					print*, "layer is ",j
				endif
				if (present(s)) then
					print*, "ngase is ",s
				endif
				print*, trim(varname1)," = ",var1
				print*, "===================================================================================="
			endif
		endif
	end subroutine print_var_real8

	subroutine print_var_int32(var1,varname1,idate,j,s,layer_ok)
		integer, intent(in) :: var1
		character(len=*), intent(in) :: varname1
		integer, intent(in) :: idate(3)
		integer, intent(in), optional :: j,s
		logical, intent(in), optional :: layer_ok
		logical :: print_all_layers
		
		! Determine if we should print all layers
		if (present(layer_ok)) then
			print_all_layers = layer_ok
		else
			print_all_layers = .false.
		endif
		
		if (idate(2)==3) then
			! Print if: j not present, OR j==0 or j==1, OR layer_ok is true
			if (.not. present(j) .or. (j==0 .or. j==1 .or. j==10)) then
				print*, "===================================================================================="
				print*, "year, day, second",idate(1),idate(2),idate(3)
				if (present(j)) then
					print*, "layer is ",j
				endif
				if (present(s)) then
					print*, "ngase is ",s
				endif
				print*, trim(varname1)," = ",var1
				print*, "===================================================================================="
			endif
		endif
	end subroutine print_var_int32
#else
	subroutine print_var_real8(var1,varname1,idate,j,s,layer_ok)
		real(r8), intent(in) :: var1
		character(len=*), intent(in) :: varname1
		integer, intent(in) :: idate(3)
		integer, intent(in), optional :: j,s
		logical, intent(in), optional :: layer_ok
		! Do nothing
	end subroutine print_var_real8

	subroutine print_var_int32(var1,varname1,idate,j,s,layer_ok)
		integer, intent(in) :: var1
		character(len=*), intent(in) :: varname1
		integer, intent(in) :: idate(3)
		integer, intent(in), optional :: j,s
		logical, intent(in), optional :: layer_ok
		! Do nothing
	end subroutine print_var_int32
#endif

END MODULE MOD_ch4
! --------- EOP ----------