#include <define.h>
#if (defined TRACER) && (defined BGC)

module MOD_Tracer_Reactive_Methane_Physics
    !=======================================================================
	! DESCRIPTION:
	! Module holding routines to calculate methane fluxes
	! The driver averages up to gridcell, weighting by finundated, and checks for balance errors.
	! Sources, sinks, "competition" for CH4 & O2, & transport are resolved in methane_tran.

	! ORIGINAL:
	! The Community Land Model version 5.0 (CLM5.0)

	! REFERENCES:
	! Lawrence, D.M., Fisher, R.A., Koven, C.D., Oleson, K.W., Swenson, S.C., Bonan, G., Collier, N.,
	! Ghimire, B., van Kampenhout, L., Kennedy, D. and Kluzek, E., 2019.
	! The Community Land Model version 5: Description of new features, benchmarking,
	! and impact of forcing uncertainty. Journal of Advances in Modeling Earth Systems, 11(12), 4245-4287.

	! REVISION:
	! Xionghui Xu, 2025, 1) Modify original CLM5 to be compatible with CoLM code structure.
    !                    2) Adapt and stabilize the process coupling for CoLM202X.
	!=======================================================================
	use MOD_Precision
#ifdef SinglePoint
   USE MOD_SingleSrfdata, only: SITE_WTD_valid
#endif
	use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
	use MOD_SPMD_Task
	use MOD_TimeManager
   use MOD_Vars_TimeInvariants, only: wetwatmax, patchlatr
	use MOD_Vars_Global, only : maxsnl,nl_soil,nl_lake,spval,PI,deg2rad
	use MOD_Const_Physical, only: denh2o, denice, tfrz, grav, vonkar
	use MOD_Tracer_Reactive_Methane_Const
   USE MOD_Namelist, only: DEF_wetland_finundation_scheme, DEF_USE_Dynamic_Wetland, &
                           USE_SITE_WTD
	USE MOD_Tracer_Reactive_Methane_GIEMS, only: giems_finundated, giems_active
   USE MOD_FloodInfiltration, only: fld_cap_p
   USE MOD_Tracer_Reactive_Methane_BgcLink, only: is_paddy_rice_live, &
      rice_days_past_planting, rice_days_since_harvest, &
      BIOME_TROPICAL_PEAT, BIOME_TROPICAL_FLOODPLAIN, BIOME_TEMPERATE_MARSH, &
      BIOME_BOREAL_FEN, BIOME_BOREAL_BOG
	USE MOD_Tracer_Reactive_Methane_VegOverride, only: get_aere_poros, get_aere_tillerC, &
	                                          get_aere_radius, get_aere_scale
   USE MOD_Tracer_Reactive_Methane_State, only: annavg_tsoi, tempavg_tsoi_int, tempavg_tsoi_cnt
   USE MOD_Tracer_Reactive_Methane_State, only: lake_icebubble_ch4_stock
   USE MOD_Tracer_Reactive_Methane_State, only: lake_cact
   USE MOD_Tracer_Reactive_Methane_State, only: methanogen_act_unsat, methanogen_act_sat
   USE MOD_Tracer_Reactive_Methane_State, only: acceptor_unsat, acceptor_sat
   USE MOD_Tracer_Reactive_Methane_WetlandVeg, only: wetveg_forest, wetveg_ombro, wetveg_bog, &
      wetveg_emerg, wetveg_dome, wetveg_moss
		USE MOD_Tracer_Reactive_Methane_State, only: f_inund_levee_patch, f_inund_flood_patch, &
      f_inund_flood_depth_patch, wetland_frac_per_patch, methane_wetland_type, &
		                                     biome_f_methane_patch, biome_redoxlag_patch, &
		                                     methane_finundated, methane_soil_finundated, &
		                                     methane_soil_zwt, methane_surf_flux_wetland, &
		                                     methane_surf_flux_soil, methane_surf_flux_lake, &
      methane_surf_flux_rice, &
      methane_prod_tot_wetland, methane_oxid_tot_wetland, &
      methane_surf_aere_wetland, methane_surf_ebul_wetland, methane_surf_diff_wetland, &
      methane_area_wetland, methane_area_soil, methane_area_rice, methane_area_lake, &
		                                     methane_surf_aere_soil, &
		                                     methane_surf_aere_rice, methane_surf_ebul_soil, &
		                                     methane_surf_ebul_rice, methane_surf_diff_soil, &
		                                     methane_surf_diff_rice, methane_prod_tot_soil, &
		                                     methane_prod_tot_rice, methane_oxid_tot_soil, &
		                                     methane_oxid_tot_rice
	!-----------------------------------------------------------------------
	implicit none
	save

	public  :: methane
	public  :: methane_host_water_reset
	public  :: methane_host_water_report
   public  :: wetland_max_wtd_of
   public  :: wetland_inflow_share
   public  :: wetland_bog_share
   public  :: wetland_peat_outflow
   public  :: wetland_moss_rss
   public  :: wetland_is_temperate_marsh
   public  :: wetland_emergent_share
   public  :: methane_tsoi_annual

   ! peatland share of each patch's wetland tile for wetland_moss_rss
   ! and wetland_f_tropical_peat, read at start-up by
   ! wetland_peat_share_read
   real(r8), allocatable, public :: wetland_tile_peat_share(:)

	private :: methane_annualupdate
	private :: methane_prod
   private :: methane_acceptor
   private :: methane_root_oxid
	private :: methane_oxid
	private :: methane_aere
	private :: methane_ebul
	private :: methane_tran
	private

	! Worst host-water/inundation disagreement seen so far, as a fraction of
	! pore volume, and whether the one-shot warning has been emitted. Module
	! state rather than an argument because the report is per run, not per
	! patch; the driver's patch loop is not OpenMP-parallel (CoLMDRIVER.F90),
	! so a plain saved scalar is safe here.
	real(r8) :: host_water_max_resid = 0._r8
	real(r8) :: host_water_reported  = 0._r8
	logical  :: host_water_warned    = .false.

contains

	!-----------------------------------------------------------------------
	!> Clear the host-water diagnostic. Called when methane state is
	!! allocated, so a model re-initialised in the same process does not
	!! inherit the previous run's worst residual or its already-warned flag.
	subroutine methane_host_water_reset ()
		host_water_max_resid = 0._r8
		host_water_reported  = 0._r8
		host_water_warned    = .false.
	end subroutine methane_host_water_reset

	!> Report the worst host-water disagreement across all worker ranks.
	!!
	!! The running maximum and the one-shot warning are per-rank state, so
	!! without this every rank warns independently and the global worst case is
	!! never printed.
	!!
	!! NOTE this runs from tracer_report, i.e. once per CoLMDRIVER call, not
	!! once per run -- an earlier comment claimed otherwise. The reduction is
	!! one double and cheap enough per timestep, but the summary line must not
	!! be printed every step, so worker 0 prints only when the global worst
	!! case grows past what it last reported.
	subroutine methane_host_water_report ()

		USE MOD_SPMD_Task
		USE MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE
		implicit none
		real(r8) :: gmax

		! Reduce over the WORKER communicator, not p_comm_glb.  This runs
		! from tracer_report, which CoLM.F90 reaches inside IF (p_is_worker):
		! master and IO ranks never arrive, so a p_comm_glb collective here
		! deadlocks every MPI run with methane enabled.  Report from worker 0
		! for the same reason -- p_is_master is not among the callers.
		gmax = host_water_max_resid
#ifdef USEMPI
		CALL mpi_reduce (host_water_max_resid, gmax, 1, MPI_REAL8, MPI_MAX, &
			0, p_comm_worker, p_err)
#endif
		if (p_iam_worker == 0 .and. gmax > host_water_reported) then
			host_water_reported = gmax
			write(6,'(A,E12.4,A,E12.4,A)') &
				' methane host-water: worst disagreement over all ranks ', gmax, &
				' of pore volume (tolerance ', DEF_METHANE%host_water_tolerance, ')'
		endif

	end subroutine methane_host_water_report

	!-----------------------------------------------------------------------
	subroutine methane (istep,ipatch,idate,patchclass,patchtype,&!input
		lb,snl,&
		dlon,dlat,&
		deltim,&
		z_soisno,dz_soisno,zi_soisno,t_soisno,t_grnd,wliq_soisno,wice_soisno,&
			forc_t,forc_pbot,forc_po2m,forc_pco2m,forc_us,forc_vs,&
		zwt,rootfr,snowdp,wat,rsur,etr,lakedepth,dz_lake,t_lake,lake_icefrac,wdsrf,wetwat,bsw,&
		smp,porsl,lai,sai,rootr,&
		annsum_npp,rr,&
		fsatmax,fsatdcf,frcsat,&
		agnpp,bgnpp,somhr,&
		crootfr,lithr,hr_vr,o_scalar,fphr,pot_f_nit_vr,pH,&
			cellorg,t_h2osfc,organic_max,&
			microbial_prod_potential_patch, microbial_oxid_potential_patch, &
		!!!! --------------------------------------------------------------------------------------------------------
		!!!!                                         sum data
		!!!! --------------------------------------------------------------------------------------------------------
		net_methane, &
		methane_prod_depth, o2_decomp_depth, co2_decomp_depth, methane_oxid_depth, o2_oxid_depth, co2_oxid_depth, &
		methane_aere_depth, methane_tran_depth, o2_aere_depth, co2_aere_depth, methane_ebul_depth, &
		o2stress, methane_stress, &
			methane_surf_flux_tot, methane_surf_flux_tot_phys, &
			methane_surf_aere, methane_surf_ebul, methane_surf_diff, methane_surf_diff_phys, &
			methane_balance_residual, methane_ch4_clip_credit, o2_cap_loss, o2_cap_gain, &
		methane_ebul_tot, methane_prod_tot, methane_oxid_tot, &
		co2_decomp_tot, co2_oxid_tot, co2_aere_tot, co2_net_tot, &
		totcol_methane, grnd_methane_cond, conc_o2, conc_methane, &
		!!!! --------------------------------------------------------------------------------------------------------
		!!!! --------------------------------------------------------------------------------------------------------
		!!!!                                         sum data (unsaturated / saturated)
		!!!! --------------------------------------------------------------------------------------------------------
		net_methane_unsat, net_methane_sat, &
		methane_prod_depth_unsat, methane_prod_depth_sat, o2_decomp_depth_unsat, o2_decomp_depth_sat, &
		co2_decomp_depth_unsat, co2_decomp_depth_sat, &
		methane_oxid_depth_unsat, methane_oxid_depth_sat, o2_oxid_depth_unsat, o2_oxid_depth_sat, &
		co2_oxid_depth_unsat, co2_oxid_depth_sat, &
		methane_aere_depth_unsat, methane_aere_depth_sat, methane_tran_depth_unsat, methane_tran_depth_sat, &
		o2_aere_depth_unsat, o2_aere_depth_sat, co2_aere_depth_unsat, co2_aere_depth_sat, &
		methane_ebul_depth_unsat, methane_ebul_depth_sat, &
		o2stress_unsat, o2stress_sat, methane_stress_unsat, methane_stress_sat, &
		methane_surf_flux_tot_unsat, methane_surf_flux_tot_sat, methane_surf_aere_unsat, methane_surf_aere_sat, &
		methane_surf_ebul_unsat, methane_surf_ebul_sat, methane_surf_diff_unsat, methane_surf_diff_sat, &
		methane_surf_diff_phys_unsat, methane_surf_diff_phys_sat, &
		methane_ebul_tot_unsat, methane_ebul_tot_sat, methane_prod_tot_unsat, methane_prod_tot_sat, &
		methane_oxid_tot_unsat, methane_oxid_tot_sat, &
		co2_decomp_tot_unsat, co2_decomp_tot_sat, &
		co2_oxid_tot_unsat, co2_oxid_tot_sat, &
		co2_net_tot_unsat, co2_net_tot_sat, &
		totcol_methane_unsat, totcol_methane_sat, grnd_methane_cond_unsat, grnd_methane_cond_sat, &
		conc_o2_unsat, conc_o2_sat, conc_methane_unsat, conc_methane_sat, &
		methane_prod_depth_lake, methane_oxid_depth_lake, methane_ebul_depth_lake, &
		co2_decomp_depth_lake, co2_oxid_depth_lake, &
		methane_surf_ebul_lake, methane_surf_diff_lake, methane_surf_flux_tot_lake, &
		methane_prod_tot_lake, methane_oxid_tot_lake, methane_ebul_tot_lake, &
		co2_decomp_tot_lake, co2_oxid_tot_lake, co2_net_tot_lake, &
		totcol_methane_lake, grnd_methane_cond_lake, conc_o2_lake, conc_methane_lake, &
		lake_water_ch4_stock, lake_water_o2_stock, lake_frozen_ch4_stock, lake_frozen_o2_stock, &
		lake_liquid_fraction_prev, lake_water_ch4_oxid, &
		lake_sed_ch4_flux, lake_sed_o2_flux, lake_air_o2_flux, &
		!!!! --------------------------------------------------------------------------------------------------------
		c_atm, forc_pmethanem, layer_sat_lag, lake_soilc, &
		annavg_agnpp, annavg_bgnpp, annavg_somhr, annavg_finrw, &
		tempavg_agnpp, tempavg_bgnpp, annsum_counter, tempavg_somhr, tempavg_finrw, &
		fsat_bef, finundated_lag, methane_dfsat_tot, f_h2osfc, &
		is_rice_paddy_in, rice_pft_frac_in, ustar_in, fq_in, raw_grnd_in, &
		store_patch_diagnostics_in, finundated_used_out, finundated_default_out)

		!=======================================================================
		! DESCRIPTION:
		! Driver for the methane emissions model
		!=======================================================================

		!===================== input ===========================================
		integer, intent(in) :: &
			istep            , &
			ipatch           , &! caller's per-patch index
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			patchclass       , &! land patch class of USGS classification or others
			patchtype        , &! land patch type (0=soil, 1=urban or built-up, 2=wetland, 3=land ice, 4=land water bodies, 99=ocean)

			lb               , &! lower bound of array   (snl+1)
			snl				     ! number of snow layers (-5~-1)

		real(r8), intent(in) :: &
			dlon                    , &! longitude (degrees)
			dlat                    , &! latitude (degrees)

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
				forc_pco2m              , &! CO2 partial pressure at observational height [Pa]
				forc_us                 , &! eastward wind speed [m/s], used as U10 proxy for lake gas exchange
				forc_vs                 , &! northward wind speed [m/s], used as U10 proxy for lake gas exchange

			zwt                     , &! the depth from ground (soil) surface to water table [m]

			rootfr   (1:nl_soil)    , &! fraction of roots in each soil layer (the sum of all layer is 1)

			snowdp                  , &! snow depth [m]
			wat                     , &! total water storage [mm] (reserved interface)
			rsur                    , &! surface runoff [mm/s] (reserved interface)
			etr                     , &! transpiration rate [mm/s]
			lakedepth               , &! static/capacity lake depth [m]
			dz_lake  (1:nl_lake)   , &! current lake layer thickness [m]
			t_lake   (1:nl_lake)   , &! current lake layer temperature [K]
			lake_icefrac(1:nl_lake) , &! lake frozen mass fraction [-]
			wdsrf                   , &! depth of surface water [mm]
			wetwat                  , &! water storage in wetland [mm]
			bsw      (1:nl_soil)   	, &! Clapp and Hornberger "b" (nlevgrnd)

			smp      (1:nl_soil)    , &! soil matrix potential [mm]
			porsl    (1:nl_soil)    , &! volumetric soil water at saturation (porosity)
			lai                     , &! leaf area index [m2/m2]
			sai                     , &! stem area index [m2/m2]

			annsum_npp              , &! annual sum NPP (g C/m2/yr)
			rr                      , &! root respiration (fine root MR + total root GR) (gC/m2/s)

			fsatmax                 , &! maximum saturated area fraction [-]
			fsatdcf                 , &! decay factor in calculation of saturated area fraction [1/m]
			frcsat                     ! fraction of saturation area
		!------------------- methane_annualupdate ------------------------------------
		real(r8), intent(in) :: &
			agnpp                   , &! aboveground NPP (gC/m2/s)
			bgnpp                   , &! belowground NPP (gC/m2/s)
			somhr                      ! (gC/m2/s) soil organic matter heterotrophic respiration

		!------------------- methane_prod --------------------------------------------
		real(r8), intent(in) :: &
			crootfr  (1:nl_soil)    , &! fraction of roots for carbon in each soil layer (the sum of all layer is 1)
			lithr                   , &! (gC/m2/s) litter heterotrophic respiration
			hr_vr    (1:nl_soil)    , &! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
			o_scalar (1:nl_soil)    , &! fraction by which decomposition is limited by DEF_METHANE%anoxia
			fphr     (1:nl_soil)    , &! fraction of potential heterotrophic respiration

			pot_f_nit_vr(1:nl_soil) , &! (gN/m3/s) potential soil nitrification flux
			pH                         ! soil water pH

		!------------------- methane_aere --------------------------------------------
		real(r8), intent(in) :: &
			rootr    (1:nl_soil)       ! effective fraction of roots in each soil layer (the sum of all layer is 1)

		!------------------- methane_tran --------------------------------------------
		real(r8), intent(in) :: &
			cellorg  (1:nl_soil)   		, &! column 3D org (kg/m^3 organic matter)
			t_h2osfc               		, &! surface water temperature
			microbial_prod_potential_patch(1:nl_soil), &! optional microbial CH4 production potential (mol/m3/s)
			microbial_oxid_potential_patch(1:nl_soil), &! optional microbial CH4 oxidation potential (mol/m3/s)
			organic_max               		! organic matter content (kg m-3) where soil is assumed to act like peat

		!!!! --------------------------------------------------------------------------------------------------------
		!!!!                                         sum data
		!!!! --------------------------------------------------------------------------------------------------------
		!------------------- ch4_flux, balance, depth variables ------------------------------
		real(r8), intent(out) :: &
			net_methane                     , & ! CH4 oxidation - production, applied to BGC CO2-HR (mol/m2/s)
			methane_prod_depth    (1:nl_soil)   , & ! production of CH4 in each soil layer (mol/m3/s)
			o2_decomp_depth   (1:nl_soil)   , & ! O2 consumption during decomposition in each soil layer (mol/m3/s)
			co2_decomp_depth  (1:nl_soil)   , & ! diagnostic CO2 from decomposition/methanogenesis before O2-stress scaling (mol/m3/s)
			methane_oxid_depth    (1:nl_soil)   , & ! CH4 consumption rate via oxidation in each soil layer (mol/m3/s)
			o2_oxid_depth     (1:nl_soil)   , & ! O2 consumption rate via oxidation in each soil layer (mol/m3/s)
			co2_oxid_depth    (1:nl_soil)   , & ! CO2 production from CH4 oxidation (mol/m3/s)
			methane_aere_depth    (1:nl_soil)   , & ! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s)
			methane_tran_depth    (1:nl_soil)   , & ! CH4 loss rate via transpiration in each soil layer (mol/m3/s)
			o2_aere_depth     (1:nl_soil)   , & ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s)
			co2_aere_depth    (1:nl_soil)   , & ! CO2 aerenchyma diagnostic flux (mol/m3/s)
			methane_ebul_depth    (1:nl_soil)   , & ! CH4 loss rate via ebullition in each soil layer (mol/m3/s)
			o2stress          (1:nl_soil)   , & ! Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs
			methane_stress         (1:nl_soil)   , & ! Ratio of methane available to total per-timestep methane sinks
			methane_surf_flux_tot               , & ! CH4 flux to atmosphere incl. numerical corrections (mol/m2/s)
			methane_surf_flux_tot_phys          , & ! CH4 physical flux before CH4 clip/residual corrections (mol/m2/s)
			methane_surf_aere                   , & ! CH4 surface flux via aerenchyma (mol/m2/s)
			methane_surf_ebul                   , & ! CH4 ebullition flux (mol/m2/s)
				methane_surf_diff                   , & ! CH4 diffusion flux plus numerical closure (mol/m2/s)
				methane_surf_diff_phys              , & ! CH4 pure-physics diffusion flux before clip/residual (mol/m2/s)
				methane_balance_residual            , & ! numerical closure flux credited to methane_surf_diff (mol/m2/s)
				methane_ch4_clip_credit             , & ! negative CH4 clip credited to methane_surf_diff (mol/m2/s)
				o2_cap_loss                        , & ! O2 removed by post-solve physical cap (mol/m2/s)
				o2_cap_gain                        , & ! O2 added by post-solve nonnegative floor (mol/m2/s)
			methane_ebul_tot                    , & ! Total CH4 ebullition (mol/m2/s)
			methane_prod_tot                    , & ! Total CH4 production (mol/m2/s)
			methane_oxid_tot                    , & ! Total CH4 oxidation (mol/m2/s)
			co2_decomp_tot                     , & ! Total CO2 decomposition/methanogenesis (mol/m2/s)
			co2_oxid_tot                       , & ! Total CO2 from CH4 oxidation (mol/m2/s)
			co2_aere_tot                       , & ! Total CO2 aerenchyma diagnostic (mol/m2/s)
			co2_net_tot                            ! Net diagnosed CO2 source (mol/m2/s)

		!------------------- total and concentration variables ------------------------------
		real(r8), intent(inout) :: &
			totcol_methane               , & ! total methane in soil column (mol/m2)
			grnd_methane_cond           , & ! tracer conductance for boundary layer [m/s]
			conc_o2  (1:nl_soil)    , & ! O2 conc in each soil layer (mol/m3)
			conc_methane (1:nl_soil)        ! CH4 conc in each soil layer (mol/m3)
		!!!! --------------------------------------------------------------------------------------------------------

		!!!! --------------------------------------------------------------------------------------------------------
		!!!!                                         sum data (unsaturated / saturated)
		!!!! --------------------------------------------------------------------------------------------------------
		!------------------- ch4_flux, balance, depth variables ------------------------------
		real(r8), intent(out) :: &
			net_methane_unsat               , & ! CH4 oxidation - production, unsaturated subcolumn (mol/m2/s)
			net_methane_sat                 , & ! CH4 oxidation - production, saturated subcolumn (mol/m2/s)
			methane_prod_depth_unsat (1:nl_soil), & ! CH4 production rate in unsaturated soil layer (mol/m3/s)
			methane_prod_depth_sat   (1:nl_soil), & ! CH4 production rate in saturated soil layer (mol/m3/s)
			o2_decomp_depth_unsat(1:nl_soil), & ! O2 consumption during decomposition (unsaturated) (mol/m3/s)
			o2_decomp_depth_sat  (1:nl_soil), & ! O2 consumption during decomposition (saturated) (mol/m3/s)
			co2_decomp_depth_unsat(1:nl_soil), & ! diagnostic CO2 decomp/methanogenesis before O2-stress scaling (unsat) (mol/m3/s)
			co2_decomp_depth_sat  (1:nl_soil), & ! diagnostic CO2 decomp/methanogenesis before O2-stress scaling (sat) (mol/m3/s)
			methane_oxid_depth_unsat (1:nl_soil), & ! CH4 oxidation rate in unsaturated soil layer (mol/m3/s)
			methane_oxid_depth_sat   (1:nl_soil), & ! CH4 oxidation rate in saturated soil layer (mol/m3/s)
			o2_oxid_depth_unsat  (1:nl_soil), & ! O2 oxidation rate in unsaturated soil layer (mol/m3/s)
			o2_oxid_depth_sat    (1:nl_soil), & ! O2 oxidation rate in saturated soil layer (mol/m3/s)
			co2_oxid_depth_unsat (1:nl_soil), & ! CO2 from CH4 oxidation (unsaturated) (mol/m3/s)
			co2_oxid_depth_sat   (1:nl_soil), & ! CO2 from CH4 oxidation (saturated) (mol/m3/s)
			methane_aere_depth_unsat (1:nl_soil), & ! CH4 loss rate via aerenchyma (unsaturated) (mol/m3/s)
			methane_aere_depth_sat   (1:nl_soil), & ! CH4 loss rate via aerenchyma (saturated) (mol/m3/s)
			methane_tran_depth_unsat (1:nl_soil), & ! CH4 loss rate via transpiration (unsaturated) (mol/m3/s)
			methane_tran_depth_sat   (1:nl_soil), & ! CH4 loss rate via transpiration (saturated) (mol/m3/s)
			o2_aere_depth_unsat  (1:nl_soil), & ! O2 gain via aerenchyma (unsaturated) (mol/m3/s)
			o2_aere_depth_sat    (1:nl_soil), & ! O2 gain via aerenchyma (saturated) (mol/m3/s)
			co2_aere_depth_unsat (1:nl_soil), & ! CO2 aerenchyma diagnostic (unsaturated) (mol/m3/s)
			co2_aere_depth_sat   (1:nl_soil), & ! CO2 aerenchyma diagnostic (saturated) (mol/m3/s)
			methane_ebul_depth_unsat (1:nl_soil), & ! CH4 ebullition loss (unsaturated) (mol/m3/s)
			methane_ebul_depth_sat   (1:nl_soil), & ! CH4 ebullition loss (saturated) (mol/m3/s)
			o2stress_unsat       (1:nl_soil), & ! O2 stress ratio (unsaturated)
			o2stress_sat         (1:nl_soil), & ! O2 stress ratio (saturated)
			methane_stress_unsat      (1:nl_soil), & ! CH4 stress ratio (unsaturated)
			methane_stress_sat        (1:nl_soil), & ! CH4 stress ratio (saturated)
			methane_surf_flux_tot_unsat         , & ! CH4 surface flux to atmosphere (unsaturated) (mol/m2/s)
			methane_surf_flux_tot_sat           , & ! CH4 surface flux to atmosphere (saturated) (mol/m2/s)
			methane_surf_aere_unsat             , & ! CH4 surface flux via aerenchyma (unsaturated) (mol/m2/s)
			methane_surf_aere_sat               , & ! CH4 surface flux via aerenchyma (saturated) (mol/m2/s)
			methane_surf_ebul_unsat             , & ! CH4 ebullition flux (unsaturated) (mol/m2/s)
			methane_surf_ebul_sat               , & ! CH4 ebullition flux (saturated) (mol/m2/s)
			methane_surf_diff_unsat             , & ! CH4 diffusion flux (unsaturated) (mol/m2/s)
			methane_surf_diff_sat               , & ! CH4 diffusion flux (saturated) (mol/m2/s)
			methane_surf_diff_phys_unsat        , & ! CH4 pure-physics diff (unsaturated) (mol/m2/s)
			methane_surf_diff_phys_sat          , & ! CH4 pure-physics diff (saturated) (mol/m2/s)
			methane_ebul_tot_unsat              , & ! Total CH4 ebullition (unsaturated) (mol/m2/s)
			methane_ebul_tot_sat                , & ! Total CH4 ebullition (saturated) (mol/m2/s)
			methane_prod_tot_unsat              , & ! Total CH4 production (unsaturated) (mol/m2/s)
			methane_prod_tot_sat                , & ! Total CH4 production (saturated) (mol/m2/s)
			methane_oxid_tot_unsat              , & ! Total CH4 oxidation (unsaturated) (mol/m2/s)
			methane_oxid_tot_sat                , & ! Total CH4 oxidation (saturated) (mol/m2/s)
			co2_decomp_tot_unsat               , & ! Total CO2 decomp/methanogenesis (unsaturated) (mol/m2/s)
			co2_decomp_tot_sat                 , & ! Total CO2 decomp/methanogenesis (saturated) (mol/m2/s)
			co2_oxid_tot_unsat                 , & ! Total CO2 from CH4 oxidation (unsaturated) (mol/m2/s)
			co2_oxid_tot_sat                   , & ! Total CO2 from CH4 oxidation (saturated) (mol/m2/s)
			co2_net_tot_unsat                  , & ! Net diagnosed CO2 source (unsaturated) (mol/m2/s)
			co2_net_tot_sat                        ! Net diagnosed CO2 source (saturated) (mol/m2/s)

		!------------------- total and concentration variables ------------------------------
		real(r8), intent(inout) :: &
			totcol_methane_unsat         , & ! total methane in soil column (unsaturated) (mol/m2)
			totcol_methane_sat           , & ! total methane in soil column (saturated) (mol/m2)
			grnd_methane_cond_unsat     , & ! tracer conductance for boundary layer (unsaturated) [m/s]
			grnd_methane_cond_sat       , & ! tracer conductance for boundary layer (saturated) [m/s]
			conc_o2_unsat (1:nl_soil), & ! O2 conc in unsaturated soil layer (mol/m3)
			conc_o2_sat   (1:nl_soil), & ! O2 conc in saturated soil layer (mol/m3)
			conc_methane_unsat(1:nl_soil), & ! CH4 conc in unsaturated soil layer (mol/m3)
			conc_methane_sat  (1:nl_soil)   ! CH4 conc in saturated soil layer (mol/m3)
		!!!! --------------------------------------------------------------------------------------------------------

		!------------------- lake methane diagnostic/state variables ------------------------------
		real(r8), intent(out) :: &
			methane_prod_depth_lake (1:nl_soil), & ! lake CH4 production rate (mol/m3/s)
			methane_oxid_depth_lake (1:nl_soil), & ! lake CH4 oxidation rate (mol/m3/s)
			methane_ebul_depth_lake (1:nl_soil), & ! lake CH4 ebullition loss rate (mol/m3/s)
			co2_decomp_depth_lake (1:nl_soil), & ! lake CO2 decomp/methanogenesis (mol/m3/s)
			co2_oxid_depth_lake   (1:nl_soil), & ! lake CO2 from CH4 oxidation (mol/m3/s)
			methane_surf_ebul_lake             , & ! lake CH4 surface ebullition flux (mol/m2/s)
			methane_surf_diff_lake             , & ! lake CH4 surface diffusion flux (mol/m2/s)
			methane_surf_flux_tot_lake         , & ! lake total CH4 surface flux (mol/m2/s)
			methane_prod_tot_lake              , & ! lake total CH4 production (mol/m2/s)
			methane_oxid_tot_lake              , & ! lake total CH4 oxidation (mol/m2/s)
			methane_ebul_tot_lake              , & ! lake total CH4 ebullition (mol/m2/s)
			co2_decomp_tot_lake                , & ! lake total CO2 decomp/methanogenesis (mol/m2/s)
			co2_oxid_tot_lake                  , & ! lake total CO2 oxidation product (mol/m2/s)
			co2_net_tot_lake                       ! lake net diagnosed CO2 source (mol/m2/s)

		real(r8), intent(inout) :: &
			totcol_methane_lake           , & ! lake sediment CH4 column stock (mol/m2)
			grnd_methane_cond_lake       , & ! lake tracer conductance (m/s)
			conc_o2_lake      (1:nl_soil), & ! lake O2 concentration (mol/m3)
			conc_methane_lake (1:nl_soil), & ! lake CH4 concentration (mol/m3)
			lake_water_ch4_stock          , & ! well-mixed water CH4 inventory (mol/m2)
			lake_water_o2_stock           , & ! well-mixed water O2 inventory (mol/m2)
			lake_frozen_ch4_stock        , & ! immobile frozen-phase CH4 inventory (mol/m2)
			lake_frozen_o2_stock         , & ! immobile frozen-phase O2 inventory (mol/m2)
			lake_liquid_fraction_prev       ! previous lake liquid fraction [-]

		real(r8), intent(out) :: &
			lake_water_ch4_oxid, & ! water-column CH4 oxidation (mol/m2/s)
			lake_sed_ch4_flux,    & ! sediment-to-water CH4 flux (mol/m2/s)
			lake_sed_o2_flux,     & ! sediment-to-water O2 flux (mol/m2/s)
			lake_air_o2_flux        ! water-to-atmosphere O2 flux (mol/m2/s)
		!!!! --------------------------------------------------------------------------------------------------------

		!------------------- atmospheric and structural variables ------------------------------
		real(r8), intent(out) :: &
			c_atm      (1:3)                 ! CH4, O2, CO2 atmospheric conc (mol/m3)

		real(r8), intent(inout) :: &
			forc_pmethanem              , & ! CH4 concentration in atmosphere (Pa)
			lake_soilc  (1:nl_soil)         ! lake sediment organic carbon per layer (gC/m3)

		! layer_sat_lag is persistent state: unsaturated updates read its previous value.
		real(r8), intent(inout) :: &
			layer_sat_lag(1:nl_soil)         ! lagged saturation ratio per layer (persists across timesteps)

		!------------------- annual accumulators ------------------------------
			real(r8), intent(inout) :: &
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
			tempavg_finrw               ! respiration-weighted temporary average of inundated zones (gC/m2/s)

		!------------------- saturated fraction ------------------------------
		real(r8), intent(inout) :: &
			fsat_bef                , & ! finundated from previous timestep
			finundated_lag          , & ! time-lagged fractional inundated area
			methane_dfsat_tot               ! retained diagnostic; area remapping itself emits no CH4 [mol/m2/s]

		real(r8), intent(in) :: &
            f_h2osfc

		! Rice component metadata is retained in the public interface for the
		! component driver; hydrology is never inferred from these fields.
		logical , intent(in), optional :: is_rice_paddy_in
		logical , intent(in), optional :: store_patch_diagnostics_in
		real(r8), intent(in), optional :: rice_pft_frac_in
		! Current Monin-Obukhov state from CoLM surface fluxes.  Keeping these
		! optional preserves the standalone/cold-start fallback until all callers
		! provide the time-varying state.
		real(r8), intent(in), optional :: ustar_in, fq_in
		! Host water vapour resistance from the ground to the reference height,
		! rd + raw under a canopy [s/m] (MOD_Vars_1DFluxes, raw_grnd).
		real(r8), intent(in), optional :: raw_grnd_in
		real(r8), intent(out), optional :: finundated_used_out, finundated_default_out

		!=================== Local Variables ============================================
		integer  :: i,j,s,l                     ! indices

		integer  :: sat                     ! 0 = unsatured, 1 = saturated
		real(r8) :: finundated              ! fractional inundated area
		logical  :: methane_cold_start      ! true only when no valid previous CH4 inundation state exists
		real(r8) :: finundated_default      ! host/routing inundation before component diagnostics
      logical  :: is_rice_paddy           ! this patch carries a flooded paddy tile
      real(r8) :: rice_pft_frac           ! its area fraction within the patch
      real(r8) :: finundated_rice         ! inundation the paddy tile alone would see
      integer  :: idpp_rice               ! days past planting, -1 when unknown
      integer  :: dsh_rice                ! days since harvest, -1 when unknown
      logical  :: anoxic_decomp_patch     ! flooded soil tile decomposes at the anoxic rate
      real(r8) :: woody_share_unsat       ! woody share of the plant conduit, unsaturated subcolumn
      real(r8) :: woody_share_sat         ! woody share of the plant conduit, saturated subcolumn
      real(r8) :: herb_share_pft          ! grass and crop share of the vegetated PFT cover
      real(r8) :: woody_share_pft         ! tree share of the patch PFT cover
      real(r8) :: ms_start, ms_end        ! midseason drying window, days past planting
      real(r8) :: drain_decay, drain_target
		logical  :: store_patch_diagnostics ! suppress direct patch writes for component solves
		real(r8) :: grnd_methane_cond_base  ! current aerodynamic tracer conductance [m/s]
		real(r8) :: totcolch4_bef
		real(r8) :: totcolch4_bef_sat
		real(r8) :: totcolch4_bef_unsat
		real(r8) :: totcolch4_bef_lake

		real(r8) :: err_methane

		logical, save :: warned_ch4_budget_methane = .false.

		real(r8) :: dfsat
		real(r8) :: redoxlags               ! redox lag time in s
		real(r8) :: redoxlags_vertical      ! Vertical redox lag time in s
			integer  :: dummyfilter(1)          ! empty filter
			integer  :: atm_month, atm_mday
			real(r8) :: atm_methane_mix
			real(r8) :: methane_balance_residual_unsat, methane_balance_residual_sat
			real(r8) :: methane_ch4_clip_credit_unsat, methane_ch4_clip_credit_sat
			real(r8) :: o2_cap_loss_unsat, o2_cap_loss_sat
			real(r8) :: o2_cap_gain_unsat, o2_cap_gain_sat

		real(r8) :: k_h_cc(0:nl_soil,ngases)! ratio of mol/m3 in liquid to mol/m3 in gas [-]

        ! liquid volumetric water content ---- water volume/all volume [m3/m3]
        real(r8) :: vol_aqu            (1:nl_soil)
        real(r8) :: vol_aqu_sat        (1:nl_soil)
        real(r8) :: vol_aqu_unsat      (1:nl_soil)

		! Mobile CH4 storage volume.  Ice-filled pore space is excluded.
        real(r8) :: vol_ch4_storage       (1:nl_soil)
        real(r8) :: vol_ch4_storage_sat   (1:nl_soil)
        real(r8) :: vol_ch4_storage_unsat (1:nl_soil)

        ! air volumetric water content ---- air volume/all volume [m3/m3]
        real(r8) :: vol_gas            (1:nl_soil)
        real(r8) :: vol_gas_sat        (1:nl_soil)
        real(r8) :: vol_gas_unsat      (1:nl_soil)

        ! water-filled proportion
        real(r8) :: f_aqu              (1:nl_soil)
        real(r8) :: f_aqu_sat          (1:nl_soil)
        real(r8) :: f_aqu_unsat        (1:nl_soil)

        ! air-filled proportion
        real(r8) :: f_gas              (1:nl_soil)
        real(r8) :: f_gas_sat          (1:nl_soil)
        real(r8) :: f_gas_unsat        (1:nl_soil)

        ! gas phase CH4 conc in each soil layer (mol/m3)
        real(r8) :: conc_ch4_gas             (1:nl_soil)
        real(r8) :: conc_ch4_gas_sat         (1:nl_soil)
        real(r8) :: conc_ch4_gas_unsat       (1:nl_soil)

        ! aqueous phase CH4 conc in each soil layer (mol/m3)
        real(r8) :: conc_ch4_aqu             (1:nl_soil)
        real(r8) :: conc_ch4_aqu_sat         (1:nl_soil)
        real(r8) :: conc_ch4_aqu_unsat       (1:nl_soil)

        ! CH4 conc in each porosity (mol/m3)
        real(r8) :: conc_ch4_porsl           (1:nl_soil)
        real(r8) :: conc_ch4_porsl_sat       (1:nl_soil)
        real(r8) :: conc_ch4_porsl_unsat     (1:nl_soil)

        ! gas phase CH4 conc in each porosity (mol/m3)
        real(r8) :: conc_ch4_gas_porsl       (1:nl_soil)
        real(r8) :: conc_ch4_gas_porsl_sat   (1:nl_soil)
        real(r8) :: conc_ch4_gas_porsl_unsat (1:nl_soil)

        ! aqueous phase CH4 conc in each porosity (mol/m3)
        real(r8) :: conc_ch4_aqu_porsl       (1:nl_soil)
        real(r8) :: conc_ch4_aqu_porsl_sat   (1:nl_soil)
        real(r8) :: conc_ch4_aqu_porsl_unsat (1:nl_soil)

        ! gas phase O2 conc in each soil layer (mol/m3)
        real(r8) :: conc_o2_gas              (1:nl_soil)
        real(r8) :: conc_o2_gas_sat          (1:nl_soil)
        real(r8) :: conc_o2_gas_unsat        (1:nl_soil)

        ! aqueous phase O2 conc in each soil layer (mol/m3)
        real(r8) :: conc_o2_aqu              (1:nl_soil)
        real(r8) :: conc_o2_aqu_sat          (1:nl_soil)
        real(r8) :: conc_o2_aqu_unsat        (1:nl_soil)

        ! O2 conc in each porosity (mol/m3)
        real(r8) :: conc_o2_porsl            (1:nl_soil)
        real(r8) :: conc_o2_porsl_sat        (1:nl_soil)
        real(r8) :: conc_o2_porsl_unsat      (1:nl_soil)

        ! gas phase O2 conc in each porosity (mol/m3)
        real(r8) :: conc_o2_gas_porsl        (1:nl_soil)
        real(r8) :: conc_o2_gas_porsl_sat    (1:nl_soil)
        real(r8) :: conc_o2_gas_porsl_unsat  (1:nl_soil)

        ! aqueous phase O2 conc in each porosity (mol/m3)
        real(r8) :: conc_o2_aqu_porsl        (1:nl_soil)
        real(r8) :: conc_o2_aqu_porsl_sat    (1:nl_soil)
        real(r8) :: conc_o2_aqu_porsl_unsat  (1:nl_soil)

		real(r8) :: err
		real(r8) :: vliq, vice, vtot, pore_volume, vliq_sat_alloc, vice_sat_alloc
		real(r8) :: host_water_tol
		real(r8) :: host_water_resid, host_water_scale
		real(r8) :: vol_liq_init, vol_ice_init, vol_gas_init
		real(r8) :: wtd_arg
      real(r8) :: hlev                ! water level of a dynamic wetland above the mean surface [m]

		integer  :: jwt                 ! index of the soil layer right above the water table (-)
        integer  :: jwt_sat            ! index of the soil layer right above the water table (-), saturated zone
        integer  :: jwt_unsat          ! index of the soil layer right above the water table (-), unsaturated zone
      real(r8) :: wt_share_unsat     ! share of layer jwt_unsat+1 below the water table (-), wt_layer_sat_share
      integer  :: jwt_tran_unsat     ! jwt_unsat as seen by transport, oxidation, ebullition and plants (-)
		real(r8) :: zwt_sat, wice_soisno_sat(maxsnl+1:nl_soil), wliq_soisno_sat(maxsnl+1:nl_soil), wdsrf_sat
		real(r8) :: zwt_unsat, wice_soisno_unsat(maxsnl+1:nl_soil), wliq_soisno_unsat(maxsnl+1:nl_soil), wdsrf_unsat
      real(r8) :: fliq_sat(1:nl_soil), fliq_unsat(1:nl_soil)   ! unfrozen water fraction wliq/(wliq+wice) per layer
		real(r8) :: routing_depth_mm
		logical  :: lake_restart_debug_print
		real(r8) :: lake_dbg_totcol_bef, lake_dbg_ch4_l1_bef, lake_dbg_o2_l1_bef
		real(r8) :: lake_dbg_ch4_l5_bef, lake_dbg_o2_l5_bef
		real(r8) :: lake_liquid_depth_init, lake_total_depth_init, lake_depth_current
		real(r8) :: lake_liquid_fraction_current, lake_phase_fraction, lake_phase_transfer
		real(r8) :: lake_loss, lake_cc          ! column sediment carbon decomposed and held [gC m-2]
		real(r8) :: lake_cdep                  ! degradable carbon supply to the sediment [gC m-2 yr-1]
		real(r8) :: lake_icebub_stock          ! CH4 of the bubbles held under the lake ice [mol m-2]
		logical  :: lake_icebub_state          ! the patch's store is read from and written to State
      logical  :: use_acceptor               ! acceptor pools apply to this patch
      real(r8) :: acc_reox_unsat(1:nl_soil)  ! potential reoxidation of the oxidised share, unsaturated [s-1]
      real(r8) :: acc_reox_sat  (1:nl_soil)  ! the same, saturated subcolumn [s-1]
      real(r8) :: high_keep                  ! low (CH4-making) share of the tile [-]
      real(r8) :: prod_tran_unsat(1:nl_soil) ! CH4 production entering the pore water, unsaturated [mol/m3/s]
      real(r8) :: prod_tran_sat  (1:nl_soil) ! the same, saturated subcolumn [mol/m3/s]
      real(r8) :: rox_prod_unsat (1:nl_soil) ! rox_on_rooted_production: rhizosphere oxidation, unsaturated [mol/m3/s]
      real(r8) :: rox_prod_sat   (1:nl_soil) ! the same, saturated subcolumn [mol/m3/s]

		!-----------------------------------------------------------------------
		! Set parameters.  When biome lookup enabled, use per-patch redoxlag
		! (set in Driver via get_biome_redoxlag).  Otherwise use global default.
		if (DEF_METHANE%use_biome_redoxlag .and. allocated(biome_redoxlag_patch) .and. &
		    ipatch >= 1 .and. ipatch <= size(biome_redoxlag_patch)) then
			redoxlags = biome_redoxlag_patch(ipatch) * secspday   ! days --> s
		else
			redoxlags = DEF_METHANE%redoxlag * secspday           ! days --> s
		endif
		! Vertical redoxlag uses same biome lookup ratio (sim ratio assumption)
		if (DEF_METHANE%use_biome_redoxlag .and. allocated(biome_redoxlag_patch) .and. &
		    ipatch >= 1 .and. ipatch <= size(biome_redoxlag_patch)) then
			redoxlags_vertical = biome_redoxlag_patch(ipatch) * secspday
		else
			redoxlags_vertical = DEF_METHANE%redoxlag_vertical * secspday
		endif
      ! electron acceptor pools on every non-lake patch
      use_acceptor = DEF_METHANE%acceptor_pool .and. patchtype /= 4 .and. &
         allocated(acceptor_sat) .and. allocated(acceptor_unsat)
      if (use_acceptor) use_acceptor = ipatch >= 1 .and. ipatch <= size(acceptor_sat, 2)
      acc_reox_unsat(:) = 0._r8
      acc_reox_sat(:)   = 0._r8
      ! at a tower, only the low share of a wetland tile makes
      ! CH4; the high share (tussocks, rims, high centres) makes none
      high_keep = 1._r8
#ifdef SinglePoint
      if (patchtype == 2 .and. DEF_METHANE%tundra_high_share_site > 0._r8) &
         high_keep = 1._r8 - min(DEF_METHANE%tundra_high_share_site, 1._r8)
#endif

		! Match CoLM's current water-vapour aerodynamic conductance only where
		! fq is the full surface-layer integral: lakes and effectively bare soil.
		! Vegetated tiles use fq-fqt plus canopy/ground resistances inside
		! MOD_LeafTemperature, so reconstructing kappa*ustar/fq there would
		! overestimate gas exchange. They keep the configured fallback unless
		! grnd_methane_cond_host takes the resolved resistance (below).
		grnd_methane_cond_base = DEF_METHANE%grnd_methane_cond_default
		if ((patchtype == 4 .or. lai + sai <= 1.e-6_r8) .and. &
		    present(ustar_in) .and. present(fq_in)) then
			if (.not. ieee_is_nan(ustar_in) .and. .not. ieee_is_nan(fq_in) .and. &
			    abs(ustar_in) < 0.5_r8*abs(spval) .and. abs(fq_in) < 0.5_r8*abs(spval) .and. &
			    ustar_in > 0._r8 .and. fq_in > 1.e-8_r8) then
				grnd_methane_cond_base = vonkar * ustar_in / fq_in
				if (grnd_methane_cond_base <= 0._r8 .or. grnd_methane_cond_base > 1._r8) &
					grnd_methane_cond_base = DEF_METHANE%grnd_methane_cond_default
			endif
		endif
		! a vegetated non-lake tile takes the host resistance from
		! the ground through the canopy air to the reference height, rd + raw,
		! as CLM passes 1/(raw(p,above_canopy)+raw(p,below_canopy)) to its
		! methane model (CanopyFluxesMod.F90; BareGroundFluxesMod.F90 uses 1/raw
		! on bare ground). Bounded to 1e-5..1 m/s (resistance 1..1e5 s/m) so a
		! near-zero u* under a stable night cannot shut the surface off, nor a
		! degenerate profile open it without limit; a missing or invalid value
		! keeps grnd_methane_cond_default.
		if (DEF_METHANE%grnd_methane_cond_host .and. patchtype /= 4 .and. &
		    lai + sai > 1.e-6_r8 .and. present(raw_grnd_in)) then
			if (.not. ieee_is_nan(raw_grnd_in) .and. abs(raw_grnd_in) < 0.5_r8*abs(spval) .and. &
			    raw_grnd_in > 0._r8) then
				grnd_methane_cond_base = min(max(1._r8/raw_grnd_in, 1.e-5_r8), 1._r8)
			endif
		endif

		! Initialize fluxes and diagnostics to zero
		lake_depth_current = sum(max(dz_lake(1:nl_lake), 0._r8))
		if (lake_depth_current <= 1.e-12_r8) lake_depth_current = max(lakedepth, 0._r8)
		methane_prod_depth        = 0._r8
		o2_decomp_depth          = 0._r8
		co2_decomp_depth         = 0._r8
		methane_oxid_depth       = 0._r8
		o2_oxid_depth            = 0._r8
		co2_oxid_depth           = 0._r8
		methane_aere_depth       = 0._r8
		methane_tran_depth       = 0._r8
		o2_aere_depth            = 0._r8
		co2_aere_depth           = 0._r8
		methane_ebul_depth       = 0._r8
		o2stress                 = 1._r8
		methane_stress           = 1._r8
		methane_ebul_tot         = 0._r8

			methane_surf_flux_tot       = 0._r8
			methane_surf_flux_tot_sat   = 0._r8
			methane_surf_flux_tot_unsat = 0._r8
				methane_balance_residual    = 0._r8
				methane_ch4_clip_credit     = 0._r8
				methane_balance_residual_unsat = 0._r8
				methane_balance_residual_sat   = 0._r8
				methane_ch4_clip_credit_unsat = 0._r8
				methane_ch4_clip_credit_sat   = 0._r8
				o2_cap_loss             = 0._r8
				o2_cap_loss_unsat       = 0._r8
				o2_cap_loss_sat         = 0._r8
				o2_cap_gain             = 0._r8
				o2_cap_gain_unsat       = 0._r8
				o2_cap_gain_sat         = 0._r8

		methane_prod_tot            = 0._r8
		methane_prod_tot_sat        = 0._r8
		methane_prod_tot_unsat      = 0._r8

		methane_oxid_tot            = 0._r8
		methane_oxid_tot_sat        = 0._r8
		methane_oxid_tot_unsat      = 0._r8

		methane_prod_tot_lake       = 0._r8
		methane_oxid_tot_lake       = 0._r8
		methane_ebul_tot_lake       = 0._r8
		methane_surf_ebul_lake      = 0._r8
		methane_surf_diff_lake      = 0._r8
		methane_surf_flux_tot_lake  = 0._r8
		lake_water_ch4_oxid         = 0._r8
		lake_sed_ch4_flux            = 0._r8
		lake_sed_o2_flux             = 0._r8
		lake_air_o2_flux             = 0._r8
		methane_prod_depth_lake     = 0._r8
		methane_oxid_depth_lake     = 0._r8
		methane_ebul_depth_lake     = 0._r8
		co2_decomp_depth_lake       = 0._r8
		co2_oxid_depth_lake         = 0._r8
		co2_decomp_tot_lake         = 0._r8
		co2_oxid_tot_lake           = 0._r8
		co2_net_tot_lake            = 0._r8

		! Adjustment to NEE for methane production - oxidation
		net_methane             = 0._r8
		net_methane_sat         = 0._r8
		net_methane_unsat       = 0._r8

			CALL julian2monthday (idate(1), idate(2), atm_month, atm_mday)
			atm_methane_mix = methane_atm_mixing_ratio(idate(1), atm_month)
			! Offline CH4 can use either the scalar atmospheric value
			! DEF_METHANE%atm_methane or an optional transient file table.
			if (DEF_METHANE%methane_offline) then
				forc_pmethanem = atm_methane_mix*forc_pbot
				! [Pa]     =[mol/mol]*[Pa]
			else
				if (ieee_is_nan(forc_pmethanem) .or. &
				    abs(forc_pmethanem) >= 0.5_r8 * abs(spval) .or. &
				    forc_pmethanem <= 0._r8) then
					forc_pmethanem = atm_methane_mix*forc_pbot
				end if
			end if
		c_atm(1) =  forc_pmethanem / rgasm / forc_t
		c_atm(2) =  forc_po2m  / rgasm / forc_t
		c_atm(3) =  max(forc_pco2m, 0._r8) / rgasm / forc_t
		! n/V = P/RT
		![mol/m3]=[Pa]         /[J/K/mol]/[K]
		![mol/m3]=[J/m3]       *[K*mol/J]*[1/K]
		totcolch4_bef = totcol_methane
		totcolch4_bef_sat = totcol_methane_sat
		totcolch4_bef_unsat = totcol_methane_unsat
		totcolch4_bef_lake = totcol_methane_lake
      ! the tile store is the low share's; rebuild it from the
      ! subcolumn stores so a restart written with another share, or with the
      ! switch off, gives the same budget reference
      if (high_keep < 1._r8 .and. fsat_bef >= 0._r8 .and. fsat_bef <= 1._r8) &
         totcolch4_bef = (totcolch4_bef_sat * fsat_bef + totcolch4_bef_unsat * (1._r8 - fsat_bef)) * high_keep
		if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
			! Generic totcol_methane is the canonical lake inventory and now
			! includes both sediment and the independent water-column stock.
			totcolch4_bef = totcol_methane
			totcolch4_bef_sat = totcolch4_bef_lake
			totcolch4_bef_unsat = 0._r8
		endif
		lake_restart_debug_print = DEF_METHANE%lake_restart_debug .and. &
			patchtype == 4 .and. DEF_METHANE%allowlakeprod .and. &
			idate(1) == DEF_METHANE%lake_restart_debug_year .and. &
			idate(2) >= DEF_METHANE%lake_restart_debug_start_doy .and. &
			idate(2) <= DEF_METHANE%lake_restart_debug_end_doy .and. &
			(DEF_METHANE%lake_restart_debug_sec < 0 .or. &
			 idate(3) == DEF_METHANE%lake_restart_debug_sec) .and. &
			totcolch4_bef_lake > 0._r8
		if (lake_restart_debug_print) then
			lake_dbg_totcol_bef = totcolch4_bef_lake
			lake_dbg_ch4_l1_bef = conc_methane_lake(1)
			lake_dbg_o2_l1_bef = conc_o2_lake(1)
			lake_dbg_ch4_l5_bef = conc_methane_lake(min(5,nl_soil))
			lake_dbg_o2_l5_bef = conc_o2_lake(min(5,nl_soil))
			write(6,'(A,1X,I6,1X,I4,1X,I7,1X,I9,1X,2F10.3,1X,24ES14.6)') &
				'CH4_LAKE_ENTRY_DEBUG', idate(1), idate(2), idate(3), ipatch, dlon, dlat, &
				lakedepth, wdsrf, t_grnd, t_h2osfc, t_soisno(1), lake_icefrac(1), &
				wliq_soisno(1), wice_soisno(1), smp(1), f_h2osfc, fsat_bef, finundated_lag, &
				layer_sat_lag(1), layer_sat_lag(min(5,nl_soil)), totcolch4_bef_lake, &
				conc_methane_lake(1), conc_o2_lake(1), conc_methane_lake(min(5,nl_soil)), &
				conc_o2_lake(min(5,nl_soil)), sum(lake_soilc(1:nl_soil)), forc_t, &
				forc_pbot, forc_po2m, forc_pmethanem
		else
			lake_dbg_totcol_bef = 0._r8
			lake_dbg_ch4_l1_bef = 0._r8
			lake_dbg_o2_l1_bef = 0._r8
			lake_dbg_ch4_l5_bef = 0._r8
			lake_dbg_o2_l5_bef = 0._r8
		endif

		totcol_methane = 0.
		totcol_methane_sat   = 0.
		totcol_methane_unsat = 0.

		finundated = 0._r8
			if (DEF_wetland_finundation_scheme == 1) then
				! Wetland tiles fully saturated (finundated=1), all other patch
				! types unsaturated (finundated=0).  The user-facing wetwat mode
				! starts from this classification and then replaces wetland tiles
				! with wetwat/wetwat_ref below.
				if (patchtype == 2) then
					finundated = 1._r8
				else
					finundated = 0._r8
				endif
			elseif (DEF_wetland_finundation_scheme == 2) then
						finundated = frcsat
			elseif (DEF_wetland_finundation_scheme == 3) then
				finundated = f_h2osfc
			elseif (DEF_wetland_finundation_scheme == 4) then
				if (patchtype == 2) then
					finundated = 1._r8
				else
					finundated = f_h2osfc
				endif
			elseif (DEF_wetland_finundation_scheme == 5) then
				! GIEMS-MC monthly wetland climatology, looked up per-patch.
				! Loaded once at init by MOD_Tracer_Reactive_Methane_GIEMS (see CoLM.F90).
			! ipatch is passed in through methane()'s arg list and forwarded
			! to giems_finundated() below.
			if (.not. giems_active) then
				write(6,*) 'ERROR: methane scheme 5 requires an active GIEMS file.'
				CALL CoLM_stop ()
				endif
				finundated = giems_finundated(ipatch, idate(1), idate(2))
				! GIEMS-MC is a grid-cell absolute fraction.  Split it once into
				! disjoint wetland-first wetland/soil patch fractions.
				IF (allocated(wetland_frac_per_patch) .and. &
				    ipatch >= 1 .and. ipatch <= size(wetland_frac_per_patch)) THEN
					finundated = methane_distribute_grid_finundation(finundated, &
						wetland_frac_per_patch(ipatch), patchtype)
				ENDIF
			elseif (DEF_wetland_finundation_scheme == 6) then
			! Water-table-driven inundation via logistic S-curve:
			!   finundated = 1 / (1 + exp((zwt - wtd_inflection) / wtd_steepness))
			! Defaults (wtd_inflection=0.30 m, wtd_steepness=0.05 m) reproduce
			! the "30 cm threshold" as the inflection point (finundated=0.5)
			! while smoothing the original step function (29 cm = 1 vs 31 cm = 0
			! creates unphysical 100x jump in CH4 production).
			! References:
			!   Walter & Heimann 2000 (sigmoid f(zwt))
			!   Sundh 2000 (CH4 ~ exp(-zwt/0.3) for boreal peatland)
			!   Bridgham 2013 GCB (review).
			if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
				finundated = 1._r8
			elseif (patchtype /= 2 .and. DEF_METHANE%use_routing_for_soil .and. &
			        allocated(f_inund_flood_patch) .and. &
			        ipatch >= 1 .and. ipatch <= size(f_inund_flood_patch)) then
				! HYBRID: soil tile uses routing-published f_inund_flood_patch
				! (captures Pantanal/Amazon overland flooding that sigmoid(zwt)
				! cannot resolve because soil zwt is too deep there).
				! Optional threshold gate: only flood when fldfrc > threshold
				! (Bates 2010, Coe 2002 hydrological connectivity concept).
				if (f_inund_flood_patch(ipatch) > DEF_METHANE%hybrid_soil_threshold) then
					finundated = f_inund_flood_patch(ipatch)
					if (allocated(wetland_frac_per_patch) .and. &
					    ipatch <= size(wetland_frac_per_patch)) then
						finundated = methane_distribute_grid_finundation(finundated, &
							wetland_frac_per_patch(ipatch), patchtype)
					endif
				else
					finundated = 0._r8
				endif
			else
				if (DEF_METHANE%wtd_steepness <= 0._r8 .or. ieee_is_nan(DEF_METHANE%wtd_steepness)) then
					write(6,*) 'ERROR: DEF_METHANE%wtd_steepness must be positive for methane scheme 6.'
					CALL CoLM_stop ()
				endif
				if (ieee_is_nan(zwt) .or. abs(zwt) >= 0.5_r8*abs(spval)) then
					write(6,*) 'ERROR: invalid zwt for methane scheme 6 logistic inundation: ', zwt
					CALL CoLM_stop ()
				endif
				! Dynamic flood extension: soil tile (patchtype==0/1) uses
				! a separate (typically deeper) sigmoid inflection.  This lets
				! seasonally-wet upland (Pantanal floodplain etc.) contribute
				! CH4 when soil zwt rises into 0.5-1.5 m range, while keeping
				! the steep 0.3 m wetland sigmoid that captures peat-style
				! peak-to-trough variation.  Disabled when wtd_inflection_soil=0.
				IF (patchtype /= 2 .and. DEF_METHANE%wtd_inflection_soil > 0._r8) then
					wtd_arg = (zwt - DEF_METHANE%wtd_inflection_soil) / &
					          max(DEF_METHANE%wtd_steepness_soil, 1.e-3_r8)
				else
					wtd_arg = (zwt - DEF_METHANE%wtd_inflection) / DEF_METHANE%wtd_steepness
				endif
				wtd_arg = min(50._r8, max(-50._r8, wtd_arg))
				finundated = 1._r8 / (1._r8 + exp(wtd_arg))
			endif
		elseif (DEF_wetland_finundation_scheme == 7) then
				! GridRiverLake fldfrc: per-patch general flood fraction
				! (levee+floodplain) from grid-based river/lake routing.
					! Published by MOD_Grid_RiverLakeFlow::publish_fldfrc_to_patches
					! at end of routing.  Methane consumes the most recently published
					! value, normally from the previous land step.  Fall back to the
					! levee variable if flood var happens to be unallocated (e.g.
				IF (allocated(f_inund_flood_patch) .and. &
				    ipatch >= 1 .and. &
				    ipatch <= size(f_inund_flood_patch)) THEN
				   finundated = f_inund_flood_patch(ipatch)
				ELSE IF (allocated(f_inund_levee_patch) .and. &
				    ipatch >= 1 .and. &
				    ipatch <= size(f_inund_levee_patch)) THEN
				   finundated = f_inund_levee_patch(ipatch)
					ELSE
						   write(6,*) 'ERROR: methane scheme 7 requires published GridRiverLake levee/floodplain fraction for patch ', &
					      ipatch
					   CALL CoLM_stop ()
					ENDIF
					! Routing is grid-relative.  Allocate it to mapped wetland first;
					! only the residual can inundate the soil patch.
					IF (allocated(wetland_frac_per_patch) .and. &
					    ipatch >= 1 .and. ipatch <= size(wetland_frac_per_patch)) THEN
						finundated = methane_distribute_grid_finundation(finundated, &
							wetland_frac_per_patch(ipatch), patchtype)
					ENDIF
      elseif (DEF_wetland_finundation_scheme == 8) then
         ! CoLM mode.  The mapped wetland tile is permanent wetland and is
         ! saturated.  Every other tile takes the routing-published flood
         ! fraction (levee field as the fallback): grid flood is allocated
         ! to the mapped wetland first so only the residual reaches the
         ! soil, gated by hybrid_soil_threshold like the hybrid soil path.
         ! Without routing the published fields stay zero.
         if (patchtype == 2) then
            if (DEF_USE_Dynamic_Wetland) then
               ! Dynamic wetland: the tile's own water table decides; the
               ! unsaturated column carries the rest of the patch.
               if (DEF_METHANE%wetland_microtopo_sigma > 0._r8) then
                  ! Share of a normally distributed surface lying below the
                  ! water level (PEAT-CLSM microtopography). Standing water
                  ! sets the level: the host takes zwt from liquid water, so
                  ! an ice-filled top layer under a pond puts it below ground.
                  if (wdsrf > 0._r8) then
                     hlev = wdsrf * 1.e-3_r8
                  else
                     hlev = -zwt
                  endif
                  finundated = 0.5_r8 * erfc(-hlev / (sqrt(2._r8) * DEF_METHANE%wetland_microtopo_sigma))
               else if (zwt <= 0._r8) then
                  finundated = 1._r8
               else
                  finundated = 0._r8
               endif
            else
               finundated = 1._r8
            endif
         else
            finundated = 0._r8
            if (allocated(f_inund_flood_patch) .and. &
               ipatch >= 1 .and. ipatch <= size(f_inund_flood_patch)) then
               finundated = f_inund_flood_patch(ipatch)
            elseif (allocated(f_inund_levee_patch) .and. &
               ipatch >= 1 .and. ipatch <= size(f_inund_levee_patch)) then
               finundated = f_inund_levee_patch(ipatch)
            endif
            if (finundated > DEF_METHANE%hybrid_soil_threshold) then
               if (allocated(wetland_frac_per_patch) .and. &
                  ipatch >= 1 .and. ipatch <= size(wetland_frac_per_patch)) then
                  finundated = methane_distribute_grid_finundation(finundated, &
                     wetland_frac_per_patch(ipatch), patchtype, &
                     debit_wetland=.false.)
               endif
               ! the soil floods only within the GLWD riverine and
               ! lacustrine floodplain of its cell (DEF_METHANE%
               ! floodplain_glwd_cap_file); fld_cap_p is 1 off the soil tiles
               if (allocated(fld_cap_p) .and. ipatch >= 1) then
                  if (ipatch <= size(fld_cap_p)) finundated = min(finundated, fld_cap_p(ipatch))
               endif
            else
               finundated = 0._r8
            endif
         endif
#ifdef SinglePoint
         ! Observed water table (single point): the table decides the
         ! classification, so one column carries the whole patch.
         if (USE_SITE_WTD .and. SITE_WTD_valid) then
            if (zwt <= 0._r8) then
               finundated = 1._r8
            else
               finundated = 0._r8
            endif
         endif
#endif
			else
			write(6,*) 'ERROR: unsupported DEF_wetland_finundation_scheme = ', &
				DEF_wetland_finundation_scheme
			CALL CoLM_stop ()
		endif
		if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
			! Lake sediment remains the saturated soil-column branch.  Its
			! diffusive CH4/O2 exchange is coupled below to a prognostic,
			! well-mixed water-column node; ebullition crosses the water column
			! directly.  Force the enabled lake gate onto the saturated branch
			! rather than silently weighting it away through the wetland scheme.
			finundated = 1._r8
		endif
		! NaN/spval defense before clamp.  Sources frcsat / f_h2osfc /
		! f_inund_levee_patch may carry sentinel values if forcing or
		! upstream module misbehaves; min/max would propagate NaN.
		if (ieee_is_nan(finundated) .or. abs(finundated) >= 0.5_r8*abs(spval)) then
			finundated = 0._r8
		endif
		if (finundated < 1.e-10_r8) finundated = 0._r8
		! Clamp to keep saturated/unsaturated weighting bounded.
		finundated = min(max(finundated, 0._r8), 1._r8)

		! CoLM's wetland hydrology hardcodes zwt=0 and frcsat=1 for
		! patchtype==2 (MOD_SoilSnowHydrology.F90:1227-1265), so the IGBP /
		! TOPMODEL / sigmoid / routing schemes all collapse to ~1.0 or near-
		! const on wetland tiles, killing the seasonal cycle.  But the WATER
		! BALANCE itself (wetwat = wdsrf + wa + wetwat + (gwat - etr +
		! qsdew + qfros - qsubl)*deltim) IS dynamic.  This override turns
		! wetwat / wetwatmax into the per-step finundated, so seasonality
		! flows from forced precip-ET balance directly into methane physics.
		if (DEF_METHANE%enable_wetwat_finundated_override .and. patchtype == 2) then
			if (wetwatmax > 0._r8) then
				finundated = max(0._r8, min(1._r8, wetwat / wetwatmax))
			endif
		endif

		! Restart sanitation normally converts NaNs to spval.  Repeat the guard at
		! the consuming boundary so a future caller cannot feed NaN into dfsat or
		! the redox-memory update (ordered comparisons do not catch NaN).
		if (ieee_is_nan(fsat_bef)) fsat_bef = spval
		methane_cold_start = (abs(fsat_bef) >= 1.e30_r8 .or. fsat_bef < 0._r8 .or. fsat_bef > 1._r8)
		if (methane_cold_start) then
			fsat_bef = finundated
		endif
		if (ieee_is_nan(finundated_lag)) finundated_lag = spval
		if (abs(finundated_lag) >= 1.e30_r8 .or. finundated_lag < 0._r8 .or. finundated_lag > 1._r8) then
			finundated_lag = finundated
		endif
		do j = 1, nl_soil
			if (ieee_is_nan(layer_sat_lag(j))) layer_sat_lag(j) = spval
			if (abs(layer_sat_lag(j)) >= 1.e30_r8 .or. &
			    layer_sat_lag(j) < 0._r8 .or. &
			    layer_sat_lag(j) > 1._r8) then
				layer_sat_lag(j) = finundated
			end if
		end do
		if (snowdp > 0._r8) then  !If snow_depth>0, keep finundated from the previous time step of snow season. (by Xiyan Xu, 05/2016)
            finundated = fsat_bef
		end if
		! Re-apply the lake force AFTER the snow override.  A snow-covered
		! lake patch must still report finundated=1 to stay on the saturated
		! sediment branch (otherwise a stale fsat_bef from a previous
		! mid-thaw step could put the lake back on the unsaturated branch
		! and zero-out CH4 production).  See lake gate ~line 658.
		if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
			finundated = 1._r8
		endif

      ! Paddy water management (removed by 68a507f8, whose comment here read
      ! "methane must not invent irrigation water"; restored from 7550e0ed).
      !
      ! It does invent it: the floor below raises the inundation the CH4 column
      ! sees without adding a drop to the host water balance. That is the known
      ! cost of this rollback and it is deliberate -- see the parameter block in
      ! MOD_Tracer_Reactive_Methane_Const. Without it a paddy is hydrologically a
      ! dry field: all seven FLUXNET-CH4 rice towers model exactly 0.0 against an
      ! observed 66.7 mg CH4 m-2 d-1, and no inundation_mode can repair that
      ! because rice has no finundated of its own -- it reads the soil one.
      !
      ! Blending by area keeps a mixed-CFT patch honest, and max() means an
      ! already-wet patch is never dried by the paddy floor.
      IF (present(is_rice_paddy_in)) THEN
         is_rice_paddy = is_rice_paddy_in
      ELSE
         is_rice_paddy = .false.
      ENDIF
      ! floodplains on the soil tile; wetland tiles keep their own BGC
      ! anoxia limit and paddy patches are left out.
      anoxic_decomp_patch = DEF_METHANE%floodplain_anoxic_decomp .and. &
         patchtype /= 2 .and. .not. is_rice_paddy
      IF (present(rice_pft_frac_in)) THEN
         rice_pft_frac = rice_pft_frac_in
      ELSE
         rice_pft_frac = 0._r8
      ENDIF
      rice_pft_frac = min(max(rice_pft_frac, 0._r8), 1._r8)

		finundated_default = finundated
      IF (is_rice_paddy .and. rice_pft_frac > 0._r8 .and. &
          (DEF_METHANE%rice_paddy_min_finundated > 0._r8 .or. DEF_METHANE%rice_host_water)) THEN
         IF (DEF_METHANE%rice_host_water) THEN
            ! the host's water on the field decides (irrigation, bunds, rain).
            IF (wdsrf >= 1._r8 .or. zwt <= 0._r8) THEN
               finundated_rice = 1._r8
            ELSE
               finundated_rice = finundated_default
            ENDIF
            finundated = rice_pft_frac * finundated_rice &
                       + (1._r8 - rice_pft_frac) * finundated_default
         ELSEIF (is_paddy_rice_live(ipatch)) THEN
            finundated_rice = max(finundated_default, &
               DEF_METHANE%rice_paddy_min_finundated)

            ! Midseason drying. croplive_p carries the season across New Year on
            ! its own, so southern-hemisphere rice needs no special case.
            idpp_rice = rice_days_past_planting(ipatch, idate(2), idate(1))
            ms_start  = DEF_METHANE%rice_midseason_start_days
            ms_end    = ms_start + DEF_METHANE%rice_midseason_drain_days
            IF (idpp_rice >= 0 .and. real(idpp_rice, r8) >= ms_start &
                .and. real(idpp_rice, r8) < ms_end) THEN
               finundated_rice = min(finundated_rice, &
                  max(DEF_METHANE%rice_midseason_drained_finundated, &
                      finundated_default))
            ENDIF
            finundated = rice_pft_frac * finundated_rice &
                       + (1._r8 - rice_pft_frac) * finundated_default
         ELSE
            ! Harvested: decay the floor to the host value rather than dropping
            ! to it, so the post-harvest emission tail is not cut off in a step.
            dsh_rice = rice_days_since_harvest(ipatch, idate(2), idate(1))
            IF (dsh_rice >= 0 .and. &
                real(dsh_rice, r8) < DEF_METHANE%rice_drain_window_days) THEN
               drain_decay  = 1._r8 - real(dsh_rice, r8) &
                            / DEF_METHANE%rice_drain_window_days
               drain_target = drain_decay * DEF_METHANE%rice_paddy_min_finundated
               finundated_rice = max(drain_target, finundated_default)
               finundated = rice_pft_frac * finundated_rice &
                          + (1._r8 - rice_pft_frac) * finundated_default
            ENDIF
         ENDIF
         finundated = min(max(finundated, 0._r8), 1._r8)
      ENDIF

		dfsat = finundated - fsat_bef

		! Update lagged finundated for redox calculation
		if (redoxlags > 0._r8) then
			finundated_lag = finundated_lag * exp(-deltim/redoxlags) &
					+ finundated * (1._r8 - exp(-deltim/redoxlags))
		else
			finundated_lag = finundated
		end if

		methane_dfsat_tot = 0._r8
		do j=1,nl_soil
			if (.not. methane_cold_start) then
				if (dfsat > 0._r8) then
					! Newly flooded area carries the dry-column concentration;
					! this convex remap conserves the area-weighted inventory.
					conc_methane_sat(j) = (fsat_bef*conc_methane_sat(j) + &
						dfsat*conc_methane_unsat(j)) / finundated
					conc_o2_sat (j) = (fsat_bef*conc_o2_sat (j) + &
						dfsat*conc_o2_unsat (j)) / finundated
				elseif (dfsat < 0._r8) then
					! The area leaving the saturated subcolumn becomes unsaturated.
					! Regrid both gases into that enlarged unsaturated area so the
					! area-weighted inventories are unchanged by bookkeeping alone.
					if (finundated < 1._r8) then
						conc_methane_unsat(j) = ((1._r8 - fsat_bef) * conc_methane_unsat(j) - &
							dfsat * conc_methane_sat(j)) / (1._r8 - finundated)
						conc_o2_unsat(j) = ((1._r8 - fsat_bef) * conc_o2_unsat(j) - &
							dfsat * conc_o2_sat(j)) / (1._r8 - finundated)
					endif
				end if
			end if
		end do

      ! the methanogen activity moves with the area between the
      ! two subcolumns the same way (area-weighted mean of what arrives)
      if (DEF_METHANE%methanogen_activity .and. patchtype == 2 .and. .not. methane_cold_start .and. &
          allocated(methanogen_act_sat) .and. allocated(methanogen_act_unsat)) then
         if (ipatch >= 1 .and. ipatch <= size(methanogen_act_sat, 2)) then
            if (dfsat > 0._r8) then
               methanogen_act_sat(:,ipatch) = (fsat_bef * methanogen_act_sat(:,ipatch) + &
                  dfsat * methanogen_act_unsat(:,ipatch)) / finundated
            elseif (dfsat < 0._r8 .and. finundated < 1._r8) then
               methanogen_act_unsat(:,ipatch) = ((1._r8 - fsat_bef) * methanogen_act_unsat(:,ipatch) - &
                  dfsat * methanogen_act_sat(:,ipatch)) / (1._r8 - finundated)
            endif
         endif
      endif

      ! the acceptor pools move with the area in the same way;
      ! both subcolumns share a layer's capacity, so the oxidised shares mix
      if (use_acceptor .and. .not. methane_cold_start) then
         if (dfsat > 0._r8) then
            acceptor_sat(:,ipatch) = (fsat_bef * acceptor_sat(:,ipatch) + &
               dfsat * acceptor_unsat(:,ipatch)) / finundated
         elseif (dfsat < 0._r8 .and. finundated < 1._r8) then
            acceptor_unsat(:,ipatch) = ((1._r8 - fsat_bef) * acceptor_unsat(:,ipatch) - &
               dfsat * acceptor_sat(:,ipatch)) / (1._r8 - finundated)
         endif
      endif

		!!!! Begin biochemistry
		! First for soil
		! Do CH4 Annual Averages
		if (patchtype /= 4) then
			call methane_annualupdate(idate, finundated, deltim, agnpp, bgnpp, somhr, &
				annavg_agnpp, annavg_bgnpp, annavg_somhr, annavg_finrw, &
				tempavg_agnpp, tempavg_bgnpp, annsum_counter, tempavg_somhr, tempavg_finrw)
		endif

		call henry_law(t_grnd,t_soisno,k_h_cc)

		! methane_tran updates a local copy of the patch's under-ice bubble
		! store, which lives in MOD_Tracer_Reactive_Methane_State with the restart
		lake_icebub_stock = 0._r8
		lake_icebub_state = .false.
		if (patchtype == 4 .and. DEF_METHANE%allowlakeprod .and. allocated(lake_icebubble_ch4_stock)) &
			lake_icebub_state = ipatch >= 1 .and. ipatch <= size(lake_icebubble_ch4_stock)
		if (lake_icebub_state) lake_icebub_stock = lake_icebubble_ch4_stock(ipatch)

		if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
			lake_total_depth_init = sum(max(dz_lake(1:nl_lake), 0._r8))
			lake_liquid_depth_init = sum(max(dz_lake(1:nl_lake), 0._r8) * &
				(1._r8 - min(max(lake_icefrac(1:nl_lake), 0._r8), 1._r8)))
			if (lake_total_depth_init <= 1.e-12_r8 .and. lakedepth > 0._r8) then
				lake_total_depth_init = lakedepth
				lake_liquid_depth_init = lakedepth * &
					(1._r8 - min(max(lake_icefrac(1), 0._r8), 1._r8))
			endif
			if (lake_total_depth_init > 1.e-12_r8) then
				lake_liquid_fraction_current = min(max( &
					lake_liquid_depth_init / lake_total_depth_init, 0._r8), 1._r8)
			else
				lake_liquid_fraction_current = 0._r8
			endif

			if (methane_cold_start) then
				lake_water_ch4_stock = max(lake_liquid_depth_init, 0._r8) * k_h_cc(0,1) * c_atm(1)
				lake_water_o2_stock = max(lake_liquid_depth_init, 0._r8) * k_h_cc(0,2) * c_atm(2)
				lake_frozen_ch4_stock = 0._r8
				lake_frozen_o2_stock = 0._r8
				lake_icebub_stock = 0._r8
			else if (ieee_is_nan(lake_liquid_fraction_prev) .or. &
			         abs(lake_liquid_fraction_prev) >= 0.5_r8 * abs(spval)) then
				! Deterministic migration for pre-schema-4 restarts: partition the
				! old undifferentiated water inventory at the current phase fraction.
				lake_phase_fraction = 1._r8 - lake_liquid_fraction_current
				lake_phase_transfer = lake_water_ch4_stock * lake_phase_fraction
				lake_water_ch4_stock = lake_water_ch4_stock - lake_phase_transfer
				lake_frozen_ch4_stock = lake_frozen_ch4_stock + lake_phase_transfer
				lake_phase_transfer = lake_water_o2_stock * lake_phase_fraction
				lake_water_o2_stock = lake_water_o2_stock - lake_phase_transfer
				lake_frozen_o2_stock = lake_frozen_o2_stock + lake_phase_transfer
			else if (lake_liquid_fraction_current < lake_liquid_fraction_prev) then
				lake_phase_fraction = (lake_liquid_fraction_prev - lake_liquid_fraction_current) / &
					lake_liquid_fraction_prev
				lake_phase_transfer = min(lake_water_ch4_stock, &
					lake_water_ch4_stock * lake_phase_fraction)
				lake_water_ch4_stock = lake_water_ch4_stock - lake_phase_transfer
				lake_frozen_ch4_stock = lake_frozen_ch4_stock + lake_phase_transfer
				lake_phase_transfer = min(lake_water_o2_stock, &
					lake_water_o2_stock * lake_phase_fraction)
				lake_water_o2_stock = lake_water_o2_stock - lake_phase_transfer
				lake_frozen_o2_stock = lake_frozen_o2_stock + lake_phase_transfer
			else if (lake_liquid_fraction_current > lake_liquid_fraction_prev) then
				lake_phase_fraction = (lake_liquid_fraction_current - lake_liquid_fraction_prev) / &
					(1._r8 - lake_liquid_fraction_prev)
				lake_phase_transfer = min(lake_frozen_ch4_stock, &
					lake_frozen_ch4_stock * lake_phase_fraction)
				lake_frozen_ch4_stock = lake_frozen_ch4_stock - lake_phase_transfer
				lake_water_ch4_stock = lake_water_ch4_stock + lake_phase_transfer
				lake_phase_transfer = min(lake_frozen_o2_stock, &
					lake_frozen_o2_stock * lake_phase_fraction)
				lake_frozen_o2_stock = lake_frozen_o2_stock - lake_phase_transfer
				lake_water_o2_stock = lake_water_o2_stock + lake_phase_transfer
			endif
			lake_liquid_fraction_prev = lake_liquid_fraction_current
		endif

		if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
			! LAKE_SOILC is initialized from the formal surface dataset field
			! lake_soilc_srf, then persisted in the Methane restart field.
			conc_o2_sat = conc_o2_lake
			conc_methane_sat = conc_methane_lake
			totcol_methane_sat = totcol_methane_lake
			grnd_methane_cond_sat = grnd_methane_cond_lake
		endif

		! Partition host soil water between conditional columns.  The saturated
		! column may borrow water only from the area represented by the unsaturated
		! column; their area-weighted liquid and ice volumes remain exactly equal to
		! the host state after the same porosity clipping the host applies before use
		! (MOD_SoilSnowHydrology.F90:312-314 -- the host clips vol_liq/vol_ice at
		! point of use and never writes the clipped values back, so the raw masses
		! can exceed pore volume once ice, at 917 kg/m3, fills a saturated layer).
		! Inundation unsupported by host water is invalid forcing, not permission
		! for the methane model to create water.
		wliq_soisno_unsat = wliq_soisno
		wice_soisno_unsat = wice_soisno
		wliq_soisno_sat = wliq_soisno
		wice_soisno_sat = wice_soisno
		if (patchtype /= 4 .and. finundated > 0._r8) then
			do j = 1, nl_soil
				pore_volume = max(porsl(j), 0._r8) * max(dz_soisno(j), 0._r8)
				if (pore_volume <= 1.e-12_r8) cycle
				vliq = max(wliq_soisno(j), 0._r8) / denh2o
				vice = max(wice_soisno(j), 0._r8) / denice
				vtot = vliq + vice
				! vtot comes from the host soil state, finundated from an
				! independent water-table S-curve, so the two disagree by more
				! than round-off at saturation. Tolerate that disagreement up to
				! DEF_METHANE%host_water_tolerance (a fraction of pore volume),
				! rescaling liquid, ice and total together so the host mass
				! balance this routine promises above is preserved; abort beyond
				! it, because a disagreement that large is invalid forcing, not
				! numerical drift. The guard this replaces compared at
				! 1.e-10*max(pore_volume,1) -- pore_volume is ~0.04 m so the max()
				! always selected 1.0, making it a flat 1e-10 m absolute
				! tolerance, i.e. round-off level between quantities that are not
				! computed from each other.
				host_water_tol = DEF_METHANE%host_water_tolerance * pore_volume

				! The two sides are NOT symmetric and must not share a tolerance.
				!
				! Deficit (vtot < finundated*pore_volume) is a legitimate
				! disagreement between two independently computed quantities, and
				! the capped allocation below reproduces the host exactly: the
				! inundated column is simply not fully saturated.

				if (vtot > pore_volume) then
					vliq = vliq * pore_volume / vtot
					vice = vice * pore_volume / vtot
					vtot = pore_volume
				endif
				if (vtot < finundated * pore_volume - host_water_tol) then
					write(6,*) 'ERROR: methane inundation exceeds host soil water: ', &
						ipatch, j, finundated, vtot, pore_volume, &
						' deficit/pore =', (finundated * pore_volume - vtot) / pore_volume, &
						' tolerance =', DEF_METHANE%host_water_tolerance
					CALL CoLM_stop ('methane inundation exceeds host soil water')
				endif

				! Inside the tolerance the host state is left untouched: vtot,
				! vliq and vice are never modified, so no host water is created
				! or discarded here.  What absorbs the disagreement is the
				! saturated-column allocation below, which is capped so the
				! area-weighted total reproduces the host exactly.
				! Record the worst residual and warn once, so that tolerating
				! the disagreement cannot hide one that is growing.
				host_water_resid = max(                                        &
					(vtot - pore_volume) / pore_volume,                         &
					(finundated * pore_volume - vtot) / pore_volume)
				if (host_water_resid > host_water_max_resid) then
					host_water_max_resid = host_water_resid
					if ((.not. host_water_warned) .and.                          &
					    host_water_resid > 0.1_r8 * DEF_METHANE%host_water_tolerance) then
						host_water_warned = .true.
						write(6,*) 'WARNING: methane host-water disagreement ',   &
							host_water_resid, ' of pore volume at patch/layer ',   &
							ipatch, j, ' (tolerance ',                            &
							DEF_METHANE%host_water_tolerance, '); tolerated'
					endif
				endif
				! Saturate the inundated column, but never beyond what the host
				! water can support.  Without the 1/finundated cap, a column
				! holding less than finundated*pore_volume gets a saturated
				! allocation the unsaturated side cannot pay for: its max(0,...)
				! floor stops the unsaturated column going negative but does NOT
				! undo the over-allocation, so the area-weighted total comes out
				! at finundated*pore_volume and the module creates exactly
				! (finundated*pore_volume - vtot) of water per layer.  That is
				! the precondition the old machine-epsilon guard enforced by
				! aborting; widening the tolerance without this cap silently
				! turned it into invented water.
				!
				! With the cap the area-weighted total is vliq (and vice)
				! exactly, in every branch: the inundated column is simply not
				! fully saturated when there is not enough water to saturate it.
				host_water_scale = min(pore_volume / max(vtot, 1.e-12_r8), &
				                       1._r8 / max(finundated, 1.e-12_r8))
				vliq_sat_alloc = vliq * host_water_scale
				vice_sat_alloc = vice * host_water_scale
				wliq_soisno_sat(j) = vliq_sat_alloc * denh2o
				wice_soisno_sat(j) = vice_sat_alloc * denice
				! a flooded soil patch takes its water from the flood, not
				! from the host soil: pores full, liquid or ice by the host's
				! phase share (by temperature when the host layer is empty), and
				! the unsaturated subcolumn keeps the host water as it is
				if (DEF_METHANE%flood_saturates_soil .and. patchtype == 0) then
					if (vliq + vice > 1.e-12_r8) then
						wliq_soisno_sat(j) = pore_volume * vliq / (vliq + vice) * denh2o
						wice_soisno_sat(j) = pore_volume * vice / (vliq + vice) * denice
					elseif (t_soisno(j) > tfrz) then
						wliq_soisno_sat(j) = pore_volume * denh2o
						wice_soisno_sat(j) = 0._r8
					else
						wliq_soisno_sat(j) = 0._r8
						wice_soisno_sat(j) = pore_volume * denice
					endif
					cycle
				endif
				if (finundated < 1._r8) then
					wliq_soisno_unsat(j) = max(0._r8, &
						(vliq - finundated * vliq_sat_alloc) / (1._r8 - finundated)) * denh2o
					wice_soisno_unsat(j) = max(0._r8, &
						(vice - finundated * vice_sat_alloc) / (1._r8 - finundated)) * denice
				endif
			enddo
		endif

      ! share of the plant conduit of each subcolumn that is woody and
      ! has the fixed cross-section woody_conduit_area: the forested share of
      ! a wetland tile in both subcolumns, all of the flooded subcolumn of a
      ! woody soil patch (flooded forest), or with pft_grass_share its
      ! tree cover; 0 when the key is off.
      woody_share_unsat = 0._r8
      woody_share_sat   = 0._r8
      if (DEF_METHANE%woody_conduit_area >= 0._r8) then
         if (patchtype == 2 .and. DEF_METHANE%wetland_veg_glwd .and. allocated(wetveg_forest)) then
            if (ipatch >= 1 .and. ipatch <= size(wetveg_forest)) &
               woody_share_unsat = min(max(wetveg_forest(ipatch), 0._r8), 1._r8)
            woody_share_sat = woody_share_unsat
         elseif (patchtype == 0 .and. .not. is_rice_paddy) then
            if (methane_patch_pft_shares(ipatch, herb_share_pft, woody_share_pft)) then
               woody_share_sat = woody_share_pft
            elseif (methane_patch_is_nongrass(patchclass)) then
               woody_share_sat = 1._r8
            endif
         endif
      endif

		!-------------------------------------------------
		! Loop
		!-------------------------------------------------
		do sat= 0, 1
			if (patchtype == 4 .and. DEF_METHANE%allowlakeprod .and. sat == 0) then
				methane_prod_depth_unsat = 0._r8
				methane_oxid_depth_unsat = 0._r8
				methane_ebul_depth_unsat = 0._r8
				methane_aere_depth_unsat = 0._r8
				methane_tran_depth_unsat = 0._r8
				o2_decomp_depth_unsat = 0._r8
				co2_decomp_depth_unsat = 0._r8
				o2_oxid_depth_unsat = 0._r8
				co2_oxid_depth_unsat = 0._r8
				o2_aere_depth_unsat = 0._r8
				co2_aere_depth_unsat = 0._r8
				o2stress_unsat = 1._r8
				methane_stress_unsat = 1._r8
					conc_methane_unsat = 0._r8
				totcol_methane_unsat = 0._r8
				methane_surf_aere_unsat = 0._r8
				methane_surf_ebul_unsat = 0._r8
				methane_surf_diff_unsat = 0._r8
				methane_surf_diff_phys_unsat = 0._r8
				methane_ebul_tot_unsat = 0._r8
				cycle
			endif
				if (sat==0) then ! unsaturated
					! Keep a prognostic, non-inundated wetland inventory and
					! disable production through methane_prod()'s j>jwt gate.
					! The host currently exposes only the wetland tile moisture,
					! so this branch retains that moisture as an explicit proxy;
					! it must not be described as a separately resolved dry soil.
					! Otherwise wetland zwt=0
					! hardcode (MOD_SoilSnowHydrology.F90:1229) gives
					! jwt_unsat=0, all layers produce, sat≈unsat → finundated
					! weighting becomes useless and ratio stuck at 1.5-1.7.
					if (DEF_METHANE%wetland_dry_unsat_branch .and. &
					    patchtype == 2) then
						zwt_unsat   = zi_soisno(nl_soil) + 1.e-3_r8   ! below deepest layer
						wdsrf_unsat = 0._r8                            ! no ponding
						jwt_unsat   = nl_soil                          ! no layer is below WT
					else
						zwt_unsat = zwt
						wdsrf_unsat = methane_wetland_water_depth(patchtype, wdsrf, wetwat)
						jwt_unsat = nl_soil
						! allow jwt to equal zero when zwt is in top layer
						do j = 1, nl_soil
							if(zwt_unsat <= zi_soisno(j)) then
								jwt_unsat = j-1
								exit
							end if
						end do
					endif
            ! (wt_layer_sat_share): layer jwt_unsat+1 holds the water
            ! table and counts as below it; only the share of its thickness
            ! under the table is saturated
            wt_share_unsat = 1._r8
            if (DEF_METHANE%wt_layer_sat_share .and. jwt_unsat < nl_soil) &
               wt_share_unsat = min(1._r8, max(0._r8, &
                  (zi_soisno(jwt_unsat+1) - zwt_unsat) / dz_soisno(jwt_unsat+1)))
            ! (wt_layer_node_transport): for gas transport, oxidation,
            ! ebullition and plant uptake the layer holding the table counts
            ! as saturated only when the table is above its node; below the
            ! node most of it is air-filled and open to the O2 from above.
            ! Production keeps jwt_unsat (with wt_layer_sat_share, on the
            ! saturated share)
            jwt_tran_unsat = jwt_unsat
            if (DEF_METHANE%wt_layer_node_transport .and. jwt_unsat < nl_soil) then
               if (zwt_unsat > z_soisno(jwt_unsat+1)) jwt_tran_unsat = jwt_unsat + 1
            endif

					if (methane_cold_start) then
						do j = 1, nl_soil
							vol_liq_init = max(0._r8, min(wliq_soisno_unsat(j) / &
								(max(dz_soisno(j), 1.e-12_r8) * denh2o), porsl(j)))
							vol_ice_init = max(0._r8, min(wice_soisno_unsat(j) / &
								(max(dz_soisno(j), 1.e-12_r8) * denice), max(porsl(j) - vol_liq_init, 0._r8)))
							vol_gas_init = max(porsl(j) - vol_liq_init - vol_ice_init, 0._r8)
							conc_methane_unsat(j) = c_atm(1) * &
								(vol_gas_init + k_h_cc(j,1) * vol_liq_init)
							conc_o2_unsat(j) = c_atm(2) * &
								(vol_gas_init + k_h_cc(j,2) * vol_liq_init)
						enddo
					endif

				do j=1,nl_soil
					if (DEF_METHANE%use_vertical_redoxlag .and. j > jwt_unsat .and. redoxlags_vertical > 0._r8) then ! saturated currently
						layer_sat_lag(j) = layer_sat_lag(j) * exp(-deltim/redoxlags_vertical) &
							+ (1._r8 - exp(-deltim/redoxlags_vertical))
					else if (DEF_METHANE%use_vertical_redoxlag .and. redoxlags_vertical > 0._r8) then
						layer_sat_lag(j) = layer_sat_lag(j) * exp(-deltim/redoxlags_vertical)
					else if (j > jwt_unsat) then  ! redoxlags_vertical = 0
						layer_sat_lag(j) = 1._r8
					else
						layer_sat_lag(j) = 0._r8
					end if
				end do

				call split_ch4_o2_phases( dz_soisno, wliq_soisno_unsat, wice_soisno_unsat, porsl, &
					conc_methane_unsat, conc_o2_unsat, k_h_cc, idate, &
					vol_aqu_unsat, vol_ch4_storage_unsat, vol_gas_unsat, f_aqu_unsat, f_gas_unsat, &
					conc_ch4_gas_unsat, conc_ch4_aqu_unsat, conc_ch4_porsl_unsat, conc_ch4_gas_porsl_unsat, conc_ch4_aqu_porsl_unsat, &
					conc_o2_gas_unsat, conc_o2_aqu_unsat, conc_o2_porsl_unsat, conc_o2_gas_porsl_unsat, conc_o2_aqu_porsl_unsat )

            ! Unfrozen water fraction per layer; a layer without water is unconstrained.
            do j = 1, nl_soil
               if (wliq_soisno_unsat(j) + wice_soisno_unsat(j) > 1.e-12_r8) then
                  fliq_unsat(j) = max(wliq_soisno_unsat(j), 0._r8) / (wliq_soisno_unsat(j) + wice_soisno_unsat(j))
               else
                  fliq_unsat(j) = 1._r8
               end if
            end do
				! Calculate CH4 production in each soil layer
				call methane_prod ( ipatch, idate, patchtype, sat, jwt_unsat, finundated, finundated_lag, rr, deltim, &
               z_soisno, dz_soisno, zi_soisno, t_soisno, fliq_unsat, &
					lai, conc_o2_unsat, rootfr, annavg_finrw, &
					crootfr, somhr, lithr, hr_vr, o_scalar, fphr, pot_f_nit_vr, pH, layer_sat_lag, lake_soilc, &
               microbial_prod_potential_patch, anoxic_decomp_patch, wt_share_unsat, &
					methane_prod_depth_unsat, o2_decomp_depth_unsat, co2_decomp_depth_unsat )
            ! alternative electron acceptors take their share of
            ! the anaerobic carbon flow below the water table
            if (use_acceptor) call methane_acceptor ( ipatch, sat, jwt_unsat, deltim, fliq_unsat, &
               cellorg, organic_max, conc_o2_aqu_porsl_unsat, &
               methane_prod_depth_unsat, o2_decomp_depth_unsat, co2_decomp_depth_unsat, acc_reox_unsat )

				! Calculate CH4 oxidation in each soil layer
            call methane_oxid ( idate, patchtype, jwt_tran_unsat, sat, t_soisno, fliq_unsat, dz_soisno, zi_soisno, smp, vol_aqu_unsat, &
					conc_o2_aqu_porsl_unsat, conc_ch4_aqu_porsl_unsat, &
					microbial_oxid_potential_patch, &
					methane_oxid_depth_unsat, o2_oxid_depth_unsat )

				! Calculate CH4 ebullition losses in each soil layer
            call methane_ebul ( idate, patchtype, jwt_tran_unsat, sat, finundated, deltim, &
					z_soisno, dz_soisno, zi_soisno, forc_pbot, lake_depth_current, lake_icefrac, &
					t_soisno, wdsrf_unsat, conc_methane_unsat, conc_ch4_gas_porsl_unsat, &
					methane_ebul_depth_unsat )

				! Calculate CH4 aerenchyma losses in each soil layer
				call methane_aere ( ipatch, idate, jwt_tran_unsat, sat, patchclass, lai, woody_share_unsat, &
               z_soisno, dz_soisno, zi_soisno, t_soisno, fliq_unsat, &
					rootfr, rootr, etr, grnd_methane_cond_base, c_atm, annsum_npp, &
					annavg_agnpp, annavg_bgnpp, &
					conc_ch4_aqu_unsat, conc_ch4_aqu_porsl_unsat, conc_ch4_gas_porsl_unsat, conc_o2_aqu_porsl_unsat, conc_o2_gas_porsl_unsat, &
					methane_aere_depth_unsat, methane_tran_depth_unsat, o2_aere_depth_unsat )
            ! a wetland tile with a fixed rhizosphere oxidation share
            ! takes no plant O2 into its soil layers (the roots use it themselves)
            IF (patchtype == 2 .and. DEF_METHANE%wetland_aereoxid >= 0._r8) o2_aere_depth_unsat = 0._r8
            IF (is_rice_paddy .and. DEF_METHANE%rice_aereoxid > 0._r8 .and. DEF_METHANE%rice_rox_only) &
               o2_aere_depth_unsat = 0._r8
            ! rox_on_rooted_production: the rhizosphere share of the CH4 made
            ! below the table in the rooted layers never enters the pore water
            prod_tran_unsat(:) = methane_prod_depth_unsat(:)
            rox_prod_unsat(:)  = 0._r8
            IF (DEF_METHANE%rox_on_rooted_production) &
               call methane_root_oxid ( patchtype, sat, jwt_unsat, is_rice_paddy, lai, dz_soisno(1:nl_soil), rootfr, &
                  methane_prod_depth_unsat, prod_tran_unsat, rox_prod_unsat )

				! Solve CH4 reaction/diffusion equation
				! Competition for oxygen will occur here.
				call methane_tran ( idate, patchtype, &
               lb, snl, jwt_tran_unsat, sat, finundated, &
					dlon, dlat, deltim, z_soisno, dz_soisno, zi_soisno, t_soisno, t_grnd, forc_us, forc_vs, &
					porsl, wliq_soisno_unsat, wice_soisno_unsat, wdsrf_unsat, lakedepth, dz_lake, t_lake, lake_icefrac, bsw, c_atm, prod_tran_unsat, o2_aere_depth_unsat, &
					cellorg, t_h2osfc, organic_max, k_h_cc, conc_ch4_gas_porsl_unsat, conc_ch4_aqu_porsl_unsat, conc_o2_gas_porsl_unsat, conc_o2_aqu_porsl_unsat, &
					vol_aqu_unsat, vol_ch4_storage_unsat, vol_gas_unsat, &
					o2stress_unsat, methane_stress_unsat, methane_surf_aere_unsat, methane_surf_ebul_unsat, methane_surf_diff_unsat, methane_ebul_tot_unsat, &
					methane_oxid_depth_unsat, methane_aere_depth_unsat, methane_tran_depth_unsat, methane_ebul_depth_unsat, &
						grnd_methane_cond_base, grnd_methane_cond_unsat, &
						o2_oxid_depth_unsat, o2_decomp_depth_unsat, &
						conc_o2_unsat, conc_methane_unsat, methane_balance_residual_unsat, methane_ch4_clip_credit_unsat, methane_surf_diff_phys_unsat, &
						o2_cap_loss_unsat, o2_cap_gain_unsat, lake_water_ch4_stock, lake_water_o2_stock, &
						lake_water_ch4_oxid, lake_sed_ch4_flux, lake_sed_o2_flux, lake_air_o2_flux, &
						lake_icebub_stock, &
						rice_rox_in = is_rice_paddy .and. DEF_METHANE%rice_aereoxid > 0._r8 )
            ! rox_on_rooted_production: counted as oxidation, with plant O2
            IF (DEF_METHANE%rox_on_rooted_production) &
               methane_oxid_depth_unsat(:) = methane_oxid_depth_unsat(:) + rox_prod_unsat(:)
            ! the reoxidation the O2 competition allowed, at the
            ! O2 stress that scaled its O2 demand
            if (use_acceptor) acceptor_unsat(:,ipatch) = min(1._r8, acceptor_unsat(:,ipatch) + &
               acc_reox_unsat(:) * o2stress_unsat(:) * deltim)

			elseif (sat==1) then ! saturated
				zwt_sat=0.

				! IF (.not.DEF_SPLIT_SOILSNOW) THEN
				! 	IF (lb >= 1) THEN
				! 		wetwat_sat = wdsrf_sat + wa_sat + wetwat_sat + (gwat - etr + qsdew + qfros - qsubl) * deltim
				! 	ELSE
				! 		wetwat_sat = wdsrf_sat + wa_sat + wetwat_sat + (gwat - etr) * deltim
				! 	ENDIF
				! ELSE
				! 	wetwat_sat = wdsrf_sat + wa_sat + wetwat_sat + (gwat - etr + qsdew_soil + qfros_soil - qsubl_soil) * deltim
				! ENDIF

				! IF (wetwat_sat > wetwatmax) THEN
				! 	wdsrf_sat  = wetwat_sat - wetwatmax
				! 	wetwat_sat = wetwatmax
	! 	wa_sat     = 0.
				! ELSEIF (wetwat_sat < 0) THEN
	            ! 	wa_sat     = wetwat_sat
				! 	wdsrf_sat  = 0.
				! 	wetwat_sat = 0.
				! ELSE
				! 	wdsrf_sat = 0.
		        !     wa_sat    = 0.
				! ENDIF
				! Lake sediment is saturated by the overlying host lake water.
				if (patchtype == 4) then
					do j = 1, nl_soil
						vliq = max(wliq_soisno(j), 0._r8) / denh2o
						vice = max(wice_soisno(j), 0._r8) / denice
						vtot = vliq + vice
						if (vtot > 1.e-12_r8) then
							wliq_soisno_sat(j) = porsl(j)*dz_soisno(j)*denh2o * vliq / vtot
							wice_soisno_sat(j) = porsl(j)*dz_soisno(j)*denice * vice / vtot
						elseif (t_soisno(j) > tfrz) then
							wliq_soisno_sat(j) = porsl(j)*dz_soisno(j)*denh2o
							wice_soisno_sat(j) = 0._r8
						else
							wliq_soisno_sat(j) = 0._r8
							wice_soisno_sat(j) = porsl(j)*dz_soisno(j)*denice
						endif
					enddo
				endif
					! Surface water depth must come from host hydrology or routing.
					! An inundated-area product alone does not supply water mass.
					wdsrf_sat = methane_wetland_water_depth(patchtype, wdsrf, wetwat)
					! Apply routing-published floodplain depth on patches that
					! also consumed routing-published fldfrc above:
					!   (a) scheme 7 routing: all patches
					!   (b) scheme 6 + use_routing_for_soil (hybrid): soil tiles only
            !   (c) scheme 8 colm: every tile except the saturated wetland
					! Wetland tiles in hybrid stay on sigmoid(zwt), so they don't
					! get the routing depth (would create area/depth mismatch).
					IF (finundated > 0._r8 .and. &
					    allocated(f_inund_flood_depth_patch) .and. &
					    ipatch >= 1 .and. ipatch <= size(f_inund_flood_depth_patch) .and. &
					    (DEF_wetland_finundation_scheme == 7 .or. &
               (DEF_wetland_finundation_scheme == 8 .and. patchtype /= 2) .or. &
					     (DEF_wetland_finundation_scheme == 6 .and. &
					      DEF_METHANE%use_routing_for_soil .and. patchtype /= 2))) THEN
					   routing_depth_mm = 1000._r8 * max(0._r8, f_inund_flood_depth_patch(ipatch))
					   wdsrf_sat = max(wdsrf_sat, max(finundated, 0.01_r8) * routing_depth_mm)
					ENDIF
					jwt_sat = 0

				if (methane_cold_start .and. patchtype /= 4) then
					do j = 1, nl_soil
						vol_liq_init = max(0._r8, min(wliq_soisno_sat(j) / &
							(max(dz_soisno(j), 1.e-12_r8) * denh2o), porsl(j)))
						vol_ice_init = max(0._r8, min(wice_soisno_sat(j) / &
							(max(dz_soisno(j), 1.e-12_r8) * denice), max(porsl(j) - vol_liq_init, 0._r8)))
						vol_gas_init = max(porsl(j) - vol_liq_init - vol_ice_init, 0._r8)
						conc_methane_sat(j) = c_atm(1) * &
							(vol_gas_init + k_h_cc(j,1) * vol_liq_init)
						conc_o2_sat(j) = c_atm(2) * &
							(vol_gas_init + k_h_cc(j,2) * vol_liq_init)
					enddo
				endif

				call split_ch4_o2_phases( dz_soisno, wliq_soisno_sat, wice_soisno_sat, porsl, &
					conc_methane_sat, conc_o2_sat, k_h_cc, idate, &
					vol_aqu_sat, vol_ch4_storage_sat, vol_gas_sat, f_aqu_sat, f_gas_sat, &
					conc_ch4_gas_sat, conc_ch4_aqu_sat, conc_ch4_porsl_sat, conc_ch4_gas_porsl_sat, conc_ch4_aqu_porsl_sat, &
					conc_o2_gas_sat, conc_o2_aqu_sat, conc_o2_porsl_sat, conc_o2_gas_porsl_sat, conc_o2_aqu_porsl_sat )

            ! Unfrozen water fraction per layer; a layer without water is unconstrained.
            do j = 1, nl_soil
               if (wliq_soisno_sat(j) + wice_soisno_sat(j) > 1.e-12_r8) then
                  fliq_sat(j) = max(wliq_soisno_sat(j), 0._r8) / (wliq_soisno_sat(j) + wice_soisno_sat(j))
               else
                  fliq_sat(j) = 1._r8
               end if
            end do
				! Calculate CH4 production in each soil layer
				call methane_prod ( ipatch, idate, patchtype, sat, jwt_sat, finundated, finundated_lag, rr, deltim, &
               z_soisno, dz_soisno, zi_soisno, t_soisno, fliq_sat, &
					lai, conc_o2_sat, rootfr, annavg_finrw, &
					crootfr, somhr, lithr, hr_vr, o_scalar, fphr, pot_f_nit_vr, pH, layer_sat_lag, lake_soilc, &
               microbial_prod_potential_patch, anoxic_decomp_patch, 1._r8, &
					methane_prod_depth_sat, o2_decomp_depth_sat, co2_decomp_depth_sat )
            ! alternative electron acceptors take their share of
            ! the anaerobic carbon flow
            if (use_acceptor) call methane_acceptor ( ipatch, sat, jwt_sat, deltim, fliq_sat, &
               cellorg, organic_max, conc_o2_aqu_porsl_sat, &
               methane_prod_depth_sat, o2_decomp_depth_sat, co2_decomp_depth_sat, acc_reox_sat )

				! Calculate CH4 oxidation in each soil layer
            call methane_oxid ( idate, patchtype, jwt_sat, sat, t_soisno, fliq_sat, dz_soisno, zi_soisno, smp, vol_aqu_sat, &
					conc_o2_aqu_porsl_sat, conc_ch4_aqu_porsl_sat, &
					microbial_oxid_potential_patch, &
					methane_oxid_depth_sat, o2_oxid_depth_sat )

				! Calculate CH4 ebullition losses in each soil layer
				call methane_ebul ( idate, patchtype, jwt_sat, sat, finundated, deltim, &
					z_soisno, dz_soisno, zi_soisno, forc_pbot, lake_depth_current, lake_icefrac, &
					t_soisno, wdsrf_sat, conc_methane_sat, conc_ch4_gas_porsl_sat, &
					methane_ebul_depth_sat )

				! Calculate CH4 aerenchyma losses in each soil layer
				if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
					methane_aere_depth_sat = 0._r8
					methane_tran_depth_sat = 0._r8
					o2_aere_depth_sat = 0._r8
				else
					call methane_aere ( ipatch, idate, jwt_sat, sat, patchclass, lai, woody_share_sat, &
                  z_soisno, dz_soisno, zi_soisno, t_soisno, fliq_sat, &
						rootfr, rootr, etr, grnd_methane_cond_base, c_atm, annsum_npp, &
						annavg_agnpp, annavg_bgnpp, &
						conc_ch4_aqu_sat, conc_ch4_aqu_porsl_sat, conc_ch4_gas_porsl_sat, conc_o2_aqu_porsl_sat, conc_o2_gas_porsl_sat, &
						methane_aere_depth_sat, methane_tran_depth_sat, o2_aere_depth_sat )
               ! a wetland tile, or the flooded subcolumn of a soil
               ! patch other than a paddy, with a fixed rhizosphere oxidation
               ! share takes no plant O2 into its soil layers (the roots use it)
               IF ((patchtype == 2 .or. (patchtype == 0 .and. .not. is_rice_paddy)) &
                   .and. DEF_METHANE%wetland_aereoxid >= 0._r8) o2_aere_depth_sat = 0._r8
               IF (is_rice_paddy .and. DEF_METHANE%rice_aereoxid > 0._r8 .and. DEF_METHANE%rice_rox_only) &
                  o2_aere_depth_sat = 0._r8
				endif
            ! rox_on_rooted_production (see the unsaturated branch)
            prod_tran_sat(:) = methane_prod_depth_sat(:)
            rox_prod_sat(:)  = 0._r8
            IF (DEF_METHANE%rox_on_rooted_production) &
               call methane_root_oxid ( patchtype, sat, jwt_sat, is_rice_paddy, lai, dz_soisno(1:nl_soil), rootfr, &
                  methane_prod_depth_sat, prod_tran_sat, rox_prod_sat )

				! Solve CH4 reaction/diffusion equation
				! Competition for oxygen will occur here.
				call methane_tran ( idate, patchtype, &
						lb, snl, jwt_sat, sat, finundated, &
						dlon, dlat, deltim, z_soisno, dz_soisno, zi_soisno, t_soisno, t_grnd, forc_us, forc_vs, &
						porsl, wliq_soisno_sat, wice_soisno_sat, wdsrf_sat, lakedepth, dz_lake, t_lake, lake_icefrac, bsw, c_atm, prod_tran_sat, o2_aere_depth_sat, &
					cellorg, t_h2osfc, organic_max, k_h_cc, conc_ch4_gas_porsl_sat, conc_ch4_aqu_porsl_sat, conc_o2_gas_porsl_sat, conc_o2_aqu_porsl_sat, &
					vol_aqu_sat, vol_ch4_storage_sat, vol_gas_sat, &
					o2stress_sat, methane_stress_sat, methane_surf_aere_sat, methane_surf_ebul_sat, methane_surf_diff_sat, methane_ebul_tot_sat, &
					methane_oxid_depth_sat, methane_aere_depth_sat, methane_tran_depth_sat, methane_ebul_depth_sat, &
						grnd_methane_cond_base, grnd_methane_cond_sat, &
						o2_oxid_depth_sat, o2_decomp_depth_sat, &
						conc_o2_sat, conc_methane_sat, methane_balance_residual_sat, methane_ch4_clip_credit_sat, methane_surf_diff_phys_sat, &
						o2_cap_loss_sat, o2_cap_gain_sat, lake_water_ch4_stock, lake_water_o2_stock, &
						lake_water_ch4_oxid, lake_sed_ch4_flux, lake_sed_o2_flux, lake_air_o2_flux, &
						lake_icebub_stock, &
						rice_rox_in = is_rice_paddy .and. DEF_METHANE%rice_aereoxid > 0._r8 )
            ! rox_on_rooted_production: counted as oxidation, with plant O2
            IF (DEF_METHANE%rox_on_rooted_production) &
               methane_oxid_depth_sat(:) = methane_oxid_depth_sat(:) + rox_prod_sat(:)
            ! the reoxidation the O2 competition allowed
            if (use_acceptor) acceptor_sat(:,ipatch) = min(1._r8, acceptor_sat(:,ipatch) + &
               acc_reox_sat(:) * o2stress_sat(:) * deltim)

			endif

		enddo

		! Diagnostic CO2 closure.  methane_prod computes decomposition/root-respiration
		! CO2 from the same base decomposition used for CH4 production, excluding
		! nitrification O2 demand and excluding any anoxia-inflated potential O2 demand.
		! co2_decomp_depth is intentionally not rescaled by the O2 competition
		! limiter in methane_tran; treat it as a carbon-source diagnostic rather
		! than a prognostic O2-limited CO2 tracer.  CH4 oxidation produces CO2 at
		! 1:1 molar ratio. Explicit CO2 aerenchyma transport is not prognosed
		! here, so it is diagnosed as zero.
		co2_oxid_depth_unsat   = methane_oxid_depth_unsat
		co2_oxid_depth_sat     = methane_oxid_depth_sat
		co2_aere_depth_unsat   = 0._r8
		co2_aere_depth_sat     = 0._r8

		methane_prod_depth = methane_prod_depth_sat * finundated + methane_prod_depth_unsat * (1.0_r8 - finundated)
		methane_oxid_depth = methane_oxid_depth_sat * finundated + methane_oxid_depth_unsat * (1.0_r8 - finundated)
		methane_aere_depth = methane_aere_depth_sat * finundated + methane_aere_depth_unsat * (1.0_r8 - finundated)
		methane_ebul_depth = methane_ebul_depth_sat * finundated + methane_ebul_depth_unsat * (1.0_r8 - finundated)
		methane_tran_depth = methane_tran_depth_sat * finundated + methane_tran_depth_unsat * (1.0_r8 - finundated)
		co2_decomp_depth = co2_decomp_depth_sat * finundated + co2_decomp_depth_unsat * (1.0_r8 - finundated)
		co2_oxid_depth   = co2_oxid_depth_sat   * finundated + co2_oxid_depth_unsat   * (1.0_r8 - finundated)
		co2_aere_depth   = co2_aere_depth_sat   * finundated + co2_aere_depth_unsat   * (1.0_r8 - finundated)

		methane_oxid_tot_unsat = sum( methane_oxid_depth_unsat(1:nl_soil) * dz_soisno(1:nl_soil) )
		methane_prod_tot_unsat = sum( methane_prod_depth_unsat(1:nl_soil) * dz_soisno(1:nl_soil) )
		co2_decomp_tot_unsat = sum( co2_decomp_depth_unsat(1:nl_soil) * dz_soisno(1:nl_soil) )
		co2_oxid_tot_unsat   = sum( co2_oxid_depth_unsat  (1:nl_soil) * dz_soisno(1:nl_soil) )
		co2_net_tot_unsat    = co2_decomp_tot_unsat + co2_oxid_tot_unsat
		net_methane_unsat = -methane_prod_tot_unsat + methane_oxid_tot_unsat
		methane_surf_flux_tot_unsat = methane_surf_diff_unsat + methane_surf_aere_unsat + methane_surf_ebul_unsat + &
			sum(methane_tran_depth_unsat(1:nl_soil) * dz_soisno(1:nl_soil))
		totcol_methane_unsat = sum( conc_methane_unsat(1:nl_soil) * dz_soisno(1:nl_soil) )

		methane_oxid_tot_sat   = sum( methane_oxid_depth_sat(1:nl_soil)   * dz_soisno(1:nl_soil) )
		methane_prod_tot_sat   = sum( methane_prod_depth_sat(1:nl_soil)   * dz_soisno(1:nl_soil) )
		co2_decomp_tot_sat = sum( co2_decomp_depth_sat(1:nl_soil) * dz_soisno(1:nl_soil) )
		co2_oxid_tot_sat   = sum( co2_oxid_depth_sat  (1:nl_soil) * dz_soisno(1:nl_soil) )
		if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
			methane_oxid_tot_sat = methane_oxid_tot_sat + lake_water_ch4_oxid
			co2_oxid_tot_sat = co2_oxid_tot_sat + lake_water_ch4_oxid
		endif
		co2_net_tot_sat    = co2_decomp_tot_sat + co2_oxid_tot_sat
		net_methane_sat   = -methane_prod_tot_sat   + methane_oxid_tot_sat
		methane_surf_flux_tot_sat = methane_surf_diff_sat  + methane_surf_aere_sat  + methane_surf_ebul_sat + &
			sum(methane_tran_depth_sat(1:nl_soil) * dz_soisno(1:nl_soil))
		totcol_methane_sat   = sum( conc_methane_sat  (1:nl_soil) * dz_soisno(1:nl_soil) )

			methane_surf_diff = methane_surf_diff_sat * finundated + methane_surf_diff_unsat * (1.0_r8 - finundated) + &
				methane_dfsat_tot
			methane_surf_diff_phys = methane_surf_diff_phys_sat * finundated + &
				methane_surf_diff_phys_unsat * (1.0_r8 - finundated)
			methane_balance_residual = methane_balance_residual_sat * finundated + &
				methane_balance_residual_unsat * (1.0_r8 - finundated)
			methane_ch4_clip_credit = methane_ch4_clip_credit_sat * finundated + &
				methane_ch4_clip_credit_unsat * (1.0_r8 - finundated)
			o2_cap_loss = o2_cap_loss_sat * finundated + o2_cap_loss_unsat * (1.0_r8 - finundated)
			o2_cap_gain = o2_cap_gain_sat * finundated + o2_cap_gain_unsat * (1.0_r8 - finundated)
			methane_surf_ebul = methane_surf_ebul_sat * finundated + methane_surf_ebul_unsat * (1.0_r8 - finundated)
		methane_surf_aere = methane_surf_aere_sat * finundated + methane_surf_aere_unsat * (1.0_r8 - finundated)
		methane_oxid_tot = methane_oxid_tot_sat * finundated + methane_oxid_tot_unsat * (1.0_r8 - finundated)
		methane_prod_tot = methane_prod_tot_sat * finundated + methane_prod_tot_unsat * (1.0_r8 - finundated)
		co2_decomp_tot = co2_decomp_tot_sat * finundated + co2_decomp_tot_unsat * (1.0_r8 - finundated)
		co2_oxid_tot   = co2_oxid_tot_sat   * finundated + co2_oxid_tot_unsat   * (1.0_r8 - finundated)
		co2_aere_tot   = sum( co2_aere_depth(1:nl_soil) * dz_soisno(1:nl_soil) )
		co2_net_tot    = co2_decomp_tot + co2_oxid_tot - co2_aere_tot
		net_methane = net_methane_sat * finundated + net_methane_unsat * (1.0_r8 - finundated)
		methane_surf_flux_tot = methane_surf_diff + methane_surf_ebul + methane_surf_aere + &
			sum(methane_tran_depth(1:nl_soil) * dz_soisno(1:nl_soil))
		methane_surf_flux_tot_phys = methane_surf_diff_phys + methane_dfsat_tot + &
			methane_surf_ebul + methane_surf_aere + &
			sum(methane_tran_depth(1:nl_soil) * dz_soisno(1:nl_soil))
		totcol_methane = totcol_methane_sat * finundated + totcol_methane_unsat * (1.0_r8 - finundated)
		grnd_methane_cond = grnd_methane_cond_sat * finundated + &
			grnd_methane_cond_unsat * (1.0_r8 - finundated)

		conc_methane = conc_methane_sat * finundated + conc_methane_unsat * (1._r8 - finundated)
		conc_o2 = conc_o2_sat * finundated + conc_o2_unsat * (1._r8 - finundated)

      ! the two subcolumns stand for the low share of a tower's
      ! wetland tile. The high share respires as usual but makes no CH4, so it
      ! has none to oxidise, store or release: the tile's CH4 fluxes and store
      ! are the low share's times high_keep, and the carbon the high share
      ! would have made into CH4 leaves as CO2 (CO2 plus CH4 production is the
      ! decomposed carbon plus root respiration in either share). The
      ! subcolumn fields and the O2 fields describe the low share.
      if (high_keep < 1._r8) then
         co2_decomp_depth   = co2_decomp_depth + (1._r8 - high_keep) * methane_prod_depth
         co2_decomp_tot     = co2_decomp_tot   + (1._r8 - high_keep) * methane_prod_tot
         methane_prod_depth = methane_prod_depth * high_keep
         methane_oxid_depth = methane_oxid_depth * high_keep
         methane_aere_depth = methane_aere_depth * high_keep
         methane_ebul_depth = methane_ebul_depth * high_keep
         methane_tran_depth = methane_tran_depth * high_keep
         co2_oxid_depth     = co2_oxid_depth     * high_keep
         co2_aere_depth     = co2_aere_depth     * high_keep
         methane_surf_diff        = methane_surf_diff        * high_keep
         methane_surf_diff_phys   = methane_surf_diff_phys   * high_keep
         methane_balance_residual = methane_balance_residual * high_keep
         methane_ch4_clip_credit  = methane_ch4_clip_credit  * high_keep
         methane_surf_ebul  = methane_surf_ebul  * high_keep
         methane_surf_aere  = methane_surf_aere  * high_keep
         methane_oxid_tot   = methane_oxid_tot   * high_keep
         methane_prod_tot   = methane_prod_tot   * high_keep
         co2_oxid_tot       = co2_oxid_tot       * high_keep
         co2_aere_tot       = co2_aere_tot       * high_keep
         co2_net_tot        = co2_decomp_tot + co2_oxid_tot - co2_aere_tot
         net_methane        = net_methane        * high_keep
         methane_surf_flux_tot      = methane_surf_flux_tot      * high_keep
         methane_surf_flux_tot_phys = methane_surf_flux_tot_phys * high_keep
         totcol_methane     = totcol_methane     * high_keep
         conc_methane       = conc_methane       * high_keep
      endif

			if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
				methane_prod_depth_lake = methane_prod_depth_sat
			methane_oxid_depth_lake = methane_oxid_depth_sat
			methane_ebul_depth_lake = methane_ebul_depth_sat
			co2_decomp_depth_lake = co2_decomp_depth_sat
			co2_oxid_depth_lake = methane_oxid_depth_lake
			methane_surf_ebul_lake = methane_surf_ebul_sat
			methane_surf_diff_lake = methane_surf_diff_sat
			methane_surf_flux_tot_lake = methane_surf_ebul_lake + methane_surf_diff_lake
			methane_prod_tot_lake = methane_prod_tot_sat
			methane_oxid_tot_lake = methane_oxid_tot_sat
			methane_ebul_tot_lake = methane_ebul_tot_sat
			co2_decomp_tot_lake = co2_decomp_tot_sat
			co2_oxid_tot_lake = co2_oxid_tot_sat
			co2_net_tot_lake = co2_net_tot_sat
			totcol_methane_lake = totcol_methane_sat
			grnd_methane_cond_lake = grnd_methane_cond_sat
			grnd_methane_cond = grnd_methane_cond_lake
			conc_o2_lake = conc_o2_sat
			conc_methane_lake = conc_methane_sat
			methane_prod_depth = methane_prod_depth_lake
			methane_oxid_depth = methane_oxid_depth_lake
			methane_ebul_depth = methane_ebul_depth_lake
			co2_decomp_depth = co2_decomp_depth_lake
			co2_oxid_depth = co2_oxid_depth_lake
			methane_surf_ebul = methane_surf_ebul_lake
			methane_surf_diff = methane_surf_diff_lake
			methane_surf_aere = 0._r8
			methane_surf_flux_tot = methane_surf_flux_tot_lake
			methane_surf_flux_tot_phys = methane_surf_diff_phys + methane_surf_ebul
			methane_prod_tot = methane_prod_tot_lake
			methane_oxid_tot = methane_oxid_tot_lake
			methane_ebul_tot = methane_ebul_tot_lake
			co2_decomp_tot = co2_decomp_tot_lake
			co2_oxid_tot = co2_oxid_tot_lake
			co2_net_tot = co2_net_tot_lake
			totcol_methane = totcol_methane_lake + lake_water_ch4_stock + lake_frozen_ch4_stock
			! CH4 held in bubbles under the ice belongs to the lake inventory
			if (lake_icebub_stock > 0._r8) totcol_methane = totcol_methane + lake_icebub_stock
			if (lake_icebub_state) lake_icebubble_ch4_stock(ipatch) = lake_icebub_stock
			net_methane = net_methane_sat
			methane_balance_residual = methane_balance_residual_sat
			methane_ch4_clip_credit = methane_ch4_clip_credit_sat
			o2_cap_loss = o2_cap_loss_sat
			o2_cap_gain = o2_cap_gain_sat
			conc_o2 = conc_o2_lake
			conc_methane = conc_methane_lake
			if (lake_restart_debug_print) then
				write(6,'(A,1X,I6,1X,I4,1X,I7,1X,I9,1X,2F10.3,1X,19ES14.6)') &
					'CH4_LAKE_RESTART_DEBUG', idate(1), idate(2), idate(3), ipatch, dlon, dlat, &
					lakedepth, wdsrf, f_h2osfc, lake_dbg_totcol_bef, totcol_methane_lake, &
					methane_prod_tot_lake, methane_oxid_tot_lake, methane_surf_diff_lake, &
					methane_surf_ebul_lake, methane_surf_flux_tot_lake, lake_dbg_ch4_l1_bef, &
					lake_dbg_o2_l1_bef, conc_methane_lake(1), conc_o2_lake(1), &
					lake_dbg_ch4_l5_bef, lake_dbg_o2_l5_bef, conc_methane_lake(min(5,nl_soil)), &
					conc_o2_lake(min(5,nl_soil)), sum(lake_soilc(1:nl_soil))
			endif
         ! the active pool, not lake_soilc, pays for the
         ! decomposition when it is on
         if (.not. DEF_METHANE%replenishlakec .and. .not. DEF_METHANE%lake_active_pool) then
            if (DEF_METHANE%lake_prod_efold > 0._r8) then
               ! the column total pays for the decomposition, spread
               ! over the layers as the carbon is
               lake_loss = sum((methane_prod_depth_sat(1:nl_soil) + co2_decomp_depth_sat(1:nl_soil)) &
                  * dz_soisno(1:nl_soil)) * deltim * catomw
               lake_cc = sum(max(lake_soilc(1:nl_soil), 0._r8) * dz_soisno(1:nl_soil))
               if (lake_cc > 0._r8) lake_soilc(1:nl_soil) = max(lake_soilc(1:nl_soil), 0._r8) &
                  * max(0._r8, 1._r8 - lake_loss / lake_cc)
            else
				do j = 1, nl_soil
					lake_soilc(j) = max(0._r8, lake_soilc(j) - &
						(methane_prod_depth_sat(j) + co2_decomp_depth_sat(j)) * deltim * catomw)
				enddo
				endif
			endif
         ! dC_act/dt = S - k(T) C_act. The pool gains the supply
         ! of the patch's latitude band and loses what methane_prod decomposed
         ! from it this step (CH4 plus CO2 carbon), so the carbon it books is
         ! the carbon the layers used.
         if (DEF_METHANE%lake_active_pool .and. allocated(lake_cact)) then
            if (ipatch >= 1 .and. ipatch <= size(lake_cact)) then
               if (dlat > 60._r8) then
                  lake_cdep = DEF_METHANE%lake_cdep_band(1)
               elseif (dlat > 45._r8) then
                  lake_cdep = DEF_METHANE%lake_cdep_band(2)
               elseif (dlat > 23.5_r8) then
                  lake_cdep = DEF_METHANE%lake_cdep_band(3)
               elseif (dlat >= -23.5_r8) then
                  lake_cdep = DEF_METHANE%lake_cdep_band(4)
               else
                  lake_cdep = DEF_METHANE%lake_cdep_band(5)
               endif
               lake_loss = sum((methane_prod_depth_sat(1:nl_soil) + co2_decomp_depth_sat(1:nl_soil)) &
                  * dz_soisno(1:nl_soil)) * deltim * catomw
               ! [gC m-2] = [gC m-2 yr-1] * [s] / [s yr-1] (365-day year)
               lake_cact(ipatch) = max(0._r8, lake_cact(ipatch) + &
                  lake_cdep * deltim / (365._r8 * secspday) - lake_loss)
            endif
         endif
      endif

			! Patch-level diagnostics for history output.  These expose the actual
			! inundation signal consumed by the methane physics and split the surface
			! flux into disjoint patchtype categories so hybrid routing/WTD behavior
			! can be diagnosed without flux back-calculation.
			store_patch_diagnostics = .true.
			IF (present(store_patch_diagnostics_in)) store_patch_diagnostics = store_patch_diagnostics_in
			IF (store_patch_diagnostics .and. allocated(methane_finundated) .and. &
			    ipatch >= 1 .and. ipatch <= size(methane_finundated)) THEN
				methane_finundated(ipatch) = finundated
				methane_soil_finundated(ipatch) = 0._r8
				methane_soil_zwt(ipatch) = spval
				methane_surf_flux_wetland(ipatch) = 0._r8
				methane_surf_flux_soil(ipatch) = 0._r8
				methane_surf_flux_lake(ipatch) = 0._r8
				methane_surf_flux_rice(ipatch) = 0._r8
				methane_surf_aere_soil(ipatch) = 0._r8
				methane_surf_aere_rice(ipatch) = 0._r8
				methane_surf_ebul_soil(ipatch) = 0._r8
				methane_surf_ebul_rice(ipatch) = 0._r8
				methane_surf_diff_soil(ipatch) = 0._r8
				methane_surf_diff_rice(ipatch) = 0._r8
				methane_prod_tot_soil(ipatch) = 0._r8
				methane_prod_tot_rice(ipatch) = 0._r8
				methane_oxid_tot_soil(ipatch) = 0._r8
				methane_oxid_tot_rice(ipatch) = 0._r8

         ! wetland process budget components: reset every step (a patch is
         ! wetland/soil/rice/lake exclusive).  soil/rice/lake process outputs
         ! come from the existing per-category vars; only the areas span all four categories.
         methane_prod_tot_wetland(ipatch) = 0._r8
         methane_oxid_tot_wetland(ipatch) = 0._r8
         methane_surf_aere_wetland(ipatch) = 0._r8
         methane_surf_ebul_wetland(ipatch) = 0._r8
         methane_surf_diff_wetland(ipatch) = 0._r8
         methane_area_wetland(ipatch) = 0._r8 ; methane_area_soil(ipatch) = 0._r8
         methane_area_rice   (ipatch) = 0._r8
         methane_area_lake   (ipatch) = 0._r8
         ! methane_area_floodplain and methane_wetland_type are published by
         ! methane_driver before this routine runs -- do NOT reset them here.

				IF (patchtype == 2) THEN
					methane_surf_flux_wetland(ipatch) = methane_surf_flux_tot
            ! Wetland is its own patch, so the split is exact.
            methane_prod_tot_wetland(ipatch) = methane_prod_tot
            methane_oxid_tot_wetland(ipatch) = methane_oxid_tot
            methane_surf_aere_wetland(ipatch) = methane_surf_aere
            methane_surf_ebul_wetland(ipatch) = methane_surf_ebul
            methane_surf_diff_wetland(ipatch) = methane_surf_diff
            methane_area_wetland(ipatch) = 1._r8
						ELSEIF (patchtype == 0) THEN
							methane_soil_finundated(ipatch) = finundated_default
							methane_soil_zwt(ipatch) = zwt

            ! soil/rice process outputs come from the existing
            ! methane_{surf_aere,surf_ebul,surf_diff,prod_tot,oxid_tot}_{soil,rice}
            ! (non-rice branch below; Driver aggregate_methane_columns for paddy).
            ! Here we only set the non-rice soil area fraction (=1 by default);
            ! Driver overrides it with (1-rice) and sets methane_area_rice for paddy.
            methane_area_soil(ipatch) = 1._r8

							methane_surf_flux_soil(ipatch) = methane_surf_flux_tot
							methane_surf_aere_soil(ipatch) = methane_surf_aere
							methane_surf_ebul_soil(ipatch) = methane_surf_ebul
							methane_surf_diff_soil(ipatch) = methane_surf_diff
							methane_prod_tot_soil(ipatch) = methane_prod_tot
							methane_oxid_tot_soil(ipatch) = methane_oxid_tot
				ELSEIF (patchtype == 4 .and. DEF_METHANE%allowlakeprod) THEN
					methane_surf_flux_lake(ipatch) = methane_surf_flux_tot_lake
            ! lake process outputs come from base methane_{prod_tot,oxid_tot,
            ! surf_ebul,surf_diff}_lake (set on the lake path above; lakes have
            ! no aerenchyma).  Here we only mark the lake area fraction.
            methane_area_lake(ipatch) = 1._r8

				ENDIF
			ENDIF

		! Column level balance
		if (.not. methane_cold_start) then
			! Check balance
			err_methane = totcol_methane - totcolch4_bef - deltim*(methane_prod_tot - methane_oxid_tot - methane_surf_flux_tot)
			! [mol/m2]    = [mol/m2] - [mol/m2] + [s]*[mol/m2/s]
			! Keep a rate-limited diagnostic for development runs.  Production
			! configurations can additionally fail fast through the common
			! single-timestep correction threshold.
			if (DEF_METHANE%numerical_correction_fatal_threshold > 0._r8 .and. &
			    abs(err_methane) > DEF_METHANE%numerical_correction_fatal_threshold) then
				write(6,*) 'ERROR: CH4 column-budget correction exceeds fatal threshold.'
				write(6,*) 'Lat,Lon,Patchtype, err, threshold = ', dlat, dlon, patchtype, &
					err_methane, DEF_METHANE%numerical_correction_fatal_threshold
				CALL CoLM_stop ('CH4 column-budget correction exceeds fatal threshold')
			endif
			if (abs(err_methane) > 1.e-7_r8 .and. .not. warned_ch4_budget_methane) then
				write(6,*)'WARNING: CH4 column budget closure > 1e-7 mol/m^2/timestep.'
				write(6,*)'Lat,Lon,Patchtype        = ', dlat,dlon, patchtype
				write(6,*)'err_methane              = ', err_methane
				write(6,*)'totcol_methane                = ', totcol_methane
				write(6,*)'totcolch4_bef            = ', totcolch4_bef
				write(6,*)'deltim*methane_prod_tot      = ', deltim*methane_prod_tot
				write(6,*)'deltim*methane_oxid_tot      = ', deltim*methane_oxid_tot
				write(6,*)'deltim*methane_surf_flux_tot = ', deltim*methane_surf_flux_tot
				write(6,*)'fsat_bef, finundated, dfsat = ', fsat_bef, finundated, dfsat
				write(6,*)'deltim*methane_dfsat_tot = ', deltim*methane_dfsat_tot
				write(6,*)'deltim*methane_balance_residual = ', deltim*methane_balance_residual
				write(6,*)'deltim*methane_ch4_clip_credit = ', deltim*methane_ch4_clip_credit
				write(6,*)'deltim*methane_surf_diff = ', deltim*methane_surf_diff
				write(6,*)'deltim*methane_surf_aere = ', deltim*methane_surf_aere
				write(6,*)'deltim*methane_surf_ebul = ', deltim*methane_surf_ebul
				write(6,*)'deltim*methane_surf_tran = ', &
					deltim*sum(methane_tran_depth(1:nl_soil)*dz_soisno(1:nl_soil))
				warned_ch4_budget_methane = .true.
			end if
		end if

		if (present(finundated_used_out)) finundated_used_out = finundated
		if (present(finundated_default_out)) finundated_default_out = finundated_default
		fsat_bef = finundated
	end subroutine methane

	pure real(r8) function methane_distribute_grid_finundation(grid_fraction, &
		wetland_fraction_in, patchtype, debit_wetland) result(patch_fraction)
		! Grid floodwater occupies mapped wetland first; soil receives only
		! the residual.  This keeps W*f_wet + (1-W)*f_soil == f_grid.
		! With debit_wetland = .false. the wetland is not charged for the flood
		! and the soil tile takes the whole grid fraction, converted to a tile
		! basis: the CoLM mode holds the mapped wetland saturated from regional
		! groundwater, a source independent of channel overbank flow.
		real(r8), intent(in) :: grid_fraction, wetland_fraction_in
		integer, intent(in) :: patchtype
		logical, intent(in), optional :: debit_wetland
		real(r8) :: flood, wetland_fraction
		logical :: debit

		debit = .true.
		if (present(debit_wetland)) debit = debit_wetland

		if (ieee_is_nan(grid_fraction) .or. abs(grid_fraction) >= 0.5_r8*abs(spval)) then
			flood = 0._r8
		else
			flood = min(max(grid_fraction, 0._r8), 1._r8)
		endif
		if (ieee_is_nan(wetland_fraction_in) .or. &
		    abs(wetland_fraction_in) >= 0.5_r8*abs(spval)) then
			wetland_fraction = 0._r8
		else
			wetland_fraction = min(max(wetland_fraction_in, 0._r8), 1._r8)
		endif

		select case (patchtype)
		case (2)
			if (wetland_fraction > 0._r8) then
				patch_fraction = min(flood / wetland_fraction, 1._r8)
			else
				patch_fraction = 0._r8
			endif
		case (0)
			if (wetland_fraction < 1._r8) then
				if (debit) then
				patch_fraction = max(flood - wetland_fraction, 0._r8) / &
					(1._r8 - wetland_fraction)
			else
					patch_fraction = min(flood / (1._r8 - wetland_fraction), 1._r8)
				endif
			else
				patch_fraction = 0._r8
			endif
		case default
			patch_fraction = flood
		end select
	end function methane_distribute_grid_finundation

	!-----------------------------------------------------------------------
	subroutine methane_annualupdate(idate, finundated, deltim,  agnpp, bgnpp, somhr, &
		annavg_agnpp, annavg_bgnpp, annavg_somhr,  annavg_finrw, &
		tempavg_agnpp,tempavg_bgnpp,annsum_counter,tempavg_somhr, tempavg_finrw)
		!-----------------------------------------------------------------------
		! DESCRIPTION:
		! Annual mean fields.
		! Only the annavg is useful, the tempavg is temp.
		!-----------------------------------------------------------------------
		use MOD_Precision
		implicit none
		!-----------------------Argument----------------------------------------
		integer, intent(in) :: &
			idate(3)                   ! model calendar at the end of this time step

		real(r8), intent(in) :: &
			finundated              , &! fractional inundated area, =sat(0 or 1)
			deltim                  , &! land model time step (sec)

			agnpp                   , &! aboveground NPP (gC/m2/s)
			bgnpp                   , &! belowground NPP (gC/m2/s)
			somhr                      ! soil organic matter heterotrophic respiration (gC/m2/s)

		real(r8), intent(inout) :: &
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
		integer :: previous_year
		real(r8):: secsperyear          ! total seconds in the endpoint year
		real(r8):: secsperyear_previous ! total seconds in the preceding year
		real(r8):: end_sec_of_year      ! endpoint seconds elapsed in endpoint year
		real(r8):: dt_previous, dt_current
		! least accumulation [s] for a calendar-year turnover (30 days)
		real(r8), parameter :: min_turnover_gap = 30._r8*86400._r8
		!-----------------------------------------------------------------------
		! CoLM advances idate before calling the driver.  An endpoint early in a
		! year can therefore belong to a step that began in the preceding year.
		! CoLM also represents an exact midnight as day N, second 86400; that
		! representation stays in the endpoint year and is handled by the normal
		! branch below.
		if ( isleapyear(idate(1)) ) then
			secsperyear = 366._r8*secspday
		else
			secsperyear = 365._r8*secspday
		endif
		end_sec_of_year = max(0._r8, real(idate(2)-1, r8)*secspday + real(idate(3), r8))

		if (end_sec_of_year < deltim) then
			previous_year = idate(1) - 1
			if ( isleapyear(previous_year) ) then
				secsperyear_previous = 366._r8*secspday
			else
				secsperyear_previous = 365._r8*secspday
			endif
			dt_current = end_sec_of_year
			dt_previous = deltim - dt_current
			CALL add_annual_contribution(dt_previous, secsperyear_previous)
			if (annsum_counter >= min_turnover_gap) CALL finish_annual_accumulator()
			CALL add_annual_contribution(dt_current, secsperyear)
		else
			CALL add_annual_contribution(deltim, secsperyear)
			! A calendar-year endpoint must turn over even after a mid-year cold
			! start, when annsum_counter has not yet reached a full year.  CoLM
			! represents that exact endpoint as the last day with sec=86400.
			! The full-year counter serves runs whose spin-up window misses
			! 1 January; a counter running a day ahead of the calendar (site
			! restarts) turned the accumulator over at 31 December 0 h and again
			! at the year end, leaving one day in the annual means, so the
			! calendar endpoint needs min_turnover_gap since the last turnover.
			if ((end_sec_of_year >= secsperyear .and. annsum_counter >= min_turnover_gap) .or. &
				 annsum_counter >= secsperyear) &
				CALL finish_annual_accumulator()
		endif

		contains

			subroutine add_annual_contribution(dt, year_seconds)
				real(r8), intent(in) :: dt, year_seconds
				if (dt <= 0._r8 .or. year_seconds <= 0._r8) return
				annsum_counter = annsum_counter + dt
				tempavg_somhr  = tempavg_somhr + dt/year_seconds * somhr
				tempavg_finrw  = tempavg_finrw + dt/year_seconds * finundated * somhr
				tempavg_agnpp  = tempavg_agnpp + dt/year_seconds * agnpp
				tempavg_bgnpp  = tempavg_bgnpp + dt/year_seconds * bgnpp
			end subroutine add_annual_contribution

			subroutine finish_annual_accumulator()
				annsum_counter = 0._r8
				annavg_somhr = tempavg_somhr
				tempavg_somhr = 0._r8
					if (annavg_somhr > 0._r8) then
						annavg_finrw = tempavg_finrw / annavg_somhr
					else
						! No respiration-weighted annual inundation estimate is available.
						! Keep it missing so the SIF block does not treat a cold/missing
						! annual state as a valid zero and suppress first-year CH4 production.
						annavg_finrw = spval
					endif
				tempavg_finrw = 0._r8
				annavg_agnpp = tempavg_agnpp
				tempavg_agnpp = 0._r8
				annavg_bgnpp = tempavg_bgnpp
				tempavg_bgnpp = 0._r8
			end subroutine finish_annual_accumulator

	end subroutine methane_annualupdate

	!-----------------------------------------------------------------------
	subroutine methane_prod (ipatch,idate,patchtype,sat,jwt,finundated,finundated_lag,rr,deltim,& !input
      z_soisno,dz_soisno,zi_soisno,t_soisno,fliq,&
		lai,conc_o2,rootfr,annavg_finrw,&
		crootfr,somhr,lithr,hr_vr,o_scalar,fphr,pot_f_nit_vr,pH,layer_sat_lag,lake_soilc,&
      microbial_prod_potential_layer, anoxic_decomp, wt_share, &
		methane_prod_depth,o2_decomp_depth,co2_decomp_depth)!output
		!-----------------------------------------------------------------------
		! DESCRIPTION:
		! Production is done below the water table, based on CN heterotrophic respiration.
		! O2 is consumed by roots & by heterotrophic aerobes.
		! Production is done separately for sat & unsat, and is adjusted for temperature, seasonal inundation,
		! pH (optional), & redox lag factor.
		!-----------------------------------------------------------------------

		!-----------------------Argument----------------------------------------
		integer , intent(in) :: &
			ipatch                      , &! local patch index (for biome_f_methane_patch lookup)
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			patchtype                   , &! land patch type (0=soil, 1=urban or built-up,
			! 2=wetland, 3=land ice, 4=land water bodies, 99=ocean)
			sat                         , &! 0 = unsaturated; 1 = saturated
			jwt                            ! index of the soil layer right above the water table (-)

		real(r8), intent(in) :: &
			finundated                  , &! fractional inundated area in soil column
			finundated_lag              , &
			rr                          , &! root respiration (fine root MR + total root GR) (gC/m2/s)

			deltim                      , &! land model time step (sec)
			z_soisno (maxsnl+1:nl_soil) , &! layer depth (m)
			dz_soisno(maxsnl+1:nl_soil) , &! layer thickness (m)
			zi_soisno(maxsnl:nl_soil)   , &! interface level below a "z" level (m)

			t_soisno (maxsnl+1:nl_soil) , &! soil temperature (K)
         fliq     (1:nl_soil)        , &! unfrozen water fraction wliq/(wliq+wice) (-)

			lai                         , &! leaf area index [m2/m2]
			conc_o2  (1:nl_soil)        , &! O2 conc in each soil layer (mol/m3) (nl_soil)
			rootfr   (1:nl_soil)        , &! fraction of roots in each soil layer (the sum of all layer is 1)

			annavg_finrw                , &! respiration-weighted annual average of finundated

			crootfr  (1:nl_soil)        , &! fraction of roots for carbon in each soil layer (the sum of all layer is 1)

			somhr                       , &! soil organic matter heterotrophic respiration (gC/m2/s)
			lithr                       , &! litter heterotrophic respiration (gC/m2/s)
			hr_vr    (1:nl_soil)        , &! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
			o_scalar (1:nl_soil)        , &! fraction by which decomposition is limited by DEF_METHANE%anoxia
			fphr     (1:nl_soil)        , &! fraction of potential heterotrophic respiration

			pot_f_nit_vr  (1:nl_soil)   , &! potential soil nitrification flux (gN/m3/s)

			pH                          , &! soil water pH
			layer_sat_lag   (1:nl_soil) , &! Lagged saturation status of soil layer in the unsaturated zone (1 = sat)
			lake_soilc      (1:nl_soil) , &! lake sediment organic carbon (gC/m3)
			microbial_prod_potential_layer(1:nl_soil) ! optional microbial CH4 production potential (mol/m3/s)

      logical, intent(in) :: anoxic_decomp  ! saturated subcolumn decomposes at the anoxic rate
      real(r8), intent(in) :: wt_share      ! share of layer jwt+1 below the water table (-); 1 unless wt_layer_sat_share

		real(r8), intent(out) :: &
			methane_prod_depth (1:nl_soil)  , &! production of CH4 in each soil layer (nl_soil) (mol/m3/s)
			o2_decomp_depth(1:nl_soil)     , &! O2 consumption during decomposition in each soil layer (nl_soil) (mol/m3/s)
			co2_decomp_depth(1:nl_soil)       ! CO2 production from decomposition/root respiration (mol/m3/s)

		!-----------------------Local Variables---------------------------------
		integer  :: j,s              ! indices
		real(r8) :: base_decomp      ! base heterotrophic respiration rate [mol C/m2/s]
		real(r8) :: partition_z
		real(r8) :: pH_fact_methane      ! pH factor in methane production
		real(r8), parameter :: pH_fact_reference = &
			10._r8**(-0.2235_r8*6.2_r8*6.2_r8 + 2.7727_r8*6.2_r8 - 8.6_r8)

		! Factors for methanogen temperature dependence being greater than soil aerobes
		real(r8) :: f_methane_adj        ! Adjusted DEF_METHANE%f_methane
		real(r8) :: t_fact_methane       ! Temperature factor calculated using additional Q10
      real(r8) :: t_base_methane       ! reference temperature of that factor [K]
      real(r8) :: f_methane_tile       ! f_methane of this tile, tropical peat share weighted

		! O2 limitation on decomposition and methanogenesis
		real(r8) :: seasonalfin      ! finundated in excess of respiration-weighted annual average

		! For calculating column average (rootfrac(j)*rr(j))
		real(r8) :: rr_vr(1:nl_soil) ! vertically resolved column-mean root respiration (g C/m2/s)

		real(r8) :: sif              ! (unitless) ratio applied to sat. prod. to account for seasonal inundation
			real(r8) :: q10lake_eff      ! lake Q10 for sediment carbon decomposition
			real(r8) :: lake_ch4_fraction
      real(r8) :: lake_c_col               ! column lake sediment carbon [gC m-2]
      real(r8) :: lake_w(1:nl_soil)        ! depth weights of lake decomposition [-]
      logical  :: use_lake_pool            ! the active pool feeds this lake patch
			real(r8) :: carbon_decomp_depth
			real(r8) :: legacy_methane_prod_depth
			real(r8) :: microbe_prod_cap
			real(r8) :: freeze_factor
			logical  :: use_microbe_override
		real(r8) :: walter_renorm        ! Walter 2001 depth attenuation renormalization factor
		real(r8) :: walter_w_raw, walter_w_atten   ! pre-pass HR weights for renorm
      ! methanogen activity (JULES-microbe, Chadburn et al. 2020)
      logical  :: use_act              ! activity applies to this call (wetland tile, switch on)
      real(r8) :: act(1:nl_soil)       ! activity m_act of this subcolumn's layers [-]
      real(r8) :: act_accl             ! acclimation factor a of k2 (eq. A6) [-]
      real(r8) :: act_min              ! activity floor alpha/f [-]
      real(r8) :: act_rate             ! activity change rate f k2 a A2 [s-1]
      real(r8) :: prod_full            ! production at full activity [mol/m3/s]
		!-----------------------------------------------------------------------
		! PATCH loop to calculate vertically resolved column-averaged root respiration
		do j=1,nl_soil
			rr_vr(j) = rr*crootfr(j)
			! [g C/m2/s]=[g C/m2/s]*[-]
		end do

		partition_z = 1._r8
		base_decomp = 0.0_r8
		use_microbe_override = DEF_METHANE%use_microbial_pools .and. &
			DEF_METHANE%use_microbial_flux_override .and. patchtype /= 4

		! ---- Walter & Heimann 2001 depth-attenuation pre-pass --------------
		! Compute a factor that preserves the unmodified HR partition-weight sum
		! while applying exp(-z/z0_methane_prod) to partition_z later.  This does
		! not claim to preserve final CH4 production after layer-specific
		! temperature, pH, and redox modifiers.  Use the same potential-
		! production domain in this pre-pass and in the application below.
		walter_renorm = 1._r8
		if (DEF_METHANE%z0_methane_prod > 0._r8 .and. patchtype /= 4) then
			walter_w_raw   = 0._r8
			walter_w_atten = 0._r8
			do j = 1, nl_soil
				if (j > jwt .or. DEF_METHANE%anoxicmicrosites) then
					walter_w_raw   = walter_w_raw   + hr_vr(j) * dz_soisno(j)
					walter_w_atten = walter_w_atten + hr_vr(j) * dz_soisno(j) * &
					                 exp(-z_soisno(j) / DEF_METHANE%z0_methane_prod)
				endif
			end do
			if (walter_w_atten > 0._r8) walter_renorm = walter_w_raw / walter_w_atten
		endif

      ! a wetland tile takes its own annual mean soil temperature as the
      ! reference of the production temperature factor once one is complete
      t_base_methane = DEF_METHANE%q10methane_base
      if (DEF_METHANE%q10methane_local_base .and. patchtype == 2 .and. allocated(annavg_tsoi)) then
         if (ipatch >= 1 .and. ipatch <= size(annavg_tsoi)) then
            if (annavg_tsoi(ipatch) > 150._r8 .and. annavg_tsoi(ipatch) < 350._r8) &
               t_base_methane = DEF_METHANE%q10methane_base + DEF_METHANE%q10methane_base_weight &
                  * (annavg_tsoi(ipatch) - DEF_METHANE%q10methane_base)
         endif
      endif

      ! the tropical peat share of a wetland tile takes its own f_methane
      f_methane_tile = DEF_METHANE%f_methane
      if (DEF_METHANE%wetland_f_tropical_peat >= 0._r8 .and. patchtype == 2 .and. &
          allocated(wetland_tile_peat_share)) then
         if (ipatch >= 1 .and. ipatch <= size(wetland_tile_peat_share)) then
            if (abs(patchlatr(ipatch)) <= 23.5_r8 * deg2rad) &
               f_methane_tile = f_methane_tile + wetland_tile_peat_share(ipatch) * &
                  (DEF_METHANE%wetland_f_tropical_peat - f_methane_tile)
         endif
      endif

      ! the forested share of a wetland tile at wetland_f_tree_ratio of f
      if (DEF_METHANE%wetland_f_tree_ratio >= 0._r8 .and. patchtype == 2 .and. &
          allocated(wetveg_forest)) then
         if (ipatch >= 1 .and. ipatch <= size(wetveg_forest)) &
            f_methane_tile = f_methane_tile * &
               (1._r8 - wetveg_forest(ipatch) * (1._r8 - DEF_METHANE%wetland_f_tree_ratio))
      endif

      ! activity of the methanogens in this subcolumn of a wetland
      ! tile; k2 acclimates to the tile's annual mean soil temperature (eq. A6)
      use_act = DEF_METHANE%methanogen_activity .and. patchtype == 2 .and. &
         allocated(methanogen_act_sat) .and. allocated(methanogen_act_unsat)
      if (use_act) use_act = ipatch >= 1 .and. ipatch <= size(methanogen_act_sat, 2)
      act_accl = 1._r8
      act_min  = 0._r8
      if (use_act) then
         if (sat == 1) then
            act(:) = methanogen_act_sat(1:nl_soil,ipatch)
         else
            act(:) = methanogen_act_unsat(1:nl_soil,ipatch)
         endif
         act_min = DEF_METHANE%methanogen_alpha / DEF_METHANE%methanogen_cue
         if (DEF_METHANE%methanogen_qev_ratio > 0._r8 .and. allocated(annavg_tsoi)) then
            if (ipatch <= size(annavg_tsoi)) then
               if (annavg_tsoi(ipatch) > 150._r8 .and. annavg_tsoi(ipatch) < 350._r8) &
                  act_accl = 1._r8 / methanogen_arrhenius(annavg_tsoi(ipatch), &
                     DEF_METHANE%methanogen_qev_ratio * DEF_METHANE%methanogen_q2)
            endif
         endif
      endif

		! column loop to partition decomposition_rate into each soil layer
		q10lake_eff = DEF_METHANE%q10methane * 1.5_r8
      if (DEF_METHANE%q10lake > 0._r8) q10lake_eff = DEF_METHANE%q10lake
      ! the reactive sediment is the fresh deposit near the lake bed;
      ! organic matter reactivity falls with age and so with depth (Middelburg
      ! 1989), so the column's decomposition is spread as exp(-z/efold) instead
      ! of following the carbon profile, with the column total unchanged
      lake_c_col = 0._r8
      lake_w(:) = 0._r8
      if (patchtype == 4 .and. DEF_METHANE%allowlakeprod .and. DEF_METHANE%lake_prod_efold > 0._r8) then
         lake_c_col = sum(max(lake_soilc(1:nl_soil), 0._r8) * dz_soisno(1:nl_soil))
         lake_w(1:nl_soil) = exp(-z_soisno(1:nl_soil) / DEF_METHANE%lake_prod_efold) * dz_soisno(1:nl_soil)
         lake_w(1:nl_soil) = lake_w(1:nl_soil) / max(sum(lake_w(1:nl_soil)), tiny(1._r8))
      endif
      ! the active pool (Const, lake_active_pool) takes the
      ! weights when lake_prod_efold is set, else the lake_soilc carbon
      ! profile the column path follows; the pool itself is updated in methane
      use_lake_pool = patchtype == 4 .and. DEF_METHANE%allowlakeprod .and. &
         DEF_METHANE%lake_active_pool .and. allocated(lake_cact)
      if (use_lake_pool) use_lake_pool = ipatch >= 1 .and. ipatch <= size(lake_cact)
      if (use_lake_pool .and. DEF_METHANE%lake_prod_efold <= 0._r8) then
         lake_w(1:nl_soil) = max(lake_soilc(1:nl_soil), 0._r8) * dz_soisno(1:nl_soil)
         if (sum(lake_w(1:nl_soil)) > 0._r8) then
            lake_w(1:nl_soil) = lake_w(1:nl_soil) / sum(lake_w(1:nl_soil))
         else
            lake_w(:) = 0._r8
            lake_w(1) = 1._r8
         endif
      endif
		do j=1,nl_soil
			if (patchtype == 4 .and. .not. DEF_METHANE%allowlakeprod) then
				methane_prod_depth(j) = 0._r8
				o2_decomp_depth(j) = 0._r8
				co2_decomp_depth(j) = 0._r8
				cycle
			endif
			if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
            if (use_lake_pool) then
               ! layer j mineralizes its share of the active
               ! pool at k(T) = lake_k20 lake_theta**(T - 20 C) (FLaMe eq. 7)
               base_decomp = DEF_METHANE%lake_k20 / secspday * &
                  DEF_METHANE%lake_theta**(t_soisno(j) - tfrz - 20._r8) * &
                  lake_cact(ipatch) * lake_w(j) / catomw
               ! [mol C m-2 s-1] = [d-1]/[s d-1] * [gC m-2] / [gC mol-1]
            elseif (DEF_METHANE%lake_prod_efold > 0._r8) then
               base_decomp = DEF_METHANE%lake_decomp_fact * DEF_METHANE%cnscalefactor * &
                  lake_c_col * lake_w(j) * &
                  q10lake_eff**((t_soisno(j) - DEF_METHANE%q10lakebase) / 10._r8) / catomw
            else
				base_decomp = DEF_METHANE%lake_decomp_fact * DEF_METHANE%cnscalefactor * &
					max(lake_soilc(j), 0._r8) * dz_soisno(j) * &
					q10lake_eff**((t_soisno(j) - DEF_METHANE%q10lakebase) / 10._r8) / catomw
            endif
				freeze_factor = min(1._r8, max(0._r8, t_soisno(j) - tfrz + 1._r8))
				base_decomp = base_decomp * freeze_factor
				if (j > jwt) then
					! Keep the anaerobic carbon split stoichiometrically bounded:
					! 2CH2O -> CH4 + CO2 implies f_ch4 <= 0.5 of decomposed C.
					lake_ch4_fraction = min(max(DEF_METHANE%f_methane, 0._r8), 0.5_r8)
               ! the active pool takes its own share, lake_f_methane
               if (use_lake_pool .and. DEF_METHANE%lake_f_methane >= 0._r8) &
                  lake_ch4_fraction = min(DEF_METHANE%lake_f_methane, 0.5_r8)
				else
					lake_ch4_fraction = 0._r8
				endif
				methane_prod_depth(j) = lake_ch4_fraction * base_decomp / dz_soisno(j)
				! Lake sediment production is treated as an anaerobic carbon split.
				! Diagnose the non-CH4 carbon as CO2 below, but do not create an
				! artificial O2 demand in the fully saturated lake sediment column.
				o2_decomp_depth(j) = 0._r8
				co2_decomp_depth(j) = max(0._r8, base_decomp / dz_soisno(j) - methane_prod_depth(j))
            ! where O2 reaches the sediment, respiration of the non-CH4
            ! carbon uses it (1 mol O2 per mol CO2); the O2 competition caps it
            ! at the O2 present, so anoxic layers are unaffected and the O2 of a
            ! cold start no longer lingers to oxidise the CH4 made below.
            if (DEF_METHANE%lake_sed_o2_demand) o2_decomp_depth(j) = co2_decomp_depth(j)
				cycle
			endif
			! Use soil heterotrophic respiration (based on Wania)
			base_decomp = (somhr+lithr) / catomw
			! [mol C/m2/s]=[g C/m2/s]   / [g C/mol C]

			! Multiply base_decomp by factor accounting for lower carbon stock in seasonally inundated areas than
			! if it were inundated all year.
			! This is to reduce emissions in seasonally inundated zones, because the eq.
			! C-flux will be less than predicted by a non-O2-lim model
			!
			! Apply SIF only when BGC decomposition is not already anoxia-limited.
			if (sat == 1 .and. DEF_METHANE%use_ch4_sif .and. &
			    .not. DEF_METHANE%bgc_anoxia_limits_decomp) then
				sif = 1._r8
				if (annavg_finrw /= spval .and. finundated > 0._r8) then
					seasonalfin = max(finundated-annavg_finrw, 0._r8)
					if (seasonalfin > 0._r8) then
						sif = (annavg_finrw + DEF_METHANE%mino2lim*seasonalfin) / finundated
						base_decomp = base_decomp * sif
					end if
				end if
			end if

         ! the soil tile's BGC carries no anoxia limit (o_scalar
         ! is 1) and flooding never builds up its carbon stock, yet CH4 is made
         ! only below a water table, in anoxic soil. In both subcolumns that soil
         ! decomposes at the anoxic rate, mino2lim times the aerobic one: the
         ! seasonal inundation factor above with no permanently inundated share
         ! (Riley et al. 2011, appendix B).
         ! (anoxic_decomp_flooded_only): the flooded subcolumn only; the
         ! non-flooded one keeps the aerobic rate (appendix B: a layer anoxic
         ! most of the year builds its carbon back up to the input rate)
         if (anoxic_decomp .and. (sat == 1 .or. .not. DEF_METHANE%anoxic_decomp_flooded_only)) then
            if (DEF_METHANE%floodplain_sif .and. sat == 1) then
               ! the flooded subcolumn of a floodplain takes the CLM4Me
               ! seasonal inundation factor (Riley et al. 2011, appendix B):
               ! the share inundated on the annual (respiration-weighted)
               ! average has carbon built up under anoxia and respires at the
               ! full rate at steady state; only the seasonal excess is slowed
               ! to mino2lim. mino2lim alone when no annual value exists.
               sif = DEF_METHANE%mino2lim
               if (annavg_finrw /= spval .and. finundated > 0._r8) then
                  seasonalfin = max(finundated - annavg_finrw, 0._r8)
                  sif = (min(annavg_finrw, finundated) + DEF_METHANE%mino2lim * seasonalfin) / finundated
               end if
               base_decomp = base_decomp * sif
            else
               base_decomp = base_decomp * DEF_METHANE%mino2lim
            end if
         end if

			! For sensitivity studies
			base_decomp = base_decomp * DEF_METHANE%cnscalefactor


			! depth dependence of production either from rootfr or decomp model
			if ( (somhr + lithr) > 0._r8) then
				partition_z = hr_vr(j) * dz_soisno(j) / (somhr + lithr)
				!   [-]     = [g C/m3/s]*[m]/[g C/m2/s]
			else
				partition_z = 1._r8
			end if

			! Walter & Heimann 2001 explicit methanogenesis depth attenuation.
			! Multiply by exp(-z/z0)*walter_renorm to preserve the pre-modifier
			! HR partition sum while sharpening the source profile (top layers up,
			! deep layers down).  Subsequent layer modifiers may change the total.
			if (DEF_METHANE%z0_methane_prod > 0._r8 .and. patchtype /= 4 .and. &
			    (j > jwt .or. DEF_METHANE%anoxicmicrosites)) then
				partition_z = partition_z * exp(-z_soisno(j) / DEF_METHANE%z0_methane_prod) * walter_renorm
			endif


			! Adjust DEF_METHANE%f_methane to account for the fact that methanogens may have a higher Q10 than aerobic decomposers.
			! Note this is crude and should ideally be applied to all anaerobic decomposition rather than just the
			! DEF_METHANE%f_methane.
			f_methane_adj = 1.0_r8

         t_fact_methane = DEF_METHANE%q10methane**((t_soisno(j) - t_base_methane)/10._r8)


			! Adjust f_methane by the temperature ratio.  When biome lookup
			! is enabled, use the per-patch value populated by methane_driver
			! via get_biome_f_methane; otherwise fall back to global scalar.
			if (DEF_METHANE%use_biome_f_methane .and. allocated(biome_f_methane_patch) .and. &
			    ipatch >= 1 .and. ipatch <= size(biome_f_methane_patch)) then
				f_methane_adj = biome_f_methane_patch(ipatch) * t_fact_methane
			else
            f_methane_adj = f_methane_tile * t_fact_methane
			endif

			! Do not allow soil methanogenesis in frozen layers.  Lake sediment
			! production already has an explicit temperature gate; apply the same
			! physical cutoff to the soil branch so winter CH4 is not produced
			! solely from residual BGC HR in ice-filled pores.
         ! With liquid_fraction_scaling the cutoff becomes the unfrozen water
         ! fraction: methanogens keep working in the liquid films of a freezing
         ! layer and stop only as the pore water actually freezes.
         ! with prod_fliq_off the decomposition moisture scalar
         ! alone limits production in a freezing layer, as in CTSM with
         ! vertically resolved soil carbon (ch4Mod: base_decomp is zeroed only
         ! for nlevdecomp == 1 or lakes)
         if (DEF_METHANE%liquid_fraction_scaling) then
            if (.not. (DEF_METHANE%prod_fliq_off .and. patchtype /= 4)) &
            f_methane_adj = f_methane_adj * fliq(j)
         else
			if (t_soisno(j) <= tfrz) f_methane_adj = 0._r8
         end if


			! Remove CN nitrogen limitation, as methanogenesis is not N limited.
			! Also remove (low) moisture limitation
			if (DEF_METHANE%methane_rmcnlim) then
				if (fphr(j) > 0._r8) then
					f_methane_adj = f_methane_adj / fphr(j)
				end if
			end if


				! If switched on, use pH factor for production based on spatial pH data defined in surface data.
         ! the flooded subcolumn of a soil tile (floodplain or
         ! paddy) takes no pH limitation: flooded mineral soils converge to
         ! near-neutral pH as reduction consumes protons
         if (patchtype /= 4 .and. DEF_METHANE%usephfact .and. &
             .not. (DEF_METHANE%flooded_soil_ph_neutral .and. patchtype == 0 .and. sat == 1))then
					if (pH <= DEF_METHANE%pHmin .or. pH >= DEF_METHANE%pHmax) then
						! Outside the empirical Dunfield production window, suppress
						! methanogenesis.  The old branch skipped the multiplier outside
						! [pHmin,pHmax], accidentally treating extreme pH as pH-neutral.
               ! With ph_factor_floor the floor is kept there as well.
               if (DEF_METHANE%ph_factor_floor > 0._r8) then
                  f_methane_adj = f_methane_adj * DEF_METHANE%ph_factor_floor
               else
						f_methane_adj = 0._r8
               end if
					else
						pH_fact_methane = 10._r8**(-0.2235_r8*pH*pH + 2.7727_r8*pH - 8.6_r8) / pH_fact_reference
						pH_fact_methane = min(1._r8, max(0._r8, pH_fact_methane))
						! fitted function using data from Dunfield et al. 1993
						! Normalized to one at pH=6.2.
						! From Lei Meng
               ! floor for communities at their native pH (ph_factor_floor)
               if (DEF_METHANE%ph_factor_floor > 0._r8) &
                  pH_fact_methane = DEF_METHANE%ph_factor_floor &
                                  + (1._r8 - DEF_METHANE%ph_factor_floor) * pH_fact_methane
						f_methane_adj = f_methane_adj * pH_fact_methane
					end if
				else
					! if no data, then no pH effects
				end if

         ! sulfate reducers outcompete methanogens in brackish and saline
         ! wetlands; tidal-marsh emission falls with porewater salinity as
         ! log10(CH4) = -0.056 S + 1.38 (Poffenbarger et al. 2011), so the
         ! yield relative to fresh water is 10**(-0.056 S)
         if (patchtype == 2 .and. DEF_METHANE%wetland_salinity_site > 0._r8) then
            f_methane_adj = f_methane_adj * 10._r8**(-0.056_r8 * DEF_METHANE%wetland_salinity_site)
         end if

         ! nitrate reducers outcompete methanogens for the same electron
         ! donors; PEPRMT-Tidal (Oikawa et al. 2024, eq. 5) scales tidal-marsh
         ! production by K / (K + [NO3])
         if (patchtype == 2 .and. DEF_METHANE%wetland_nitrate_site > 0._r8) then
            f_methane_adj = f_methane_adj * DEF_METHANE%wetland_nitrate_ki &
                          / (DEF_METHANE%wetland_nitrate_ki + DEF_METHANE%wetland_nitrate_site)
         end if

			! Redox factor
			if ((patchtype /= 4) .and. sat==1 .and. finundated_lag < finundated)  then
				f_methane_adj = f_methane_adj * finundated_lag / finundated
			elseif (sat == 0 .and. j > jwt) then ! Assume lag in decay of alternative electron acceptors vertically
				f_methane_adj = f_methane_adj * layer_sat_lag(j)
			end if
			! Alternative electron acceptors will be consumed first after soil is inundated.

			f_methane_adj = min(f_methane_adj, 0.5_r8)
			! Must be less than 0.5 because otherwise the actual implied aerobic respiration would be negative.
			! The total of aer. respiration + methanogenesis must remain equal to the SOMHR calculated in CN,
			! so that the NEE is sensible. Even perfectly anaerobic conditions with no alternative
			! electron acceptors would predict no more than 0.5 b/c some oxygen is present in organic matter.
			! e.g. 2CH2O --> CH4 + CO2.

         ! (wt_layer_sat_share): the layer holding the water table makes
         ! CH4 on its saturated share; the share above the table is aerobic
         ! soil, producing only through anoxic microsites when those are on.
         ! The carbon not turned into CH4 stays in the aerobic O2 demand and
         ! the CO2 below.
         if (j == jwt + 1 .and. wt_share < 1._r8) then
            if (DEF_METHANE%anoxicmicrosites) then
               f_methane_adj = f_methane_adj * (wt_share + (1._r8 - wt_share) &
                  / (1._r8 + DEF_METHANE%oxinhib * conc_o2(j)))
            else
               f_methane_adj = f_methane_adj * wt_share
            endif
         endif

			carbon_decomp_depth = max(0._r8, base_decomp * partition_z / dz_soisno(j))
			legacy_methane_prod_depth = 0._r8
			if (j > jwt) then
				legacy_methane_prod_depth = &
					max(0._r8, f_methane_adj * base_decomp * partition_z / dz_soisno(j))
			elseif (DEF_METHANE%anoxicmicrosites) then
				legacy_methane_prod_depth = max(0._r8, &
					f_methane_adj * base_decomp * partition_z / dz_soisno(j) / &
					(1._r8 + DEF_METHANE%oxinhib * conc_o2(j)))
			endif
			microbe_prod_cap = 0.5_r8 * carbon_decomp_depth
			if (legacy_methane_prod_depth > 0._r8) then
				microbe_prod_cap = min(microbe_prod_cap, &
					DEF_METHANE%max_microbe_prod_multiplier * legacy_methane_prod_depth)
			elseif (use_microbe_override) then
				! If the legacy pathway forbids production (e.g. frozen layers or
				! above-WT layers without anoxic microsites), keep the microbial
				! override at zero even if a stale/future potential is nonzero.
				microbe_prod_cap = 0._r8
			endif
			if (use_microbe_override .and. &
			    abs(microbial_prod_potential_layer(j)) < 0.5_r8 * abs(spval)) then
				if (j > jwt .or. DEF_METHANE%anoxicmicrosites) then
					methane_prod_depth(j) = min(max(0._r8, microbial_prod_potential_layer(j)), &
						microbe_prod_cap)
				else
					methane_prod_depth(j) = 0._r8
				endif
			elseif (j  >  jwt) then ! Below the water table so anaerobic CH4 production can occur
				! partition decomposition to layer
				! turn into per volume-total by dz
				methane_prod_depth(j) = f_methane_adj * base_decomp * partition_z / dz_soisno(j)
				! [mol C/m3/s]    = [-]       * [mol C/m2/s]* [-]         / [m]
			else ! Above the WT
				if (DEF_METHANE%anoxicmicrosites) then
					methane_prod_depth(j) = f_methane_adj * base_decomp * partition_z / dz_soisno(j) &
						/ (1._r8 + DEF_METHANE%oxinhib*conc_o2(j))
				else
					methane_prod_depth(j) = 0._r8 ! [mol/m3 total/s]
				endif ! DEF_METHANE%anoxicmicrosites
			endif ! WT

         ! biomass and substrate held at their full-activity
         ! steady state, production is the full-activity production times the
         ! activity at the start of the step (Chadburn et al. 2020, eq. A7).
         ! The activity then rises as f k2 a A2 m_act while the specific growth
         ! rate f k2 a A2 (phi = 1) exceeds mu in a layer below the water table
         ! with substrate supply, and falls at that rate otherwise (eq. A5),
         ! within [alpha/f, 1].
         prod_full = methane_prod_depth(j)
         if (use_act) then
            methane_prod_depth(j) = prod_full * act(j)
            act_rate = DEF_METHANE%methanogen_cue * DEF_METHANE%methanogen_k2 / 3600._r8 * &
               act_accl * methanogen_arrhenius(t_soisno(j), DEF_METHANE%methanogen_q2)
            if (j > jwt .and. prod_full > 0._r8 .and. &
                act_rate * 3600._r8 > DEF_METHANE%methanogen_mu) then
               act(j) = min(1._r8, act(j) * exp(act_rate * deltim))
            else
               act(j) = max(act_min, act(j) * exp(-act_rate * deltim))
            endif
         endif

			! Heterotrophic O2 demand is the aerobic remainder, not the total
			! decomposed carbon.  Each mole of CH4 represents two moles of carbon
			! decomposed anaerobically by the reaction documented above.
			! the full-activity production is used, as the carbon
			! dormant methanogens leave unconverted is not respired aerobically.
			o2_decomp_depth(j) = max(0._r8, carbon_decomp_depth - &
				2._r8 * prod_full)
			! Recover potential aerobic demand only after the anaerobic share has
			! been removed.  Dividing total decomposition first would charge the
			! methanogenic carbon an O2 demand again.
			if (DEF_METHANE%anoxia .and. .not. DEF_METHANE%bgc_anoxia_limits_decomp) then
				if (o_scalar(j) > 0._r8) &
					o2_decomp_depth(j) = o2_decomp_depth(j) / max(o_scalar(j), 0.01_r8)
			end if

			! Root respiration contributes CO2 in every rooted layer, so its O2
			! demand must use the same vertical support, including saturated soil.
			if (patchtype /= 4) then
				o2_decomp_depth(j) = o2_decomp_depth(j) + rr_vr(j)/catomw/dz_soisno(j)
				! [mol/m3/s]       = [mol/m3/s]         + [g C/m2/s]/[g C/mol C]/[m]
			end if

			! Add oxygen demand for nitrification.
			if (DEF_METHANE%use_nitrif_denitrif) then
				o2_decomp_depth(j) = o2_decomp_depth(j) + pot_f_nit_vr(j) * 2.0_r8/14.0_r8
				! [mol/m3/s]       = [mol/m3/s]         + [g N/m3/s]/[g N/mol N]
			end if

			co2_decomp_depth(j) = max(0._r8, &
				carbon_decomp_depth - methane_prod_depth(j) + &
				rr_vr(j) / catomw / dz_soisno(j))

		end do ! nl_soil

      if (use_act) then
         if (sat == 1) then
            methanogen_act_sat(1:nl_soil,ipatch) = act(:)
         else
            methanogen_act_unsat(1:nl_soil,ipatch) = act(:)
         endif
      endif

	end subroutine methane_prod

	!---------------------------------------------------------------------------
   subroutine methane_acceptor (ipatch, sat, jwt, deltim, fliq, cellorg, organic_max, conc_o2_aqu, &
      methane_prod_depth, o2_decomp_depth, co2_decomp_depth, reox)
      !-----------------------------------------------------------------------
      ! DESCRIPTION:
      ! alternative electron acceptors of one subcolumn of a
      ! non-lake patch (Segers and Leffelaar 2001, JGR 106, eqs. 18, 20-30,
      ! p.3515; parameters in MOD_Tracer_Reactive_Methane_Const). Below the
      ! water table the acceptors take the share 1 - zeta of the anaerobic
      ! carbon flow 2 P, whose carbon leaves as CO2 instead of CH4, and are
      ! reduced by 8 mol e- per mol of CH4 not made; the pool is updated here.
      ! The reoxidation the aqueous O2 would drive is returned as a rate of
      ! the oxidised share and its O2 demand added to the decomposition
      ! demand; the caller applies it after the O2 competition of
      ! methane_tran, scaled by the O2 stress that scaled that demand, so the
      ! O2 taken and the acceptors oxidised stay in step.
      !-----------------------------------------------------------------------
      integer,  intent(in) :: ipatch                    ! patch index
      integer,  intent(in) :: sat                       ! 0 = unsaturated; 1 = saturated
      integer,  intent(in) :: jwt                       ! index of the soil layer right above the water table (-)
      real(r8), intent(in) :: deltim                    ! land model time step (s)
      real(r8), intent(in) :: fliq(1:nl_soil)           ! unfrozen water fraction wliq/(wliq+wice) (-)
      real(r8), intent(in) :: cellorg(1:nl_soil)        ! organic matter density (kg OM m-3)
      real(r8), intent(in) :: organic_max               ! organic matter density where soil acts as peat (kg OM m-3)
      real(r8), intent(in) :: conc_o2_aqu(1:nl_soil)    ! aqueous O2 concentration (mol m-3 water)
      real(r8), intent(inout) :: methane_prod_depth(1:nl_soil) ! CH4 production (mol m-3 s-1)
      real(r8), intent(inout) :: o2_decomp_depth(1:nl_soil)    ! O2 demand of decomposition (mol m-3 s-1)
      real(r8), intent(inout) :: co2_decomp_depth(1:nl_soil)   ! CO2 from decomposition (mol m-3 s-1)
      real(r8), intent(out)   :: reox(1:nl_soil)               ! potential reoxidation of the oxidised share (s-1)

      integer  :: j
      real(r8) :: f_om       ! organic share of the layer (-)
      real(r8) :: cap        ! acceptor capacity A_max (mol e- m-3)
      real(r8) :: share      ! oxidised share A / A_max (-)
      real(r8) :: pool       ! oxidised acceptors A (mol e- m-3)
      real(r8) :: zeta       ! share of the anaerobic carbon flow going to methanogenesis (-)
      real(r8) :: divert     ! CH4 production taken over by the acceptors (mol CH4 m-3 s-1)
      real(r8) :: f_aer      ! O2 dependence of the reoxidation (-)
      real(r8), parameter :: e_per_ch4 = 8._r8     ! mol e- per mol CH4 not made: 2 mol C at nu_rd = 4 (Table 1, p.3512)
      real(r8), parameter :: o2_per_e  = 0.25_r8   ! nu_ro, mol O2 per mol e- reoxidised (Table 1, p.3512)

      reox(:) = 0._r8
      do j = 1, nl_soil
         f_om = 0._r8
         if (organic_max > 0._r8 .and. cellorg(j) > 0._r8 .and. cellorg(j) < 1.e29_r8) &
            f_om = min(cellorg(j) / organic_max, 1._r8)
         cap = (1._r8 - f_om) * DEF_METHANE%acceptor_cap_min + f_om * DEF_METHANE%acceptor_cap_org
         if (cap <= 0._r8) cycle
         if (sat == 1) then
            share = acceptor_sat(j,ipatch)
         else
            share = acceptor_unsat(j,ipatch)
         endif
         share = min(max(share, 0._r8), 1._r8)
         pool = share * cap

         ! Reduction below the water table: methanogenesis keeps the share
         ! zeta (eqs. 24-27); the acceptors take the rest while they last
         ! (eq. 22). Production above the table (anoxic microsites) keeps its
         ! own O2 inhibition and is left alone.
         if (j > jwt .and. methane_prod_depth(j) > 0._r8 .and. pool > 0._r8) then
            zeta = 1._r8 / (1._r8 + DEF_METHANE%acceptor_eta0 * pool / (DEF_METHANE%acceptor_k_half + pool))
            divert = min((1._r8 - zeta) * methane_prod_depth(j), pool / (e_per_ch4 * deltim))
            methane_prod_depth(j) = methane_prod_depth(j) - divert
            ! 2 CH2O -> CH4 + CO2 becomes 2 CO2: one more CO2 per CH4 not made
            co2_decomp_depth(j) = co2_decomp_depth(j) + divert
            pool = max(0._r8, pool - e_per_ch4 * divert * deltim)
            share = pool / cap
         endif

         ! Reoxidation of the reduced part where O2 and liquid water are
         ! present (eqs. 18, 28), relaxing exactly over the step; O2 at nu_ro
         ! per electron (eq. 30)
         f_aer = max(conc_o2_aqu(j), 0._r8)
         f_aer = f_aer / (f_aer + DEF_METHANE%acceptor_k_o2)
         reox(j) = (1._r8 - share) * (1._r8 - exp(-f_aer * min(max(fliq(j), 0._r8), 1._r8) * &
            DEF_METHANE%acceptor_k_reox * deltim)) / deltim
         o2_decomp_depth(j) = o2_decomp_depth(j) + o2_per_e * reox(j) * cap

         if (sat == 1) then
            acceptor_sat(j,ipatch) = share
         else
            acceptor_unsat(j,ipatch) = share
         endif
      enddo

   end subroutine methane_acceptor

   !---------------------------------------------------------------------------
   subroutine methane_root_oxid (patchtype, sat, jwt, is_rice_paddy, lai, dz_soi, rootfr, &
      methane_prod_depth, prod_tran, rox_depth)
      !-----------------------------------------------------------------------
      ! DESCRIPTION:
      ! with rox_on_rooted_production: the fixed rhizosphere
      ! oxidation share x of one subcolumn acts on the CH4 produced below the
      ! water table in the rooted layers. x is rice_aereoxid on a paddy
      ! column (> 0 is on), wetland_aereoxid on a wetland tile and on the
      ! flooded subcolumn of another soil patch (>= 0 is on), as for the
      ! plant O2 cut-off in methane. Layer j below the table keeps
      ! P_j - R_j for the transport, R_j = x w_j P_j, where
      !    w_j = rho_j / max_k rho_k,   rho_j = rootfr_j / dz_j
      ! is the layer's root density relative to the most densely rooted
      ! layer of the profile. The caller books R as CH4 oxidation after the
      ! O2 competition of methane_tran; its O2 comes from the plant, so as
      ! methane_aere it needs a live canopy (lai > 0): a harvested paddy or
      ! a leafless tile oxidises nothing at the roots.
      !-----------------------------------------------------------------------
      integer,  intent(in)  :: patchtype                     ! land patch type
      integer,  intent(in)  :: sat                           ! 0 = unsaturated; 1 = saturated
      integer,  intent(in)  :: jwt                           ! index of the soil layer right above the water table (-)
      logical,  intent(in)  :: is_rice_paddy                 ! paddy rice column
      real(r8), intent(in)  :: lai                           ! leaf area index of the column (m2 m-2)
      real(r8), intent(in)  :: dz_soi(1:nl_soil)             ! layer thickness (m)
      real(r8), intent(in)  :: rootfr(1:nl_soil)             ! fraction of roots in each soil layer (-)
      real(r8), intent(in)  :: methane_prod_depth(1:nl_soil) ! CH4 production (mol m-3 s-1)
      real(r8), intent(out) :: prod_tran(1:nl_soil)          ! production left to the pore water (mol m-3 s-1)
      real(r8), intent(out) :: rox_depth(1:nl_soil)          ! rhizosphere oxidation R (mol m-3 s-1)

      integer  :: j
      real(r8) :: x          ! rhizosphere oxidation share of the subcolumn (-)
      real(r8) :: rho_max    ! largest root density of the profile (m-1)
      real(r8) :: w          ! relative root density of the layer (-)

      prod_tran(:) = methane_prod_depth(:)
      rox_depth(:) = 0._r8

      x = -1._r8
      if (is_rice_paddy) then
         if (DEF_METHANE%rice_aereoxid > 0._r8) x = DEF_METHANE%rice_aereoxid
      else if (patchtype == 2 .or. (patchtype == 0 .and. sat == 1)) then
         if (DEF_METHANE%wetland_aereoxid >= 0._r8) x = DEF_METHANE%wetland_aereoxid
      end if
      if (x <= 0._r8 .or. .not. (lai > 0._r8)) return

      rho_max = 0._r8
      do j = 1, nl_soil
         if (dz_soi(j) > 0._r8) rho_max = max(rho_max, max(rootfr(j), 0._r8) / dz_soi(j))
      end do
      if (rho_max <= 0._r8) return

      do j = max(jwt, 0) + 1, nl_soil
         if (methane_prod_depth(j) > 0._r8 .and. dz_soi(j) > 0._r8) then
            w = min(max(rootfr(j), 0._r8) / dz_soi(j) / rho_max, 1._r8)
            rox_depth(j) = x * w * methane_prod_depth(j)
            prod_tran(j) = methane_prod_depth(j) - rox_depth(j)
         end if
      end do

   end subroutine methane_root_oxid

   !---------------------------------------------------------------------------
   subroutine methane_oxid (idate, patchtype, jwt,  sat, t_soisno, fliq, dz_soisno, zi_soisno, smp, vol_aqu, &
		conc_o2_aqu_porsl, conc_ch4_aqu_porsl, &
		microbial_oxid_potential_layer, &
		methane_oxid_depth, o2_oxid_depth)
		!-----------------------------------------------------------------------
		! DESCRIPTION:
		! Oxidation is based on double Michaelis-Mentin kinetics, and is adjusted for low soil moisture.
		! Oxidation participates in the shared O2 competition in methane_tran.
		! The microbial override potential may already include local microbial
		! substrate/O2 limits; methane_tran can still reduce it when roots,
		! decomposition, and methanotrophs compete for the same timestep O2.
		!-----------------------------------------------------------------------

		!-----------------------Argument---------- -----------------------------
		integer , intent(in) :: &
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			patchtype        , &! land patch type
			jwt                    , &! index of the soil layer right above the water table (-)
			sat                       ! 0 = unsaturated; 1 = saturated

		real(r8), intent(in) :: &
			t_soisno (maxsnl+1:nl_soil)    , &! soil temperature (Kelvin)
         fliq     (1:nl_soil)           , &! unfrozen water fraction wliq/(wliq+wice) (-)
			dz_soisno(maxsnl+1:nl_soil)    , &! soil layer thickness (m)
			zi_soisno(maxsnl:nl_soil)      , &! interface level below a "z" level (m)
			smp      (1:nl_soil)   , &! soil matrix potential [mm]
			vol_aqu  (1:nl_soil)   , &! liquid volumetric water content

			conc_o2_aqu_porsl  (1:nl_soil)   , &! aqueous phase O2 conc in each porosity (mol/m3)
			conc_ch4_aqu_porsl (1:nl_soil)   , &! aqueous phase CH4 conc in each porosity (mol/m3)
			microbial_oxid_potential_layer(1:nl_soil) ! optional microbial CH4 oxidation potential (mol/m3/s)

		real(r8), intent(out) :: &
			methane_oxid_depth (1:nl_soil)   , &! CH4 consumption rate via oxidation in each soil layer (mol/m3/s)
			o2_oxid_depth  (1:nl_soil)      ! O2 consumption rate via oxidation in each soil layer (mol/m3/s)

		!-----------------------Local Variables---------------------------------
		integer :: j                              ! indices
		real(r8):: t0                             ! Base temperature for Q10
		real(r8):: oxid_a                         ! Oxidation predicted by method A (temperature & enzyme limited) (mol CH4/m3/s)
		real(r8):: smp_fact                       ! factor for reduction based on soil moisture (unitless)
		real(r8):: k_m_eff                        ! effective DEF_METHANE%k_m
		real(r8):: k_m_o2_eff                     ! effective DEF_METHANE%k_m_o2
		real(r8):: vmax_eff                       ! effective vmax
		real(r8):: oxid_scale                     ! optional landunit-specific oxidation multiplier
		real(r8):: lake_oxid_layer_factor         ! fraction of layer inside lake oxic sediment cap
		real(r8):: lake_oxid_layer_top            ! layer top depth for lake oxic cap (m)
		real(r8):: lake_oxid_layer_bot            ! layer bottom depth for lake oxic cap (m)
		real(r8):: lake_oxid_layer_overlap        ! layer thickness inside lake oxic cap (m)
		logical :: use_microbe_override
		!-----------------------------------------------------------------------
		t0 = tfrz + 12._r8 ! Walter, for Michigan site where the 45 M/h comes from
		use_microbe_override = DEF_METHANE%use_microbial_pools .and. &
			DEF_METHANE%use_microbial_flux_override .and. patchtype /= 4

		! Loop to determine oxidation in each layer
		do j=1,nl_soil
			if (sat == 1 .or. j > jwt) then ! Below the water table
				! Literature (e.g. Bender & Conrad, 1992) suggests lower DEF_METHANE%k_m and vmax for high-CH4-affinity methanotrophs in
				! upland soils consuming ambient methane.
				k_m_eff = DEF_METHANE%k_m
				vmax_eff = DEF_METHANE%vmax_methane_oxid
			else
				k_m_eff = DEF_METHANE%k_m_unsat
				vmax_eff = DEF_METHANE%vmax_oxid_unsat
            if (patchtype == 2 .and. DEF_METHANE%wetland_oxic_cap_kinetics) then
               ! The layers above a wetland's water table are its oxic cap, not
               ! upland soil.  CH4 diffusing up from the saturated zone keeps
               ! them orders of magnitude above ambient air, so the low-affinity
               ! kinetics of the saturated zone apply there as well; the
               ! high-affinity values stay with upland patches, which is the
               ! case Bender & Conrad (1992) describe.  Same choice the lean
               ! line makes through its per-WFT table (B-5-2).
               k_m_eff = DEF_METHANE%k_m
               vmax_eff = DEF_METHANE%vmax_methane_oxid
            end if
			end if
			k_m_o2_eff = DEF_METHANE%k_m_o2
			oxid_scale = 1._r8
			lake_oxid_layer_factor = 1._r8
			if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
				if (DEF_METHANE%lake_vmax_methane_oxid >= 0._r8) &
					vmax_eff = DEF_METHANE%lake_vmax_methane_oxid
				if (DEF_METHANE%lake_k_m_o2 > 0._r8) k_m_o2_eff = DEF_METHANE%lake_k_m_o2
				oxid_scale = DEF_METHANE%lake_oxid_scale
				if (DEF_METHANE%lake_oxic_sediment_depth > 0._r8) then
					if (j == 1) then
						lake_oxid_layer_top = 0._r8
					else
						lake_oxid_layer_top = zi_soisno(j-1)
					endif
					lake_oxid_layer_bot = zi_soisno(j)
					lake_oxid_layer_overlap = max(0._r8, &
						min(DEF_METHANE%lake_oxic_sediment_depth, lake_oxid_layer_bot) - lake_oxid_layer_top)
					lake_oxid_layer_factor = min(1._r8, &
						lake_oxid_layer_overlap / max(dz_soisno(j), 1.e-12_r8))
				endif
			endif

			! Continuous soil-water-potential limitation follows CLM5.
			! Smoothly reduce oxidation under dry matrix-potential stress.
			! Reference: Schnell & King 1996, Figure 3; CLM5 Tech Note 2.25.11
			if (j <= jwt .and. smp(j) < 0._r8) then
				smp_fact = exp(-smp(j)/DEF_METHANE%smp_crit)
			else
				smp_fact = 1._r8
			end if

			oxid_a              = oxid_scale * vmax_eff * vol_aqu(j)* conc_ch4_aqu_porsl(j) / (k_m_eff + conc_ch4_aqu_porsl(j)) &
			! [mol/m3/s]        = [mol/m3/s]   * [-]     [mol/m3-w]    [mol/m3-w]  [mol/m3-w]
				* conc_o2_aqu_porsl(j) / (k_m_o2_eff + conc_o2_aqu_porsl(j)) &
				* DEF_METHANE%q10_methane_oxid ** ((t_soisno(j) - t0) / 10._r8) * smp_fact &
				* lake_oxid_layer_factor

			! For all landunits / levels, prevent oxidation if at or below freezing
         ! (with liquid_fraction_scaling, scale it by the unfrozen water fraction)
         if (DEF_METHANE%liquid_fraction_scaling) then
            oxid_a = oxid_a * fliq(j)
         else
			if (t_soisno(j) <= tfrz) oxid_a = 0._r8
         end if

			if (use_microbe_override .and. &
			    abs(microbial_oxid_potential_layer(j)) < 0.5_r8 * abs(spval)) then
				methane_oxid_depth(j) = max(0._r8, microbial_oxid_potential_layer(j))
			else
				methane_oxid_depth(j) = oxid_a
			endif
			o2_oxid_depth(j) = methane_oxid_depth(j) * 2._r8
		end do
	end subroutine methane_oxid

		!---------------------------------------------------------------------------
	subroutine methane_aere (ipatch, idate,jwt, sat, patchclass, lai, woody_share,&
      z_soisno, dz_soisno, zi_soisno, t_soisno, fliq,&
		rootfr, rootr, etr, grnd_methane_cond_base, c_atm, annsum_npp,&
		annavg_agnpp, annavg_bgnpp, &
		conc_ch4_aqu, conc_ch4_aqu_porsl,conc_ch4_gas_porsl,conc_o2_aqu_porsl,conc_o2_gas_porsl,&
		methane_aere_depth, methane_tran_depth, o2_aere_depth)
		!-----------------------------------------------------------------------
		! DESCRIPTION:
		! Arctic c3 grass (which is often present in fens) and all vegetation in inundated areas is assumed to have
		! some root porosity. Currently, root porosity is allowed to be different for grasses & non-grasses.
		! CH4 diffuses out and O2 diffuses into the soil.  CH4 is also lossed via transpiration, which is both
		! included in the "aere" variables and output separately.  In practice this value is small.
		! By default upland veg. has small 5% porosity but this can be switched to be equal to inundated porosity.
		!-----------------------------------------------------------------------

		!-----------------------Argument----------------------------------------
		integer , intent(in) :: &
			ipatch            , &! caller's per-patch index for vegetation overrides
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			jwt                    , &! index of the soil layer right above the water table (-)
			sat                    , &! 0 = unsatured, 1 = saturated
			patchclass                ! land patch class of USGS classification or others


		real(r8), intent(in) :: &
			lai                    , &! adjusted leaf area index for seasonal variation [-]
         woody_share            , &! share of the conduit at woody_conduit_area, 0 when off [-]

			z_soisno (maxsnl+1:nl_soil)    , &! layer depth (m)
			dz_soisno(maxsnl+1:nl_soil)    , &! layer thickness (m)
			zi_soisno(maxsnl:nl_soil)      , &! interface level below a "z" level (m)

			t_soisno (maxsnl+1:nl_soil)    , &! soil temperature (Kelvin)
         fliq     (1:nl_soil)   , &! unfrozen water fraction wliq/(wliq+wice) (-)
			rootfr   (1:nl_soil)   , &! fraction of roots in each soil layer
			rootr    (1:nl_soil)   , &! effective fraction of roots in each soil layer (SMS method only)
			! rootr here for effective per-layer transpiration, which may not be the same as rootfr
			etr                    , &! transpiration rate [mm/s]
			grnd_methane_cond_base     , &! base tracer conductance before snow/pond resistance [m/s]

			c_atm(3)               , &! CH4, O2, CO2 atmospheric conc  (mol/m3)

			! These variables help us swap between big-leaf and fates boundary conditions
			annsum_npp             , &! annual sum NPP (gC/m2/yr)
			annavg_agnpp           , &! annual avg aboveground NPP (gC/m2/s)
			annavg_bgnpp           , &! annual avg belowground NPP (gC/m2/s)

      ! conc_methane / methane_{prod,oxid,ebul}_depth removed from methane_aere
      ! args: aerenchyma competition is resolved in methane_tran, not here.

			conc_ch4_aqu       (1:nl_soil) , &! aqueous phase CH4 concentration [mol/m3 water]
			conc_ch4_aqu_porsl (1:nl_soil) , &! aqueous phase CH4 conc in each porosity [mol/m3]
			conc_ch4_gas_porsl (1:nl_soil) , &! gas phase CH4 conc in each porosity [mol/m3]
			conc_o2_aqu_porsl  (1:nl_soil) , &! aqueous phase O2 conc in each porosity [mol/m3]
			conc_o2_gas_porsl  (1:nl_soil)    ! gas phase O2 conc in each porosity [mol/m3]


		real(r8), intent(out) :: &
			methane_aere_depth  (1:nl_soil)  , &! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s)
			methane_tran_depth  (1:nl_soil)  , &! CH4 loss rate via transpiration in each soil layer (mol/m3/s)
			o2_aere_depth   (1:nl_soil)     ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s)

		! methane aerenchyma parameters
		real(r8) :: tranloss(1:nl_soil)    ! loss due to transpiration (mol / m3 /s)
		real(r8) :: aere    (1:nl_soil)
		real(r8) :: oxaere  (1:nl_soil)    ! (mol / m3 /s)

		real(r8) :: poros, poros_tiller_real
		real(r8) :: poros_default
      real(r8) :: herb_share, woody_share_pft   ! PFT make-up of a soil patch

		!-----------------------------------------------------------------------
      if (methane_patch_pft_shares(ipatch, herb_share, woody_share_pft)) then
         ! area-weighted porosity of the herbs and the woody plants;
         ! where the trees have a conduit of their own the tillers
         ! keep the grass porosity, as on a wetland tile
         if (woody_share > 0._r8) then
            poros_default = DEF_METHANE%poros_tiller
         else
            poros_default = DEF_METHANE%poros_tiller * &
               (herb_share + (1._r8 - herb_share) * DEF_METHANE%nongrassporosratio)
         endif
      elseif (methane_patch_is_nongrass(patchclass)) then
			   poros_default = DEF_METHANE%poros_tiller * DEF_METHANE%nongrassporosratio
			else
			   poros_default = DEF_METHANE%poros_tiller
			endif
			poros = get_aere_poros(ipatch, poros_default)
		if (sat == 0) then
			poros_tiller_real = poros * DEF_METHANE%unsat_aere_ratio ! unsaturate
			! Apply a proportional floor on the unsat path: the absolute
			! porosmin=0.05 used to silently lift the non-grass branch back
			! up to the grass value (both ended at 0.05 because porosmin
			! caught poros*ratio = 0.3*1/3*0.05/0.3 = 0.0167 -> 0.05).
			! That defeats nongrassporosratio for forest wetlands.  Keep
			! the species contrast by scaling porosmin by the same
			! unsat_aere_ratio.
			poros_tiller_real = max(poros_tiller_real, &
				DEF_METHANE%porosmin * DEF_METHANE%unsat_aere_ratio)
		else
			poros_tiller_real = poros
			poros_tiller_real = max(poros_tiller_real, DEF_METHANE%porosmin)
		endif

      call SiteOxAere(ipatch, idate, jwt,  sat, poros_tiller_real, lai, woody_share, &
         z_soisno, dz_soisno,  zi_soisno,  t_soisno,  fliq,  &
		rootfr, rootr, grnd_methane_cond_base, etr, &
			annsum_npp, annavg_agnpp, annavg_bgnpp, c_atm, conc_ch4_aqu, conc_ch4_aqu_porsl, conc_ch4_gas_porsl, conc_o2_aqu_porsl,conc_o2_gas_porsl,&
			tranloss, aere, oxaere)

		methane_aere_depth = aere
		methane_tran_depth = tranloss
		o2_aere_depth = oxaere

	end subroutine methane_aere


	!---------------------------------------------------------------------------
   subroutine SiteOxAere(ipatch, idate,jwt,  sat, poros_tiller_real, lai, woody_share, &
      z_soisno, dz_soisno,  zi_soisno,  t_soisno,  fliq,  &
			rootfr, rootr, grnd_methane_cond_base, etr, &
		annsum_npp, annavg_agnpp, annavg_bgnpp, c_atm, conc_ch4_aqu, conc_ch4_aqu_porsl, conc_ch4_gas_porsl, conc_o2_aqu_porsl,conc_o2_gas_porsl,&
		tranloss, aere, oxaere)
		!-----------------------------------------------------------------------
		! DESCRIPTION:
		! Site(column) level fluxes for O2 gain rate via
		! aerenchyma and methane losss rates from transpiration
		!-----------------------------------------------------------------------

		!-----------------------Argument----------------------------------------
		integer , intent(in) :: &
			ipatch            , &! caller's per-patch index for vegetation overrides
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			jwt                    , &! index of the soil layer right above the water table (-)
			sat                       ! 0 = unsatured, 1 = saturated

		real(r8), intent(in) :: &
			poros_tiller_real      , &
			lai                    , &! leaf area index [m2/m2]
         woody_share            , &! share of the conduit at woody_conduit_area, 0 when off [-]

			z_soisno (maxsnl+1:nl_soil)    , &! layer depth (m)
			dz_soisno(maxsnl+1:nl_soil)    , &! layer thickness (m)
			zi_soisno(maxsnl:nl_soil)      , &! interface level below a "z" level (m)

			t_soisno (maxsnl+1:nl_soil)    , &! soil temperature (Kelvin)
         fliq     (1:nl_soil)   , &! unfrozen water fraction wliq/(wliq+wice) (-)
			rootfr   (1:nl_soil)   , &! fraction of roots in each soil layer
			rootr    (1:nl_soil)   , &! root resistance of a layer, all layers sum to 1
			grnd_methane_cond_base     , &! base tracer conductance before snow/pond resistance [m/s]
			etr                    , &! transpiration rate [mm/s]

			annsum_npp             , &! annual sum NPP (g C/m2/yr)
			annavg_agnpp           , &! annual average aboveground NPP (g C/m2/s)
			annavg_bgnpp           , &! annual average belowground NPP (g C/m2/s)

			c_atm(3)               , &! CH4, O2, CO2 atmospheric conc  (mol/m3)

			conc_ch4_aqu       (1:nl_soil) , &! aqueous phase CH4 concentration [mol/m3 water]
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
			real(r8) :: secsperyear ! seconds in current model year
			real(r8) :: aere_ch4_resis    ! aerenchyma resistance [s/m]
		real(r8) :: grnd_ch4_resis    ! boundary layer resistance [s/m]
		real(r8) :: aere_o2_resis    ! aerenchyma resistance [s/m]
		real(r8) :: grnd_o2_resis    ! boundary layer resistance [s/m]
		real(r8) :: aerecond    ! aerenchyma conductance [m/s]
		! Per-patch aerenchyma overrides for wetland (patchtype==2):
		! get_aere_* return DEF_METHANE%* defaults for non-wetland patches.
		real(r8) :: tillerC_eff      ! gC per tiller for current patch
		real(r8) :: radius_eff       ! tiller radius for current patch (m)
		real(r8) :: scale_eff        ! scale_factor_aere for current patch
		real(r8), parameter :: smallnumber = 1.e-12_r8
		!-----------------------------------------------------------------------
		! This parameter is poorly constrained and should be done on a patch-specific basis...

		! Attn EK: This calculation of aerenchyma properties is very uncertain. Let's check in once all
		! the new components are in; if there is any tuning to be done to get a realistic global flux,
		! this would probably be the place.  We will have to document clearly in the Tech Note
		! any major changes from the Riley et al. 2011 version. (There are a few other minor ones.)

			tillerC_eff = get_aere_tillerC(ipatch, DEF_METHANE%tiller_C)
			radius_eff  = get_aere_radius (ipatch, DEF_METHANE%aere_radius)
			scale_eff   = get_aere_scale  (ipatch, DEF_METHANE%scale_factor_aere)

			anpp = annsum_npp
			if (ieee_is_nan(anpp) .or. abs(anpp) >= 0.5_r8*abs(spval)) anpp = 0._r8
			if (anpp <= 0._r8 .and. .not. ieee_is_nan(annavg_agnpp) .and. .not. ieee_is_nan(annavg_bgnpp) .and. &
				abs(annavg_agnpp) < 0.5_r8*abs(spval) .and. abs(annavg_bgnpp) < 0.5_r8*abs(spval) .and. &
				(annavg_agnpp + annavg_bgnpp) > 0._r8) then
				! If annual-sum NPP is unavailable, derive annual total NPP from annual-mean rates.
				if (isleapyear(idate(1))) then
					secsperyear = 366._r8 * secspday
				else
					secsperyear = 365._r8 * secspday
				end if
				anpp = (annavg_agnpp + annavg_bgnpp) * secsperyear
				! [gC/m2] = [gC/m2/s] * [s/yr]
			end if
			anpp = max(anpp, 0._r8) ! NPP can be negative b/c of consumption of storage pools

		if (.not. ieee_is_nan(annavg_agnpp) .and. .not. ieee_is_nan(annavg_bgnpp) .and. &
			abs(annavg_agnpp) < 0.5_r8*abs(spval) .and. abs(annavg_bgnpp) < 0.5_r8*abs(spval) .and. &
			annavg_agnpp > 0._r8 .and. annavg_bgnpp > 0._r8) then
			nppratio = annavg_bgnpp / (annavg_agnpp + annavg_bgnpp)
		else
			nppratio = 0.5_r8
		end if

		do j=1,nl_soil
			! Calculate transpiration loss
			if (DEF_METHANE%transpirationloss .and. lai > 0) then
				tranloss(j) = conc_ch4_aqu_porsl(j) * rootr(j)*etr / dz_soisno(j) / 1000._r8
				! [mol/m3/s]= [mol/m3 water] * [-]     *[mm/s]/ [m]        /   [mm/m]
				! Use rootr here for effective per-layer transpiration, which may not be the same as rootfr
				tranloss(j) = max(tranloss(j), 0._r8) ! in case transpiration is pathological
			else
				tranloss(j) = 0._r8
			end if

			! Calculate aerenchyma diffusion.  With liquid_fraction_scaling a frozen
         ! layer is not switched off at tfrz; its conductance is scaled by the
         ! unfrozen water fraction fliq(j) below instead.
         if (j > jwt .and. (t_soisno(j) > tfrz .or. DEF_METHANE%liquid_fraction_scaling) .and. &
             (lai > 0 .or. woody_share * DEF_METHANE%woody_conduit_area > 0._r8)) then ! Below water table
            if (lai > 0) then
				! Estimate area of tillers (see Wania thesis)
				!m_tiller = anpp * r_leaf_root * lai ! (4.17 Wania)
				!m_tiller = 600._r8 * 0.5_r8 * 2._r8  ! used to be 300
				! Note: this calculation is based on Arctic graminoids, and should be refined for woody plants, if not
				! done on a PFT-specific basis.

				! Use the current effective LAI for aerenchyma area, following CLM5.
				m_tiller = anpp * nppratio * lai
				! [g C/m2] = [g C/m2/yr aliased as g C/m2] * [-] * [m2/m2]

				n_tiller = m_tiller / tillerC_eff
				! [tiller/m2] = [g C/m2]/ [g C/tiller]

				area_tiller = scale_eff * n_tiller * poros_tiller_real * PI * radius_eff**2
				! [m2/m2]   = [-]               * [tiller/m2] * [-]       * [-]* [m2/tiller]
            else
               area_tiller = 0._r8
            end if
            ! the woody share of the conduit has the fixed cross-section
            ! woody_conduit_area, whatever the NPP, LAI and phenology: stems
            ! vent pore water CH4 by gas diffusion, at night and without leaves
            ! (Pangala et al. 2013, 2014, 2017). The tillers keep the rest.
            if (woody_share > 0._r8) area_tiller = (1._r8 - woody_share) * area_tiller &
               + woody_share * DEF_METHANE%woody_conduit_area
            ! [m2/m2]   = [-] * [m2/m2] + [-] * [m2/m2]

				aere_ch4_resis = 1._r8/((area_tiller * rootfr(j) * d_con_g(1,1) * 1e-4_r8 / (z_soisno(j)*DEF_METHANE%rob))+smallnumber)
				! [s/m]         = /([m2/m2]          * [-]       * [m2/s]                 / [m]         /[-])
				! Add in boundary layer resistance
				grnd_ch4_resis = 1._r8/(grnd_methane_cond_base+smallnumber)
				aerecond = 1._r8/(aere_ch4_resis + grnd_ch4_resis)
            ! Ice in the pore space blocks the root-soil interface: scale the
            ! conductance by the layer's unfrozen water fraction, the same
            ! quantity methane_prod uses (exactly 1 in an unfrozen layer).
            if (DEF_METHANE%liquid_fraction_scaling) aerecond = aerecond * fliq(j)
				! aerecond = max(aerecond,1.e-8_r8)
				aere(j) = aerecond*(conc_ch4_gas_porsl(j) - c_atm(1)) / dz_soisno(j)
				! [mol/m3/s] = [mol/m3]                      / [m]          / [s/m]
					! Allow negative CH4 aerenchyma fluxes: when pore-air CH4 is
					! below atmospheric CH4, aerenchyma can be an uptake pathway.
					! Downstream methane competition only limits positive sinks.

				! Do oxygen diffusion into layer
				aere_o2_resis = 1._r8/((area_tiller * rootfr(j) * d_con_g(2,1) * 1e-4_r8 / (z_soisno(j)*DEF_METHANE%rob)) + smallnumber)
				grnd_o2_resis = 1._r8/(grnd_methane_cond_base+smallnumber)

				oxaere(j) = -(conc_o2_gas_porsl(j) - c_atm(2)) / (dz_soisno(j)*(aere_o2_resis + grnd_o2_resis)) ![mol/m3-total/s]
				oxaere(j) = max(oxaere(j), 0._r8)
				! Diffusion in is positive; prevent backwards diffusion
            ! O2 rides the same aerenchyma, so it sees the same ice throttle.
            if (DEF_METHANE%liquid_fraction_scaling) oxaere(j) = oxaere(j) * fliq(j)
				if ( .not. DEF_METHANE%use_aereoxid_prog ) then ! fixed aere oxid proportion; will be done in methane_tran
					oxaere(j) = 0._r8
				end if
			else
				aere(j) = 0._r8
				oxaere(j) = 0._r8
			end if ! veg type, below water table, & above freezing
		end do

	end subroutine SiteOxAere

	!---------------------------------------------------------------------------
	subroutine methane_ebul (idate, patchtype, jwt, sat, finundated, deltim, &
		z_soisno, dz_soisno, zi_soisno, forc_pbot, lakedepth, lake_icefrac, &
		t_soisno, wdsrf, conc_methane, conc_ch4_gas_porsl,&
		methane_ebul_depth)
		!-----------------------------------------------------------------------
		! DESCRIPTION:
		! Bubbling is based on temperature & pressure dependent solubility (k_h_cc),
		! with assumed proportion of bubbles
		! which are CH4, and assumed early nucleation at DEF_METHANE%vgc_max sat (Wania).
		! Bubbles are released to the water table surface in methane_tran.
		!-----------------------------------------------------------------------

		!-----------------------Argument---------- -----------------------------
		integer , intent(in) :: &
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			patchtype                  , &! land patch type; 4 = lake/water body
			jwt                        , &! index of the soil layer right above the water table (-)
			sat                        ! 0 = unsaturated; 1 = saturated

		real(r8), intent(in) :: &
			finundated                 , &
			deltim                     , &! land model time step (sec)
			z_soisno (maxsnl+1:nl_soil), &! layer depth (m)
			dz_soisno(maxsnl+1:nl_soil), &! layer thickness (m)
			zi_soisno(maxsnl:nl_soil)  , &! interface level below a "z" level (m)

			forc_pbot                  , &! atm bottom level pressure (or reference height) (Pa)
			lakedepth                  , &! lake depth (m), CTSM hydrostatic pressure over sediments
			lake_icefrac(1:nl_lake)    , &! lake frozen mass fraction [-]
			t_soisno (maxsnl+1:nl_soil), &! soil temperature (Kelvin)
			wdsrf                      , &! depth of surface water [mm]
			conc_methane       (1:nl_soil) , &! CH4 conc in each soil layer (mol/m3)

			conc_ch4_gas_porsl (1:nl_soil)! gas phase CH4 conc in each porosity (mol/m3)

		real(r8), intent(out) :: &
			methane_ebul_depth (1:nl_soil)    ! CH4 loss rate via ebullition in each soil layer (mol/m3/s)

		!-----------------------Local Variables---------------------------------
		integer :: j      ! indices

			real(r8) :: vgc     ! gas phase volumetric CH4 content (m3 CH4/m3 pore air)
			real(r8) :: pressure! sum atmospheric and hydrostatic pressure
			real(r8) :: ebul_timescale
			real(r8), parameter :: smallnumber = 1.e-12_r8
		!-----------------------------------------------------------------------
		! Ebullition follows the CLM5/Riley trigger and post-event target gas fraction.
		ebul_timescale = deltim ! Allow fast bubbling

      ! with the under-ice bubble store on, the bed keeps bubbling under
      ! ice and methane_tran holds the bubbles below the ice
      if (patchtype == 4 .and. DEF_METHANE%allowlakeprod .and. lake_icefrac(1) > 0.1_r8 .and. &
          DEF_METHANE%lake_icebubble_release <= 0._r8) then
			methane_ebul_depth = 0._r8
			return
		endif

		! column loop to estimate ebullition CH4 flux from each soil layer
		do j=1,nl_soil
			if (j  >  jwt .and. t_soisno(j) > tfrz) then ! Ebullition occurs only below the water table
				if (patchtype == 4 .and. DEF_METHANE%allowlakeprod) then
					! CTSM lake path: sediment bubbles must overcome the full
					! lake water-column head above the sediment surface.
					pressure = forc_pbot + denh2o * grav * (z_soisno(j) + max(lakedepth, 0._r8))
				else
					pressure = forc_pbot + denh2o * grav * (z_soisno(j)-zi_soisno(jwt))
					! [Pa]   = [Pa]      + [kg/m3]* [m/s2]* [m]
					! [Pa]   = [N/m2] = [kg]*[m/s2]/[m2] = [kg/m/s2]
				endif
				if (patchtype /= 4 .and. sat == 1 .and. finundated>0._r8) then ! Add ponding pressure head
					! Convert gridcell-mean ponded depth to local inundated depth with a 1% floor.
					pressure = pressure + denh2o * grav * wdsrf/1000._r8/max(finundated, 0.01_r8)
					! [Pa]   = [Pa]     + [kg/m3]* [m/s2]* [mm]/[mm/m]
				end if

				! Compare partial pressure to ambient pressure.
				vgc = conc_ch4_gas_porsl(j) * rgasm * t_soisno(j) / pressure
				! [-]= [mol/m3]           * [Pa*m3/K/mol]*[K]   / [Pa]

					! Trigger bubbles at vgc_max, then relax the layer toward the
					! lower post-event target vgc_max*bubble_f.  This implements
					! the intended two-threshold release; the old single threshold
					! vgc_max*bubble_f fired too early and had no hysteresis window.
					if (vgc > DEF_METHANE%vgc_max) then
						methane_ebul_depth(j) = (vgc - DEF_METHANE%vgc_max * DEF_METHANE%bubble_f) / &
							max(vgc, smallnumber) * conc_methane(j) / ebul_timescale
						! [mol/m3/s]      = [-]                                       * [mol/m3]    / [s]
					else
					methane_ebul_depth(j) = 0._r8
				endif

			else ! above the water table or freezing
				methane_ebul_depth(j) = 0._r8
			endif ! below the water table and not freezing
		end do ! j

	end subroutine methane_ebul


		!---------------------------------------------------------------------------
	subroutine methane_tran (idate,patchtype, &
		lb, snl, jwt, sat, finundated,&
		dlon, dlat, deltim, z_soisno, dz_soisno, zi_soisno,  t_soisno, t_grnd, forc_us, forc_vs, &
		porsl, wliq_soisno, wice_soisno, wdsrf, lakedepth, dz_lake, t_lake, lake_icefrac, bsw, c_atm, methane_prod_depth, o2_aere_depth,&
		cellorg,t_h2osfc, organic_max, k_h_cc, conc_ch4_gas_porsl,conc_ch4_aqu_porsl,conc_o2_gas_porsl,conc_o2_aqu_porsl, &
		vol_aqu,vol_ch4_storage,vol_gas,&
		o2stress, methane_stress, methane_surf_aere, methane_surf_ebul, methane_surf_diff, methane_ebul_tot, &
		methane_oxid_depth, methane_aere_depth, methane_tran_depth, methane_ebul_depth, &
			grnd_methane_cond_base, grnd_methane_cond_effective, o2_oxid_depth, o2_decomp_depth, &
			conc_o2, conc_methane, methane_balance_residual, methane_ch4_clip_credit, methane_surf_diff_phys, &
		o2_cap_loss, o2_cap_gain, lake_water_ch4_stock, lake_water_o2_stock, &
		lake_water_ch4_oxid, lake_sed_ch4_flux, lake_sed_o2_flux, lake_air_o2_flux, lake_icebub_stock, &
		rice_rox_in )
		!-----------------------------------------------------------------------
		! DESCRIPTION:
		! Solves the reaction & diffusion equation for the timestep.
		! 1  "Competition" between processes for CH4 & O2 demand is done.
		! 2  Concentrations are apportioned into gas & liquid fractions;
		!    only the gas fraction is considered for diffusion in unsat.
		! 3  Snow/pond resistance modifies the soil-atmosphere boundary.  Lake
		!    sediment is instead coupled to a prognostic, well-mixed water node;
		!    the water node exchanges independently with the atmosphere.
		! 4  Diffusivity is set based on soil texture and organic matter fraction.
		!    A backward-Euler M-matrix solution is used for positivity and stability.
		! 5  CH4 diffusive flux is calculated and consistency is checked.
		!-----------------------------------------------------------------------

		!-----------------------Argument----------------------------------------
		integer, intent(in) :: &
			idate(3)         , &! current date (year, days of the year, seconds of the day)
			patchtype        	! land patch type (0=soil, 1=urban or built-up, 2=wetland,
										! 3=land ice, 4=land water bodies, 99=ocean


		integer , intent(in) :: &
			lb                , &! lower bound of array (snl+1)
			snl				  , &!  number of snow layers     (-5~-1)
			jwt               , &! index of the soil layer right above the water table (-)
			sat                  ! 0 = unsaturated; 1 = saturated

		real(r8), intent(in) :: &
			finundated                      , &
			dlon   	   				        , &! logitude
			dlat     	   			        , &! latitude

			deltim                  	    , &! land model time step (sec)
			z_soisno (maxsnl+1:nl_soil)   	, &! layer depth (m)
			dz_soisno(maxsnl+1:nl_soil)   	, &! layer thickness (m)
			zi_soisno(maxsnl:nl_soil)   	, &! interface level below a "z" level (m)


			t_soisno (maxsnl+1:nl_soil)    	, &! soil temperature (Kelvin)
			t_grnd                 		    , &! ground surface temperature [k]
			forc_us                         , &! eastward wind speed [m/s], U10 proxy for lake gas exchange
			forc_vs                         , &! northward wind speed [m/s], U10 proxy for lake gas exchange

			porsl             (1:nl_soil)   , &! volumetric soil water at saturation (porosity)
			wliq_soisno(maxsnl+1:nl_soil)	, &! liquid water in layers [kg/m2]
			wice_soisno(maxsnl+1:nl_soil) 	, &! ice lens in layers [kg/m2]
			wdsrf                  		    , &! depth of surface water [mm]
			lakedepth                       , &! static/capacity lake water depth [m]
			dz_lake           (1:nl_lake)   , &! current lake layer thickness [m]
			t_lake            (1:nl_lake)   , &! current lake layer temperature [K]
			lake_icefrac      (1:nl_lake)   , &! lake-layer ice fraction
			bsw               (1:nl_soil)   , &! Clapp and Hornberger "b" (nlevgrnd)

			c_atm(3)               		    , &! CH4, O2, CO2 atmospheric conc  (mol/m3)

			methane_prod_depth    (1:nl_soil)   , &! production of CH4 in each soil layer (mol/m3/s)
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
			vol_ch4_storage   (1:nl_soil)   , &! liquid/ice volume retaining CH4 [m3/m3]
			vol_gas           (1:nl_soil)   , &
			grnd_methane_cond_base              ! base conductance before snow/pond resistance [m/s]

		real(r8), intent(out) :: &
			o2stress          (1:nl_soil)   , &! Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs
			methane_stress         (1:nl_soil)   , &! Ratio of methane available to the total per-timestep methane sinks
			methane_surf_aere                   , &! Total column CH4 aerenchyma (mol/m2/s)
			methane_surf_ebul                   , &! CH4 ebullition to atmosphere (mol/m2/s)
				methane_surf_diff                   , &! CH4 surface flux plus numerical closure (mol/m2/s)
				methane_ebul_tot                    , &! Total column CH4 ebullition (mol/m2/s)
				methane_balance_residual            , &! numerical closure flux credited to methane_surf_diff (mol/m2/s)
				methane_ch4_clip_credit             , &! negative CH4 clip credited to methane_surf_diff (mol/m2/s)
				methane_surf_diff_phys              , &! pure diffusive flux before clip/residual (mol/m2/s)
				o2_cap_loss                        , &! retained diagnostic; no upper O2 cap is applied (mol/m2/s)
				o2_cap_gain                         , &! O2 added by post-solve nonnegative floor (mol/m2/s)
				lake_water_ch4_oxid                , &! water-column CH4 oxidation (mol/m2/s)
				lake_sed_ch4_flux                   , &! sediment-to-water CH4 flux (mol/m2/s)
				lake_sed_o2_flux                    , &! sediment-to-water O2 flux (mol/m2/s)
				lake_air_o2_flux                    , &! water-to-atmosphere O2 flux (mol/m2/s)
				grnd_methane_cond_effective            ! effective CH4 conductance including snow/pond resistance [m/s]

		real(r8), intent(inout) :: &
			methane_oxid_depth    (1:nl_soil)   , &! InOut: CH4 consumption rate via oxidation in each soil layer (mol/m3/s)
			methane_aere_depth    (1:nl_soil)   , &! InOut: CH4 loss rate via aerenchyma in each soil layer (mol/m3/s)
			methane_tran_depth    (1:nl_soil)   , &! InOut: CH4 loss rate via transpiration in each soil layer (mol/m3/s)
			methane_ebul_depth    (1:nl_soil)   , &! InOut: CH4 loss rate via ebullition in each soil layer (mol/m3/s)
			o2_oxid_depth     (1:nl_soil)   , &! InOut: O2 loss rate via ebullition in each soil layer (mol/m3/s)
			o2_decomp_depth   (1:nl_soil)   , &! InOut: O2 consumption during decomposition in each soil layer (mol/m3/s)

			conc_o2           (1:nl_soil)   , &! InOut: O2 conc in each soil layer (mol/m3)
			conc_methane      (1:nl_soil)   , &! InOut: CH4 conc in each soil layer (mol/m3)
			lake_water_ch4_stock             , &! well-mixed water CH4 inventory (mol/m2)
			lake_water_o2_stock              , &! well-mixed water O2 inventory (mol/m2)
			lake_icebub_stock                  ! CH4 of bubbles held under the lake ice (mol/m2)

		!-----------------------Local Variables------------------------------
		integer :: j,s,i			                                               ! indices
		integer :: jtop                                                        ! top level at each column
		real(r8) :: at (0:nl_soil)                     ! "a" vector for tridiagonal matrix
		real(r8) :: bt (0:nl_soil)                     ! "b" vector for tridiagonal matrix
		real(r8) :: ct (0:nl_soil)                     ! "c" vector for tridiagonal matrix
		real(r8) :: rt (0:nl_soil)                     ! "r" vector for tridiagonal solution
		real(r8) :: f_a                                                        ! air-filled fraction of available pore space
		real(r8) :: diffus (0:nl_soil)                 ! diffusivity (m2/s)
		real(r8) :: dzj
		real(r8) :: dp1_zp1 (0:nl_soil)                ! diffusivity/delta_z for next j
		real(r8) :: dm1_zm1 (0:nl_soil)                ! diffusivity/delta_z for previous j
		real(r8) :: t_soisno_c                                                 ! soil temperature   (maxsnl+1:nl_soil)
			real(r8) :: deficit                                                    ! mol CH4 /m^2 that must be subtracted from diffusive flux to atm. to make up
			real(r8) :: ch4_clip_deficit_col                                       ! column-integrated nonnegative correction (mol CH4/m2/timestep)
		! for keeping concentrations always above zero
		logical, save :: warned_ch4_clip = .false.
		logical, save :: warned_o2_cap   = .false.
		logical, save :: warned_ch4_competition_demand = .false.
		logical, save :: warned_ch4_competition_left   = .false.
		logical, save :: warned_o2_competition_demand  = .false.
		logical, save :: warned_o2_competition_left    = .false.
		logical, save :: warned_ch4_diffusion_residual = .false.
		logical, save :: warned_lake_zero_depth        = .false.
		real(r8) :: conc_ch4_bef(1:nl_soil)            ! concentration at the beginning of the timestep
		real(r8) :: err_methane                            ! Error (Mol CH4 /m^2) [+ = too much CH4]
		real(r8) :: conc_ch4_rel(0:nl_soil)            ! Concentration per volume of air or water
		real(r8) :: conc_o2_rel(0:nl_soil)             ! Concentration per volume of air or water
		real(r8) :: conc_ch4_rel_old(0:nl_soil)        ! Concentration during last Crank-Nich. loop
			real(r8), parameter :: smallnumber = 1.e-12_r8
			! Use one positivity-preserving transport discretization for both
			! reactive gases.  A Crank-Nicolson O2 step could overshoot sharp
			! gradients and then rely on a mass-creating nonnegative clip.
			logical, parameter :: backward_euler_transport = .true.
			real(r8) :: snowdiff                                                   ! snow diffusivity (m^2/s)
			real(r8) :: snow_resis                           ! Cumulative Snow resistance (s/m). Also includes
			real(r8) :: pond_resis                                                    ! Additional resistance from ponding, up to pondmx water on top of top soil layer (s/m)
			real(r8) :: pondz                                                      ! Depth of ponding (m)
			real(r8) :: ponddiff                                                   ! Pondwater diffusivity (m^2/s)
			real(r8) :: lake_exchange_vel                                          ! Lake piston velocity [m/s]
			real(r8) :: lake_k600_cmhr                                             ! Cole & Caraco k600 [cm/hr]
			real(r8) :: lake_schmidt_ch4                                           ! CH4 Schmidt number in freshwater [-]
			real(r8) :: lake_wind                                                  ! Wind speed [m/s], using forcing wind as U10 proxy
		real(r8) :: liquid_volfrac, ice_volfrac, open_pore_volfrac
		real(r8) :: spec_grnd_cond(1:ngases)            ! species grnd conductance (s/m)
		real(r8) :: airfrac                                                    ! air fraction in snow
		real(r8) :: waterfrac                                                  ! water fraction in snow
		real(r8) :: icefrac                                                    ! ice fraction in snow
		real(r8) :: epsilon_t (1:nl_soil,1:ngases)
		real(r8) :: source (1:nl_soil,1:ngases)        ! source
		real(r8) :: om_frac                                                    ! organic matter fraction
			real(r8) :: o2demand, methane_demand, methane_supply                        ! mol/m^3/s
			real(r8) :: o2_before_cap, o2_cap_loss_col, o2_cap_gain_col
				real(r8) :: aere_oxid_flux
      logical , intent(in), optional :: rice_rox_in   ! paddy rice column with rhizosphere oxidation
      logical  :: rice_rox
      real(r8) :: rice_rox_flux
		logical :: is_lake_water
		logical :: lake_stratified          ! water column stably stratified
		real(r8) :: lake_rho_top, lake_rho_bot
		integer  :: lake_jbot, jj_lake
		real(r8) :: lake_total_depth, lake_liquid_depth, lake_water_storage_depth
		real(r8) :: lake_water_temp, lake_temp_weight, lake_temp_weight_sum
		real(r8) :: lake_sed_water_cond, lake_water_stock_old
		real(r8) :: lake_ch4_stock_bef, lake_water_oxid_potential
		real(r8) :: lake_vmax_oxid, lake_k_m_o2_eff
		real(r8) :: lake_ch4_conc, lake_o2_conc
		real(r8) :: lake_k_m_eff            ! CH4 constant of water-node oxidation (mol/m3)
		real(r8) :: lake_water_sod          ! sediment O2 demand on the water node (mol/m2/s)
		real(r8) :: lake_t_bot              ! temperature of the deepest lake layer (K)
		real(r8) :: lake_ebul_dissol        ! ebullition dissolved into the water node (mol/m2/s)
		real(r8) :: lake_icebub_stock_bef   ! under-ice bubble store at the start of the step (mol/m2)
		real(r8) :: lake_icebub_trap        ! ebullition held under the ice (mol/m2/s)
		real(r8) :: lake_icebub_out         ! CH4 leaving the store this step (mol/m2)
		real(r8) :: lake_icebub_release     ! release from the store (mol/m2/s)
		real(r8) :: lake_icebub_dissol      ! part of the release dissolved into the water node (mol/m2/s)
		logical  :: lake_icebub_trapping    ! bubbles reaching the surface go under the ice this step
		real(r8) :: lake_hyps_dz, lake_hyps_a_top, lake_hyps_a_bot
		real(r8) :: lake_hyps_reach, lake_hyps_area
		integer  :: k_hyps
		integer , parameter :: n_lake_hyps = 20          ! depth bands of the lake basin
		real(r8), parameter :: lake_hyps_shape = 2._r8   ! FLaMe valley-shape parameter s
		real(r8), parameter :: o2_molar_mass = 31.998_r8 ! g O2 per mol
				!-----------------------------------------------------------------------
				methane_balance_residual = 0._r8
				methane_ch4_clip_credit = 0._r8
      rice_rox = .false.
      IF (present(rice_rox_in)) rice_rox = rice_rox_in
				o2_cap_loss = 0._r8
				o2_cap_gain = 0._r8
		lake_water_ch4_oxid = 0._r8
		lake_sed_ch4_flux = 0._r8
		lake_sed_o2_flux = 0._r8
		lake_air_o2_flux = 0._r8
		lake_water_sod = 0._r8
		lake_ebul_dissol = 0._r8
		lake_icebub_trap = 0._r8
		lake_icebub_out = 0._r8
		lake_icebub_release = 0._r8
		lake_icebub_dissol = 0._r8
		lake_icebub_trapping = .false.
		lake_icebub_stock_bef = lake_icebub_stock

		is_lake_water = patchtype == 4 .and. DEF_METHANE%allowlakeprod
		lake_total_depth = 0._r8
		lake_liquid_depth = 0._r8
		lake_water_storage_depth = 0._r8
		lake_water_temp = 0._r8
		lake_temp_weight = 0._r8
		lake_temp_weight_sum = 0._r8
		lake_ch4_stock_bef = lake_water_ch4_stock
		if (is_lake_water) then
			lake_total_depth = sum(max(dz_lake(1:nl_lake), 0._r8))
			lake_liquid_depth = sum(max(dz_lake(1:nl_lake), 0._r8) * &
				(1._r8 - min(max(lake_icefrac(1:nl_lake), 0._r8), 1._r8)))
			do j = 1, nl_lake
				lake_temp_weight = max(dz_lake(j), 0._r8) * &
					(1._r8 - min(max(lake_icefrac(j), 0._r8), 1._r8))
				if (.not. ieee_is_nan(t_lake(j)) .and. t_lake(j) > 150._r8 .and. &
				    t_lake(j) < 350._r8) then
					lake_water_temp = lake_water_temp + lake_temp_weight * t_lake(j)
					lake_temp_weight_sum = lake_temp_weight_sum + lake_temp_weight
				endif
			enddo
			if (lake_temp_weight_sum > smallnumber) then
				lake_water_temp = lake_water_temp / lake_temp_weight_sum
			else
				lake_water_temp = t_h2osfc
			endif
			if (lake_total_depth <= smallnumber .and. lakedepth > 0._r8) then
				lake_total_depth = lakedepth
				lake_liquid_depth = lakedepth * &
					(1._r8 - min(max(lake_icefrac(1), 0._r8), 1._r8))
			endif
			! Dissolved concentration is referenced only to liquid water.  At full
			! freeze the conserved water-column inventory is immobile because all
			! exchange conductances are zero; the small numerical capacity keeps the
			! implicit row nonsingular without treating ice as dissolved storage.
			lake_water_storage_depth = max(lake_liquid_depth, smallnumber)

			lake_ch4_conc = max(lake_water_ch4_stock, 0._r8) / lake_water_storage_depth
			lake_o2_conc = max(lake_water_o2_stock, 0._r8) / lake_water_storage_depth
			lake_vmax_oxid = DEF_METHANE%vmax_methane_oxid
			if (DEF_METHANE%lake_vmax_methane_oxid >= 0._r8) &
				lake_vmax_oxid = DEF_METHANE%lake_vmax_methane_oxid
			lake_k_m_o2_eff = DEF_METHANE%k_m_o2
			if (DEF_METHANE%lake_k_m_o2 > 0._r8) lake_k_m_o2_eff = DEF_METHANE%lake_k_m_o2
			! lake water may take its own CH4 constant
			lake_k_m_eff = DEF_METHANE%k_m
			if (DEF_METHANE%lake_k_m > 0._r8) lake_k_m_eff = DEF_METHANE%lake_k_m
			lake_water_oxid_potential = 0._r8
			if (lake_liquid_depth > smallnumber .and. lake_water_temp > tfrz) then
				lake_water_oxid_potential = DEF_METHANE%lake_oxid_scale * lake_vmax_oxid * &
					lake_ch4_conc / (lake_k_m_eff + lake_ch4_conc) * &
					lake_o2_conc / (lake_k_m_o2_eff + lake_o2_conc) * &
					DEF_METHANE%q10_methane_oxid ** ((lake_water_temp - (tfrz + 12._r8)) / 10._r8) * &
					lake_liquid_depth
			endif
			if (deltim > 0._r8) then
				lake_water_ch4_oxid = min(max(lake_water_oxid_potential, 0._r8), &
					max(lake_water_ch4_stock, 0._r8) / deltim, &
					max(lake_water_o2_stock, 0._r8) / (2._r8 * deltim))
			endif
			lake_water_ch4_stock = max(0._r8, lake_water_ch4_stock - lake_water_ch4_oxid * deltim)
			lake_water_o2_stock = max(0._r8, lake_water_o2_stock - 2._r8 * lake_water_ch4_oxid * deltim)

			! the lake bed draws O2 from the water node,
			! sod20 * 1.065**(T_bot - 20 C) * O2 / (K_O2 + O2), with T_bot the
			! deepest lake layer and K_O2 that of the water CH4 oxidation.  It
			! takes what the CH4 oxidation left, is counted with the sediment-
			! water O2 exchange (lake_sed_o2_flux) after the transport solve, and
			! keeps acting under ice, where the bottom water stays liquid.
			if (DEF_METHANE%lake_sod20 > 0._r8 .and. lake_liquid_depth > smallnumber .and. &
			    deltim > 0._r8) then
				lake_jbot = 1
				do jj_lake = 1, nl_lake
					if (dz_lake(jj_lake) > smallnumber) lake_jbot = jj_lake
				end do
				lake_t_bot = t_lake(lake_jbot)
				if (ieee_is_nan(lake_t_bot) .or. lake_t_bot <= 150._r8 .or. lake_t_bot >= 350._r8) &
					lake_t_bot = lake_water_temp
				lake_o2_conc = max(lake_water_o2_stock, 0._r8) / lake_water_storage_depth
				! [mol/m2/s] = [g O2/m2/d] / ([g O2/mol] * [s/d])
				lake_water_sod = DEF_METHANE%lake_sod20 / (o2_molar_mass * secspday) * &
					1.065_r8 ** (lake_t_bot - 293.15_r8) * &
					lake_o2_conc / (lake_k_m_o2_eff + lake_o2_conc)
				lake_water_sod = min(lake_water_sod, max(lake_water_o2_stock, 0._r8) / deltim)
				lake_water_o2_stock = max(0._r8, lake_water_o2_stock - lake_water_sod * deltim)
			endif
		endif

				! When use_aereoxid_prog=.false. (CTSM legacy fixed-fraction root-zone
			! aerenchyma oxidation) we MUST fold aere_oxid_flux into the
			! methanotrophic oxidation budget BEFORE the O2/CH4 competition loop.
			! Previously this block sat after the competition (just before
			! `source(j,1)`), which let the additional CH4 oxidation and 2:1 O2
			! demand bypass the stress caps — pushing total O2 sink beyond what
			! is available in low-O2 soils.  Now the augmented demands enter the
			! competition like any other sink and get scaled by o2stress/methane_stress.
			do j = 1,nl_soil
				if ( .not. DEF_METHANE%use_aereoxid_prog .and. methane_aere_depth(j) > 0._r8 ) then
					aere_oxid_flux = DEF_METHANE%aereoxid * methane_aere_depth(j)
					methane_oxid_depth(j)  = methane_oxid_depth(j)  + aere_oxid_flux
					o2_oxid_depth(j)       = o2_oxid_depth(j)       + 2._r8 * aere_oxid_flux
					methane_aere_depth(j)  = methane_aere_depth(j)  - aere_oxid_flux
				end if
			end do

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

			! Signed plant uptake is a same-layer source; only positive fluxes compete as sinks.
			methane_supply = conc_methane(j) / deltim + methane_prod_depth(j) - &
				min(methane_aere_depth(j), 0._r8) - min(methane_tran_depth(j), 0._r8)
			methane_demand = methane_oxid_depth(j) + max(methane_aere_depth(j), 0._r8) + &
				max(methane_tran_depth(j), 0._r8) + methane_ebul_depth(j)
			if (methane_demand > 0._r8) then
				methane_stress(j) = min(methane_supply / methane_demand, 1._r8)
			else
				methane_stress(j) = 1._r8
			end if

			! Resolve methane oxidation
			if (o2stress(j) < 1._r8 .or. methane_stress(j) < 1._r8) then
				if (methane_stress(j) <= o2stress(j)) then
					! methane limited
					if (o2stress(j) < 1._r8) then
						! Recalculate oxygen limitation
						o2demand = o2_decomp_depth(j)
						if (o2demand > 0._r8) then
							o2stress(j) = min((conc_o2(j)/deltim + o2_aere_depth(j) - methane_stress(j)*o2_oxid_depth(j))/o2demand, 1._r8)
						else
							o2stress(j) = 1._r8
						end if
					end if
					! Reset oxidation
					methane_oxid_depth(j) = methane_oxid_depth(j) * methane_stress(j)
					o2_oxid_depth(j)  = o2_oxid_depth(j) * methane_stress(j)
				else
					! oxygen limited
					if (methane_stress(j) < 1._r8) then
						! Recalculate methane limitation
						methane_demand = max(methane_aere_depth(j), 0._r8) + &
							max(methane_tran_depth(j), 0._r8) + methane_ebul_depth(j)
						if (methane_demand > 0._r8) then
							methane_stress(j) = min((methane_supply - &
								o2stress(j)*methane_oxid_depth(j)) / methane_demand, 1._r8)
						else
							methane_stress(j) = 1._r8
						end if
					end if
					! Reset oxidation
					methane_oxid_depth(j) = methane_oxid_depth(j) * o2stress(j)
					o2_oxid_depth(j) = o2_oxid_depth(j) * o2stress(j)
				end if
			end if

			! Reset non-methanotroph demands
				IF (methane_aere_depth(j) > 0._r8) methane_aere_depth(j) = methane_aere_depth(j) * methane_stress(j)
				IF (methane_tran_depth(j) > 0._r8) methane_tran_depth(j) = methane_tran_depth(j) * methane_stress(j)
			methane_ebul_depth(j) = methane_ebul_depth(j) * methane_stress(j)
			o2_decomp_depth(j) = o2_decomp_depth(j) * o2stress(j)

         ! rhizosphere oxidation on paddy rice. A fixed share
         ! of the CH4 entering the aerenchyma is oxidised at the roots by the O2
         ! the roots release (CLM4Me's aereoxid meaning): that O2 comes
         ! from inside the plant, so it is not drawn from the layer's O2, which
         ! root respiration and aerobic decomposition have used up by now in a
         ! flooded paddy (taking it from there oxidises next to nothing).
         ! With rox_on_rooted_production the share acts on the production of
         ! the rooted layers instead (methane_root_oxid), not here.
         IF (rice_rox .and. methane_aere_depth(j) > 0._r8 .and. &
             .not. DEF_METHANE%rox_on_rooted_production) THEN
            rice_rox_flux = DEF_METHANE%rice_aereoxid * methane_aere_depth(j)
            methane_aere_depth(j) = methane_aere_depth(j) - rice_rox_flux
            methane_oxid_depth(j) = methane_oxid_depth(j) + rice_rox_flux
         ENDIF
         ! the same at the roots of wetland plants and of the
         ! plants of a flooded soil subcolumn
         IF ((patchtype == 2 .or. (patchtype == 0 .and. sat == 1 .and. .not. rice_rox)) &
             .and. DEF_METHANE%wetland_aereoxid >= 0._r8 .and. methane_aere_depth(j) > 0._r8 &
             .and. .not. DEF_METHANE%rox_on_rooted_production) THEN
            rice_rox_flux = DEF_METHANE%wetland_aereoxid * methane_aere_depth(j)
            methane_aere_depth(j) = methane_aere_depth(j) - rice_rox_flux
            methane_oxid_depth(j) = methane_oxid_depth(j) + rice_rox_flux
         ENDIF


		end do !j


		! Accumulate ebullition to place in first layer above water table, or directly to atmosphere
		methane_ebul_tot = 0._r8
		do j = 1,nl_soil
			methane_ebul_tot = methane_ebul_tot + methane_ebul_depth(j) * dz_soisno(j)
		end do

		! bubbles rising from the lake bed lose CH4 to the water
		! on the way up.  The bed area follows the valley-shaped basin of FLaMe,
		! A(z)/A0 = 1/2 + atanh((1 - 2z/zmax) tanh s)/(2s), s = 2, zmax twice the
		! lake depth, cut into n_lake_hyps equal depth bands; the CH4 of a band's
		! bubbles reaches the air as exp(-z/h_b) averaged over the band.  The rest
		! enters the water node as a source of the CH4 solve below and is taken
		! off the ebullition to the air.  The sediment column loses the full
		! ebullition either way.
		if (is_lake_water .and. jwt == 0 .and. DEF_METHANE%lake_bubble_dissol_depth > 0._r8 .and. &
		    lake_liquid_depth > smallnumber .and. methane_ebul_tot > 0._r8) then
			lake_hyps_dz = 2._r8 * lake_total_depth / real(n_lake_hyps, r8)
			lake_hyps_reach = 0._r8
			lake_hyps_area = 0._r8
			lake_hyps_a_top = 1._r8
			do k_hyps = 1, n_lake_hyps
				lake_hyps_a_bot = 0.5_r8 + atanh((1._r8 - 2._r8 * real(k_hyps, r8) / &
					real(n_lake_hyps, r8)) * tanh(lake_hyps_shape)) / (2._r8 * lake_hyps_shape)
				! band mean of exp(-z/h_b) between the band's top and bottom depths
				lake_hyps_reach = lake_hyps_reach + (lake_hyps_a_top - lake_hyps_a_bot) * &
					DEF_METHANE%lake_bubble_dissol_depth / lake_hyps_dz * &
					(exp(-real(k_hyps - 1, r8) * lake_hyps_dz / DEF_METHANE%lake_bubble_dissol_depth) - &
					 exp(-real(k_hyps, r8) * lake_hyps_dz / DEF_METHANE%lake_bubble_dissol_depth))
				lake_hyps_area = lake_hyps_area + (lake_hyps_a_top - lake_hyps_a_bot)
				lake_hyps_a_top = lake_hyps_a_bot
			enddo
			if (lake_hyps_area > smallnumber) then
				lake_ebul_dissol = methane_ebul_tot * &
					(1._r8 - min(max(lake_hyps_reach / lake_hyps_area, 0._r8), 1._r8))
			endif
		endif

		! while the top lake layer is more than 0.1 ice, the
		! bubbles that reach the surface (what the dissolution on the rise leaves of them) are held under
		! the ice in lake_icebub_stock instead of entering the air.  Once the ice
		! fraction is 0.1 or below the store empties at lake_icebubble_release
		! (d-1), integrated exactly over the step; the share lake_icebubble_dissol
		! of what leaves enters the water node as a source of the CH4 solve below
		! (then water oxidation and gas exchange), the rest joins the ebullition to
		! the air.  The store is counted in the CH4 balance below.
		if (is_lake_water .and. jwt == 0 .and. DEF_METHANE%lake_icebubble_release > 0._r8 .and. &
		    deltim > 0._r8) then
			if (lake_icefrac(1) > 0.1_r8) then
				lake_icebub_trapping = .true.
				lake_icebub_trap = methane_ebul_tot - lake_ebul_dissol
				lake_icebub_stock = lake_icebub_stock + lake_icebub_trap * deltim
			elseif (lake_icebub_stock > 0._r8) then
				lake_icebub_out = lake_icebub_stock * &
					(1._r8 - exp(-DEF_METHANE%lake_icebubble_release * deltim / secspday))
				lake_icebub_stock = lake_icebub_stock - lake_icebub_out
				lake_icebub_release = lake_icebub_out / deltim
				if (lake_liquid_depth > smallnumber) &
					lake_icebub_dissol = DEF_METHANE%lake_icebubble_dissol * lake_icebub_release
			endif
		endif

		! Set the source term for each species (no need to do j=0, since epsilon_t and source not used there)
		! Note that because of the semi-implicit diffusion and the 30 min timestep combined with explicit
		! sources, occasionally negative concentration will result. In this case it is brought to zero and the
		! surface flux is adjusted to conserve. This results in some inaccuracy as compared to a shorter timestep
		! or iterative solution.
		do j = 1,nl_soil

			! NOTE: The fixed-fraction aerenchyma oxidation that used to live
			! here is now applied BEFORE the O2/CH4 competition loop above so
			! it participates in the stress budget.

			source(j,1) = methane_prod_depth(j) - methane_oxid_depth(j) - &
					methane_aere_depth(j) - methane_tran_depth(j) - &
					methane_ebul_depth(j) ! [mol/m3-total/s]

			! aerenchyma added to surface flux below
			! ebul added to soil depth just above WT
			! Original 1e-12 mol/m3/s strict abort tripped on roundoff at multiple
			! patches when only_wetland=.false. exposes wider patch population
			! Replace fatal with rate-limited warning (one per rank
			! per check); source/conc are recomputed each step from physics so a
			! small floating-point imbalance does not corrupt the long-term budget.
			if (source(j,1) + conc_methane(j) / deltim < -1.e-6_r8 .and. .not. warned_ch4_competition_demand)then
				write(6,*) 'WARNING: methane competition demand > 1e-6 mol/m3/s, j:', &
						source(j,1) + conc_methane(j) / deltim, j
				write(6,*)'Lat,Lon=',dlat,dlon
				warned_ch4_competition_demand = .true.
			else if (methane_stress(j) < 1._r8 .and. source(j,1) + conc_methane(j) / deltim > 1.e-6_r8 &
			         .and. .not. warned_ch4_competition_left) then
				write(6,*) 'WARNING: methane limited yet leftover > 1e-6 mol/m3/s, j:', &
						source(j,1) + conc_methane(j) / deltim, j
				write(6,*)'Lat,Lon=',dlat,dlon
				warned_ch4_competition_left = .true.
			end if

			source(j,2) = -o2_oxid_depth(j) - o2_decomp_depth(j) + o2_aere_depth(j) ! O2 [mol/m3/s]

			if (source(j,2) + conc_o2(j) / deltim < -1.e-6_r8 .and. .not. warned_o2_competition_demand) then
				write(6,*) 'WARNING: oxygen competition demand > 1e-6 mol/m3/s, j:', &
						source(j,2) + conc_o2(j) / deltim, j
				write(6,*)'Lat,Lon=',dlat,dlon
				warned_o2_competition_demand = .true.
			else if (o2stress(j) < 1._r8 .and. source(j,2) + conc_o2(j) / deltim > 1.e-6_r8 &
			         .and. .not. warned_o2_competition_left) then
				write(6,*) 'WARNING: oxygen limited yet leftover > 1e-6 mol/m3/s, j:', &
						source(j,2) + conc_o2(j) / deltim, j
				write(6,*)'Lat,Lon=',dlat,dlon
				warned_o2_competition_left = .true.
			end if
			conc_ch4_bef(j) = conc_methane(j) !For Balance Check

		enddo ! j


		! Accumulate aerenchyma to add directly to atmospheric flux
		methane_surf_aere = 0._r8
		do j = 1,nl_soil
			methane_surf_aere = methane_surf_aere + methane_aere_depth(j) * dz_soisno(j)
		enddo

		! Bubbles entering an unsaturated column dissolve immediately above
		! the water table, then remain available for oxidation and diffusion.
		! With the water table at the surface they pass directly to the air.
		if (jwt /= 0) then
			source(jwt,1) = source(jwt,1) + methane_ebul_tot/dz_soisno(jwt)
		endif

		! Calculate concentration relative to m^3 of air or water: needed for the diffusion
		do j = 0,nl_soil
			if (j == 0) then
				if (is_lake_water) then
					conc_ch4_rel(j) = lake_water_ch4_stock / lake_water_storage_depth
					conc_o2_rel(j)  = lake_water_o2_stock / lake_water_storage_depth
				else
					conc_ch4_rel(j) = c_atm(1)
					conc_o2_rel(j)  = c_atm(2)
				endif
			else
					if (porsl(j) <= smallnumber) then
						do s = 1, 2
							epsilon_t(j,s) = smallnumber
						end do
						conc_ch4_rel(j) = 0._r8
						conc_o2_rel(j)  = 0._r8
					else if (j <= jwt) then  ! Above the WT
						! Both gases use mobile gas/liquid pore space; ice stores neither.
						epsilon_t(j,1) = max(vol_gas(j) + k_h_cc(j,1)*vol_ch4_storage(j), smallnumber)
						epsilon_t(j,2) = max(vol_gas(j) + k_h_cc(j,2)*vol_aqu(j), smallnumber)
					conc_ch4_rel(j) = conc_methane(j)/epsilon_t(j,1)
					conc_o2_rel(j)  = conc_o2_gas_porsl(j)
					! Partition between the liquid and gas phases. The gas phase will drive the diffusion.
				else ! Below the WT
					! Below the water table, mobile liquid controls storage/diffusion.
					epsilon_t(j,1) = max(vol_ch4_storage(j), smallnumber)
					epsilon_t(j,2) = max(vol_aqu(j), smallnumber)
					conc_ch4_rel(j) = conc_methane(j)/epsilon_t(j,1)
					conc_o2_rel(j)  = conc_o2(j) /epsilon_t(j,2)
				end if
			end if
		end do


		! Loop over species
		do s = 1, 2 ! transport CH4 and O2; CO2 is not prognosed here
			lake_exchange_vel = 0._r8
			! Apply snow and pond resistance to the immutable base conductance.
			! The effective result is diagnostic output only and never feeds the next step.
			do j = maxsnl + 1,0
				if (j == maxsnl + 1) then
					! Needed to prevent overflow when ground is frozen, e.g. for lakes
					snow_resis = 0._r8
				end if

				! Add snow resistance
				if (j >= snl + 1) then
					! For the ice layer, all = ice + water + air, no soil
					t_soisno_c = t_soisno(j) - tfrz
					icefrac = wice_soisno(j)/denice/dz_soisno(j)
					! [-]   = [kg/m2]       / [kg/m3] / [m]
					waterfrac = wliq_soisno(j)/denh2o/dz_soisno(j)
					airfrac = max(1._r8 - icefrac - waterfrac, 0._r8)

					! Calculate snow diffusivity
					if (airfrac > 0.05_r8) then
						! Use Millington-Quirk Expression, as hydraulic properties (bsw) not available
						snowdiff = (d_con_g(s,1) + d_con_g(s,2)*t_soisno_c) * 1.e-4_r8 * &
								airfrac**(10._r8/3._r8) / (airfrac+waterfrac)**2 &
								* DEF_METHANE%scale_factor_gasdiff
					else !solute diffusion in water only, airfrac -> 0
						snowdiff = (airfrac+waterfrac)**DEF_METHANE%satpow * (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
								* DEF_METHANE%scale_factor_liqdiff
					end if

					snowdiff = max(snowdiff, smallnumber)
					snow_resis = snow_resis + dz_soisno(j)/snowdiff
				end if

				if (j == 0) then ! final loop
					! Add pond resistance
					pond_resis = 0._r8

					! First old pond formulation up to pondmx
					liquid_volfrac = wliq_soisno(1)/(dz_soisno(1)*denh2o)
					ice_volfrac = wice_soisno(1)/(dz_soisno(1)*denice)
					open_pore_volfrac = max(porsl(1) - ice_volfrac, 0._r8)
					if (patchtype /= 4 .and. snl == 0 .and. liquid_volfrac > open_pore_volfrac) then
						t_soisno_c = t_soisno(1) - tfrz
						if (t_soisno(1) <= tfrz) then
							ponddiff = (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
									* (wliq_soisno(1)/denh2o+smallnumber)/ &
									(wliq_soisno(1)/denh2o+wice_soisno(1)/denice+smallnumber) &
									* DEF_METHANE%scale_factor_liqdiff
						else ! Unfrozen
							ponddiff = (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
									* DEF_METHANE%scale_factor_liqdiff
						end if
						pondz = dz_soisno(1) * (liquid_volfrac - open_pore_volfrac)
						pond_resis = pondz / ponddiff
					end if

					! Now add new wdsrf form
					if (patchtype /= 4 .and. sat == 1 .and. finundated>0._r8) then
						if (t_h2osfc >= tfrz) then
							t_soisno_c = t_h2osfc - tfrz
							ponddiff = (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
									* DEF_METHANE%scale_factor_liqdiff
							! Guard against fill-and-spill transients with very small inundated area.
							pondz = wdsrf / 1000._r8/max(finundated, 0.01_r8) ! Assume all wdsrf corresponds to sat area
							! [m] = [mm]  /  [mm/m]
							pond_resis = pond_resis + pondz / ponddiff
						else if (wdsrf/max(finundated, 0.01_r8) > DEF_METHANE%capthick) then
							! assume surface ice is impermeable
							pond_resis = pond_resis + 1._r8/smallnumber
						end if
					end if

					! Lake air exchange is a separate water-side piston boundary;
					! it must not be put in series with sediment diffusion.  The
					! prognostic j=0 water node below supplies the residence time.
					if (is_lake_water .and. lake_total_depth <= smallnumber) then
						if (DEF_METHANE%lake_zero_depth_fatal) then
							write(6,*) 'ERROR: lake methane enabled but current dz_lake/lakedepth are non-positive.'
							CALL CoLM_stop ()
						else if (.not. warned_lake_zero_depth) then
							write(6,*) 'WARNING: lake methane enabled with zero water depth; exchange is disabled.'
							warned_lake_zero_depth = .true.
                  end if
               end if
					! while the lake is stably stratified (bottom water denser than
					! surface water by more than lake_strat_drho, density as in MOD_Lake),
					! the water that meets the sediment does not reach the air; the CH4
					! it gathers leaves at overturn. The CH4 node stands for that water.
					lake_stratified = .false.
					if (is_lake_water .and. DEF_METHANE%lake_strat_drho > 0._r8) then
						lake_jbot = 1
						do jj_lake = 1, nl_lake
							if (dz_lake(jj_lake) > smallnumber) lake_jbot = jj_lake
						end do
						lake_rho_top = 1000._r8 * (1._r8 - 1.9549e-5_r8 * abs(t_lake(1) - 277._r8)**1.68_r8)
						lake_rho_bot = 1000._r8 * (1._r8 - 1.9549e-5_r8 * abs(t_lake(lake_jbot) - 277._r8)**1.68_r8)
						lake_stratified = (lake_rho_bot - lake_rho_top) > DEF_METHANE%lake_strat_drho
					end if
					if (is_lake_water .and. lake_liquid_depth > smallnumber .and. &
					    t_h2osfc >= tfrz .and. lake_icefrac(1) <= 0.1_r8 .and. .not. lake_stratified) then
						t_soisno_c = t_h2osfc - tfrz
						if (ieee_is_nan(forc_us) .or. ieee_is_nan(forc_vs) .or. &
						    abs(forc_us) > 1.e30_r8 .or. abs(forc_vs) > 1.e30_r8) then
							lake_wind = 0._r8
						else
							lake_wind = sqrt(max(forc_us*forc_us + forc_vs*forc_vs, 0._r8))
						endif
						lake_k600_cmhr = 2.07_r8 + 0.215_r8 * max(lake_wind, 0._r8)**1.7_r8
						t_soisno_c = max(0._r8, min(30._r8, t_soisno_c))
						lake_schmidt_ch4 = max(300._r8, s_con(s,1) + s_con(s,2)*t_soisno_c + &
							s_con(s,3)*t_soisno_c**2 + s_con(s,4)*t_soisno_c**3)
						lake_exchange_vel = lake_k600_cmhr * 1.e-2_r8 / 3600._r8 * &
							(lake_schmidt_ch4 / 600._r8)**(-2._r8/3._r8)
					endif

					spec_grnd_cond(s) = 1._r8/(1._r8/max(grnd_methane_cond_base, smallnumber) + &
						snow_resis + pond_resis)
					if (is_lake_water) spec_grnd_cond(s) = k_h_cc(0,s) * lake_exchange_vel
				end if
			end do ! j

			! Determine gas diffusion and fraction of open pore (f_a)
			do j = 1,nl_soil
				t_soisno_c = t_soisno(j) - tfrz

				if (porsl(j) <= smallnumber) then
					diffus(j) = smallnumber
				else if (j <= jwt) then  ! Above the WT
					if (organic_max > 0._r8) then
						om_frac = min(DEF_METHANE%om_frac_sf*cellorg(j)/organic_max, 1._r8)
						! Use first power, not square as in iniTimeConst
					else
						om_frac = 1._r8
					end if
					diffus (j) = (d_con_g(s,1) + d_con_g(s,2)*t_soisno_c) * 1.e-4_r8 * &
								(om_frac * vol_gas(j)**(10._r8/3._r8) / porsl(j)**2._r8 + &
								(1._r8-om_frac) * vol_gas(j)**2._r8 * (vol_gas(j)/porsl(j))**(3._r8 / max(bsw(j), 0.5_r8)) ) &
								* DEF_METHANE%scale_factor_gasdiff
					else ! Below the WT use liquid-water diffusivity; ice-filled pore space is excluded via vol_aqu.
						diffus (j) = max(vol_aqu(j), smallnumber)**DEF_METHANE%satpow * &
							(d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
							* DEF_METHANE%scale_factor_liqdiff
					end if ! Above/below the WT
					if (patchtype == 4 .and. DEF_METHANE%allowlakeprod .and. j > jwt) then
						if (s == 1) diffus(j) = diffus(j) * DEF_METHANE%lake_liqdiff_scale
						if (s == 2) diffus(j) = diffus(j) * DEF_METHANE%lake_o2_liqdiff_scale
					end if
					diffus(j) = max(diffus(j), smallnumber) ! Prevent overflow
				enddo ! j

			lake_sed_water_cond = 0._r8
			if (is_lake_water .and. lake_liquid_depth > smallnumber) then
				! Half of the top sediment layer is the only unresolved distance
				! between its centre and the well-mixed overlying water.
				lake_sed_water_cond = 2._r8 * diffus(1) / max(dz_soisno(1), smallnumber)
			endif

			! Zero conductance arrays so end-rows that intentionally skip one
			! direction stay at a defined value under -finit-real=nan builds.
			dp1_zp1(:) = 0._r8
			dm1_zm1(:) = 0._r8

			do j = 1,nl_soil

				! Set up coefficients for tridiagonal solver.
				if (is_lake_water .and. j == 1) then
					dm1_zm1(j) = lake_sed_water_cond
					if (j < nl_soil) &
						dp1_zp1(j) = 2._r8/(dz_soisno(j)/diffus(j)+dz_soisno(j+1)/diffus(j+1))
				else if (j == 1 .and. j /= jwt .and. j /= jwt+1) then
					dm1_zm1(j) = 1._r8/(1._r8/spec_grnd_cond(s)+dz_soisno(j)/(diffus(j)*2._r8))
					! replace Diffusivity / Delta_z by conductance (grnd_methane_cond) for top layer
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
			end do ! j; nl_soil

			lake_water_stock_old = 0._r8
			if (is_lake_water) then
				if (s == 1) then
					lake_water_stock_old = lake_water_ch4_stock
				else
					lake_water_stock_old = lake_water_o2_stock
				endif
			endif

			! Perform a second loop for the tridiagonal coefficients since need dp1_zp1 and dm1_z1 at each depth
				do j = 0,nl_soil
					conc_ch4_rel_old(j) = conc_ch4_rel(j)

					if (j > 0) dzj = dz_soisno(j)
					if (j == 0) then
						at(j) = 0._r8
						if (is_lake_water) then
							bt(j) = lake_water_storage_depth / deltim + &
								lake_sed_water_cond + lake_exchange_vel
							ct(j) = -lake_sed_water_cond
							rt(j) = lake_water_stock_old / deltim + &
								lake_exchange_vel * k_h_cc(0,s) * c_atm(s)
						else
							bt(j) = 1._r8
							ct(j) = 0._r8
							rt(j) = c_atm(s) ! fixed atmospheric boundary
						endif
					elseif (is_lake_water) then
						! Lake CH4 and O2 both use a backward-Euler M-matrix.
						! The top sediment row shares Gsw with the water row above,
						! so their internal exchange cancels exactly in the budget.
						at(j) = -dm1_zm1(j) / dzj
						if (j < nl_soil) then
							bt(j) = epsilon_t(j,s) / deltim + &
								(dp1_zp1(j) + dm1_zm1(j)) / dzj
							ct(j) = -dp1_zp1(j) / dzj
						else
							bt(j) = epsilon_t(j,s) / deltim + dm1_zm1(j) / dzj
							ct(j) = 0._r8
						endif
					elseif (backward_euler_transport) then
						! CH4 and O2 use a fully implicit backward-Euler diffusion
						! matrix.  The old Crank-Nicolson half-explicit
						! diffusion term could overshoot sharp CH4 gradients in
						! cold / low-storage layers and then rely on the
						! negative-concentration clip branch.  With the
						! competition limiter above, the CH4 RHS is non-negative;
						! this M-matrix transport step is therefore
						! positivity-preserving while keeping the same
						! Henry-law interface coefficients.
						if (j < nl_soil .and. j == jwt) then
							at(j) = -1._r8 / dzj * dm1_zm1(j)
							bt(j) = epsilon_t(j,s) / deltim + 1._r8 / dzj * &
								(dp1_zp1(j)*k_h_cc(j,s) + dm1_zm1(j))
							ct(j) = -1._r8 / dzj * dp1_zp1(j)
						elseif (j < nl_soil .and. j == jwt+1) then
							at(j) = -1._r8 / dzj * dm1_zm1(j) * k_h_cc(j-1,s)
							bt(j) = epsilon_t(j,s) / deltim + 1._r8 / dzj * &
								(dp1_zp1(j) + dm1_zm1(j))
							ct(j) = -1._r8 / dzj * dp1_zp1(j)
						elseif (j < nl_soil) then
							at(j) = -1._r8 / dzj * dm1_zm1(j)
							bt(j) = epsilon_t(j,s) / deltim + 1._r8 / dzj * &
								(dp1_zp1(j) + dm1_zm1(j))
							ct(j) = -1._r8 / dzj * dp1_zp1(j)
						else if (j == nl_soil .and. j== jwt+1) then
							at(j) = -1._r8 / dzj * dm1_zm1(j) * k_h_cc(j-1,s)
							bt(j) = epsilon_t(j,s) / deltim + 1._r8 / dzj * dm1_zm1(j)
							ct(j) = 0._r8
						else
							at(j) = -1._r8 / dzj * dm1_zm1(j)
							bt(j) = epsilon_t(j,s) / deltim + 1._r8 / dzj * dm1_zm1(j)
							ct(j) = 0._r8
						endif
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

						! Source and epsilon_t are constant over the implicit step,
						! so a single evaluation is sufficient (formerly a CN-style
						! average against a same-step "old" copy, which was a no-op).
						dzj = dz_soisno(j)
						if (backward_euler_transport) then
							! Backward-Euler CH4 diffusion has no explicit
							! diffusion contribution on the RHS.  With
							! source(j,1) limited so source + storage/dt is
							! non-negative, the implicit matrix above prevents
							! large negative post-solve concentrations.
							rt(j) = epsilon_t(j,s) / deltim * conc_ch4_rel(j) + source(j,s)
						elseif (j < nl_soil .and. j == jwt) then ! concentration inside needs to be mult. by k_h_cc for dp1_zp1 term
							rt(j) = epsilon_t(j,s) / deltim * conc_ch4_rel(j) +           &
								0.5_r8 / dzj * (dp1_zp1(j) * (conc_ch4_rel(j+1)-conc_ch4_rel(j)*k_h_cc(j,s)) - &
								dm1_zm1(j) * (conc_ch4_rel(j)  -conc_ch4_rel(j-1))) + &
							source(j,s)
					elseif (j < nl_soil .and. j == jwt+1) then
						! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
						rt(j) = epsilon_t(j,s) / deltim * conc_ch4_rel(j) +           &
							0.5_r8 / dzj * (dp1_zp1(j) * (conc_ch4_rel(j+1)-conc_ch4_rel(j)) - &
							dm1_zm1(j) * (conc_ch4_rel(j) -conc_ch4_rel(j-1)*k_h_cc(j-1,s))) + &
							source(j,s)
					elseif (j < nl_soil) then
						rt(j) = epsilon_t(j,s) / deltim * conc_ch4_rel(j) +           &
							0.5_r8 / dzj * (dp1_zp1(j) * (conc_ch4_rel(j+1)-conc_ch4_rel(j)) - &
							dm1_zm1(j) * (conc_ch4_rel(j)  -conc_ch4_rel(j-1))) + &
							source(j,s)
					else if (j == nl_soil .and. j== jwt+1) then
						! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
						rt(j) = epsilon_t(j,s) / deltim * conc_ch4_rel(j) +           &
							0.5_r8 / dzj * ( - dm1_zm1(j) * (conc_ch4_rel(j) -conc_ch4_rel(j-1)*k_h_cc(j-1,s))) + &
							source(j,s)
					else  !j==nl_soil
						rt(j) = epsilon_t(j,s) / deltim * conc_ch4_rel(j) +           &
							0.5_r8 / dzj * ( - dm1_zm1(j) * (conc_ch4_rel(j)  -conc_ch4_rel(j-1))) + &
							source(j,s)
					endif

				enddo ! j; nl_soil

					! bubble CH4 dissolved on the way up enters the lake water node
					if (lake_ebul_dissol > 0._r8) rt(0) = rt(0) + lake_ebul_dissol
					! the dissolved share of the under-ice bubble release
					if (lake_icebub_dissol > 0._r8) rt(0) = rt(0) + lake_icebub_dissol

				call Tridiagonal(0, nl_soil, &
					jtop, &
					at(:), &
					bt(:), &
					ct(:), &
					rt(:), &
					conc_ch4_rel(0:nl_soil))



					! Calculate net methane flux to the atmosphere from the surface (+ to atm).
					if (is_lake_water) then
						lake_water_ch4_stock = lake_water_storage_depth * conc_ch4_rel(0)
						methane_surf_diff = lake_exchange_vel * &
							(conc_ch4_rel(0) - k_h_cc(0,s) * c_atm(s))
						lake_sed_ch4_flux = lake_sed_water_cond * &
							(conc_ch4_rel(1) - conc_ch4_rel(0))
					elseif (jwt /= 0) then ! WT not at the surface
						if (backward_euler_transport) then
							methane_surf_diff = dm1_zm1(1) * (conc_ch4_rel(1) - c_atm(s)) ! [mol/m2/s]
						else
							methane_surf_diff = dm1_zm1(1) * ( (conc_ch4_rel(1)+conc_ch4_rel_old(1))/2._r8 &
								- c_atm(s)) ! [mol/m2/s]
						endif
					else ! WT at the surface; i.e., jwt==0
						if (backward_euler_transport) then
							methane_surf_diff = dm1_zm1(1) * (conc_ch4_rel(1) &
								- c_atm(s)*k_h_cc(0,s)) ! [mol/m2/s]
						else
							methane_surf_diff = dm1_zm1(1) * ( (conc_ch4_rel(1)+conc_ch4_rel_old(1))/2._r8 &
								- c_atm(s)*k_h_cc(0,s)) ! [mol/m2/s]
						endif
						! atmospheric concentration gets mult. by k_h_cc as above
					endif
				if (jwt /= 0) then
					methane_surf_ebul = 0._r8
				else
					methane_surf_ebul = methane_ebul_tot ! [mol/m2/s]
					! less the part dissolved into the lake water node
					if (lake_ebul_dissol > 0._r8) methane_surf_ebul = methane_ebul_tot - lake_ebul_dissol
					! held under the ice, or plus the release not dissolved
					if (lake_icebub_trapping) then
						methane_surf_ebul = 0._r8
					elseif (lake_icebub_release > 0._r8) then
						methane_surf_ebul = methane_surf_ebul + (lake_icebub_release - lake_icebub_dissol)
					endif
				endif

				! Snapshot the raw diffusive flux BEFORE the negative-concentration
				! clip-credit and the column-closure residual fold-in below, so
				! history can expose the pure-physics value separately from the
				! numerical correction (methane_balance_residual already captures
				! the residual; the clip credit is its own contribution).
				methane_surf_diff_phys = methane_surf_diff

				! Ensure that concentrations stay above 0
				! This should be done after the flux, so that the flux calculation is consistent.
				! Original strict 1e-12 abort tripped on implicit-solver roundoff at
				! cold start; relaxing to 1e-6 still aborted because cold-start sink
				! pull on initially-empty CH4 columns drives deficits on the order of
				! the per-timestep flux scale.
				! The clip+credit branch conserves column mass.  Diagnose small clips
				! and let production configurations reject material corrections.
					ch4_clip_deficit_col = 0._r8
					do j = 1,nl_soil
						if (conc_ch4_rel(j) < 0._r8) then
							deficit = - conc_ch4_rel(j)*epsilon_t(j,1)*dz_soisno(j)  ! Mol/m^2 added
							ch4_clip_deficit_col = ch4_clip_deficit_col + deficit
							! Always clip and credit; column mass is conserved by the surface flux update.
							methane_surf_diff = methane_surf_diff - deficit / deltim
							methane_ch4_clip_credit = methane_ch4_clip_credit - deficit / deltim
							conc_ch4_rel(j) = 0._r8
					end if
				enddo
					if (DEF_METHANE%numerical_correction_fatal_threshold > 0._r8 .and. &
					    ch4_clip_deficit_col > DEF_METHANE%numerical_correction_fatal_threshold) then
						write(6,*) 'ERROR: CH4 nonnegative column clip exceeds fatal threshold.'
						write(6,*) 'Lat,Lon, column deficit, threshold = ', dlat, dlon, &
							ch4_clip_deficit_col, DEF_METHANE%numerical_correction_fatal_threshold
						CALL CoLM_stop ('CH4 nonnegative column clip exceeds fatal threshold')
					endif
					if (ch4_clip_deficit_col > 1.e-3_r8 .and. .not. warned_ch4_clip) then
						write(6,*)'WARNING: methane_tran implicit column clip > 1e-3 mol/m2 timestep.'
						write(6,*)'Lat,Lon, column deficit=',dlat,dlon,ch4_clip_deficit_col
						warned_ch4_clip = .true.
					endif


			elseif (s == 2) then  ! O2

				! Set rt, since it depends on conc
				do j = 1,nl_soil

					! Source and epsilon_t are constant over the implicit step
					! (see the CH4 branch above for the same simplification).
					dzj     = dz_soisno(j)
					if (backward_euler_transport) then
						! The competition limiter makes storage/dt + source
						! nonnegative.  With the M-matrix assembled above this
						! prevents the O2 overshoot that previously required a
						! mass-creating post-solve floor.
						rt(j) = epsilon_t(j,s) / deltim * conc_o2_rel(j) + source(j,s)
					elseif (is_lake_water) then
						rt(j) = epsilon_t(j,s) / deltim * conc_o2_rel(j) + source(j,s)
					elseif (j < nl_soil .and. j == jwt) then ! concentration inside needs to be mult. by k_h_cc for dp1_zp1 term
						rt(j) = epsilon_t(j,s) / deltim * conc_o2_rel(j) +           &
							0.5_r8 / dzj * (dp1_zp1(j) * (conc_o2_rel(j+1)-conc_o2_rel(j)*k_h_cc(j,s)) - &
							dm1_zm1(j) * (conc_o2_rel(j)  -conc_o2_rel(j-1))) + &
							source(j,s)
					elseif (j < nl_soil .and. j == jwt+1) then
						! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
						rt(j) = epsilon_t(j,s) / deltim * conc_o2_rel(j) +           &
							0.5_r8 / dzj * (dp1_zp1(j) * (conc_o2_rel(j+1)-conc_o2_rel(j)) - &
							dm1_zm1(j) * (conc_o2_rel(j) -conc_o2_rel(j-1)*k_h_cc(j-1,s))) + &
							source(j,s)
					elseif (j < nl_soil) then
						rt(j) = epsilon_t(j,s) / deltim * conc_o2_rel(j) +           &
							0.5_r8 / dzj * (dp1_zp1(j) * (conc_o2_rel(j+1)-conc_o2_rel(j)) - &
							dm1_zm1(j) * (conc_o2_rel(j)  -conc_o2_rel(j-1))) + &
							source(j,s)
					else if (j == nl_soil .and. j== jwt+1) then
						! concentration above needs to be mult. by k_h_cc for dm1_zm1 term
						rt(j) = epsilon_t(j,s) / deltim * conc_o2_rel(j) +           &
							0.5_r8 / dzj * ( - dm1_zm1(j) * (conc_o2_rel(j) -conc_o2_rel(j-1)*k_h_cc(j-1,s))) + &
							source(j,s)
					else  !j==nl_soil
						rt(j) = epsilon_t(j,s) / deltim * conc_o2_rel(j) +           &
							0.5_r8 / dzj * ( - dm1_zm1(j) * (conc_o2_rel(j)  -conc_o2_rel(j-1))) + &
							source(j,s)
					endif
				enddo ! j; nl_soil

				call Tridiagonal(0, nl_soil, jtop, &
					at(:), &
					bt(:), &
					ct(:), &
					rt(:), &
					conc_o2_rel(0:nl_soil))

					if (is_lake_water) then
						lake_water_o2_stock = lake_water_storage_depth * conc_o2_rel(0)
						lake_sed_o2_flux = lake_sed_water_cond * &
							(conc_o2_rel(1) - conc_o2_rel(0))
						! the bed's O2 demand is O2 the water gives to the sediment
						if (lake_water_sod > 0._r8) lake_sed_o2_flux = lake_sed_o2_flux - lake_water_sod
						lake_air_o2_flux = lake_exchange_vel * &
							(conc_o2_rel(0) - k_h_cc(0,s) * c_atm(s))
					endif

					! Allow physical O2 supersaturation.  Backward Euler should
					! preserve non-negativity; a material negative value is an
					! invariant violation, not a quantity to hide with a clip.
					o2_cap_loss_col = 0._r8
					o2_cap_gain_col = 0._r8
					do j = 1,nl_soil

						o2_before_cap = conc_o2_rel(j)
						! No per-layer abort on a negative ratio.  conc_o2_rel is
						! conc_o2 / epsilon_t, and epsilon_t falls to the liquid content
						! of the layer, so the division amplifies solver round-off by
						! 1 / epsilon_t.  That factor varies continuously with how frozen
						! the layer is - 8.3e3 at vol_aqu = 1.2e-4, 1e12 once vol_aqu hits
						! the floor - so no fixed ratio threshold can separate round-off
						! from divergence.  Observed artefacts reached a ratio of -1902 at
						! a stock of -1.9e-9 mol m-3, and -4.1e-3 at -4.9e-7 mol m-3; both
						! are physically zero.
						! The column-integrated guard below is the same mechanism CH4 uses
						! (ch4_clip_deficit_col): it accumulates the correction as a stock
						! in mol m-2, which is immune to the amplification, and aborts on
						! DEF_METHANE%numerical_correction_fatal_threshold.  Round-off
						! artefacts contribute ~5e-8 mol m-2 against a 1e-3 threshold,
						! while a real divergence saturates it at once.
						conc_o2_rel(j) = max (conc_o2_rel(j), 0._r8)
						if (conc_o2_rel(j) > o2_before_cap) then
							o2_cap_gain_col = o2_cap_gain_col + &
								(conc_o2_rel(j) - o2_before_cap) * epsilon_t(j,2) * dz_soisno(j)
						endif
						enddo
					if (deltim > 0._r8) then
						o2_cap_loss = o2_cap_loss_col / deltim
						o2_cap_gain = o2_cap_gain_col / deltim
					endif
					if (DEF_METHANE%numerical_correction_fatal_threshold > 0._r8 .and. &
					    o2_cap_gain_col > DEF_METHANE%numerical_correction_fatal_threshold) then
						write(6,*) 'ERROR: O2 nonnegative floor exceeds fatal threshold.'
						write(6,*) 'Lat,Lon, correction, threshold = ', dlat, dlon, o2_cap_gain_col, &
							DEF_METHANE%numerical_correction_fatal_threshold
						CALL CoLM_stop ('O2 nonnegative floor exceeds fatal threshold')
					endif
				! Retain the diagnostic for roundoff-scale corrections only.
					if (o2_cap_gain_col > 1.e-12_r8 &
					    .and. .not. warned_o2_cap) then
						write(6,*)'WARNING: O2 roundoff adjustment > 1e-12 mol/m2 timestep.'
						write(6,*)'Lat,Lon, cap_loss, cap_gain=',dlat,dlon,o2_cap_loss_col,o2_cap_gain_col
						warned_o2_cap = .true.
					endif

		   endif  ! species

		enddo  ! species

		! Update absolute concentrations per unit volume
		do j = 1,nl_soil ! No need to update the atm. level concentrations
			conc_methane(j) = conc_ch4_rel(j)*epsilon_t(j,1)
			conc_o2(j)  = conc_o2_rel(j) *epsilon_t(j,2)
		end do

		! Do Balance Check and absorb small
		!    discrepancy into surface flux.
		err_methane = 0._r8
		do j = 1,nl_soil
			err_methane = err_methane + (conc_methane(j) - conc_ch4_bef(j))*dz_soisno(j)
			err_methane = err_methane - methane_prod_depth(j)*dz_soisno(j)*deltim
			err_methane = err_methane + methane_oxid_depth(j)*dz_soisno(j)*deltim
			err_methane = err_methane + methane_tran_depth(j)*dz_soisno(j)*deltim
		end do
		if (is_lake_water) then
			err_methane = err_methane + lake_water_ch4_stock - lake_ch4_stock_bef + &
				lake_water_ch4_oxid * deltim
			! change of the under-ice bubble store
			if (lake_icebub_stock /= lake_icebub_stock_bef) &
				err_methane = err_methane + (lake_icebub_stock - lake_icebub_stock_bef)
		endif

		! History output uses the effective CH4 conductance, including snow and pond resistance.
		grnd_methane_cond_effective = spec_grnd_cond(1)

		err_methane = err_methane + (methane_surf_aere + methane_surf_ebul + methane_surf_diff)*deltim


		! Always clip-credit the closure error into methane_surf_diff so the
		! column inventory stays conserved.  Production runs reject material
		! residuals through the same single-timestep threshold.
		methane_balance_residual = -err_methane/deltim
		methane_surf_diff = methane_surf_diff + methane_balance_residual
		if (DEF_METHANE%numerical_correction_fatal_threshold > 0._r8 .and. &
		    abs(err_methane) > DEF_METHANE%numerical_correction_fatal_threshold) then
			write(6,*) 'ERROR: CH4 transport closure residual exceeds fatal threshold.'
			write(6,*) 'Date, Lat, Lon, residual, threshold = ', idate, dlat, dlon, err_methane, &
				DEF_METHANE%numerical_correction_fatal_threshold
			CALL CoLM_stop ('CH4 transport closure residual exceeds fatal threshold')
		endif
		if (abs(err_methane) > 1.e-3_r8 .and. .not. warned_ch4_diffusion_residual) then
			write(6,*)'WARNING: CH4 diffusion residual > 1e-3 mol/m^2/timestep, istep, err=', idate, err_methane
			write(6,*)'Lat,Lon=',dlat,dlon
			warned_ch4_diffusion_residual = .true.
		end if

			end subroutine methane_tran

			!---------------------------------------------------------------------------
			real(r8) function methane_wetland_water_depth(patchtype, wdsrf, wetwat)
				! Use one surface-water-depth definition for all inundation schemes.
				! In the non-dynamic wetland hydrology branch, wdsrf/wa/wetwat are
				! merged into wetwat, then overflow is stored back in wdsrf. Those two
				! post-step reservoirs are disjoint and both contribute to ponded
				! depth. Dynamic wetland tiles follow the soil-ground branch, so wetwat
				! is not part of the active ponded-water state there.
				integer, intent(in) :: patchtype
				real(r8), intent(in) :: wdsrf, wetwat

				if (patchtype == 2 .and. .not. DEF_USE_Dynamic_Wetland) then
					methane_wetland_water_depth = max(wdsrf, 0._r8) + max(wetwat, 0._r8)
				else
					methane_wetland_water_depth = max(wdsrf, 0._r8)
				endif
			end function methane_wetland_water_depth

			!---------------------------------------------------------------------------
   real(r8) function wetland_max_wtd_of(ipatch)
      ! Water-table floor [m] of a dynamic wetland tile: one global value, or
      ! one per wetland class when use_biome_wetland_max_wtd is set. The class
      ! is published by the driver after the host hydrology of the same step,
      ! so hydrology reads the previous step's class; it follows latitude and
      ! the modelled top-layer carbon, and the zero it starts from falls back
      ! to the global value. With wetland_dim_class the per-class
      ! values go by the tile's shares instead of its class: the temperate
      ! marsh floor on the emergent marsh share s_e, the tropical peat floor
      ! on the tropical peat dome share s_d and the global floor on the rest,
      ! s_e z_marsh + s_d z_peat + (1 - s_e - s_d) z. With
      ! wetland_max_wtd_bog_share the open-bog share of the
      ! part taking lateral inflow, w_b = wetland_bog_share /
      ! wetland_inflow_share, holds its floor at that depth, and the floor is
      ! the area mean (1 - w_b) z_class + w_b z_bog.
      integer, intent(in) :: ipatch
      real(r8) :: fbog, wb, se, sd

      wetland_max_wtd_of = DEF_METHANE%wetland_max_wtd
      if (DEF_METHANE%use_biome_wetland_max_wtd .and. DEF_METHANE%wetland_dim_class) then
         if (allocated(wetveg_emerg) .and. allocated(wetveg_dome)) then
            if (ipatch >= 1 .and. ipatch <= min(size(wetveg_emerg), size(wetveg_dome))) then
               se = min(max(wetveg_emerg(ipatch), 0._r8), 1._r8)
               sd = min(max(wetveg_dome(ipatch), 0._r8), 1._r8 - se)
               wetland_max_wtd_of = se * DEF_METHANE%wetland_max_wtd_temperate_marsh &
                  + sd * DEF_METHANE%wetland_max_wtd_tropical_peat &
                  + (1._r8 - se - sd) * DEF_METHANE%wetland_max_wtd
            endif
         endif
      elseif (DEF_METHANE%use_biome_wetland_max_wtd .and. allocated(methane_wetland_type)) then
         if (ipatch >= 1 .and. ipatch <= size(methane_wetland_type)) then
      select case (nint(methane_wetland_type(ipatch)))
      case (BIOME_TROPICAL_PEAT)
         wetland_max_wtd_of = DEF_METHANE%wetland_max_wtd_tropical_peat
      case (BIOME_TROPICAL_FLOODPLAIN)
         wetland_max_wtd_of = DEF_METHANE%wetland_max_wtd_tropical_floodplain
      case (BIOME_TEMPERATE_MARSH)
         wetland_max_wtd_of = DEF_METHANE%wetland_max_wtd_temperate_marsh
      case (BIOME_BOREAL_FEN)
         wetland_max_wtd_of = DEF_METHANE%wetland_max_wtd_boreal_fen
      case (BIOME_BOREAL_BOG)
         wetland_max_wtd_of = DEF_METHANE%wetland_max_wtd_boreal_bog
      end select
         endif
      endif

      ! A class without a floor keeps none.
      if (wetland_max_wtd_of < 0._r8) return
      fbog = wetland_bog_share(ipatch)
      if (fbog > 0._r8) then
         wb = fbog / wetland_inflow_share(ipatch)
         wetland_max_wtd_of = (1._r8 - wb) * wetland_max_wtd_of &
            + wb * DEF_METHANE%wetland_max_wtd_bog_share
      endif
   end function wetland_max_wtd_of

   !---------------------------------------------------------------------------
   real(r8) function wetland_inflow_share(ipatch)
      ! share of a wetland tile that takes lateral inflow when its water
      ! table falls below the floor. A rain-fed (ombrotrophic) bog is fed by
      ! precipitation alone, so only the minerotrophic rest is held up;
      ! 1 when no rain-fed share is known. The rain-fed share is the
      ! rainfed_share of wetland_veg_file, which v2/scripts/mk_rainfed_share.py
      ! fills with the permafrost-bog (PEB) share of the BAWLD wetland classes
      ! (Olefeldt et al. 2021), open bogs keeping the floor,
      ! or wetland_ombro_share_site at a tower. With wetland_dim_class
      ! it is the permafrost plateau part s_r - s_b of the water-source
      ! dimension's rain-fed share, from wetland_share_rain_site and
      ! wetland_share_bog_site at a tower. The host scales the refill by
      ! it, or with DEF_WETLAND_INFLOW_AREA_SPLIT splits the tile by it.
      integer, intent(in) :: ipatch

      wetland_inflow_share = 1._r8
      if (.not. allocated(wetveg_ombro)) return
      if (ipatch < 1 .or. ipatch > size(wetveg_ombro)) return
      wetland_inflow_share = 1._r8 - min(max(wetveg_ombro(ipatch), 0._r8), 1._r8)
   end function wetland_inflow_share

   !---------------------------------------------------------------------------
   real(r8) function wetland_bog_share(ipatch)
      ! open-bog share of a wetland tile, the bog_share of
      ! wetland_veg_file, which v2/scripts/mk_bog_share.py fills with the bog
      ! (BOG) share of the BAWLD wetland classes (Olefeldt et al. 2021), or
      ! wetland_bog_share_site at a tower (wetland_share_bog_site with
      ! wetland_dim_class). Bogs sit within the part that
      ! takes lateral inflow (the rain-fed rest is permafrost bog), so
      ! the share is capped at wetland_inflow_share. 0 unless
      ! wetland_max_wtd_bog_share is set, and without a bog share.
      integer, intent(in) :: ipatch

      wetland_bog_share = 0._r8
      if (DEF_METHANE%wetland_max_wtd_bog_share < 0._r8) return
      if (.not. allocated(wetveg_bog)) return
      if (ipatch < 1 .or. ipatch > size(wetveg_bog)) return
      wetland_bog_share = min(max(wetveg_bog(ipatch), 0._r8), wetland_inflow_share(ipatch))
   end function wetland_bog_share

   !---------------------------------------------------------------------------
   real(r8) function wetland_peat_outflow(nl_soil, zi, hksat, om_density, zwt)
      ! Lateral outflow [mm s-1] of a dynamic wetland through its saturated
      ! zone (PEAT-CLSM, Bechtold et al. 2019, eqs 6-9): Q = c T, with T the
      ! transmissivity below the water table. In each layer the peat part
      ! integrates K0 (1 + 100 z)^(-m) over the saturated depth and the
      ! mineral part takes the layer's saturated conductivity; the organic
      ! fraction OM_density / organic_max weighs the two. A ponded tile
      ! drains at T(0): the macro-scale conductivity includes the flow
      ! through hollows over the surface.
      integer,  intent(in) :: nl_soil
      real(r8), intent(in) :: zi(0:nl_soil)          ! layer interfaces [m]
      real(r8), intent(in) :: hksat(1:nl_soil)       ! saturated conductivity [mm s-1]
      real(r8), intent(in) :: om_density(1:nl_soil)  ! organic matter density [kg m-3]
      real(r8), intent(in) :: zwt                    ! water-table depth [m]

      integer  :: j
      real(r8) :: ztop, zbot, fom, tpeat, trans, a, e

      wetland_peat_outflow = 0._r8
      if (.not. DEF_METHANE%wetland_peat_drainage) return

      e = 1._r8 - DEF_METHANE%peat_m
      a = DEF_METHANE%peat_K0 / (100._r8 * (DEF_METHANE%peat_m - 1._r8))
      trans = 0._r8
      do j = 1, nl_soil
         ztop = max(zi(j-1), max(zwt, 0._r8))
         zbot = zi(j)
         if (zbot <= ztop) cycle
         fom = 0._r8
         if (om_density(j) > 0._r8 .and. om_density(j) < 1.e29_r8) &
            fom = min(om_density(j) / DEF_METHANE%organic_max, 1._r8)
         tpeat = a * ((1._r8 + 100._r8 * ztop)**e - (1._r8 + 100._r8 * zbot)**e)
         trans = trans + fom * tpeat + (1._r8 - fom) * max(hksat(j), 0._r8) * 1.e-3_r8 * (zbot - ztop)
      enddo
      wetland_peat_outflow = DEF_METHANE%peat_c * trans * 1.e3_r8
   end function wetland_peat_outflow

			!---------------------------------------------------------------------------
   logical function wetland_is_temperate_marsh(ipatch)
      ! Whether a patch's wetland tile is of the temperate marsh class, the
      ! emergent marsh whose litter and standing dead lie on its ponded water
      ! (Goulden et al. 2007); false before the class is first published.
      integer, intent(in) :: ipatch

      wetland_is_temperate_marsh = .false.
      if (.not. allocated(methane_wetland_type)) return
      if (ipatch < 1 .or. ipatch > size(methane_wetland_type)) return
      wetland_is_temperate_marsh = nint(methane_wetland_type(ipatch)) == BIOME_TEMPERATE_MARSH
   end function wetland_is_temperate_marsh

   !---------------------------------------------------------------------------
   real(r8) function wetland_emergent_share(ipatch)
      ! Share of a patch's wetland tile under emergent marsh, whose litter
      ! and standing dead lie on its ponded water (Goulden et al. 2007), for
      ! the host's pond litter resistance. With wetland_dim_class the
      ! tile's emergent marsh share (GLWD 17, or wetland_share_emerg_site at
      ! a tower); otherwise 1 for the temperate marsh class and 0 for any
      ! other class and before the class is first published.
      integer, intent(in) :: ipatch

      wetland_emergent_share = 0._r8
      if (DEF_METHANE%wetland_dim_class) then
         if (.not. allocated(wetveg_emerg)) return
         if (ipatch < 1 .or. ipatch > size(wetveg_emerg)) return
         wetland_emergent_share = min(max(wetveg_emerg(ipatch), 0._r8), 1._r8)
      elseif (wetland_is_temperate_marsh(ipatch)) then
         wetland_emergent_share = 1._r8
      endif
   end function wetland_emergent_share

   !---------------------------------------------------------------------------
   real(r8) function wetland_moss_rss(ipatch, zwt, wdsrf)
      ! Floor [s m-1] of the snow-free soil surface resistance of a dynamic
      ! wetland tile: the evaporation resistance of its moss and peat
      ! surface, 100 s m-1 while the water table is within 0.1 m of the
      ! surface, rising below that with an e-folding depth of 0.155 m (fitted
      ! to the bulk resistance by water-table class of Alekseychik et al.
      ! 2018) up to 1000 s m-1, times the tile's peatland share. Tropical
      ! classes, whose peat carries no moss carpet, and ponded tiles, which
      ! evaporate as open water, take none; so does a tile before its class
      ! is first published. With wetland_dim_class the weight is the
      ! tile's moss-carpeted peatland share (GLWD 22-25, or
      ! wetland_share_moss_site at a tower) in place of the class test and
      ! the peatland share.
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: zwt     ! water-table depth [m]
      real(r8), intent(in) :: wdsrf   ! ponded water depth [mm]

      real(r8), parameter :: r0 = 100._r8, z0 = 0.10_r8, zs = 0.155_r8, rmax = 1000._r8
      real(r8) :: smoss

      wetland_moss_rss = 0._r8
      if (.not. DEF_METHANE%wetland_moss_rss) return
      if (DEF_METHANE%wetland_dim_class) then
         if (.not. allocated(wetveg_moss)) return
         if (ipatch < 1 .or. ipatch > size(wetveg_moss)) return
         smoss = min(max(wetveg_moss(ipatch), 0._r8), 1._r8)
      else
         if (.not. allocated(wetland_tile_peat_share)) return
         if (.not. allocated(methane_wetland_type)) return
         if (ipatch < 1 .or. ipatch > size(wetland_tile_peat_share)) return
         select case (nint(methane_wetland_type(ipatch)))
         case (BIOME_TEMPERATE_MARSH, BIOME_BOREAL_FEN, BIOME_BOREAL_BOG)
            smoss = wetland_tile_peat_share(ipatch)
         case default
            return
         end select
      endif

      if (wdsrf > 1._r8) return
      wetland_moss_rss = r0
      if (zwt > z0) wetland_moss_rss = min(rmax, r0 * exp((zwt - z0) / zs))
      wetland_moss_rss = smoss * wetland_moss_rss
   end function wetland_moss_rss

   !---------------------------------------------------------------------------
   subroutine methane_tsoi_annual(ipatch, deltim, z_soi, dz_soi, t_soi, year_end)
      ! Annual mean soil temperature of a wetland tile over the top metre, the
      ! reference of its CH4 production temperature factor: the
      ! thickness-weighted mean of the layers whose node lies within 1 m is
      ! integrated in time, and every full 365 days the mean replaces
      ! annavg_tsoi and the integral starts again. With
      ! q10methane_base_unfrozen only the steps with that mean above freezing
      ! are integrated, and the turnover is at the end of each calendar year.
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim               ! time step [s]
      real(r8), intent(in) :: z_soi(1:nl_soil)     ! layer node depth [m]
      real(r8), intent(in) :: dz_soi(1:nl_soil)    ! layer thickness [m]
      real(r8), intent(in) :: t_soi(1:nl_soil)     ! layer temperature [K]
      logical,  intent(in) :: year_end             ! last step of a calendar year

      real(r8), parameter :: zmax = 1._r8, year_seconds = 365._r8 * 86400._r8
      real(r8) :: tcol, wsum
      integer  :: j

      if (.not. allocated(annavg_tsoi)) return
      if (ipatch < 1 .or. ipatch > size(annavg_tsoi)) return
      tcol = 0._r8; wsum = 0._r8
      do j = 1, nl_soil
         if (z_soi(j) > zmax) exit
         tcol = tcol + t_soi(j) * dz_soi(j)
         wsum = wsum + dz_soi(j)
      enddo
      if (wsum <= 0._r8) return
      tcol = tcol / wsum
      if (DEF_METHANE%q10methane_base_unfrozen) then
         ! mean over the thawed part of the calendar year
         if (tcol > tfrz) then
            tempavg_tsoi_int(ipatch) = tempavg_tsoi_int(ipatch) + tcol * deltim
            tempavg_tsoi_cnt(ipatch) = tempavg_tsoi_cnt(ipatch) + deltim
         endif
         if (year_end) then
            if (tempavg_tsoi_cnt(ipatch) > 0._r8) &
               annavg_tsoi(ipatch) = tempavg_tsoi_int(ipatch) / tempavg_tsoi_cnt(ipatch)
            tempavg_tsoi_int(ipatch) = 0._r8
            tempavg_tsoi_cnt(ipatch) = 0._r8
         endif
         return
      endif
      tempavg_tsoi_int(ipatch) = tempavg_tsoi_int(ipatch) + tcol * deltim
      tempavg_tsoi_cnt(ipatch) = tempavg_tsoi_cnt(ipatch) + deltim
      if (tempavg_tsoi_cnt(ipatch) >= year_seconds) then
         annavg_tsoi(ipatch) = tempavg_tsoi_int(ipatch) / tempavg_tsoi_cnt(ipatch)
         tempavg_tsoi_int(ipatch) = 0._r8
         tempavg_tsoi_cnt(ipatch) = 0._r8
      endif
   end subroutine methane_tsoi_annual

   !---------------------------------------------------------------------------
   pure real(r8) function methanogen_arrhenius(t, q)
      ! Arrhenius function of JULES-microbe (Chadburn et al.
      ! 2020, eq. A4, p.18 of 28): q**(0.1 Tc / (1 - Tc/T0)) with Tc in C and
      ! T0 = -273.15 C, i.e. q**(27.315 (T - 273.15) / T) with T in K; 1 at 0 C.
      real(r8), intent(in) :: t   ! temperature [K]
      real(r8), intent(in) :: q   ! Q10-like temperature sensitivity [-]

      methanogen_arrhenius = q ** (27.315_r8 * (t - 273.15_r8) / max(t, 1._r8))
   end function methanogen_arrhenius

   !---------------------------------------------------------------------------
			logical function methane_patch_is_nongrass(patchclass)
				integer, intent(in) :: patchclass

#ifdef LULC_USGS
				! USGS: woody/non-grass classes that should use the reduced
				! non-grass aerenchyma porosity.  Include mixed cropland/woodland,
				! shrubland and mixed shrub/grass classes in addition to
				! savanna/forest classes.
				methane_patch_is_nongrass = patchclass == 6 .or. patchclass == 8 .or. &
				                            patchclass == 9 .or. &
				                            (patchclass >= 10 .and. patchclass <= 15) .or. &
				                            patchclass == 18 .or. patchclass == 21
#else
				! IGBP/PFT-class builds: tree/shrub/savanna classes are non-grass.
				! Grasslands/croplands/wetlands keep the grass-like tiller porosity
				! unless a future PFT-level parameter is wired in.
				methane_patch_is_nongrass = (patchclass >= 1 .and. patchclass <= 9) .or. patchclass == 14
#endif
		end function methane_patch_is_nongrass

		!---------------------------------------------------------------------------
   ! pft_grass_share: PFT make-up of a patch that carries PFTs, for the
   ! gridded PFT and PC builds, where every natural soil patch has
   ! patchclass 1 and the class test above calls it non-grass.
   !   herb_share : grasses (PFT 12-14) and crops (PFT >= 15) over the
   !                vegetated cover; bare ground (PFT 0) has no tillers.
   !   woody_share: trees (PFT 1-8) over the whole patch, for the woody
   !                conduit of woody_conduit_area: its evidence is
   !                stem emission of flooded trees (Pangala et al. 2017), shrub
   !                stems are not measured; bare ground holds no conduit.
   ! A patch without vegetated cover gets herb_share 1 and woody_share 0.
   ! Returns .false. with the key off, in single-point builds (patchclass
   ! is the tower's IGBP class there) and for patches without PFTs
   ! (wetland, lake, urban, ice); the callers then keep the class test.
   logical function methane_patch_pft_shares(ipatch, herb_share, woody_share) result(found)
      USE MOD_LandPFT, only: patch_pft_s, patch_pft_e
      USE MOD_Vars_PFTimeInvariants, only: pftclass, pftfrac
      integer,  intent(in)  :: ipatch
      real(r8), intent(out) :: herb_share, woody_share
      integer  :: m, ps, pe
      real(r8) :: frac, f_all, f_herb, f_woody, f_tree

      found = .false.
      herb_share  = 1._r8
      woody_share = 0._r8
#ifndef SinglePoint
      if (.not. DEF_METHANE%pft_grass_share) return
      if (.not. allocated(patch_pft_s) .or. .not. allocated(patch_pft_e)) return
      if (.not. allocated(pftclass) .or. .not. allocated(pftfrac)) return
      if (ipatch < 1 .or. ipatch > size(patch_pft_s) .or. ipatch > size(patch_pft_e)) return
      ps = patch_pft_s(ipatch)
      pe = patch_pft_e(ipatch)
      if (ps < 1 .or. pe < ps .or. pe > size(pftclass) .or. pe > size(pftfrac)) return

      f_all   = 0._r8
      f_herb  = 0._r8
      f_woody = 0._r8
      f_tree  = 0._r8
      do m = ps, pe
         frac = max(pftfrac(m), 0._r8)
         f_all = f_all + frac
         if (pftclass(m) >= 1 .and. pftclass(m) <= 11) then
            f_woody = f_woody + frac
            if (pftclass(m) <= 8) f_tree = f_tree + frac
         elseif (pftclass(m) >= 12) then
            f_herb = f_herb + frac
         endif
      enddo
      if (f_herb + f_woody > 0._r8) herb_share = f_herb / (f_herb + f_woody)
      if (f_all > 0._r8) woody_share = min(f_tree / f_all, 1._r8)
      found = .true.
#endif
   end function methane_patch_pft_shares

   !---------------------------------------------------------------------------
		subroutine Tridiagonal (lbj, ubj, jtop, a, b, c, r, u)
		!-----------------------------------------------------------------------
		! DESCRIPTION:
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
		real(r8):: t_eff                  ! bounded effective temperature [K]

			do j = 0,nl_soil
				do s=1,ngases
					if (j==0) then
					t_eff = max(t_grnd, 200._r8)
					k_h = kh_theta(s)*exp(c_h(s) * (1._r8 / t_eff - 1._r8 / kh_tbase))
					! [mol/L/atm] = [mol/L/atm]*e**([K]*[1/K])
					k_h_cc(j,s) = k_h * rgasLatm * t_eff
					! [-]  = [mol/L/atm]*[L*atm/mol/K]*[K]
				else
					t_eff = max(t_soisno(j), 200._r8)
					k_h = kh_theta(s)*exp(c_h(s) * (1._r8 / t_eff - 1._r8 / kh_tbase))
					! [mol/L/atm] = [mol/L/atm]*e**([K]*[1/K])
					k_h_cc(j,s) = k_h * rgasLatm * t_eff
					! [-]  = [mol/L/atm]*[L*atm/mol/K]*[K]
				end if
			end do
		end do
	end subroutine henry_law

	subroutine split_ch4_o2_phases( dz_soisno, wliq_soisno, wice_soisno, porsl, &
									conc_methane, conc_o2, k_h_cc, idate,&
									vol_aqu,vol_ch4_storage,vol_gas,f_aqu,f_gas,&
									conc_ch4_gas,conc_ch4_aqu,conc_ch4_porsl,conc_ch4_gas_porsl,conc_ch4_aqu_porsl,&
									conc_o2_gas,conc_o2_aqu,conc_o2_porsl,conc_o2_gas_porsl,conc_o2_aqu_porsl)

		implicit none
		real(r8), intent(in) :: dz_soisno(maxsnl+1:nl_soil)     ! layer thickness [m]
		real(r8), intent(in) :: wliq_soisno(maxsnl+1:nl_soil)   ! liquid water in layers [kg/m2]
		real(r8), intent(in) :: wice_soisno(maxsnl+1:nl_soil)   ! ice in layers [kg/m2]
		real(r8), intent(in) :: porsl(1:nl_soil)                ! volumetric soil water at saturation (porosity)
		real(r8), intent(in) :: conc_methane(1:nl_soil)          ! CH4 concentration in each soil layer [mol/m3]
		real(r8), intent(in) :: conc_o2(1:nl_soil)           ! O2 concentration in each soil layer [mol/m3]
		real(r8), intent(in) :: k_h_cc(0:nl_soil,ngases)        ! ratio of mol/m3 in liquid to mol/m3 in gas [-]
		integer, intent(in) :: idate(3)

		integer :: j,s
		real(r8), parameter :: smallnumber = 1.e-12_r8
		real(r8) :: vol_ice, mobile_pore, f_ch4_storage
		real(r8), intent(out) :: vol_aqu(1:nl_soil)          ! liquid volumetric water content [m3/m3]
		real(r8), intent(out) :: vol_ch4_storage(1:nl_soil)  ! liquid/ice volume retaining CH4 [m3/m3]
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
			if (porsl(j) <= smallnumber .or. dz_soisno(j) <= smallnumber) then
				vol_aqu(j) = 0._r8
				vol_ch4_storage(j) = 0._r8
				vol_gas(j) = 0._r8
				f_aqu(j) = 0._r8
				f_gas(j) = 0._r8
				conc_ch4_gas(j) = 0._r8
				conc_ch4_aqu(j) = 0._r8
				conc_ch4_porsl(j) = 0._r8
				conc_ch4_gas_porsl(j) = 0._r8
				conc_ch4_aqu_porsl(j) = 0._r8
				conc_o2_gas(j) = 0._r8
				conc_o2_aqu(j) = 0._r8
				conc_o2_porsl(j) = 0._r8
				conc_o2_gas_porsl(j) = 0._r8
				conc_o2_aqu_porsl(j) = 0._r8
				cycle
			end if
			! ---- Compute volumetric liquid/ice/gas contents ----
			vol_aqu(j) = max(0._r8, min(wliq_soisno(j)/(dz_soisno(j)*denh2o), porsl(j)))
			! [m3/m3] = [kg/m2] / ([m] * [kg/m3])
			vol_ice = max(0._r8, min(wice_soisno(j)/(dz_soisno(j)*denice), max(porsl(j)-vol_aqu(j), 0._r8)))

			! ---- Compute volumetric gas content, excluding ice-filled pores ----
			vol_gas(j) = max(porsl(j) - vol_aqu(j) - vol_ice, 0._r8)
			mobile_pore = vol_aqu(j) + vol_gas(j)
			! Ice-filled pore space is neither gas nor liquid Henry-law storage.
			! When liquid freezes, the shrinking mobile capacity expels CH4 through
			! the conservative transport solve instead of hiding it in an ice pool.
			vol_ch4_storage(j) = vol_aqu(j)

			! ---- Compute filled proportions ----
			f_aqu(j) = vol_aqu(j)/porsl(j)
			f_gas(j) = vol_gas(j)/porsl(j)
			f_ch4_storage = vol_ch4_storage(j)/porsl(j)

			! ---- CH4 partitioning between gas and aqueous phases ----
			if (mobile_pore <= smallnumber) then
				conc_ch4_aqu(j) = 0._r8
				conc_ch4_gas(j) = 0._r8
				conc_o2_aqu(j) = 0._r8
				conc_o2_gas(j) = 0._r8
			else
				conc_ch4_aqu(j) = conc_methane(j)/(f_ch4_storage+f_gas(j)/k_h_cc(j,1))
				conc_ch4_gas(j) = conc_methane(j)/(k_h_cc(j,1)*f_ch4_storage+f_gas(j))

				! ---- O2 partitioning between gas and aqueous phases ----
				conc_o2_aqu(j) = conc_o2(j)/(f_aqu(j)+f_gas(j)/k_h_cc(j,2))
				conc_o2_gas(j) = conc_o2(j)/(k_h_cc(j,2)*f_aqu(j)+f_gas(j))
			endif

			! ---- Concentrations normalized by porosity ----
			conc_ch4_porsl(j)      = conc_methane(j)/porsl(j)
			conc_ch4_aqu_porsl(j)  = conc_ch4_aqu(j)/porsl(j)
			conc_ch4_gas_porsl(j)  = conc_ch4_gas(j)/porsl(j)
			conc_o2_porsl(j)       = conc_o2(j)/porsl(j)
			conc_o2_aqu_porsl(j)   = conc_o2_aqu(j)/porsl(j)
			conc_o2_gas_porsl(j)   = conc_o2_gas(j)/porsl(j)


		end do
	end subroutine split_ch4_o2_phases
END MODULE MOD_Tracer_Reactive_Methane_Physics
#endif
! --------- EOP ----------
