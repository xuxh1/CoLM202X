#include "define.h"

module MOD_ch4_test
!-----------------------------------------------------------------------
! !DESCRIPTION:
! Module holding routines to calculate methane fluxes
! The driver averages up to gridcell, weighting by finundated, and checks for balance errors.
! Sources, sinks, "competition" for CH4 & O2, & transport are resolved in ch4_tran.
!-----------------------------------------------------------------------
!USES:
	use MOD_Precision
	use MOD_Const_Physical, only: rgas, denh2o, denice, tfrz, grav
	use MOD_Const_ch4, only : ngases,catomw, s_con, d_con_w, d_con_g, c_h_inv, kh_theta, kh_tbase
	use MOD_Const_ch4, only : secspday, rgasm, rgasLatm
	use MOD_ch4varcon, only : replenishlakec, allowlakeprod, ch4offline
	use MOD_ch4varcon, only : usephfact, anoxicmicrosites, ch4rmcnlim 
	use MOD_ch4varcon, only : iulog, use_cn, use_nitrif_denitrif, use_lch4, use_fates_bgc, anoxia
	use MOD_ch4varcon, only : transpirationloss, use_aereoxid_prog
	use MOD_ch4varcon, only : ch4frzout
	use MOD_Namelist, only : DEF_USE_VariablySaturatedFlow
	use MOD_Vars_Global, only : maxsnl,nl_soil,spval,PI,deg2rad
	use MOD_SPMD_Task
	! 
	implicit none
	save
	! !PUBLIC MEMBER FUNCTIONS:
	! public  :: readParams
	! public  :: ch4_init_column_balance_check
	! public  :: ch4_init_gridcell_balance_check
	public  :: ch4
	
	! !PRIVATE MEMBER FUNCTIONS:
	private :: ch4_annualupdate
	private :: ch4_prod
	private :: ch4_oxid
	private :: ch4_aere
	private :: SiteOxAere
	private :: ch4_ebul
	private :: ch4_tran
	private :: Tridiagonal



	type, private :: params_type
		! ch4 production constants
		real(r8) :: q10ch4 =2              ! additional Q10 for methane production ABOVE the soil decomposition temperature relationship (2)
		real(r8) :: q10ch4base = 295._r8   ! temperature at which the effective f_ch4 actually equals the constant f_ch4 (295)
		real(r8) :: f_ch4 = 0.2            ! ratio of CH4 production to total C mineralization (0.2)
		! ! real(r8) :: rootlitfrac        ! Fraction of soil organic matter associated with roots
		real(r8) :: cnscalefactor=1        ! scale factor on CN decomposition for assigning methane flux (?)
		real(r8) :: redoxlag =30           ! Number of days to lag in the calculation of finundated_lag (30)
		real(r8) :: lake_decomp_fact =1    ! Base decomposition rate (1/s) at 25C (1)
		real(r8) :: redoxlag_vertical=30   ! time lag (days) to inhibit production for newly unsaturated layers (30)
		real(r8) :: pHmax = 9._r8          ! maximum pH for methane production(= 9._r8)
		real(r8) :: pHmin = 2.2_r8         ! minimum pH for methane production(= 2.2_r8)
		real(r8) :: oxinhib = 400          ! inhibition of methane production by oxygen (m^3/mol) (400)

		real(r8) :: mino2lim = 0.2         ! minimum anaerobic decomposition rate as a fraction of potential aerobic rate (0.2)

		! ! ch4 oxidation constants
		real(r8) :: vmax_ch4_oxid = 45.e-6_r8 * 1000._r8 / 3600._r8       ! oxidation rate constant (= 45.e-6_r8 * 1000._r8 / 3600._r8) [mol/m3-w/s];
		real(r8) :: k_m = 5.e-6_r8 * 1000._r8                 ! Michaelis-Menten oxidation rate constant for CH4 concentration 
		real(r8) :: q10_ch4_oxid = 1         ! Q10 oxidation constant (?)
		real(r8) :: smp_crit =-2.4e5_r8            ! Critical soil moisture potential(mm)
		real(r8) :: k_m_o2 =20.e-6_r8 * 1000._r8              ! Michaelis-Menten oxidation rate constant for O2 concentration
		real(r8) :: k_m_unsat = 5.e-6_r8 * 1000._r8 / 10._r8           ! Michaelis-Menten oxidation rate constant for CH4 concentration
		real(r8) :: vmax_oxid_unsat = 45.e-6_r8 * 1000._r8 / 3600._r8 / 10._r8     ! (= 45.e-6_r8 * 1000._r8 / 3600._r8 / 10._r8) [mol/m3-w/s]

		! ! ch4 aerenchyma constants
		real(r8) :: aereoxid =1            ! fraction of methane flux entering aerenchyma rhizosphere that will be(?)

		! ! oxidized rather than emitted
		real(r8) :: scale_factor_aere = 1   ! scale factor on the aerenchyma area for sensitivity tests (1)
		! real(r8) :: nongrassporosratio = 1/3   ! Ratio of root porosity in non-grass to grass, used for aerenchyma transport
		real(r8) :: unsat_aere_ratio = 0.05_r8 / 0.3_r8    ! Ratio to multiply upland vegetation aerenchyma porosity by compared to inundated systems (= 0.05_r8 / 0.3_r8)
		real(r8) :: porosmin = 0.05_r8            ! minimum aerenchyma porosity (unitless)(= 0.05_r8) 

		! ! ch4 ebbulition constants
		real(r8) :: vgc_max  =0.15            ! ratio of saturation pressure triggering ebullition (0.15)

		! ! ch4 transport constants
		real(r8) :: satpow  =2             ! exponent on watsat for saturated soil solute diffusion (2?)
		real(r8) :: scale_factor_gasdiff = 2! For sensitivity tests; convection would allow this to be > 1(?)
		real(r8) :: scale_factor_liqdiff = 2! For sensitivity tests; convection would allow this to be > 1(?)
		real(r8) :: capthick = 100._r8            ! min thickness before assuming h2osfc is impermeable (mm) (= 100._r8)

		! ! additional constants
		! real(r8) :: f_sat =0.95               ! volumetric soil water defining top of water table or where production is allowed (=0.95)
		real(r8) :: qflxlagd  = 30._r8          ! days to lag qflx_surf_lag in the tropics (days) ( = 30._r8)
		real(r8) :: highlatfact = 2._r8         ! multiple of qflxlagd for high latitudes	(= 2._r8)	
		real(r8) :: q10lakebase = 298._r8         ! (K) base temperature for lake CH4 production (= 298._r8)
		real(r8) :: atmch4  = 1.7e-6_r8           ! Atmospheric CH4 mixing ratio to prescribe if not provided by the atmospheric model (= 1.7e-6_r8) (mol/mol)
		real(r8) :: rob = 3._r8                 ! ratio of root length to vertical depth ("root obliquity") (= 3._r8)
		real(r8) :: om_frac_sf = 1           ! Scale factor for organic matter fraction (unitless)(?)
	end type params_type
	type(params_type), private ::  params_inst

contains

	!-----------------------------------------------------------------------
	subroutine ch4 (ipatch,patchtype,&!input
		patchlonr,patchlatr,&
		lb,nl_soil,maxsnl,snl,&
		deltim,&
		z_soisno,dz_soisno,zi_soisno,t_soisno,t_grnd,wliq_soisno,wice_soisno,&
		forc_t,forc_pbot,forc_po2m,forc_pco2m,&
		zwt,rootfr,snowdp,wat,rsur,etr,lakedepth,lake_icefrac,wdsrf,bsw,&
		smp,porsl,lai,&
		annsum_npp,rr,&
		idate,agnpp,bgnpp,somhr,&
		crootfr,lithr,hr_vr,o_scalar,fphr,pot_f_nit_vr,pH,&
		rootr,&
		cellorg,t_h2osfc,organic_max,&
		c_atm,ch4_surf_flux_tot,net_methane,&!output
		annavg_agnpp,annavg_bgnpp,annavg_somhr,annavg_finrw,&
		ch4_prod_depth,o2_decomp_depth,&
		ch4_oxid_depth,o2_oxid_depth,&
		ch4_aere_depth,ch4_tran_depth,o2_aere_depth,&
		ch4_ebul_depth,&
		o2stress,ch4stress,ch4_surf_aere,ch4_surf_ebul,ch4_surf_diff,ch4_ebul_total,&
		ch4_first_time,totcolch4,forc_pch4m,grnd_ch4_cond,conc_o2,conc_ch4,layer_sat_lag,lake_soilc,&!inout
		tempavg_agnpp,tempavg_bgnpp,annsum_counter,&
		tempavg_somhr,tempavg_finrw)
!=======================================================================
! !DESCRIPTION:
! Driver for the methane emissions model
!=======================================================================
		use MOD_Precision

		implicit none

!-----------------------------------------------------------------------

		integer, INTENT(in) :: &
			ipatch           , &! patch index
			patchtype           ! land patch type (0=soil, 1=urban or built-up, 2=wetland,
									  ! 3=land ice, 4=land water bodies, 99=ocean
			! istep            , &! the i time step

		real(r8), intent(in) :: &
			patchlonr   ,&! logitude in radians
			patchlatr     ! latitude in radians

		integer, INTENT(in) :: &
			lb               , &! lower bound of array (snl+1)
			nl_soil          , &! upper bound of array (10)
			maxsnl			 , &!  max number of snow layers (-5)
			snl				    !  number of snow layers     (-5~-1)

		real(r8), INTENT(in) :: &
			deltim                  , &! land model time step (sec)
			z_soisno (lb:nl_soil)    , &! layer depth (m)
			dz_soisno(lb:nl_soil)    , &! layer thickness (m)
			zi_soisno(lb-1:nl_soil)    , &! interface level below a "z" level (m)

			t_soisno (lb:nl_soil)    , &! soil temperature (Kelvin)
			t_grnd                 		, &! ground surface temperature [k]
			wliq_soisno(lb:nl_soil)	, &! liquid water in layers [kg/m2]
			wice_soisno(lb:nl_soil) , &! ice lens in layers [kg/m2]

			forc_t                  , &! temperature at reference height [kelvin]
			forc_pbot               , &! atm bottom level pressure (or reference height) (pa)
			forc_po2m               , &! O2 concentration in atmos. (pascals)
			forc_pco2m              , &! CO2 concentration in atmos. (pascals)

			zwt                     , &! the depth from ground (soil) surface to water table [m]

			rootfr   (1:nl_soil)    , &! fraction of roots in each soil layer

			snowdp                  , &! snow depth [m]
			wat                     , &! total water storage [mm]
			rsur                    , &! surface runoff (mm h2o/s)
			etr                     , &! transpiration rate [mm/s]
			lakedepth               , &! lake depth
			lake_icefrac (1:nl_soil), &! lake mass fraction of lake layer that is frozen
			wdsrf                   , &! depth of surface water [mm]
			bsw      (1:nl_soil)   	, &! Clapp and Hornberger "b" (nlevgrnd)             

			smp      (1:nl_soil)    , &! soil matrix potential [mm]
			porsl    (1:nl_soil)    , &! volumetric soil water at saturation (porosity)
			lai                        ! leaf area index [m2/m2]

		real(r8), INTENT(in) :: &
			annsum_npp              , &! annual sum NPP (gC/m2/yr)
			rr                         ! root respiration (fine root MR + total root GR) (gC/m2/s)
			! froot_xsmr               ! fine root maintenance respiration storage C due to available C deficit (gC m-2 s-1)

		!-------------------ch4_annualupdate------------------------------
		integer, INTENT(in) :: &
			idate(3)             	   ! current date (year, days of the year, seconds of the day)

		real(r8), INTENT(in) :: &			
			agnpp                   , &! aboveground NPP (gC/m2/s)
			bgnpp                   , &! belowground NPP (gC/m2/s)
			somhr                      ! (gC/m2/s) soil organic matter heterotrophic respiration

		!-------------------ch4_prod------------------------------
		real(r8), INTENT(in) :: &
			crootfr  (1:nl_soil)    , &! fraction of roots for carbon in each soil layer
			lithr                   , &! (gC/m2/s) litter heterotrophic respiration        
			hr_vr    (1:nl_soil)    , &! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
			o_scalar (1:nl_soil)    , &! fraction by which decomposition is limited by anoxia
			fphr     (1:nl_soil)    , &! fraction of potential heterotrophic respiration 

			pot_f_nit_vr(1:nl_soil) , &! (gN/m3/s) potential soil nitrification flux 
			pH                         ! soil water pH                                     

		!-------------------ch4_aere------------------------------
		real(r8), INTENT(in) :: &
			rootr    (1:nl_soil)       ! effective fraction of roots in each soil layer (SMS method only)
		
		!-------------------ch4_aere------------------------------
		real(r8), INTENT(in) :: &
			cellorg  (1:nl_soil)   		, &! column 3D org (kg/m^3 organic matter)
			t_h2osfc               		, &! surface water temperature               
			organic_max               		! organic matter content (kg m-3) where soil is assumed to act like peat


		!-------------------------------------------------------------------------------------
		real(r8), INTENT(out) :: &
		! 	! ch4_oxid_depth (1:nl_soil)   , &! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) 
		! 	! o2_oxid_depth  (1:nl_soil)   , &! O2 consumption rate via oxidation in each soil layer (mol/m3/s) 
		! 	! co2_oxid_depth (1:nl_soil)   , &! CO2 production rate via oxidation in each soil layer (mol/m3/s) 
		! 	! o2_decomp_depth(1:nl_soil)      ! O2 consumption during decomposition in each soil layer (mol/m3/s)
		

		! 	ch4_surf_diff                , &! CH4 surface flux (mol/m2/s)                            
		! 	ch4_surf_ebul                , &! CH4 ebullition to atmosphere (mol/m2/s)                 
		! 	ch4_surf_aere                , &! Total column CH4 aerenchyma (mol/m2/s)            
		
		! 	ch4_oxid_depth (1:nl_soil)   , &! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) 
		! 	ch4_prod_depth (1:nl_soil)   , &! production of CH4 in each soil layer (mol/m3/s)

		! 	soilc          (1:nl_soil)   , &! total soil organic matter found in level (g C / m3) 
		! 	! conc_o2_lake   (1:nl_soil)   , &! O2 conc  in each soil layer (mol/m3)
		! 	! ch4_dfsat_flux               , &! CH4 flux to atm due to decreasing finundated (kg C/m^2/s) [+]

		! 	rsur_lag                     , &! time-lagged surface runoff (mm H2O /s)
		! 	layer_lag  (1:nl_soil)       , &! Lagged saturation status of soil layer (1 = sat)
			c_atm      (1:3)             , &! CH4, O2, CO2 atmospheric conc  (mol/m3)         
		! 	ch4co2                       , &! CO2 production from CH4 oxidation (g C/m**2/s)
		! 	ch4prod                      , &! average CH4 production (g C/m^2/s)       
			ch4_surf_flux_tot            , &! CH4 flux to atm. (kg C/m**2/s)
			net_methane                     ! average net methane correction to CO2 flux (g C/m^2/s)

		!-------------------ch4_annualupdate------------------------------
		real(r8), INTENT(out) :: &
			annavg_agnpp            , &! annual average above-ground NPP (gC/m2/s)         
			annavg_bgnpp            , &! annual average below-ground NPP (gC/m2/s)         
			annavg_somhr            , &! annual average SOM heterotrophic resp. (gC/m2/s)  
			annavg_finrw               ! respiration-weighted annual average of finundated 
		
		!-------------------ch4_prod------------------------------
		real(r8), INTENT(out) :: &            
			ch4_prod_depth    (1:nl_soil)         , &! production of CH4 in each soil layer  (mol/m3/s)
			o2_decomp_depth   (1:nl_soil)            ! O2 consumption during decomposition in each soil layer (mol/m3/s)

		!-------------------ch4_oxid------------------------------
		real(r8), INTENT(out) :: &
			ch4_oxid_depth (1:nl_soil)       , &! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) 
			o2_oxid_depth  (1:nl_soil)          ! O2 consumption rate via oxidation in each soil layer (mol/m3/s) 
		
		!-------------------ch4_aere------------------------------
		real(r8), INTENT(out) :: &
			ch4_aere_depth  (1:nl_soil)  , &! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) 
			ch4_tran_depth  (1:nl_soil)  , &! CH4 loss rate via transpiration in each soil layer (mol/m3/s) 
			o2_aere_depth   (1:nl_soil)     ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s) 

		!-------------------ch4_ebul------------------------------
		real(r8), INTENT(out) :: &
			ch4_ebul_depth (1:nl_soil)      ! CH4 loss rate via ebullition in each soil layer (mol/m3/s)

		!-------------------ch4_tran------------------------------
		real(r8), INTENT(out) :: &
			o2stress          (1:nl_soil)  , &! Output: Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs
			ch4stress         (1:nl_soil)  , &! Output: Ratio of methane available to the total per-timestep methane sinks 
			ch4_surf_aere                  , &! Output: Total column CH4 aerenchyma (mol/m2/s)
			ch4_surf_ebul                  , &! Output: CH4 ebullition to atmosphere (mol/m2/s)
			ch4_surf_diff                  , &! Output: CH4 surface flux (mol/m2/s)
			ch4_ebul_total                    ! Output: Total column CH4 ebullition (mol/m2/s)


		!-------------------------------------------------------------------------------------
		logical, INTENT(inout) ::&
			ch4_first_time
		
		real(r8), INTENT(inout) :: &
			totcolch4               , &! total methane in soil column, start of timestep (g C / m^2)
			forc_pch4m              , &! CH4 concentration in atmos. (pascals)
			grnd_ch4_cond           , &! tracer conductance for boundary layer [m/s]
			conc_o2  (1:nl_soil)    , &! O2 conc in each soil layer (mol/m3) 
			conc_ch4   (1:nl_soil)  , &! CH4 conc in each soil layer (mol/m3) 
			layer_sat_lag(1:nl_soil), &
			lake_soilc  (1:nl_soil)    ! total soil organic matter found in level (g C / m^3) (nl_soil)


		!-------------------ch4_annualupdate------------------------------
		real(r8), INTENT(inout) :: &
			tempavg_agnpp           , &! temporary average above-ground NPP (gC/m2/s)      
			tempavg_bgnpp           , &! temporary average below-ground NPP (gC/m2/s)      
			annsum_counter          , &! seconds since last annual accumulator turnover    
			tempavg_somhr           , &! temporary average SOM heterotrophic resp. (gC/m2/s)
			tempavg_finrw              ! respiration-weighted annual average of finundated 
!-----------------------Local Variables------------------------------
		integer  :: j,s      ! indices

		integer  :: sat                     ! 0 = unsatured, 1 = saturated
		integer  :: finundated              ! fractional inundated area, =sat(0 or 1)
		integer  :: jwt                     ! index of the soil layer right above the water table (-)
		real(r8) :: vol_liq  (1:nl_soil)    ! liquid volumetric water content
		real(r8) :: dzmm(1:nl_soil)         ! dz m to mm

		real(r8) :: lon,lat                 ! lon,lat

		real(r8) :: ch4_prod_tot            ! CH4 production for column (g C/m**2/s)
		real(r8) :: ch4_oxid_tot            ! CH4 oxidation for column (g C/m**2/s)

		real(r8) :: total                   ! diff + aere + ebul
		real(r8) :: dfsat
		real(r8) :: fsat_bef                ! finundated from previous timestep
		real(r8) :: errch4                  ! g C / m^2
		!real(r8) :: zwt_actual
		! real(r8) :: qflxlags              ! Time to lag qflx_surf_lag (s)
		real(r8) :: redoxlag                ! Redox time lag 
		real(r8) :: redoxlag_vertical       ! Vertical redox lag time 
		real(r8) :: atmch4                  ! Atmospheric CH4 mixing ratio to
											! prescribe if not provided by the atmospheric model (= 1.7e-6_r8) (mol/mol)
		real(r8) :: redoxlags               ! Redox time lag in s
		real(r8) :: redoxlags_vertical      ! Vertical redox lag time in s
		real(r8) :: qflxlagd                ! days to lag qflx_surf_lag in the tropics (days)
		real(r8) :: highlatfact             ! multiple of qflxlagd for high latitudes
		integer  :: dummyfilter(1)          ! empty filter
		character(len=32) :: subname='ch4'  ! subroutine name

		real(r8) :: totcolch4_bef           ! total methane in soil column, start of timestep (g C / m^2)
!-----------------------------------------------------------------------

		! Set parameters
		redoxlag          = params_inst%redoxlag
		redoxlag_vertical = params_inst%redoxlag_vertical
		atmch4            = params_inst%atmch4
		qflxlagd          = params_inst%qflxlagd
		highlatfact       = params_inst%highlatfact

		redoxlags = redoxlag*secspday ! days --> s
		redoxlags_vertical = redoxlag_vertical*secspday ! days --> s
		
		totcolch4_bef = totcolch4
		totcolch4 = 0
!-----------------------------------------------------------------------
        ! compute jwt index
		! The layer index of the first unsaturated layer,
		! i.e., the layer right above the water table
		jwt = nl_soil
		! allow jwt to equal zero when zwt is in top layer
		do j = 1, nl_soil
			if(zwt <= zi_soisno(j)) then
				jwt = j-1
				exit
			end if
		end do

		if(jwt==0) then
			sat = 1
		else
			sat = 0
		end if
	
		finundated = sat

		do j= 1, nl_soil
			dzmm(j) = dz_soisno(j)*denh2o
			vol_liq(j) = wliq_soisno(j)/dzmm(j)
        end do
        
        lon = patchlonr/deg2rad
        lat = patchlatr/deg2rad
        print*, lon,lat
!-----------------------------------------------------------------------
		! Initialize local fluxes to zero: necessary for columns outside the filters because averaging up to gridcell will be done
		ch4_surf_flux_tot     = 0._r8
		ch4_prod_tot          = 0._r8
		ch4_oxid_tot          = 0._r8

		! Adjustment to NEE for methane production - oxidation
		net_methane           = 0._r8
	
		if (ch4offline) then
			forc_pch4m = atmch4*forc_pbot
		else
			if (forc_pch4m == 0._r8) then
				write(6,*) 'not using ch4offline, but methane concentration not passed from the atmosphere', &
				'to land model! CLM Model is stopping.'
				CALL CoLM_stop ()
			end if
		end if

		c_atm(1) =  forc_pch4m / rgasm / forc_t ! [mol/m3 air]
		c_atm(2) =  forc_po2m  / rgasm / forc_t ! [mol/m3 air]
		c_atm(3) =  forc_pco2m / rgasm / forc_t ! [mol/m3 air]
	
		!!!! Begin biochemistry
	
		! First for soil

		! Do CH4 Annual Averages
		call ch4_annualupdate(idate, finundated, deltim,  agnpp, bgnpp, somhr, &
			annavg_agnpp, annavg_bgnpp, annavg_somhr,  annavg_finrw, &
			tempavg_agnpp,tempavg_bgnpp,annsum_counter,tempavg_somhr, tempavg_finrw)
 
!-------------------------------------------------
! Loop
!-------------------------------------------------
		! Update lagged saturation status of layer
		if (sat == 0) then ! unsaturated
			do j=1,nl_soil
				if (j > jwt .and. redoxlags_vertical > 0._r8) then ! saturated currently
					layer_sat_lag(j) = layer_sat_lag(j) * exp(-deltim/redoxlags_vertical) &
						+ (1._r8 - exp(-deltim/redoxlags_vertical))
				else if (redoxlags_vertical > 0._r8) then
					layer_sat_lag(j) = layer_sat_lag(j) * exp(-deltim/redoxlags_vertical)
				else if (j > jwt) then  ! redoxlags_vertical = 0
					layer_sat_lag(j) = 1._r8
				else
					layer_sat_lag(j) = 0._r8
				end if
			end do
		end if ! saturated no change

		! calculate CH4 production in each soil layer
		call ch4_prod (patchtype,sat,finundated,jwt,nl_soil,rr,deltim,& !input
			z_soisno,dz_soisno,zi_soisno,t_soisno,&
			lai,conc_o2,rootfr,annavg_finrw,&
			crootfr,somhr,lithr,hr_vr,o_scalar,fphr,pot_f_nit_vr,&
			pH,lake_soilc,&
			ch4_prod_depth,o2_decomp_depth,&!output
			layer_sat_lag)!inout

		! calculate CH4 oxidation in each soil layer
		call ch4_oxid (nl_soil,  jwt,  sat,  deltim,  z_soisno,  dz_soisno,  zi_soisno, &
			t_soisno,      smp,      vol_liq,    porsl,   conc_o2,   conc_ch4,              &
			ch4_oxid_depth,          o2_oxid_depth) 

		! calculate CH4 aerenchyma losses in each soil layer
		call ch4_aere (patchtype, nl_soil, jwt, sat, lai,     deltim, &
			z_soisno, dz_soisno,  zi_soisno,     t_soisno, vol_liq, porsl,  &
			rootfr,   rootr, etr, grnd_ch4_cond, c_atm,    annsum_npp,      &
			annavg_agnpp,    annavg_bgnpp, conc_o2, conc_ch4, ch4_prod_depth,&
			ch4_aere_depth, ch4_tran_depth, o2_aere_depth)

		! calculate CH4 ebullition losses in each soil layer
		call ch4_ebul (patchtype, nl_soil, jwt, sat, deltim, &
			z_soisno, dz_soisno, zi_soisno, lakedepth, forc_pbot,&
			t_soisno, lake_icefrac, porsl, wdsrf, conc_ch4,&
			ch4_ebul_depth)

		! Solve CH4 reaction/diffusion equation 
		! Competition for oxygen will occur here.
		call ch4_tran (patchtype,&
			lb, nl_soil, snl, maxsnl,jwt, sat,&
			lon,lat,deltim, z_soisno, dz_soisno, zi_soisno,  t_soisno, t_grnd, &
			vol_liq, porsl, wliq_soisno, wice_soisno, wdsrf,bsw, c_atm, ch4_prod_depth, o2_aere_depth,&
			cellorg,t_h2osfc, organic_max, &
			o2stress, ch4stress, ch4_surf_aere, ch4_surf_ebul, ch4_surf_diff, ch4_ebul_total, &
			ch4_oxid_depth, ch4_aere_depth, ch4_ebul_depth, &
			grnd_ch4_cond, o2_oxid_depth, o2_decomp_depth, conc_o2, conc_ch4 )
	
!-------------------------------------------------
! Now do over lakes
!-------------------------------------------------
 
		if (patchtype == 4) then
			sat = 1
			jwt = 0
	
			! calculate CH4 production in each soil layer
			call ch4_prod (patchtype,sat,finundated,jwt,nl_soil,rr,deltim,& !input
				z_soisno,dz_soisno,zi_soisno,t_soisno,&
				lai,conc_o2,rootfr,annavg_finrw,&
				crootfr,somhr,lithr,hr_vr,o_scalar,fphr,pot_f_nit_vr,&
				pH,lake_soilc,&
				ch4_prod_depth,o2_decomp_depth,&!output
				layer_sat_lag)!inout

			! calculate CH4 oxidation in each soil layer
			call ch4_oxid (nl_soil,  jwt,  sat,  deltim,  z_soisno,  dz_soisno,  zi_soisno, &
				t_soisno,      smp,      vol_liq,    porsl,   conc_o2,   conc_ch4,              &
				ch4_oxid_depth,          o2_oxid_depth) 

			! calculate CH4 aerenchyma losses in each soil layer
			call ch4_aere (patchtype, nl_soil, jwt, sat, lai,     deltim, &
				z_soisno, dz_soisno,  zi_soisno,     t_soisno, vol_liq, porsl,  &
				rootfr,   rootr, etr, grnd_ch4_cond, c_atm,    annsum_npp,      &
				annavg_agnpp,    annavg_bgnpp, conc_o2, conc_ch4, ch4_prod_depth,&
				ch4_aere_depth, ch4_tran_depth, o2_aere_depth)

			! calculate CH4 ebullition losses in each soil layer
			call ch4_ebul (patchtype, nl_soil, jwt, sat, deltim, &
				z_soisno, dz_soisno, zi_soisno, lakedepth, forc_pbot,&
				t_soisno, lake_icefrac, porsl, wdsrf, conc_ch4,&
				ch4_ebul_depth)

			! Solve CH4 reaction/diffusion equation 
			! Competition for oxygen will occur here.
			call ch4_tran (patchtype,&
				lb, nl_soil, snl, maxsnl,jwt, sat,&
				lon,lat,deltim, z_soisno, dz_soisno, zi_soisno,  t_soisno, t_grnd, &
				vol_liq, porsl, wliq_soisno, wice_soisno, wdsrf,bsw, c_atm, ch4_prod_depth, o2_aere_depth,&
				cellorg,t_h2osfc, organic_max, &
				o2stress, ch4stress, ch4_surf_aere, ch4_surf_ebul, ch4_surf_diff, ch4_ebul_total, &
				ch4_oxid_depth, ch4_aere_depth, ch4_ebul_depth, &
				grnd_ch4_cond, o2_oxid_depth, o2_decomp_depth, conc_o2, conc_ch4 )
	
		end if
 
!-------------------------------------------------------------------------------
! Average up to gridcell flux and column oxidation and production rate.
!-------------------------------------------------------------------------------

		! First weight the soil columns by finundated.
		do j=1,nl_soil
			if (j == 1) then
				total = ch4_surf_diff + ch4_surf_aere + ch4_surf_ebul
				ch4_surf_flux_tot = total* catomw / 1000._r8
				!Convert from mol to kg C
				! ch4_oxid_tot and ch4_prod_tot are initialized to zero above
			end if

			ch4_oxid_tot = ch4_oxid_tot + ch4_oxid_depth(j) * dz_soisno(j) * catomw
			!Convert from mol to g C
			ch4_prod_tot = ch4_prod_tot + ch4_prod_depth(j) * dz_soisno(j) * catomw
			!Convert from mol to g C
			if (j == nl_soil) then
				! Adjustment to NEE flux to atm. for methane production
				net_methane = net_methane - ch4_prod_tot
				! Adjustment to NEE flux to atm. for methane oxidation
				net_methane = net_methane + ch4_oxid_tot
			end if
		end do
	
		! ! Correct for discrepancies in CH4 concentration from changing finundated
		! ch4_surf_flux_tot = ch4_surf_flux_tot + ch4_dfsat_flux

	
		if (patchtype==4) then
			do j=1,nl_soil
				if (j == 1) then
					! ch4_oxid_tot and ch4_prod_tot are initialized to zero above
					total = ch4_surf_diff + ch4_surf_aere + ch4_surf_ebul
					ch4_surf_flux_tot = total*catomw / 1000._r8
				end if
	
				ch4_oxid_tot = ch4_oxid_tot + ch4_oxid_depth(j)*dz_soisno(j)*catomw
				ch4_prod_tot = ch4_prod_tot + ch4_prod_depth(j)*dz_soisno(j)*catomw
	
				if (.not. replenishlakec) then
					!Adjust lake_soilc for production.
					lake_soilc(j) = lake_soilc(j) - 2._r8*ch4_prod_depth(j)*deltim*catomw
					! Factor of 2 is for CO2 that comes off with CH4 because of stoichiometry
				end if
	
				if (j == nl_soil) then
				! Adjustment to NEE flux to atm. for methane production
				if (.not. replenishlakec) then
					net_methane = net_methane + ch4_prod_tot
					! Here this is positive because it is actually the CO2 that comes off with the methane
					! NOTE THIS MODE ASSUMES TRANSIENT CARBON SUPPLY FROM LAKES; COUPLED MODEL WILL NOT CONSERVE CARBON
					! IN THIS MODE.
				else ! replenishlakec
					net_methane = net_methane - ch4_prod_tot
					! Keep total C constant, just shift from CO2 to methane
				end if

				! Adjustment to NEE flux to atm. for methane oxidation
				net_methane = net_methane + ch4_oxid_tot

				end if
			end do
		end if  ! ch4_surf_flux_tot, ch4_oxid_tot, and ch4_prod_tot should be initialized to 0 above if .not. allowlakeprod
	
		! Finalize CH4 balance and check for errors
	
		do j = 1, nl_soil
			totcolch4 = totcolch4 + conc_ch4(j)*dz_soisno(j)*catomw
			! mol CH4 --> g C
		end do
	
		! Column level balance
	
		if (.not. ch4_first_time) then
			! Check balance
			errch4 = totcolch4 - totcolch4_bef - deltim*(ch4_prod_tot - ch4_oxid_tot - ch4_surf_flux_tot*1000._r8) 
			! kg C --> g C
			if (abs(errch4) > 1.e-7_r8) then ! g C / m^2 / timestep
				! write(6,*)'Patch-level CH4 Conservation Error in CH4Mod driver, istep,  errch4 (gC /m^2.timestep)',istep,errch4
				write(6,*)'Latrad,Lonrad,Patchtype=',patchlatr,patchlonr,patchtype
				write(6,*)'totcolch4                 = ', totcolch4
				write(6,*)'totcolch4_bef             = ', totcolch4_bef
				write(6,*)'deltim*ch4_prod_tot           = ', deltim*ch4_prod_tot
				write(6,*)'deltim*ch4_oxid_tot           = ', deltim*ch4_oxid_tot
				write(6,*)'deltim*ch4_surf_flux_tot*1000 = ', deltim*ch4_surf_flux_tot*1000._r8
				CALL CoLM_stop ()
			end if
		end if
	
		if (allowlakeprod) then
			if (.not. ch4_first_time) then
				! Check balance
				errch4 = totcolch4 - totcolch4_bef - deltim*(ch4_prod_tot - ch4_oxid_tot - ch4_surf_flux_tot*1000._r8) 
				! kg C --> g C
				if (abs(errch4) > 1.e-7_r8) then ! g C / m^2 / timestep
					! write(6,*)'Column-level CH4 Conservation Error in CH4Mod driver for lake column, istep, errch4 (gC/m^2.timestep)',istep,errch4
					write(6,*)'Latrad,Lonrad,Patchtype=',patchlatr,patchlonr,patchtype
					write(6,*)'totcolch4                = ', totcolch4
					write(6,*)'totcolch4_bef            = ', totcolch4_bef
					write(6,*)'deltim*ch4_prod_tot           = ', deltim*ch4_prod_tot
					write(6,*)'deltim*ch4_oxid_tot           = ', deltim*ch4_oxid_tot
					write(6,*)'deltim*ch4_surf_flux_tot*1000 = ', deltim*ch4_surf_flux_tot*1000._r8
					CALL CoLM_stop ()
				end if
			end if
		end if
	
		!
		! Gricell level balance
		!
	
		! Skip the check if dynamic lakes are on and it's
		! - the beginning of a new year (ok for restart runs) OR
		! - the beginning of a simulation (needed for hybrid/startup runs)
		! See (https://github.com/ESCOMP/CTSM/issues/43#issuecomment-1282609233)
		! 
		! if ( is_beg_curr_year() .and. get_do_transient_lakes() .or. &
		! 	is_first_step() .and. get_do_transient_lakes() )then
		! 	ch4_first_time = .true.
		! end if
	

		! if (.not. ch4_first_time) then
		! 	! Check balance
		! 	errch4 = totcolch4 - totcolch4_bef + deltim * &
		! 	(nem_grc + ch4_surf_flux_tot * 1000._r8)  ! kg C --> g C

		! 	if (abs(errch4) > 1.e-7_r8) then  ! g C / m^2 / timestep
		! 		! write(6,*)'Gridcell-level CH4 Conservation Error in CH4Mod driver, istep, errch4 (gC /m^2.timestep)', &
		! 		! istep, errch4
		! 		write(6,*)'Latrad,Lonrad=',patchlatr,patchlonr
		! 		write(6,*)'totcolch4     =', totcolch4
		! 		write(6,*)'totcolch4_bef =', totcolch4_bef
		! 		write(6,*)'deltim * nem_grc   =', deltim * nem_grc
		! 		write(6,*)'deltim * ch4_surf_flux_tot * 1000 =', deltim * ch4_surf_flux_tot * 1000._r8
		! 		CALL CoLM_stop ()
		! 	end if
		! end if
	
		ch4_first_time = .false.

	end subroutine ch4


!-----------------------------------------------------------------------
	subroutine ch4_annualupdate(idate, finundated, deltim,  agnpp, bgnpp, somhr, &
		annavg_agnpp, annavg_bgnpp, annavg_somhr,  annavg_finrw, &
		tempavg_agnpp,tempavg_bgnpp,annsum_counter,tempavg_somhr, tempavg_finrw)
!
! !DESCRIPTION: Annual mean fields.
!
        use MOD_Precision
		implicit none

!-----------------------Argument----------------------------------------
		integer, INTENT(in) :: &
			idate(3)             , &! current date (year, days of the year, seconds of the day)
			finundated              ! fractional inundated area, =sat(0 or 1)

		real(r8), INTENT(in) :: &
			deltim                  , &! land model time step (sec)

			agnpp                   , &! aboveground NPP (gC/m2/s)
			bgnpp                   , &! belowground NPP (gC/m2/s)
			somhr                      ! (gC/m2/s) soil organic matter heterotrophic respiration

		real(r8), INTENT(out) :: &
			annavg_agnpp            , &! output: annual average above-ground NPP (gC/m2/s)         
			annavg_bgnpp            , &! output: annual average below-ground NPP (gC/m2/s)         
			annavg_somhr            , &! output: annual average SOM heterotrophic resp. (gC/m2/s)  
			annavg_finrw               ! output: respiration-weighted annual average of finundated 
		
		real(r8), INTENT(inout) :: &
			tempavg_agnpp           , &! inout: temporary average above-ground NPP (gC/m2/s)      
			tempavg_bgnpp           , &! inout: temporary average below-ground NPP (gC/m2/s)      
			annsum_counter          , &! inout: seconds since last annual accumulator turnover    
			tempavg_somhr           , &! inout: temporary average SOM heterotrophic resp. (gC/m2/s)
			tempavg_finrw              ! inout: respiration-weighted annual average of finundated 
!-----------------------Local Variables------------------------------         
		real(r8):: secsperyear         ! seconds since this year
!-----------------------------------------------------------------------
		! set time steps
		secsperyear = idate(2)*secspday + idate(3)

		annsum_counter = annsum_counter + deltim

		if (annsum_counter >= secsperyear) then

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
		else
			tempavg_somhr  = tempavg_somhr + deltim/secsperyear * somhr
			tempavg_finrw  = tempavg_finrw + deltim/secsperyear * finundated * somhr
		end if
		

		if (annsum_counter >= secsperyear) then
			
			annavg_agnpp = tempavg_agnpp
			tempavg_agnpp = 0._r8
			
			annavg_bgnpp = tempavg_bgnpp
			tempavg_bgnpp = 0._r8
			
		else
			tempavg_agnpp = tempavg_agnpp + deltim/secsperyear * agnpp
			tempavg_bgnpp = tempavg_bgnpp + deltim/secsperyear * bgnpp
		end if
		
		! column loop
		if (annsum_counter >= secsperyear) annsum_counter = 0._r8

	end subroutine ch4_annualupdate


!-----------------------------------------------------------------------
	subroutine ch4_prod (patchtype,sat,finundated,jwt,nl_soil,rr,deltim,& !input
		z_soisno,dz_soisno,zi_soisno,t_soisno,&
		lai,conc_o2,rootfr,annavg_finrw,&
		crootfr,somhr,lithr,hr_vr,o_scalar,fphr,pot_f_nit_vr,&
		pH,lake_soilc,&
		ch4_prod_depth,o2_decomp_depth,&
		layer_sat_lag)!output
  !-----------------------DESCRIPTION------------------------------
  ! Production is done below the water table, based on CN heterotrophic respiration.
  ! O2 is consumed by roots & by heterotrophic aerobes.
  ! Production is done separately for sat & unsat, and is adjusted for temperature, seasonal inundation,
  ! pH (optional), & redox lag factor.
  !-----------------------Argument---------- ------------------------------
		integer, INTENT(in) :: &
			patchtype            ! land patch type (0=soil, 1=urban or built-up, 2=wetland,
									! 3=land ice, 4=land water bodies, 99=ocean
  
		integer , INTENT(in) :: &
			sat                     , &! 0 = unsaturated; 1 = saturated 
			finundated              , &! fractional inundated area in soil column 
			jwt                     , &! index of the soil layer right above the water table (-) 
		  	nl_soil                    ! number of soil layers
  
		real(r8), INTENT(in) :: &
			rr                      , &! root respiration (fine root MR + total root GR) (gC/m2/s)

	
			deltim                  , &! land model time step (sec)
			z_soisno (1:nl_soil)    , &! layer depth (m)
			dz_soisno(1:nl_soil)    , &! layer thickness (m)
			zi_soisno(0:nl_soil)    , &! interface level below a "z" level (m)
	
			t_soisno (1:nl_soil)    , &! soil temperature (Kelvin)
	
			lai                     , &! leaf area index [m2/m2]
			conc_o2   (1:nl_soil)   , &! Input: O2 conc in each soil layer (mol/m3) (nl_soil)   
			rootfr   (1:nl_soil)    , &! fraction of roots in each soil layer
	
			annavg_finrw            ! Input: respiration-weighted annual average of finundated 
	
		real(r8), INTENT(in) :: &
			crootfr  (1:nl_soil)    , &! fraction of roots for carbon in each soil layer
	
			somhr                   , &! (gC/m2/s) soil organic matter heterotrophic respiration
			lithr                   , &! (gC/m2/s) litter heterotrophic respiration        
			hr_vr    (1:nl_soil)    , &! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
			o_scalar (1:nl_soil)    , &! fraction by which decomposition is limited by anoxia
			fphr     (1:nl_soil)    , &! fraction of potential heterotrophic respiration 
	
			pot_f_nit_vr(1:nl_soil) , &! (gN/m3/s) potential soil nitrification flux 
	
			pH                      , &! Input: soil water pH                                     
			lake_soilc   (1:nl_soil)   ! Input: total soil organic matter found in level (g C / m^3) (nl_soil)
				
  
		real(r8), INTENT(out) :: &            
			ch4_prod_depth    (1:nl_soil)         , &! Output: production of CH4 in each soil layer (nl_soil) (mol/m3/s)
			o2_decomp_depth   (1:nl_soil)            ! Output: O2 consumption during decomposition in each soil layer (nl_soil) (mol/m3/s)
		
		real(r8), INTENT(inout) :: &            
			layer_sat_lag   (1:nl_soil)           ! Lagged saturation status of soil layer in the unsaturated zone (1 = sat)
	
!-----------------------Local Variables------------------------------
		integer  :: j,s        ! indices
		real(r8) :: base_decomp      ! base rate (mol/m2/s)
		real(r8) :: partition_z
  
		real(r8) :: q10lake          ! For now, take to be the same as q10ch4 * 1.5.
		real(r8) :: q10lakebase      ! (K) base temperature for lake CH4 production
		real(r8) :: mino2lim         ! minimum anaerobic decomposition rate as a fraction of potential aerobic rate
  
		real(r8) :: q10ch4           ! additional Q10 for methane production ABOVE the soil decomposition temperature relationship  
		real(r8) :: q10ch4base       ! temperature at which the effective f_ch4 actually equals the constant f_ch4
		real(r8) :: f_ch4            ! ratio of CH4 production to total C mineralization
		! real(r8) :: rootlitfrac      ! Fraction of soil organic matter associated with roots
		real(r8) :: cnscalefactor    ! scale factor on CN decomposition for assigning methane flux
		real(r8) :: lake_decomp_fact ! Base decomposition rate (1/s) at 25C
  
		! added by Lei Meng to account for pH influence of CH4 production 
		real(r8) :: pHmax 
		real(r8) :: pHmin 
		real(r8) :: pH_fact_ch4      ! pH factor in methane production
  
		! Factors for methanogen temperature dependence being greater than soil aerobes
		real(r8) :: f_ch4_adj         ! Adjusted f_ch4
		real(r8) :: t_fact_ch4        ! Temperature factor calculated using additional Q10
  
		! O2 limitation on decomposition and methanogenesis
		real(r8) :: seasonalfin       ! finundated in excess of respiration-weighted annual average
		real(r8) :: oxinhib           ! inhibition of methane production by oxygen (m^3/mol)
  
		! For calculating column average (rootfrac(j)*rr(j))
		real(r8) :: rr_vr(1:nl_soil)  ! vertically resolved column-mean root respiration (g C/m^2/s)
  
		! 
		real(r8) :: sif   ! (unitless) ratio applied to sat. prod. to account for seasonal inundation
  
  
		character(len=32) :: subname='ch4_prod' ! subroutine name
  !-----------------------------------------------------------------------
  
		! Set transport parameters
		q10ch4           = params_inst%q10ch4
		q10ch4base       = params_inst%q10ch4base
		f_ch4            = params_inst%f_ch4
		cnscalefactor    = params_inst%cnscalefactor
		lake_decomp_fact = params_inst%lake_decomp_fact
		pHmax            = params_inst%pHmax
		pHmin            = params_inst%pHmin
		oxinhib          = params_inst%oxinhib
		q10lakebase      = params_inst%q10lakebase
  
		mino2lim         = params_inst%mino2lim
  
		q10lake = q10ch4 * 1.5_r8
  
		print*, 'q10lake',q10lake
		! PATCH loop to calculate vertically resolved column-averaged root respiration
		if (patchtype /= 4) then
		   rr_vr(:) = 0.0_r8
  
			if (lai > 0._r8) then
				do j=1,nl_soil
					rr_vr(j) = rr_vr(j) + rr*crootfr(j)
				end do
			end if
		  
		end if
		print*, rr_vr
  
		partition_z = 1._r8
		base_decomp = 0.0_r8
  
		! column loop to partition decomposition_rate into each soil layer
		do j=1,nl_soil
			print*, '               j==              ',j
			if (patchtype /= 4) then
				! Use soil heterotrophic respiration (based on Wania)
				base_decomp = (somhr+lithr) / catomw
				! Convert from gC to molC
				! Multiply base_decomp by factor accounting for lower carbon stock in seasonally inundated areas than
				! if it were inundated all year.
				! This is to reduce emissions in seasonally inundated zones, because the eq.
				! C-flux will be less than predicted by a non-O2-lim model
				if (sat == 1) then
					sif = 1._r8
					if (.not. anoxia) then
						if (annavg_finrw /= spval) then
						seasonalfin = max(finundated-annavg_finrw, 0._r8)
						if (seasonalfin > 0._r8) then
							sif = (annavg_finrw + mino2lim*seasonalfin) / finundated
							base_decomp = base_decomp * sif
						end if
						end if
					end if ! anoxia
				end if
	
				! For sensitivity studies
				base_decomp = base_decomp * cnscalefactor
	
			else !lake
  
				base_decomp = lake_decomp_fact * lake_soilc(j) * dz_soisno(j) * &
					q10lake**((t_soisno(j)-q10lakebase)/10._r8) / catomw
				! convert from g C to mol C
			end if
			print*, 'base_decomp',base_decomp
  
			! For all landunits, prevent production or oxygen consumption when soil is at or below freezing.
			! If using VERTSOILC, it is OK to use base_decomp as given because liquid water stress will limit decomp.
			if (t_soisno(j) <= tfrz .and. (patchtype == 4)) base_decomp = 0._r8
  
			! depth dependence of production either from rootfr or decomp model
			if (patchtype /= 4) then ! use default rootfr, averaged to the column level in the ch4 driver, or vert HR
				if ( (somhr + lithr) > 0._r8) then
					partition_z = hr_vr(j) * dz_soisno(j) / (somhr + lithr)
				else
					partition_z = 1._r8
				end if
			else ! lake
				partition_z = 1._r8
			endif
			print*, 'partition_z',partition_z
	
			! Adjust f_ch4 to account for the fact that methanogens may have a higher Q10 than aerobic decomposers.
			! Note this is crude and should ideally be applied to all anaerobic decomposition rather than just the
			! f_ch4.
			f_ch4_adj = 1.0_r8
			if (patchtype /= 4) then
				t_fact_ch4 = q10ch4**((t_soisno(j) - q10ch4base)/10._r8)
				! Adjust f_ch4 by the ratio
				f_ch4_adj = f_ch4 * t_fact_ch4
	
				! Remove CN nitrogen limitation, as methanogenesis is not N limited.
				! Also remove (low) moisture limitation
				if (ch4rmcnlim) then
					if (fphr(j) > 0._r8) then
						f_ch4_adj = f_ch4_adj / fphr(j)
					end if
				end if
	
			else ! lake
				f_ch4_adj = 0.5_r8 ! For lakes assume no redox limitation. Production only depends on temp, soil C, and
				! lifetime parameter.
			end if
			print*, 'f_ch4_adj',f_ch4_adj
	
			! If switched on, use pH factor for production based on spatial pH data defined in surface data.
			if ((patchtype /= 4) .and. usephfact )then 
				if (  pH >  pHmin .and.pH <  pHmax) then
					pH_fact_ch4 = 10._r8**(-0.2235_r8*pH*pH + 2.7727_r8*pH - 8.6_r8)
					! fitted function using data from Dunfield et al. 1993  
					! Strictly less than one, with optimum at 6.5
					! From Lei Meng
					f_ch4_adj = f_ch4_adj * pH_fact_ch4
				end if
			else
				! if no data, then no pH effects
			end if
	
			! Redox factor                                ??????????????????????????????????
			if (j > jwt) then ! Assume lag in decay of alternative electron acceptors vertically
				f_ch4_adj = f_ch4_adj * layer_sat_lag(j)
			end if
			! Alternative electron acceptors will be consumed first after soil is inundated.
	
			f_ch4_adj = min(f_ch4_adj, 0.5_r8)
			! Must be less than 0.5 because otherwise the actual implied aerobic respiration would be negative.
			! The total of aer. respiration + methanogenesis must remain equal to the SOMHR calculated in CN,
			! so that the NEE is sensible. Even perfectly anaerobic conditions with no alternative
			! electron acceptors would predict no more than 0.5 b/c some oxygen is present in organic matter.
			! e.g. 2CH2O --> CH4 + CO2.
			print*, 'f_ch4_adj',f_ch4_adj
	
	
			! Decomposition uses 1 mol O2 per mol CO2 produced (happens below WT also, to deplete O2 below WT)
			! o2_decomp_depth is the demand in the absense of O2 supply limitation, in addition to autotrophic respiration.
			! Competition will be done in ch4_prod
	
			o2_decomp_depth(j) = base_decomp * partition_z / dz_soisno(j)
			if (anoxia) then
				! Divide off o_scalar to use potential O2-unlimited HR to represent aerobe demand for oxygen competition
				if (patchtype /= 4) then
					if (o_scalar(j) > 0._r8) then
						o2_decomp_depth(j) = o2_decomp_depth(j) / o_scalar(j)
					end if
				end if
			end if ! anoxia
	
			! Add root respiration
			if (patchtype /= 4) then
				o2_decomp_depth(j) = o2_decomp_depth(j) + rr_vr(j)/catomw/dz_soisno(j) ! mol/m^3/s
				! g C/m2/s ! gC/mol O2 ! m
			end if
	
			! Add oxygen demand for nitrification
			if (use_nitrif_denitrif) then
				if (patchtype /= 4) then
					o2_decomp_depth(j) = o2_decomp_depth(j) + pot_f_nit_vr(j) * 2.0_r8/14.0_r8
					! g N/m^3/s           mol O2 / g N
				end if
			end if
			print*, 'o2_decomp_depth',j,'==',o2_decomp_depth(j)
	
			if (j  >  jwt) then ! Below the water table so anaerobic CH4 production can occur
				! partition decomposition to layer
				! turn into per volume-total by dz
				ch4_prod_depth(j) = f_ch4_adj * base_decomp * partition_z / dz_soisno(j)! [mol/m3-total/s]
			else ! Above the WT
				if (anoxicmicrosites) then
					ch4_prod_depth(j) = f_ch4_adj * base_decomp * partition_z / dz_soisno(j) &
						/ (1._r8 + oxinhib*conc_o2(j))
				else
					ch4_prod_depth(j) = 0._r8 ! [mol/m3 total/s]
				endif ! anoxicmicrosites
			endif ! WT
			print*, 'ch4_prod_depth',j,'==',ch4_prod_depth(j)
	
		end do ! nl_soil
	
	end subroutine ch4_prod

!-----------------------------------------------------------------------
	subroutine ch4_oxid (nl_soil,  jwt,  sat,  deltim,  z_soisno,  dz_soisno,  zi_soisno, &
		t_soisno,      smp,      vol_liq,    porsl,   conc_o2,   conc_ch4,              &
		ch4_oxid_depth,          o2_oxid_depth) 
!-----------------------DESCRIPTION------------------------------------------------------------------
! Oxidation is based on double Michaelis-Mentin kinetics, and is adjusted for low soil moisture.
! Oxidation will be limited by available oxygen in ch4_tran.
!-----------------------Argument---------- ----------------------------------------------------------
		integer , INTENT(in) :: &
			nl_soil                , &! number of soil layers
			jwt                    , &! index of the soil layer right above the water table (-) 
			sat                       ! 0 = unsaturated; 1 = saturated 

		real(r8), INTENT(in) :: &
			deltim                 , &! land model time step (sec)
			z_soisno (1:nl_soil)   , &! layer depth (m)
			dz_soisno(1:nl_soil)   , &! layer thickness (m)
			zi_soisno(0:nl_soil)   , &! interface level below a "z" level (m)

			t_soisno (1:nl_soil)   , &! soil temperature (Kelvin)
			smp      (1:nl_soil)   , &! soil matrix potential [mm]
			vol_liq  (1:nl_soil)   , &! liquid volumetric water content
			porsl    (1:nl_soil)   , &! volumetric soil water at saturation (porosity)

			conc_o2  (1:nl_soil)   , &! O2 conc in each soil layer (mol/m3) 
			conc_ch4 (1:nl_soil)      ! CH4 conc in each soil layer (mol/m3) 

		real(r8), INTENT(out) :: &
			ch4_oxid_depth (1:nl_soil)   , &! CH4 consumption rate via oxidation in each soil layer (mol/m3/s) 
			o2_oxid_depth  (1:nl_soil)      ! O2 consumption rate via oxidation in each soil layer (mol/m3/s) 
!-----------------------Local Variables------------------------------
		integer :: j                              ! indices
		real(r8):: t0                             ! Base temperature for Q10
		real(r8):: porevol                        ! air-filled volume ratio to total soil volume
		real(r8):: vol_liq_min                    ! vol_liq restricted to be below porsl
		real(r8):: conc_ch4_rel                   ! concentration with respect to water volume (mol/m^3 water)
		real(r8):: conc_o2_rel                    ! concentration with respect to water volume (mol/m^3 water)
		real(r8):: oxid_a                         ! Oxidation predicted by method A (temperature & enzyme limited) (mol CH4/m3/s)
		real(r8):: smp_fact                       ! factor for reduction based on soil moisture (unitless)
		real(r8):: k_h_cc, k_h_inv                ! see functions below for description
		real(r8):: k_m_eff                        ! effective k_m
		real(r8):: vmax_eff                       ! effective vmax

		! ch4 oxidation parameters
		real(r8) :: vmax_ch4_oxid                 ! oxidation rate constant (= 45.e-6_r8 * 1000._r8 / 3600._r8) [mol/m3-w/s];
		real(r8) :: k_m 			      ! Michaelis-Menten oxidation rate constant for CH4 concentration 
		real(r8) :: q10_ch4_oxid                   ! Q10 oxidation constant
		real(r8) :: smp_crit                      ! Critical soil moisture potential
		real(r8) :: k_m_o2                        ! Michaelis-Menten oxidation rate constant for O2 concentration
		real(r8) :: k_m_unsat                     ! Michaelis-Menten oxidation rate constant for CH4 concentration
		real(r8) :: vmax_oxid_unsat               ! (= 45.e-6_r8 * 1000._r8 / 3600._r8 / 10._r8) [mol/m3-w/s]   
!-----------------------------------------------------------------------
		! Set oxidation parameters
		vmax_ch4_oxid   = params_inst%vmax_ch4_oxid
		k_m             = params_inst%k_m
		q10_ch4_oxid    = params_inst%q10_ch4_oxid
		smp_crit        = params_inst%smp_crit
		k_m_o2          = params_inst%k_m_o2
		k_m_unsat       = params_inst%k_m_unsat
		vmax_oxid_unsat = params_inst%vmax_oxid_unsat

		t0 = tfrz + 12._r8 ! Walter, for Michigan site where the 45 M/h comes from

		! Loop to determine oxidation in each layer
		do j=1,nl_soil
			if (sat == 1 .or. j > jwt) then
				! Literature (e.g. Bender & Conrad, 1992) suggests lower k_m and vmax for high-CH4-affinity methanotrophs in
				! upland soils consuming ambient methane.
				k_m_eff = k_m
				vmax_eff = vmax_ch4_oxid
			else
				k_m_eff = k_m_unsat
				vmax_eff = vmax_oxid_unsat
			end if

			porevol = max(porsl(j) - vol_liq(j), 0._r8)
			vol_liq_min = min(porsl(j), vol_liq(j))
			if (j <= jwt .and. smp(j) < 0._r8) then
				smp_fact = exp(-smp(j)/smp_crit)
				! Schnell & King, 1996, Figure 3
			else
				smp_fact = 1._r8
			end if

			if (j  <=  jwt) then ! Above the water table
				k_h_inv = exp(-c_h_inv(1) * (1._r8 / t_soisno(j) - 1._r8 / kh_tbase) + log (kh_theta(1)))
				k_h_cc = t_soisno(j) / k_h_inv * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)]
				conc_ch4_rel = conc_ch4(j) / (vol_liq_min + porevol/k_h_cc)

				k_h_inv = exp(-c_h_inv(2) * (1._r8 / t_soisno(j) - 1._r8 / kh_tbase) + log (kh_theta(2)))
				k_h_cc = t_soisno(j) / k_h_inv * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)]
				conc_o2_rel  = conc_o2(j) / (vol_liq_min + porevol/k_h_cc)
			else
				conc_ch4_rel = conc_ch4(j) / porsl(j)
				conc_o2_rel  = conc_o2(j) / porsl(j)
			endif

			oxid_a              = vmax_eff     * vol_liq_min* conc_ch4_rel / (k_m_eff + conc_ch4_rel) &
							![mol/m3-t/s]         [mol/m3-w/s]    [m3-w/m3-t]     [mol/m3-w]    [mol/m3-w]  [mol/m3-w]
				* conc_o2_rel / (k_m_o2 + conc_o2_rel) &
				* q10_ch4_oxid ** ((t_soisno(j) - t0) / 10._r8) * smp_fact

			! For all landunits / levels, prevent oxidation if at or below freezing
			if (t_soisno(j) <= tfrz) oxid_a = 0._r8

			ch4_oxid_depth(j) = oxid_a
			o2_oxid_depth(j) = ch4_oxid_depth(j) * 2._r8

		end do
	end subroutine ch4_oxid


!--------------------------------------------------------------------------------------  
	subroutine ch4_aere (patchtype, nl_soil, jwt, sat, lai,     deltim, &
		z_soisno, dz_soisno,  zi_soisno,     t_soisno, vol_liq, porsl,  &
		rootfr,   rootr, etr, grnd_ch4_cond, c_atm,    annsum_npp,      &
		annavg_agnpp,    annavg_bgnpp, conc_o2, conc_ch4, ch4_prod_depth,&
		ch4_aere_depth, ch4_tran_depth, o2_aere_depth)
!-------------------------------------------------------------------------------------------------------------
! !DESCRIPTION:
! Arctic c3 grass (which is often present in fens) and all vegetation in inundated areas is assumed to have
! some root porosity. Currently, root porosity is allowed to be different for grasses & non-grasses.
! CH4 diffuses out and O2 diffuses into the soil.  CH4 is also lossed via transpiration, which is both
! included in the "aere" variables and output separately.  In practice this value is small.
! By default upland veg. has small 5% porosity but this can be switched to be equal to inundated porosity.

!-----------------------Argument-------------------------------------------------------------------------------
	integer, INTENT(in) :: &
		patchtype       ! land patch type (0=soil, 1=urban or built-up, 2=wetland,
						! 3=land ice, 4=land water bodies, 99=ocean

	integer , INTENT(in) :: &
		nl_soil                , &! number of soil layers
		jwt                    , &! index of the soil layer right above the water table (-) 
		sat                       ! 0 = unsatured, 1 = saturated 

	real(r8), INTENT(in) :: &
		lai                    , &! adjusted leaf area index for seasonal variation [-]

		deltim                 , &! land model time step (sec)
		z_soisno (1:nl_soil)   , &! layer depth (m)
		dz_soisno(1:nl_soil)   , &! layer thickness (m)
		zi_soisno(0:nl_soil)   , &! interface level below a "z" level (m)

		t_soisno (1:nl_soil)   , &! soil temperature (Kelvin)
		vol_liq  (1:nl_soil)   , &! liquid volumetric water content
		porsl    (1:nl_soil)   , &! volumetric soil water at saturation (porosity)
		rootfr   (1:nl_soil)   , &! fraction of roots in each soil layer
		
		rootr    (1:nl_soil)   , &! effective fraction of roots in each soil layer (SMS method only)
		! rootr here for effective per-layer transpiration, which may not be the same as rootfr
		etr                    , &! transpiration rate [mm/s]
		grnd_ch4_cond          , &! tracer conductance for boundary layer [m/s]

		c_atm(3)               , &! CH4, O2, CO2 atmospheric conc  (mol/m3)

		! These variables help us swap between big-leaf and fates boundary conditions
		annsum_npp             , &! annual sum NPP (molC m-2 s-1)
		annavg_agnpp           , &! annual avg aboveground NPP (gC/m2/s)
		annavg_bgnpp           , &! annual avg belowground NPP (gC/m2/s)

		! These variables help us swap between saturated and unsaturated boundary conditions
		conc_o2  (1:nl_soil)        , &! O2 conc in each soil layer (mol/m3) 
		conc_ch4 (1:nl_soil)        , &! CH4 conc in each soil layer (mol/m3) 
		ch4_prod_depth (1:nl_soil)     ! production of CH4 in each soil layer (mol/m3/s) 

	real(r8), INTENT(out) :: &
		ch4_aere_depth  (1:nl_soil)  , &! CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) 
		ch4_tran_depth  (1:nl_soil)  , &! CH4 loss rate via transpiration in each soil layer (mol/m3/s) 
		o2_aere_depth   (1:nl_soil)     ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s) 

!-----------------------Local Variables------------------------------
	integer  :: j              ! indices

	! ch4 aerenchyma parameters
	real(r8) :: tranloss(1:nl_soil)    ! loss due to transpiration (mol / m3 /s)
	real(r8) :: aere    (1:nl_soil) 
	real(r8) :: oxaere  (1:nl_soil)    ! (mol / m3 /s)

	real(r8) :: aeretran
	real(r8) :: poros_tiller

	! These pointers help us swap between saturated and unsaturated boundary conditions
	real(r8), parameter :: smallnumber = 1.e-12_r8   
!-----------------------------------------------------------------------

		! Initialize ch4_aere_depth
		do j=1,nl_soil
			ch4_aere_depth(j) = 0._r8
			ch4_tran_depth(j) = 0._r8
			o2_aere_depth(j) = 0._r8
		end do


		! point loop to partition aerenchyma flux into each soil layer
		if (patchtype /= 4) then

			poros_tiller = 0.3

			call SiteOxAere(nl_soil,  jwt,  sat,lai,    z_soisno, dz_soisno,  zi_soisno,  t_soisno,  &
			vol_liq,  porsl,  rootfr,   rootr,  poros_tiller,grnd_ch4_cond, etr,   &
			annsum_npp, annavg_agnpp,   annavg_bgnpp,  c_atm,      conc_o2, conc_ch4,        &
			tranloss, aere,   oxaere)

			do j = 1,nl_soil
				! Impose limitation based on available methane during timestep
				! By imposing the limitation here, don't allow aerenchyma access to methane from other Patches.
				aeretran = min(aere(j)+tranloss(j), conc_ch4(j)/deltim + ch4_prod_depth(j))
				ch4_aere_depth (j) = ch4_aere_depth(j) + aeretran
				ch4_tran_depth (j) = ch4_tran_depth(j) + min(tranloss(j), aeretran)
				o2_aere_depth  (j) = o2_aere_depth (j) + oxaere(j)
			end do ! over levels
		end if ! not lake

	end subroutine ch4_aere


!--------------------------------------------------------------------------------------  
	subroutine SiteOxAere(nl_soil,  jwt,  sat, lai,    z_soisno, dz_soisno,  zi_soisno,  t_soisno,  &
		vol_liq,  porsl,  rootfr,   rootr,  poros_tiller,grnd_ch4_cond, etr,   &
		annsum_npp, annavg_agnpp,   annavg_bgnpp,  c_atm,      conc_o2, conc_ch4,        &
		tranloss, aere,   oxaere)
!-----------------------DESCRIPTION--------------------------------------
! Site(column) level fluxes for O2 gain rate via
! aerenchyma and ch4 losss rates from transpiration
!-----------------------Argument---------- ------------------------------
		integer , INTENT(in) :: &
			nl_soil                , &! number of soil layers
			jwt                    , &! index of the soil layer right above the water table (-) 
			sat                       ! 0 = unsatured, 1 = saturated 

		real(r8), INTENT(in) :: &
			lai                    , &! leaf area index [m2/m2]

			z_soisno (1:nl_soil)   , &! layer depth (m)
			dz_soisno(1:nl_soil)   , &! layer thickness (m)
			zi_soisno(0:nl_soil)   , &! interface level below a "z" level (m)

			t_soisno (1:nl_soil)   , &! soil temperature (Kelvin)
			vol_liq  (1:nl_soil)   , &! liquid volumetric water content
			porsl    (1:nl_soil)   , &! volumetric soil water at saturation (porosity)
			rootfr   (1:nl_soil)   , &! fraction of roots in each soil layer
			rootr    (1:nl_soil)   , &! root resistance of a layer, all layers add to 1
			grnd_ch4_cond          , &! tracer conductance for boundary layer [m/s] 
			etr                    , &! transpiration rate [mm/s]



			annsum_npp             , &! annual sum NPP (molC m-2 yr-1)
			annavg_agnpp           , &! annual average aboveground NPP (gC/m2/s)
			annavg_bgnpp           , &! annual average belowground NPP (gC/m2/s)

			c_atm(3)               , &! CH4, O2, CO2 atmospheric conc  (mol/m3)
			conc_o2  (1:nl_soil)   , &! O2 conc in each soil layer (mol/m3) 
			conc_ch4 (1:nl_soil)      ! CH4 conc in each soil layer (mol/m3) 

		real(r8), INTENT(out) :: &
			tranloss        (1:nl_soil)  , &! CH4 in soil water tran rate via plant transpiration in each soil layer (mol/m3/s) 
			aere            (1:nl_soil)  , &! CH4 tran rate via aerenchyma in each soil layer (mol/m3/s) 
			oxaere          (1:nl_soil)     ! O2 gain rate via aerenchyma in each soil layer (mol/m3/s) 

		real(r8), INTENT(inout) :: &
			poros_tiller           

!-----------------------Local Variables------------------------------         
		integer  :: j
		real(r8) :: oxdiffus
		real(r8) :: area_tiller ! cross-sectional area of tillers (m^2/m^2)
		real(r8) :: diffus_aere ! gas diffusivity through aerenchyma (m^2/s)
		real(r8) :: m_tiller 
		real(r8) :: n_tiller 
		real(r8) :: vol_liq_min
		real(r8) :: k_h_cc, k_h_inv
		real(r8) :: anpp, nppratio
		real(r8) :: conc_ch4_wat
		real(r8) :: aerecond    ! aerenchyma conductance (m/s)
		real(r8), parameter :: smallnumber = 1.e-12_r8

		real(r8) :: scale_factor_aere
		real(r8) :: unsat_aere_ratio
		real(r8) :: porosmin
		real(r8) :: rob
!-----------------------------------------------------------------------

		scale_factor_aere = params_inst%scale_factor_aere
		unsat_aere_ratio = params_inst%unsat_aere_ratio
		porosmin = params_inst%porosmin

		rob = params_inst%rob

		! This parameter is poorly constrained and should be done on a patch-specific basis...
		diffus_aere = d_con_g(1,1)*1.e-4_r8  ! for CH4: m^2/s

		do j=1,nl_soil
			! Calculate transpiration loss
			if (transpirationloss .and. lai > 0) then
				! Calculate water concentration
				vol_liq_min = min(porsl(j), vol_liq(j))
				k_h_inv = exp(-c_h_inv(1) * (1._r8 / t_soisno(j) - 1._r8 / kh_tbase) + log (kh_theta(1)))
				k_h_cc = t_soisno(j) / k_h_inv * rgasLatm
				conc_ch4_wat = conc_ch4(j) / ( (porsl(j)-vol_liq_min)/k_h_cc + vol_liq_min)

				tranloss(j) = conc_ch4_wat * rootr(j)*etr / dz_soisno(j) / 1000._r8
				! mol/m3/s    mol/m3                  mm / s  m           mm/m
				! Use rootr here for effective per-layer transpiration, which may not be the same as rootfr
				tranloss(j) = max(tranloss(j), 0._r8) ! in case transpiration is pathological
			else
				tranloss(j) = 0._r8
			end if

			! Calculate aerenchyma diffusion
			if (j > jwt .and. t_soisno(j) > tfrz .and. lai > 0) then
				! Attn EK: This calculation of aerenchyma properties is very uncertain. Let's check in once all
				! the new components are in; if there is any tuning to be done to get a realistic global flux,
				! this would probably be the place.  We will have to document clearly in the Tech Note
				! any major changes from the Riley et al. 2011 version. (There are a few other minor ones.)

				anpp = annsum_npp ! g C / m^2/yr
				anpp = max(anpp, 0._r8) ! NPP can be negative b/c of consumption of storage pools

				if (annavg_agnpp /= spval .and. annavg_bgnpp /= spval .and. &
					annavg_agnpp > 0._r8 .and. annavg_bgnpp > 0._r8) then
					nppratio = annavg_bgnpp / (annavg_agnpp + annavg_bgnpp)
				else
					nppratio = 0.5_r8
				end if

				! Estimate area of tillers (see Wania thesis)
				!m_tiller = anpp * r_leaf_root * lai ! (4.17 Wania)
				!m_tiller = 600._r8 * 0.5_r8 * 2._r8  ! used to be 300
				! Note: this calculation is based on Arctic graminoids, and should be refined for woody plants, if not
				! done on a PFT-specific basis.

				m_tiller = anpp * nppratio * 4._r8  !replace the elai(p) by constant 4 (by Xiyan Xu, 05/2016)


				n_tiller = m_tiller / 0.22_r8

				if (sat == 0) then ! unsaturate
					poros_tiller = poros_tiller * unsat_aere_ratio
				end if

				poros_tiller = max(poros_tiller, porosmin)

				area_tiller = scale_factor_aere * n_tiller * poros_tiller * PI * 2.9e-3_r8**2._r8 ! (m2/m2)

				k_h_inv = exp(-c_h_inv(1) * (1._r8 / t_soisno(j) - 1._r8 / kh_tbase) + log (kh_theta(1))) ! (4.12) Wania (L atm/mol)
				k_h_cc = t_soisno(j) / k_h_inv * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)]
				aerecond = area_tiller * rootfr(j) * diffus_aere / (z_soisno(j)*rob)
				! Add in boundary layer resistance
				aerecond = 1._r8 / (1._r8/(aerecond+smallnumber) + 1._r8/(grnd_ch4_cond+smallnumber))

				aere(j) = aerecond * (conc_ch4(j)/porsl(j)/k_h_cc - c_atm(1)) / dz_soisno(j) ![mol/m3-total/s]
				!ZS: Added porsl & Henry's const.
				aere(j) = max(aere(j), 0._r8) ! prevent backwards diffusion

				! Do oxygen diffusion into layer
				k_h_inv = exp(-c_h_inv(2) * (1._r8 / t_soisno(j) - 1._r8 / kh_tbase) + log (kh_theta(2)))
				k_h_cc = t_soisno(j) / k_h_inv * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)]
				oxdiffus = diffus_aere * d_con_g(2,1) / d_con_g(1,1) ! adjust for O2:CH4 molecular diffusion
				aerecond = area_tiller * rootfr(j) * oxdiffus / (z_soisno(j)*rob)
				aerecond = 1._r8 / (1._r8/(aerecond+smallnumber) + 1._r8/(grnd_ch4_cond+smallnumber))
				oxaere(j) = -aerecond *(conc_o2(j)/porsl(j)/k_h_cc - c_atm(2)) / dz_soisno(j) ![mol/m3-total/s]
				oxaere(j) = max(oxaere(j), 0._r8)
				! Diffusion in is positive; prevent backwards diffusion
				if ( .not. use_aereoxid_prog ) then ! fixed aere oxid proportion; will be done in ch4_tran
					oxaere(j) = 0._r8
				end if
			else
				aere(j) = 0._r8
				oxaere(j) = 0._r8
			end if ! veg type, below water table, & above freezing
		end do

  	end subroutine SiteOxAere

	
!-----------------------------------------------------------------------
	subroutine ch4_ebul (patchtype, nl_soil, jwt, sat, deltim, &
		z_soisno, dz_soisno, zi_soisno, lakedepth, forc_pbot,&
		t_soisno, lake_icefrac, porsl, wdsrf, conc_ch4,&
		ch4_ebul_depth)

!-----------------------DESCRIPTION------------------------------
! Bubbling is based on temperature & pressure dependent solubility (k_h_cc), 
! with assumed proportion of bubbles
! which are CH4, and assumed early nucleation at vgc_max sat (Wania).
! Bubbles are released to the water table surface in ch4_tran.  
!-----------------------Argument---------- ------------------------------
		integer, INTENT(in) :: &
			patchtype           ! land patch type (0=soil, 1=urban or built-up, 2=wetland,
									  ! 3=land ice, 4=land water bodies, 99=ocean

		integer , INTENT(in) :: &
			nl_soil                 , &! number of soil layers
			jwt                     , &! index of the soil layer right above the water table (-) 
         	sat                        ! 0 = unsaturated; 1 = saturated 

		real(r8), INTENT(in) :: &
			deltim                     , &! land model time step (sec)
			z_soisno (1:nl_soil)       , &! layer depth (m)
			dz_soisno(1:nl_soil)       , &! layer thickness (m)
			zi_soisno(0:nl_soil)       , &! interface level below a "z" level (m)

			lakedepth                  , &! lake depth
			forc_pbot                  , &! atm bottom level pressure (or reference height) (pa)
			t_soisno (1:nl_soil)       , &! soil temperature (Kelvin)
			lake_icefrac (1:nl_soil)   , &! lake mass fraction of lake layer that is frozen
			porsl    (1:nl_soil)       , &! volumetric soil water at saturation (porosity)
			! vol_liq  (1:nl_soil)       , &! liquid volumetric water content
			wdsrf                      , &! depth of surface water [mm]
			conc_ch4       (1:nl_soil)    ! Output: CH4 conc in each soil layer (mol/m3) 

		real(r8), INTENT(out) :: &
			ch4_ebul_depth (1:nl_soil)   ! CH4 loss rate via ebullition in each soil layer (mol/m3/s)

!-----------------------Local Variables------------------------------
		integer :: j      ! indices

		real(r8) :: vgc     ! volumetric CH4 content (m3 CH4/m3 pore air)
		real(r8) :: vgc_min ! minimum aqueous CH4 content when ebullition ceases
		real(r8) :: k_h_inv ! 
		real(r8) :: k_h     ! 
		real(r8) :: k_h_cc  ! 
		real(r8) :: pressure! sum atmospheric and hydrostatic pressure
		real(r8) :: bubble_f! CH4 content in gas bubbles (Kellner et al. 2006)
		real(r8) :: ebul_timescale

		real(r8) :: vgc_max   ! ratio of saturation pressure triggering ebullition
!-----------------------------------------------------------------------
		! Set ebullution parameter
		vgc_max = params_inst%vgc_max

		bubble_f = 0.57_r8 ! CH4 content in gas bubbles (Kellner et al. 2006)
		vgc_min = vgc_max
		ebul_timescale = deltim ! Allow fast bubbling

		! IF (.not. DEF_USE_VARIABLY_SATURATED_FLOW) THEN
		! 	wdsrf = 0
		! END IF

		! column loop to estimate ebullition CH4 flux from each soil layer
		do j=1,nl_soil
			if (j  >  jwt .and. t_soisno(j) > tfrz) then ! Ebullition occurs only below the water table
				k_h_inv = exp(-c_h_inv(1) * (1._r8 / t_soisno(j) - 1._r8 / kh_tbase) + log (kh_theta(1))) ! (4.12 Wania) (atm.L/mol)
				k_h = 1._r8 / k_h_inv ! (mol/L.atm)
				k_h_cc = t_soisno(j) * k_h * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)] 

				if (patchtype/=4) then
					pressure = forc_pbot + denh2o * grav * (z_soisno(j)-zi_soisno(jwt)) ! (Pa)
					if (sat == 1 ) then ! Add ponding pressure head
						pressure = pressure + denh2o * grav * wdsrf/1000._r8
						! mm     / mm/m
					end if
				else
					pressure = forc_pbot + denh2o * grav * (z_soisno(j) + lakedepth)
				end if

				! Compare partial pressure to ambient pressure.
				vgc = conc_ch4(j) / porsl(j) / k_h_cc * rgasm * t_soisno(j) / pressure
				! [mol/m3t]      [m3w/m3t]   [m3g/m3w]  [Pa/(mol/m3g)]          [Pa]

				if (vgc > vgc_max * bubble_f) then ! If greater than max value, remove amount down to vgc_min
					ch4_ebul_depth(j) = (vgc - vgc_min * bubble_f) * conc_ch4(j) / ebul_timescale
					! [mol/m3t/s]                                       [mol/m3t]         [s]
				else
					ch4_ebul_depth(j) = 0._r8
				endif

			else ! above the water table or freezing
				ch4_ebul_depth(j) = 0._r8
			endif ! below the water table and not freezing

			! Prevent ebullition from reaching the surface for frozen lakes
			if (patchtype==4 .and. lake_icefrac(1) > 0.1_r8) ch4_ebul_depth(j) = 0._r8
		end do ! j

	end subroutine ch4_ebul

!-----------------------------------------------------------------------
	subroutine ch4_tran (patchtype, &
		lb, nl_soil, snl, maxsnl,jwt, sat,&
		lon, lat, deltim, z_soisno, dz_soisno, zi_soisno,  t_soisno, t_grnd, &
		vol_liq, porsl, wliq_soisno, wice_soisno, wdsrf, bsw, c_atm, ch4_prod_depth, o2_aere_depth,&
		cellorg,t_h2osfc, organic_max, &
		o2stress, ch4stress, ch4_surf_aere, ch4_surf_ebul, ch4_surf_diff, ch4_ebul_total, &
		ch4_oxid_depth, ch4_aere_depth, ch4_ebul_depth, &
		grnd_ch4_cond, o2_oxid_depth, o2_decomp_depth, conc_o2, conc_ch4 )
!------------------------------------------------------------------------------
! !DESCRIPTION:
! Solves the reaction & diffusion equation for the timestep.  First "competition" between processes for
! CH4 & O2 demand is done.  Then concentrations are apportioned into gas & liquid fractions; only the gas
! fraction is considered for diffusion in unsat.  Snow and lake water resistance to diffusion is added as
! a bulk term in the ground conductance (which is really a surface layer conductance), but concentrations
! are not tracked and oxidation is not allowed inside snow and lake water.
! Diffusivity is set based on soil texture and organic matter fraction. A Crank-Nicholson solution is used.
! Then CH4 diffusive flux is calculated and consistency is checked.

!-----------------------Argument----------------------------------------
		integer, INTENT(in) :: &
			patchtype        	! land patch type (0=soil, 1=urban or built-up, 2=wetland,
										! 3=land ice, 4=land water bodies, 99=ocean
			! istep             , &! the i time step


		integer , INTENT(in) :: &
			lb               , &! lower bound of array (snl+1)
			nl_soil          , &! upper bound of array (10)
			snl				  , &!  number of snow layers     (-5~-1)
			maxsnl			  , &!  max number of snow layers (-5)
			jwt              , &! index of the soil layer right above the water table (-) 
			sat                 ! 0 = unsaturated; 1 = saturated 

		real(r8), INTENT(in) :: &
			lon   	   				   , &! logitude 
			lat     	   					, &! latitude 

			deltim                  	, &! land model time step (sec)
			z_soisno (lb:nl_soil)   	, &! layer depth (m)
			dz_soisno(lb:nl_soil)   	, &! layer thickness (m)
			zi_soisno(lb-1:nl_soil)   	, &! interface level below a "z" level (m)


			t_soisno (lb:nl_soil)    	, &! soil temperature (Kelvin)
			t_grnd                 		, &! ground surface temperature [k]

			vol_liq  (1:nl_soil)   		, &! liquid volumetric water content
			porsl    (1:nl_soil)   		, &! volumetric soil water at saturation (porosity)
			wliq_soisno(lb:nl_soil)		, &! liquid water in layers [kg/m2]
			wice_soisno(lb:nl_soil) 	, &! ice lens in layers [kg/m2]
			wdsrf                  		, &! depth of surface water [mm]
			bsw      (1:nl_soil)   		, &! Clapp and Hornberger "b" (nlevgrnd)             

			c_atm(3)               		, &! CH4, O2, CO2 atmospheric conc  (mol/m3)

			ch4_prod_depth    (1:nl_soil)  , &! production of CH4 in each soil layer (mol/m3/s) 
			o2_aere_depth     (1:nl_soil)  , &! O2 gain rate via aerenchyma in each soil layer (mol/m3/s) 


			cellorg  (1:nl_soil)   		, &! column 3D org (kg/m^3 organic matter)
			t_h2osfc               		, &! surface water temperature               
			organic_max               		! organic matter content (kg m-3) where soil is assumed to act like peat

		real(r8), INTENT(out) :: &
			o2stress          (1:nl_soil)  , &! Output: Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs
			ch4stress         (1:nl_soil)  , &! Output: Ratio of methane available to the total per-timestep methane sinks 
			ch4_surf_aere                  , &! Output: Total column CH4 aerenchyma (mol/m2/s)
			ch4_surf_ebul                  , &! Output: CH4 ebullition to atmosphere (mol/m2/s)
			ch4_surf_diff                  , &! Output: CH4 surface flux (mol/m2/s)
			ch4_ebul_total                    ! Output: Total column CH4 ebullition (mol/m2/s)

		real(r8), INTENT(inout) :: &
			ch4_oxid_depth    (1:nl_soil)  , &! InOut: CH4 consumption rate via oxidation in each soil layer (mol/m3/s) 
			ch4_aere_depth    (1:nl_soil)  , &! InOut: CH4 loss rate via aerenchyma in each soil layer (mol/m3/s) 
			ch4_ebul_depth    (1:nl_soil)  , &! InOut: CH4 loss rate via ebullition in each soil layer (mol/m3/s)
			o2_oxid_depth     (1:nl_soil)  , &! InOut: O2 loss rate via ebullition in each soil layer (mol/m3/s) 
			o2_decomp_depth   (1:nl_soil)  , &! InOut: O2 consumption during decomposition in each soil layer (mol/m3/s)

			grnd_ch4_cond                  , &! InOut: tracer conductance for boundary layer [m/s]  
			conc_o2           (1:nl_soil)  , &! InOut: O2 conc in each soil layer (mol/m3) 
			conc_ch4          (1:nl_soil)     ! InOut: CH4 conc in each soil layer (mol/m3) 
  !-----------------------Local Variables------------------------------
		integer :: j,s,i			                                               ! indices
		integer :: jtop                                                        ! top level at each column
		real(r8) :: at (0:nl_soil)                     ! "a" vector for tridiagonal matrix
		real(r8) :: bt (0:nl_soil)                     ! "b" vector for tridiagonal matrix
		real(r8) :: ct (0:nl_soil)                     ! "c" vector for tridiagonal matrix
		real(r8) :: rt (0:nl_soil)                     ! "r" vector for tridiagonal solution
		real(r8) :: f_a                                                        ! air-filled fraction of available pore space
		real(r8) :: diffus (0:nl_soil)                 ! diffusivity (m2/s)
		real(r8) :: k_h_inv                                                    ! 1/Henry's Law Constant in Latm/mol
		real(r8) :: k_h_cc(0:nl_soil,ngases)           ! ratio of mol/m3 in liquid to mol/m3 in gas
		real(r8) :: dzj                                                        ! 
		real(r8) :: dp1_zp1 (0:nl_soil)                ! diffusivity/delta_z for next j
		real(r8) :: dm1_zm1 (0:nl_soil)                ! diffusivity/delta_z for previous j
		real(r8) :: t_soisno_c                                                 ! soil temperature   (maxsnl+1:nl_soil)
		real(r8) :: eps                                                        ! either epsilon_a or epsilon_w, depending on where in soil, wrt WT
		real(r8) :: deficit                                                    ! mol CH4 /m^2 that must be subtracted from diffusive flux to atm. to make up
		! for keeping concentrations always above zero
		real(r8) :: conc_ch4_bef(1:nl_soil)            ! concentration at the beginning of the timestep
		real(r8) :: errch4                            ! Error (Mol CH4 /m^2) [+ = too much CH4]
		real(r8) :: conc_ch4_rel(0:nl_soil)            ! Concentration per volume of air or water
		real(r8) :: conc_o2_rel(0:nl_soil)             ! Concentration per volume of air or water
		real(r8) :: conc_ch4_rel_old(0:nl_soil)        ! Concentration during last Crank-Nich. loop
		real(r8) :: vol_liq_min(1:nl_soil)          ! vol_liq restricted to be <= porsl
		real(r8), parameter :: smallnumber = 1.e-12_r8
		real(r8) :: snowdiff                                                   ! snow diffusivity (m^2/s)
		real(r8) :: snowres                           ! Cumulative Snow resistance (s/m). Also includes
		real(r8) :: pondres                                                    ! Additional resistance from ponding, up to pondmx water on top of top soil layer (s/m)
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
		real(r8) :: liqfrac(1:nl_soil)

		real(r8) :: satpow                                                     ! exponent on porsl for saturated soil solute diffusion
		real(r8) :: scale_factor_gasdiff                                       ! For sensitivity tests; convection would allow this to be > 1
		real(r8) :: scale_factor_liqdiff                                       ! For sensitivity tests; convection would allow this to be > 1
		real(r8) :: capthick                                                   ! (mm) min thickness before assuming wdsrf is impermeable
		real(r8) :: aereoxid                                                   ! fraction of methane flux entering aerenchyma rhizosphere 

		real(r8) :: om_frac_sf                                                 ! Scale factor for organic matter fraction (unitless)
  !-----------------------------------------------------------------------
  
		! Set transport parameters
		satpow               = params_inst%satpow
		scale_factor_gasdiff = params_inst%scale_factor_gasdiff
		scale_factor_liqdiff = params_inst%scale_factor_liqdiff
		capthick             = params_inst%capthick 
		aereoxid             = params_inst%aereoxid
  
		om_frac_sf           = params_inst%om_frac_sf
		! Perform competition for oxygen and methane in each soil layer if demands over the course of the timestep
		! exceed that available. Assign to each process in proportion to the quantity demanded in the absense of
		! the limitation.
		do j = 1,nl_soil
			o2demand = o2_decomp_depth(j) + o2_oxid_depth(j) ! o2_decomp_depth includes autotrophic root respiration
			if (o2demand > 0._r8) then
				if ( (conc_o2(j) / deltim + o2_aere_depth(j)) > o2demand )then
					o2stress(j) = 1._r8
				else
					o2stress(j) = (conc_o2(j) / deltim + o2_aere_depth(j)) / o2demand
				end if
			else
				o2stress(j) = 1._r8
			end if
	
			ch4demand = ch4_oxid_depth(j) + ch4_aere_depth(j) + ch4_ebul_depth(j)
			if (ch4demand > 0._r8) then
				ch4stress(j) = min((conc_ch4(j) / deltim + ch4_prod_depth(j)) / ch4demand, 1._r8)
			else
				ch4stress(j) = 1._r8
			end if
	
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
					o2_oxid_depth(j) = o2_oxid_depth(j) * ch4stress(j)
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
		end do !j
  
  
		! Accumulate ebullition to place in first layer above water table, or directly to atmosphere
		do j = 1,nl_soil
			if (j == 1) ch4_ebul_total = 0._r8
			ch4_ebul_total = ch4_ebul_total + ch4_ebul_depth(j) * dz_soisno(j)
		end do
  
		! Set the Henry's Law coefficients
		do j = 0,nl_soil
			do s=1,2         
				if (j == 0) then
					k_h_inv = exp(-c_h_inv(s) * (1._r8 / t_grnd - 1._r8 / kh_tbase) + log (kh_theta(s)))
					! (4.12) Wania (L atm/mol)
					k_h_cc(j,s) = t_grnd / k_h_inv * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)]
				else
					k_h_inv = exp(-c_h_inv(s) * (1._r8 / t_soisno(j) - 1._r8 / kh_tbase) + log (kh_theta(s)))
					! (4.12) Wania (L atm/mol)
					k_h_cc(j,s) = t_soisno(j) / k_h_inv * rgasLatm ! (4.21) Wania [(mol/m3w) / (mol/m3g)]
				end if
			end do
		end do
  
  
		! Set the source term for each species (no need to do j=0, since epsilon_t and source not used there)
		! Note that because of the semi-implicit diffusion and the 30 min timestep combined with explicit
		! sources, occasionally negative concentration will result. In this case it is brought to zero and the
		! surface flux is adjusted to conserve. This results in some inaccuracy as compared to a shorter timestep
		! or iterative solution.
		do j = 1,nl_soil
	
			if ( .not. use_aereoxid_prog ) then
				! First remove the CH4 oxidation that occurs at the base of root tissues (aere), and add to oxidation
				ch4_oxid_depth(j) = ch4_oxid_depth(j) + aereoxid * ch4_aere_depth(j)
				ch4_aere_depth(j) = ch4_aere_depth(j) - aereoxid * ch4_aere_depth(j)
			end if ! else oxygen is allowed to diffuse in via aerenchyma
	
			source(j,1) = ch4_prod_depth(j) - ch4_oxid_depth(j) - &
					ch4_aere_depth(j) - ch4_ebul_depth(j) ! [mol/m3-total/s]
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
	
			conc_ch4_bef(j) = conc_ch4(j) !For Balance Check
		enddo ! j
  
  
		! Accumulate aerenchyma to add directly to atmospheric flux
		do j = 1,nl_soil
			if (j==1) ch4_surf_aere = 0._r8
			ch4_surf_aere = ch4_surf_aere + ch4_aere_depth(j) * dz_soisno(j)
		enddo
  
		! Add in ebullition to source at depth just above WT
		if (jwt /= 0) then
		   	source(jwt,1) = source(jwt,1) + ch4_ebul_total/dz_soisno(jwt)
		endif
  
		! Calculate concentration relative to m^3 of air or water: needed for the diffusion
		do j = 0,nl_soil
			if (j == 0) then
				conc_ch4_rel(j) = c_atm(1)
				conc_o2_rel(j)  = c_atm(2)
			else
				vol_liq_min(j) = min(porsl(j), vol_liq(j))
				liqfrac(j) = 1._r8
				if (ch4frzout) then
					liqfrac(j) = max(0.05_r8, (wliq_soisno(j)/denh2o+smallnumber)/ &
						(wliq_soisno(j)/denh2o+wice_soisno(j)/denice+smallnumber))
				else
					liqfrac(j) = 1._r8
				end if
				if (j <= jwt) then  ! Above the WT
					do s=1,2
						epsilon_t(j,s) = porsl(j)- (1._r8-k_h_cc(j,s))*vol_liq_min(j)*liqfrac(j)
					end do
					! Partition between the liquid and gas phases. The gas phase will drive the diffusion.
				else ! Below the WT
					do s=1,2
						epsilon_t(j,s) = porsl(j)*liqfrac(j)
					end do
				end if
				conc_ch4_rel(j) = conc_ch4(j)/epsilon_t(j,1)
				conc_o2_rel(j)  = conc_o2(j) /epsilon_t(j,2)
			end if
		end do
  
  
		! Loop over species
		do s = 1, 2 ! 1=CH4; 2=O2; 3=CO2
			! Adjust the grnd_ch4_cond to keep it positive, and add the snow resistance & pond resistance
			do j = maxsnl + 1,0
				if (j == maxsnl + 1) then
					if (grnd_ch4_cond < smallnumber .and. s==1) grnd_ch4_cond = smallnumber
					! Needed to prevent overflow when ground is frozen, e.g. for lakes
					snowres = 0._r8
				end if
	
				! Add snow resistance
				if (j >= snl + 1) then
					t_soisno_c = t_soisno(j) - tfrz
					icefrac = wice_soisno(j)/denice/dz_soisno(j)
					waterfrac = wliq_soisno(j)/denh2o/dz_soisno(j)
					airfrac = max(1._r8 - icefrac - waterfrac, 0._r8)
					! Calculate snow diffusivity
					if (airfrac > 0.05_r8) then
						f_a = airfrac / (airfrac + waterfrac)
						eps = airfrac ! Air-filled fraction of total snow volume
						! Use Millington-Quirk Expression, as hydraulic properties (bsw) not available
						snowdiff = (d_con_g(s,1) + d_con_g(s,2)*t_soisno_c) * 1.e-4_r8 * &
								f_a**(10._r8/3._r8) / (airfrac+waterfrac)**2 &
								* scale_factor_gasdiff
					else !solute diffusion in water only
						eps = waterfrac  ! Water-filled fraction of total soil volume
						snowdiff = eps**satpow * (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
								* scale_factor_liqdiff
					end if
					snowdiff = max(snowdiff, smallnumber)
					snowres = snowres + dz_soisno(j)/snowdiff
				end if
	
				if (j == 0) then ! final loop
					! Add pond resistance
					pondres = 0._r8

					! First old pond formulation up to pondmx
					if (patchtype /= 4 .and. snl == 0 .and. vol_liq(1) > porsl(1)) then
						t_soisno_c = t_soisno(1) - tfrz
						if (t_soisno(1) <= tfrz) then
							ponddiff = (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
									* (wliq_soisno(1)/denh2o+smallnumber)/ &
									(wliq_soisno(1)/denh2o+wice_soisno(1)/denice+smallnumber) &
									* scale_factor_liqdiff
						else ! Unfrozen
							ponddiff = (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
									* scale_factor_liqdiff
						end if
						pondz = dz_soisno(1) * (vol_liq(1) - porsl(1))
						pondres = pondz / ponddiff
					end if

					! Now add new wdsrf form
					if (patchtype /= 4 .and. sat == 1) then
						if (t_h2osfc >= tfrz) then
							t_soisno_c = t_h2osfc - tfrz
							ponddiff = (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
									* scale_factor_liqdiff
							pondz = wdsrf / 1000._r8 ! Assume all wdsrf corresponds to sat area
							! mm      /  mm/m
							pondres = pondres + pondz / ponddiff
						else if (wdsrf > capthick) then
							! assume surface ice is impermeable
							pondres = 1/smallnumber
						end if
					end if

					spec_grnd_cond(s) = 1._r8/(1._r8/grnd_ch4_cond + snowres + pondres)
				end if
			end do ! j
	
			! Determine gas diffusion and fraction of open pore (f_a)
			do j = 1,nl_soil
				t_soisno_c = t_soisno(j) - tfrz
	
				if (j <= jwt) then  ! Above the WT
					f_a = 1._r8 - vol_liq_min(j) / porsl(j)
					! Provisionally calculate diffusivity as linear combination of the Millington-Quirk 
					! expression in Wania (for peat) & Moldrup (for mineral soil)
					eps =  porsl(j)-vol_liq_min(j) ! Air-filled fraction of total soil volume
					if (organic_max > 0._r8) then
						om_frac = min(om_frac_sf*cellorg(j)/organic_max, 1._r8)
						! Use first power, not square as in iniTimeConst
					else
						om_frac = 1._r8
					end if
					diffus (j) = (d_con_g(s,1) + d_con_g(s,2)*t_soisno_c) * 1.e-4_r8 * &
								(om_frac * f_a**(10._r8/3._r8) / porsl(j)**2._r8 + &
								(1._r8-om_frac) * eps**2._r8 * f_a**(3._r8 / bsw(j)) ) &
								* scale_factor_gasdiff
				else ! Below the WT use saturated diffusivity and only water in epsilon_t
					! Note the following is not currently corrected for the effect on diffusivity of excess ice in soil under
					! lakes (which is currently experimental only).
					eps = porsl(j)  ! Water-filled fraction of total soil volume
					diffus (j) = eps**satpow * (d_con_w(s,1) + d_con_w(s,2)*t_soisno_c + d_con_w(s,3)*t_soisno_c**2) * 1.e-9_r8 &
						* scale_factor_liqdiff
					if (t_soisno(j)<=tfrz) then
						diffus(j) = diffus(j)*(wliq_soisno(j)/denh2o+smallnumber)/ &
							(wliq_soisno(j)/denh2o+wice_soisno(j)/denice+smallnumber)
					end if
				end if ! Above/below the WT
				diffus(j) = max(diffus(j), smallnumber) ! Prevent overflow
	
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
					!ch4_surf_ebul = ch4_ebul_total ! [mol/m2/s]
				else ! WT at the surface; i.e., jwt==0
					ch4_surf_diff = dm1_zm1(1) * ( (conc_ch4_rel(1)+conc_ch4_rel_old(1))/2._r8 &
						- c_atm(s)*k_h_cc(0,s)) ! [mol/m2/s]
					! atmospheric concentration gets mult. by k_h_cc as above
					ch4_surf_ebul = ch4_ebul_total ! [mol/m2/s]
				endif
	
				! Ensure that concentrations stay above 0
				! This should be done after the flux, so that the flux calculation is consistent.
				do j = 1,nl_soil
					if (conc_ch4_rel(j) < 0._r8) then
						deficit = - conc_ch4_rel(j)*epsilon_t(j,1)*dz_soisno(j)  ! Mol/m^2 added
						if (deficit > 1.e-3_r8 * scale_factor_gasdiff) then
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
	
				enddo ! j; nl_soil
  
				call Tridiagonal(0, nl_soil, jtop, &
							at(:), &
							bt(:), &
							ct(:), &
							rt(:), &
							conc_o2_rel(0:nl_soil))
  
				! Ensure that concentrations stay above 0
				do j = 1,nl_soil
					conc_o2_rel(j) = max (conc_o2_rel(j), 1.e-12_r8)
					! In case of pathologically large aerenchyma conductance. Should be OK in general but
					! this will maintain stability even if a PATCH with very small weight somehow has an absurd NPP or LAI.
					! Also, oxygen above ambient will probably bubble.
					conc_o2_rel(j) = min (conc_o2_rel(j), c_atm(2)/epsilon_t(j,2))
				enddo
  
		   endif  ! species
  
		enddo  ! species
  
		! Update absolute concentrations per unit volume
		do j = 1,nl_soil ! No need to update the atm. level concentrations
		   conc_ch4(j) = conc_ch4_rel(j)*epsilon_t(j,1)
		   conc_o2(j)  = conc_o2_rel(j) *epsilon_t(j,2)
		end do
  
		! Do Balance Check and absorb small
		!    discrepancy into surface flux.
		do j = 1,nl_soil
		   if (j == 1) errch4 = 0._r8
		   errch4 = errch4 + (conc_ch4(j) - conc_ch4_bef(j))*dz_soisno(j)
		   errch4 = errch4 - ch4_prod_depth(j)*dz_soisno(j)*deltim
		   errch4 = errch4 + ch4_oxid_depth(j)*dz_soisno(j)*deltim
		end do
  
		
		! For history make sure that grnd_ch4_cond includes snow, for methane diffusivity
		grnd_ch4_cond = spec_grnd_cond(1)
  
		errch4 = errch4 + (ch4_surf_aere + ch4_surf_ebul + ch4_surf_diff)*deltim
  
		if (abs(errch4) < 1.e-8_r8) then 
		   ch4_surf_diff = ch4_surf_diff - errch4/deltim
		else ! errch4 > 1e-8 mol / m^2 / timestep
		   ! write(6,*)'CH4 Conservation Error in CH4Mod during diffusion, istep, errch4 (mol /m^2.timestep)', &
			  ! 	istep,errch4
			  write(6,*)'Lat,Lon=',lat,lon
			  CALL CoLM_stop ()
		end if
  
	end subroutine ch4_tran

!------------------------------------------------------------------------
	subroutine Tridiagonal (lbj, ubj, jtop, a, b, c, r, u)
!-----------------------DESCRIPTION--------------------------------------
! Tridiagonal matrix solution

!-----------------------Argument---------- ------------------------------
		implicit none         
		integer , intent(in)    :: lbj, ubj   ! lbinning and ubing level indices
		integer , intent(in)    :: jtop       ! top level for each column [col]
		real(r8), intent(in)    :: a(lbj:ubj) ! "a" left off diagonal of tridiagonal matrix [col, j]
		real(r8), intent(in)    :: b(lbj:ubj) ! "b" diagonal column for tridiagonal matrix [col, j]
		real(r8), intent(in)    :: c(lbj:ubj) ! "c" right off diagonal tridiagonal matrix [col, j]
		real(r8), intent(in)    :: r(lbj:ubj) ! "r" forcing term of tridiagonal matrix [col, j]
		real(r8), intent(inout) :: u(lbj:ubj) ! solution [col, j]
!-----------------------Local Variables------------------------------
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

END MODULE MOD_ch4_test
! --------- EOP ----------

