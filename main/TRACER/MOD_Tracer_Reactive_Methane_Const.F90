#include <define.h>
#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Reactive_Methane_Const
!=======================================================================
! methane constants
!=======================================================================
!
! PRIMARY REFERENCES (full bibliographic detail for citations used in
! parameter comments and biome lookup function headers below):
!
!   [Bridgham 2013] Bridgham, S. D., Cadillo-Quiroz, H., Keller, J. K., &
!     Zhuang, Q. (2013). Methane emissions from wetlands: biogeochemical,
!     microbial, and modeling perspectives from local to global scales.
!     Global Change Biology, 19(5), 1325-1346.  doi:10.1111/gcb.12131
!
!   [Hamilton 1996] Hamilton, S. K., Sippel, S. J., & Melack, J. M. (1996).
!     Inundation patterns in the Pantanal wetland of South America
!     determined from passive microwave remote sensing.
!     Archiv fur Hydrobiologie, 137(1), 1-23.
!
!   [Holzapfel-Pschorn 1985] Holzapfel-Pschorn, A., Conrad, R., & Seiler, W.
!     (1985). Production, oxidation and emission of methane in rice paddies.
!     FEMS Microbiology Ecology, 1(6), 343-351.
!     doi:10.1111/j.1574-6968.1985.tb01605.x
!
!   [Le Mer & Roger 2001] Le Mer, J., & Roger, P. (2001).  Production,
!     oxidation, emission and consumption of methane by soils: A review.
!     European Journal of Soil Biology, 37(1), 25-50.
!     doi:10.1016/S1164-5563(01)01067-6
!
!   [Marani & Alvalá 2007] Marani, L., & Alvalá, P. C. (2007).  Methane
!     emissions from lakes and floodplains in Pantanal, Brazil.
!     Atmospheric Environment, 41(8), 1627-1633.
!     doi:10.1016/j.atmosenv.2006.10.046
!
!   [Pangala 2017] Pangala, S. R., Enrich-Prast, A., Basso, L. S., Peixoto,
!     R. B., Bastviken, D., Hornibrook, E. R. C., Gatti, L. V., Marotta, H.,
!     Calazans, L. S. B., Sakuragui, C. M., Bastos, W. R., Malm, O., Gloor,
!     E., Miller, J. B., & Gauci, V. (2017).  Large emissions from
!     floodplain trees close the Amazon methane budget.
!     Nature, 552(7684), 230-234.  doi:10.1038/nature24639
!
!   [Walter 2001] Walter, B. P., Heimann, M., & Matthews, E. (2001).
!     Modeling modern methane emissions from natural wetlands: 1. Model
!     description and results.
!     J. Geophys. Res. Atmos., 106(D24), 34189-34206.
!     doi:10.1029/2001JD900165
!
!   [Wania 2010] Wania, R., Ross, I., & Prentice, I. C. (2010).
!     Implementation and evaluation of a new methane model within a
!     dynamic global vegetation model: LPJ-WHyMe v1.3.1.
!     Geoscientific Model Development, 3(2), 565-584.
!     doi:10.5194/gmd-3-565-2010
!
!   [Whalen & Reeburgh 1990] Whalen, S. C., & Reeburgh, W. S. (1990).
!     Consumption of atmospheric methane by tundra soils.
!     Nature, 346(6280), 160-162.  doi:10.1038/346160a0
!
! IMPORTANT: Specific numeric values in DEF_METHANE biome lookup arrays
! (f_methane_*, redoxlag_*, z0_methane_prod, hybrid_soil_threshold) are
! AUTHOR-SELECTED midpoints or order-of-magnitude estimates within the
! ranges discussed in these references, NOT direct numerical quotations.
! Default 'hybrid' mode is tuned against Pantanal in-situ flux
! [Marani & Alvalá 2007].
!=======================================================================
	USE MOD_Precision
	USE MOD_SPMD_Task, only: CoLM_Stop, p_is_master
	USE MOD_Tracer_Defs, only: tracer_lower
	USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite

	IMPLICIT NONE

	PRIVATE
	PUBLIC :: DEF_METHANE, DEF_METHANE_hydrology
	PUBLIC :: read_methane_namelist, configure_methane_inundation_mode
	PUBLIC :: methane_atm_mixing_ratio, methane_history_enabled
	PUBLIC :: methane_history_accumulation_mode

	!------------------------------------------------------------------
	! Constants for the Methane reactive tracer
	!------------------------------------------------------------------
	! Note some of these constants are also used in CNNitrifDenitrifMod

	integer :: iloop  ! loop index

	integer, public, parameter :: ngases      =   3     ! CH4, O2, & CO2
	integer, public, parameter :: METHANE_COMP_SOIL = 1
	integer, public, parameter :: METHANE_COMP_RICE = 2
	integer, public, parameter :: N_METHANE_COMP    = 2

	!------------------------------------------------------------------

	real(r8), public, parameter :: catomw = 12.011_r8 ! molar mass of C atoms (g/mol)
	! Carbon content of soil organic matter, gC per kg OM. Converts between
	! cellorg [kg OM m-3] and organic carbon: BgcLink divides by it going
	! C -> OM, the microbial substrate pools multiply by it going OM -> C.
	! Named because a bare 580 repeated across files is one edit away from
	! the two directions disagreeing, and nothing would catch that.
	real(r8), public, parameter :: gc_per_kg_om = 580._r8
	real(r8), public, parameter :: methane_atomw = 16.04_r8 ! molar mass of CH4 atoms (g/mol)

	real(r8), public :: s_con(ngases,4)    ! Schmidt # calculation constants (spp, #)
	data (s_con(1,iloop),iloop=1,4) /1898.0_r8, -110.1_r8, 2.834_r8, -0.02791_r8/ ! CH4
	data (s_con(2,iloop),iloop=1,4) /1801.0_r8, -120.1_r8, 3.7818_r8, -0.047608_r8/ ! O2
	data (s_con(3,iloop),iloop=1,4) /1911.0_r8, -113.7_r8, 2.967_r8, -0.02943_r8/ ! CO2

	real(r8), public :: d_con_w(ngases,3)    ! water diffusivity constants (spp, #)  (*10^-9 m2/s)
	data (d_con_w(1,iloop),iloop=1,3) /0.9798_r8, 0.02986_r8, 0.0004381_r8/ ! CH4
	data (d_con_w(2,iloop),iloop=1,3) /1.172_r8, 0.03443_r8, 0.0005048_r8/ ! O2
	data (d_con_w(3,iloop),iloop=1,3) /0.939_r8, 0.02671_r8, 0.0004095_r8/ ! CO2

	real(r8), public :: d_con_g(ngases,2)    ! gas diffusivity constants (spp, #) (*10^-4 m2/s)
	data (d_con_g(1,iloop),iloop=1,2) /0.1875_r8, 0.0013_r8/ ! CH4
	data (d_con_g(2,iloop),iloop=1,2) /0.1759_r8, 0.00117_r8/ ! O2
	data (d_con_g(3,iloop),iloop=1,2) /0.1325_r8, 0.0009_r8/ ! CO2

	real(r8), public :: c_h(ngases)    ! constant (K) for Henry's law (4.12, Wania)
	data c_h(1:3) /1600._r8, 1500._r8, 2400._r8/ ! CH4, O2, CO2

	real(r8), public :: kh_theta(ngases)    ! Henry's constant (mol/L/atm) at standard temperature (298K)
   data kh_theta(1:3) /1.4e-3_r8, 1.3e-3_r8, 3.4e-2_r8/ ! CH4, O2, CO2

	real(r8), public :: kh_tbase = 298.15_r8 ! base temperature for calculation of Henry's constant (K)

!------------------------------------------------------------------

	real(r8), public, parameter :: rgasm = 8.31446261815324_r8 ! Universal gas constant [J mol-1 K-1]
                                 ![J/K/mol]=[N*m/K/mol]=[Pa*m3/K/mol]
	!!! rgas Different from CoLM
	!!! rgas in CoLM is gas constant for dry air [J/kg/K]
	!!! not Universal gas constant
	real(r8), public, parameter :: rgasLatm = 0.08206_r8 ! L*atm/mol/K
   ! rgasLatm      = rgasm         / 101325 / 0.001
   ! [L*atm/mol/K] = [Pa*m3/K/mol] / [Pa/atm] / [m3/L]

	real(r8), public, parameter :: secspday = 86400._r8 ! Seconds per day

   type Methane_type
      ! ---------------------------------------------------- CH4 species controls ----------------------------------------------------
      ! These controls are private to the CH4 reactive tracer and are read from
      ! standard_ch4_parameter.nml with the rest of DEF_METHANE.  Keep the main
      ! nl_colm namelist species-neutral; it should only map reactive tracer
      ! names to parameter files via DEF_TRACER_PARAM_FILES.
      character(len=32) :: inundation_mode = 'hybrid'
      logical :: enable_rice_paddy = .false.
      logical :: only_wetland = .false.
      logical :: use_spatial_ph = .false.

      ! ---------------------------------------------------- Variable parameter ----------------------------------------------------
      ! methane production constants
      real(r8) :: q10methane =2._r8             ! additional Q10 for methane production ABOVE the soil decomposition temperature relationship (doc:Q10 Baseline:2 Range:1.5~4) (params:2.0)
      real(r8) :: f_methane = 0.2_r8            ! ratio of CH4 production to total C mineralization (Baseline:0.2 Range:NA params:0.2)

      ! Biome-specific f_methane lookup (CH4 yield = CH4 / total anaerobic decomp).
      ! When use_biome_f_methane=.true., per-patch f_methane is selected by
      ! climate zone (mirrors get_wetland_veg_proxy classification) instead
      ! of the global default 0.20.  Default .false. preserves backwards compat.
      !
      ! Values below are AUTHOR-SELECTED midpoints within published literature
      ! ranges, NOT direct numerical quotes from the cited papers.  The cited
      ! references provide the underlying f_methane / yield ranges (typically
      ! 0.03-0.50 across biomes per Bridgham 2013), from which the specific
      ! values below were chosen to match Pantanal in-situ flux (Marani &
      ! Alvalá 2007).  Re-tuning against other in-situ sites is encouraged.
      logical  :: use_biome_f_methane            = .false.
      real(r8) :: f_methane_tropical_peat        = 0.12_r8  ! Inferred from Bridgham 2013 GCB tropical range; Saunders et al. papyrus productivity work
      real(r8) :: f_methane_tropical_floodplain  = 0.07_r8  ! Tuned to match Marani & Alvalá 2007 Pantanal in-situ (200-300 / 30-80 mg/m²/d); within Bridgham 2013 0.03-0.15 tropical range
      real(r8) :: f_methane_floodplain           = 0.10_r8  ! Routing-activated floodplain / seasonally flooded non-rice soil CH4 yield
      real(r8) :: f_methane_temperate_marsh      = 0.15_r8  ! Inferred from Bridgham 2013 GCB Typha/Carex range 0.10-0.25
      real(r8) :: f_methane_boreal_fen           = 0.12_r8  ! Inferred from Wania et al. 2010 GMD LPJ-WHyMe fen range 0.10-0.20; lower than CTSM global 0.20
      real(r8) :: f_methane_boreal_bog           = 0.08_r8  ! Inferred from Bridgham 2013 Sphagnum range 0.05-0.15 (fermenter-dominated)
      real(r8) :: f_methane_rice_paddy           = 0.10_r8  ! Rice paddy CH4/CO2 ratio; conservative value within CLM4Me sensitivity range (0.1-0.3)
      real(r8) :: f_methane_upland_soil          = 0.05_r8  ! Inferred from Le Mer & Roger 2001 review (upland anaerobic microsites)

      ! Biome-specific redoxlag lookup (days).  When enabled, methane_prod
      ! consumes biome_redoxlag_patch(ipatch) instead of global DEF_METHANE%redoxlag.
      ! Concept: faster microbial response in warm tropical wetlands, slower
      ! in cold boreal peat.  Default .false. preserves backwards compatibility.
      !
      ! Values below are AUTHOR-CHOSEN order-of-magnitude estimates based on
      ! the qualitative conclusions in the cited papers, NOT direct numerical
      ! quotes.  Boreal slow-response from Whalen & Reeburgh 1990 (cold-microbe
      ! consumption); tropical rapid-response from Pangala et al. 2017 (Amazon
      ! tree stem CH4); CTSM legacy ~30d from Walter & Heimann 2001.
      ! Less constrained than f_methane lookup — wider uncertainty (50-100%).
      logical  :: use_biome_redoxlag             = .false.
      real(r8) :: redoxlag_tropical_peat         = 12._r8   ! Tropical mid-range
      real(r8) :: redoxlag_tropical_floodplain   = 7._r8    ! Tropical rapid (Pangala et al. 2017 Nature 552:230 qualitative)
      real(r8) :: redoxlag_temperate_marsh       = 20._r8   ! Temperate intermediate
      real(r8) :: redoxlag_boreal_fen            = 30._r8   ! CTSM legacy default (Walter & Heimann 2001 JGR ~14-30d transition)
      real(r8) :: redoxlag_boreal_bog            = 45._r8   ! Boreal slowest (cold Sphagnum, Whalen & Reeburgh 1990 Nature)
      real(r8) :: redoxlag_rice_paddy            = 25._r8   ! Managed flood/drain redox response; faster than wetland 30d but avoids unrealistically immediate production
      real(r8) :: redoxlag_upland_soil           = 30._r8   ! Default (only when Plan A activates upland CH4)
      ! redoxlag_wetland_dim: with wetland_dim_class and
      !   use_biome_redoxlag every wetland tile takes this one lag in place of
      !   the class values above, for the newly inundated area
      !   (finundated_lag) and the newly saturated layers (layer_sat_lag)
      !   alike; rice and upland soil keep theirs. Anoxic incubations of
      !   northern soils give a median lag of 27 d for samples from the zone
      !   of a fluctuating water table (mean 20 d, n 23; 10 d when permanently
      !   inundated), and the lags did not differ among incubation
      !   temperatures (Treat et al. 2015, GCB 21, 2787-2803, Table 2).
      real(r8) :: redoxlag_wetland_dim           = 27._r8   ! [days]

      ! Explicit methanogenesis depth attenuation.  CTSM BGC SOC profile
      ! decays with z0_BGC ~ 0.5 m, but boreal peatland observations
      ! (see Walter & Heimann 2001 JGR 106:34189) indicate CH4 production
      ! concentrates in the top 30 cm.  When z0_methane_prod > 0, partition_z
      ! is multiplied by exp(-z/z0_methane_prod) and the column re-normalized.
      ! NOTE: the normalisation fixes the sum of the BASE weights only.  Each
      ! layer is then multiplied by its own temperature, pH and redox factor,
      ! so the column total is NOT preserved -- see the explicit statement in
      ! methane_prod (Physics.F90, "does not claim to preserve final CH4
      ! production after layer-specific ... modifiers").  This comment claimed
      ! conservation until 2026-07-30; it was wrong, and the practical
      ! consequence is that z0_methane_prod moves column totals, not just the
      ! shape of the profile.  Default 0 = disabled, falls back to CTSM BGC
      ! profile alone.
      !
      ! Recommended 0.30 m: chosen by author within the 0.2-0.5 m range
      ! discussed in Walter & Heimann 2001 (specific value 0.30 not directly
      ! quoted from paper).  Tropical wetlands may require shorter z0 (~0.15 m,
      ! Pangala 2017) — biome-specific z0 not yet implemented.
      real(r8) :: z0_methane_prod = 0._r8

      ! methane oxidation constants
	      ! Oxidation vmax is per aqueous/active-water volume; methane_oxid
	      ! multiplies by vol_aqu to produce bulk-soil mol m-3 s-1 rates.
	      real(r8) :: vmax_methane_oxid = 1.25e-5_r8       ! Unit: mol m-3-aqueous s-1
	      real(r8) :: vmax_oxid_unsat = 1.25e-6_r8         ! Unit: mol m-3-aqueous s-1
      real(r8) :: k_m = 5.e-3_r8                ! Michaelis-Menten oxidation rate constant for CH4 concentration (params:5e-3 code:5.e-6_r8 * 1000._r8) (doc:KCH4 Baseline:5e-3 Range:5e-4~5e-2)
      real(r8) :: k_m_unsat = 5.e-4_r8           ! Michaelis-Menten oxidation rate constant for CH4 concentration (params:5e-4 code:5.e-6_r8 * 1000._r8 / 10._r8) (doc:KCH4 Baseline:5e-3 Range:5e-4~5e-2)
      real(r8) :: k_m_o2 =2.e-2_r8             ! Michaelis-Menten oxidation rate constant for O2 concentration (params:2e-2 code:20.e-6_r8 * 1000._r8) (doc:KO2 Baseline:2e-2 Range:2e-3~2e-1)
      real(r8) :: q10_methane_oxid = 1.9_r8         ! Q10 oxidation constant (? params:1.9)
      ! Lake-only oxidation controls.  Defaults preserve the legacy/global
      ! methane oxidation kinetics; use these instead of changing k_m_o2 or
      ! vmax_methane_oxid globally when diagnosing lake CH4.
      real(r8) :: lake_oxid_scale = 1.0_r8          ! multiplicative factor on lake CH4 oxidation only
      real(r8) :: lake_k_m_o2 = -1.0_r8             ! >0 overrides k_m_o2 for lake oxidation only
      real(r8) :: lake_vmax_methane_oxid = -1.0_r8  ! >=0 overrides vmax_methane_oxid for lake oxidation only
      real(r8) :: lake_oxic_sediment_depth = -1.0_r8 ! m; >0 limits lake sediment oxidation to this top depth, -1 disables

      ! Experimental microbial-pool state.  Runtime validation rejects
      ! enabling it until biomass growth/loss is coupled to donor/sink carbon.
      logical  :: use_microbial_pools = .false.
      logical  :: use_microbial_flux_override = .false.
      logical  :: use_microbial_dormancy = .false.
      real(r8) :: B_init_methanogen = 1.0_r8
      real(r8) :: B_init_methanotroph = 1.0_r8
      real(r8) :: B_min_methanogen = 1.0e-3_r8
      real(r8) :: B_min_methanotroph = 1.0e-3_r8
      ! Biomass caps as fractions of local layer organic carbon
      ! (cellorg [kg OM m-3] * 580 gC kgOM-1).  Set <=0 to disable only
      ! this stricter fractional limit; the shared physical C ceiling remains.
      real(r8) :: B_max_fraction_methanogen = -1.0_r8
      real(r8) :: B_max_fraction_methanotroph = -1.0_r8
      real(r8) :: mu_max_methanogen = 0.2_r8
      real(r8) :: mu_max_methanotroph = 0.5_r8
      real(r8) :: gamma_methanogen = 0.05_r8
      real(r8) :: gamma_methanotroph = 0.10_r8
      real(r8) :: gamma_microbial_dormant = 0.005_r8
      real(r8) :: gamma_microbial_freeze = 0.001_r8
      real(r8) :: K_substrate_methanogen_pool = 0.04_r8   ! substrate-pool half-saturation [mol C m-3]
      real(r8) :: K_inh_O2_methanogen = 1.0e-3_r8
      real(r8) :: kappa_m_methanogen = 1.0e-2_r8
      real(r8) :: kappa_m_methanotroph = 1.0e-2_r8
      ! When microbial flux override is enabled, do not allow microbial
      ! production to exceed this multiple of the legacy production tendency.
      real(r8) :: max_microbe_prod_multiplier = 3.0_r8
      real(r8) :: q10_microbe_growth = 3.0_r8
      real(r8) :: T_ref_microbe = 298.15_r8
      real(r8) :: dormancy_rate_active = 0.1_r8
      real(r8) :: dormancy_rate_revive = 1.0_r8
      real(r8) :: dormancy_threshold_methanogen_fS = 0.1_r8
      real(r8) :: dormancy_threshold_methanogen_fO2 = 0.1_r8
      real(r8) :: dormancy_threshold_methanotroph_fS = 0.1_r8
      real(r8) :: dormancy_threshold_methanotroph_fO2 = 0.1_r8

      ! methane ebullition constants
      real(r8) :: vgc_max  =0.15_r8            ! gas-volume fraction scale for ebullition threshold/target (unitless)

      ! methane aerenchyma constants
      real(r8) :: nongrassporosratio = 1._r8/3._r8   ! Ratio of root porosity in non-grass to grass, used for aerenchyma transport (params:0.33)
      real(r8) :: poros_tiller = 0.3_r8    ! Porosity for grass tiller
      real(r8) :: unsat_aere_ratio = 0.05_r8/0.3_r8    ! Ratio to multiply upland vegetation aerenchyma porosity by compared to inundated systems (params:0.1666666667 code:0.05_r8 / 0.3_r8)
      real(r8) :: porosmin = 0.05_r8            ! minimum aerenchyma porosity (unitless)(params:0.05 code:0.05_r8)
      real(r8) :: aere_radius = 2.9e-3_r8 ! Aerenchyma radius
      real(r8) :: rob = 3._r8                 ! ratio of root length to vertical depth ("root obliquity") (params:3. code:3._r8)
      real(r8) :: scale_factor_aere = 1._r8   ! scale factor on the aerenchyma area for sensitivity tests (1 params:1.) (doc:Fa Baseline:1 Range:0.5~1.5)

      ! methane transport constants
      real(r8) :: scale_factor_gasdiff = 1._r8   ! For sensitivity tests; convection would allow this to be > 1(? params:1.) (doc:fD0 Basline:1 Range:1,10 Unit:m2 s-1)
      real(r8) :: scale_factor_liqdiff = 1._r8   ! For sensitivity tests; convection would allow this to be > 1(? params:1.) (doc:fD0 Basline:1 Range:1,10 Unit:m2 s-1)
      real(r8) :: lake_liqdiff_scale = 1._r8      ! Lake-only multiplier on saturated/liquid CH4 diffusivity
      real(r8) :: lake_o2_liqdiff_scale = 1._r8   ! Lake-only multiplier on saturated/liquid O2 diffusivity
      ! Fallback atmospheric CH4/O2 boundary conductance [m/s], used only
      ! when current CoLM Monin-Obukhov state is absent or invalid.
      real(r8) :: grnd_methane_cond_default = 1.e-6_r8
      ! grnd_methane_cond_host: the boundary-layer
      !   conductance of a vegetated non-lake tile is 1/(rd + raw), the host
      !   water vapour resistance from the ground through the canopy air to
      !   the reference height (MOD_Vars_1DFluxes, raw_grnd), instead of
      !   grnd_methane_cond_default. CLM does the same: grnd_ch4_cond =
      !   1/(raw(p,above_canopy)+raw(p,below_canopy)) under a canopy
      !   (CanopyFluxesMod.F90) and 1/raw on bare ground
      !   (BareGroundFluxesMod.F90). Off by default.
      logical  :: grnd_methane_cond_host = .false.

      ! -------------------------------------------------- Invariant parameter ----------------------------------------------------
      ! methane production constants
      real(r8) :: mino2lim = 0.2_r8         ! minimum anaerobic decomposition rate as a fraction of potential aerobic rate (0.2+ params:0.2)
      real(r8) :: q10methane_base = 295._r8 ! temperature at which the effective f_methane actually equals the constant f_methane (295+ params:295)
      ! q10methane_local_base: the reference temperature of
      !   the production temperature factor of a wetland tile is the tile's own
      !   annual mean soil temperature over the top metre instead of
      !   q10methane_base: Q10 describes the seasonal response about the
      !   site's annual mean (Walter and Heimann 2000, eq. 5), so f_methane is
      !   the anaerobic CH4 share at the site's mean temperature. The previous
      !   full 365-day mean is used; q10methane_base until the first one is
      !   complete. Off by default.
      logical  :: q10methane_local_base = .false.
      ! q10methane_base_unfrozen: with q10methane_local_base,
      !   the reference is the mean top-metre soil temperature over the
      !   unfrozen part of the last calendar year (steps whose column mean is
      !   above freezing) instead of the full 365-day mean. The within-site
      !   seasonal temperature dependence of CH4 emission holds about the
      !   average temperature of the measured, thawed season, and emission at a
      !   fixed temperature falls with that site average across wetlands
      !   (Yvon-Durocher et al. 2014, Extended Data Fig. 2e); frozen months,
      !   when methanogens are inactive, do not set it. Off by default.
      logical  :: q10methane_base_unfrozen = .false.
      ! q10methane_base_weight: with q10methane_local_base,
      !   the reference is q10methane_base + weight x (tile mean -
      !   q10methane_base): 1 takes the tile's own mean (full acclimation),
      !   0 the fixed base; in between, CH4 emission at a fixed temperature
      !   still falls with the site mean but less than fully (the geographic
      !   response is weaker than the seasonal one, Yvon-Durocher et al.
      !   2014). Between 0 and 1.
      real(r8) :: q10methane_base_weight = 1._r8
      real(r8) :: q10lakebase = 298._r8       ! (K) base temperature for lake CH4 production (params:298. code:298._r8)
      ! q10lake: Q10 of lake CH4 production; <= 0 keeps the
      !   CLM4Me tie to 1.5 x q10methane, so a change of the wetland Q10 does
      !   not carry over to lakes when this is set.
      real(r8) :: q10lake = -1._r8
      ! methanogen_activity: CH4 production of a
      !   wetland tile (patchtype 2) is multiplied by the methanogen activity
      !   level m_act of each layer of each subcolumn, the dormancy variable of
      !   JULES-microbe (Chadburn et al. 2020, GBC 34, e2020GB006678, eq. A5,
      !   p.18 of 28). Methanogen biomass and dissolved substrate are held at
      !   the steady state they reach at full activity, where production is
      !   half the substrate supply (eqs. A2, A3, A7), i.e. the present
      !   production; production is then that value times m_act, and no
      !   microbial carbon pool is carried. m_act rises as f k2 a A2 m_act
      !   while the specific growth rate f k2 a A2 exceeds the dormancy
      !   threshold mu in a layer below the water table that has substrate
      !   supply (growth stops above the table, p.14), and falls at the same
      !   rate otherwise, between alpha/f and 1 (eq. A5, p.18); substrate is
      !   taken as not limiting there (phi = 1, p.19). A2 = q2**(0.1 T / (1 +
      !   T/273.15)), T in C (eq. A4, p.18); a = A(T_ann, qev_ratio q2)**-1 is
      !   the acclimation of eq. A6 (p.19) with T_ann the tile's last annual
      !   mean soil temperature over the top metre in place of the 5-yr
      !   relaxation (a = 1 until one exists; qev_ratio <= 0: no acclimation).
      !   The carbon a dormant community leaves unconverted goes to CO2, so
      !   closure against the BGC pools holds, and it adds no O2 demand.
      !   m_act starts at 1 (the present scheme); rice, floodplain and lake
      !   are unchanged. An alternative to q10methane_local_base and its
      !   options, not combined with it. Off by default.
      logical  :: methanogen_activity  = .false.
      real(r8) :: methanogen_k2        = 0.01_r8     ! [hr-1] k2, maximum respiration per unit biomass at 0 C (Table A1, p.17; Price and Sowers 2004, p.20)
      real(r8) :: methanogen_q2        = 4.3_r8      ! [-] Q2, temperature sensitivity of methanogens (Table A1, p.17; Yvon-Durocher et al. 2014, p.20)
      real(r8) :: methanogen_cue       = 0.03_r8     ! [-] f, carbon use efficiency, 0.6/19 (Table A1, p.17; p.20)
      real(r8) :: methanogen_alpha     = 0.001_r8    ! [-] alpha, maintenance respiration as a fraction of k2; activity floor alpha/f (Table A1, p.17; p.18)
      real(r8) :: methanogen_mu        = 0.00042_r8  ! [hr-1] mu, growth-rate threshold for dormancy and reactivation (Table A1, p.17; p.20)
      real(r8) :: methanogen_qev_ratio = 0.55_r8     ! [-] Q_ev / Q2 of the acclimation (Table A1, p.17; eq. A6, p.19; Bradford et al. 2019, p.20)
      ! acceptor_pool: a dynamic pool of alternative
      !   electron acceptors (NO3-, Mn, Fe(III), SO4--, humic quinones lumped
      !   as one, Segers and Leffelaar 2001, JGR 106, 3511-3528, p.3515) in
      !   each layer of each subcolumn of every non-lake patch (wetland,
      !   flooded soil, paddy). Its capacity is
      !     A_max = (1 - f_om) acceptor_cap_min + f_om acceptor_cap_org,
      !   f_om = min(cellorg / organic_max, 1) the layer's organic share as in
      !   the gas diffusivity of methane_tran: the humic capacity scales with
      !   organic matter (Wu and Blodau 2013, GMD 6, PEATBOG, p.1181; Keller
      !   and Bridgham 2007, abstract), the mineral part carries Fe(III).
      !   Below the water table the anaerobic carbon flow of the present
      !   scheme, 2 P (P the CH4 production, 2 CH2O -> CH4 + CO2), is split by
      !   eqs. 20, 22 and 24-27 (p.3515): methanogenesis takes the share
      !     zeta = 1 / (1 + eta0 A / (K + A)),
      !   A the oxidised acceptors, eta0 = Vm_rd K_mg / (K_rd Vm_mg) = 100 and
      !   K = K_rd,eo = 10 mol m-3 (Table 1, p.3512, with theta = 1 as there),
      !   so CH4 production becomes zeta P and the carbon of the other
      !   (1 - zeta) 2 P leaves as CO2, taking nu_rd = 4 mol e- per mol C,
      !   i.e. 8 mol e- per mol of CH4 not made (PEATBOG electron balance,
      !   eqs. A114, A120, A121, p.1181 and p.1200). Present acceptors all
      !   but stop methanogenesis (zeta << 1, p.3515); the pool empties at
      !   the rate the anaerobic carbon flow sets, "several days to several
      !   weeks" (p.3517), faster when warm. Where O2 is present the reduced
      !   part is reoxidised, eqs. 18, 28-30 (p.3515):
      !     dA/dt = f_aer f_liq k_ro (A_max - A),  f_aer = c_O2 / (c_O2 + K_ae),
      !   c_O2 the aqueous O2, f_liq the unfrozen water share (no change
      !   while frozen), at nu_ro = 0.25 mol O2 per mol e-, added to the O2
      !   demand and scaled by the layer's O2 stress like the other
      !   consumers. A water table that falls and rises again, a flood
      !   over dry soil (the pools move with the area between the
      !   subcolumns) and a layer that went into winter oxidised and thaws
      !   saturated all start with acceptors to use up first: the dynamic
      !   form of the redox lags (redoxlag, redoxlag_vertical, biome lags),
      !   which must be off with it. The CO2 made by the acceptor route is
      !   carbon the BGC already respired (net_methane only falls), so the
      !   carbon closure holds. Pools start full; the state is the oxidised
      !   share A / A_max. Lakes unchanged. Off by default; bit-for-bit
      !   unchanged when off.
      logical  :: acceptor_pool    = .false.
      real(r8) :: acceptor_cap_min = 5._r8       ! [mol e- m-3 soil] capacity of mineral soil; no global data, Segers typical total acceptors c_etot 5 (1-10) (Table 1, p.3512)
      real(r8) :: acceptor_cap_org = 5._r8       ! [mol e- m-3 soil] capacity of peat (f_om = 1); c_etot 5 (1-10) (Table 1, p.3512); PEATBOG 2000-4000 mmol m-2 over the top 0.6 m, 3.3-6.7 (p.1181)
      real(r8) :: acceptor_k_half  = 10._r8      ! [mol e- m-3] K_rd,eo, half saturation of acceptor reduction, 10 (1-100) mM (Table 1, p.3512; eq. 27)
      real(r8) :: acceptor_eta0    = 100._r8     ! [-] eta0 = Vm_rd K_mg / (K_rd Vm_mg) = 1e-4 x 0.1 / (0.01 x 1e-5) (eq. 26, p.3515; Table 1, p.3512)
      real(r8) :: acceptor_k_reox  = 1.e-5_r8    ! [s-1] k_ro, reoxidation rate constant, 1e-5 (1e-6 to 1e-4) (Table 1, p.3512; eq. 28)
      real(r8) :: acceptor_k_o2    = 0.02_r8     ! [mol O2 m-3 water] K_ae,O2, 20 (0.3-40) uM (Table 1, p.3512; eq. 18; "about 0.01 mol m-3", p.3516)
      real(r8) :: cnscalefactor=1._r8        ! scale factor on CN decomposition for assigning methane flux (?- params:1.)

      real(r8) :: redoxlag =30._r8           ! Number of days to lag in the calculation of finundated_lag (30+ params:30.)
      real(r8) :: lake_decomp_fact =9.e-11_r8    ! Base decomposition rate (1/s) at 25C (1 params:9e-11)
      real(r8) :: redoxlag_vertical=30._r8   ! time lag (days) to inhibit production for newly unsaturated layers (30+ params:0.)
      real(r8) :: pHmax = 9._r8          ! maximum pH for methane production(params:9. code:9.)
      real(r8) :: pHmin = 2.2_r8         ! minimum pH for methane production(params:2.2 code:2.2)
      ! ph_factor_floor: floor phi of the pH factor, phi + (1 - phi) D(pH).
      ! D is fitted to peat pushed off its own pH by buffers (Dunfield et al. 1993);
      ! at their native pH acidic peats and bogs keep 0.23-0.8 of their maximum
      ! (Dunfield et al. 1993 Fig. 5; Ye et al. 2012; Kotsyurbenko et al. 2007),
      ! so D overstates the difference between acidic and neutral wetlands.
      ! 0 (default) keeps D; outside [pHmin, pHmax] production takes phi.
      real(r8) :: ph_factor_floor = 0._r8   ! [-]
      real(r8) :: oxinhib = 400._r8          ! inhibition of methane production by oxygen (m^3/mol) (400+? params:400.)

      ! methane oxidation constants
      real(r8) :: smp_crit =-2.4e5_r8            ! Critical soil moisture potential (mm) (params:-2.4e5)

      ! methane ebullition constants
      real(r8) :: bubble_f =0.57_r8            ! CH4 content in gas bubbles (Kellner et al. 2006)

      ! methane aerenchyma constants
      real(r8) :: aereoxid =0._r8            ! fraction of methane flux entering aerenchyma rhizosphere that will be(? params:0.)
      real(r8) :: tiller_C = 0.22_r8 ! Per tiller 0.22 g C [g C/tiller]

      ! methane transport constants
      real(r8) :: satpow  =2._r8             ! exponent on watsat for saturated soil solute diffusion (2? params:2.)
      real(r8) :: capthick = 100._r8         ! min thickness before assuming h2osfc is impermeable (mm) (params:100.code:100._r8)

      ! additional constants
      real(r8) :: atm_methane  = 1.7e-6_r8         ! Atmospheric CH4 mixing ratio fallback (mol/mol)
      logical :: use_transient_atm_methane = .false. ! true: read time-varying atmospheric CH4 from atm_methane_file
      character(len=256) :: atm_methane_file = 'null' ! ASCII table: year value, or year month value
      character(len=16) :: atm_methane_file_units = 'auto' ! auto, mol/mol, ppmv, or ppbv
      real(r8) :: om_frac_sf = 1._r8          ! Scale factor for organic matter fraction (unitless)(? params:NA)

      ! -------------------------------------------------- Switch ----------------------------------------------------
      logical :: use_aereoxid_prog = .true. ! if false then aereoxid is read off of
      ! the parameter file and may be modifed by the user (default aereoxid on the
      ! file is 0.0).

      logical :: transpirationloss = .true. ! switch for activating CH4 loss from transpiration
                                    ! Transpiration loss assumes that the methane concentration in dissolved soil
                                    ! water remains constant through the plant and is released when the water evaporates
                                    ! from the stomata.
                                    ! Currently hard-wired to true; impact is < 1 Tg CH4/yr

      logical :: allowlakeprod = .false. ! Switch to allow production under lakes based on soil carbon dataset
                                 ! (Methane can be produced, and CO2 produced from methane oxidation,
                                 ! which will slowly reduce the available carbon stock, if ! replenishlakec, but no other biogeochem is done.)
                                 ! Note: switching this off turns off ALL lake methane biogeochem. However, 0 values
                                 ! will still be averaged into the concentration _sat history fields.

	      logical :: usephfact = .true.  ! Switch to use pH factor in methane production; spatial pH input is required when use_spatial_ph is true.

      ! Smooth WTD->finundated transition for scheme 6 (logistic S-curve).
      ! Replaces the original step function (zwt<=0.30 -> 1, else 0).
      ! Reference: Walter & Heimann 2000 (sigmoid f(zwt)), Sundh 2000
      ! (exp(-zwt/0.3)), Bridgham 2013 GCB (review).
      !   finundated = 1 / (1 + exp((zwt - wtd_inflection)/wtd_steepness))
      ! Default wtd_inflection=0.30m preserves the original threshold as the
      ! S-curve's symmetric center (finundated = 0.5 at zwt = 30 cm).
      real(r8) :: wtd_inflection = 0.30_r8   ! [m] zwt at which finundated = 0.5 (wetland tile, patchtype==2)
      real(r8) :: wtd_steepness  = 0.05_r8   ! [m] S-curve width; smaller = sharper

      ! Soil-tile sigmoid params (patchtype==0/1): deeper inflection so seasonally
      ! wet upland (Pantanal floodplain etc.) can also produce CH4 when soil zwt
      ! rises into the 0.5-1.5 m range during wet season.  Disabled if =0
      ! (default 0 keeps backwards compat: soil patches use wetland's sigmoid).
      real(r8) :: wtd_inflection_soil = 0._r8 ! [m] >0 enables soil-tile dynamic flood extension
      real(r8) :: wtd_steepness_soil  = 0.3_r8

      ! Hybrid mode: on soil tiles (patchtype != 2) take finundated from the
      ! routing-published f_inund_flood_patch instead of sigmoid(zwt).  This
      ! captures tropical floodplain overland flooding (Pantanal, Amazon) that
      ! sigmoid(zwt) cannot resolve because soil zwt is too deep there.
      ! Wetland tiles still use sigmoid(zwt) to keep dyn_wtd seasonality.
      ! Requires GridRiverLakeFlow active + DEF_USE_Dynamic_Wetland=.true.
      logical  :: use_routing_for_soil = .false.
      ! colm_floodplain_class: in the colm mode
      !   (scheme 8), which resets use_routing_for_soil although every soil
      !   tile takes the routing flood fraction, a soil tile flooded beyond the
      !   static wetland share is treated as a floodplain (floodplain redox lag
      !   and biome class, methane_area_floodplain diagnostic), as in the
      !   hybrid mode. Off by default.
      logical  :: colm_floodplain_class = .false.
      ! flooded_soil_ph_neutral: the flooded
      !   (saturated) subcolumn of a soil tile, floodplain or paddy, takes no pH
      !   factor. Flooded mineral soils converge to pH 6.7-7.2 as reduction
      !   consumes protons (Ponnamperuma, after Sahrawat 2015), near the
      !   methanogenesis optimum; rice models (CH4MOD, DNDC-Rice) carry no pH
      !   term and CLM4Me keeps the factor at 1. The dryland pH of the soil map
      !   otherwise cuts paddy and floodplain production by about a quarter.
      !   Peat wetland tiles keep their pH factor. Off by default.
      logical  :: flooded_soil_ph_neutral = .false.
      ! prod_fliq_off: methanogenesis in a freezing
      !   layer of a non-lake tile is not multiplied by the unfrozen water
      !   fraction on top of the decomposition moisture scalar, which already
      !   drops as the pore water freezes (double count). CTSM zeroes the
      !   production of a frozen layer only without vertically resolved soil
      !   carbon or in lakes (ch4Mod.F90, base_decomp). Cold-season emission is
      !   45-50% of the annual total in observations but 27% in the GCP
      !   models (Ito et al. 2023; Zona et al. 2016). Off by default.
      logical  :: prod_fliq_off = .false.

      ! Hybrid soil-tile threshold gate: when use_routing_for_soil=.true.,
      ! soil tile only produces CH4 when its routing fldfrc exceeds this
      ! threshold.  Below threshold, soil tile is treated as fully dry
      ! (finundated=0).  Conceptually represents hydrological connectivity:
      ! river overbank does not flood soil until stage exceeds bankfull.
      ! Default 0 disables the gate (back compat).
      !
      ! Recommended 0.05-0.10 for tropical 2-deg grid (empirically tuned
      ! against Pantanal P/T ratio match; not derived from any single paper).
      ! Pantanal routing fldfrc ranges 0.05-0.20 wet→dry, so threshold 0.05
      ! gates out only the driest months.
      real(r8) :: hybrid_soil_threshold    = 0._r8

      logical :: replenishlakec = .false. ! Finite lake-sediment carbon stock by default.
                                    ! Enabling replenishment imposes an external carbon source
                                    ! and is therefore suitable only for explicit sensitivity tests.

      logical :: wetland_fixed_substrate = .false. ! Hold the patchtype 2 decomposition pools at
                                    ! their initial values: heterotrophic respiration is still
                                    ! computed from them and consumed as CH4/CO2, but the pools are
                                    ! never debited.  The permanent wetland tile carries no PFT, so
                                    ! it receives no litter and can only drain; this switch makes the
                                    ! substrate a prescribed boundary condition instead.  Like
                                    ! replenishlakec it imposes an external carbon source and breaks
                                    ! closure against the BGC pools -- sensitivity tests only.

      logical :: methane_offline = .true.    ! Only offline land CH4 is implemented in this repository.
                                 ! Setting false is rejected during validation until a host atmosphere
                                 ! flux publisher and NEM-to-NEE coupling are available.

      logical :: methane_rmcnlim = .false.   ! Remove the N and low moisture limitations on SOM HR when calculating
                                 ! methanogenesis.
                                 ! Note: this option has not been extensively tested.
                                 ! Currently hardwired off.

      logical :: anoxicmicrosites = .false. ! Use Arah & Stephen 1998 expression to allow production above the water table
                                    ! Currently hardwired off; expression is crude.

      logical :: methane_frzout = .false.    ! Retired compatibility flag; must remain false.
                                 ! Ice is always excluded from mobile CH4 storage.

      ! public :: methane_conrd ! Read and initialize CH4 constants

      logical :: use_nitrif_denitrif = .true.

      logical :: anoxia  = .true. ! true => anoxia is applied to heterotrophic respiration also considered in CH4 model
                                    ! default value reset in controlMod
                                    ! Whether to enable the anoxia for the seasonally induatded zones
      logical :: use_vertical_redoxlag = .true. ! Whether to enable the vertical redox lag effect

      ! SIF is enabled only when coupled BGC does not already apply anoxia limits.
      logical :: bgc_anoxia_limits_decomp = .false.  ! true: BGC o_scalar already limits decomposition
      logical :: use_ch4_sif              = .true.   ! true: apply CH4 seasonal inundation factor
      ! floodplain_anoxic_decomp: on non-wetland, non-paddy
      !   patches CH4 is produced (below the water table, in both the flooded
      !   and the non-flooded subcolumn) from decomposition at the anoxic
      !   rate, mino2lim times the aerobic rate their BGC carries (o_scalar is
      !   1 there, and flooding never builds up their carbon stock). This is
      !   the CLM4Me seasonal inundation factor with no permanently inundated
      !   share (Riley et al. 2011, appendix B). Off by default.
      logical :: floodplain_anoxic_decomp = .false.
      ! anoxic_decomp_flooded_only: with floodplain_anoxic_decomp on, only the
      !   flooded subcolumn decomposes at the anoxic rate (seasonally scaled when
      !   floodplain_sif is on); the non-flooded subcolumn keeps the aerobic
      !   rate in every layer. Its layers below the water table are saturated
      !   for most of the year, and appendix B of Riley et al. (2011) has the
      !   carbon stock of a layer anoxic a share phi of the year build up to
      !   I tau / (phi beta + 1 - phi), so production returns to the input
      !   rate as phi goes to 1 instead of staying at beta (mino2lim) times
      !   it; the layers above the table are aerobic soil. Off by default
      !   (floodplain_anoxic_decomp in both subcolumns).
      logical :: anoxic_decomp_flooded_only = .false.
      ! wt_layer_sat_share: in the non-flooded subcolumn the
      !   layer holding the water table makes CH4 on the share of its
      !   thickness below the table only; the share above it is unsaturated
      !   soil (CH4 only from anoxic microsites when anoxicmicrosites is on).
      !   Off by default: the whole layer counts as below the table.
      logical :: wt_layer_sat_share = .false.
      ! wt_layer_node_transport: in the non-flooded
      !   subcolumn the layer holding the water table counts as saturated for
      !   gas and liquid transport, oxidation kinetics, ebullition and plant
      !   uptake only when the table is above its node. The off setting counts
      !   it saturated whatever the depth of the table inside it, which cuts
      !   the O2 of its air-filled part and ends the oxic zone at its top, on
      !   average half a layer above the table (0.17 m instead of 0.25 m for a
      !   table at 0.25 m); the node rule puts the end of the oxic zone at the
      !   layer interface nearest the table. Production keeps the interface
      !   rule (wt_layer_sat_share sets its share). Off by default.
      logical :: wt_layer_node_transport = .false.
      ! If BGC later applies real o_scalar limits, set (true, false).
      ! frozen_anoxic_decomp: a frozen wetland layer keeps
      !   the anoxic limit mino2lim on decomposition instead of 1. Ice brings
      !   no oxygen, and CLM applies o_scalar in every layer whatever its
      !   temperature. With 1, a layer held at the freezing point, where
      !   w_scalar is still near 1, decomposed five times faster than a
      !   thawed one. Off by default.
      logical :: frozen_anoxic_decomp = .false.
      ! liquid_fraction_scaling: CH4 production, oxidation and
      !   aerenchyma transport in a soil layer scale with its unfrozen water
      !   fraction wliq/(wliq+wice) instead of stopping at tfrz. Off by default.
      logical :: liquid_fraction_scaling = .false.
      ! wetland_oxic_cap_kinetics: the layers above a wetland's
      !   water table oxidise CH4 with the saturated-zone kinetics (k_m,
      !   vmax_methane_oxid) instead of the upland values. Off by default.
      logical :: wetland_oxic_cap_kinetics = .false.

      ! Global-run diagnostics / guard rails.
      ! write_ch4_history=false suppresses all CH4 history variables.
      ! ch4_history_vars accepts:
      !   'core'       : compact global-run set (default): surface flux,
      !                  production, oxidation, and total column CH4;
      !   'diagnostic' : previous broad core set for debugging/global audits;
      !   'all'        : legacy behavior, emit every CH4 diagnostic;
      !   'none'       : suppress all CH4 diagnostics;
      !   comma-separated NetCDF variable names for an exact custom set.
      logical :: write_ch4_history = .true.
      character(len=4096) :: ch4_history_vars = 'core'
      ! A positive value aborts the whole MPI job when any single-timestep
      ! CH4/O2 numerical correction exceeds this mol/m2 threshold.  A
      ! non-positive value keeps the existing diagnostic-only behavior.
      real(r8) :: numerical_correction_fatal_threshold = -1._r8
      ! Tolerance, as a FRACTION OF PORE VOLUME, on the disagreement between the
      ! host soil water (wliq/denh2o + wice/denice) and what the inundated
      ! fraction implies, when methane partitions the column into its saturated
      ! and unsaturated halves.
      !
      ! The two quantities are computed independently -- vtot from the host soil
      ! state, finundated from a water-table S-curve -- so they cannot be
      ! compared at machine epsilon. The guard this replaces used
      ! 1.e-10*max(pore_volume,1), i.e. a flat 1e-10 m absolute tolerance
      ! (pore_volume is ~0.04 m, so the max() always selected 1.0), which is
      ! round-off level and aborted global runs on their first timestep.
      !
      ! Within the tolerance the host water is left exactly as it is; the
      ! saturated-column allocation is capped instead, so the area-weighted
      ! total still reproduces the host state and no water is created. The
      ! largest residual seen is reported once. Beyond the tolerance the run
      ! still stops: a disagreement that large is invalid forcing, not
      ! numerical drift.
      real(r8) :: host_water_tolerance = 0.05_r8
      ! By default a missing/non-positive lake depth emits one warning and
      ! continues.  Set true for production global runs to fail fast when
      ! lake CH4 is enabled but landdata lacks valid lakedepth.
      logical :: lake_zero_depth_fatal = .false.
      ! Targeted restart-continuity diagnostics for lake CH4/O2 state.
      ! Default is off.  When enabled, lake patches print one compact line
      ! for selected days/times so restart-vs-continuous divergence can be
      ! diagnosed without changing physics.
      logical :: lake_restart_debug = .false.
      integer :: lake_restart_debug_year = 1996
      integer :: lake_restart_debug_start_doy = 1
      integer :: lake_restart_debug_end_doy = 3
      integer :: lake_restart_debug_sec = 1800

      ! Water-balance finundated override for wetland tiles.
      ! Replaces scheme-computed finundated on patchtype==2 with
      ! clamp(wetwat / wetwatmax, 0, 1).  Soil patches keep scheme finundated.
      logical :: enable_wetwat_finundated_override = .false.

      ! Force the wetland unsaturated branch to use a dry soil-column state.
      ! This keeps saturated and unsaturated branches distinct when wetland
      ! finundated is supplied by a seasonal area signal.
      logical :: wetland_dry_unsat_branch = .false.

      ! Paddy water management (removed by 68a507f8, restored from 7550e0ed).
      !
      ! These are methane-only: they set the inundation the CH4 column sees and
      ! do NOT irrigate the host. That is a known inconsistency -- the water is
      ! invented here and the host water balance never sees it -- accepted
      ! deliberately, because without it a paddy is hydrologically a dry field
      ! and the seven FLUXNET-CH4 rice towers all model exactly 0.0 against
      ! observed 66.7 mg CH4 m-2 d-1. Record it in the calibration archive: any
      ! parameter tuned on rice under this scheme carries the inconsistency,
      ! and a later coupled path through MOD_Irrigation must re-tune.
      !
      ! rice_paddy_min_finundated: floor on finundated while CN reports the crop
      !   alive. Blended with the scheme value by max(), so an already-wet patch
      !   -- a wetland tile carrying a rice CFT, or scheme 6 with the water table
      !   at the surface -- is never dried by it. 0 (default) leaves the paddy
      !   on the scheme value; run/paper_ch4_parameter.nml uses 0.85.
      real(r8) :: rice_paddy_min_finundated     = 0._r8
      ! Midseason drying: the standard Asian practice of draining for 7-10 days
      ! around 30-40 days after planting. Timing and depth are tunable; the
      ! defaults are the mid-range of that practice, not a fitted value.
      real(r8) :: rice_midseason_start_days     = 35._r8
      real(r8) :: rice_midseason_drain_days     = 10._r8
      real(r8) :: rice_midseason_drained_finundated = 0.30_r8

      ! Rice physiology/aerenchyma can remain active briefly after harvest.
      ! Also the window over which the paddy drains back to the host value.
      real(r8) :: rice_drain_window_days        = 30._r8

      ! wetland_max_wtd: deepest water table [m below surface] a wetland tile
      !   may reach under DEF_USE_Dynamic_Wetland. When the table falls
      !   deeper, WATER_VSF feeds the aquifer from the side (negative
      !   subsurface runoff, booked in rnof) until the table is back at this
      !   depth, standing in for the lateral inflow the closed wetland bucket
      !   never receives. A negative value switches it off. With
      !   use_biome_wetland_max_wtd the value is taken per wetland class
      !   (BIOME_* in BgcLink) instead of the single global value.
      real(r8) :: wetland_max_wtd                     = -1._r8
      logical  :: use_biome_wetland_max_wtd           = .false.
      real(r8) :: wetland_max_wtd_tropical_peat       = -1._r8
      real(r8) :: wetland_max_wtd_tropical_floodplain = -1._r8
      real(r8) :: wetland_max_wtd_temperate_marsh     = -1._r8
      real(r8) :: wetland_max_wtd_boreal_fen          = -1._r8
      real(r8) :: wetland_max_wtd_boreal_bog          = -1._r8
      ! wetland_max_wtd_bog_share: water-table floor
      !   z_bog [m below surface] of the open-bog share of a dynamic wetland
      !   tile. A bog is fed by precipitation alone and draws down in summer
      !   until its acrotelm is drained; the catotelm below stays saturated,
      !   and evaporation and outflow die away there. Models of rain-fed
      !   peatlands put that depth at 0.3 m: a measured maximum water-table
      !   depth z_b of 0.30 m, with actual evapotranspiration falling to zero
      !   at 0.28 m (Granberg et al. 1999, Table 1 p. 3776, eq. 8 p. 3774);
      !   a 0.3 m acrotelm over a saturated catotelm (LPJ-WHyMe, Wania et al.
      !   2009a, sect. 2.3.2 p. 6); saturation below 0.30 m (TEM, Zhuang et
      !   al. 2004, appendix D6 pp. 18-19). Boreal bogs reach a median deepest
      !   growing-season water table of 20 cm (quartiles 16-27.5 cm, 19
      !   sites) against 8 cm in fens (BAWLD-CH4 field WTMin, Kuhn et al.
      !   2021). The bog share b is the bog_share of wetland_veg_file (BAWLD
      !   BOG / (PEB + WTU + MAR + BOG + FEN), Olefeldt et al. 2021, from
      !   v2/scripts/mk_bog_share.py) or wetland_bog_share_site at a tower.
      !   Within the share that takes lateral inflow, 1 - s with s the
      !   rain-fed share, it weighs w_b = min(max(b / (1 - s), 0), 1),
      !   and that share's floor becomes the area mean
      !   (1 - w_b) z_class + w_b z_bog of the class floor z_class above
      !   (left off where z_class is off) and z_bog. With
      !   DEF_WETLAND_LATERAL_INFLOW the bog share takes none of the r R_up
      !   either, the inflow falls to (1 - s - b) r R_up, and the host refills
      !   the bog share up to z_bog instead. 0.30 is the Granberg depth,
      !   0.25-0.35 the sensitivity range; negative is off (default).
      real(r8) :: wetland_max_wtd_bog_share           = -1._r8  ! [m]

      ! wetland_peat_drainage: lateral outflow of a dynamic wetland through
      !   its saturated zone, Q = peat_c * T(zwt) (PEAT-CLSM, Bechtold et al.
      !   2019, JAMES 11, eqs 6-9). T integrates below the water table a
      !   macro-scale conductivity that falls with depth as
      !   peat_K0 / (1 + 100 z)^peat_m in peat and equals the layer's
      !   saturated conductivity in mineral soil; each layer mixes the two by
      !   its organic fraction OM_density / organic_max (Lawrence and Slater
      !   2008; organic_max 130 kg m-3 as in the CLM parameter file). peat_c
      !   is hydraulic gradient over flow length. The water leaves as positive
      !   subsurface runoff. Off by default.
      logical  :: wetland_peat_drainage = .false.
      real(r8) :: peat_K0               = 10._r8      ! [m s-1]
      real(r8) :: peat_m                = 3._r8       ! [-]
      real(r8) :: peat_c                = 1.5e-5_r8   ! [m-1]
      real(r8) :: organic_max           = 130._r8     ! [kg OM m-3]
      ! wetland_microtopo_sigma: standard deviation [m] of a normally
      !   distributed wetland surface (PEAT-CLSM microtopography). The
      !   inundated fraction of a dynamic wetland is then the share of the
      !   surface below the water level (ponded depth, else minus zwt); a
      !   value <= 0 keeps the step (inundated when zwt <= 0).
      real(r8) :: wetland_microtopo_sigma = -1._r8
      ! wetland_plant_input: plant carbon input of the permanent-wetland tile.
      !   The tile's own canopy assimilation times
      !   wetland_npp_frac (NPP/GPP) enters the litter pools every step along
      !   the tile's root profile, split labile/cellulose/lignin as CoLM's
      !   grass litter (0.25/0.5/0.25), with nitrogen at wetland_litter_cn.
      !   At steady state heterotrophic respiration then follows productivity
      !   instead of the initial stock. Needs wetland_fixed_substrate off.
      logical  :: wetland_plant_input = .false.
      real(r8) :: wetland_npp_frac    = 0.5_r8     ! [-]
      real(r8) :: wetland_litter_cn   = 46._r8     ! [g C / g N]
      ! wetland_bg_frac: share of that input that follows the
      !   root profile; the rest is aboveground litter laid on the surface
      !   along CoLM's leaf-litter profile. 1 keeps everything on the roots.
      real(r8) :: wetland_bg_frac     = 1._r8      ! [-]
      ! wetland_veg_glwd: replace the five-zone wetland
      !   vegetation proxy by the tile's GLWD make-up read from wetland_veg_file
      !   (forested share and class areas per grid cell). The forested share
      !   takes wetland_bg_frac_forest as its belowground input share (tropical
      !   peat swamp forest 0.07-0.23) and nongrassporosratio of the grass
      !   aerenchyma porosity; the non-forested share keeps its remote-sensing
      !   LAI only up to the measured peak, wetland_lai_open_peat for open
      !   peatland (GLWD 23, 25) and wetland_lai_marsh for marsh (17, 19, 27).
      logical  :: wetland_veg_glwd       = .false.
      character(len=256) :: wetland_veg_file = 'null'
      ! site_flood_file: single-point runs have no river
      !   routing, so a river-fed tower never floods. The file gives monthly
      !   flood_frac and flood_depth [m] per site (dims site, time; with lat,
      !   lon, year, month); the soil patch of the nearest site within 0.05
      !   degree takes them as the routing floodplain fraction and depth,
      !   linearly between mid-months (monthly climatology outside the years).
      character(len=256) :: site_flood_file = 'null'
      ! floodplain_glwd_cap_file: colm mode only. The routing
      !   flood of the soil tiles is bounded by the riverine and lacustrine
      !   floodplain of the cell in GLWD v2 (maximum extent over 1984-2020;
      !   Lehner et al. 2025, ESSD 17, 2277-2329, table 2): classes 8-15 and
      !   the large river deltas (30), which supersede the riverine classes
      !   inside the delta outlines (ibid., sect. 3.4). The flooded fraction is
      !   at most A(8-15, 30) / A(soil), A(soil) being the area of the cell's
      !   patchtype-0 patches. The same bound goes to the flood the land sees
      !   (the site flood series of DEF_FLOODPLAIN_INFILTRATION, and in gridded
      !   runs the flood fraction of DEF_GridRiverLake_FloodFeedback), so the
      !   methane and the soil water see one flooded area; the routing's
      !   storage and discharge are left alone. The file holds lat,
      !   lon, area_class_08 ... area_class_15 and area_class_30 [km2] on the
      !   model grid (v2/data/glwd33_2deg.nc). 'null': no bound.
      character(len=256) :: floodplain_glwd_cap_file = 'null'
      real(r8) :: wetland_bg_frac_forest = 0.15_r8   ! [-]
      real(r8) :: wetland_lai_open_peat  = 0.6_r8    ! [m2 m-2]
      real(r8) :: wetland_lai_marsh      = 3.0_r8    ! [m2 m-2]
      ! Single-point overrides of the gridded make-up (a tower sits in one
      ! wetland, the grid cell holds a mix); < 0 keeps the file value.
      real(r8) :: wetland_forest_share_site = -1._r8  ! [-]
      ! wetland_ombro_share_site: rain-fed (ombrotrophic) share
      !   of a tower's wetland, which receives no lateral inflow; < 0 keeps the
      !   rainfed_share of wetland_veg_file (0 when the file has none).
      real(r8) :: wetland_ombro_share_site  = -1._r8  ! [-]
      ! wetland_bog_share_site: open-bog share of a
      !   tower's wetland for wetland_max_wtd_bog_share; < 0 keeps the
      !   bog_share of wetland_veg_file (0 when the file has none).
      real(r8) :: wetland_bog_share_site    = -1._r8  ! [-]
      ! wetland_dim_class: the wetland class tree
      !   (latitude and modelled top-layer carbon, get_biome_f_methane) no
      !   longer decides the water-table floor, the pond litter resistance
      !   (DEF_WETLAND_POND_LITTER_RSS) and the moss surface resistance
      !   (wetland_moss_rss); each follows the share of the tile's wetland
      !   area under the process dimension that controls it, taken from the
      !   GLWD v2 classes of wetland_veg_file as shares of the tile classes
      !   16-19 and 22-27 (Lehner et al. 2025): emergent marsh s_e (17,
      !   isolated regularly flooded non-forested), tropical peat dome s_d
      !   (26-27) and moss-carpeted peatland s_m (boreal and temperate peat,
      !   22-25). With use_biome_wetland_max_wtd the floor is
      !     s_e wetland_max_wtd_temperate_marsh + s_d wetland_max_wtd_tropical_peat
      !       + (1 - s_e - s_d) wetland_max_wtd,
      !   the pond litter resistance is scaled by s_e in place of the
      !   temperate marsh class, and the moss resistance by s_m in place of
      !   the temperate marsh and boreal classes times the peatland share.
      !   The values are the class values; only which tile takes them
      !   changes (v2/results/design_wetland_class_260930.txt, sect. 5).
      !   Step 2: with use_biome_redoxlag every wetland tile takes
      !   redoxlag_wetland_dim, and with wetland_forest_input_herb the
      !   moss-carpeted forested share keeps wetland_moss_input_frac of the
      !   stand input. Step 3, the water-source dimension: the
      !   rain-fed share s_r, fed by precipitation alone (Bridgham et al.
      !   2013), is the permafrost peat plateau share s_p (BAWLD PEB,
      !   rainfed_share of wetland_veg_file) plus the open-bog share s_b
      !   (BAWLD BOG, bog_share), both of the BAWLD wetland area (Olefeldt et
      !   al. 2021). The plateaus, raised by ground ice, take no lateral
      !   inflow (as with DEF_WETLAND_INFLOW_AREA_SPLIT), so the fed
      !   share is 1 - (s_r - s_b); open bogs keep the floor refill, at
      !   wetland_max_wtd_bog_share when set. Only the tower keys
      !   below are new; with them unset the rain-fed and open-bog shares come
      !   from the wetland vegetation file. Needs wetland_veg_glwd. Off by default.
      logical  :: wetland_dim_class         = .false.
      ! A tower's own shares for wetland_dim_class, from its wetland type
      ! (a tower sits in one wetland, not in the cell's mix); < 0 keeps the
      ! share of wetland_veg_file. For the water source, the rain-fed share
      ! s_r and the open-bog share s_b <= s_r within it replace
      ! wetland_ombro_share_site (then s_r - s_b) and wetland_bog_share_site.
      real(r8) :: wetland_share_emerg_site  = -1._r8  ! [-]
      real(r8) :: wetland_share_dome_site   = -1._r8  ! [-]
      real(r8) :: wetland_share_moss_site   = -1._r8  ! [-]
      real(r8) :: wetland_share_rain_site   = -1._r8  ! [-]
      real(r8) :: wetland_share_bog_site    = -1._r8  ! [-]
      real(r8) :: wetland_lai_cap_site      = -1._r8  ! [m2 m-2]
      ! wetland_forest_htop_trop: canopy top height [m] of the
      !   forested share f of a tropical wetland tile (wetland_veg_glwd; within
      !   23.5 degrees of the equator, the tropical band of get_biome_f_methane).
      !   The host gives the tile the canopy height of its land class (IGBP 11,
      !   htop0 0.5 m, HTOP_readin), so a peat swamp forest of 23-35 m (Sakabe
      !   et al. 2018; Wong et al. 2018) gets the roughness length and
      !   displacement height of a meadow (LeafTemperature), a canopy weakly
      !   coupled to the air above it. With the key the tile's height is
      !     htop = f wetland_forest_htop_trop + (1 - f) htop0,
      !   averaged by share as the PFT heights of a soil patch (HTOP_readin);
      !   nothing else of the canopy changes. < 0 keeps the land-class height
      !   (default).
      real(r8) :: wetland_forest_htop_trop  = -1._r8  ! [m]
      ! wetland_cover_frac_site (marsh towers): vascular plant cover
      !   of a tower's footprint. The plant carbon input of the wetland tile
      !   (wetland_plant_input), root exudates included, is scaled by it, for
      !   a rewetted or managed marsh whose emergent plants leave part of the
      !   footprint as open water. < 0 keeps the whole input (default).
      real(r8) :: wetland_cover_frac_site   = -1._r8  ! [-]
      ! wetland_burial_frac_site (managed marsh towers): share of the
      !   wetland tile's plant litter input (wetland_plant_input, what is left
      !   after the root exudates) buried as new peat instead of entering the
      !   litter pools. Rewetted Delta marshes bury 280-350 g C m-2 yr-1, about
      !   four tenths of their NPP (US-Myb, US-Tw1; Arias-Ortiz et al. 2021),
      !   which a steady-state litter input would decompose. The buried carbon
      !   enters the passive pool (soil3) along the root profile with its
      !   litter N; soil3 turns over in millennia on the anoxic tile and
      !   wetland_bgc_sasu holds it. < 0 buries nothing (default).
      real(r8) :: wetland_burial_frac_site  = -1._r8  ! [-]
      ! wetland_lai_shape (shoulder seasons): with wetland_veg_glwd,
      !   scale the non-forested share's remote-sensing LAI by cap / annual
      !   peak instead of clipping it at the cap each month. The clip holds
      !   an open fen at its July LAI from April to November wherever the
      !   remote-sensing peak exceeds the cap, and with it the aerenchyma
      !   cross-section (area_tiller ~ LAI) and the tile's assimilation; the
      !   scaling keeps the measured peak and the remote-sensing seasonal
      !   shape. Off by default (clip).
      logical  :: wetland_lai_shape         = .false.
      ! wetland_forest_input_herb (boreal wetlands): with
      !   wetland_plant_input and wetland_veg_glwd, the forested share of a
      !   wetland tile feeds the methanogens as a herb layer at the LAI of the
      !   non-forested share instead of through its remote-sensing canopy. A
      !   GLWD tile carries the LAI of the dryland forest around it; tree
      !   litter is woody and falls on the surface, and tree and
      !   shrub fine roots sit mostly above the water table while sedge roots
      !   reach below it (Moore et al. 2002), so grass litter along the roots
      !   overstates the substrate the trees give the methanogens. The tile's assimilation is split between the two shares by LAI and
      !   the forested share rescaled to the non-forested LAI, so the input is
      !   the assimilation times LAI(non-forested) / LAI(tile), and all of it
      !   takes the herb belowground share wetland_bg_frac instead of the
      !   forest-weighted one; litter quality, wetland_litter_cn and
      !   wetland_exudate_frac stay as for grass. Dropping the forested share's
      !   input instead would leave a fully forested tile without fresh input,
      !   whereas forested and open peatlands of one BAWLD class emit alike
      !   (BAWLD, Olefeldt et al. 2021; BAWLD-CH4, Kuhn et al. 2021). It acts
      !   wherever a tile has a forested share, tropical peat swamp forest
      !   (GLWD 26) included, whose input falls to the cell's marsh LAI
      !   (wetland_lai_marsh) as well. The aerenchyma keeps the tile's NPP.
      !   Off by default.
      logical  :: wetland_forest_input_herb = .false.
      ! wetland_moss_input_frac: with wetland_forest_input_herb
      !   and wetland_dim_class, the forested share of a tile keeps at least
      !   its moss layer, this share of the stand input, where the ground is
      !   moss-carpeted peat (s_m, GLWD 22-25). The input ratio becomes
      !     (1 - f) r + f max(r, wetland_moss_input_frac s_m),
      !   r the herb layer ratio and f the forested share, taking the moss
      !   share inside the forested share as the tile's s_m (exact for a
      !   tower, whose shares are 0 or 1). The herb layer alone leaves a
      !   black spruce peatland, whose understory is Sphagnum and feather
      !   moss, almost without input. Mosses give 48 % of the productivity
      !   of boreal and tundra wetlands, moss over moss plus aboveground
      !   vascular NPP (53 % in fens, 58 % in bogs), against 20 % in uplands
      !   (Turetsky et al. 2010, Can. J. For. Res. 40, 1237-1264, p. 1237 and
      !   p. 1242). The share applies to the tile's input, which counts no
      !   moss (biased low) but includes the belowground NPP that the
      !   published share leaves out (biased high). The input keeps the herb
      !   belowground share.
      real(r8) :: wetland_moss_input_frac   = 0.48_r8  ! [-]
      ! wetland_anoxia_catotelm: below the water table the
      !   wetland anoxia scalar falls with depth d under the table from
      !   mino2lim to this floor, floor + (mino2lim - floor) exp(-d/efold),
      !   as in HPM (Frolking et al. 2010, eq. 9, c2 = 0.3 m). Anoxic/oxic
      !   CO2 production is 0.06-0.14 in 5-45 cm peat cores (Scanlon and
      !   Moore 2000). Negative keeps mino2lim at every depth (default).
      real(r8) :: wetland_anoxia_catotelm   = -1._r8  ! [-]
      real(r8) :: wetland_anoxia_efold      = 0.3_r8  ! [m]
      ! wetland_tau_s3: base turnover time of the wetland
      !   tile's passive pool in place of CENTURY's 222 yr, at the same
      !   reference (t_scalar = 1 near 30 C). ORCHIDEE-PEAT v2 uses a passive
      !   rate of 0.0006 per year at 30 C, i.e. 1667 yr (Qiu et al. 2019,
      !   Table 1); LPJ-WHy a slow pool of 1000 yr at 10 C. Negative keeps
      !   tau_s3 (default).
      real(r8) :: wetland_tau_s3            = -1._r8  ! [yr]
      ! Wetland N closure and burial, all off by default.
      ! wetland_n_unlimited: litter-to-SOM immobilization on the wetland tile
      !   runs at its potential rate and the N a layer lacks is supplied, as in
      !   CLM's carbon-only supplemental N (suplnitro = 'ALL'). The tile's
      !   input is prescribed from its canopy, so nothing matches it to the N
      !   supply, and the closed cascade frees the N litter needs only once
      !   soil2 has built up. f_fpi keeps the share the soil supplied.
      ! wetland_n_uptake: the N the plant input returns in its litter is taken
      !   up from the tile's mineral N left after immobilization, each layer in
      !   proportion to its residual; what the soil cannot give counts as
      !   fixed N. f_fpg is the share the soil supplied. Needs
      !   wetland_plant_input.
      ! wetland_vert_transp: CLM's SOM mixing (SoilBiogeochemLittVertTransp,
      !   bioturbation 1 cm2 yr-1) on the wetland tile, with a downward
      !   advection wetland_burial_velocity for peat accretion burying the
      !   acrotelm into the catotelm: 12 g C m-2 yr-1 (Wania et al. 2009,
      !   after Clymo 1984) over an acrotelm of 25 kg C m-3 is 4.8e-4 m yr-1.
      !   Carbon and N crossing the bottom interface leave as buried peat.
      logical  :: wetland_n_unlimited       = .false.
      logical  :: wetland_n_uptake          = .false.
      logical  :: wetland_vert_transp       = .false.
      real(r8) :: wetland_burial_velocity   = 0._r8   ! [m yr-1]
      ! wetland_bgc_sasu: semi-analytic spin-up of the wetland
      !   tile's decomposition pools. At the end of every spin-up cycle but
      !   the last, litter, CWD, soil1 and soil2 of each layer jump to the
      !   steady state of that cycle's flux-weighted rates before the N limit
      !   (fpi), -(A K)^-1 I; soil3 keeps the seeded peat. CoLM's CNSASU idea
      !   (Lu et al. 2020) for the patchtype 2 shim, which CNSASU never
      !   reaches. Acts only while DEF_simulation_time%spinup_repeat > 1.
      !   Solves each layer on its own, so it excludes wetland_vert_transp.
      !   Off by default.
      logical  :: wetland_bgc_sasu          = .false.
      ! wetland_exudate_frac: share of the wetland plant input
      !   released as root exudates. It is respired within the step in the
      !   layer where it is released, along the land class's root profile,
      !   and counts as litter heterotrophic respiration for CH4 production;
      !   it never enters the pools. Sedges feed methanogens below the water
      !   table this way (Whiting and Chanton 1992; Strom et al. 2003);
      !   LPJ-WHyMe releases 0.15 of NPP as exudates (Wania et al. 2010).
      !   0 keeps all of the input as litter (default).
      !   Retired: the exudate is taken off the litter input and also booked
      !   on decomp_hr_vr before CDecompStateUpdate, which debits it from the
      !   metabolic litter a second time (MOD_BGC_CNCStateUpdate1.F90), so
      !   each step's exudate leaves the carbon budget twice and the deep
      !   metabolic litter turns negative. root_exudate_frac books the same
      !   exudate after the pool update; a value above 0 here stops the run.
      real(r8) :: wetland_exudate_frac      = 0._r8   ! [-]
      ! root_exudate_frac: share of the net primary
      !   production of every non-lake tile released by the roots as
      !   exudates, the live plant carbon that feeds methanogens without
      !   passing through the litter first (fexu, Wania et al. 2010, Table 4
      !   p. 573: 0.15; Table 1 p. 568 lists 0.175; LPJ-WHyMe takes it from
      !   NPP each time step, p. 567). It is released along the root profile
      !   and respired within the step, counted as litter heterotrophic
      !   respiration: below the water table and in flooded subcolumns the
      !   same f_methane, temperature, pH and redox chain as the rest of the
      !   respiration turns a share of it into CH4, and the remainder leaves
      !   as CO2. The same carbon leaves the litter input, so the tile's
      !   carbon is conserved:
      !   - soil tile: NPP is the host's running mean lag_npp (e-folding
      !     nfix_timeconst), not the step's gpp - ar, whose positive part
      !     would count the day and drop the night's respiration. The host
      !     has already put this step's litterfall into the pools, so the
      !     exudate is taken from the metabolic, cellulose and lignin
      !     litter of each layer in the fine-root litter split
      !     (fr_flab/fr_fcel/fr_flig), at most half a pool per step; at
      !     steady state that is the litter input less the exudate. Litter N
      !     stays.
      !   - wetland tile: the exudate leaves the plant input
      !     (wetland_plant_input) as with wetland_exudate_frac; only where
      !     wetland_exudate_frac is 0, since it releases the same exudate there.
      !   - paddy rice with DEF_RICE_ROOT_EXUDATE > 0 keeps its own
      !     exudate; the rice share of the tile is left out.
      !   <= 0 is off (default).
      real(r8) :: root_exudate_frac         = 0._r8   ! [-]
      ! rice_aere_override: .true. gives live paddy rice the
      !   tiller geometry of get_rice_veg_proxy (porosity 0.40, radius 0.75 mm,
      !   1.0 gC per tiller), whose aerenchyma cross-section per unit tiller
      !   carbon is 1/51 of the CLM4Me default; ebullition then carries most
      !   of the paddy CH4. .false. keeps the CLM4Me crop geometry (porosity
      !   0.3, radius 2.9 mm, 0.22 gC per tiller; Riley et al. 2011, Wania et
      !   al. 2010), and the plants become the main pathway, as observed in
      !   paddies (Cicerone and Shetter 1981). .true. by default.
      logical  :: rice_aere_override        = .true.
      ! rice_aereoxid: share of the CH4 entering the
      !   aerenchyma of paddy rice that is oxidised in the rhizosphere by the
      !   O2 the roots release from inside the plant (not drawn from the soil
      !   layer's O2, as CLM4Me's aereoxid). The CLM4Me legacy switch (aereoxid with
      !   use_aereoxid_prog = .false.) acts on every column and shuts the plant
      !   O2 off, so it cannot serve here. Paddies oxidise 40% rising to 90% of
      !   the CH4 produced over the season (Cao et al. 1995, after Schutz et al.
      !   1989), 0.65 as the season mean emitted fraction 0.35 (Huang et al.
      !   1998). 0 (off) by default.
      real(r8) :: rice_aereoxid             = 0._r8   ! [-]
      ! wetland_aereoxid: share of the CH4 entering the
      !   aerenchyma of a wetland tile, or of the flooded subcolumn of a soil
      !   patch other than a paddy, that is oxidised at the roots, with the
      !   O2 the roots release (as rice_aereoxid); the tile then takes no plant
      !   O2 into its soil layers, which in the prognostic scheme fuelled the
      !   oxidation of 58-90% of production in flooded subcolumns.
      !   In situ rhizosphere oxidation of wetland plants is mostly 10-50% of
      !   the potential emission (King 1996; Calhoun and King 1997; Kankaala
      !   and Bergstrom 2004), lower for sedges (Turner et al. 2020).
      !   < 0 (off) by default.
      real(r8) :: wetland_aereoxid          = -1._r8  ! [-]
      ! rice_rox_only: with rice_aereoxid, the paddy
      !   column takes no plant O2 into its soil layers either, so root O2 is
      !   counted once, as the rhizosphere oxidation share. Off by default
      !   (off keeps both).
      logical  :: rice_rox_only             = .false.
      ! rox_on_rooted_production: what the fixed rhizosphere
      !   oxidation share x of a tile acts on (wetland_aereoxid on a wetland
      !   tile and on the flooded subcolumn of a soil patch other than a
      !   paddy, rice_aereoxid on a paddy column; x off leaves the tile alone).
      !   .false. (default): the CH4 entering the aerenchyma of each layer
      !   below the water table, as CLM4Me's aereoxid (Riley et al. 2011).
      !   Where the plant conduit is small (woody tiles, flooded forest) the
      !   share then has almost nothing to act on, and with the plant O2 kept
      !   out of the soil the roots oxidise no CH4 outside the conduit.
      !   .true.: the CH4 produced below the water table in the rooted layers,
      !   before it enters the pore water, whatever route it would leave by.
      !   Layer j loses R_j = x w_j P_j, with P_j its production and
      !     w_j = rho_j / max_k rho_k,   rho_j = rootfr_j / dz_j,
      !   the root density of the layer relative to the most densely rooted
      !   layer of the tile's root profile (the profile the aerenchyma uses).
      !   Root O2 leaks out mainly near the root tips (Laanbroek 2010, Annals
      !   of Botany 105, 141-153), so the share of a layer's pore water within
      !   its reach grows with the root density of the layer; rootfr_j itself
      !   is a share per layer and would tie the oxidation to the layer
      !   thickness, and rootfr_j renormalised over the layers below the table
      !   would raise it as the table falls. w does not
      !   depend on where the water table is, and is 0 without roots. The
      !   measured shares count the CH4 made in the root zone, or its
      !   potential emission, that the rhizosphere oxidises whatever path it
      !   takes out: in situ 10-50% for wetland plants (King 1996, mean 27%;
      !   Laanbroek 2010, Table 1), 40-90% of the CH4 produced in paddies over
      !   the season (Cao et al. 1995, after Schutz et al. 1989); excluding
      !   the roots of a tropical peat swamp forest raised the CH4 flux across
      !   the peat surface by 85-92%, by the authors' reading through the loss
      !   of root O2 (Girkin et al. 2020, ERL 15, 064013), a flux no
      !   aerenchyma carries. x is reached in
      !   the most densely rooted layer; the tile mean is x times the
      !   production-weighted w. R counts as CH4 oxidation (methane_oxid_depth,
      !   and CO2), its O2 comes from inside the plant and is not taken from
      !   the layer, so it needs a live canopy as the aerenchyma does (lai >
      !   0; none in a harvested paddy), and the plant O2 stays out of the
      !   soil layers as under wetland_aereoxid and rice_rox_only. Needs
      !   wetland_aereoxid >= 0 or rice_aereoxid > 0. Bit-for-bit unchanged
      !   when off.
      logical  :: rox_on_rooted_production  = .false.
      ! rice_host_water: .true. floods the paddy column when
      !   the host holds water on the field (surface water >= 1 mm or the
      !   water table at the surface), i.e. from irrigation and the bunds of
      !   DEF_USE_IRRIGATION / DEF_PADDY_RICE_BUND, instead of the imposed
      !   floor rice_paddy_min_finundated with its midseason drainage and
      !   post-harvest decay. Rainfed rice then floods only when rain fills
      !   the bunds. .false. by default.
      logical  :: rice_host_water           = .false.
      ! lake_sed_o2_demand: .true. lets respiration of the
      !   non-CH4 carbon in lake sediment consume O2 where O2 is present, as in
      !   soils; the sediment below the oxic surface film then turns anoxic.
      !   .false. keeps no O2 demand, so the 1 mol m-3 O2 every lake layer is
      !   allocated with stays for years and oxidises the CH4 made at depth
      !   (DE-Dgw: oxidation about equal to production after nine years).
      logical  :: lake_sed_o2_demand        = .false.
      ! lake_prod_efold: e-folding depth [m] of lake sediment
      !   decomposition below the bed, the fresh deposit being the reactive
      !   part (Middelburg 1989, doi:10.1016/0016-7037(89)90239-1); column
      !   total unchanged. <= 0 keeps the carbon profile (CLM4Me).
      real(r8) :: lake_prod_efold           = -1._r8  ! [m]
      ! floodplain_sif: .true. gives the flooded subcolumn of
      !   a floodplain the CLM4Me seasonal inundation factor instead of
      !   mino2lim: full rate for the annually inundated share, mino2lim for
      !   the seasonal excess (Riley et al. 2011, doi:10.5194/bg-8-1925-2011).
      logical  :: floodplain_sif            = .false.
      ! flood_saturates_soil: .true. fills the pores of the
      !   flooded subcolumn of a soil patch (routing or site flood) with
      !   water for the gas diffusion and oxygen of the CH4 column, instead
      !   of lending it only what the host soil holds: the host soil does not
      !   know it is flooded, the flood water sits in the river storage, and
      !   a column saturated to a third lets air-filled pores carry oxygen
      !   down (BW-Nxr: O2 1-1.6 mol m-3 under flood, wetland tiles 0.02-0.09).
      !   The host water balance is untouched.
      logical  :: flood_saturates_soil      = .false.
      ! lake_strat_drho: bottom-minus-surface water density
      !   [kg m-3] above which the lake counts as stratified and its CH4 node
      !   exchanges no gas with the air (stored CH4 leaves at overturn).
      !   <= 0: always exchanging when ice-free.
      real(r8) :: lake_strat_drho           = -1._r8  ! [kg m-3]
      ! lake_bubble_dissol_depth: e-folding rise distance h_b
      !   [m] over which a bubble from the lake bed loses its CH4 to the water.
      !   The lake's ebullition reaches the air as the share
      !   sum_k dA_k exp(-z_k/h_b) / sum_k dA_k of the bed area dA_k at depth
      !   z_k, the bed following the valley-shaped basin of FLaMe (shape 2,
      !   maximum depth twice the lake depth; Maisonnier et al. 2025); the
      !   rest dissolves into the water node, where it is oxidised or leaves
      !   by gas exchange like the CH4 the sediment releases by diffusion.
      !   LAKE 2.0 brings 68-70% of the bubble CH4 leaving 12.5 m to the
      !   surface (Stepanenko et al. 2016), h_b about 34 m; a 6 mm bubble
      !   keeps about 30% of its CH4 over 23 m (McGinnis et al. 2006), h_b
      !   about 19 m. The ebullition threshold keeps the full water head.
      !   <= 0: all lake ebullition goes to the air (default).
      real(r8) :: lake_bubble_dissol_depth  = -1._r8  ! [m]
      ! lake_k_m: Michaelis-Menten constant [mol m-3] for CH4
      !   in the oxidation of the lake water node, which otherwise takes the
      !   soil value k_m (5e-3). LAKE 2.0 calibrates 3.75e-2 for lake water
      !   (Stepanenko et al. 2016); hypolimnion water of two Alaskan lakes
      !   gives 4.5e-3 and 1.1e-2 (Lofton et al. 2014). The maximum rate is
      !   set by lake_vmax_methane_oxid. <= 0 keeps k_m (default).
      real(r8) :: lake_k_m                  = -1._r8  ! [mol m-3]
      ! lake_sod20: sediment O2 demand at 20 C [g O2 m-2 d-1]
      !   drawn from the lake water node, lake_sod20 * 1.065**(T_bot - 20 C)
      !   * O2 / (K_O2 + O2), with T_bot the temperature of the deepest lake
      !   layer and K_O2 the O2 constant of the water-node CH4 oxidation
      !   (lake_k_m_o2, else k_m_o2). Guo et al. (2023) take 1.5-9.9 g O2
      !   m-2 d-1 and 1.065 for a tropical floodplain lake; LAKE 2.0 applies
      !   the demand after Walker and Snodgrass (1986). With the stratification
      !   gate (lake_strat_drho) it lets the water cut off from the air
      !   turn anoxic. The O2 taken is counted in lake_sed_o2_flux.
      !   <= 0: no sediment O2 demand on the water (default).
      real(r8) :: lake_sod20                = -1._r8  ! [g O2 m-2 d-1]
      ! lake_icebubble_release: rate [d-1] at which the CH4 of
      !   the bubbles held under lake ice leaves once the ice fraction of the
      !   top lake layer falls to 0.1 or below. While it is above 0.1 the lake
      !   bed keeps bubbling and the bubbles that reach the surface go into the
      !   patch store lake_icebubble_ch4_stock instead of the air. ALBM
      !   releases ice-trapped bubbles at 0.14-1 d-1 (Tan et al. 2024);
      !   bLake4Me releases 60% of them at ice-out (Tan et al. 2015), LAKE2.6
      !   traps 90% and releases 67.5% at ice break (Li et al. 2026).
      !   <= 0: the lake bubbles only when the ice fraction is 0.1 or below,
      !   nothing is stored (default).
      real(r8) :: lake_icebubble_release    = -1._r8  ! [d-1]
      ! lake_icebubble_dissol: share [0-1] of the CH4 released
      !   from the under-ice bubble store that dissolves into the lake water
      !   node, where it is oxidised or leaves by gas exchange like the CH4
      !   the sediment releases by diffusion; the rest joins the ebullition to
      !   the air. 0.4 is the share bLake4Me does not release (Tan et al.
      !   2015); Johnson et al. (2022) oxidise 75% of the CH4 stored under ice
      !   at ice-out. Used only with lake_icebubble_release > 0.
      real(r8) :: lake_icebubble_dissol     = 0.4_r8  ! [-]
      ! lake_active_pool: lake CH4 production is fed
      !   by an active pool of degradable organic carbon C_act [g C m-2] in
      !   each lake patch, supplied at a constant rate S and mineralized at
      !   first order,
      !     dC_act/dt = S - k(T) C_act,   k(T) = lake_k20 lake_theta**(T - 20),
      !   T in C: the autochthonous pool of FLaMe v1.0 (Maisonnier et al.
      !   2025, ESD 16, 1779-1808, doi:10.5194/esd-16-1779-2025, eqs. 2, 7
      !   and 10). S is the degradable deposition, the sediment supply net of
      !   burial, so FLaMe's burial term (eq. 8, k_bur = k20/2) is left out.
      !   The pool is spread over the sediment layers as the present
      !   production is: exp(-z/lake_prod_efold) when that is set,
      !   else the lake_soilc carbon profile. Layer j mineralizes its share at
      !   k(T_j), times the 1 K freezing ramp of lake sediment; CH4 is
      !   f_methane (lake_f_methane when set) of that (at most 0.5; a
      !   share of the whole mineralization,
      !   as FLaMe f_mm 1/4, 1/6-1/2) and the rest leaves as CO2 on the present path,
      !   with its O2 demand under lake_sed_o2_demand. At steady state
      !   C_act = S / k and production is f_methane S whatever the
      !   temperature, which sets only the pool size and the seasonal timing.
      !   The whole-column lake_soilc path (lake_decomp_fact, q10lake,
      !   q10lakebase, cnscalefactor) then makes no CH4 and lake_soilc is not
      !   debited. The pool starts empty (restart field ch4_lake_cact, 0 in a
      !   restart written before it); 1/k is 107 d at 28 C and 171 d at 4 C.
      !   Needs allowlakeprod. Off by default; bit-for-bit unchanged when off.
      logical  :: lake_active_pool          = .false.
      ! lake_cdep_band: S, the degradable organic
      !   carbon supply to the sediment of a lake patch [g C m-2 yr-1], by the
      !   patch latitude: (1) north of 60 N, (2) 45-60 N, (3) 23.5-45 N,
      !   (4) 23.5 S-23.5 N, (5) south of 23.5 S. S is the sediment
      !   mineralization that balances the organic carbon burial B of natural
      !   lakes at burial efficiency BE = 1/3, S = B (1 - BE) / BE = 2 B (the
      !   BE implied by FLaMe; measured BE mean 48% over 27 sediment sites,
      !   22% where autochthonous carbon dominates, Sobek et al. 2009, L&O 54,
      !   2243-2254, doi:10.4319/lo.2009.54.6.2243). B is the median burial
      !   of natural lakes by latitude band (Mendonca et al. 2017, Nat.
      !   Commun. 8, 1694, doi:10.1038/s41467-017-01789-6, Supplementary
      !   Data 1): tropics 29.6, 23.3-40 degrees 17.6 (taken for bands 3 and
      !   5), north of 60 N 10 between 2.1 north of 66.3 N and 27.9 at
      !   55-66.3 N; 45-60 N takes 15, undisturbed boreal lakes (Heathcote et
      !   al. 2015, Nat. Commun. 6, 10016, doi:10.1038/ncomms10016), not the
      !   median 34 of the mostly agricultural 40-55 degree lakes. The
      !   measured sediment mineralization has median 40 and range 5.3-196
      !   g C m-2 yr-1 (Sobek et al. 2009).
      real(r8) :: lake_cdep_band(5)         = [20._r8, 30._r8, 35._r8, 59._r8, 35._r8]  ! [g C m-2 yr-1]
      ! lake_k20: first-order mineralization rate of
      !   the active pool at 20 C [d-1]; FLaMe k20 0.008 (0.003-0.015) d-1
      !   (Maisonnier et al. 2025, Table 1). Used with lake_active_pool.
      real(r8) :: lake_k20                  = 0.008_r8  ! [d-1]
      ! lake_theta: temperature coefficient of k;
      !   FLaMe 1.02 (Maisonnier et al. 2025, Table 1), a Q10 of 1.22; the
      !   measured Q10 of lake sediment mineralization, 2.1-2.3 (Gudasz et
      !   al. 2010, Nature 466, 478-481, doi:10.1038/nature09186), is theta
      !   1.077-1.087. Used with lake_active_pool.
      real(r8) :: lake_theta                = 1.02_r8   ! [-]
      ! lake_f_methane: share of the mineralization of the
      !   active pool (lake_active_pool) that becomes CH4 below the lake water
      !   table, in place of f_methane for the lake patches the pool feeds, so
      !   that steady-state production is lake_f_methane S (lake_cdep_band).
      !   The pool is debited by the whole mineralization whatever the split;
      !   the rest leaves as CO2 as before. CH4 is 20-56% of whole-lake carbon
      !   mineralization (Bastviken 2009, cited in Schenk et al. 2021); fresh
      !   algal and leaf matter incubated anoxic at 20-22 C gives CH4/(CH4 +
      !   CO2) of 0.32-0.39 (Grasset et al. 2018, L&O 63, 1488-1501,
      !   doi:10.1002/lno.10786); anoxic stratified lakes release 47 +- 9% of
      !   their organic input as CH4 (Kelly and Chynoweth 1981, L&O 26,
      !   891-897, doi:10.4319/lo.1981.26.5.0891), eutrophic Lake 227 55%
      !   (Rudd and Hamilton 1978, L&O 23, 337-348,
      !   doi:10.4319/lo.1978.23.2.0337); FLaMe takes 1/4 (1/6-1/2;
      !   Maisonnier et al. 2025). At most 0.5 (2 CH2O -> CH4 + CO2), as
      !   f_methane. Needs lake_active_pool; < 0 keeps f_methane (default).
      real(r8) :: lake_f_methane            = -1._r8  ! [-]
      ! wetland_salinity_site: porewater salinity [psu] of a
      !   tower's wetland; the yield is scaled by 10**(-0.056 S) (Poffenbarger
      !   et al. 2011, doi:10.1007/s13157-011-0197-0). <= 0: fresh, no scaling.
      real(r8) :: wetland_salinity_site     = -1._r8  ! [psu]
      ! wetland_nitrate_site: nitrate of the water feeding a
      !   tower's tidal wetland [mg L-1], the creek or river annual mean.
      !   Nitrate reducers take the electron donors of methanogens; the yield
      !   is scaled by K / (K + NO3), K = wetland_nitrate_ki, as in PEPRMT-Tidal
      !   (Oikawa et al. 2024, eq. 5, doi:10.1029/2023JG007943; K of the v1.0
      !   code, doi:10.5281/zenodo.10278505). <= 0: no scaling.
      real(r8) :: wetland_nitrate_site      = -1._r8  ! [mg L-1]
      real(r8) :: wetland_nitrate_ki        = 0.102_r8  ! [mg L-1]
      ! tundra_high_share_site (site diagnosis only):
      !   share of a tower's wetland tile taken by high microsites (tussocks,
      !   polygon rims, high centres). They respire as usual but make almost no
      !   CH4: ecosys ridges emit 0.0-0.1 against 4-6 g C m-2 yr-1 in centres
      !   and troughs (abstract) while their heterotrophic respiration,
      !   132-151, is no lower than there, 103-147 g C m-2 yr-1 (Table 4;
      !   Grant et al. 2017, JGR Biogeosci. 122, 3174-3187,
      !   doi:10.1002/2017JG004037), and
      !   moss-lichen plots emit 0.00 against 1.68 mg CH4-C m-2 h-1 in wet
      !   sedge (Davidson et al. 2016, Ecosystems 19, 1116-1132,
      !   doi:10.1007/s10021-016-9991-0). The high share makes no CH4, so it has
      !   none to oxidise or release; the tile's two subcolumns stand for the
      !   low remainder. Single point only (no global map of high microsites);
      !   <= 0 is off (default).
      real(r8) :: tundra_high_share_site    = -1._r8  ! [-]
      ! wetland_moss_rss: the moss and peat surface of a
      !   dynamic wetland tile resists evaporation, so its snow-free soil
      !   surface resistance is held at least at that of the surface. Sphagnum
      !   carpets keep about 100 s m-1 with the water table near the surface
      !   (Kim and Verma 1996; Kettridge et al. 2013), rising as it falls below
      !   about 0.1 m, to near 500 s m-1 at 0.35 m in Finnish and Swedish mires
      !   (Alekseychik et al. 2018). The resistance is weighted by the tile's
      !   peatland share; tropical classes and ponded tiles take none.
      !   Off by default.
      logical  :: wetland_moss_rss          = .false.
      ! wetland_root_efold: e-folding depth [m] of the root
      !   profile of a wetland tile, used for its plant carbon input (root
      !   litter, exudates) and for the CH4 root fraction (aerenchyma, root
      !   respiration); root water uptake keeps the land-class profile. IGBP
      !   wetland takes the grassland profile of Zeng (2001), 52% of roots in
      !   the top 9 cm; sedge fens and marshes root deeper (Jackson beta
      !   0.955-0.976 and 0.94-0.97 in measured profiles), and LPJ-WHyMe
      !   gives flood-tolerant graminoids exp(-z/0.2517 m) (Wania et al.
      !   2010). A value <= 0 keeps the land-class profile (default).
      real(r8) :: wetland_root_efold        = -1._r8  ! [m]
      ! wetland_f_tropical_peat: f_methane of the peatland
      !   share of a wetland tile within 23.5 degrees of the equator.
      !   Anaerobic incubations of tropical forest peat give CH4 as 0.2-0.3%
      !   of the decomposed carbon at 25 C (Girkin et al. 2020), and the
      !   Maludam tower emits about 0.5% of half its ecosystem respiration
      !   (Wong et al. 2018), against f_methane 0.2, which has no measured
      !   source (Riley et al. 2011, Table 1). A negative value keeps
      !   f_methane (default).
      real(r8) :: wetland_f_tropical_peat   = -1._r8  ! [-]
      ! wetland_f_tree_ratio: the forested share of a
      !   wetland tile makes CH4 at this fraction of f_methane. Tree-
      !   covered towers emit a median 1.9% of half their ecosystem
      !   respiration against 6.4% at open ones (FLUXNET-CH4 TREE flag), while
      !   the model's respiration there is about right; cutting their plant
      !   transport instead raised emission (less plant O2, more ebullition).
      !   Anaerobic incubations give a maximum CH4 production of 2.4 in
      !   tree-dominated against 19 ug C gC-1 d-1 in graminoid-dominated
      !   soils (Treat et al. 2015, Fig. 3a). A negative value keeps f_methane
      !   (default).
      real(r8) :: wetland_f_tree_ratio      = -1._r8  ! [-]
      ! woody_conduit_area: cross-section [m2 m-2] of the gas
      !   conduit of woody plants standing in water, for the forested share
      !   f_w of a wetland tile (wetland_veg_glwd), in both subcolumns, and for
      !   all of the flooded subcolumn of a woody soil patch (the non-grass
      !   classes of methane_patch_is_nongrass: flooded forest). Below
      !   the water table the plant conduit of the subcolumn is then
      !     A = (1 - f_w) A_tiller + f_w A_w,
      !   with A_tiller the CLM4Me tiller cross-section of the herbs at the
      !   grass porosity (NPP x belowground share x LAI / tiller_C x
      !   poros_tiller x pi aere_radius**2) and A_w this key, and CH4 and O2
      !   diffuse through A as through the tillers (Riley et al. 2011). A_w is
      !   fixed: no NPP, LAI, phenology, scale_factor_aere or unsat_aere_ratio,
      !   so the conduit stays open at night and without leaves. Wetland tree
      !   stems vent the CH4 of the soil pore water, mostly by gas diffusion
      !   out of the lenticels of the lower stem, unrelated to leaf area and
      !   transpiration, in winter as in summer (Pangala et al. 2013, 2014,
      !   2015, 2017; Jeffrey et al. 2024), while sedge tillers put on a peat
      !   swamp forest give it a ground conductance 2-4 orders of magnitude
      !   above the stem fluxes measured there (Pangala et al. 2013). Of the
      !   global models only JULES gives trees a conduit of their own (Gedney
      !   et al. 2019). Stem fluxes over pore water CH4 put A_w at 1e-6 to
      !   4e-5 for tropical peat and temperate swamp forest and 2e-5 to 1e-3
      !   for Amazon floodplain forest. The rhizosphere oxidation share
      !   (wetland_aereoxid) acts on this conduit as on the tillers.
      !   < 0 (off) keeps the tillers on the forested share (default).
      real(r8) :: woody_conduit_area        = -1._r8  ! [m2 m-2]
      ! pft_grass_share: a gridded PFT or PC build merges
      !   every natural soil type into patchclass 1 (MOD_LandPatch.F90), so
      !   methane_patch_is_nongrass calls every natural soil patch non-grass:
      !   all of them get the non-grass aerenchyma porosity and, with
      !   woody_conduit_area, the woody conduit on all of their flooded
      !   subcolumn. With this key a soil patch takes both from its own PFT
      !   make-up instead: the tiller porosity is
      !     poros_tiller x (g + (1 - g) x nongrassporosratio),
      !   with g the share of grasses (PFT 12-14) and crops (CLM4Me's grass
      !   test, Riley et al. 2011) in the vegetated cover (bare ground left
      !   out: it has no tillers), and the woody share of the flooded
      !   subcolumn (woody_conduit_area) is the tree cover (PFT 1-8) over
      !   the whole patch (bare ground included: the conduit is a cross-
      !   section per unit ground under trees; the stem emission behind
      !   it is measured on flooded trees, Pangala et al. 2017, not shrubs).
      !   There the tillers are the herbs and keep the grass porosity, as on
      !   the non-forested share of a wetland tile. Single-point builds keep
      !   the tower's IGBP class. Off (default) keeps the class test.
      logical  :: pft_grass_share           = .false.
      ! pft_root_profile: the methane module takes a
      !   patch's root profile from its land class (MOD_Const_LC rootfr), so
      !   in a gridded PFT or PC build every natural soil patch (patchclass 1)
      !   gets the evergreen needleleaf forest profile. With this key a patch
      !   with PFTs takes the pftfrac-weighted mean, over its vegetated PFTs,
      !   of the host's PFT profiles (MOD_Const_PFT rootfr_p, those of the
      !   host's root water uptake). The profile weights the layers of the
      !   aerenchyma conductance and of the rhizosphere oxidation.
      !   Single-point builds keep the tower's class. Off (default) keeps the
      !   class profile.
      logical  :: pft_root_profile          = .false.

      ! R2 short-term SOC fix (methane-only): paddy soils accumulate SOC
      ! ~2-3x faster than upland soils under long flooding (Pan 2010 GCB,
      ! Inubushi 2003), but CoLM reads soil C profile as a grid-mean from
      ! cnsteadystate.nc with no paddy-specific value.  This systematically
      ! under-estimates the methanogenic substrate available on rice patches.
      ! At the same time CoLM assumes 100% straw return (stems all to litter
      ! at harvest), which slightly over-estimates substrate input; net of
      ! the two is ~ -30% to -50% on rice CH4 production.
      !
      ! Reserved compatibility parameter.  The former methane-only production
      ! multiplier was not carbon conservative because it did not debit the
      ! BGC substrate pool.  Validation therefore requires exactly 1 until a
      ! coupled carbon-debit path is implemented.
      real(r8) :: rice_substrate_boost          = 1.0_r8
   END type Methane_type

   type Methane_hydrology_type
      ! vdcf and pc remain in the namelist only for backward compatibility.
      ! The current CLM fill-and-spill formulation uses slopebeta/slopemax;
      ! non-default values for the retired knobs are rejected below.
      real(r8) :: vdcf = 2._r8
      real(r8) :: slopebeta = -3._r8
      real(r8) :: slopemax = 0.4_r8
      real(r8) :: pc = 0.4_r8
   END type Methane_hydrology_type

   type (Methane_type) :: DEF_METHANE
   type (Methane_hydrology_type) :: DEF_METHANE_hydrology

   real(r8), save :: atm_ch4_file_molmol(0:3000,12)
   logical,  save :: atm_ch4_file_loaded = .false.

CONTAINS

	   SUBROUTINE read_methane_namelist (nlfile)

   USE MOD_Namelist
   IMPLICIT NONE

   character(len=*), intent(in) :: nlfile
   ! Local variables
   logical :: fexists
   integer :: ierr
   character(len=512) :: iomsg
   integer :: unit_nml

   namelist /nl_colm_methane_parameter/ DEF_METHANE,DEF_METHANE_hydrology

      ! A model instance may be finalized and initialized again in the same
      ! process.  Reconstruct both parameter objects so fields omitted by the
      ! next namelist cannot inherit overrides from the previous instance.
      DEF_METHANE = Methane_type()
      DEF_METHANE_hydrology = Methane_hydrology_type()
      ! The forcing filename and contents are namelist state, so a new read
      ! must also invalidate the saved file cache even when the path is
      ! unchanged.
      atm_ch4_file_loaded = .false.
      atm_ch4_file_molmol(:,:) = -1._r8

      ! Read on every rank so all workers use the same methane constants.
      ! The previous master-only read left non-master ranks at default values
      ! whenever users overrode DEF_METHANE in the parameter namelist.
      INQUIRE (file=trim(nlfile), exist=fexists)
      IF (.not. fexists) THEN
         CALL CoLM_Stop (' ***** ERROR: methane parameter file does not exist: '// trim(nlfile))
      ENDIF
      open(newunit=unit_nml, status='OLD', file=trim(nlfile), form="FORMATTED")
      iomsg = ''
      read(unit_nml, nml=nl_colm_methane_parameter, iostat=ierr, iomsg=iomsg)
      IF (ierr /= 0) THEN
         close(unit_nml)
         IF (p_is_master) write(*,'(A,A,A)') &
            'ERROR read_methane_namelist: invalid &nl_colm_methane_parameter in ', trim(nlfile), ': '
         IF (p_is_master) write(*,'(A)') trim(iomsg)
         CALL CoLM_Stop (' ***** ERROR: Problem reading namelist: '// trim(nlfile))
      ENDIF
      close(unit_nml)

      CALL validate_methane_namelist ()

	   END SUBROUTINE read_methane_namelist

	   SUBROUTINE configure_methane_inundation_mode ()
	      ! Resolve the user-facing four-option CH4 inundation mode into the
	      ! internal scheme integer plus paired methane switches.  The old
	      ! integer scheme remains only as the physics dispatch key; users
	      ! should set DEF_METHANE%inundation_mode to one of:
      !   wetwat, satellite/giems, routing, dynamic_wtd, hybrid, colm.
	      USE MOD_Namelist, only: DEF_wetland_finundation_scheme, &
	                              DEF_USE_Dynamic_Wetland
	      IMPLICIT NONE

	      character(len=32) :: mode

	      mode = tracer_lower(adjustl(trim(DEF_METHANE%inundation_mode)))

	      ! Reset only hydrology-mode-owned switches before dispatch so repeated
	      ! calls cannot retain a previous mode.  Biome production, redox,
	      ! vertical profile and numeric flood thresholds remain independent
	      ! namelist choices.
	      DEF_METHANE%enable_wetwat_finundated_override = .false.
	      DEF_METHANE%wetland_dry_unsat_branch = .true.
	      DEF_METHANE%use_routing_for_soil = .false.

	      SELECT CASE (trim(mode))
       CASE ('saturated')
         ! Permanently saturated wetland: scheme 1's base assigns finundated=1
         ! to every wetland tile (the original CoLM "wetland at surface water
         ! table" behaviour). Unlike 'wetwat' the wetwat/wetwatmax override is
         ! left off, so a single-point wetland with no lateral inflow (where the
         ! surface store drains under ET) stays inundated instead of collapsing
         ! to finundated=0. Requires DEF_USE_Dynamic_Wetland=.false. so the soil
         ! column is held saturated (zwt=0) consistently.
         DEF_wetland_finundation_scheme = 1
         DEF_METHANE%enable_wetwat_finundated_override = .false.
         DEF_METHANE%wetland_dry_unsat_branch = .true.
         IF (DEF_USE_Dynamic_Wetland) THEN
            IF (p_is_master) write(6,*) &
               '***** ERROR: saturated methane inundation mode requires DEF_USE_Dynamic_Wetland = .false.'
            CALL CoLM_Stop (' ***** ERROR: invalid methane inundation mode / dynamic wetland combination')
         ENDIF

	      CASE ('wetwat')
	         DEF_wetland_finundation_scheme = 1
	         DEF_METHANE%enable_wetwat_finundated_override = .true.
	         DEF_METHANE%wetland_dry_unsat_branch = .true.
	         IF (DEF_USE_Dynamic_Wetland) THEN
	            IF (p_is_master) write(6,*) &
	               '***** ERROR: wetwat methane inundation mode requires DEF_USE_Dynamic_Wetland = .false.'
	            CALL CoLM_Stop (' ***** ERROR: invalid methane inundation mode / dynamic wetland combination')
	         ENDIF

	      CASE ('satellite','giems')
	         DEF_wetland_finundation_scheme = 5
	         DEF_METHANE%enable_wetwat_finundated_override = .false.
	         DEF_METHANE%wetland_dry_unsat_branch = .true.
	         IF (DEF_USE_Dynamic_Wetland) THEN
	            IF (p_is_master) write(6,*) &
	               '***** ERROR: satellite methane inundation mode requires DEF_USE_Dynamic_Wetland = .false.'
	            CALL CoLM_Stop (' ***** ERROR: invalid methane inundation mode / dynamic wetland combination')
	         ENDIF

       CASE ('colm')
         ! CoLM mode.  The mapped wetland tile is permanent
         ! wetland (finundated 1, host wetland bucket); every other tile
         ! takes the routing flood fraction and depth, wetland first, gated
         ! by hybrid_soil_threshold.  Without routing -- a single point --
         ! the published fields stay zero, so the soil tile is a plain soil
         ! column and only a rain-fed wetland can be represented.  With
         ! USE_SITE_WTD the observed table takes over: finundated is 1 when
         ! it is at or above the surface and 0 otherwise, so one column
         ! carries the patch.  The unsaturated branch is live (no forced
         ! dry column) and no sigmoid is applied.
         DEF_wetland_finundation_scheme = 8
         DEF_METHANE%enable_wetwat_finundated_override = .false.
         DEF_METHANE%wetland_dry_unsat_branch = .false.
         ! Soil tiles flooded by routing are floodplains (biome parameters
         ! and the floodplain area diagnostic), as in the hybrid mode.
         DEF_METHANE%use_routing_for_soil = .true.
         ! With DEF_USE_Dynamic_Wetland the wetland tile keeps its own water
         ! table and the scheme 8 branch follows it (Physics).

	      CASE ('routing')
	         DEF_wetland_finundation_scheme = 7
	         DEF_METHANE%enable_wetwat_finundated_override = .false.
	         DEF_METHANE%wetland_dry_unsat_branch = .true.
	         IF (DEF_USE_Dynamic_Wetland) THEN
	            IF (p_is_master) write(6,*) &
	               '***** ERROR: routing methane inundation mode requires DEF_USE_Dynamic_Wetland = .false.'
	            CALL CoLM_Stop (' ***** ERROR: invalid methane inundation mode / dynamic wetland combination')
	         ENDIF

	      CASE ('dynamic_wtd','dynamic-wtd')
	         DEF_wetland_finundation_scheme = 6
	         DEF_METHANE%enable_wetwat_finundated_override = .false.
	         ! Dynamic wetland hydrology supplies the WTD forcing; keep the
	         ! dry unsaturated branch active for wetland tiles.
	         DEF_METHANE%wetland_dry_unsat_branch = .true.
	         DEF_METHANE%use_routing_for_soil = .false.
	         IF (.not. DEF_USE_Dynamic_Wetland) THEN
	            IF (p_is_master) write(6,*) &
	               '***** ERROR: dynamic_wtd requires DEF_USE_Dynamic_Wetland = .true.'
	            CALL CoLM_Stop (' ***** ERROR: invalid methane inundation mode / dynamic wetland combination')
	         ENDIF

	      CASE ('hybrid','dh_all_thr05','dyn_routing_hybrid')
		         ! Site-calibrated hybrid mode; not a globally validated default.
	         ! Combines routing and dynamic-WTD hydrology.  Biome yield,
	         ! redox lag and vertical source attenuation are independent
	         ! namelist controls; the standard parameter file enables their
	         ! recommended values explicitly.
	         ! Numeric values below are author-selected midpoints, not direct quotes:
	         !   - dyn_routing_hybrid: wetland sigmoid(zwt) + soil routing fldfrc
	         !   - biome f_methane lookup (range from Bridgham 2013 GCB review)
	         !   - biome redoxlag lookup (Pangala 2017 / Whalen 1990 qualitative)
	         !   - hybrid soil threshold 0.05 (empirical, tuned to Pantanal P/T)
	         !   - depth attenuation z0=0.30m (within Walter & Heimann 2001 range)
	         ! Tuned to match Pantanal in-situ flux (Marani & Alvalá 2007).
	         ! Names dyn_routing_hybrid/hybrid kept as backwards-compatible aliases.
	         DEF_wetland_finundation_scheme              = 6
	         DEF_METHANE%enable_wetwat_finundated_override = .false.
	         DEF_METHANE%wetland_dry_unsat_branch        = .true.
	         DEF_METHANE%use_routing_for_soil            = .true.
	         IF (.not. DEF_USE_Dynamic_Wetland) THEN
	            IF (p_is_master) write(6,*) &
	               '***** ERROR: hybrid mode requires DEF_USE_Dynamic_Wetland = .true.'
	            CALL CoLM_Stop (' ***** ERROR: invalid methane inundation mode / dynamic wetland combination')
	         ENDIF

	      CASE DEFAULT
	         IF (p_is_master) write(6,*) &
	            '***** ERROR: unsupported DEF_METHANE%inundation_mode = ', trim(DEF_METHANE%inundation_mode), &
            '; expected saturated, wetwat, satellite, routing, dynamic_wtd, hybrid, or colm.'
	         CALL CoLM_Stop (' ***** ERROR: unsupported methane inundation mode')
	      END SELECT

	      IF (p_is_master) write(6,'(A,A,A,I0,A,L1,A,L1)') &
	         ' CH4 inundation mode: ', trim(mode), &
	         ' -> scheme=', DEF_wetland_finundation_scheme, &
	         ' wetwat_override=', DEF_METHANE%enable_wetwat_finundated_override, &
	         ' dry_unsat=', DEF_METHANE%wetland_dry_unsat_branch
	   END SUBROUTINE configure_methane_inundation_mode

	   SUBROUTINE validate_methane_namelist ()
      ! Range-check user-overridable methane parameters.  Catches negative
      ! production rates, non-positive Q10 / Michaelis-Menten constants, and
      ! pH window inversions that would otherwise propagate as silent NaNs or
      ! negative fluxes.
	      IMPLICIT NONE
      logical :: bad

      bad = .false.

      ! Ordered comparisons do not reject NaN: every comparison with NaN is
      ! false.  Check the complete real-valued namelist surface once before
      ! the semantic range checks below so no user override can silently
      ! inject a non-finite rate, scale, threshold, or temperature.
      IF (.not. all(ieee_is_finite([ &
         DEF_METHANE%q10methane, DEF_METHANE%f_methane, DEF_METHANE%f_methane_tropical_peat, &
         DEF_METHANE%f_methane_tropical_floodplain, DEF_METHANE%f_methane_floodplain, &
         DEF_METHANE%f_methane_temperate_marsh, DEF_METHANE%f_methane_boreal_fen, &
         DEF_METHANE%f_methane_boreal_bog, DEF_METHANE%f_methane_rice_paddy, &
         DEF_METHANE%f_methane_upland_soil, DEF_METHANE%redoxlag_tropical_peat, &
         DEF_METHANE%redoxlag_tropical_floodplain, DEF_METHANE%redoxlag_temperate_marsh, &
         DEF_METHANE%redoxlag_boreal_fen, DEF_METHANE%redoxlag_boreal_bog, &
         DEF_METHANE%redoxlag_rice_paddy, DEF_METHANE%redoxlag_upland_soil, &
         DEF_METHANE%redoxlag_wetland_dim, DEF_METHANE%wetland_moss_input_frac, &
         DEF_METHANE%z0_methane_prod, DEF_METHANE%vmax_methane_oxid, &
         DEF_METHANE%vmax_oxid_unsat, DEF_METHANE%k_m, DEF_METHANE%k_m_unsat, &
         DEF_METHANE%k_m_o2, DEF_METHANE%q10_methane_oxid, DEF_METHANE%lake_oxid_scale, &
         DEF_METHANE%lake_k_m_o2, DEF_METHANE%lake_vmax_methane_oxid, &
         DEF_METHANE%lake_oxic_sediment_depth, DEF_METHANE%lake_bubble_dissol_depth, &
         DEF_METHANE%lake_k_m, DEF_METHANE%lake_sod20, DEF_METHANE%lake_icebubble_release, &
         DEF_METHANE%lake_icebubble_dissol, DEF_METHANE%lake_cdep_band, &
         DEF_METHANE%lake_k20, DEF_METHANE%lake_theta, DEF_METHANE%lake_f_methane, &
         DEF_METHANE%B_init_methanogen, &
         DEF_METHANE%B_init_methanotroph, DEF_METHANE%B_min_methanogen, &
         DEF_METHANE%B_min_methanotroph, DEF_METHANE%B_max_fraction_methanogen, &
         DEF_METHANE%B_max_fraction_methanotroph, DEF_METHANE%mu_max_methanogen, &
         DEF_METHANE%mu_max_methanotroph, DEF_METHANE%gamma_methanogen, &
         DEF_METHANE%gamma_methanotroph, DEF_METHANE%gamma_microbial_dormant, &
         DEF_METHANE%gamma_microbial_freeze, DEF_METHANE%K_substrate_methanogen_pool, &
         DEF_METHANE%K_inh_O2_methanogen, DEF_METHANE%kappa_m_methanogen, &
         DEF_METHANE%kappa_m_methanotroph, DEF_METHANE%max_microbe_prod_multiplier, &
         DEF_METHANE%q10_microbe_growth, DEF_METHANE%T_ref_microbe, &
         DEF_METHANE%dormancy_rate_active, DEF_METHANE%dormancy_rate_revive, &
         DEF_METHANE%dormancy_threshold_methanogen_fS, &
         DEF_METHANE%dormancy_threshold_methanogen_fO2, &
         DEF_METHANE%dormancy_threshold_methanotroph_fS, &
         DEF_METHANE%dormancy_threshold_methanotroph_fO2, DEF_METHANE%vgc_max, &
         DEF_METHANE%nongrassporosratio, DEF_METHANE%poros_tiller, &
         DEF_METHANE%unsat_aere_ratio, DEF_METHANE%porosmin, DEF_METHANE%aere_radius, &
         DEF_METHANE%rob, DEF_METHANE%scale_factor_aere, DEF_METHANE%scale_factor_gasdiff, &
         DEF_METHANE%scale_factor_liqdiff, DEF_METHANE%lake_liqdiff_scale, &
         DEF_METHANE%lake_o2_liqdiff_scale, DEF_METHANE%grnd_methane_cond_default, &
         DEF_METHANE%mino2lim, DEF_METHANE%q10methane_base, DEF_METHANE%q10lakebase, &
         DEF_METHANE%cnscalefactor, DEF_METHANE%redoxlag, DEF_METHANE%lake_decomp_fact, &
         DEF_METHANE%redoxlag_vertical, DEF_METHANE%pHmax, DEF_METHANE%pHmin, &
         DEF_METHANE%ph_factor_floor, &
         DEF_METHANE%oxinhib, DEF_METHANE%smp_crit, DEF_METHANE%bubble_f, &
         DEF_METHANE%aereoxid, DEF_METHANE%tiller_C, DEF_METHANE%satpow, &
         DEF_METHANE%capthick, DEF_METHANE%atm_methane, DEF_METHANE%om_frac_sf, &
         DEF_METHANE%wtd_inflection, DEF_METHANE%wtd_steepness, &
         DEF_METHANE%wtd_inflection_soil, DEF_METHANE%wtd_steepness_soil, &
         DEF_METHANE%hybrid_soil_threshold, DEF_METHANE%rice_drain_window_days, &
         DEF_METHANE%rice_paddy_min_finundated, DEF_METHANE%rice_midseason_start_days, &
         DEF_METHANE%wetland_max_wtd, DEF_METHANE%wetland_max_wtd_tropical_peat, &
         DEF_METHANE%wetland_max_wtd_tropical_floodplain, &
         DEF_METHANE%wetland_max_wtd_temperate_marsh, &
         DEF_METHANE%wetland_max_wtd_boreal_fen, DEF_METHANE%wetland_max_wtd_boreal_bog, &
         DEF_METHANE%wetland_max_wtd_bog_share, DEF_METHANE%wetland_bog_share_site, &
         DEF_METHANE%wetland_share_emerg_site, DEF_METHANE%wetland_share_dome_site, &
         DEF_METHANE%wetland_share_moss_site, DEF_METHANE%wetland_share_rain_site, &
         DEF_METHANE%wetland_share_bog_site, &
         DEF_METHANE%rice_midseason_drain_days, DEF_METHANE%rice_midseason_drained_finundated, &
         DEF_METHANE%rice_substrate_boost, DEF_METHANE%numerical_correction_fatal_threshold, &
         DEF_METHANE%host_water_tolerance, &
         DEF_METHANE%peat_K0, DEF_METHANE%peat_m, DEF_METHANE%peat_c, DEF_METHANE%organic_max, &
         DEF_METHANE%wetland_microtopo_sigma, DEF_METHANE%wetland_npp_frac, DEF_METHANE%wetland_litter_cn, &
         DEF_METHANE%wetland_bg_frac, DEF_METHANE%wetland_bg_frac_forest, &
         DEF_METHANE%wetland_lai_open_peat, DEF_METHANE%wetland_lai_marsh, &
         DEF_METHANE%wetland_forest_share_site, DEF_METHANE%wetland_lai_cap_site, &
         DEF_METHANE%wetland_cover_frac_site, DEF_METHANE%wetland_burial_frac_site, &
         DEF_METHANE%wetland_anoxia_catotelm, DEF_METHANE%wetland_anoxia_efold, &
         DEF_METHANE%wetland_tau_s3, DEF_METHANE%wetland_burial_velocity, &
         DEF_METHANE%wetland_exudate_frac, DEF_METHANE%root_exudate_frac, &
         DEF_METHANE%rice_aereoxid, DEF_METHANE%wetland_root_efold, &
         DEF_METHANE%wetland_f_tropical_peat, DEF_METHANE%wetland_f_tree_ratio, &
         DEF_METHANE%woody_conduit_area, DEF_METHANE%wetland_forest_htop_trop, &
         DEF_METHANE%tundra_high_share_site, &
         DEF_METHANE%methanogen_k2, DEF_METHANE%methanogen_q2, DEF_METHANE%methanogen_cue, &
         DEF_METHANE%methanogen_alpha, DEF_METHANE%methanogen_mu, DEF_METHANE%methanogen_qev_ratio, &
         DEF_METHANE%acceptor_cap_min, DEF_METHANE%acceptor_cap_org, DEF_METHANE%acceptor_k_half, &
         DEF_METHANE%acceptor_eta0, DEF_METHANE%acceptor_k_reox, DEF_METHANE%acceptor_k_o2, &
         DEF_METHANE_hydrology%vdcf, &
         DEF_METHANE_hydrology%slopebeta, DEF_METHANE_hydrology%slopemax, &
         DEF_METHANE_hydrology%pc]))) THEN
         IF (p_is_master) write(6,*) '***** ERROR: methane namelist contains a NaN or infinite real parameter.'
         bad = .true.
      ENDIF

      IF (.not. DEF_METHANE%methane_offline) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: methane_offline=.false. requires an online atmosphere/NEE coupling that is not implemented.'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%methane_frzout) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: methane_frzout is retired; ice is always excluded from mobile CH4 storage.'
         bad = .true.
      ENDIF

      IF (DEF_METHANE%f_methane < 0._r8 .or. DEF_METHANE%f_methane > 0.5_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: f_methane out of [0,0.5]: ', DEF_METHANE%f_methane
         bad = .true.
      ENDIF
      IF (DEF_METHANE%f_methane_tropical_peat < 0._r8 .or. &
          DEF_METHANE%f_methane_tropical_peat > 0.5_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: f_methane_tropical_peat out of [0,0.5]: ', &
            DEF_METHANE%f_methane_tropical_peat
         bad = .true.
      ENDIF
      IF (DEF_METHANE%f_methane_tropical_floodplain < 0._r8 .or. &
          DEF_METHANE%f_methane_tropical_floodplain > 0.5_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: f_methane_tropical_floodplain out of [0,0.5]: ', &
            DEF_METHANE%f_methane_tropical_floodplain
         bad = .true.
      ENDIF
      IF (DEF_METHANE%f_methane_floodplain < 0._r8 .or. &
          DEF_METHANE%f_methane_floodplain > 0.5_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: f_methane_floodplain out of [0,0.5]: ', &
            DEF_METHANE%f_methane_floodplain
         bad = .true.
      ENDIF
      IF (DEF_METHANE%f_methane_temperate_marsh < 0._r8 .or. &
          DEF_METHANE%f_methane_temperate_marsh > 0.5_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: f_methane_temperate_marsh out of [0,0.5]: ', &
            DEF_METHANE%f_methane_temperate_marsh
         bad = .true.
      ENDIF
      IF (DEF_METHANE%f_methane_boreal_fen < 0._r8 .or. &
          DEF_METHANE%f_methane_boreal_fen > 0.5_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: f_methane_boreal_fen out of [0,0.5]: ', &
            DEF_METHANE%f_methane_boreal_fen
         bad = .true.
      ENDIF
      IF (DEF_METHANE%f_methane_boreal_bog < 0._r8 .or. &
          DEF_METHANE%f_methane_boreal_bog > 0.5_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: f_methane_boreal_bog out of [0,0.5]: ', &
            DEF_METHANE%f_methane_boreal_bog
         bad = .true.
      ENDIF
      IF (DEF_METHANE%f_methane_rice_paddy < 0._r8 .or. &
          DEF_METHANE%f_methane_rice_paddy > 0.5_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: f_methane_rice_paddy out of [0,0.5]: ', &
            DEF_METHANE%f_methane_rice_paddy
         bad = .true.
      ENDIF
      IF (DEF_METHANE%f_methane_upland_soil < 0._r8 .or. &
          DEF_METHANE%f_methane_upland_soil > 0.5_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: f_methane_upland_soil out of [0,0.5]: ', &
            DEF_METHANE%f_methane_upland_soil
         bad = .true.
      ENDIF
      IF (DEF_METHANE%q10methane <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: q10methane must be > 0: ', DEF_METHANE%q10methane
         bad = .true.
      ENDIF
      IF (DEF_METHANE%q10_methane_oxid <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: q10_methane_oxid must be > 0: ', DEF_METHANE%q10_methane_oxid
         bad = .true.
      ENDIF
      IF (DEF_METHANE%vmax_methane_oxid < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: vmax_methane_oxid must be >= 0: ', DEF_METHANE%vmax_methane_oxid
         bad = .true.
      ENDIF
      IF (DEF_METHANE%vmax_oxid_unsat < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: vmax_oxid_unsat must be >= 0: ', DEF_METHANE%vmax_oxid_unsat
         bad = .true.
      ENDIF
      IF (DEF_METHANE%k_m       <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: k_m must be > 0: ', DEF_METHANE%k_m
         bad = .true.
      ENDIF
      IF (DEF_METHANE%k_m_unsat <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: k_m_unsat must be > 0: ', DEF_METHANE%k_m_unsat
         bad = .true.
      ENDIF
      IF (DEF_METHANE%k_m_o2    <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: k_m_o2 must be > 0: ', DEF_METHANE%k_m_o2
         bad = .true.
      ENDIF
      IF (DEF_METHANE%lake_oxid_scale < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: lake_oxid_scale must be >= 0: ', &
            DEF_METHANE%lake_oxid_scale
         bad = .true.
      ENDIF
      IF (DEF_METHANE%lake_k_m_o2 /= -1._r8 .and. DEF_METHANE%lake_k_m_o2 <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: lake_k_m_o2 must be -1 or > 0: ', &
            DEF_METHANE%lake_k_m_o2
         bad = .true.
      ENDIF
      IF (DEF_METHANE%lake_vmax_methane_oxid /= -1._r8 .and. &
          DEF_METHANE%lake_vmax_methane_oxid < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: lake_vmax_methane_oxid must be -1 or >= 0: ', &
            DEF_METHANE%lake_vmax_methane_oxid
         bad = .true.
      ENDIF
      IF (DEF_METHANE%lake_oxic_sediment_depth /= -1._r8 .and. &
          DEF_METHANE%lake_oxic_sediment_depth <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: lake_oxic_sediment_depth must be -1 or > 0: ', &
            DEF_METHANE%lake_oxic_sediment_depth
         bad = .true.
      ENDIF
      ! above 100 d-1 the store empties within about an hour of ice-out,
      ! two orders of magnitude beyond the published 0.14-1 d-1
      IF (DEF_METHANE%lake_icebubble_release > 100._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: lake_icebubble_release must be <= 100 d-1 (<= 0 off): ', &
            DEF_METHANE%lake_icebubble_release
         bad = .true.
      ENDIF
      IF (DEF_METHANE%lake_icebubble_dissol < 0._r8 .or. DEF_METHANE%lake_icebubble_dissol > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: lake_icebubble_dissol out of [0,1]: ', &
            DEF_METHANE%lake_icebubble_dissol
         bad = .true.
      ENDIF
      ! the supply stays within five times the largest measured
      ! sediment mineralization (196 g C m-2 yr-1, Sobek et al. 2009); k20 within
      ! about seven times FLaMe's upper bound; theta from none to a Q10 of 6
      IF (any(DEF_METHANE%lake_cdep_band < 0._r8) .or. any(DEF_METHANE%lake_cdep_band > 1000._r8)) THEN
         IF (p_is_master) write(6,*) '***** ERROR: lake_cdep_band out of [0,1000] g C m-2 yr-1: ', &
            DEF_METHANE%lake_cdep_band
         bad = .true.
      ENDIF
      IF (DEF_METHANE%lake_k20 <= 0._r8 .or. DEF_METHANE%lake_k20 > 0.1_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: lake_k20 out of (0,0.1] d-1: ', DEF_METHANE%lake_k20
         bad = .true.
      ENDIF
      IF (DEF_METHANE%lake_theta < 1._r8 .or. DEF_METHANE%lake_theta > 1.2_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: lake_theta out of [1,1.2]: ', DEF_METHANE%lake_theta
         bad = .true.
      ENDIF
      IF (DEF_METHANE%lake_active_pool .and. .not. DEF_METHANE%allowlakeprod) THEN
         IF (p_is_master) write(6,*) '***** ERROR: lake_active_pool needs allowlakeprod'
         bad = .true.
      ENDIF
      ! same stoichiometric bound as f_methane; set without the pool it
      ! would be read and never used
      IF (DEF_METHANE%lake_f_methane > 0.5_r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: lake_f_methane must be negative (off) ', &
            'or lie in [0, 0.5]: ', DEF_METHANE%lake_f_methane
         bad = .true.
      ENDIF
      IF (DEF_METHANE%lake_f_methane >= 0._r8 .and. .not. DEF_METHANE%lake_active_pool) THEN
         IF (p_is_master) write(6,*) '***** ERROR: lake_f_methane >= 0 needs lake_active_pool'
         bad = .true.
      ENDIF
      ! a nitrate set with a non-positive half-inhibition constant
      ! would divide by zero or turn the factor negative
      IF (DEF_METHANE%wetland_nitrate_site > 0._r8 .and. DEF_METHANE%wetland_nitrate_ki <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_nitrate_site needs wetland_nitrate_ki > 0: ', &
            DEF_METHANE%wetland_nitrate_ki
         bad = .true.
      ENDIF
      IF (DEF_METHANE%aereoxid < 0._r8 .or. DEF_METHANE%aereoxid > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: aereoxid out of [0,1]: ', DEF_METHANE%aereoxid
         bad = .true.
      ENDIF
      IF (DEF_METHANE%bubble_f <= 0._r8 .or. DEF_METHANE%bubble_f > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: bubble_f out of (0,1]: ', DEF_METHANE%bubble_f
         bad = .true.
      ENDIF
      IF (DEF_METHANE%pHmin >= DEF_METHANE%pHmax) THEN
         IF (p_is_master) write(6,*) '***** ERROR: pHmin >= pHmax: ', DEF_METHANE%pHmin, DEF_METHANE%pHmax
         bad = .true.
      ENDIF
      IF (DEF_METHANE%ph_factor_floor < 0._r8 .or. DEF_METHANE%ph_factor_floor > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: ph_factor_floor out of [0,1]: ', DEF_METHANE%ph_factor_floor
         bad = .true.
      ENDIF
      IF (DEF_METHANE%mino2lim < 0._r8 .or. DEF_METHANE%mino2lim > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: mino2lim out of [0,1]: ', DEF_METHANE%mino2lim
         bad = .true.
      ENDIF
      IF (DEF_METHANE%atm_methane < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: atm_methane must be >= 0: ', DEF_METHANE%atm_methane
         bad = .true.
      ENDIF
      ! A fraction of pore volume: 0 restores the machine-epsilon behaviour that
      ! aborted on its first timestep, and >=1 would accept any disagreement at
      ! all, which defeats the guard.
      IF (DEF_METHANE%host_water_tolerance < 0._r8 .or. &
          DEF_METHANE%host_water_tolerance >= 1._r8) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: host_water_tolerance out of [0,1): ', DEF_METHANE%host_water_tolerance
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wtd_steepness <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wtd_steepness must be > 0: ', DEF_METHANE%wtd_steepness
         bad = .true.
      ENDIF
      IF (.not. ieee_is_finite(DEF_METHANE%wtd_inflection_soil) .or. &
          DEF_METHANE%wtd_inflection_soil < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wtd_inflection_soil must be finite and >= 0 m: ', &
            DEF_METHANE%wtd_inflection_soil
         bad = .true.
      ENDIF
      IF (.not. ieee_is_finite(DEF_METHANE%wtd_steepness_soil) .or. &
          DEF_METHANE%wtd_steepness_soil <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wtd_steepness_soil must be finite and > 0 m: ', &
            DEF_METHANE%wtd_steepness_soil
         bad = .true.
      ENDIF
      IF (.not. ieee_is_finite(DEF_METHANE%hybrid_soil_threshold) .or. &
          DEF_METHANE%hybrid_soil_threshold < 0._r8 .or. &
          DEF_METHANE%hybrid_soil_threshold > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: hybrid_soil_threshold must be finite and in [0,1]: ', &
            DEF_METHANE%hybrid_soil_threshold
         bad = .true.
      ENDIF
      IF (.not. ieee_is_finite(DEF_METHANE%z0_methane_prod) .or. &
          DEF_METHANE%z0_methane_prod < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: z0_methane_prod must be finite and >= 0 m: ', &
            DEF_METHANE%z0_methane_prod
         bad = .true.
      ENDIF
      IF (DEF_METHANE%vgc_max <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: vgc_max must be > 0: ', DEF_METHANE%vgc_max
         bad = .true.
      ENDIF
      IF (DEF_METHANE%poros_tiller < 0._r8 .or. DEF_METHANE%poros_tiller > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: poros_tiller out of [0,1]: ', DEF_METHANE%poros_tiller
         bad = .true.
      ENDIF
      IF (DEF_METHANE%nongrassporosratio < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: nongrassporosratio must be >= 0: ', DEF_METHANE%nongrassporosratio
         bad = .true.
      ENDIF
      IF (DEF_METHANE%unsat_aere_ratio < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: unsat_aere_ratio must be >= 0: ', DEF_METHANE%unsat_aere_ratio
         bad = .true.
      ENDIF
      IF (DEF_METHANE%porosmin < 0._r8 .or. DEF_METHANE%porosmin > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: porosmin out of [0,1]: ', DEF_METHANE%porosmin
         bad = .true.
      ENDIF
      IF (DEF_METHANE%aere_radius <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: aere_radius must be > 0: ', DEF_METHANE%aere_radius
         bad = .true.
      ENDIF
      IF (DEF_METHANE%rob <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: rob must be > 0: ', DEF_METHANE%rob
         bad = .true.
      ENDIF
      IF (DEF_METHANE%scale_factor_aere < 0._r8 .or. &
          DEF_METHANE%scale_factor_gasdiff < 0._r8 .or. &
          DEF_METHANE%scale_factor_liqdiff < 0._r8 .or. &
          DEF_METHANE%lake_liqdiff_scale < 0._r8 .or. &
          DEF_METHANE%lake_o2_liqdiff_scale < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: methane transport/aerenchyma scale factors must be >= 0: ', &
            DEF_METHANE%scale_factor_aere, DEF_METHANE%scale_factor_gasdiff, &
            DEF_METHANE%scale_factor_liqdiff, DEF_METHANE%lake_liqdiff_scale, &
            DEF_METHANE%lake_o2_liqdiff_scale
         bad = .true.
      ENDIF
      IF (DEF_METHANE%grnd_methane_cond_default <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: grnd_methane_cond_default must be > 0: ', &
            DEF_METHANE%grnd_methane_cond_default
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wetland_aereoxid > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_aereoxid must not exceed 1: ', DEF_METHANE%wetland_aereoxid
         bad = .true.
      ENDIF
      IF (DEF_METHANE%woody_conduit_area > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: woody_conduit_area must be negative (off) ', &
            'or lie in [0, 1] m2 m-2: ', DEF_METHANE%woody_conduit_area
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wetland_forest_htop_trop >= 0._r8 .and. &
          (DEF_METHANE%wetland_forest_htop_trop < 1._r8 .or. DEF_METHANE%wetland_forest_htop_trop > 100._r8 .or. &
           .not. DEF_METHANE%wetland_veg_glwd)) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_forest_htop_trop must be negative (off) ', &
            'or lie in [1, 100] m, and needs wetland_veg_glwd: ', DEF_METHANE%wetland_forest_htop_trop
         bad = .true.
      ENDIF
      IF (DEF_METHANE%rox_on_rooted_production .and. DEF_METHANE%wetland_aereoxid < 0._r8 .and. &
          DEF_METHANE%rice_aereoxid <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: rox_on_rooted_production needs a rhizosphere ', &
            'oxidation share: wetland_aereoxid >= 0 or rice_aereoxid > 0'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%q10methane_base_weight < 0._r8 .or. DEF_METHANE%q10methane_base_weight > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: q10methane_base_weight must lie in [0,1]: ', &
            DEF_METHANE%q10methane_base_weight
         bad = .true.
      ENDIF
      IF (DEF_METHANE%q10methane_base <= 0._r8 .or. DEF_METHANE%q10lakebase <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: q10 base temperatures must be > 0 K: ', &
            DEF_METHANE%q10methane_base, DEF_METHANE%q10lakebase
         bad = .true.
      ENDIF
      IF (DEF_METHANE%methanogen_activity) THEN
         IF (DEF_METHANE%methanogen_k2 <= 0._r8 .or. DEF_METHANE%methanogen_q2 <= 0._r8 .or. &
             DEF_METHANE%methanogen_mu <= 0._r8 .or. DEF_METHANE%methanogen_alpha <= 0._r8 .or. &
             DEF_METHANE%methanogen_cue > 1._r8 .or. &
             DEF_METHANE%methanogen_alpha >= DEF_METHANE%methanogen_cue) THEN
            IF (p_is_master) write(6,*) '***** ERROR: methanogen_activity needs methanogen_k2, ', &
               'methanogen_q2, methanogen_mu > 0 and 0 < methanogen_alpha < methanogen_cue <= 1: ', &
               DEF_METHANE%methanogen_k2, DEF_METHANE%methanogen_q2, DEF_METHANE%methanogen_mu, &
               DEF_METHANE%methanogen_alpha, DEF_METHANE%methanogen_cue
            bad = .true.
         ENDIF
         IF (DEF_METHANE%q10methane_local_base .or. DEF_METHANE%q10methane_base_unfrozen) THEN
            IF (p_is_master) write(6,*) '***** ERROR: methanogen_activity (candidate 30) is an ', &
               'alternative to q10methane_local_base / q10methane_base_unfrozen (C-27, C-51); ', &
               'switch one of them off'
            bad = .true.
         ENDIF
      ENDIF
      IF (DEF_METHANE%acceptor_pool) THEN
         IF (DEF_METHANE%acceptor_cap_min < 0._r8 .or. DEF_METHANE%acceptor_cap_org < 0._r8 .or. &
             DEF_METHANE%acceptor_k_half <= 0._r8 .or. DEF_METHANE%acceptor_eta0 < 0._r8 .or. &
             DEF_METHANE%acceptor_k_reox < 0._r8 .or. DEF_METHANE%acceptor_k_o2 <= 0._r8) THEN
            IF (p_is_master) write(6,*) '***** ERROR: acceptor_pool needs acceptor_cap_min, ', &
               'acceptor_cap_org, acceptor_eta0, acceptor_k_reox >= 0 and acceptor_k_half, ', &
               'acceptor_k_o2 > 0: ', DEF_METHANE%acceptor_cap_min, DEF_METHANE%acceptor_cap_org, &
               DEF_METHANE%acceptor_eta0, DEF_METHANE%acceptor_k_reox, DEF_METHANE%acceptor_k_half, &
               DEF_METHANE%acceptor_k_o2
            bad = .true.
         ENDIF
         ! The pool is the dynamic form of the redox lags; both at once would
         ! count the same acceptors twice.
         IF (DEF_METHANE%use_biome_redoxlag .or. DEF_METHANE%redoxlag > 0._r8 .or. &
             (DEF_METHANE%use_vertical_redoxlag .and. DEF_METHANE%redoxlag_vertical > 0._r8)) THEN
            IF (p_is_master) write(6,*) '***** ERROR: acceptor_pool (candidate 16) replaces the ', &
               'redox lags; set use_biome_redoxlag = .false., redoxlag = 0 and ', &
               'use_vertical_redoxlag = .false. (or redoxlag_vertical = 0): ', &
               DEF_METHANE%use_biome_redoxlag, DEF_METHANE%redoxlag, &
               DEF_METHANE%use_vertical_redoxlag, DEF_METHANE%redoxlag_vertical
            bad = .true.
         ENDIF
      ENDIF
      IF (DEF_METHANE%cnscalefactor < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: cnscalefactor must be >= 0: ', DEF_METHANE%cnscalefactor
         bad = .true.
      ENDIF
      IF (DEF_METHANE%redoxlag < 0._r8 .or. DEF_METHANE%redoxlag_vertical < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: redox lags must be >= 0 days: ', &
            DEF_METHANE%redoxlag, DEF_METHANE%redoxlag_vertical
         bad = .true.
      ENDIF
      IF (min(DEF_METHANE%redoxlag_tropical_peat, &
              DEF_METHANE%redoxlag_tropical_floodplain, &
              DEF_METHANE%redoxlag_temperate_marsh, DEF_METHANE%redoxlag_boreal_fen, &
              DEF_METHANE%redoxlag_boreal_bog, DEF_METHANE%redoxlag_rice_paddy, &
              DEF_METHANE%redoxlag_upland_soil, DEF_METHANE%redoxlag_wetland_dim) < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: biome redox lags must be >= 0 days.'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%oxinhib < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: oxinhib must be >= 0 m3 mol-1: ', DEF_METHANE%oxinhib
         bad = .true.
      ENDIF
      IF (DEF_METHANE%B_max_fraction_methanogen > 1._r8 .or. &
          DEF_METHANE%B_max_fraction_methanotroph > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: enabled microbial biomass fractions must be <= 1: ', &
            DEF_METHANE%B_max_fraction_methanogen, DEF_METHANE%B_max_fraction_methanotroph
         bad = .true.
      ENDIF
      IF (DEF_METHANE%lake_decomp_fact < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: lake_decomp_fact must be >= 0: ', DEF_METHANE%lake_decomp_fact
         bad = .true.
      ENDIF
      IF (DEF_METHANE%smp_crit >= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: smp_crit must be negative [mm]: ', DEF_METHANE%smp_crit
         bad = .true.
      ENDIF
      IF (DEF_METHANE%satpow <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: satpow must be > 0: ', DEF_METHANE%satpow
         bad = .true.
      ENDIF
      IF (DEF_METHANE%capthick < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: capthick must be >= 0: ', DEF_METHANE%capthick
         bad = .true.
      ENDIF
      ! Rice physiology and reserved substrate-control range checks.
      IF (DEF_METHANE%rice_drain_window_days <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: rice_drain_window_days must be > 0: ', &
            DEF_METHANE%rice_drain_window_days
         bad = .true.
      ENDIF
      IF (DEF_METHANE%rice_paddy_min_finundated < 0._r8 .or. &
          DEF_METHANE%rice_paddy_min_finundated > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: rice_paddy_min_finundated must be in [0,1]: ', &
            DEF_METHANE%rice_paddy_min_finundated
         bad = .true.
      ENDIF
      IF (DEF_METHANE%rice_midseason_drained_finundated < 0._r8 .or. &
          DEF_METHANE%rice_midseason_drained_finundated > 1._r8) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: rice_midseason_drained_finundated must be in [0,1]: ', &
            DEF_METHANE%rice_midseason_drained_finundated
         bad = .true.
      ENDIF
      IF (DEF_METHANE%rice_midseason_start_days < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: rice_midseason_start_days must be >= 0: ', &
            DEF_METHANE%rice_midseason_start_days
         bad = .true.
      ENDIF
      IF (DEF_METHANE%rice_midseason_drain_days < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: rice_midseason_drain_days must be >= 0: ', &
            DEF_METHANE%rice_midseason_drain_days
         bad = .true.
      ENDIF
      IF (abs(DEF_METHANE%rice_substrate_boost - 1._r8) > 10._r8 * epsilon(1._r8)) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: rice_substrate_boost must remain 1 until methane production debits BGC carbon: ', &
            DEF_METHANE%rice_substrate_boost
         bad = .true.
      ENDIF
      IF (len_trim(DEF_METHANE%ch4_history_vars) == 0) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: ch4_history_vars must be core/diagnostic/all/none or a comma-separated list'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%tiller_C <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: tiller_C must be > 0: ', DEF_METHANE%tiller_C
         bad = .true.
      ENDIF
      IF (DEF_METHANE%om_frac_sf < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: om_frac_sf must be >= 0: ', DEF_METHANE%om_frac_sf
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wetland_veg_glwd .and. (trim(DEF_METHANE%wetland_veg_file) == 'null' .or. &
         DEF_METHANE%wetland_bg_frac_forest < 0._r8 .or. DEF_METHANE%wetland_bg_frac_forest > 1._r8 .or. &
         DEF_METHANE%wetland_lai_open_peat <= 0._r8 .or. DEF_METHANE%wetland_lai_marsh <= 0._r8)) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_veg_glwd needs wetland_veg_file, ', &
            '0 <= wetland_bg_frac_forest <= 1 and positive wetland_lai_open_peat, wetland_lai_marsh'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wetland_max_wtd_bog_share > 1._r8 .or. DEF_METHANE%wetland_bog_share_site > 1._r8 .or. &
         (DEF_METHANE%wetland_max_wtd_bog_share >= 0._r8 .and. .not. DEF_METHANE%wetland_veg_glwd)) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_max_wtd_bog_share must be negative (off) or ', &
            'lie in [0, 1] m and needs wetland_veg_glwd (bog_share of wetland_veg_file, or ', &
            'wetland_bog_share_site), and wetland_bog_share_site must not exceed 1: ', &
            DEF_METHANE%wetland_max_wtd_bog_share, DEF_METHANE%wetland_bog_share_site
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wetland_share_emerg_site > 1._r8 .or. DEF_METHANE%wetland_share_dome_site > 1._r8 .or. &
          DEF_METHANE%wetland_share_moss_site > 1._r8 .or. &
          (DEF_METHANE%wetland_share_emerg_site >= 0._r8 .and. DEF_METHANE%wetland_share_dome_site >= 0._r8 .and. &
           DEF_METHANE%wetland_share_emerg_site + DEF_METHANE%wetland_share_dome_site > 1._r8)) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_share_emerg_site, wetland_share_dome_site and ', &
            'wetland_share_moss_site must be negative (file) or lie in [0, 1], and the emergent and ', &
            'dome shares must not sum above 1: ', DEF_METHANE%wetland_share_emerg_site, &
            DEF_METHANE%wetland_share_dome_site, DEF_METHANE%wetland_share_moss_site
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wetland_share_rain_site > 1._r8 .or. DEF_METHANE%wetland_share_bog_site > 1._r8 .or. &
          (DEF_METHANE%wetland_share_rain_site >= 0._r8 .and. &
           DEF_METHANE%wetland_share_bog_site > DEF_METHANE%wetland_share_rain_site)) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_share_rain_site and wetland_share_bog_site ', &
            'must be negative (file) or lie in [0, 1], and the open-bog share must not exceed the ', &
            'rain-fed share: ', DEF_METHANE%wetland_share_rain_site, DEF_METHANE%wetland_share_bog_site
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wetland_dim_class .and. (.not. DEF_METHANE%wetland_veg_glwd .or. &
          (DEF_METHANE%use_biome_wetland_max_wtd .and. min(DEF_METHANE%wetland_max_wtd, &
           DEF_METHANE%wetland_max_wtd_temperate_marsh, DEF_METHANE%wetland_max_wtd_tropical_peat) < 0._r8))) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_dim_class needs wetland_veg_glwd (the GLWD ', &
            'classes of wetland_veg_file) and, with use_biome_wetland_max_wtd, non-negative ', &
            'wetland_max_wtd, wetland_max_wtd_temperate_marsh and wetland_max_wtd_tropical_peat'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wetland_cover_frac_site > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_cover_frac_site must be negative (off) ', &
            'or lie in [0, 1]: ', DEF_METHANE%wetland_cover_frac_site
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wetland_burial_frac_site > 1._r8 .or. &
         (DEF_METHANE%wetland_burial_frac_site > 0._r8 .and. .not. DEF_METHANE%wetland_plant_input)) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_burial_frac_site must be negative (off) ', &
            'or lie in [0, 1], and needs wetland_plant_input: ', DEF_METHANE%wetland_burial_frac_site
         bad = .true.
      ENDIF
      IF (DEF_METHANE%tundra_high_share_site > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: tundra_high_share_site must not exceed 1: ', &
            DEF_METHANE%tundra_high_share_site
         bad = .true.
      ENDIF
#ifndef SinglePoint
      IF (DEF_METHANE%tundra_high_share_site > 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: tundra_high_share_site is a single-point key ', &
            '(no global map of high microsites): ', DEF_METHANE%tundra_high_share_site
         bad = .true.
      ENDIF
#endif
      IF (DEF_METHANE%wetland_forest_input_herb .and. &
         .not. (DEF_METHANE%wetland_plant_input .and. DEF_METHANE%wetland_veg_glwd)) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_forest_input_herb needs wetland_plant_input ', &
            'and wetland_veg_glwd (the forested share and LAI cap of the tile)'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wetland_moss_input_frac < 0._r8 .or. DEF_METHANE%wetland_moss_input_frac > 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_moss_input_frac must lie in [0, 1]: ', &
            DEF_METHANE%wetland_moss_input_frac
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wetland_tau_s3 == 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_tau_s3 must be positive, or negative to keep tau_s3'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wetland_bgc_sasu .and. DEF_METHANE%wetland_fixed_substrate) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_bgc_sasu needs wetland_fixed_substrate off: ', &
            'held pools have no steady state to jump to'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%rice_aereoxid < 0._r8 .or. DEF_METHANE%rice_aereoxid >= 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: rice_aereoxid must lie in [0, 1): ', DEF_METHANE%rice_aereoxid
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wetland_exudate_frac /= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_exudate_frac (C-23) is retired, it ', &
            'debits each exudate from the metabolic litter twice; set it to 0 and use ', &
            'root_exudate_frac: ', DEF_METHANE%wetland_exudate_frac
         bad = .true.
      ENDIF
      IF (DEF_METHANE%root_exudate_frac >= 1._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: root_exudate_frac must lie below 1: ', &
            DEF_METHANE%root_exudate_frac
         bad = .true.
      ELSEIF (DEF_METHANE%root_exudate_frac > 0._r8 .and. p_is_master) THEN
         write(6,'(A,F6.3)') ' root_exudate_frac (candidate 6): share of NPP exuded ', &
            DEF_METHANE%root_exudate_frac
         IF (DEF_METHANE%wetland_exudate_frac > 0._r8) write(6,'(A)') &
            '   wetland tiles keep wetland_exudate_frac (C-23); candidate 6 acts on soil tiles'
      ENDIF
      IF (DEF_METHANE%wetland_bgc_sasu .and. DEF_METHANE%wetland_vert_transp) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_bgc_sasu solves each layer on its own ', &
            'and cannot be combined with wetland_vert_transp'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wetland_anoxia_catotelm >= 0._r8 .and. &
         (DEF_METHANE%wetland_anoxia_catotelm > DEF_METHANE%mino2lim .or. &
          DEF_METHANE%wetland_anoxia_efold <= 0._r8)) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_anoxia_catotelm must not exceed mino2lim ', &
            'and wetland_anoxia_efold must be positive: ', &
            DEF_METHANE%wetland_anoxia_catotelm, DEF_METHANE%mino2lim, DEF_METHANE%wetland_anoxia_efold
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wetland_plant_input .and. (DEF_METHANE%wetland_fixed_substrate .or. &
         DEF_METHANE%wetland_npp_frac <= 0._r8 .or. DEF_METHANE%wetland_npp_frac > 1._r8 .or. &
         DEF_METHANE%wetland_litter_cn <= 0._r8 .or. &
         DEF_METHANE%wetland_bg_frac < 0._r8 .or. DEF_METHANE%wetland_bg_frac > 1._r8)) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_plant_input needs wetland_fixed_substrate off, ', &
            '0 < wetland_npp_frac <= 1, wetland_litter_cn > 0 and 0 <= wetland_bg_frac <= 1: ', &
            DEF_METHANE%wetland_fixed_substrate, DEF_METHANE%wetland_npp_frac, DEF_METHANE%wetland_litter_cn, &
            DEF_METHANE%wetland_bg_frac
         bad = .true.
      ENDIF
      IF ((DEF_METHANE%wetland_n_unlimited .or. DEF_METHANE%wetland_n_uptake .or. &
           DEF_METHANE%wetland_vert_transp) .and. DEF_METHANE%wetland_fixed_substrate) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_n_unlimited, wetland_n_uptake and ', &
            'wetland_vert_transp need wetland_fixed_substrate off'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wetland_n_uptake .and. .not. DEF_METHANE%wetland_plant_input) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_n_uptake needs wetland_plant_input'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wetland_burial_velocity < 0._r8 .or. DEF_METHANE%wetland_burial_velocity > 0.01_r8 .or. &
         (DEF_METHANE%wetland_burial_velocity > 0._r8 .and. .not. DEF_METHANE%wetland_vert_transp)) THEN
         IF (p_is_master) write(6,*) '***** ERROR: wetland_burial_velocity must lie in [0, 0.01] m yr-1 ', &
            'and needs wetland_vert_transp: ', DEF_METHANE%wetland_burial_velocity
         bad = .true.
      ENDIF
      IF (DEF_METHANE%wetland_peat_drainage .and. &
         (DEF_METHANE%peat_K0 <= 0._r8 .or. DEF_METHANE%peat_m <= 1._r8 .or. &
          DEF_METHANE%peat_c <= 0._r8 .or. DEF_METHANE%organic_max <= 0._r8)) THEN
         IF (p_is_master) write(6,*) '***** ERROR: peat drainage needs peat_K0, peat_c, organic_max > 0 and peat_m > 1: ', &
            DEF_METHANE%peat_K0, DEF_METHANE%peat_m, DEF_METHANE%peat_c, DEF_METHANE%organic_max
         bad = .true.
      ENDIF
      IF (DEF_METHANE%K_substrate_methanogen_pool <= 0._r8 .or. &
          DEF_METHANE%K_inh_O2_methanogen <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: microbial half-saturation/inhibition constants must be > 0: ', &
            DEF_METHANE%K_substrate_methanogen_pool, DEF_METHANE%K_inh_O2_methanogen
         bad = .true.
      ENDIF
      IF (DEF_METHANE%B_init_methanogen < 0._r8 .or. DEF_METHANE%B_init_methanotroph < 0._r8 .or. &
          DEF_METHANE%B_min_methanogen < 0._r8 .or. DEF_METHANE%B_min_methanotroph < 0._r8 .or. &
          DEF_METHANE%mu_max_methanogen < 0._r8 .or. DEF_METHANE%mu_max_methanotroph < 0._r8 .or. &
          DEF_METHANE%gamma_methanogen < 0._r8 .or. DEF_METHANE%gamma_methanotroph < 0._r8 .or. &
          DEF_METHANE%gamma_microbial_dormant < 0._r8 .or. DEF_METHANE%gamma_microbial_freeze < 0._r8 .or. &
          DEF_METHANE%kappa_m_methanogen < 0._r8 .or. DEF_METHANE%kappa_m_methanotroph < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: microbial biomass/rate/loss parameters must be >= 0'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%use_microbial_flux_override .and. &
          .not. DEF_METHANE%use_microbial_pools) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: use_microbial_flux_override requires use_microbial_pools.'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%use_microbial_dormancy .and. &
          .not. DEF_METHANE%use_microbial_pools) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: use_microbial_dormancy requires use_microbial_pools.'
         bad = .true.
      ENDIF
      IF (DEF_METHANE%use_microbial_pools) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: microbial pools are disabled until biomass growth/loss has donor/sink carbon coupling.'
         bad = .true.
         IF (.not. ieee_is_finite(DEF_METHANE%B_max_fraction_methanogen) .or. &
             .not. ieee_is_finite(DEF_METHANE%B_max_fraction_methanotroph) .or. &
             DEF_METHANE%B_max_fraction_methanogen <= 0._r8 .or. &
             DEF_METHANE%B_max_fraction_methanotroph <= 0._r8 .or. &
             DEF_METHANE%B_max_fraction_methanogen > 1._r8 .or. &
             DEF_METHANE%B_max_fraction_methanotroph > 1._r8) THEN
            IF (p_is_master) write(6,*) &
               '***** ERROR: microbial pools require finite B_max fractions in (0,1]: ', &
               DEF_METHANE%B_max_fraction_methanogen, DEF_METHANE%B_max_fraction_methanotroph
            bad = .true.
         ENDIF
      ENDIF
      IF (DEF_METHANE%max_microbe_prod_multiplier <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: max_microbe_prod_multiplier must be > 0: ', &
            DEF_METHANE%max_microbe_prod_multiplier
         bad = .true.
      ENDIF
      IF (DEF_METHANE%q10_microbe_growth <= 0._r8 .or. DEF_METHANE%T_ref_microbe <= 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: microbial Q10 and reference temperature must be > 0: ', &
            DEF_METHANE%q10_microbe_growth, DEF_METHANE%T_ref_microbe
         bad = .true.
      ENDIF
      IF (DEF_METHANE%dormancy_rate_active < 0._r8 .or. DEF_METHANE%dormancy_rate_revive < 0._r8 .or. &
          DEF_METHANE%dormancy_threshold_methanogen_fS < 0._r8 .or. &
          DEF_METHANE%dormancy_threshold_methanogen_fO2 < 0._r8 .or. &
          DEF_METHANE%dormancy_threshold_methanotroph_fS < 0._r8 .or. &
          DEF_METHANE%dormancy_threshold_methanotroph_fO2 < 0._r8) THEN
         IF (p_is_master) write(6,*) '***** ERROR: microbial dormancy rates/thresholds must be >= 0'
         bad = .true.
      ENDIF
      IF (DEF_METHANE_hydrology%slopemax <= 0._r8 .or. &
          DEF_METHANE_hydrology%slopebeta >= 0._r8) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: methane hydrology requires slopemax > 0 and slopebeta < 0: ', &
            DEF_METHANE_hydrology%slopemax, DEF_METHANE_hydrology%slopebeta
         bad = .true.
      ENDIF
      IF (abs(DEF_METHANE_hydrology%vdcf - 2._r8) > 64._r8 * epsilon(1._r8) .or. &
          abs(DEF_METHANE_hydrology%pc - 0.4_r8) > 64._r8 * epsilon(1._r8)) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: retired methane hydrology knobs vdcf/pc must retain defaults 2.0/0.4: ', &
            DEF_METHANE_hydrology%vdcf, DEF_METHANE_hydrology%pc
         bad = .true.
      ENDIF
      IF (DEF_METHANE%use_transient_atm_methane .and. &
          (len_trim(DEF_METHANE%atm_methane_file) == 0 .or. &
           trim(tracer_lower(adjustl(DEF_METHANE%atm_methane_file))) == 'null')) THEN
         IF (p_is_master) write(6,*) &
            '***** ERROR: transient atmospheric CH4 requires atm_methane_file.'
         bad = .true.
      ENDIF

      IF (bad) CALL CoLM_Stop (' ***** ERROR: methane namelist validation failed')
   END SUBROUTINE validate_methane_namelist

   real(r8) FUNCTION methane_atm_mixing_ratio (year, month)
      integer, intent(in) :: year, month
      integer :: iy, im

      methane_atm_mixing_ratio = DEF_METHANE%atm_methane
      IF (.not. DEF_METHANE%use_transient_atm_methane) RETURN

      IF (year < 0 .or. year > 3000 .or. month < 1 .or. month > 12) THEN
         write(6,*) 'ERROR: atmospheric CH4 request is outside file-table bounds: ', year, month
         CALL CoLM_Stop ('invalid transient atmospheric CH4 date')
      ENDIF

      CALL load_methane_atm_file ()
      iy = year
      im = month
      IF (.not. ieee_is_finite(atm_ch4_file_molmol(iy,im)) .or. &
          atm_ch4_file_molmol(iy,im) <= 0._r8) THEN
         write(6,*) 'ERROR: transient atmospheric CH4 table has no valid value for ', iy, im
         CALL CoLM_Stop ('incomplete transient atmospheric CH4 table')
      ENDIF
      methane_atm_mixing_ratio = atm_ch4_file_molmol(iy,im)
   END FUNCTION methane_atm_mixing_ratio

   integer FUNCTION methane_history_accumulation_mode ()
      ! AccFlux has three coarse storage modes.  The mode marker distinguishes
      ! off/core/custom, while AccFlux also stores an exact selector fingerprint
      ! because distinct custom selectors can accumulate different families.
      character(len=4096) :: selector

      methane_history_accumulation_mode = 0
      IF (.not. DEF_METHANE%write_ch4_history) RETURN

      selector = tracer_lower(adjustl(trim(DEF_METHANE%ch4_history_vars)))
      SELECT CASE (trim(selector))
      CASE ('none', 'off', 'false', '.false.')
         methane_history_accumulation_mode = 0
      CASE ('core', 'default', 'minimal', 'fast')
         methane_history_accumulation_mode = 1
      CASE DEFAULT
         methane_history_accumulation_mode = 2
      END SELECT
   END FUNCTION methane_history_accumulation_mode

   logical FUNCTION methane_history_enabled (varname)
      ! Runtime per-variable gate for methane history output.
      ! See DEF_METHANE%ch4_history_vars in Methane_type.
      character(len=*), intent(in) :: varname
      character(len=4096), save :: cached_raw = ' '
      character(len=4096), save :: cached_list = ' '
      integer, save :: cached_mode = -1
      character(len=4096) :: list
      character(len=256)  :: v

      methane_history_enabled = .false.
      IF (.not. DEF_METHANE%write_ch4_history) RETURN

      list = tracer_lower(adjustl(trim(DEF_METHANE%ch4_history_vars)))
      v    = tracer_lower(adjustl(trim(varname)))

      ! History gates are queried dozens of times per history write.  Cache
      ! the normalized selector/list so custom comma lists do not get reparsed
      ! for every variable.  Rebuild automatically if the namelist value is
      ! changed between calls.
      IF (trim(list) /= trim(cached_raw)) THEN
         cached_raw  = list
         cached_list = compact_commas(trim(list))
         SELECT CASE (trim(cached_list))
         CASE ('all','*')
            cached_mode = 1
         CASE ('none','off','false','.false.')
            cached_mode = 0
         CASE ('core','default','minimal','fast')
            cached_mode = 2
         CASE ('diagnostic','extended','debug')
            cached_mode = 3
         CASE DEFAULT
            cached_mode = 4
         END SELECT
      ENDIF

      SELECT CASE (cached_mode)
      CASE (1)
         methane_history_enabled = .true.
      CASE (2)
         methane_history_enabled = methane_history_is_core(trim(v))
      CASE (3)
         methane_history_enabled = methane_history_is_diagnostic(trim(v))
      CASE (4)
         methane_history_enabled = index(','//trim(cached_list)//',', ','//trim(v)//',') > 0
      CASE DEFAULT
         methane_history_enabled = .false.
      END SELECT
   END FUNCTION methane_history_enabled

   logical FUNCTION methane_history_is_core (v)
      character(len=*), intent(in) :: v

      ! Keep the default history set intentionally small for global CH4 runs.
      ! Use ch4_history_vars='diagnostic' to recover the previous broad core.
      SELECT CASE (trim(v))
      CASE ( &
         'f_methane_surf_flux_tot', &
         'f_methane_surf_flux_tot_active', &
         'f_methane_surf_flux_global_total_with_lake', &
	     'f_methane_surf_flux_global_phys_with_lake', &
	     'f_methane_balance_residual_global_with_lake', &
	     'f_methane_ch4_clip_credit_global_with_lake', &
         'f_methane_surf_flux_tot_phys', &
         'f_methane_balance_residual', &
         'f_methane_ch4_clip_credit', &
         'f_o2_cap_loss', &
         'f_o2_cap_gain', &
         'f_methane_surf_flux_wetland', &
         'f_methane_surf_flux_soil', &
         'f_methane_surf_flux_lake', &
         'f_methane_surf_flux_rice', &
         'f_methane_prod_tot', &
         'f_methane_oxid_tot', &
         'f_totcol_methane')
         methane_history_is_core = .true.
      CASE DEFAULT
         methane_history_is_core = .false.
      END SELECT
   END FUNCTION methane_history_is_core

   logical FUNCTION methane_history_is_diagnostic (v)
      character(len=*), intent(in) :: v

      ! Broad diagnostic set preserved from the former 'core' selector.
      SELECT CASE (trim(v))
      CASE ( &
         'f_net_methane', &
         'f_methane_surf_flux_tot', &
         'f_methane_surf_flux_tot_active', &
         'f_methane_surf_flux_active_total_without_lake', &
         'f_methane_surf_flux_global_total_with_lake', &
	     'f_methane_surf_flux_global_phys_with_lake', &
	     'f_methane_balance_residual_global_with_lake', &
	     'f_methane_ch4_clip_credit_global_with_lake', &
         'f_methane_surf_flux_tot_phys', &
         'f_methane_surf_aere', &
         'f_methane_surf_aere_soil', &
         'f_methane_surf_aere_rice', &
         'f_methane_surf_ebul', &
         'f_methane_surf_ebul_soil', &
         'f_methane_surf_ebul_rice', &
         'f_methane_surf_diff', &
         'f_methane_surf_diff_soil', &
         'f_methane_surf_diff_rice', &
         'f_methane_surf_diff_phys', &
         'f_methane_balance_residual', &
         'f_methane_ch4_clip_credit', &
         'f_o2_cap_loss', &
         'f_o2_cap_gain', &
         'f_methane_prod_tot', &
         'f_methane_prod_tot_soil', &
         'f_methane_prod_tot_rice', &
         'f_methane_oxid_tot', &
         'f_methane_oxid_tot_soil', &
         'f_methane_oxid_tot_rice', &
         'f_co2_decomp_tot', &
         'f_co2_oxid_tot', &
         'f_co2_net_tot', &
         'f_totcol_methane', &
         'f_grnd_methane_cond', &
         'f_methane_surf_flux_tot_lake', &
         'f_methane_surf_ebul_lake', &
         'f_methane_surf_diff_lake', &
         'f_methane_prod_tot_lake', &
         'f_methane_oxid_tot_lake', &
         'f_co2_net_tot_lake', &
         'f_totcol_methane_lake', &
         'f_lake_water_ch4_stock', &
         'f_lake_water_o2_stock', &
         'f_lake_water_ch4_oxid', &
         'f_lake_sed_ch4_flux', &
         'f_lake_sed_o2_flux', &
         'f_lake_air_o2_flux', &
         'f_forc_pmethanem', &
         'f_layer_sat_lag', &
         'f_annavg_finrw', &
         'f_methane_dfsat_tot', &
         'f_f_h2osfc', &
         'f_methane_finundated', &
         'f_methane_soil_finundated', &
         'f_methane_soil_zwt', &
         'f_inund_flood_patch', &
         'f_inund_flood_depth_patch', &
         'f_wetland_frac_patch', &
         'f_methane_surf_flux_wetland', &
         'f_methane_surf_flux_soil', &
         'f_methane_surf_flux_lake', &
         'f_methane_surf_flux_lake_intensive', &
         'f_methane_surf_flux_rice', &
          'f_methane_surf_flux_rice_intensive', &
       ! Category-split CH4 budget components (wetland/soil/lake/rice).
          'f_methane_prod_tot_wetland', &
          'f_methane_oxid_tot_wetland', &
          'f_methane_surf_aere_wetland', &
          'f_methane_surf_ebul_wetland', &
          'f_methane_surf_diff_wetland', &
          'f_methane_area_wetland', &
          'f_methane_area_soil', &
          'f_methane_area_rice', &
          'f_methane_area_lake', &
          'f_methane_floodplain_frac', &
          'f_methane_wetland_type')
         methane_history_is_diagnostic = .true.
      CASE DEFAULT
         methane_history_is_diagnostic = .false.
      END SELECT
   END FUNCTION methane_history_is_diagnostic

   character(len=4096) FUNCTION compact_commas (s)
      ! Lower-level namelist convenience: allow users to write
      ! "f_a, f_b" with spaces after commas.
      character(len=*), intent(in) :: s
      integer :: i, n

      compact_commas = ' '
      n = 0
      DO i = 1, len_trim(s)
         IF (s(i:i) == ' ' .or. s(i:i) == char(9)) CYCLE
         n = n + 1
         IF (n <= len(compact_commas)) compact_commas(n:n) = s(i:i)
      ENDDO
   END FUNCTION compact_commas

   SUBROUTINE load_methane_atm_file ()
      character(len=512) :: line
      integer :: iu, ios, ios3, ios2, line_number
      integer :: yr, mon
      real(r8) :: val, val_molmol

      IF (atm_ch4_file_loaded) RETURN
      atm_ch4_file_loaded = .true.
      atm_ch4_file_molmol(:,:) = -1._r8

      open(newunit=iu, file=trim(DEF_METHANE%atm_methane_file), &
         status='old', action='read', iostat=ios)
      IF (ios /= 0) THEN
         write(6,*) 'ERROR: cannot open transient CH4 atmospheric file: ', &
            trim(DEF_METHANE%atm_methane_file)
         CALL CoLM_Stop ('cannot open transient atmospheric CH4 file')
      ENDIF

      line_number = 0
      DO
         read(iu,'(A)',iostat=ios) line
         IF (ios /= 0) EXIT
         line_number = line_number + 1
         line = adjustl(line)
         IF (len_trim(line) == 0) CYCLE
         IF (line(1:1) == '#' .or. line(1:1) == '!') CYCLE

         read(line,*,iostat=ios3) yr, mon, val
         IF (ios3 == 0) THEN
            IF (yr < 0 .or. yr > 3000 .or. mon < 1 .or. mon > 12) THEN
               write(6,*) 'ERROR: invalid atmospheric CH4 date at table line ', line_number
               CALL CoLM_Stop ('invalid transient atmospheric CH4 table date')
            ENDIF
            val_molmol = methane_atm_to_molmol(val)
            IF (.not. ieee_is_finite(val_molmol) .or. val_molmol <= 0._r8) THEN
               write(6,*) 'ERROR: invalid atmospheric CH4 value at table line ', line_number
               CALL CoLM_Stop ('invalid transient atmospheric CH4 value')
            ENDIF
            atm_ch4_file_molmol(yr,mon) = val_molmol
         ELSE
            read(line,*,iostat=ios2) yr, val
            IF (ios2 == 0 .and. yr >= 0 .and. yr <= 3000) THEN
               val_molmol = methane_atm_to_molmol(val)
               IF (.not. ieee_is_finite(val_molmol) .or. val_molmol <= 0._r8) THEN
                  write(6,*) 'ERROR: invalid atmospheric CH4 value at table line ', line_number
                  CALL CoLM_Stop ('invalid transient atmospheric CH4 value')
               ENDIF
               atm_ch4_file_molmol(yr,:) = val_molmol
            ELSE
               write(6,*) 'ERROR: malformed atmospheric CH4 table line ', line_number, ': ', trim(line)
               CALL CoLM_Stop ('malformed transient atmospheric CH4 table')
            ENDIF
         ENDIF
      ENDDO
      IF (ios > 0) THEN
         write(6,*) 'ERROR: failed while reading atmospheric CH4 table after line ', line_number
         CALL CoLM_Stop ('transient atmospheric CH4 table read failure')
      ENDIF
      close(iu)
   END SUBROUTINE load_methane_atm_file

   real(r8) FUNCTION methane_atm_to_molmol (val)
      real(r8), intent(in) :: val
      character(len=16) :: units

      units = adjustl(trim(DEF_METHANE%atm_methane_file_units))
      SELECT CASE (units)
      CASE ('mol/mol','molmol','vmr','MOL/MOL','MOLMOL','VMR')
         methane_atm_to_molmol = val
      CASE ('ppmv','ppm','PPMV','PPM')
         methane_atm_to_molmol = val * 1.e-6_r8
      CASE ('ppbv','ppb','PPBV','PPB')
         methane_atm_to_molmol = val * 1.e-9_r8
      CASE DEFAULT
         ! auto: accept common CH4 conventions: 1700 ppbv, 1.7 ppmv, or 1.7e-6 mol/mol.
         IF (val > 100._r8) THEN
            methane_atm_to_molmol = val * 1.e-9_r8
         ELSEIF (val > 0.1_r8) THEN
            methane_atm_to_molmol = val * 1.e-6_r8
         ELSE
            methane_atm_to_molmol = val
         ENDIF
      END SELECT
   END FUNCTION methane_atm_to_molmol
END MODULE MOD_Tracer_Reactive_Methane_Const
#endif
