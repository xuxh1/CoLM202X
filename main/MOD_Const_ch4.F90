MODULE MOD_Const_ch4
!=======================================================================
! ch4 constants
!=======================================================================
	USE MOD_Precision
	USE MOD_ForcingDownscaling, only: SHR_CONST_RGAS
	IMPLICIT NONE
	
	PUBLIC

	!------------------------------------------------------------------
	! (Non-tunable) Constants for the CH4 submodel (Tuneable constants in ch4varcon)
	!------------------------------------------------------------------
	! Note some of these constants are also used in CNNitrifDenitrifMod

	integer, private :: i  ! loop index

	integer, public, parameter :: ngases      =   3     ! CH4, O2, & CO2

	!------------------------------------------------------------------

	real(r8), public, parameter :: catomw = 12.011_r8 ! molar mass of C atoms (g/mol)

	real(r8), public :: s_con(ngases,4)    ! Schmidt # calculation constants (spp, #)
	data (s_con(1,i),i=1,4) /1898_r8, -110.1_r8, 2.834_r8, -0.02791_r8/ ! CH4
	data (s_con(2,i),i=1,4) /1801_r8, -120.1_r8, 3.7818_r8, -0.047608_r8/ ! O2
	data (s_con(3,i),i=1,4) /1911_r8, -113.7_r8, 2.967_r8, -0.02943_r8/ ! CO2
	
	real(r8), public :: d_con_w(ngases,3)    ! water diffusivity constants (spp, #)  (mult. by 10^-4)
	data (d_con_w(1,i),i=1,3) /0.9798_r8, 0.02986_r8, 0.0004381_r8/ ! CH4
	data (d_con_w(2,i),i=1,3) /1.172_r8, 0.03443_r8, 0.0005048_r8/ ! O2
	data (d_con_w(3,i),i=1,3) /0.939_r8, 0.02671_r8, 0.0004095_r8/ ! CO2
	
	real(r8), public :: d_con_g(ngases,2)    ! gas diffusivity constants (spp, #) (cm^2/s) (mult. by 10^-9)
	data (d_con_g(1,i),i=1,2) /0.1875_r8, 0.0013_r8/ ! CH4
	data (d_con_g(2,i),i=1,2) /0.1759_r8, 0.00117_r8/ ! O2
	data (d_con_g(3,i),i=1,2) /0.1325_r8, 0.0009_r8/ ! CO2
	
	real(r8), public :: c_h_inv(ngases)    ! constant (K) for Henry's law (4.12, Wania)
	data c_h_inv(1:3) /1600._r8, 1500._r8, 2400._r8/ ! CH4, O2, CO2
	
	real(r8), public :: kh_theta(ngases)    ! Henry's constant (L.atm/mol) at standard temperature (298K)
	data kh_theta(1:3) /714.29_r8, 769.23_r8, 29.4_r8/ ! CH4, O2, CO2
	
	real(r8), public :: kh_tbase = 298._r8 ! base temperature for calculation of Henry's constant (K)

!------------------------------------------------------------------

	real(r8), public, parameter :: rgasm = SHR_CONST_RGAS/1000._r8 ! J mol-1 K-1; Universal gas constant 
                                 ![J/K/mol]=[J/K/kmol]/[kmol/mol]
	!!! rgas Different from CoLM
	!!! rgas in CoLM is gas constant for dry air [J/kg/K]
	!!! not Universal gas constant
	real(r8), public, parameter :: rgasLatm = 0.0821_r8 ! L.atm/mol.K

	real(r8), public, parameter :: secspday = 86400._r8 ! Seconds per day

!------------------------------------------------------------------
   !!!!!!!!!!!!!!!!!!!!!!! The parameters here are inaccurate and need to be reconfirmed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! params_type
   ! ch4 production constants
   real(r8), public, parameter :: q10ch4 =1.33              ! additional Q10 for methane production ABOVE the soil decomposition temperature relationship (doc:Q10 Baseline:2 Range:1.5~4) (params:1.33)
   real(r8), public, parameter :: q10lake =q10ch4*1.5_r8    ! For now, take to be the same as q10ch4 * 1.5.

   real(r8), public, parameter :: q10ch4base = 295._r8 ! temperature at which the effective f_ch4 actually equals the constant f_ch4 (295+ params:295)
   real(r8), public, parameter :: f_ch4 = 0.2            ! ratio of CH4 production to total C mineralization (Baseline:0.2 Range:NA params:0.2) 
                                                         ! -------- defination differ with documentation
   ! ! real(r8) :: rootlitfrac        ! Fraction of soil organic matter associated with roots (params:0.5)
   real(r8), public, parameter :: cnscalefactor=1.        ! scale factor on CN decomposition for assigning methane flux (?- params:1.)
   real(r8), public, parameter :: redoxlag =30.           ! Number of days to lag in the calculation of finundated_lag (30+ params:30.)
   real(r8), public, parameter :: lake_decomp_fact =9e-11    ! Base decomposition rate (1/s) at 25C (1 params:9e-11)
   real(r8), public, parameter :: redoxlag_vertical=30._r8   ! time lag (days) to inhibit production for newly unsaturated layers (30+ params:0.)
   real(r8), public, parameter :: pHmax = 9._r8          ! maximum pH for methane production(params:9. code:9.)
   real(r8), public, parameter :: pHmin = 2.2_r8         ! minimum pH for methane production(params:2.2 code:2.2)
   real(r8), public, parameter :: oxinhib = 400._r8          ! inhibition of methane production by oxygen (m^3/mol) (400+? params:400.)

   real(r8), public, parameter :: mino2lim = 0.2_r8         ! minimum anaerobic decomposition rate as a fraction of potential aerobic rate (0.2+ params:0.2)

   ! ! ch4 oxidation constants
   real(r8), public, parameter :: vmax_ch4_oxid = 45.e-6_r8 * 1000._r8 / 3600._r8       ! oxidation rate constant (params:1.25e-5 code:45.e-6_r8 * 1000._r8 / 3600._r8) (doc:Ro,max Baseline:1.25e-5 Range:1.25e-6~1.25e-4 Unit:mol m-3 s-1)
   real(r8), public, parameter :: k_m = 5.e-6_r8 * 1000._r8                 ! Michaelis-Menten oxidation rate constant for CH4 concentration (params:5e-3 code:5.e-6_r8 * 1000._r8) (doc:KCH4 Baseline:5e-3 Range:5e-4~5e-2)
   real(r8), public, parameter :: q10_ch4_oxid = 1.9_r8         ! Q10 oxidation constant (? params:1.9)
   real(r8), public, parameter :: smp_crit =-2.4e5_r8            ! Critical soil moisture potential (mm) (params:-2.4e5)
   real(r8), public, parameter :: k_m_o2 =20.e-6_r8 * 1000._r8              ! Michaelis-Menten oxidation rate constant for O2 concentration (params:2e-2 code:20.e-6_r8 * 1000._r8) (doc:KO2 Baseline:2e-2 Range:2e-3~2e-1)
   real(r8), public, parameter :: k_m_unsat = 5.e-6_r8 * 1000._r8 / 10._r8           ! Michaelis-Menten oxidation rate constant for CH4 concentration (params:5e-4 code:5.e-6_r8 * 1000._r8 / 10._r8) (doc:KCH4 Baseline:5e-3 Range:5e-4~5e-2)
   real(r8), public, parameter :: vmax_oxid_unsat = 45.e-6_r8 * 1000._r8 / 3600._r8 / 10._r8     ! (params:1.25e-6 code:45.e-6_r8 * 1000._r8 / 3600._r8 / 10._r8) (doc:Ro,max Baseline:1.25e-5 Range:1.25e-6~1.25e-4 Unit:mol m-3 s-1)

   ! ! ch4 aerenchyma constants
   real(r8), public, parameter :: aereoxid =0._r8            ! fraction of methane flux entering aerenchyma rhizosphere that will be(? params:0.)

   ! ! oxidized rather than emitted
   real(r8), public, parameter :: scale_factor_aere = 1._r8   ! scale factor on the aerenchyma area for sensitivity tests (1 params:1.) (doc:Fa Baseline:1 Range:0.5~1.5)
   ! real(r8) :: nongrassporosratio = 1/3   ! Ratio of root porosity in non-grass to grass, used for aerenchyma transport (params:0.33)
   real(r8), public, parameter :: unsat_aere_ratio = 0.05_r8 / 0.3_r8    ! Ratio to multiply upland vegetation aerenchyma porosity by compared to inundated systems (params:0.1666666667 code:0.05_r8 / 0.3_r8)
   real(r8), public, parameter :: porosmin = 0.05_r8            ! minimum aerenchyma porosity (unitless)(params:0.05 code:0.05_r8) 

   ! ! ch4 ebbulition constants
   real(r8), public, parameter :: vgc_max  =0.15_r8            ! ratio of saturation pressure triggering ebullition (params:0.15) (doc:Ce,max Ce,min Baseline:0.15 Unit:mol m-3 ?)

   ! ! ch4 transport constants
   real(r8), public, parameter :: satpow  =2._r8             ! exponent on watsat for saturated soil solute diffusion (2? params:2.)
   real(r8), public, parameter :: scale_factor_gasdiff = 1   ! For sensitivity tests; convection would allow this to be > 1(? params:1.) (doc:fD0 Basline:1 Range:1,10 Unit:m2 s-1)
   real(r8), public, parameter :: scale_factor_liqdiff = 1   ! For sensitivity tests; convection would allow this to be > 1(? params:1.) (doc:fD0 Basline:1 Range:1,10 Unit:m2 s-1)
   real(r8), public, parameter :: capthick = 100._r8         ! min thickness before assuming h2osfc is impermeable (mm) (params:100.code:100._r8)

   ! ! additional constants
   ! real(r8) :: f_sat =0.95               ! volumetric soil water defining top of water table or where production is allowed (params:0.95 code:0.95)
   real(r8), public, parameter :: qflxlagd  = 30._r8          ! days to lag qflx_surf_lag in the tropics (days) (params:30 code:30._r8)
   real(r8), public, parameter :: highlatfact = 2._r8         ! multiple of qflxlagd for high latitudes	(params:2. code:2._r8)	
   real(r8), public, parameter :: q10lakebase = 298._r8       ! (K) base temperature for lake CH4 production (params:298. code:298._r8)
   real(r8), public, parameter :: atmch4  = 1.7e-6_r8         ! Atmospheric CH4 mixing ratio to prescribe if not provided by the atmospheric model (params:1.7e-6 code:1.7e-6_r8 search:1.9e-6) (mol/mol) could change with year
   real(r8), public, parameter :: rob = 3._r8                 ! ratio of root length to vertical depth ("root obliquity") (params:3. code:3._r8)
   real(r8), public, parameter :: om_frac_sf = 1._r8          ! Scale factor for organic matter fraction (unitless)(? params:NA)

END MODULE MOD_Const_ch4
    