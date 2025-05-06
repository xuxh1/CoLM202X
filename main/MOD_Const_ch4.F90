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

	real(r8), public, parameter :: rgasm = SHR_CONST_RGAS/1000._r8 ! J/mol.K; Universal gas constant 
	!!! rgas Different from CoLM
	!!! rgas in CoLM is gas constant for dry air [J/kg/K]
	!!! not Universal gas constant
	real(r8), parameter :: rgasLatm = 0.0821_r8 ! L.atm/mol.K

	real(r8), public, parameter :: secspday = 86400._r8 ! Seconds per day

   ! should set in const_ch4
   type, public :: params_type
      ! ch4 production constants
      real(r8) :: q10ch4 =1.33              ! additional Q10 for methane production ABOVE the soil decomposition temperature relationship (2+)
      real(r8) :: q10ch4base = 295._r8 ! temperature at which the effective f_ch4 actually equals the constant f_ch4 (295+)
      real(r8) :: f_ch4 = 0.2            ! ratio of CH4 production to total C mineralization (0.2+?) -------- defination differ with documentation
      ! ! real(r8) :: rootlitfrac        ! Fraction of soil organic matter associated with roots
      real(r8) :: cnscalefactor=1        ! scale factor on CN decomposition for assigning methane flux (?-)
      real(r8) :: redoxlag =30           ! Number of days to lag in the calculation of finundated_lag (30+)
      real(r8) :: lake_decomp_fact =9e-11    ! Base decomposition rate (1/s) at 25C (1)
      real(r8) :: redoxlag_vertical=0   ! time lag (days) to inhibit production for newly unsaturated layers (30+)
      real(r8) :: pHmax = 9._r8          ! maximum pH for methane production(= 9._r8)
      real(r8) :: pHmin = 2.2_r8         ! minimum pH for methane production(= 2.2_r8)
      real(r8) :: oxinhib = 400          ! inhibition of methane production by oxygen (m^3/mol) (400+?)

      real(r8) :: mino2lim = 0.2         ! minimum anaerobic decomposition rate as a fraction of potential aerobic rate (0.2+)

      ! ! ch4 oxidation constants
      real(r8) :: vmax_ch4_oxid = 45.e-6_r8 * 1000._r8 / 3600._r8       ! oxidation rate constant (= 45.e-6_r8 * 1000._r8 / 3600._r8) [mol/m3-w/s];
      real(r8) :: k_m = 5.e-6_r8 * 1000._r8                 ! Michaelis-Menten oxidation rate constant for CH4 concentration 
      real(r8) :: q10_ch4_oxid = 1         ! Q10 oxidation constant (?)
      real(r8) :: smp_crit =-2.4e5_r8            ! Critical soil moisture potential(mm)
      real(r8) :: k_m_o2 =20.e-6_r8 * 1000._r8              ! Michaelis-Menten oxidation rate constant for O2 concentration
      real(r8) :: k_m_unsat = 5.e-6_r8 * 1000._r8 / 10._r8           ! Michaelis-Menten oxidation rate constant for CH4 concentration
      real(r8) :: vmax_oxid_unsat = 45.e-6_r8 * 1000._r8 / 3600._r8 / 10._r8     ! (= 45.e-6_r8 * 1000._r8 / 3600._r8 / 10._r8) [mol/m3-w/s]

      ! ! ch4 aerenchyma constants
      real(r8) :: aereoxid =0            ! fraction of methane flux entering aerenchyma rhizosphere that will be(?)

      ! ! oxidized rather than emitted
      real(r8) :: scale_factor_aere = 1   ! scale factor on the aerenchyma area for sensitivity tests (1)
      ! real(r8) :: nongrassporosratio = 1/3   ! Ratio of root porosity in non-grass to grass, used for aerenchyma transport
      real(r8) :: unsat_aere_ratio = 0.05_r8 / 0.3_r8    ! Ratio to multiply upland vegetation aerenchyma porosity by compared to inundated systems (= 0.05_r8 / 0.3_r8)
      real(r8) :: porosmin = 0.05_r8            ! minimum aerenchyma porosity (unitless)(= 0.05_r8) 

      ! ! ch4 ebbulition constants
      real(r8) :: vgc_max  =0.15            ! ratio of saturation pressure triggering ebullition (0.15)

      ! ! ch4 transport constants
      real(r8) :: satpow  =2             ! exponent on watsat for saturated soil solute diffusion (2?)
      real(r8) :: scale_factor_gasdiff = 1! For sensitivity tests; convection would allow this to be > 1(?)
      real(r8) :: scale_factor_liqdiff = 1! For sensitivity tests; convection would allow this to be > 1(?)
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

END MODULE MOD_Const_ch4
    