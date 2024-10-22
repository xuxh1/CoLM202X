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

END MODULE MOD_Const_ch4
    