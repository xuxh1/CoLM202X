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
	use MOD_Vars_Global, only : maxsnl,nl_soil,nl_lake,spval,PI,deg2rad
	use MOD_Const_Physical, only: rgas, denh2o, denice, tfrz, grav
	use MOD_Const_ch4
    use MOD_ch4varcon
	!-----------------------------------------------------------------------
	implicit none
	save
	public  :: ch4
	
	private :: ch4_annualupdate
	private :: ch4_prod
	private :: ch4_oxid
	private :: ch4_aere
	private :: ch4_ebul
	private :: ch4_tran

contains