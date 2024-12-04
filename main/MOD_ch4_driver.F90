#include <define.h>
#ifdef CH4
	SUBROUTINE ch4_driver (istep,i,idate,patchtype,deltim,lb,snl,dlon,dlat,&!input
		z_soisno,dz_soisno,zi_soisno,t_soisno,t_grnd,wliq_soisno,wice_soisno,&
		forc_t,forc_pbot,forc_po2m,forc_pco2m,&
		zwt,rootfr,snowdp,wat,rsur,etr,lakedepth,lake_icefrac,wdsrf,bsw,&
		smp,porsl,lai,rootr)

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
		use MOD_Vars_Global, only : maxsnl,nl_soil,nl_lake,spval,PI,deg2rad,z_soi,zi_soi,dz_soi
		use MOD_ch4_test
		use MOD_SPMD_Task
		USE MOD_LandPFT, only: patch_pft_s, patch_pft_e
		USE MOD_Vars_PFTimeInvariants,  only: pftfrac

		USE MOD_BGC_Vars_1DFluxes, only: decomp_hr, decomp_hr_vr, pot_f_nit_vr,&
		froot_mr, cpool_froot_gr, cpool_livecroot_gr, cpool_deadcroot_gr, &
		cpool_froot_storage_gr, cpool_livecroot_storage_gr, cpool_deadcroot_storage_gr, &
		transfer_froot_gr, transfer_livecroot_gr, transfer_deadcroot_gr, &
		somhr, lithr, hr_vr, rr, agnpp, bgnpp, annsum_npp, fphr

		USE MOD_BGC_Vars_TimeVariables, only: decomp_cpools_vr
	
		! USE MOD_BGC_Vars_TimeVariables, only: froot_mr, cpool_froot_gr, cpool_froot_storage_gr, transfer_froot_gr
		USE MOD_BGC_Vars_TimeVariables, only: c_atm, ch4_surf_flux_tot, net_methane,annavg_agnpp,annavg_bgnpp,&
		annavg_somhr,annavg_finrw,ch4_prod_depth,o2_decomp_depth,ch4_oxid_depth,o2_oxid_depth,&
		ch4_aere_depth,ch4_tran_depth,o2_aere_depth,ch4_ebul_depth,o2stress,ch4stress,ch4_surf_aere,&
		ch4_surf_ebul,ch4_surf_diff,ch4_ebul_total
	
		USE MOD_BGC_Vars_TimeVariables, only: totcolch4,forc_pch4m,grnd_ch4_cond,conc_o2,conc_ch4,&
		layer_sat_lag,lake_soilc,tempavg_agnpp,tempavg_bgnpp,annsum_counter,tempavg_somhr,tempavg_finrw, &
		o_scalar

		USE MOD_BGC_Vars_TimeInvariants, only: organic_max

		! USE MOD_BGC_Vars_PFTimeVariables, only: annsum_npp_p
		! USE MOD_BGC_Vars_1DPFTFluxes, only: froot_mr_p, cpool_froot_gr_p, cpool_froot_storage_gr_p, transfer_froot_gr_p
		IMPLICIT NONE
		integer ,intent(inout) :: istep  

		integer ,intent(in) :: i         ! patch index
		integer ,intent(in) :: idate(1:3)  ! current date (year, day of the year, seconds of the day)
		integer ,intent(in) :: patchtype ! land patch type (0=soil, 1=urban and built-up,
								         ! 2=wetland, 3=land ice, 4=land water bodies, 99 = ocean)
		real(r8),intent(in) :: deltim    ! time step in seconds
		integer ,intent(in) :: &
				lb,&
				snl
		real(r8),intent(in) :: dlon     ! longitude (degrees)
		real(r8),intent(in) :: dlat     ! latitude (degrees)

		real(r8),intent(in) :: &
				z_soisno   (maxsnl+1:nl_soil) , &! layer depth (m)
				dz_soisno  (maxsnl+1:nl_soil) , &! layer thickness (m)
				zi_soisno  (maxsnl:nl_soil)   , &! interface level below a "z" level (m)
				t_soisno   (maxsnl+1:nl_soil) , &! soil + snow layer temperature [K]
				t_grnd                 		  , &! ground surface temperature [k]
				wliq_soisno(maxsnl+1:nl_soil) , &! liquid water (kg/m2)
				wice_soisno(maxsnl+1:nl_soil) , &! ice lens (kg/m2)
				forc_t                        , &! temperature at agcm reference height [kelvin]
				forc_pbot                     , &! atmosphere pressure at the bottom of the atmos. model level [pa]
				forc_po2m                     , &! partial pressure of O2 at observational height [pa]
				forc_pco2m                    , &! partial pressure of CO2 at observational height [pa]
				zwt                           , &! the depth to water table [m]
				rootfr     (1:nl_soil)          , &! fraction of roots in each soil layer
				snowdp                        , &! snow depth (m)
				wat                           , &! total water storage
				rsur                          , &! surface runoff (mm h2o/s)
				etr                           , &! transpiration rate [mm/s]
				lakedepth                     , &! lake depth (m)
				lake_icefrac(1:nl_lake)         , &! lake mass fraction of lake layer that is frozen
				wdsrf                         , &! depth of surface water [mm]
				bsw         (1:nl_soil)         , &! clapp and hornbereger "b" parameter [-]
				smp         (1:nl_soil)         , &! soil matrix potential [mm]
				porsl       (1:nl_soil)         , &! fraction of soil that is voids [-]
				lai                           , &! leaf area index
				rootr       (1:nl_soil)            ! water exchange between soil and root. Positive: soil->root [?]

		integer :: ps, pe
		integer j
		logical :: ch4_first_time = .false.
		! real(r8):: &
				! annsum_npp_tmp
				! froot_mr,&
				! cpool_froot_gr,&
				! cpool_froot_storage_gr,&
				! transfer_froot_gr
		real(r8):: &
				! agnpp                   , &! aboveground NPP (gC/m2/s)
				! bgnpp                   , &! belowground NPP (gC/m2/s)
				! rr                      , &! root respiration (fine root MR + total root GR) (gC/m2/s)
				! somhr                   , &! (gC/m2/s) soil organic matter heterotrophic respiration
				! lithr                   , &! (gC/m2/s) litter heterotrophic respiration        
				! hr_vr    (1:nl_soil)      , &! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
				crootfr  (1:nl_soil)      , &! fraction of roots for carbon in each soil layer
				! o_scalar (1:nl_soil)      , &! fraction by which decomposition is limited by anoxia
				! fphr     (1:nl_soil)      , &! fraction of potential heterotrophic respiration 
				pH                      , &! soil water pH                                     
				cellorg  (1:nl_soil)      , &! column 3D org (kg/m3 organic matter)
				t_h2osfc               	     ! surface water temperature               
				! organic_max                ! organic matter content (kg/m3) where soil is assumed to act like peat

		ps = patch_pft_s(i)      
		pe = patch_pft_e(i)

		! agnpp = annsum_npp(i)/365/86400/2
		! bgnpp = annsum_npp(i)/365/86400/2

		! annsum_npp_tmp = 600
		! agnpp = annsum_npp_tmp/365/86400/2
		! bgnpp = annsum_npp_tmp/365/86400/2

		! rr = froot_mr(i) + cpool_froot_gr(i) + cpool_livecroot_gr(i) + cpool_deadcroot_gr(i) + &
		!  cpool_froot_storage_gr(i) + cpool_livecroot_storage_gr(i) + cpool_deadcroot_storage_gr(i) + &
		!  transfer_froot_gr(i) + transfer_livecroot_gr(i) + transfer_deadcroot_gr(i)

		! somhr = decomp_hr(i)/2
		! lithr = decomp_hr(i)/2
		! hr_vr(1:10) = sum(decomp_hr_vr(1:10,1:10,i),dim=2)
		crootfr(:) = rootfr(:)
		! o_scalar(:) = 1
		! fphr(:) = 1
		pH = 7

		cellorg = 0.
		cellorg(:) = cellorg(:) + sum(decomp_cpools_vr(1:10, 1:7, i), dim=2)
		t_h2osfc = (t_grnd+forc_t)/2
		! organic_max = cellorg(1)

		if (istep == 1) then
			ch4_first_time = .true.
		endif

		CALL ch4 (i,idate(1:3),patchtype,lb,nl_soil,maxsnl,snl,dlon,dlat,deltim,&
		z_soisno(maxsnl+1:),dz_soisno(maxsnl+1:),zi_soisno(maxsnl:),t_soisno(maxsnl+1:),&
		t_grnd,wliq_soisno(maxsnl+1:),wice_soisno(maxsnl+1:),&
		forc_t,forc_pbot,forc_po2m,forc_pco2m,&
		zwt,rootfr,snowdp,wat,rsur,etr,lakedepth,lake_icefrac,wdsrf,bsw,&
		smp,porsl,lai,rootr,&
		annsum_npp(i),rr(i),agnpp(i),bgnpp(i),somhr(i),&
		crootfr(1:nl_soil),lithr(i),hr_vr(1:nl_soil,i),o_scalar(1:nl_soil,i),fphr(1:nl_soil,i),pot_f_nit_vr(1:nl_soil,i),pH,&
		cellorg(1:nl_soil),t_h2osfc,organic_max,&
		c_atm(1:3,i),ch4_surf_flux_tot(i),net_methane(i),annavg_agnpp(i),annavg_bgnpp(i),annavg_somhr(i),annavg_finrw(i),&
		ch4_prod_depth(1:nl_soil,i),o2_decomp_depth(1:nl_soil,i),&
		ch4_oxid_depth(1:nl_soil,i),o2_oxid_depth(1:nl_soil,i),&
		ch4_aere_depth(1:nl_soil,i),ch4_tran_depth(1:nl_soil,i),o2_aere_depth(1:nl_soil,i),&
		ch4_ebul_depth(1:nl_soil,i),&
		o2stress(1:nl_soil,i),ch4stress(1:nl_soil,i),ch4_surf_aere(i),ch4_surf_ebul(i),ch4_surf_diff(i),ch4_ebul_total(i),&
		ch4_first_time,totcolch4(i),forc_pch4m(i),grnd_ch4_cond(i),conc_o2(1:nl_soil,i),conc_ch4(1:nl_soil,i),&
		layer_sat_lag(1:nl_soil,i),lake_soilc(1:nl_soil,i),&
		tempavg_agnpp(i),tempavg_bgnpp(i),annsum_counter(i),&
		tempavg_somhr(i),tempavg_finrw(i))
	
	END SUBROUTINE ch4_driver
#endif