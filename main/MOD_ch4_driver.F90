#include <define.h>
#ifdef CH4
	SUBROUTINE ch4_driver (istep,i,idate,patchtype,deltim,lb,snl,dlon,dlat,&!input
		z_soisno,dz_soisno,zi_soisno,t_soisno,t_grnd,wliq_soisno,wice_soisno,&
		forc_t,forc_pbot,forc_po2m,forc_pco2m,&
		zwt,rootfr,snowdp,wat,rsur,etr,lakedepth,lake_icefrac,wdsrf,bsw,&
		smp,porsl,lai,rootr,fsatmax,fsatdcf,frcsat)

		use MOD_Precision
		use MOD_Const_Physical, only: rgas, denh2o, denice, tfrz, grav
		use MOD_Const_ch4
		! use MOD_ch4varcon
		use MOD_Namelist, only : DEF_USE_VariablySaturatedFlow
		use MOD_Vars_Global, only : maxsnl,nl_soil,nl_lake,spval,PI,deg2rad,z_soi,zi_soi,dz_soi
		use MOD_ch4
		use MOD_SPMD_Task
		USE MOD_LandPFT, only: patch_pft_s, patch_pft_e
		USE MOD_Vars_PFTimeInvariants,  only: pftfrac

		USE MOD_BGC_Vars_1DFluxes, only: decomp_hr, decomp_hr_vr, pot_f_nit_vr,&
		froot_mr, cpool_froot_gr, cpool_livecroot_gr, cpool_deadcroot_gr, &
		cpool_froot_storage_gr, cpool_livecroot_storage_gr, cpool_deadcroot_storage_gr, &
		transfer_froot_gr, transfer_livecroot_gr, transfer_deadcroot_gr, &
		somhr, lithr, hr_vr, rr, agnpp, bgnpp, annsum_npp, fphr

		USE MOD_BGC_Vars_TimeVariables, only: decomp_cpools_vr, o_scalar
	
		USE MOD_BGC_Vars_TimeVariables, only: &
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
		tempavg_agnpp, tempavg_bgnpp, annsum_counter, tempavg_somhr, tempavg_finrw, &
		fsat_bef, finundated_lag, ch4_dfsat_tot

		

		USE MOD_BGC_Vars_TimeInvariants, only: organic_max

		! USE MOD_BGC_Vars_PFTimeVariables, only: annsum_npp_p
		! USE MOD_BGC_Vars_1DPFTFluxes, only: froot_mr_p, cpool_froot_gr_p, cpool_froot_storage_gr_p, transfer_froot_gr_p
		IMPLICIT NONE
		integer ,intent(in) :: istep  

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
				rootfr     (1:nl_soil)        , &! fraction of roots in each soil layer
				snowdp                        , &! snow depth (m)
				wat                           , &! total water storage
				rsur                          , &! surface runoff (mm h2o/s)
				etr                           , &! transpiration rate [mm/s]
				lakedepth                     , &! lake depth (m)
				lake_icefrac(1:nl_lake)       , &! lake mass fraction of lake layer that is frozen
				wdsrf                         , &! depth of surface water [mm]
				bsw         (1:nl_soil)       , &! clapp and hornbereger "b" parameter [-]
				smp         (1:nl_soil)       , &! soil matrix potential [mm]
				porsl       (1:nl_soil)       , &! fraction of soil that is voids [-]
				lai                           , &! leaf area index
				rootr       (1:nl_soil)       , &! water exchange between soil and root. Positive: soil->root [?]

				fsatmax                       , &! maximum saturated area fraction [-]
				fsatdcf                       , &! decay factor in calculation of saturated area fraction [1/m]
        		frcsat                           ! fraction of saturation area

		integer :: ps, pe
		integer j
		real(r8):: &
				crootfr  (1:nl_soil)     , &! fraction of roots for carbon in each soil layer
				pH                       , &! soil water pH                                     
				cellorg  (1:nl_soil)     , &! column 3D org (kg/m3 organic matter)
				t_h2osfc             	    ! surface water temperature               

		ps = patch_pft_s(i)      
		pe = patch_pft_e(i)

		crootfr(:) = rootfr(:)
		pH = 7
		cellorg = 0.
		cellorg(:) = (cellorg(:) + sum(decomp_cpools_vr(1:10, 1:7, i), dim=2))*1000
		t_h2osfc = t_grnd

		CALL ch4 (istep,idate(1:3),patchtype,lb,snl,dlon,dlat,deltim,&
		z_soisno(maxsnl+1:),dz_soisno(maxsnl+1:),zi_soisno(maxsnl:),t_soisno(maxsnl+1:),&
		t_grnd,wliq_soisno(maxsnl+1:),wice_soisno(maxsnl+1:),&
		forc_t,forc_pbot,forc_po2m,forc_pco2m,&
		zwt,rootfr,snowdp,wat,rsur,etr,wdsrf,bsw,&
		smp,porsl,lai,rootr,&
		annsum_npp(i),rr(i),&
		fsatmax,fsatdcf,frcsat,&
		agnpp(i),bgnpp(i),somhr(i),&
		crootfr(1:nl_soil),lithr(i),hr_vr(1:nl_soil,i),o_scalar(1:nl_soil,i),fphr(1:nl_soil,i),pot_f_nit_vr(1:nl_soil,i),pH,&
		cellorg(1:nl_soil),t_h2osfc,organic_max,&
		!!!! --------------------------------------------------------------------------------------------------------
		!!!!                                         sum data   
		!!!! --------------------------------------------------------------------------------------------------------
		net_methane(i), &
		ch4_prod_depth(1:nl_soil,i), o2_decomp_depth(1:nl_soil,i), &
		ch4_oxid_depth(1:nl_soil,i), o2_oxid_depth(1:nl_soil,i), &
		ch4_aere_depth(1:nl_soil,i), ch4_tran_depth(1:nl_soil,i), o2_aere_depth(1:nl_soil,i), &
		ch4_ebul_depth(1:nl_soil,i), &
		o2stress(1:nl_soil,i), ch4stress(1:nl_soil,i), &
		ch4_surf_flux_tot(i), ch4_surf_aere(i), ch4_surf_ebul(i), ch4_surf_diff(i), &
		ch4_ebul_tot(i), ch4_prod_tot(i), ch4_oxid_tot(i), &
		totcolch4(i), grnd_ch4_cond(i), conc_o2(1:nl_soil,i), conc_ch4(1:nl_soil,i), &
		!!!! --------------------------------------------------------------------------------------------------------
		!!!! --------------------------------------------------------------------------------------------------------
		!!!!                                         sum data (unsaturated / saturated)
		!!!! --------------------------------------------------------------------------------------------------------
		net_methane_unsat(i), net_methane_sat(i), &
		ch4_prod_depth_unsat(1:nl_soil,i), ch4_prod_depth_sat(1:nl_soil,i), &
		o2_decomp_depth_unsat(1:nl_soil,i), o2_decomp_depth_sat(1:nl_soil,i), &
		ch4_oxid_depth_unsat(1:nl_soil,i), ch4_oxid_depth_sat(1:nl_soil,i), &
		o2_oxid_depth_unsat(1:nl_soil,i), o2_oxid_depth_sat(1:nl_soil,i), &
		ch4_aere_depth_unsat(1:nl_soil,i), ch4_aere_depth_sat(1:nl_soil,i), &
		ch4_tran_depth_unsat(1:nl_soil,i), ch4_tran_depth_sat(1:nl_soil,i), &
		o2_aere_depth_unsat(1:nl_soil,i), o2_aere_depth_sat(1:nl_soil,i), &
		ch4_ebul_depth_unsat(1:nl_soil,i), ch4_ebul_depth_sat(1:nl_soil,i), &
		o2stress_unsat(1:nl_soil,i), o2stress_sat(1:nl_soil,i), &
		ch4stress_unsat(1:nl_soil,i), ch4stress_sat(1:nl_soil,i), &
		ch4_surf_flux_tot_unsat(i), ch4_surf_flux_tot_sat(i), &
		ch4_surf_aere_unsat(i), ch4_surf_aere_sat(i), &
		ch4_surf_ebul_unsat(i), ch4_surf_ebul_sat(i), &
		ch4_surf_diff_unsat(i), ch4_surf_diff_sat(i), &
		ch4_ebul_tot_unsat(i), ch4_ebul_tot_sat(i), &
		ch4_prod_tot_unsat(i), ch4_prod_tot_sat(i), &
		ch4_oxid_tot_unsat(i), ch4_oxid_tot_sat(i), &
		totcolch4_unsat(i), totcolch4_sat(i), &
		grnd_ch4_cond_unsat(i), grnd_ch4_cond_sat(i), &
		conc_o2_unsat(1:nl_soil,i), conc_o2_sat(1:nl_soil,i), &
		conc_ch4_unsat(1:nl_soil,i), conc_ch4_sat(1:nl_soil,i), &
		!!!! --------------------------------------------------------------------------------------------------------
		c_atm(1:3,i), forc_pch4m(i), layer_sat_lag(1:nl_soil,i), lake_soilc(1:nl_soil,i), &
		annavg_agnpp(i), annavg_bgnpp(i), annavg_somhr(i), annavg_finrw(i), &
		tempavg_agnpp(i), tempavg_bgnpp(i), annsum_counter(i), tempavg_somhr(i), tempavg_finrw(i), fsat_bef(i), finundated_lag(i), ch4_dfsat_tot(i))

	END SUBROUTINE ch4_driver
#endif