#include <define.h>
#ifdef BGC

   SUBROUTINE bgc_driver (i,idate,deltim,dlat,dlon)

!-----------------------------------------------------------------------------------------------------------
! !DESCRIPTION:
! The trunk subroutine of the CoLM biogeochemistry module. The bgc_driver link different bgc processes, and
! sequentially run each process step by step. bgc_driver was called by CoLMDRIVER includes vegetation
! and soil CN cycle processes.
! 
!
! !ORIGINAL:
! The Community Land Model version 5.0 (CLM5.0)
!
! !REFERENCE:
! Lawrence, D.M., Fisher, R.A., Koven, C.D., Oleson, K.W., Swenson, S.C., Bonan, G., Collier, N., 
! Ghimire, B., van Kampenhout, L., Kennedy, D. and Kluzek, E., 2019. 
! The Community Land Model version 5: Description of new features, benchmarking,
! and impact of forcing uncertainty. Journal of Advances in Modeling Earth Systems, 11(12), 4245-4287.
!
! !REVISION:
! Xingjie Lu, 2022, modify original CLM5 to be compatible with CoLM code structure.


   USE MOD_Precision
   USE MOD_Namelist, only : DEF_USE_SASU, DEF_USE_NITRIF, DEF_USE_CNSOYFIXN, DEF_USE_FIRE, DEF_USE_IRRIGATION
   USE MOD_Const_Physical, only : tfrz, denh2o, denice
   USE MOD_Vars_PFTimeInvariants, only: pftfrac
   USE MOD_LandPFT, only: patch_pft_s, patch_pft_e
   USE MOD_BGC_Vars_1DFluxes, only: plant_ndemand, ndep_to_sminn
   USE MOD_BGC_Vars_1DPFTFluxes, only: plant_ndemand_p, cpool_to_leafc_p, crop_seedc_to_leaf_p
   USE MOD_BGC_Veg_CNMResp, only: CNMResp
   USE MOD_BGC_Soil_BiogeochemDecompCascadeBGC, only: decomp_rate_constants_bgc
   USE MOD_BGC_Soil_BiogeochemPotential, only: SoilBiogeochemPotential
   USE MOD_BGC_Soil_BiogeochemVerticalProfile, only: SoilBiogeochemVerticalProfile
   USE MOD_BGC_Veg_NutrientCompetition, only: calc_plant_nutrient_demand_CLM45_default,&
                                           calc_plant_nutrient_competition_CLM45_default
   USE MOD_BGC_Soil_BiogeochemNitrifDenitrif, only: SoilBiogeochemNitrifDenitrif

   USE MOD_BGC_Soil_BiogeochemCompetition, only: SoilBiogeochemCompetition
   USE MOD_BGC_Soil_BiogeochemDecomp, only: SoilBiogeochemDecomp
   USE MOD_BGC_Veg_CNPhenology, only: CNPhenology
   USE MOD_BGC_Veg_CNGResp, only: CNGResp
   USE MOD_BGC_CNCStateUpdate1, only: CStateUpdate1
   USE MOD_BGC_CNNStateUpdate1, only: NStateUpdate1
   USE MOD_BGC_Soil_BiogeochemNStateUpdate1, only: SoilBiogeochemNStateUpdate1
   USE MOD_BGC_Soil_BiogeochemLittVertTransp, only: SoilBiogeochemLittVertTransp
   USE MOD_BGC_Veg_CNGapMortality, only: CNGapMortality
   USE MOD_BGC_CNCStateUpdate2, only: CStateUpdate2
   USE MOD_BGC_CNNStateUpdate2, only: NStateUpdate2
   USE MOD_BGC_CNCStateUpdate3, only: CStateUpdate3
   USE MOD_BGC_Soil_BiogeochemNLeaching, only: SoilBiogeochemNLeaching
   USE MOD_BGC_CNNStateUpdate3, only: NstateUpdate3
   USE MOD_BGC_CNSummary, only: CNDriverSummarizeStates, CNDriverSummarizeFluxes
   USE MOD_BGC_CNAnnualUpdate, only: CNAnnualUpdate
   USE MOD_BGC_CNZeroFluxes, only: CNZeroFluxes
   USE MOD_BGC_Veg_CNVegStructUpdate, only: CNVegStructUpdate
   USE MOD_BGC_CNBalanceCheck, only: BeginCNBalance, CBalanceCheck, NBalanceCheck
   USE MOD_BGC_CNSASU, only: CNSASU
   USE MOD_BGC_Veg_CNNDynamics, only: CNNFixation
#ifdef CROP
   USE MOD_BGC_Veg_CNNDynamics, only: CNNFert, CNSoyfix
   USE MOD_Irrigation, only: CalIrrigationNeeded
#endif
   USE MOD_TimeManager
   USE MOD_Vars_Global, only: nl_soil, nl_soil_full, ndecomp_pools, ndecomp_pools_vr, ndecomp_transitions, npcropmin, &
                       z_soi,dz_soi,zi_soi,nbedrock,zmin_bedrock

   USE MOD_BGC_Vars_TimeVariables, only: sminn_vr, col_begnb, skip_balance_check, decomp_cpools_vr
   USE MOD_BGC_Veg_CNFireBase, only: CNFireFluxes
   USE MOD_BGC_Veg_CNFireLi2016, only: CNFireArea

   IMPLICIT NONE

   integer ,intent(in) :: i        ! patch index
   real(r8),intent(in) :: deltim   ! time step in seconds
   integer ,intent(in) :: idate(3) ! current date (year, day of the year, seconds of the day)
   real(r8),intent(in) :: dlat     ! latitude (degrees)
   real(r8),intent(in) :: dlon     ! longitude (degrees)

   integer :: ps, pe
   integer j

   !   --------------------------------------------------------------ch4----------------------------------------------
!   ---------------------------------------------------------------------------------------------------------------
#ifdef CH4
! set input
real(r8) :: &
annsum_npp              , &! annual sum NPP (gC/m2/yr)
rr                      , &! root respiration (fine root MR + total root GR) (gC/m2/s)
agnpp                   , &! aboveground NPP (gC/m2/s)
bgnpp                   , &! belowground NPP (gC/m2/s)
somhr,&
crootfr(1:nl_soil),lithr,hr_vr(1:nl_soil),o_scalar(1:nl_soil),fphr(1:nl_soil),pot_f_nit_vr(1:nl_soil),pH,&
cellorg(1:nl_soil),t_h2osfc,organic_max
! set output
real(r8), INTENT(out) :: &
c_atm(1:3),&
ch4_surf_flux_tot            , &! CH4 flux to atm. (kg C/m**2/s)
net_methane                  , &! average net methane correction to CO2 flux (g C/m^2/s)
annavg_agnpp,&
annavg_bgnpp,&
annavg_somhr,&
annavg_finrw,&
ch4_prod_depth(1:nl_soil),o2_decomp_depth(1:nl_soil),&
ch4_oxid_depth(1:nl_soil),o2_oxid_depth(1:nl_soil),&
ch4_aere_depth(1:nl_soil),ch4_tran_depth(1:nl_soil),o2_aere_depth(1:nl_soil),&
ch4_ebul_depth(1:nl_soil),&
o2stress(1:nl_soil),ch4stress(1:nl_soil),ch4_surf_aere,ch4_surf_ebul,ch4_surf_diff,ch4_ebul_total
! set inout
logical ::&
ch4_first_time
real(r8) :: &
totcolch4,&
forc_pch4m              , &! CH4 concentration in atmos. (pascals)
grnd_ch4_cond,&
conc_o2  (1:nl_soil)    , &! O2 conc in each soil layer (mol/m3) 
conc_ch4(1:nl_soil),&
layer_sat_lag(1:nl_soil),&
lake_soilc(1:nl_soil),&
tempavg_agnpp,&
tempavg_bgnpp,&
annsum_counter,&
tempavg_somhr,&
tempavg_finrw
#endif   

      ps = patch_pft_s(i)      
      pe = patch_pft_e(i)
      CALL BeginCNBalance(i)
      CALL CNZeroFluxes(i, ps, pe, nl_soil, ndecomp_pools, ndecomp_transitions)
      CALL CNNFixation(i,idate)
      CALL CNMResp(i, ps, pe, nl_soil, npcropmin)
      CALL decomp_rate_constants_bgc(i,nl_soil,z_soi)
      CALL SoilBiogeochemPotential(i,nl_soil,ndecomp_pools,ndecomp_transitions)
      CALL SoilBiogeochemVerticalProfile(i,ps,pe,nl_soil,nl_soil_full,nbedrock,zmin_bedrock,z_soi,dz_soi)
      IF(DEF_USE_NITRIF)THEN
         CALL SoilBiogeochemNitrifDenitrif(i,nl_soil,dz_soi)
      ENDIF
      CALL calc_plant_nutrient_demand_CLM45_default(i,ps,pe,deltim,npcropmin)
    
      plant_ndemand(i) = sum( plant_ndemand_p(ps:pe)*pftfrac(ps:pe) )
  
      CALL SoilBiogeochemCompetition(i,deltim,nl_soil,dz_soi)
      CALL calc_plant_nutrient_competition_CLM45_default(i,ps,pe,npcropmin)
#ifdef CROP
      IF(DEF_USE_CNSOYFIXN)THEN
         CALL CNSoyfix (i, ps, pe, nl_soil)
      ENDIF
#endif
  
      CALL SoilBiogeochemDecomp(i,nl_soil,ndecomp_pools,ndecomp_transitions, dz_soi)
      CALL CNPhenology(i,ps,pe,nl_soil,idate(1:3),dz_soi,deltim,dlat,npcropmin,phase=1)
      CALL CNPhenology(i,ps,pe,nl_soil,idate(1:3),dz_soi,deltim,dlat,npcropmin,phase=2)
#ifdef CROP
      CALL CNNFert(i, ps, pe)
#endif
      CALL CNGResp(i, ps, pe, npcropmin)
#ifdef CROP
      IF(DEF_USE_IRRIGATION)THEN
         CALL CalIrrigationNeeded(i,ps,pe,idate,nl_soil,nbedrock,z_soi,dz_soi,deltim,dlon,npcropmin)
      ENDIF
#endif
      ! update vegetation pools from phenology, allocation and nitrogen uptake
      ! update soil N pools from decomposition and nitrogen competition
      CALL CStateUpdate1(i, ps, pe, deltim, nl_soil, ndecomp_transitions, npcropmin)
      CALL NStateUpdate1(i, ps, pe, deltim, nl_soil, ndecomp_transitions, npcropmin,dz_soi)
      CALL SoilBiogeochemNStateUpdate1(i,deltim,nl_soil,ndecomp_transitions,dz_soi)
  
      ! update soil pools from soil vertical mixing
      CALL SoilBiogeochemLittVertTransp(i,deltim,nl_soil,nl_soil_full,ndecomp_pools,nbedrock,z_soi,zi_soi,dz_soi)
  
      ! update vegetation pools from gap mortality
      CALL CNGapMortality(i, ps, pe, nl_soil,npcropmin)
      CALL CStateUpdate2(i, ps, pe, deltim, nl_soil)
      CALL NStateUpdate2(i, ps, pe, deltim, nl_soil, dz_soi)
  
      IF(DEF_USE_FIRE)THEN
         ! update vegetation and fire pools from fire
         CALL CNFireArea(i,ps,pe,dlat,nl_soil,idate,dz_soi)
         CALL CNFireFluxes(i,ps,pe,dlat,nl_soil,ndecomp_pools)
      ENDIF   
      CALL CStateUpdate3(i,ps,pe,deltim, nl_soil, ndecomp_pools)
      CALL CNAnnualUpdate(i,ps,pe,deltim,idate(1:3))
  
  ! update soil mineral nitrogen pools leaching
      CALL SoilBiogeochemNLeaching(i,deltim,nl_soil,zi_soi,dz_soi)
      CALL NstateUpdate3(i, ps, pe, deltim, nl_soil, ndecomp_pools,dz_soi)
  
      IF(DEF_USE_SASU)THEN
         CALL CNSASU(i,ps,pe,deltim,idate(1:3),nl_soil,ndecomp_transitions,ndecomp_pools,ndecomp_pools_vr)! only for spin up
      ENDIF
  
      CALL CNDriverSummarizeStates(i,ps,pe,nl_soil,dz_soi,ndecomp_pools, .false.)
      CALL CNDriverSummarizeFluxes(i,ps,pe,nl_soil,dz_soi,ndecomp_transitions,ndecomp_pools,deltim)
  
      IF( .not. skip_balance_check(i) )THEN
         CALL CBalanceCheck(i,ps,pe,deltim,dlat,dlon)
         CALL NBalanceCheck(i,deltim,dlat,dlon)
  
  
      ELSE
         skip_balance_check(i) = .false.
      ENDIF
  
      CALL CNVegStructUpdate(i,ps,pe,deltim,npcropmin)


#ifdef CH4
! set input
annsum_npp = 1
rr = 1

agnpp  = 1
bgnpp  = 1
somhr  = 1

crootfr = 1
lithr = 1
hr_vr = 1
o_scalar = 1
fphr = 1
pot_f_nit_vr =1
pH=1

cellorg =1
t_h2osfc = 1
organic_max = 1

!  set inout
ch4_first_time = .true.
totcolch4 = 1
forc_pch4m = 1
grnd_ch4_cond = 1
conc_o2  = 1 
conc_ch4 = 1
layer_sat_lag = 1
lake_soilc=1

tempavg_agnpp = 1
tempavg_bgnpp = 1
annsum_counter = 1
tempavg_somhr =1 
tempavg_finrw = 1
CALL ch4 (ipatch,patchtype,&!input
patchlonr,patchlatr,&
lb,nl_soil,maxsnl,snl,&
deltim,&
z_soisno(lb:nl_soil),dz_soisno(lb:nl_soil),zi_soisno(lb-1:nl_soil),t_soisno(lb:nl_soil),t_grnd,wliq_soisno(lb:nl_soil),wice_soisno(lb:nl_soil),&
forc_t,forc_pbot,forc_po2m,forc_pco2m,&
zwt,rootfr(1:nl_soil),snowdp,wat,rsur,etr,lakedepth,lake_icefrac(1:nl_soil),wdsrf,bsw(1:nl_soil),&
smp(1:nl_soil),porsl(1:nl_soil),lai,&
annsum_npp,rr,&
idate,agnpp,bgnpp,somhr,&
crootfr(1:nl_soil),lithr,hr_vr(1:nl_soil),o_scalar(1:nl_soil),fphr(1:nl_soil),pot_f_nit_vr(1:nl_soil),pH,&
rootr(1:nl_soil),&
cellorg(1:nl_soil),t_h2osfc,organic_max,&
c_atm(1:3),ch4_surf_flux_tot,net_methane,&!output
annavg_agnpp,annavg_bgnpp,annavg_somhr,annavg_finrw,&
ch4_prod_depth(1:nl_soil),o2_decomp_depth(1:nl_soil),&
ch4_oxid_depth(1:nl_soil),o2_oxid_depth(1:nl_soil),&
ch4_aere_depth(1:nl_soil),ch4_tran_depth(1:nl_soil),o2_aere_depth(1:nl_soil),&
ch4_ebul_depth(1:nl_soil),&
o2stress(1:nl_soil),ch4stress(1:nl_soil),ch4_surf_aere,ch4_surf_ebul,ch4_surf_diff,ch4_ebul_total,&
ch4_first_time,totcolch4,forc_pch4m,grnd_ch4_cond,conc_o2(1:nl_soil),conc_ch4(1:nl_soil),layer_sat_lag(1:nl_soil),lake_soilc(1:nl_soil),&!inout
tempavg_agnpp,tempavg_bgnpp,annsum_counter,&
tempavg_somhr,tempavg_finrw)
!  print*, ch4_surf_flux_tot,net_methane
!  print*, totcolch4,grnd_ch4_cond
!  print*, forc_pch4m,layer_sat_lag
#endif

   END SUBROUTINE bgc_driver

#endif
