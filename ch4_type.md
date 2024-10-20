#input:
    ipatch,patchtype,&
i,patchtype,
    patchlonr,patchlatr,&
dlat,dlon,
    lb,nl_soil,maxsnl,snl,&
,nl_soil,maxsnl,snl,
    deltim,&
deltim,
    z_soisno,dz_soisno,zi_soisno,t_soisno,t_grnd,wliq_soisno,wice_soisno,&
z_soi,dz_soi,zi_soi,t_soi,t_grnd...
    forc_t,forc_pbot,forc_po2m,forc_pco2m,&
    zwt,rootfr,snowdp,wat,rsur,etr,lakedepth,lake_icefrac,wdsrf,bsw,&
    smp,porsl,lai,&


    annsum_npp,rr,&
    idate,agnpp,bgnpp,somhr,&
    
    crootfr,lithr,hr_vr,o_scalar,fphr,pot_f_nit_vr,pH,&
    rootr,&
    cellorg,t_h2osfc,organic_max,&

#inout:
    ch4_first_time,totcolch4,forc_pch4m,grnd_ch4_cond,conc_o2,conc_ch4,layer_sat_lag,lake_soilc,&
    tempavg_agnpp,tempavg_bgnpp,annsum_counter,&
    tempavg_somhr,tempavg_finrw

#output:
    c_atm,ch4_surf_flux_tot,net_methane,&
    annavg_agnpp,annavg_bgnpp,annavg_somhr,annavg_finrw,&
    ch4_prod_depth,o2_decomp_depth,&
    ch4_oxid_depth,o2_oxid_depth,&
    ch4_aere_depth,ch4_tran_depth,o2_aere_depth,&
    ch4_ebul_depth,&
    o2stress,ch4stress,ch4_surf_aere,ch4_surf_ebul,ch4_surf_diff,ch4_ebul_total,&

