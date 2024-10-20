	real(r8), INTENT(in) :: &
        annsum_npp              , &! annual sum NPP (gC/m2/yr)
annsum_npp_p
        rr                      , &! root respiration (fine root MR + total root GR) (gC/m2/s)
需要具体区分

        agnpp                   , &! aboveground NPP (gC/m2/s)
        bgnpp                   , &! belowground NPP (gC/m2/s)
自己加
        crootfr  (1:nl_soil)    , &! fraction of roots for carbon in each soil layer

        somhr                   , &! (gC/m2/s) soil organic matter heterotrophic respiration
        lithr                   , &! (gC/m2/s) litter heterotrophic respiration        
        hr_vr    (1:nl_soil)    , &! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
        o_scalar (1:nl_soil)    , &! fraction by which decomposition is limited by anoxia
        fphr     (1:nl_soil)    , &! fraction of potential heterotrophic respiration 

        pot_f_nit_vr(1:nl_soil) , &! (gN/m3/s) potential soil nitrification flux 
        pH                      , &! soil water pH 
        rootr    (1:nl_soil)    , &! effective fraction of roots in each soil layer (SMS method only)
        ! rootr here for effective per-layer transpiration, which may not be the same as rootfr

        cellorg  (1:nl_soil)    , &! column 3D org (kg/m^3 organic matter)
        organic_max                ! organic matter content (kg m-3) where soil is assumed to act like peat


    real(r8), INTENT(inout) :: &
        forc_pch4m              , &! CH4 concentration in atmos. (pascals)
        conc_o2  (1:nl_soil)    , &! O2 conc in each soil layer (mol/m3) 
        conc_ch4   (1:nl_soil)  , &! CH4 conc in each soil layer (mol/m3) 
        lake_soilc  (1:nl_soil) , &! total soil organic matter found in level (g C / m^3) (nl_soil)

        grnd_ch4_cond           , &! tracer conductance for boundary layer [m/s]

        totcolch4                  ! total methane in soil column, start of timestep (g C / m^2)
