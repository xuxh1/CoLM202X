module MOD_ch4varcon

    !-----------------------------------------------------------------------
    ! Module containing CH4 parameters and logical switches and routine to read constants from CLM namelist.
    !
    ! use shr_kind_mod, only : r8 => shr_kind_r8
    ! use abortutils  , only : endrun
    ! use clm_varctl  , only : iulog
    ! use clm_varctl  , only : NLFileName_in
    !
    implicit none
    !
    ! Methane Model Parameters
    !
    private
  
    logical, public :: use_aereoxid_prog = .true. ! if false then aereoxid is read off of
    ! the parameter file and may be modifed by the user (default aereoxid on the
    ! file is 0.0).
  
    logical, public :: transpirationloss = .true. ! switch for activating CH4 loss from transpiration
                                        ! Transpiration loss assumes that the methane concentration in dissolved soil
                                        ! water remains constant through the plant and is released when the water evaporates
                                        ! from the stomata.
                                        ! Currently hard-wired to true; impact is < 1 Tg CH4/yr
  
    logical, public :: allowlakeprod = .false. ! Switch to allow production under lakes based on soil carbon dataset
                                       ! (Methane can be produced, and CO2 produced from methane oxidation,
                                       ! which will slowly reduce the available carbon stock, if ! replenishlakec, but no other biogeochem is done.)
                                       ! Note: switching this off turns off ALL lake methane biogeochem. However, 0 values
                                       ! will still be averaged into the concentration _sat history fields.
  
    logical, public :: usephfact = .false. ! Switch to use pH factor in methane production
  
    logical, public :: replenishlakec = .true. ! Switch for keeping carbon storage under lakes constant
                                        ! so that lakes do not affect the carbon balance
                                        ! Good for long term rather than transient warming experiments
                 ! NOTE SWITCHING THIS OFF ASSUMES TRANSIENT CARBON SUPPLY FROM LAKES; COUPLED MODEL WILL NOT CONSERVE CARBON
                 ! IN THIS MODE.
  
    ! inundatrion fraction -- which is used in methane code and potentially soil code
    integer, public :: finundation_mtd        ! Finundation method type to use, one of the following
    integer, public, parameter :: finundation_mtd_h2osfc        = 0 ! Use prognostic fsat h2osfc
    integer, public, parameter :: finundation_mtd_ZWT_inversion = 1 ! Use inversion of ZWT to Prigent satellite inundation obs. data
    integer, public, parameter :: finundation_mtd_TWS_inversion = 2 ! Use inversion of TWS to Prigent satellite inundation obs. data
  
    logical, public :: usefrootc = .false.    ! Use CLMCN fine root C rather than ann NPP & LAI based parameterization to
                                      ! calculate tiller C for aerenchyma area calculation.
                                      ! The NPP & LAI param. was based on Wania for Arctic sedges and may not be
                                      ! appropriate for woody Patches, although nongrassporosratio above partly adjusts
                                      ! for this.  However, using fine root C reduces the aerenchyma area by a large
                                      ! factor.
  
    logical, public :: ch4offline = .true.    ! true --> Methane is not passed between the land & atmosphere.
                                      ! NEM is not added to NEE flux to atm. to correct for methane production,
                                      ! and ambient CH4 is set to constant 2009 value.
  
    logical, public :: ch4rmcnlim = .false.   ! Remove the N and low moisture limitations on SOM HR when calculating
                                      ! methanogenesis.
                                      ! Note: this option has not been extensively tested.
                                      ! Currently hardwired off.
  
    logical, public :: anoxicmicrosites = .false. ! Use Arah & Stephen 1998 expression to allow production above the water table
                                          ! Currently hardwired off; expression is crude.
  
    logical, public :: ch4frzout = .false.    ! Exclude CH4 from frozen fraction of soil pore H2O, to simulate "freeze-out" pulse
                                      ! as in Mastepanov 2008.
                                      ! Causes slight increase in emissions in the fall and decrease in the spring.
                                      ! Currently hardwired off; small impact.
  
    ! public :: ch4conrd ! Read and initialize CH4 constants



!-----------------------------------------------------------------------

    integer, public :: iulog = 6        ! "stdout" log file unit number, default is 6
    logical, public :: use_cn              = .false.
    logical, public :: use_lch4            = .true.
    logical, public :: use_nitrif_denitrif = .true.
    logical, public :: use_fates_bgc = .false.                 ! true => use FATES along with CLM soil biogeochemistry

    ! true => anoxia is applied to heterotrophic respiration also considered in CH4 model
    ! default value reset in controlMod
    logical, public :: anoxia  = .true. 

end module MOD_ch4varcon