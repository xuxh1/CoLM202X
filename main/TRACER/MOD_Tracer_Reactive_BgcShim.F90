#include <define.h>

#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Reactive_BgcShim
!=======================================================================
! Reactive tracer wetland/BGC coupling shim.
!
! This module is the reactive-tracer boundary for invoking the BGC
! decomposition cascade needed by wetland CH4.  Reactive_Methane_Impl
! should orchestrate methane driver calls, not reach directly into
! MOD_BGC_* internals.
!=======================================================================

   USE MOD_Precision
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   USE MOD_SPMD_Task, only: CoLM_stop, p_is_master
   USE MOD_Vars_Global, only: nl_soil, z_soi, dz_soi, &
      ndecomp_pools, ndecomp_transitions
   USE MOD_LandPatch, only: numpatch
   USE MOD_BGC_Soil_BiogeochemDecompCascadeBGC, only: decomp_rate_constants_bgc
   USE MOD_BGC_Soil_BiogeochemPotential,        only: SoilBiogeochemPotential
   USE MOD_BGC_Soil_BiogeochemCompetition,      only: SoilBiogeochemCompetitionNoPlant
   USE MOD_BGC_Soil_BiogeochemDecomp,           only: SoilBiogeochemDecomp
   USE MOD_BGC_Vars_1DFluxes, only: decomp_cpools_sourcesink, decomp_npools_sourcesink, &
      decomp_hr_vr, decomp_ctransfer_vr, &
      decomp_ntransfer_vr, decomp_sminn_flux_vr, sminn_to_denit_decomp_vr, &
      pmnf_decomp, p_decomp_cpool_loss, net_nmin_vr, gross_nmin_vr, &
      net_nmin, gross_nmin, potential_immob_vr, phr_vr, pot_f_nit_vr, &
      decomp_hr, somc_fire, som_c_leached, som_n_leached, denit, f_n2o_nit, &
      smin_no3_leached, smin_no3_runoff, sminn_leached, sminn_to_plant
   USE MOD_BGC_Vars_TimeVariables, only: fpi_vr, o_scalar, t_scalar, w_scalar, &
      depth_scalar, decomp_k, decomp_cpools_vr, decomp_npools_vr
   USE MOD_Namelist, only: DEF_BGC_DEBUG_SCALARS, DEF_WETLAND_PEAT_C_PASSIVE
   USE MOD_Vars_TimeVariables, only: t_soisno, smp, wliq_soisno, wice_soisno, zwt
   USE MOD_Vars_TimeInvariants, only: porsl, patchclass, patchtype
   USE MOD_BGC_Vars_TimeInvariants, only: i_met_lit, i_cel_lit, i_lig_lit, i_soil3, tau_s3, &
      i_cwd, i_soil1, i_soil2, donor_pool, receiver_pool, rf_decomp, pathfrac_decomp, &
      floating_cn_ratio, initial_cn_ratio, is_litter
   USE MOD_Vars_1DFluxes, only: assim
   USE MOD_Const_LC, only: rootfr
   USE MOD_Const_Physical, only: tfrz
   USE MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE
   USE MOD_BGC_Soil_BiogeochemVerticalProfile, only: surfprof_exp
   USE MOD_Tracer_Reactive_Methane_WetlandVeg, only: wetveg_bg_frac, wetland_root_profile, wetveg_forest_input_ratio
   USE MOD_Tracer_Reactive_Methane_State, only: methane_finundated, root_exudate_vr

   IMPLICIT NONE
   PRIVATE

   ! Call counter, used only to throttle the decomposition-scalar dump.
   integer :: dbg_calls = 0

   ! Wetland semi-analytic spin-up (DEF_METHANE%wetland_bgc_sasu). The sums
   ! cover one spin-up cycle and are cleared at every jump. They need no
   ! restart: CoLM replays every spin-up cycle from the spin-up start date
   ! inside one run, so a cycle is always summed whole.
   logical  :: sasu_active = .false.
   real(r8), allocatable :: sasu_time    (:)       ! summed time [s]
   real(r8), allocatable :: sasu_c_acc   (:,:,:)   ! sum of C dt [gC m-3 s]
   real(r8), allocatable :: sasu_kc_acc  (:,:,:)   ! sum of k C dt, k before fpi [gC m-3]
   real(r8), allocatable :: sasu_k_acc   (:,:,:)   ! sum of k dt [-]
   real(r8), allocatable :: sasu_ic_acc  (:,:,:)   ! plant C input [gC m-3]
   real(r8), allocatable :: sasu_in_acc  (:,:,:)   ! plant N input [gN m-3]
   real(r8), allocatable :: sasu_lit_pot (:)       ! column litter C loss before fpi [gC m-2]
   real(r8), allocatable :: sasu_lit_act (:)       ! column litter C loss after fpi [gC m-2]
   ! A jumped pool denser than this [gC m-3], 17 times pure organic matter
   ! (1 g cm-3 at 580 gC/kg), marks a failed solve; the layer is left alone.
   real(r8), parameter :: SASU_C_MAX = 1.e7_r8

   PUBLIC :: reactive_bgc_run_wetland_decomp
   PUBLIC :: wetland_bgc_sasu_init
   PUBLIC :: wetland_bgc_sasu_cycle_end


CONTAINS

   SUBROUTINE reactive_bgc_set_wetland_anoxia (ipatch)

! !DESCRIPTION:
!  Anoxia limitation on decomposition, which CoLM202X otherwise lacks.
!
!  The wetland hydrology holds every thawed layer at saturation, so those layers
!  are anaerobic and decompose at mino2lim of the aerobic rate -- the parameter
!  already carries exactly that definition. Without it a saturated tile has no
!  limiter left at all: t_scalar tracks temperature, depth_scalar is fixed, and
!  w_scalar is 1 precisely because the tile is waterlogged. Frozen layers keep 1
!  so the suppression is not counted twice against t_scalar and w_scalar --
!  unless frozen_anoxic_decomp is on: just below freezing w_scalar is
!  still near 1, and 1 there made a freezing layer decompose five times faster
!  than a thawed one.
!
!  Set BEFORE decomposition runs -- decomp_rate_constants_bgc folds o_scalar into
!  decomp_k inside itself, and its patchtype 2 exemption is what preserves this.
!
!  CALLERS MUST RESTRICT THIS TO patchtype == 2. It lived inside the wetland
!  decomposition shim until 2026-08-02, and that shim is called for every patch,
!  so an upland soil was silently limited to 20% of its aerobic rate too. No
!  single-point run could show it -- all 44 towers are patchtype 2 -- but a
!  global run decomposes its entire land surface through this.
!
!  With wetland_anoxia_catotelm >= 0 the scalar also deepens below the
!  water table: mino2lim at the table, falling towards the catotelm floor
!  over wetland_anoxia_efold (HPM). The inundated share has its table at the
!  surface, the rest at the host water table; above the table it stays at
!  mino2lim, which with CENTURY's litter rates matches peatland litter bags.

      IMPLICIT NONE
      integer, intent(in) :: ipatch
      integer :: j
      real(r8) :: o_top, fsat

      IF (.not. allocated(o_scalar)) RETURN
      o_top = max(DEF_METHANE%mino2lim, 1.e-6_r8)
      fsat = 0._r8
      IF (allocated(methane_finundated)) fsat = min(max(methane_finundated(ipatch), 0._r8), 1._r8)
      DO j = 1, nl_soil
         IF (t_soisno(j,ipatch) > tfrz .or. DEF_METHANE%frozen_anoxic_decomp) THEN
            IF (DEF_METHANE%wetland_anoxia_catotelm >= 0._r8) THEN
               o_scalar(j,ipatch) = fsat * anoxia_below_table(z_soi(j)) &
                  + (1._r8 - fsat) * anoxia_below_table(max(z_soi(j) - zwt(ipatch), 0._r8))
            ELSE
               o_scalar(j,ipatch) = o_top
            ENDIF
         ELSE
            o_scalar(j,ipatch) = 1._r8
         ENDIF
      ENDDO

   CONTAINS

      real(r8) FUNCTION anoxia_below_table (d)
         real(r8), intent(in) :: d    ! depth below the water table [m]
         anoxia_below_table = DEF_METHANE%wetland_anoxia_catotelm &
            + (o_top - DEF_METHANE%wetland_anoxia_catotelm) * exp(-d / DEF_METHANE%wetland_anoxia_efold)
      END FUNCTION anoxia_below_table

   END SUBROUTINE reactive_bgc_set_wetland_anoxia

   SUBROUTINE reactive_bgc_run_wetland_decomp (ipatch, deltim)

      IMPLICIT NONE
      integer, intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      integer :: j, k
      real(r8) :: litter_nrate, plant_ndemand_col
      real(r8) :: exu_vr(1:nl_soil)   ! root exudate respiration [gC m-3 s-1]

      IF (.not. ieee_is_finite(deltim) .or. deltim <= 0._r8) THEN
         CALL CoLM_stop(' ***** ERROR: wetland CH4/BGC coupling requires a finite positive timestep')
      ENDIF

      ! nothing is exuded unless the plant input below sets it
      IF (DEF_METHANE%root_exudate_frac > 0._r8 .and. allocated(root_exudate_vr)) &
         root_exudate_vr(:,ipatch) = 0._r8

      ! Start from the same clean per-patch flux state as the full BGC driver.
      IF (allocated(decomp_cpools_sourcesink))   decomp_cpools_sourcesink  (1:nl_soil,:,ipatch) = 0._r8
      IF (allocated(decomp_npools_sourcesink))   decomp_npools_sourcesink  (1:nl_soil,:,ipatch) = 0._r8
      IF (allocated(decomp_hr_vr))              decomp_hr_vr             (1:nl_soil,:,ipatch) = 0._r8
      IF (allocated(decomp_ctransfer_vr))       decomp_ctransfer_vr      (1:nl_soil,:,ipatch) = 0._r8
      IF (allocated(decomp_ntransfer_vr))       decomp_ntransfer_vr      (1:nl_soil,:,ipatch) = 0._r8
      IF (allocated(decomp_sminn_flux_vr))      decomp_sminn_flux_vr     (1:nl_soil,:,ipatch) = 0._r8
      IF (allocated(sminn_to_denit_decomp_vr))  sminn_to_denit_decomp_vr (1:nl_soil,:,ipatch) = 0._r8
      IF (allocated(pmnf_decomp))               pmnf_decomp              (1:nl_soil,:,ipatch) = 0._r8
      IF (allocated(p_decomp_cpool_loss))       p_decomp_cpool_loss      (1:nl_soil,:,ipatch) = 0._r8
      IF (allocated(net_nmin_vr))               net_nmin_vr              (1:nl_soil,ipatch)   = 0._r8
      IF (allocated(gross_nmin_vr))             gross_nmin_vr            (1:nl_soil,ipatch)   = 0._r8
      IF (allocated(potential_immob_vr))        potential_immob_vr       (1:nl_soil,ipatch)   = 0._r8
      IF (allocated(phr_vr))                    phr_vr                   (1:nl_soil,ipatch)   = 0._r8
      IF (allocated(pot_f_nit_vr))              pot_f_nit_vr             (1:nl_soil,ipatch)   = 0._r8
      IF (allocated(o_scalar))                  o_scalar                 (1:nl_soil,ipatch)   = 1._r8
      IF (allocated(fpi_vr))                    fpi_vr                   (1:nl_soil,ipatch)   = 1._r8
      IF (allocated(net_nmin))                  net_nmin                 (ipatch)             = 0._r8
      IF (allocated(gross_nmin))                gross_nmin               (ipatch)             = 0._r8
      IF (allocated(decomp_hr))                 decomp_hr                (ipatch)             = 0._r8
      IF (allocated(somc_fire))                 somc_fire                (ipatch)             = 0._r8
      IF (allocated(som_c_leached))             som_c_leached            (ipatch)             = 0._r8
      IF (allocated(som_n_leached))             som_n_leached            (ipatch)             = 0._r8
      IF (allocated(denit))                     denit                    (ipatch)             = 0._r8
      IF (allocated(f_n2o_nit))                 f_n2o_nit                (ipatch)             = 0._r8
      IF (allocated(smin_no3_leached))          smin_no3_leached         (ipatch)             = 0._r8
      IF (allocated(smin_no3_runoff))           smin_no3_runoff          (ipatch)             = 0._r8
      IF (allocated(sminn_leached))             sminn_leached            (ipatch)             = 0._r8
      IF (allocated(sminn_to_plant))            sminn_to_plant           (ipatch)             = 0._r8

      ! Plant carbon input of the tile, on the source/sink
      ! that CDecompStateUpdate adds to the pools once they are not fixed.
      litter_nrate = 0._r8
      exu_vr(:) = 0._r8
      IF (DEF_METHANE%wetland_plant_input) CALL wetland_plant_litter_input (ipatch, deltim, litter_nrate, exu_vr)

      ! Anoxia limit on the waterlogged tile, only when BGC owns it under the
      ! CH4/BGC contract (bgc_anoxia_limits_decomp); otherwise o_scalar stays 1.
      IF (DEF_METHANE%bgc_anoxia_limits_decomp) CALL reactive_bgc_set_wetland_anoxia (ipatch)

      CALL decomp_rate_constants_bgc (ipatch, nl_soil, z_soi)
      ! the peat passive pool turns over on its own base time.
      IF (DEF_METHANE%wetland_tau_s3 > 0._r8) &
         decomp_k(1:nl_soil,i_soil3,ipatch) = decomp_k(1:nl_soil,i_soil3,ipatch) &
            * tau_s3 / DEF_METHANE%wetland_tau_s3
      CALL SoilBiogeochemPotential   (ipatch, nl_soil, ndecomp_pools, ndecomp_transitions)
      ! Track B N closure: with wetland_n_uptake the litter N above is taken up
      ! from the tile's mineral N instead of appearing from nowhere; with
      ! wetland_n_unlimited immobilization is not held back by the layer's
      ! mineral N. Both off leave the call as before (demand 0, limited).
      IF (DEF_METHANE%wetland_n_uptake) THEN
         plant_ndemand_col = litter_nrate
         CALL SoilBiogeochemCompetitionNoPlant (ipatch, deltim, nl_soil, dz_soi, &
            DEF_METHANE%wetland_n_unlimited, plant_ndemand_col)
      ELSE
         CALL SoilBiogeochemCompetitionNoPlant (ipatch, deltim, nl_soil, dz_soi, &
            DEF_METHANE%wetland_n_unlimited)
      ENDIF
      CALL SoilBiogeochemDecomp      (ipatch, nl_soil, ndecomp_pools, ndecomp_transitions, dz_soi)

      ! root exudates are respired where released; booked on the
      ! metabolic-litter transition so the CH4 side counts them as litter HR.
      ! SoilBiogeochemDecomp assigns decomp_hr_vr, so this comes after it.
      IF (DEF_METHANE%wetland_exudate_frac > 0._r8) THEN
         DO k = 1, ndecomp_transitions
            IF (donor_pool(k) == i_met_lit) THEN
               decomp_hr_vr(1:nl_soil,k,ipatch) = decomp_hr_vr(1:nl_soil,k,ipatch) + exu_vr(:)
               EXIT
            ENDIF
         ENDDO
      ENDIF

      ! Semi-analytic spin-up: sum this step for the jump at the cycle's end.
      IF (sasu_active) CALL wetland_bgc_sasu_accumulate (ipatch, deltim)

      ! Decomposition-scalar dump.  decomp_k is the product of these five
      ! multipliers with a per-pool base rate, so printing the multipliers
      ! says which one differs between a site that emits and one that does
      ! not, without having to reason about the pool index.
      IF (DEF_BGC_DEBUG_SCALARS) CALL dump_decomp_scalars (ipatch, deltim)

   END SUBROUTINE reactive_bgc_run_wetland_decomp

   SUBROUTINE wetland_plant_litter_input (ipatch, deltim, nrate, exu_vr)

      ! Canopy assimilation [mol CO2 m-2 s-1] times NPP/GPP gives the carbon
      ! the tile's plants return to the soil at steady state; it enters the
      ! litter pools along the land class's root profile with CoLM's grass
      ! litter split and a grass litter C:N (mean of leaf litter 50 and fine
      ! roots 42 in MOD_Const_PFT). With wetland_bg_frac < 1 only that
      ! share follows the roots; the aboveground rest is laid on the surface
      ! along CoLM's leaf-litter profile, exp(-surfprof_exp z).
      ! With wetland_forest_input_herb the forested share enters as a herb
      ! layer at the non-forested share's LAI, with the herb belowground share;
      ! with wetland_dim_class moss-carpeted peat keeps at least its
      ! moss layer, wetland_moss_input_frac of the stand input.
      ! With wetland_burial_frac_site that share of the litter goes to soil3.
      ! nrate returns the litter N added this step [gN m-2 s-1], buried litter
      ! included, which wetland_n_uptake asks the soil to supply.
      IMPLICIT NONE
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      real(r8), intent(out) :: nrate
      real(r8), intent(out) :: exu_vr(1:nl_soil)   ! exudate respiration [gC m-3 s-1]
      integer  :: j
      real(r8) :: cin, cbur, s, bg, prof(1:nl_soil), sprof(1:nl_soil)

      nrate = 0._r8
      exu_vr(:) = 0._r8
      IF (.not. allocated(decomp_cpools_sourcesink)) RETURN
      cin = assim(ipatch)
      IF (.not. ieee_is_finite(cin) .or. cin <= 0._r8 .or. cin > 1.e-2_r8) RETURN
      cin = cin * 12.011_r8 * DEF_METHANE%wetland_npp_frac * deltim       ! g C m-2 per step
      IF (DEF_METHANE%wetland_forest_input_herb) cin = cin * wetveg_forest_input_ratio(ipatch)
      ! Vascular cover of a marsh tower: the open water of the
      ! footprint gives no plant input; the exudate share below follows.
      IF (DEF_METHANE%wetland_cover_frac_site >= 0._r8) cin = cin * DEF_METHANE%wetland_cover_frac_site

      CALL wetland_root_profile (ipatch, prof)

      ! the exudate share is respired where the roots release it and
      ! leaves the litter input.
      IF (DEF_METHANE%wetland_exudate_frac > 0._r8) THEN
         exu_vr(:) = cin * DEF_METHANE%wetland_exudate_frac * prof(:) / dz_soi(1:nl_soil) / deltim
         cin = cin * (1._r8 - DEF_METHANE%wetland_exudate_frac)
      ! where wetland_exudate_frac is 0: the same share of the input leaves the
      ! litter input, but it is booked on decomp_hr_vr only after the pools
      ! are updated (tracer_ch4_bgc_finalize_step); CDecompStateUpdate would
      ! otherwise debit the metabolic litter for carbon it never received.
      ELSEIF (DEF_METHANE%root_exudate_frac > 0._r8 .and. allocated(root_exudate_vr)) THEN
         root_exudate_vr(:,ipatch) = cin * DEF_METHANE%root_exudate_frac * prof(:) &
            / dz_soi(1:nl_soil) / deltim
         cin = cin * (1._r8 - DEF_METHANE%root_exudate_frac)
      ENDIF
      IF (allocated(decomp_npools_sourcesink)) nrate = cin / DEF_METHANE%wetland_litter_cn / deltim

      ! Managed marshes: the buried share of the litter becomes new
      ! peat in the passive pool, along the root profile and with its litter
      ! N, so nrate above still covers all the organic N added. Exudates are
      ! not litter and keep their share of the whole input.
      IF (DEF_METHANE%wetland_burial_frac_site > 0._r8) THEN
         cbur = cin * DEF_METHANE%wetland_burial_frac_site
         cin  = cin - cbur
         DO j = 1, nl_soil
            decomp_cpools_sourcesink(j,i_soil3,ipatch) = decomp_cpools_sourcesink(j,i_soil3,ipatch) &
               + cbur * prof(j) / dz_soi(j)
            IF (allocated(decomp_npools_sourcesink)) &
               decomp_npools_sourcesink(j,i_soil3,ipatch) = decomp_npools_sourcesink(j,i_soil3,ipatch) &
                  + cbur * prof(j) / dz_soi(j) / DEF_METHANE%wetland_litter_cn
         ENDDO
      ENDIF

      bg = wetveg_bg_frac(ipatch)      ! wetland_bg_frac, forest-weighted under wetland_veg_glwd
      IF (DEF_METHANE%wetland_forest_input_herb) bg = DEF_METHANE%wetland_bg_frac
      IF (bg < 1._r8) THEN
         sprof(:) = exp(-surfprof_exp * z_soi(1:nl_soil)) * dz_soi(1:nl_soil)
         sprof(:) = sprof(:) / sum(sprof)
         prof(:) = bg * prof(:) + (1._r8 - bg) * sprof(:)
      ENDIF

      DO j = 1, nl_soil
         decomp_cpools_sourcesink(j,i_met_lit,ipatch) = decomp_cpools_sourcesink(j,i_met_lit,ipatch) &
            + cin * prof(j) / dz_soi(j) * 0.25_r8
         decomp_cpools_sourcesink(j,i_cel_lit,ipatch) = decomp_cpools_sourcesink(j,i_cel_lit,ipatch) &
            + cin * prof(j) / dz_soi(j) * 0.50_r8
         decomp_cpools_sourcesink(j,i_lig_lit,ipatch) = decomp_cpools_sourcesink(j,i_lig_lit,ipatch) &
            + cin * prof(j) / dz_soi(j) * 0.25_r8
         IF (allocated(decomp_npools_sourcesink)) THEN
            decomp_npools_sourcesink(j,i_met_lit,ipatch) = decomp_npools_sourcesink(j,i_met_lit,ipatch) &
               + cin * prof(j) / dz_soi(j) * 0.25_r8 / DEF_METHANE%wetland_litter_cn
            decomp_npools_sourcesink(j,i_cel_lit,ipatch) = decomp_npools_sourcesink(j,i_cel_lit,ipatch) &
               + cin * prof(j) / dz_soi(j) * 0.50_r8 / DEF_METHANE%wetland_litter_cn
            decomp_npools_sourcesink(j,i_lig_lit,ipatch) = decomp_npools_sourcesink(j,i_lig_lit,ipatch) &
               + cin * prof(j) / dz_soi(j) * 0.25_r8 / DEF_METHANE%wetland_litter_cn
         ENDIF
      ENDDO

   END SUBROUTINE wetland_plant_litter_input

   SUBROUTINE dump_decomp_scalars (ipatch, deltim)

      IMPLICIT NONE
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      integer :: j, every

      ! roughly monthly, whatever the timestep
      every = max(1, nint(30._r8 * 86400._r8 / deltim))
      dbg_calls = dbg_calls + 1
      IF (mod(dbg_calls, every) /= 1) RETURN

      write(6,'(A,I8,A,I6)') 'BGCSCAL step=', dbg_calls, ' patch=', ipatch
      ! smp and the liquid/ice split are here because w_scalar turned out to be
      ! the multiplier that differs, and the floor it sits on is a test on smp:
      ! knowing w_scalar is 0.001 does not say whether the soil is genuinely at
      ! -200 m or whether smp itself is wrong.
      write(6,'(A)') 'BGCSCAL  lyr    t_soisno    t_scalar    w_scalar    o_scalar'&
                  // '  depth_scal      fpi_vr       hr_vr      smp_mm        wliq'&
                  // '        wice       porsl'
      DO j = 1, nl_soil
         write(6,'(A,I5,11E12.4)') 'BGCSCAL', j, &
            t_soisno(j,ipatch), t_scalar(j,ipatch), w_scalar(j,ipatch), &
            o_scalar(j,ipatch), depth_scalar(j,ipatch), fpi_vr(j,ipatch), &
            sum(decomp_hr_vr(j,:,ipatch)), smp(j,ipatch), &
            wliq_soisno(j,ipatch), wice_soisno(j,ipatch), porsl(j,ipatch)
      ENDDO

   END SUBROUTINE dump_decomp_scalars

   SUBROUTINE wetland_bgc_sasu_init (spinup_cycles)

! !DESCRIPTION:
!  Arm the wetland semi-analytic spin-up (DEF_METHANE%wetland_bgc_sasu) for
!  this run. It acts only when the run starts with two or more spin-up
!  cycles: the pools jump at the end of every cycle but the last, and the
!  last one runs plain, so the fast pools and mineral N settle on the jumped
!  state before history begins. A run without spin-up (stage 2 of the site
!  runs) sums and jumps nothing.

      IMPLICIT NONE
      integer, intent(in) :: spinup_cycles   ! spin-up cycles of this run, 0 without spin-up

      sasu_active = DEF_METHANE%wetland_bgc_sasu .and. spinup_cycles > 1
#ifdef WETLAND_PFT
      ! bgc_driver decomposes the wetland patch here, and CNSASU belongs to it.
      sasu_active = .false.
#endif
      IF (.not. (DEF_METHANE%wetland_bgc_sasu .and. p_is_master)) RETURN
      IF (sasu_active) THEN
         write(6,'(A,I4)') ' wetland_bgc_sasu: wetland pools jump at the end of spin-up cycles 1 to', &
            spinup_cycles - 1
         IF (.not. DEF_WETLAND_PEAT_C_PASSIVE) write(6,'(A)') &
            ' WARNING: wetland_bgc_sasu without DEF_WETLAND_PEAT_C_PASSIVE: the seeded peat in soil2', &
            ' and CWD jumps to the steady state of the plant input; only soil3 is held.'
      ELSE
         write(6,'(A)') ' wetland_bgc_sasu: this run makes no spin-up cycle to end, no jump'
      ENDIF

   END SUBROUTINE wetland_bgc_sasu_init

   SUBROUTINE wetland_bgc_sasu_accumulate (ipatch, deltim)

! !DESCRIPTION:
!  Sum this step's pools, rates and plant input for the jump at the end of
!  the spin-up cycle. Called after the decomposition fluxes are set and
!  before tracer_ch4_bgc_finalize_step updates the pools, so the pools are
!  the ones the fluxes came from and decomp_cpools_sourcesink still holds the
!  plant input alone. The rate is decomp_k, i.e. before fpi (see
!  wetland_bgc_sasu_jump); the realized litter loss is summed only to report
!  the cycle's mean fpi.

      IMPLICIT NONE
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      integer  :: j, k
      real(r8) :: c

      IF (.not. allocated(decomp_cpools_sourcesink) .or. &
          .not. allocated(decomp_npools_sourcesink)) RETURN

      IF (.not. allocated(sasu_time)) THEN
         allocate (sasu_time    (numpatch))                        ; sasu_time    (:)     = 0._r8
         allocate (sasu_c_acc   (nl_soil,ndecomp_pools,numpatch)) ; sasu_c_acc   (:,:,:) = 0._r8
         allocate (sasu_kc_acc  (nl_soil,ndecomp_pools,numpatch)) ; sasu_kc_acc  (:,:,:) = 0._r8
         allocate (sasu_k_acc   (nl_soil,ndecomp_pools,numpatch)) ; sasu_k_acc   (:,:,:) = 0._r8
         allocate (sasu_ic_acc  (nl_soil,ndecomp_pools,numpatch)) ; sasu_ic_acc  (:,:,:) = 0._r8
         allocate (sasu_in_acc  (nl_soil,ndecomp_pools,numpatch)) ; sasu_in_acc  (:,:,:) = 0._r8
         allocate (sasu_lit_pot (numpatch))                        ; sasu_lit_pot (:)     = 0._r8
         allocate (sasu_lit_act (numpatch))                        ; sasu_lit_act (:)     = 0._r8
      ENDIF

      sasu_time(ipatch) = sasu_time(ipatch) + deltim
      DO k = 1, ndecomp_pools
         DO j = 1, nl_soil
            c = max(decomp_cpools_vr(j,k,ipatch), 0._r8)
            sasu_c_acc (j,k,ipatch) = sasu_c_acc (j,k,ipatch) + c * deltim
            sasu_kc_acc(j,k,ipatch) = sasu_kc_acc(j,k,ipatch) + decomp_k(j,k,ipatch) * c * deltim
            sasu_k_acc (j,k,ipatch) = sasu_k_acc (j,k,ipatch) + decomp_k(j,k,ipatch) * deltim
            sasu_ic_acc(j,k,ipatch) = sasu_ic_acc(j,k,ipatch) + decomp_cpools_sourcesink(j,k,ipatch)
            sasu_in_acc(j,k,ipatch) = sasu_in_acc(j,k,ipatch) + decomp_npools_sourcesink(j,k,ipatch)
            IF (is_litter(k)) sasu_lit_pot(ipatch) = sasu_lit_pot(ipatch) &
               + decomp_k(j,k,ipatch) * c * dz_soi(j) * deltim
         ENDDO
      ENDDO
      DO k = 1, ndecomp_transitions
         IF (.not. is_litter(donor_pool(k))) CYCLE
         DO j = 1, nl_soil
            sasu_lit_act(ipatch) = sasu_lit_act(ipatch) + p_decomp_cpool_loss(j,k,ipatch) * dz_soi(j) * deltim
         ENDDO
      ENDDO

   END SUBROUTINE wetland_bgc_sasu_accumulate

   SUBROUTINE wetland_bgc_sasu_cycle_end (last_jump)

! !DESCRIPTION:
!  End of a spin-up cycle that another one follows: every wetland patch's
!  pools jump to the steady state of the cycle just run, and the sums are
!  cleared. With last_jump the next cycle is the last of the spin-up and
!  runs plain, so summing stops here. Called on every worker at the same
!  step; the summary line is reduced over the workers.

#ifdef USEMPI
      USE MOD_SPMD_Task, only: p_comm_worker, p_iam_worker, p_err, MPI_IN_PLACE, MPI_REAL8, MPI_SUM
#else
      USE MOD_SPMD_Task, only: p_iam_worker
#endif
      IMPLICIT NONE
      logical, intent(in) :: last_jump
      integer  :: i, nskip
      real(r8) :: before(4), after(4), fpi_mean
      ! patches jumped, layers skipped, column litter/CWD/soil1/soil2 before, after
      real(r8) :: tot(10)

      IF (.not. sasu_active) RETURN

      tot(:) = 0._r8
      IF (allocated(sasu_time)) THEN
         DO i = 1, numpatch
            IF (patchtype(i) /= 2 .or. sasu_time(i) <= 0._r8) CYCLE
            CALL wetland_bgc_sasu_jump (i, nskip, before, after)
            fpi_mean = 1._r8
            IF (sasu_lit_pot(i) > 0._r8) fpi_mean = sasu_lit_act(i) / sasu_lit_pot(i)
            tot(1)    = tot(1) + 1._r8
            tot(2)    = tot(2) + real(nskip, r8)
            tot(3:6)  = tot(3:6)  + before
            tot(7:10) = tot(7:10) + after
#ifdef SinglePoint
            write(6,'(A,I6,A,4ES11.3,A,4ES11.3,A,F6.3,A,I3)') ' WETSASU patch', i, &
               ' litter/CWD/soil1/soil2 [gC m-2]', before, ' ->', after, &
               '  cycle fpi', fpi_mean, '  layers kept', nskip
#endif
         ENDDO
      ENDIF

#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, tot, size(tot), MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
#endif
      IF (p_iam_worker == 0 .and. tot(1) > 0._r8) THEN
         write(6,'(A,I8,A,I6,A,4ES11.3,A,4ES11.3)') ' WETSASU jump: patches', nint(tot(1)), &
            '  layers kept', nint(tot(2)), '  mean litter/CWD/soil1/soil2 [gC m-2]', &
            tot(3:6) / tot(1), ' ->', tot(7:10) / tot(1)
      ENDIF

      IF (last_jump) THEN
         sasu_active = .false.
         IF (allocated(sasu_time)) THEN
            deallocate (sasu_time, sasu_c_acc, sasu_kc_acc, sasu_k_acc, sasu_ic_acc, sasu_in_acc, &
               sasu_lit_pot, sasu_lit_act)
         ENDIF
      ELSEIF (allocated(sasu_time)) THEN
         sasu_time   (:)     = 0._r8
         sasu_c_acc  (:,:,:) = 0._r8
         sasu_kc_acc (:,:,:) = 0._r8
         sasu_k_acc  (:,:,:) = 0._r8
         sasu_ic_acc (:,:,:) = 0._r8
         sasu_in_acc (:,:,:) = 0._r8
         sasu_lit_pot(:)     = 0._r8
         sasu_lit_act(:)     = 0._r8
      ENDIF

   END SUBROUTINE wetland_bgc_sasu_cycle_end

   SUBROUTINE wetland_bgc_sasu_jump (ipatch, nskip, before, after)

! !DESCRIPTION:
!  Steady state of one wetland patch's pools for the spin-up cycle just run,
!  layer by layer: the wetland tile has no vertical transport (bgc_driver,
!  which calls SoilBiogeochemLittVertTransp, runs for patchtype 0 only).
!
!  In a layer the pools follow dC/dt = I(t) + A K(t) C, with A the CENTURY
!  cascade (pathfrac_decomp, 1 - rf_decomp) and K the diagonal of rates.
!  Summed over a cycle, the cycle-mean pools Cm satisfy 0 = Im + A Ke Cm
!  exactly once the cycle repeats itself, Ke being the flux-weighted rate
!  sum(k C dt) / sum(C dt). The pools are set to Cm = -(A Ke)^-1 Im, the
!  steady state for this cycle's rates; the rates shift little from cycle to
!  cycle and the next jump corrects them, as in CoLM's CNSASU (Lu et al.
!  2020), which divides the same summed fluxes by the pools at the start of
!  the year instead.
!
!  Rates are taken before fpi. This tile's organic N leaves only as net
!  mineralization, and there is no N loss downstream, so at steady state net
!  mineralization equals the litter N input and mineral N builds up: fpi is
!  1 there. A cycle that starts with the fresh pools empty has fpi far
!  below 1 at the productive towers (0.04-0.2), because the soil pools whose
!  mineralization would feed the litter are not built yet. Jumping with those
!  rates lifts litter by 1/fpi and keeps the tile N-locked for decades.
!
!  soil3 is held: it carries the seeded old peat and the buried litter
!  (wetland_burial_frac_site), which are not in equilibrium with today's
!  input and turn over in millennia; its own input is left out. Its outflow to
!  soil1 enters as an input, and the flows into it are losses. Floating-C:N
!  pools get their N from the same system on N (litter N leaves with its C,
!  so litter settles at the C:N of its input); fixed pools take C over C:N.

      IMPLICIT NONE
      integer,  intent(in)  :: ipatch
      integer,  intent(out) :: nskip                   ! layers left alone
      real(r8), intent(out) :: before(4), after(4)     ! column litter, CWD, soil1, soil2 C [gC m-2]

      integer  :: j, k, l, m, d, r, nj, nf
      integer  :: pj(ndecomp_pools), pf(ndecomp_pools), pos(ndecomp_pools), posf(ndecomp_pools)
      real(r8) :: tsum, fr
      real(r8) :: keff(ndecomp_pools)
      real(r8) :: a (ndecomp_pools,ndecomp_pools), b (ndecomp_pools)
      real(r8) :: an(ndecomp_pools,ndecomp_pools), bn(ndecomp_pools)
      logical  :: ok

      ! Pools that jump (all but soil3), and the floating-C:N ones among them.
      nj = 0; nf = 0; pos(:) = 0; posf(:) = 0
      DO l = 1, ndecomp_pools
         IF (l == i_soil3) CYCLE
         nj = nj + 1; pj(nj) = l; pos(l) = nj
         IF (floating_cn_ratio(l)) THEN
            nf = nf + 1; pf(nf) = l; posf(l) = nf
         ENDIF
      ENDDO

      tsum  = sasu_time(ipatch)
      nskip = 0
      CALL sasu_column_c (ipatch, before)

      DO j = 1, nl_soil

         ! A pool that stayed empty all cycle has no flux-weighted rate; its
         ! time-mean rate serves, and with no input it stays empty.
         DO l = 1, ndecomp_pools
            IF (sasu_c_acc(j,l,ipatch) > 1.e-12_r8 * tsum) THEN
               keff(l) = sasu_kc_acc(j,l,ipatch) / sasu_c_acc(j,l,ipatch)
            ELSE
               keff(l) = sasu_k_acc(j,l,ipatch) / tsum
            ENDIF
         ENDDO

         a (:,:) = 0._r8; b (:) = 0._r8
         an(:,:) = 0._r8; bn(:) = 0._r8
         DO m = 1, nj
            a(m,m) = - keff(pj(m))
            b(m)   = - sasu_ic_acc(j,pj(m),ipatch) / tsum
         ENDDO
         DO m = 1, nf
            an(m,m) = - keff(pf(m))
            bn(m)   = - sasu_in_acc(j,pf(m),ipatch) / tsum
         ENDDO
         DO k = 1, ndecomp_transitions
            d  = donor_pool(k)
            r  = receiver_pool(k)
            IF (r < 1) CYCLE                  ! respired whole: on the diagonal already
            IF (pos(r) == 0) CYCLE            ! into the held soil3: a loss, on the diagonal
            fr = pathfrac_decomp(j,k,ipatch)
            IF (pos(d) > 0) THEN
               a(pos(r),pos(d)) = a(pos(r),pos(d)) + (1._r8 - rf_decomp(j,k,ipatch)) * fr * keff(d)
               ! N moves with the donor's C:N, none of it respired (ntransfer)
               IF (posf(d) > 0 .and. posf(r) > 0) &
                  an(posf(r),posf(d)) = an(posf(r),posf(d)) + fr * keff(d)
            ELSE                               ! out of the held soil3: an input
               b(pos(r)) = b(pos(r)) - (1._r8 - rf_decomp(j,k,ipatch)) * fr &
                  * sasu_kc_acc(j,d,ipatch) / tsum
            ENDIF
         ENDDO

         CALL sasu_solve (nj, a(1:nj,1:nj), b(1:nj), ok)
         IF (ok .and. nf > 0) CALL sasu_solve (nf, an(1:nf,1:nf), bn(1:nf), ok)
         IF (ok) ok = all(ieee_is_finite(b(1:nj))) .and. all(ieee_is_finite(bn(1:nf))) &
            .and. all(b(1:nj) < SASU_C_MAX)
         IF (.not. ok) THEN
            nskip = nskip + 1
            CYCLE
         ENDIF

         DO m = 1, nj
            l = pj(m)
            decomp_cpools_vr(j,l,ipatch) = max(b(m), 0._r8)
            IF (floating_cn_ratio(l)) THEN
               decomp_npools_vr(j,l,ipatch) = max(bn(posf(l)), 0._r8)
            ELSE
               decomp_npools_vr(j,l,ipatch) = decomp_cpools_vr(j,l,ipatch) / initial_cn_ratio(l)
            ENDIF
         ENDDO
      ENDDO

      CALL sasu_column_c (ipatch, after)

   END SUBROUTINE wetland_bgc_sasu_jump

   SUBROUTINE sasu_column_c (ipatch, col)

      ! Column litter, CWD, soil1 and soil2 carbon [gC m-2], for the log.
      IMPLICIT NONE
      integer,  intent(in)  :: ipatch
      real(r8), intent(out) :: col(4)

      col(1) = sum((decomp_cpools_vr(1:nl_soil,i_met_lit,ipatch) + decomp_cpools_vr(1:nl_soil,i_cel_lit,ipatch) &
         + decomp_cpools_vr(1:nl_soil,i_lig_lit,ipatch)) * dz_soi(1:nl_soil))
      col(2) = sum(decomp_cpools_vr(1:nl_soil,i_cwd  ,ipatch) * dz_soi(1:nl_soil))
      col(3) = sum(decomp_cpools_vr(1:nl_soil,i_soil1,ipatch) * dz_soi(1:nl_soil))
      col(4) = sum(decomp_cpools_vr(1:nl_soil,i_soil2,ipatch) * dz_soi(1:nl_soil))

   END SUBROUTINE sasu_column_c

   SUBROUTINE sasu_solve (n, a, b, ok)

      ! Gaussian elimination with partial pivoting for a x = b; x replaces b.
      ! n is at most ndecomp_pools, so no library call is worth it.
      IMPLICIT NONE
      integer,  intent(in)    :: n
      real(r8), intent(inout) :: a(n,n), b(n)
      logical,  intent(out)   :: ok
      integer  :: i, k, p
      real(r8) :: f, row(n)

      ok = .false.
      DO k = 1, n
         p = k - 1 + maxloc(abs(a(k:n,k)), dim=1)
         IF (.not. (abs(a(p,k)) > 0._r8)) RETURN
         IF (p /= k) THEN
            row(:) = a(k,:); a(k,:) = a(p,:); a(p,:) = row(:)
            f = b(k); b(k) = b(p); b(p) = f
         ENDIF
         DO i = k + 1, n
            f = a(i,k) / a(k,k)
            a(i,k:n) = a(i,k:n) - f * a(k,k:n)
            b(i) = b(i) - f * b(k)
         ENDDO
      ENDDO
      DO k = n, 1, -1
         b(k) = (b(k) - sum(a(k,k+1:n) * b(k+1:n))) / a(k,k)
      ENDDO
      ok = .true.

   END SUBROUTINE sasu_solve

END MODULE MOD_Tracer_Reactive_BgcShim
#endif
