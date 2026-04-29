! ============================================================================
! MODULE_GLACIER
!
! Purpose:
!   Land-ice (glacier and ice sheet) thermal and hydrological processes for
!   the CSSPv2-GLC land surface model. Provides coupled energy-balance and
!   mass-balance calculations.
!
! Original author: Shuwen Li and Xing Yuan
! Modified by: Shuwen Li (swli), 2023 - 2026
!     * Integrated a physically based glacier energy-balance module into
!       the CSSPv2 land surface model (CSSPv2-GLC); see Li and Yuan (2026)
!       for the full model description and evaluation
!     * Rewrote groundfluxes_glacier with Andreas (1987) scalar roughness
!       and Smeets and Van den Broeke (2008) stability-corrected fluxes
!     * Enhanced groundtem_glacier with shortwave radiation penetration
!       (Beer-Lambert law) and a constant geothermal heat-flux bottom
!       boundary condition
!     * Added meltf_refreeze_glc for refreezing-aware phase change with
!       explicit separation of melt and refreeze mass accounting
!     * Added depth-weighted snow/ice albedo blending in GLC_netsolar

! References:
!   [1] Yuan, X., C. Li, P. Ji, X.Z. Xi, S.W. Li, Y. Jiao, H. Yi, X.Z. Liang,
!       and Q.C. Zeng (2026): Development and evaluation of Conjunctive
!       Surface-Subsurface Process Model version 3 (CSSPv3) at global scale.
!       Journal of Advances in Modeling Earth Systems. (Revised)
!   [2] Li, S., and X. Yuan (2026): Evaluation of a Land Surface-Glacier
!       Coupled Model over the Three-River Headwaters Region. Water, 18, 1030.
! ============================================================================
module MODULE_GLACIER
#include "configs.inc"
      USE shr_kind_mod, ONLY: r8 => shr_kind_rn
      USE module_sf_clm_constant
      IMPLICIT NONE

      ! PUBLIC MEMBER FUNCTIONS:
      PUBLIC :: glacier_temp
      PUBLIC :: glacier_water
      PUBLIC :: rain_snow_temp
      PUBLIC :: GLC_netsolar

      ! PRIVATE MEMBER FUNCTIONS:
      private :: groundfluxes_glacier
      private :: glacier_surface_temperature
      private :: groundtem_glacier
      PRIVATE :: combo_glc
      private :: winddriftcompaction_glc
      private :: snowcompaction_glc
      private :: snowlayerscombine_glc
      private :: snowlayersdivide_glc
      private :: meltf_glc
      private :: qsadv_glc
      private :: moninobuk_glc
      private :: moninobukini_glc
      private :: snowwater
      private :: tridia
      private :: calc_dewpoint
      private :: meltf_refreeze_glc

      contains


      ! ======================================================================
      ! Purpose: Compute net solar radiation absorbed by the glacier surface,
      ! accounting for snow albedo, ice albedo, and snow cover fraction.
      !
      ! Authors: Shuwen Li and Xing Yuan
      ! Modified by: Shuwen Li (swli), 2023 - 2026
      !   * Added depth-weighted albedo blending between overlying snow and
      !     underlying glacier ice using an exponential attenuation with
      !     snow depth (threshold snowdph = 0.057 m)
      !   * Attempted dewpoint-based albedo parameterization for glacier ice
      !
      ! References:
      ! Li, S., and X. Yuan (2026): Evaluation of a Land Surface-Glacier
!     ! Coupled Model over the Three-River Headwaters Region. Water, 18, 1030.
      ! ======================================================================
      SUBROUTINE GLC_netsolar(tm, psrf, qm, cosza, scv, snowdp, fsno, sag, &
                              sabg, solis, alb_glacier, ipr, jpr)
      USE module_sf_clm_constant, ONLY: tfrz, vegmin, &
                                        alblis, alblil, albfls, albfll, albsns, albsnl
      IMPLICIT NONE

      ! Input variables
      REAL(r8), INTENT(IN)  :: snowdp        ! Snow depth [m]
      REAL(r8), INTENT(IN)  :: fsno          ! Snow cover fraction [-]
      REAL(r8), INTENT(IN)  :: sag           ! Snow age factor
      REAL(r8), INTENT(IN)  :: scv           ! Snow water equivalent [kg/m2]
      REAL(r8), INTENT(IN)  :: solis(2,2)    ! Downward solar radiation [W/m2]
      REAL(r8), INTENT(IN)  :: tm            ! Surface temperature [K]
      REAL(r8), INTENT(IN)  :: psrf          ! Surface pressure [Pa]
      REAL(r8), INTENT(IN)  :: qm            ! Specific humidity [kg/kg]
      REAL(r8), INTENT(IN)  :: cosza         ! Cosine of solar zenith angle [-]
      integer,  INTENT(IN)  :: ipr, jpr

      ! Output variables
      REAL(r8), INTENT(OUT) :: sabg          ! Solar absorption by glacier surface [W/m2]
      REAL(r8)              :: alb_glacier   ! Glacier/ice albedo [-]

      ! Local variables
      REAL(r8) :: alb_snow      ! Average snow albedo
      REAL(r8) :: alb_ice       ! Final surface albedo
      REAL(r8) :: solis_sum     ! Total downward solar radiation [W/m2]
      REAL(r8) :: albis         ! Depth-weighted blended albedo
      REAL(r8) :: albsno_s(2,2) ! Snow albedo by band/type
      REAL(r8) :: czen, czf, age, cff
      REAL(r8) :: dfalbl, dfalbs, dralbl, dralbs

      real(r8), save :: &
            sl   = 2.0,  &  ! Zenith angle albedo dependence factor
            cons = 0.2,  &  ! NIR snow albedo constant
            conn = 0.5      ! Visible snow albedo constant

      REAL(r8), PARAMETER :: snowdph = 0.057  ! Snow depth threshold for albedo blending [m]


      albsno_s(:,:) = 0.0_r8
      sabg = 0.0_r8

      if (cosza <= 0.) return

      czen = max(cosza, 0.001)
      if (scv > 0.) then
            ! Correction for snow age
            age    = 1. - 1./(1. + sag)
            dfalbs = albsns*(1. - cons*age)
            ! Zenith angle correction for new snow
            cff    = max(0., (1. + 1./sl)/(1. + czen*2.*sl) - 1./sl)
            czf    = 0.4*cff*(1. - dfalbs)
            dralbs = dfalbs + czf
            dfalbl = albsnl*(1. - conn*age)
            czf    = 0.4*cff*(1. - dfalbl)
            dralbl = dfalbl + czf

            albsno_s(1,1) = dralbs
            albsno_s(2,1) = dralbl
            albsno_s(1,2) = dfalbs
            albsno_s(2,2) = dfalbl
      endif

      ! Compute average snow albedo
      alb_snow = (albsno_s(1,1) + albsno_s(1,2) + albsno_s(2,1) + albsno_s(2,2)) / 4.0

      ! Depth-weighted blended albedo: thin snow allows underlying ice to show through
      albis = alb_snow + (alb_glacier - alb_snow) * EXP(-snowdp / snowdph)

      ! Final surface albedo based on snow cover fraction
      IF (fsno > 0.0) THEN
            alb_ice = albis * fsno + alb_glacier * (1.0 - fsno)
      ELSE
            alb_ice = alb_glacier
      END IF
      alb_ice = min(0.90_r8, max(0.05_r8, alb_ice))

      solis_sum = SUM(solis(:,:))
      sabg = solis_sum * (1.0 - alb_ice)

      END SUBROUTINE GLC_netsolar


      ! ======================================================================
      ! GLACIER_TEMP
      ! Purpose: Main subroutine for thermal processes and surface fluxes of land ice.
      !
      ! Authors: Shuwen Li and Xing Yuan
      ! Modified by Shuwen Li (swli), 2023 - 2026:
      !   * Reorganized surface flux correction and energy-balance error
      !     accounting for consistency with updated phase-change subroutine
      !   * Added precipitation advected heat to the energy balance
      !
      ! References:
      !   [1] Yuan, X., C. Li, P. Ji, X.Z. Xi, S.W. Li, Y. Jiao, H. Yi, X.Z. Liang,
      !       and Q.C. Zeng (2026): Development and evaluation of Conjunctive
      !       Surface-Subsurface Process Model version 3 (CSSPv3) at global scale.
      !       Journal of Advances in Modeling Earth Systems. (Revised)
      !   [2] Li, S., and X. Yuan (2026): Evaluation of a Land Surface-Glacier
      !       Coupled Model over the Three-River Headwaters Region. Water, 18, 1030.
      ! ======================================================================
      subroutine GLACIER_TEMP (scwater, lb, nl_ice, dtime, &
                        hw, ht, hq, &
                        us, vs, tm, qm, &
                        rhoair, psrf, nl_snow, sabg, &
                        frl, fsno, dz, z, &
                        zi, tss, wice, wliq, &
                        scv, snowdp, imelt, taux, &
                        tauy, fsena, fevpa, lfevpa, &
                        fseng, fevpg, olrg, fgrnd, &
                        qseva, qsdew, qsubl, qfros, &
                        sm, tref, qref, &
                        errore, emis, z0ma, zol, &
                        rib, ustar, qstar, tstar, &
                        fm, fh, fq, pg_rain, &
                        pg_snow, t_precip, sabvg, xmf, melt, smelt, refreeze, &
                        msub, &
                        ipr, jpr)

      USE shr_kind_mod, ONLY: r8 => shr_kind_rn
      use module_sf_clm_constant, only : hvap, hsub, rair, cpair, stefnc, tfrz, cpliq, cpice, cnfac, capr, denliq

      IMPLICIT NONE

      integer, INTENT(in) :: &
            scwater,    &  ! Land patch type (3 = land ice)
            lb,         &  ! Lower array bound
            nl_ice,     &  ! Upper array bound (ice layers)
            nl_snow
      integer ipr, jpr

      real(r8), INTENT(in) :: &
            dtime,      &  ! Model time step [s]
            hw,         &  ! Wind observation height [m]
            ht,         &  ! Temperature observation height [m]
            hq,         &  ! Humidity observation height [m]
            us,         &  ! Eastward wind component [m/s]
            vs,         &  ! Northward wind component [m/s]
            tm,         &  ! Air temperature at reference height [K]
            qm,         &  ! Specific humidity at reference height [kg/kg]
            rhoair,     &  ! Air density [kg/m3]
            psrf,       &  ! Surface pressure [Pa]
            t_precip,   &  ! Precipitation temperature [K]
            pg_rain,    &  ! Rainfall rate [kg/(m2 s)]
            pg_snow,    &  ! Snowfall rate [kg/(m2 s)]
            sabg,       &  ! Solar radiation absorbed by ground [W/m2]
            frl,        &  ! Atmospheric longwave radiation [W/m2]
            fsno,       &  ! Snow cover fraction [-]
            dz(lb:nl_ice), &  ! Layer thickness [m]
            z(lb:nl_ice),  &  ! Node depth [m]
            zi(lb-1:nl_ice)   ! Interface depth [m]

      REAL(r8), intent(in) :: sabvg  ! Solar absorbed by ground + vegetation [W/m2]

      real(r8), INTENT(inout) :: &
            tss(lb:nl_ice),  &  ! Snow/ice temperature [K]
            wice(lb:nl_ice), &  ! Ice content [kg/m2]
            wliq(lb:nl_ice), &  ! Liquid water content [kg/m2]
            scv,             &  ! Snow water equivalent [kg/m2]
            snowdp              ! Snow depth [m]

      integer, INTENT(out) :: &
            imelt(nl_snow+1:nl_ice)  ! Phase change flag: 1=melt, 2=freeze, 0=none

      real(r8), INTENT(out) :: &
            taux,    &  ! East-west wind stress [kg/m/s2]
            tauy,    &  ! North-south wind stress [kg/m/s2]
            fsena,   &  ! Sensible heat flux to atmosphere [W/m2]
            lfevpa,  &  ! Latent heat flux to atmosphere [W/m2]
            fseng,   &  ! Sensible heat flux from ground [W/m2]
            fevpg,   &  ! Evaporation flux from ground [mm/s]
            olrg,    &  ! Outgoing longwave radiation [W/m2]
            fgrnd,   &  ! Ground heat flux [W/m2]
            fevpa,   &  ! Evapotranspiration to atmosphere [mm H2O/s]
            qseva,   &  ! Ground surface evaporation rate [mm H2O/s]
            qsdew,   &  ! Ground surface dew formation [mm H2O/s]
            qsubl,   &  ! Sublimation from snow pack [mm H2O/s]
            qfros,   &  ! Surface dew added to snow pack [mm H2O/s]
            sm,      &  ! Snowmelt rate [kg/(m2 s)]
            tref,    &  ! 2m air temperature [K]
            qref,    &  ! 2m specific humidity [kg/kg]
            emis,    &  ! Bulk surface emissivity [-]
            z0ma,    &  ! Effective roughness length [m]
            zol,     &  ! Dimensionless height (z/L) in Monin-Obukhov theory
            rib,     &  ! Bulk Richardson number
            ustar,   &  ! Friction velocity [m/s]
            qstar,   &  ! Humidity scale parameter [kg/kg]
            tstar,   &  ! Temperature scale parameter [K]
            fm,      &  ! Momentum profile function integral
            fh,      &  ! Heat profile function integral
            fq          ! Moisture profile function integral

      ! Local variables
      integer j

      real(r8) :: &
            cgrnd,        &  ! d(ice energy flux)/d(ice temp) [W/m2/K]
            cgrndl,       &  ! d(ice sensible heat)/d(ice temp) [W/m2/K]
            cgrnds,       &  ! d(ice latent heat)/d(ice temp) [W/m2/K]
            degdT,        &  ! d(eg)/dT
            dqgdT,        &  ! d(qg)/dT
            eg,           &  ! Vapor pressure at temperature T [Pa]
            egsmax,       &  ! Max evaporation available in one time step
            egidif,       &  ! Excess evaporation over egsmax
            emg,          &  ! Ground emissivity (0.97)
            errore,       &  ! Energy balance error [W/m2]
            fact(lb:nl_ice), & ! Tridiagonal matrix factor
            htvp,         &  ! Latent heat of vaporization or sublimation [J/kg]
            qg,           &  ! Ground specific humidity [kg/kg]
            qsatg,        &  ! Saturation humidity [kg/kg]
            qsatgdT,      &  ! d(qsatg)/dT
            qred,         &  ! Ice surface relative humidity
            thm,          &  ! Intermediate variable (tm + 0.0098*ht)
            th,           &  ! Potential temperature [K]
            thv,          &  ! Virtual potential temperature [K]
            tg,           &  ! Ground surface temperature [K]
            tss_bef(lb:nl_ice), & ! Snow/ice temperature before update
            tinc,         &  ! Temperature increment between time steps
            ur,           &  ! Wind speed at reference height [m/s]
            xmf,          &  ! Total latent heat of phase change [W/m2]
            melt,         &  ! Total melt mass [kg/m2]
            smelt,        &  ! Surface melt mass [kg/m2]
            refreeze,     &  ! Total refreeze mass [kg/m2]
            msub

      ! ======================================================================
      ! [1] Initialization
      ! ======================================================================
      taux   = 0.;  tauy   = 0.
      fsena  = 0.;  fevpa  = 0.
      lfevpa = 0.;  tref   = tm
      fseng  = 0.;  fevpg  = 0.
      qref   = qm
      emis   = 0.;  z0ma   = 0.
      zol    = 0.;  rib    = 0.
      ustar  = 0.;  qstar  = 0.;  tstar  = 0.
      msub   = 0.;  fgrnd  = 0.
      errore = 0.

      tg = tss(lb)
      tss_bef(lb:nl_ice) = tss(lb:nl_ice)

      emg = 0.97

      ! Latent heat: sublimation if below freezing, evaporation otherwise
      if (tg >= tfrz) then
            htvp = hvap
      else
            htvp = hsub
      endif

      thm = tm + 0.0098*ht
      th  = tm*(100000./psrf)**(rair/cpair)
      thv = th*(1. + 0.61*qm)
      ur  = max(0.1, sqrt(us*us + vs*vs))

      ! ======================================================================
      ! [2] Specific humidity and derivative at ground surface
      ! ======================================================================
      qred = 1.
      call qsadv_glc(tm, psrf, eg, degdT, qsatg, qsatgdT)
      qg    = qred*qsatg
      dqgdT = qred*qsatgdT

      ! ======================================================================
      ! [3] Sensible and latent fluxes and their derivatives
      ! ======================================================================
      call groundfluxes_glacier(hw, ht, hq, us, vs, tm, qm, rhoair, ur, thm, th, thv, tg, qg, dqgdT, htvp, &
                                fsno, cgrnd, cgrndl, cgrnds, taux, tauy, fsena, fevpa, fseng, fevpg, tref, qref, &
                                z0ma, zol, rib, ustar, qstar, tstar, fm, fh, fq, ipr, jpr)

      ! ======================================================================
      ! [4] Update surface temperature via thermal diffusion solver
      ! ======================================================================
      call groundtem_glacier(scwater, lb, nl_ice, dtime, tm, ht, capr, cnfac, dz, z, zi, nl_snow, tss, tss_bef, wice, wliq, &
                             scv, snowdp, frl, sabg, sabvg, fseng, fevpg, cgrnd, cgrnds, cgrndl, htvp, emg, tg, &
                             imelt, sm, xmf, melt, smelt, refreeze, fact, pg_rain, pg_snow, t_precip, ipr, jpr)

      ! ======================================================================
      ! [5] Correct fluxes to updated ice temperature
      ! ======================================================================
      tg   = tss(lb)
      tinc = tss(lb) - tss_bef(lb)
      fseng = fseng + tinc*cgrnds
      fevpg = fevpg + tinc*cgrndl

      ! Limit evaporation to available water
      egsmax = (wice(lb) + wliq(lb)) / dtime
      egidif = max(0., fevpg - egsmax)
      fevpg  = min(fevpg, egsmax)
      fseng  = fseng + htvp*egidif

      ! Total fluxes to atmosphere
      fsena  = fseng
      fevpa  = fevpg
      lfevpa = htvp*fevpg

      qseva = 0.;  qsubl = 0.;  qfros = 0.;  qsdew = 0.

      if (tg < tfrz) then
            qsubl = max(fevpg, 0.)
            qsdew = min(fevpg, 0.)
      else
            qseva = max(fevpg, 0.)
            qfros = min(fevpg, 0.)
      endif

      msub = qfros + qsdew + qseva + qsubl

      ! Ground heat flux
      fgrnd = sabg + emg*frl &
            - emg*stefnc*tss_bef(lb)**3*(tss_bef(lb) + 4.*tinc) &
            - (fseng + fevpg*htvp) &
            + cpliq * pg_rain * (t_precip - tss(lb)) &
            + cpice * pg_snow * (t_precip - tss(lb))

      ! Outgoing longwave radiation
      olrg = (1. - emg)*frl + emg*stefnc * tss_bef(lb)**4 &
           + 4.*emg*stefnc*tss_bef(lb)**3*tinc

      emis = emg

      ! ======================================================================
      ! [6] Energy balance error
      ! NOTE: the +0.06 term is the geothermal heat flux [W/m2] applied as a
      ! hardcoded constant, following Davies and Davies (2010) [Ref. 8].
      ! This should match the bottom_heat_flux in groundtem_glacier for
      ! physical consistency.
      ! ======================================================================
      errore = (sabg + frl - olrg) - (fsena + lfevpa) - xmf + 0.06 &
             + (cpliq * pg_rain * (t_precip - tfrz) + cpice * pg_snow * (t_precip - tfrz))

      do j = lb, nl_ice
            errore = errore - ((wice(j)*cpice + wliq(j)*cpliq) * (tss(j) - tss_bef(j)) / dtime)
      enddo

      end subroutine GLACIER_TEMP


      ! ======================================================================
      ! groundfluxes_glacier
      ! Purpose: Compute surface turbulent heat and momentum fluxes over glacier ice.
      ! Authors: Shuwen Li and Xing Yuan
      !
      ! References:
      ! Li, S., and X. Yuan (2026): Evaluation of a Land Surface-Glacier
!     ! Coupled Model over the Three-River Headwaters Region. Water, 18, 1030.
      ! ======================================================================
      subroutine groundfluxes_glacier(hu, ht, hq, &
                                      us, vs, tm, qm, rhoair, &
                                      ur, thm, th, thv, tg, qg, dqgdT, htvp, &
                                      fsno, cgrnd, cgrndl, cgrnds, &
                                      taux, tauy, fsena, fevpa, fseng, fevpg, tref, qref, &
                                      z0ma, zol, rib, ustar, qstar, tstar, fm, fh, fq, ipr, jpr)

      USE shr_kind_mod, ONLY: r8 => shr_kind_rn
      use module_sf_clm_constant, only : cpair, vonkar, grav, z0sno, z0ice, hsub
      implicit none

      integer, intent(in) :: ipr, jpr

      real(r8), INTENT(in) :: &
            hu,      &  ! Wind observation height [m]
            ht,      &  ! Temperature observation height [m]
            hq,      &  ! Humidity observation height [m]
            us,      &  ! Eastward wind component [m/s]
            vs,      &  ! Northward wind component [m/s]
            tm,      &  ! Air temperature at reference height [K]
            qm,      &  ! Specific humidity at reference height [kg/kg]
            rhoair,  &  ! Air density [kg/m3]
            fsno,    &  ! Snow cover fraction [-]
            ur,      &  ! Wind speed at reference height [m/s]
            thm,     &  ! Intermediate variable (tm + 0.0098*ht)
            th,      &  ! Potential temperature [K]
            thv,     &  ! Virtual potential temperature [K]
            tg,      &  ! Ground surface temperature [K]
            qg,      &  ! Ground specific humidity [kg/kg]
            dqgdT,   &  ! d(qg)/dT
            htvp        ! Latent heat of vaporization or sublimation [J/kg]

      real(r8), INTENT(out) :: &
            taux,    &  ! East-west wind stress [kg/m/s2]
            tauy,    &  ! North-south wind stress [kg/m/s2]
            fsena,   &  ! Sensible heat flux to atmosphere [W/m2]
            fevpa,   &  ! Evaporation flux to atmosphere [mm/s]
            fseng,   &  ! Sensible heat flux from ground [W/m2]
            fevpg,   &  ! Evaporation flux from ground [mm/s]
            cgrnd,   &  ! d(ice energy flux)/d(ice temp) [W/m2/K]
            cgrndl,  &  ! d(ice sensible heat)/d(ice temp) [W/m2/K]
            cgrnds,  &  ! d(ice latent heat)/d(ice temp) [W/m2/K]
            tref,    &  ! 2m air temperature [K]
            qref,    &  ! 2m specific humidity [kg/kg]
            z0ma,    &  ! Effective roughness length [m]
            zol,     &  ! Dimensionless height z/L
            rib,     &  ! Bulk Richardson number
            ustar,   &  ! Friction velocity [m/s]
            tstar,   &  ! Temperature scale parameter [K]
            qstar,   &  ! Humidity scale parameter [kg/kg]
            fm,      &  ! Momentum profile function integral
            fh,      &  ! Heat profile function integral
            fq          ! Moisture profile function integral

      ! Local constants
      real(r8), parameter :: &
            nu_air         = 1.5e-5,  &  ! Kinematic viscosity of air [m2/s]
            Pr             = 0.71,    &  ! Prandtl number
            Sc             = 0.60,    &  ! Schmidt number
            CONVERGENCE_TOL = 0.1,   &  ! Convergence tolerance
            MIN_WIND       = 0.1,    &  ! Minimum wind speed [m/s]
            MAX_ZETA       = 10.0,   &  ! Maximum stability parameter
            MIN_ZETA       = -10.0,  &  ! Minimum stability parameter
            DAMPING_FACTOR = 0.7        ! Damping factor to prevent oscillations

      integer :: niters = 10, iter, nmozsgn, convergence_count

      real(r8) :: &
            beta,      &  ! Convective velocity coefficient
            displax,   &  ! Zero displacement height [m]
            dth,       &  ! Temperature difference: reference height - surface
            dqh,       &  ! Humidity difference: reference height - surface
            dthv,      &  ! Virtual potential temperature difference
            obu,       &  ! Monin-Obukhov length [m]
            obuold,    &  ! Monin-Obukhov length from previous iteration
            ram,       &  ! Aerodynamic resistance [s/m]
            rah,       &  ! Heat resistance [s/m]
            raw,       &  ! Moisture resistance [s/m]
            raih,      &  ! Intermediate variable [kg/m2/s]
            raiw,      &  ! Intermediate variable [kg/m2/s]
            fh2m,      &  ! Temperature profile relation at 2m
            fq2m,      &  ! Humidity profile relation at 2m
            fm10m,     &  ! Momentum profile function at 10m
            thvstar,   &  ! Virtual potential temperature scale
            um,        &  ! Wind speed including stability effects [m/s]
            zeta,      &  ! Dimensionless height z/L
            zeta_old,  &  ! Previous iteration zeta
            zii,       &  ! Convective boundary layer height [m]
            zldis,     &  ! Reference height minus zero displacement height [m]
            z0mg,      &  ! Momentum roughness length [m]
            z0hg,      &  ! Heat roughness length [m]
            z0qg,      &  ! Moisture roughness length [m]
            relative_change  ! Relative change in obu between iterations

      logical :: converged

      ! -----------------------------------------------------------------------
      ! Determine initial momentum roughness length
      ! -----------------------------------------------------------------------
      if (fsno > 0.5) then
            z0mg = z0sno   ! Snow-covered surface
      else
            z0mg = z0ice   ! Bare ice surface
      endif
      z0mg = max(z0mg, 1.0e-6)
      z0mg = min(z0mg, 0.1)

      ! Initialize variables
      beta    = 1.0
      zii     = 1000.0
      z0ma    = z0mg
      displax = 0.0

      dth   = thm - tg
      dqh   = qm  - qg
      dthv  = dth*(1.0 + 0.61*qm) + 0.61*th*dqh
      zldis = hu - displax

      nmozsgn         = 0
      convergence_count = 0
      converged       = .false.
      obuold          = 0.0

      call moninobukini_glc(ur, th, thm, thv, dth, dqh, dthv, zldis, z0mg, um, obu)

      ! Initialize zeta before iteration
      if (abs(obu) > 1.e-6) then
            zeta = zldis / obu
      else
            zeta = 0.0
      endif
      zeta = max(MIN_ZETA, min(MAX_ZETA, zeta))

      ! -----------------------------------------------------------------------
      ! Stability iteration
      ! -----------------------------------------------------------------------
      ITERATION : do iter = 1, niters

            zeta_old = zeta

            call calculate_scalar_roughness(z0mg, ustar, z0hg, z0qg)

            call moninobuk_glc(hu, ht, hq, displax, z0mg, z0hg, z0qg, obu, um, &
                               ustar, fh2m, fq2m, fm10m, fm, fh, fq)

            if (abs(fh) > 1.e-10) then
                  tstar = vonkar/fh*dth
            else
                  tstar = 0.0
            endif

            if (abs(fq) > 1.e-10) then
                  qstar = vonkar/fq*dqh
            else
                  qstar = 0.0
            endif

            thvstar = tstar*(1.0 + 0.61*qm) + 0.61*th*qstar

            if (abs(ustar) > 1.e-10 .and. abs(thv) > 1.e-10) then
                  zeta = zldis*vonkar*grav*thvstar/(ustar**2*thv)
            else
                  zeta = 0.0
            endif
            zeta = max(MIN_ZETA, min(MAX_ZETA, zeta))

            if (abs(zeta) > 1.e-10) then
                  obu = zldis/zeta
            else
                  obu = 1.e6
            endif

            call adjust_wind_speed(ur, ustar, thvstar, thv, zeta, zii, beta, um)

            if (iter > 1) then
                  relative_change = abs((obu - obuold) / max(abs(obu), 1.0))
                  if (relative_change < CONVERGENCE_TOL) then
                        convergence_count = convergence_count + 1
                        if (convergence_count >= 2) then
                              converged = .true.
                              EXIT
                        endif
                  else
                        convergence_count = 0
                  endif
            endif

            if (iter > 1 .and. obuold*obu < 0.0) then
                  nmozsgn = nmozsgn + 1
                  if (nmozsgn >= 4) EXIT
            endif

            if (iter > 2) then
                  obu = DAMPING_FACTOR * obu + (1.0 - DAMPING_FACTOR) * obuold
            endif

            obuold = obu

      enddo ITERATION

      ! -----------------------------------------------------------------------
      ! Compute final resistances and derivatives
      ! -----------------------------------------------------------------------
      if (abs(ustar) > 1.e-10) then
            ram = 1.0/(ustar*ustar/um)
            rah = 1.0/(vonkar/fh*ustar)
            raw = 1.0/(vonkar/fq*ustar)
      else
            ram = 1.e6
            rah = 1.e6
            raw = 1.e6
      endif

      raih = rhoair*cpair/rah
      raiw = rhoair/raw

      cgrnds = raih
      cgrndl = raiw*dqgdT
      cgrnd  = cgrnds + htvp*cgrndl

      zol = zeta
      if (abs(um) > 1.e-10 .and. abs(fh) > 1.e-10) then
            rib = min(5.0, zol*ustar**2/(vonkar**2/fh*um**2))
      else
            rib = 0.0
      endif

      ! Surface fluxes
      taux  = -rhoair*us/ram
      tauy  = -rhoair*vs/ram
      fseng = -raih*dth
      fevpg = -raiw*dqh
      fsena = fseng
      fevpa = fevpg

      ! 2m temperature and humidity
      if (abs(fh) > 1.e-10 .and. abs(fq) > 1.e-10) then
            tref = thm + vonkar/fh*dth * (fh2m/vonkar - fh/vonkar)
            qref = qm  + vonkar/fq*dqh * (fq2m/vonkar - fq/vonkar)
      else
            tref = thm
            qref = qm
      endif

      tref = max(200.0, min(350.0, tref))
      qref = max(0.0,   min(0.1,   qref))

      contains

      ! Internal subroutine: compute scalar roughness lengths
      ! Andreas (1987) [Ref. 3] / Smeets and Van den Broeke (2008) [Ref. 4]
      subroutine calculate_scalar_roughness(z0m, ustar_in, z0h_out, z0q_out)
      implicit none
      real(r8), intent(in)  :: z0m, ustar_in
      real(r8), intent(out) :: z0h_out, z0q_out
      real(r8) :: Re_star_local, ln_z0h_z0m_local

      Re_star_local = ustar_in * z0m / nu_air

      if (Re_star_local < 0.135) then
            ! Smooth surface regime - Andreas (1987) [Ref. 3]
            z0h_out = z0m * exp(-vonkar * (Pr**(2.0/3.0) - 1.0) / Pr**(1.0/3.0))
            z0q_out = z0m * exp(-vonkar * (Sc**(2.0/3.0) - 1.0) / Sc**(1.0/3.0))
      else
            ! Rough surface regime - Smeets and Van den Broeke (2008) [Ref. 4]
            ln_z0h_z0m_local = -2.0 - 1.29 * sqrt(Re_star_local)
            ln_z0h_z0m_local = max(ln_z0h_z0m_local, -10.0)
            z0h_out = z0m * exp(ln_z0h_z0m_local)
            z0q_out = z0h_out
      endif

      z0h_out = max(1.0e-9, min(z0m, z0h_out))
      z0q_out = max(1.0e-9, min(z0m, z0q_out))

      end subroutine calculate_scalar_roughness

      ! Internal subroutine: adjust wind speed for stability
      ! Convective velocity scale - Deardorff (1970) [Ref. 6]
      subroutine adjust_wind_speed(ur_in, ustar_in, thvstar_in, thv_in, zeta_in, zii_in, beta_in, um_out)
      implicit none
      real(r8), intent(in)  :: ur_in, ustar_in, thvstar_in, thv_in, zeta_in, zii_in, beta_in
      real(r8), intent(out) :: um_out
      real(r8) :: w_star_cubed_local, wc_local, wc2_local

      if (zeta_in >= 0.0) then
            ! Stable: mechanical turbulence only
            um_out = max(ur_in, MIN_WIND)
      else
            ! Unstable: include convective contribution - Deardorff (1970) [Ref. 6]
            w_star_cubed_local = grav * ustar_in * (-thvstar_in) * zii_in / thv_in
            if (w_star_cubed_local > 0.0) then
                  wc_local  = w_star_cubed_local**(1.0/3.0)
                  wc2_local = (beta_in * wc_local)**2
                  um_out = sqrt(ur_in*ur_in + wc2_local)
            else
                  um_out = max(ur_in, MIN_WIND)
            endif
      endif
      um_out = max(MIN_WIND, min(100.0, um_out))

      end subroutine adjust_wind_speed

      end subroutine groundfluxes_glacier


      ! ======================================================================
      ! groundtem_glacier
      ! Purpose: Snow and land ice temperature calculation using Crank-Nicolson scheme.
      !
      ! Authors: Shuwen Li and Xing Yuan
      ! Modified by Shuwen Li (swli)
      !   * 04/2025 - Replaced fixed bottom-boundary temperature with a
      !               constant geothermal heat-flux boundary condition
      !               (default 0.06 W/m2 from Davies and Davies, 2010
      !               [Ref. 8])
      !   * 04/2025 - Fixed right-hand-side vector construction in the
      !               tridiagonal Crank-Nicolson system to use tss_bef(j)
      !               as the base state at every layer
      !
      ! References:
      ! Li, S., and X. Yuan (2026): Evaluation of a Land Surface-Glacier
!     ! Coupled Model over the Three-River Headwaters Region. Water, 18, 1030.
      ! ======================================================================
      SUBROUTINE groundtem_glacier(scwater, lb, nl_ice, dtime, tm, ht, &
                              capr, cnfac, dz, z, zi, nl_snow, &
                              tss, tss_bef, wice, wliq, scv, snowdp, &
                              frl, sabg, sabvg, fseng, fevpg, cgrnd, cgrnds, cgrndl, htvp, emg, tg, &
                              imelt, sm, xmf, melt, smelt, refreeze, fact, pg_rain, pg_snow, t_precip, &
                              ipr, jpr, geo_heat)

      USE shr_kind_mod, ONLY: r8 => shr_kind_rn
      use module_sf_clm_constant, only : stefnc, cpice, cpliq, denliq, denice, tfrz, tkliq, tkice, tkair, mu_i, mu_s

      IMPLICIT NONE

      ! Input variables
      integer, INTENT(in)  :: scwater     ! Land patch type (3 = land ice)
      integer              :: ipr, jpr, nl_snow
      integer, INTENT(in)  :: lb          ! Lower array bound
      integer, INTENT(in)  :: nl_ice      ! Upper array bound
      real(r8), INTENT(in) :: dtime       ! Time step [s]
      real(r8), INTENT(in) :: capr        ! Surface temperature tuning factor
      real(r8), INTENT(in) :: cnfac       ! Crank-Nicolson factor (0 to 1)
      real(r8), INTENT(in) :: ht          ! Temperature observation height [m]
      real(r8), INTENT(in) :: tm          ! Air temperature at reference height [K]
      real(r8), INTENT(in) :: dz(lb:nl_ice)     ! Layer thickness [m]
      real(r8), INTENT(in) :: z (lb:nl_ice)     ! Node depth [m]
      real(r8), INTENT(in) :: zi(lb-1:nl_ice)   ! Interface depth [m]
      real(r8), INTENT(inout) :: tg       ! Temperature of the ice or snow surface [K]
      real(r8), INTENT(in) :: sabg        ! Solar absorbed by ground [W/m2]
      real(r8), INTENT(in) :: frl         ! Atmospheric longwave radiation [W/m2]
      real(r8), INTENT(in) :: fseng       ! Sensible heat flux from ground [W/m2]
      real(r8), INTENT(in) :: fevpg       ! Evaporation flux from ground [mm/s]
      real(r8), INTENT(in) :: cgrnd, cgrnds, cgrndl
      real(r8), INTENT(in) :: htvp        ! Latent heat [J/kg]
      real(r8), INTENT(in) :: emg         ! Ground emissivity
      real(r8), INTENT(in) :: t_precip    ! Precipitation temperature [K]
      real(r8), INTENT(in) :: pg_rain     ! Rainfall [kg/(m2 s)]
      real(r8), INTENT(in) :: pg_snow     ! Snowfall [kg/(m2 s)]
      real(r8), INTENT(in) :: sabvg       ! Solar absorbed by vegetation [W/m2]
      real(r8), INTENT(in), OPTIONAL :: geo_heat  ! Geothermal heat flux [W/m2]

      real(r8), INTENT(inout) :: tss_bef(lb:nl_ice)   ! Temperature at time step before
      real(r8), INTENT(inout) :: tss (lb:nl_ice)       ! Temperature at time step
      real(r8), INTENT(inout) :: wice(lb:nl_ice)       ! Ice water content [kg/m²]
      real(r8), INTENT(inout) :: wliq(lb:nl_ice)       ! Liquid water content [kg/m²]
      real(r8), INTENT(inout) :: scv                  ! Snow water equivalent [kg/m²]
      real(r8), INTENT(inout) :: snowdp               ! Snow depth [m]

      ! Output variables
      real(r8), INTENT(out) :: sm      ! Snow melt [mm/s]
      real(r8), INTENT(out) :: xmf     ! Snow mass flux [kg/(m² s)]
      real(r8), INTENT(out) :: melt    ! Total melt [mm]
      real(r8), INTENT(out) :: smelt   ! Snowmelt [mm/s]
      real(r8), INTENT(out) :: refreeze  ! Refreezing rate [mm/s]
      real(r8), INTENT(out) :: fact(lb:nl_ice)   ! Temperature correction factor
      integer, INTENT(out) :: imelt(nl_snow+1:nl_ice)   ! Melt flag for each snow layer

      ! Local variables
      real(r8) :: rhosnow(lb:0)       ! Snow density [kg/m³]
      real(r8) :: cv(lb:nl_ice)       ! Heat capacity of ice or water [J/kg K]
      real(r8) :: thk(lb:nl_ice)      ! Thickness of layers [m]
      real(r8) :: tk(lb:nl_ice)       ! Temperature at each layer [K]
      real(r8) :: at(lb:nl_ice), bt(lb:nl_ice), ct(lb:nl_ice), rt(lb:nl_ice)  ! Coefficients for Crank-Nicolson scheme
      real(r8) :: fn(lb:nl_ice), fn1(lb:nl_ice)   ! Crank-Nicolson variables
      real(r8) :: dzm, dzp            ! Layer thickness changes
      real(r8) :: tss_bef_loc(lb:nl_ice)   ! Temporary storage for temperature before change
      real(r8) :: hs, dhsdT           ! Snow or ice enthalpy and temperature gradient
      real(r8) :: brr(lb:nl_ice)      ! Heat flux coefficient
      real(r8) :: d_tss(lb:nl_ice)    ! Change in temperature
      real(r8) :: thk_1              ! Thickness of previous layer
      real(r8) :: T_base             ! Base temperature for melt/refreeze
      real(r8) :: bottom_heat_flux   ! Heat flux at the bottom of the ice layer
      real(r8) :: I_itf              ! Ice temperature factor
      real(r8) :: sabg_absorbed      ! Absorbed solar radiation
      real(r8), dimension(lb:nl_ice) :: dIdz   ! Solar radiation gradient across depth
      real(r8) :: z_depth, z_depth_prev   ! Depth and previous depth for gradient calculations
      real(r8) :: I_bottom, I_top     ! Solar radiation at the bottom and top of the layer
      real(r8) :: DZ_INTER           ! Layer interface depth for Crank-Nicolson
      integer  :: i, j                ! Loop indices
      logical  :: use_legacy_melt = .false.  ! Flag for legacy melt computation


      tss_bef_loc(:) = tss_bef(:)

      ! ======================================================================
      ! 1. Heat capacity initialization
      ! ======================================================================
      cv(lb:) = wice(lb:)*cpice + wliq(lb:)*cpliq
      rhosnow(lb:0) = 0.0_r8
      if (lb == 1 .and. scv > 0.0_r8) cv(1) = cv(1) + cpice*scv
      if (lb <= 0) cv(:0) = cpliq*wliq(:0) + cpice*wice(:0)
      do j = lb, nl_ice
            if (cv(j) == 0.) then
                  wliq(j) = denliq*dz(j)*0.03
                  wice(j) = denice*dz(j)*0.9
                  cv(j)   = wice(j)*cpice + wliq(j)*cpliq
            endif
      enddo

      DZ_INTER = 0.05_r8
      T_base   = tfrz - 6.5
      dIdz(:)  = 0.0_r8

      ! Set bottom geothermal heat flux
      IF (PRESENT(geo_heat)) THEN
            bottom_heat_flux = geo_heat
      ELSE
            ! Default value: Davies and Davies (2010) [Ref. 8]
            bottom_heat_flux = 0.06_r8
      ENDIF

      ! ======================================================================
      ! 2. Thermal conductivity (ice layers)
      ! ======================================================================
      do j = lb, nl_ice
            if (tss(j) <= tfrz) then
                  thk(j) = tkice
            else
                  thk(j) = tkliq
            end if
      end do

      ! Snow layer thermal conductivity (density-dependent, Sturm et al. 1997 [Ref. 7])
      if (lb <= 0) then
            do j = lb, 0
                  rhosnow(j) = (wice(j) + wliq(j)) / max(dz(j), 1.0e-6_r8)
                  rhosnow(j) = min(rhosnow(j), 917.0_r8)
                  thk(j) = 0.029_r8 * (1.0_r8 + 1.0e-4_r8 * rhosnow(j)**2)
                  thk(j) = min(thk(j), 2.5_r8)  ! Cap at 2.5 to prevent thermal runaway
            end do
      end if

      ! ======================================================================
      ! 3. Interface thermal conductivity
      ! ======================================================================
      do j = lb, nl_ice-1
            if ((j==0) .and. (z(j+1)-zi(j) < zi(j)-z(j))) then
                  tk(j) = 2.0_r8 * thk(j) * thk(j+1) / max(thk(j)+thk(j+1), 1.0e-6_r8)
                  tk(j) = max(0.5_r8*thk(j+1), tk(j))
            else
                  tk(j) = thk(j)*thk(j+1)*(z(j+1)-z(j)) / &
                        max(thk(j)*(z(j+1)-zi(j)) + thk(j+1)*(zi(j)-z(j)), 1.0e-6_r8)
            end if
      end do
      tk(nl_ice) = thk(nl_ice)

      ! ======================================================================
      ! 4. Radiation penetration (Beer-Lambert law; see Li and Yuan, 2026 [Ref. 2])
      ! ======================================================================
      if (lb <= 0) then
            I_itf = sabg * exp(-mu_s * snowdp)
      else
            I_itf = sabg
      end if

      sabg_absorbed = 0.0_r8

      do j = lb, nl_ice
            if (j <= 0) then  ! Snow layers
                  z_depth = 0.0_r8
                  do i = lb, j
                        z_depth = z_depth + dz(i)
                  end do
                  z_depth_prev = z_depth - dz(j)
                  if (dz(j) < 1.0e-8_r8) then
                        dIdz(j) = 0.0_r8
                  else
                        I_top    = sabg * exp(-mu_s * z_depth_prev)
                        I_bottom = sabg * exp(-mu_s * z_depth)
                        dIdz(j)  = (I_top - I_bottom) / dz(j)
                  end if
            else              ! Ice layers
                  z_depth = 0.0_r8
                  if (lb <= 0) then
                        do i = 1, j
                              z_depth = z_depth + dz(i)
                        end do
                  else
                        do i = lb, j
                              z_depth = z_depth + dz(i)
                        end do
                  end if
                  z_depth_prev = z_depth - dz(j)
                  if (dz(j) < 1.0e-8_r8) then
                        dIdz(j) = 0.0_r8
                  else
                        I_top    = I_itf * exp(-mu_i * z_depth_prev)
                        I_bottom = I_itf * exp(-mu_i * z_depth)
                        dIdz(j)  = (I_top - I_bottom) / dz(j)
                  end if
            end if
      end do

      sabg_absorbed = dIdz(lb) * dz(lb)

      ! ======================================================================
      ! 5. Factorial factors for numerical solution
      ! ======================================================================
      j = lb
      thk_1 = 0.5*(z(j)-zi(j-1) + capr*(z(j+1)-zi(j-1)))
      thk_1 = max(thk_1, 1.0e-2_r8)

      do j = lb, nl_ice
            fact(j) = dtime / cv(j)
      end do

      ! ======================================================================
      ! 6. Initial heat fluxes
      ! ======================================================================
      do j = lb, nl_ice-1
            d_tss(j) = tss(j+1) - tss(j)
            dzp = max(z(j+1) - z(j), 1.0e-6_r8)
            fn(j) = tk(j) * d_tss(j) / dzp
      end do
      fn(nl_ice) = bottom_heat_flux

      ! ======================================================================
      ! 7. Surface energy balance
      ! ======================================================================
      hs = sabg_absorbed + emg*frl - emg*stefnc*tg**4 - (fseng + fevpg*htvp) &
         + cpliq*pg_rain*(t_precip - tg) + cpice*pg_snow*(t_precip - tg)

      dhsdT = -cgrnd - 4.0_r8*emg*stefnc*tg**3 - cpliq*pg_rain - cpice*pg_snow

      ! ======================================================================
      ! 8. Tridiagonal matrix setup (Crank-Nicolson)
      ! ======================================================================
      ! Surface layer
      j   = lb
      dzp = z(j+1) - z(j)
      at(j) = 0.
      bt(j) = 1. + (1.-cnfac)*fact(j)*tk(j)/dzp - cnfac*fact(j)*dhsdT
      ct(j) =    -(1.-cnfac)*fact(j)*tk(j)/dzp
      rt(j) = tss_bef(j) + fact(j)*(hs - cnfac*dhsdT*tss_bef(j) + cnfac*fn(j))

      ! Internal layers
      do j = lb+1, nl_ice-1
            dzp = max(z(j+1)-z(j), 1.0e-6_r8)
            dzm = max(z(j)-z(j-1), 1.0e-6_r8)
            at(j) = -(1.0_r8 - cnfac)*fact(j)*tk(j-1)/dzm
            bt(j) =  1.0_r8 + (1.0_r8 - cnfac)*fact(j)*(tk(j)/dzp + tk(j-1)/dzm)
            ct(j) = -(1.0_r8 - cnfac)*fact(j)*tk(j)/dzp
            rt(j) = tss_bef(j) + cnfac*fact(j)*(fn(j)-fn(j-1)) + fact(j)*dIdz(j)*dz(j)
      end do

      ! Bottom layer (constant heat flux boundary)
      j   = nl_ice
      dzm = max(z(j)-z(j-1), 1.0e-6_r8)
      at(j) = -(1.0_r8 - cnfac)*fact(j)*tk(j-1)/dzm
      bt(j) =  1.0_r8 + (1.0_r8 - cnfac)*fact(j)*tk(j-1)/dzm
      ct(j) =  0.0_r8
      rt(j) = tss_bef(j) + fact(j)*(cnfac*(fn(j)-fn(j-1))) + fact(j)*dIdz(j)*dz(j)

      ! ======================================================================
      ! 9. Solve tridiagonal system
      ! ======================================================================
      i = size(at)
      call tridia(i, at, bt, ct, rt, tss)

      if (tss(lb) /= tss(lb) .or. sabg /= sabg) then
            print *, 'INPUT ERROR: tg or sabg is NaN'
            stop
      end if

      ! ======================================================================
      ! 10. Final heat fluxes
      ! ======================================================================
      do j = lb, nl_ice-1
            dzp = max(z(j+1)-z(j), 1.0e-6_r8)
            d_tss(j) = tss(j+1) - tss(j)
            if (tss(j) > 400.0_r8 .or. tss(j) < 0.0_r8) then
                  print *, 'CRASH PREVENTION: tss value is physically impossible!'
                  print *, 'tss(',j,')=', tss(j)
                  print *, 'Location: ipr=', ipr, ' jpr=', jpr
            end if
            fn1(j) = tk(j)*d_tss(j)/dzp
      end do
      fn1(nl_ice) = bottom_heat_flux

      j = lb
      brr(j) = cnfac*fn(j) + (1.0-cnfac)*fn1(j)
      do j = lb+1, nl_ice
            brr(j) = cnfac*(fn(j)-fn(j-1)) + (1.0-cnfac)*(fn1(j)-fn1(j-1))
      end do

      ! ======================================================================
      ! 11. Phase change (melt and refreeze)
      ! ======================================================================
      if (use_legacy_melt) then
            call meltf_glc(scwater, lb, nl_ice, nl_snow, dtime, tm, ht, &
                        fact(lb:), brr(lb:), hs, dhsdT, tg, &
                        tss_bef(lb:), tss(lb:), wliq(lb:), wice(lb:), imelt(nl_snow+1:), &
                        scv, snowdp, sm, xmf, melt, smelt, refreeze, dz(1:), ipr, jpr)
      else
            call meltf_refreeze_glc(scwater, lb, nl_ice, nl_snow, dtime, tm, ht, &
                        fact(lb:), brr(lb:), hs, dhsdT, tg, &
                        tss_bef(lb:), tss(lb:), wliq(lb:), wice(lb:), imelt(nl_snow+1:), &
                        scv, snowdp, sm, xmf, melt, smelt, refreeze, dz(lb:), ipr, jpr)
      end if

      end subroutine groundtem_glacier


      ! ======================================================================
      ! meltf_refreeze_glc
      ! Phase change calculations (melt and refreeze) for glacier/snow layers.
      ! Authors: Shuwen Li and Xing Yuan
      ! Written by: Shuwen Li (swli), 2023 - 2026
      !   * New replacement for the legacy meltf_glc subroutine
      !   * Separates melt and refreeze mass accounting so that each is
      !     reported as a non-negative quantity independently
      !   * Refreeze energy is bounded by the layer cold content using the
      !     actual heat capacity (wice*cpice + wliq*cpliq), with an upper
      !     limit of 50% of the available liquid water per time step to
      !     prevent unphysical over-refreezing
      !
      ! References:
      ! Li, S., and X. Yuan (2026): Evaluation of a Land Surface-Glacier
!     ! Coupled Model over the Three-River Headwaters Region. Water, 18, 1030.
      ! ======================================================================
      SUBROUTINE meltf_refreeze_glc(scwater, lb, nl_ice, nl_snow, dtime, tm, ht, &
                  fact, brr, hs, dhsdT, tg, &
                  tss_bef, tss, wliq, wice, imelt, &
                  scv, snowdp, sm, xmf, melt, smelt, refreeze, dz, ipr, jpr)

      USE shr_kind_mod, ONLY: r8 => shr_kind_rn
      USE module_sf_clm_constant, ONLY : tfrz, hfus, denice, denliq, cpice, grav

      IMPLICIT NONE

      ! Input variables
      INTEGER,  INTENT(IN)    :: scwater         ! Land patch type (3 = land ice)
      INTEGER,  INTENT(IN)    :: lb              ! Lower array bound
      INTEGER,  INTENT(IN)    :: nl_ice           ! Upper array bound
      INTEGER,  INTENT(IN)    :: nl_snow          ! Number of snow layers
      REAL(r8), INTENT(IN)    :: dtime            ! Time step [s]
      REAL(r8), INTENT(IN)    :: tm               ! Temperature at reference height [K]
      REAL(r8), INTENT(IN)    :: ht               ! Observation height for temperature [m]
      REAL(r8), INTENT(IN)    :: fact(lb:)         ! Temperature correction factor [K]
      REAL(r8), INTENT(IN)    :: brr(lb:)          ! Heat flux coefficient [W/m²K]

      ! In-out variables
      REAL(r8), INTENT(INOUT) :: hs               ! Heat storage [J/m²]
      REAL(r8), INTENT(IN)    :: dhsdT            ! Temperature change in snow/ice layer [K/s]
      REAL(r8), INTENT(IN)    :: tg               ! Surface temperature of glacier [K]
      REAL(r8), INTENT(IN)    :: tss_bef(lb:)     ! Temperature before current timestep [K]
      REAL(r8), INTENT(INOUT) :: tss(lb:)         ! Current temperature in each layer [K]
      REAL(r8), INTENT(INOUT) :: wice(lb:)        ! Ice water content [kg/m²]
      REAL(r8), INTENT(INOUT) :: wliq(lb:)        ! Liquid water content [kg/m²]
      REAL(r8), INTENT(INOUT) :: scv              ! Snow water equivalent [kg/m²]
      REAL(r8), INTENT(INOUT) :: snowdp           ! Snow depth [m]
      REAL(r8), INTENT(OUT)   :: sm               ! Snow melt [mm/s]
      REAL(r8), INTENT(OUT)   :: xmf              ! Snow mass flux [kg/(m² s)]
      REAL(r8), INTENT(OUT)   :: melt             ! Total melt [mm]
      REAL(r8), INTENT(OUT)   :: smelt            ! Snow melt rate [mm/s]
      REAL(r8), INTENT(OUT)   :: refreeze         ! Refreeze rate [mm/s]
      INTEGER, INTENT(OUT)    :: imelt(nl_snow+1:) ! Melt flag for each snow layer

      ! Local variables
      REAL(r8) :: wmass0(lb:nl_ice)           ! Initial water mass in each layer [kg/m²]
      REAL(r8) :: wice0(lb:nl_ice)            ! Initial ice water content [kg/m²]
      REAL(r8) :: wliq0(lb:nl_ice)            ! Initial liquid water content [kg/m²]
      REAL(r8) :: hm(lb:nl_ice)               ! Heat transfer coefficient [W/m²K]
      REAL(r8) :: xm(lb:nl_ice)               ! Melt flux [kg/(m² s)]
      REAL(r8) :: heatr                       ! Heat energy [J]
      REAL(r8) :: temp1                       ! Temporary temperature for calculations [K]
      REAL(r8) :: propor                      ! Proportional factor for heat flux
      REAL(r8) :: tinc                        ! Time increment [s]
      REAL(r8) :: scvold                      ! Previous snow water equivalent [kg/m²]
      REAL(r8) :: we                          ! Energy available for melt [J/m²]
      REAL(r8) :: we_old                      ! Previous energy available [J/m²]
      REAL(r8) :: total_mass                  ! Total mass of the snow/ice layer [kg/m²]
      INTEGER  :: j                           ! Loop index
      REAL(r8) :: layer_heat_capacity         ! Heat capacity of the layer [J/m²K]
      REAL(r8) :: max_refreeze_energy        ! Maximum energy available for refreeze [J/m²]
      REAL(r8) :: max_refreeze_mass          ! Maximum mass that can refreeze [kg/m²]
      REAL(r8) :: layer_refreeze             ! Refreeze energy in the layer [J/m²]
      REAL(r8) :: layer_melt                 ! Melt energy in the layer [J/m²]
      REAL(r8) :: delta_temp                 ! Change in temperature due to melt or refreeze [K]
      REAL(r8) :: energy_available           ! Energy available for phase change [J]
      REAL(r8), PARAMETER :: retain_cap = 0.02_r8 ! Pore water retention coefficient
      REAL(r8) :: total_melt                 ! Total melt for the time step [mm]
      REAL(r8) :: total_refreeze             ! Total refreeze for the time step [mm]


      ! ======================================================================
      ! Initialization
      ! ======================================================================
      sm            = 0.0_r8
      xmf           = 0.0_r8
      total_melt    = 0.0_r8
      total_refreeze = 0.0_r8

      do j = lb, nl_ice
            imelt(j)  = 0
            hm(j)     = 0.0_r8
            xm(j)     = 0.0_r8
            wice0(j)  = wice(j)
            wliq0(j)  = wliq(j)
            wmass0(j) = wice(j) + wliq(j)
      enddo

      scvold = scv
      we     = 0.0_r8
      if (lb <= 0) we = sum(wice(lb:0) + wliq(lb:0))

      ! ======================================================================
      ! Phase change condition detection
      ! ======================================================================
      do j = lb, nl_ice
            ! Melt: ice above freezing
            if (tss(j) > tfrz .and. wice(j) > 0.0_r8) then
                  imelt(j) = 1
                  tss(j)   = tfrz
            endif
            ! Freeze: liquid below freezing
            if (imelt(j) == 0) then
                  if (wliq(j) > 0.0_r8 .and. tss(j) < tfrz) imelt(j) = 2
            endif
      enddo

      ! Thin snow layer special case
      if (lb == 1 .and. scv > 0.0_r8 .and. tss(lb) > tfrz) then
            imelt(1) = 1
            tss(1)   = tfrz
      endif

      ! ======================================================================
      ! Residual energy calculation
      ! ======================================================================
      do j = lb, nl_ice
            if (imelt(j) > 0) then
                  tinc = tss(j) - tss_bef(j)
                  if (j > lb) then
                        hm(j) = brr(j) - tinc / fact(j)
                  else
                        hm(j) = hs + dhsdT*tinc + brr(j) - tinc / fact(j)
                  endif
            endif
      enddo

      ! Clear invalid phase change flags
      do j = lb, nl_ice
            if (imelt(j) == 1 .and. hm(j) < 0.0_r8) then
                  hm(j) = 0.0_r8;  imelt(j) = 0
            else if (imelt(j) == 2 .and. hm(j) > 0.0_r8) then
                  hm(j) = 0.0_r8;  imelt(j) = 0
            endif
      enddo

      ! ======================================================================
      ! Phase change mass calculation
      ! ======================================================================
      heatr = 0.0_r8
      smelt = 0.0_r8

      do j = lb, nl_ice
            if (imelt(j) > 0 .and. abs(hm(j)) > 0.0_r8) then

                  xm(j) = hm(j) * dtime / hfus  ! Phase change mass [kg/m2]

                  ! Surface snow melt
                  if (j == lb .and. scv > 0.0_r8 .and. xm(j) > 0.0_r8) then
                        temp1 = scv
                        scv   = max(0.0_r8, temp1 - xm(j))
                        propor = scv / max(temp1, 1.e-10_r8)
                        snowdp = propor * snowdp

                        layer_melt  = temp1 - scv
                        total_melt  = total_melt + layer_melt
                        sm          = layer_melt / dtime

                        heatr = hm(j) - hfus * layer_melt / dtime
                        if (heatr > 0.0_r8) then
                              xm(j) = heatr * dtime / hfus
                              hm(j) = heatr
                        else
                              xm(j) = 0.0_r8
                              hm(j) = 0.0_r8
                        endif
                  endif

                  ! Refreeze (imelt=2)
                  if (imelt(j) == 2) then
                        layer_heat_capacity  = wice(j)*cpice + wliq(j)*cpliq
                        energy_available     = layer_heat_capacity * (tfrz - tss(j))
                        max_refreeze_mass    = energy_available / hfus
                        layer_refreeze       = min(max_refreeze_mass, wliq(j))
                        layer_refreeze       = max(0.0_r8, layer_refreeze)
                        layer_refreeze       = min(layer_refreeze, 0.5_r8 * wliq(j))

                        if (layer_refreeze > 1.e-10_r8) then
                              wice(j) = wice(j) + layer_refreeze
                              wliq(j) = wliq(j) - layer_refreeze
                              delta_temp = (layer_refreeze * hfus) / max(layer_heat_capacity, 1.e-10_r8)
                              tss(j)  = min(tfrz, tss(j) + delta_temp)
                              total_refreeze = total_refreeze + layer_refreeze
                        endif

                  ! Regular melt (imelt=1)
                  else if (imelt(j) == 1) then
                        if (xm(j) > 1.e-10_r8) then
                              layer_melt  = min(xm(j), wice0(j))
                              wice(j)     = wice0(j) - layer_melt
                              wliq(j)     = wliq0(j) + layer_melt
                              total_melt  = total_melt + layer_melt
                        endif
                  endif

            endif
      enddo

      ! ======================================================================
      ! Final output
      ! ======================================================================
      melt     = max(0.0_r8, total_melt)
      refreeze = max(0.0_r8, total_refreeze)
      xmf      = (total_melt - total_refreeze) * hfus / dtime

      if (tg > tfrz - 2.0_r8 .and. hs > 0.0_r8) then
            smelt = max(0.0_r8, sm * dtime)
      else
            smelt = 0.0_r8
      endif

      ! Mass conservation check for snow layers
      if (lb <= 0) then
            we_old = we
            we     = sum(wice(lb:0) + wliq(lb:0))
            total_mass = we
            if (abs(we - we_old) > max(1.e-3_r8, 1.e-6_r8 * total_mass)) then
                  print*, 'meltf_glc mass conservation error:', we - we_old, ' at:', ipr, jpr
            endif
      endif

      end subroutine meltf_refreeze_glc


      ! ======================================================================
      ! GLACIER_WATER
      ! Prupose: Update liquid water within snow layers and water flux onto ice surface.
      ! Snow meltwater percolates downward. 
      ! Authors: Shuwen Li and Xing Yuan
      ! Modified by: Shuwen Li 2025
      ! Once underlying snow is saturated,excess meltwater runs off. Meltwater on glacier 
      ! ice stays in place until refreezing.
      !
      ! References:
      !   [1] Yuan, X., C. Li, P. Ji, X.Z. Xi, S.W. Li, Y. Jiao, H. Yi, X.Z. Liang,
      !       and Q.C. Zeng (2026): Development and evaluation of Conjunctive
      !       Surface-Subsurface Process Model version 3 (CSSPv3) at global scale.
      !       Journal of Advances in Modeling Earth Systems. (Revised)
      !   [2] Li, S., and X. Yuan (2026): Evaluation of a Land Surface-Glacier
      !       Coupled Model over the Three-River Headwaters Region. Water, 18, 1030.
      ! ======================================================================
      SUBROUTINE GLACIER_WATER(nl_ice, nl_snow, dtime, &
                              z, dz, zi, tss, &
                              wliq, wice, pg_rain, pg_snow, &
                              sm, scv, snowdp, imelt, &
                              fiold, snl, qseva, qsdew, &
                              qsubl, qfros, gwat, &
                              us, vs, melt, ipr, jpr)

      USE shr_kind_mod, ONLY: r8 => shr_kind_rn
      USE module_sf_clm_constant, ONLY : denice, denliq, tfrz, wimp, wssi

      IMPLICIT NONE

      ! Input variables
      INTEGER,  INTENT(IN)    :: nl_ice        ! Upper bound of ice/snow column
      INTEGER,  INTENT(IN)    :: nl_snow       ! Number/index of snow layers
      INTEGER,  INTENT(IN)    :: ipr           ! Grid index in x-direction
      INTEGER,  INTENT(IN)    :: jpr           ! Grid index in y-direction

      REAL(r8), INTENT(IN)    :: dtime         ! Time step [s]
      REAL(r8), INTENT(IN)    :: pg_rain       ! Rainfall rate [kg m-2 s-1]
      REAL(r8), INTENT(IN)    :: pg_snow       ! Snowfall rate [kg m-2 s-1]
      REAL(r8), INTENT(IN)    :: sm            ! Surface melt rate or snow melt source term
      REAL(r8), INTENT(IN)    :: qseva         ! Ground evaporation flux [kg m-2 s-1]
      REAL(r8), INTENT(IN)    :: qsdew         ! Dew deposition flux [kg m-2 s-1]
      REAL(r8), INTENT(IN)    :: qsubl         ! Sublimation flux [kg m-2 s-1]
      REAL(r8), INTENT(IN)    :: qfros         ! Frost deposition flux [kg m-2 s-1]
      REAL(r8), INTENT(IN)    :: melt          ! Meltwater production
      REAL(r8), INTENT(IN)    :: us            ! Zonal wind speed [m s-1]
      REAL(r8), INTENT(IN)    :: vs            ! Meridional wind speed [m s-1]

      REAL(r8), INTENT(IN)    :: fiold(nl_snow+1:nl_ice)   ! Previous ice fraction or layer ice state
      INTEGER,  INTENT(IN)    :: imelt(nl_snow+1:nl_ice)   ! Melt flag for each layer

      ! In-out variables
      INTEGER,  INTENT(INOUT) :: snl           ! Number of active snow layers

      REAL(r8), INTENT(INOUT) :: z(nl_snow+1:nl_ice)       ! Layer midpoint depth [m]
      REAL(r8), INTENT(INOUT) :: dz(nl_snow+1:nl_ice)      ! Layer thickness [m]
      REAL(r8), INTENT(INOUT) :: zi(nl_snow:nl_ice)        ! Layer interface depth [m]
      REAL(r8), INTENT(INOUT) :: tss(nl_snow+1:nl_ice)     ! Snow/ice layer temperature [K]
      REAL(r8), INTENT(INOUT) :: wice(nl_snow+1:nl_ice)    ! Ice mass in each layer [kg m-2]
      REAL(r8), INTENT(INOUT) :: wliq(nl_snow+1:nl_ice)    ! Liquid water mass in each layer [kg m-2]
      REAL(r8), INTENT(INOUT) :: scv                       ! Snow water equivalent [kg m-2]
      REAL(r8), INTENT(INOUT) :: snowdp                    ! Snow depth [m]

      ! Output variables
      REAL(r8), INTENT(OUT)   :: gwat          ! Net water input from the surface [mm s-1]

      ! Local variables
      INTEGER                 :: lb            ! Lower bound of active snow/ice column

      ! ======================================================================
      ! [1] Update snow layer water content and compute surface water flux
      ! NOTE: gwat is reassigned in the if/else block below; the first
      ! assignment here is a fallback for the bare-ice case only.
      ! ======================================================================
      gwat = pg_rain - qseva - qsubl + melt/dtime
      lb   = snl + 1

      if (lb == 1) then
            ! Bare ice surface
            gwat = pg_rain + melt/dtime - (qseva + qsubl)
      else
            call snowwater(lb, dtime, &
                           pg_rain, qseva, qsdew, qsubl, qfros, &
                           dz(lb:0), wice(lb:0), wliq(lb:0), gwat)
      endif

      ! ======================================================================
      ! [2] Snow compaction, layer combination and division
      ! ======================================================================
      if (snl <= 0) then
            lb = snl + 1
            call snowcompaction_glc(lb, dtime, &
                                    imelt(lb:0), fiold(lb:0), tss(lb:0), &
                                    wliq(lb:0), wice(lb:0), us, vs, dz(lb:0))

            lb = nl_snow + 1
            call snowlayerscombine_glc(lb, snl, &
                                       z(lb:1), dz(lb:1), zi(lb-1:1), &
                                       wliq(lb:1), wice(lb:1), tss(lb:1), scv, snowdp, ipr, jpr)

            if (snl < 0) &
                  call snowlayersdivide_glc(lb, snl, &
                                            z(lb:0), dz(lb:0), zi(lb-1:0), &
                                            wliq(lb:0), wice(lb:0), tss(lb:0))
      endif

      if (snl == 0) then
            wice(nl_snow+1:snl) = 0.
            wliq(nl_snow+1:snl) = 0.
            tss (nl_snow+1:snl) = 0
            z   (nl_snow+1:snl) = 0.
            dz  (nl_snow+1:snl) = 0.
      endif

      if (lb >= 1) then
            wliq(1) = max(1.e-8, wliq(1) + qsdew * dtime)
            wice(1) = max(1.e-8, wice(1) + (qfros - qsubl) * dtime)
      end if
      gwat   = gwat + wliq(1) / dtime
      wliq(1) = 0

      end subroutine GLACIER_WATER


      ! ======================================================================
      ! snowwater
      ! Purpose: Renew ice and liquid water mass in snow layers.
      ! Meltwater percolates with bypass flow to prevent density explosion.
      !
      ! Authors: Shuwen Li and Xing Yuan
      ! Modified by Shuwen Li (swli), 2023 - 2026:
      !   * Added a density safety valve that removes excess liquid water
      !     when a layer's total mass exceeds 1000 kg/m3 after inflow
      !   * Introduced preferential bypass flow: meltwater unable to enter
      !     a saturated layer below is routed directly to runoff instead of
      !     accumulating in the current layer
      !
      ! References:
      ! Li, S., and X. Yuan (2026): Evaluation of a Land Surface-Glacier
!     ! Coupled Model over the Three-River Headwaters Region. Water, 18, 1030.
      ! ======================================================================
      SUBROUTINE GLACIER_WATER(nl_ice, nl_snow, dtime, &
                              z, dz, zi, tss, &
                              wliq, wice, pg_rain, pg_snow, &
                              sm, scv, snowdp, imelt, &
                              fiold, snl, qseva, qsdew, &
                              qsubl, qfros, gwat, &
                              us, vs, melt, ipr, jpr)

      USE shr_kind_mod, ONLY: r8 => shr_kind_rn
      USE module_sf_clm_constant, ONLY : denice, denliq, tfrz, wimp, wssi

      IMPLICIT NONE

      ! Input variables
      INTEGER,  INTENT(IN)    :: nl_ice        ! Upper bound of ice/snow column
      INTEGER,  INTENT(IN)    :: nl_snow       ! Number/index of snow layers
      INTEGER,  INTENT(IN)    :: ipr           ! Grid index in x-direction
      INTEGER,  INTENT(IN)    :: jpr           ! Grid index in y-direction

      REAL(r8), INTENT(IN)    :: dtime         ! Time step [s]
      REAL(r8), INTENT(IN)    :: pg_rain       ! Rainfall rate [kg m-2 s-1]
      REAL(r8), INTENT(IN)    :: pg_snow       ! Snowfall rate [kg m-2 s-1]
      REAL(r8), INTENT(IN)    :: sm            ! Surface melt rate or snow melt source term
      REAL(r8), INTENT(IN)    :: qseva         ! Ground evaporation flux [kg m-2 s-1]
      REAL(r8), INTENT(IN)    :: qsdew         ! Dew deposition flux [kg m-2 s-1]
      REAL(r8), INTENT(IN)    :: qsubl         ! Sublimation flux [kg m-2 s-1]
      REAL(r8), INTENT(IN)    :: qfros         ! Frost deposition flux [kg m-2 s-1]
      REAL(r8), INTENT(IN)    :: melt          ! Meltwater production
      REAL(r8), INTENT(IN)    :: us            ! Zonal wind speed [m s-1]
      REAL(r8), INTENT(IN)    :: vs            ! Meridional wind speed [m s-1]

      REAL(r8), INTENT(IN)    :: fiold(nl_snow+1:nl_ice)   ! Previous ice fraction or layer ice state
      INTEGER,  INTENT(IN)    :: imelt(nl_snow+1:nl_ice)   ! Melt flag for each layer

      ! In-out variables
      INTEGER,  INTENT(INOUT) :: snl           ! Number of active snow layers

      REAL(r8), INTENT(INOUT) :: z(nl_snow+1:nl_ice)       ! Layer midpoint depth [m]
      REAL(r8), INTENT(INOUT) :: dz(nl_snow+1:nl_ice)      ! Layer thickness [m]
      REAL(r8), INTENT(INOUT) :: zi(nl_snow:nl_ice)        ! Layer interface depth [m]
      REAL(r8), INTENT(INOUT) :: tss(nl_snow+1:nl_ice)     ! Snow/ice layer temperature [K]
      REAL(r8), INTENT(INOUT) :: wice(nl_snow+1:nl_ice)    ! Ice mass in each layer [kg m-2]
      REAL(r8), INTENT(INOUT) :: wliq(nl_snow+1:nl_ice)    ! Liquid water mass in each layer [kg m-2]
      REAL(r8), INTENT(INOUT) :: scv                       ! Snow water equivalent [kg m-2]
      REAL(r8), INTENT(INOUT) :: snowdp                    ! Snow depth [m]

      ! Output variables
      REAL(r8), INTENT(OUT)   :: gwat          ! Net water input from the surface [mm s-1]

      ! Local variables
      INTEGER                 :: lb            ! Lower bound of active snow/ice column
      ! Update ice and liquid mass at top layer
      wgdif = wice(lb) + (qfros - qsubl)*dtime
      wice(lb) = wgdif
      if (wgdif < 0.) then
            wice(lb) = 0.
            wliq(lb) = wliq(lb) + wgdif
      endif
      wliq(lb) = max(0., wliq(lb) + (pg_rain + qsdew - qseva)*dtime)

      ! Compute porosity and volume fractions (based on pre-update values)
      do j = lb, 0
            vol_ice(j)      = min(1., wice(j)/(dz(j)*denice))
            eff_porosity(j) = max(0.01, 1. - vol_ice(j))
            vol_liq(j)      = min(eff_porosity(j), wliq(j)/(dz(j)*denliq))
      enddo

      qin          = 0.0_r8
      total_runoff = 0.0_r8

      do j = lb, 0
            wliq(j) = wliq(j) + qin

            ! Safety valve: prevent density explosion from incoming flow
            max_allowed_mass = 1000.0_r8 * max(dz(j), 1.0e-6_r8)
            if ((wice(j) + wliq(j)) > max_allowed_mass) then
                  excess_mass  = (wice(j) + wliq(j)) - max_allowed_mass
                  removed_mass = min(wliq(j), excess_mass)
                  wliq(j)      = wliq(j) - removed_mass
                  total_runoff = total_runoff + removed_mass
            end if

            if (j <= -1) then  ! Internal snow layers
                  if (eff_porosity(j) < wimp .OR. eff_porosity(j+1) < wimp) then
                        qout = 0.
                  else
                        ! Potential outflow based on layer holding capacity
                        q_pot_m  = max(0._r8, (vol_liq(j) - wssi*eff_porosity(j))*dz(j))
                        ! Space available in the layer below
                        space_next_m = max(0._r8, (1._r8 - vol_ice(j+1) - vol_liq(j+1))*dz(j+1))
                        ! Water that can flow to the layer below
                        q_real_m     = min(q_pot_m, space_next_m)
                        ! Water that bypasses (large-pore preferential flow)
                        q_bypass_m   = q_pot_m - q_real_m
                        qout         = q_real_m * 1000.0_r8
                        removed_mass = q_bypass_m * 1000.0_r8
                        wliq(j)      = wliq(j) - qout - removed_mass
                        total_runoff = total_runoff + removed_mass
                        qin          = qout
                  endif
            else  ! Bottom snow layer (j=0): outflow to glacier surface
                  qout    = max(0., (vol_liq(j) - wssi*eff_porosity(j))*dz(j))
                  qout    = qout * 1000.
                  wliq(j) = wliq(j) - qout
                  qin     = qout
            endif
      enddo

      qout_snowb = (qin + total_runoff) / dtime

      end subroutine snowwater


      ! ======================================================================
      ! snowcompaction_glc
      ! Purpose: Compute snow compaction from four mechanisms:
      !   1. Destructive metamorphism (SNTHERM.89/99)
      !   2. Overburden (Vionnet et al., 2012 [Ref. 10])
      !   3. Melt metamorphism
      !   4. Wind drift (van Kampenhout et al., 2017 [Ref. 11])
      !
      ! Original author: Shuwen Li and Xing Yuan
      ! ======================================================================
      subroutine snowcompaction_glc(lb, dtime, imelt, fiold, &
                  tss, wliq, wice, us, vs, dz)

      USE shr_kind_mod, ONLY: r8 => shr_kind_rn
      use module_sf_clm_constant, only : denice, denliq, tfrz
      implicit none

      integer,  INTENT(in)    :: lb
      real(r8), INTENT(in)    :: dtime
      integer,  INTENT(in)    :: imelt(lb:0)
      real(r8), INTENT(in)    :: fiold(lb:0)
      real(r8), INTENT(in)    :: tss(lb:0)
      real(r8), INTENT(in)    :: wice(lb:0), wliq(lb:0)
      real(r8), INTENT(in)    :: us, vs
      real(r8), INTENT(inout) :: dz(lb:0)

      integer j

      real(r8), parameter :: c2  = 23.0e-3    ! [m3/kg]
      real(r8), parameter :: c3  = 2.777e-6   ! [1/s]
      real(r8), parameter :: c4  = 0.04       ! [1/K]
      real(r8), parameter :: c5  = 2.0        ! Liquid water enhancement factor
      real(r8), parameter :: dm  = 100.0      ! Destructive metamorphism density limit [kg/m3]

      real(r8) :: burden, ddz1, ddz2, ddz3, ddz4
      real(r8) :: dexpf, fi, td, pdzdtc, void, wx, bi
      real(r8) :: zpseudo, f1, f2, eta, forc_wind
      logical  :: mobile

      burden  = 0.0
      zpseudo = 0.0
      mobile  = .true.

      do j = lb, 0
            wx   = wice(j) + wliq(j)
            void = 1.0 - (wice(j)/denice + wliq(j)/denliq)/dz(j)

            if (void <= 0.001 .or. wice(j) <= .1) then
                  burden = burden + wx
                  mobile = .false.
                  CYCLE
            endif

            bi    = wice(j) / dz(j)
            fi    = wice(j) / wx
            td    = tfrz - tss(j)
            dexpf = exp(-c4*td)

            ! Destructive metamorphism
            ddz1 = -c3*dexpf
            if (bi > dm) ddz1 = ddz1*exp(-46.0e-3*(bi-dm))
            if (wliq(j) > 0.01*dz(j)) ddz1 = ddz1*c5

            ! Overburden compaction - Vionnet et al. (2012) [Ref. 10]
            f1  = 1.0/(1.0 + 60.0*wliq(j)/(denliq*dz(j)))
            f2  = 4.0
            eta = f1*f2*(bi/450.0)*exp(0.1*td + c2*bi)*7.62237e6
            ddz2 = -(burden + wx/2.0) / eta

            ! Melt compaction
            if (imelt(j) == 1 .and. fiold(j) /= 0.) then
                  ddz3 = -1.0/dtime * max(0.0, (fiold(j) - fi)/fiold(j))
            else
                  ddz3 = 0.0
            endif

            ! Wind drift compaction - van Kampenhout et al. (2017) [Ref. 11]
            forc_wind = sqrt(us**2 + vs**2)
            call winddriftcompaction_glc(bi, forc_wind, dz(j), zpseudo, mobile, ddz4)

            pdzdtc = ddz1 + ddz2 + ddz3 + ddz4
            dz(j)  = dz(j)*(1.0 + pdzdtc*dtime)
            dz(j)  = max(dz(j), (wice(j)/denice + wliq(j)/denliq))

            burden = burden + wx
      end do

      end subroutine snowcompaction_glc


      ! ======================================================================
      ! snowlayerscombine_glc
      ! Purpose: Combine thin snow layers with neighbors when below minimum thickness.
      ! Authors: Shuwen Li and Xing Yuan
      ! ======================================================================
      subroutine snowlayerscombine_glc(lb, snl, &
                  z, dz, zi, wliq, wice, tss, scv, snowdp, ipr, jpr)

      USE shr_kind_mod, ONLY: r8 => shr_kind_rn
      implicit none

      integer,  INTENT(in)    :: lb
      real(r8), INTENT(inout) :: wice(lb:1), wliq(lb:1), tss(lb:1)
      real(r8), INTENT(inout) :: dz(lb:1), z(lb:1), zi(lb-1:1)
      real(r8), INTENT(inout) :: snowdp, scv
      integer,  INTENT(inout) :: snl
      integer,  INTENT(in)    :: ipr, jpr

      real(r8) :: drr, zwice, zwliq
      real(r8) :: dzmin(5)
      integer  :: i, j, k, l, msn_old, mssi, neibor

      data dzmin /0.010, 0.015, 0.025, 0.055, 0.115/

      ! Combine layers with negligible ice mass with the layer below
      msn_old = snl
      do j = msn_old+1, 0
            if (wice(j) <= .1) then
                  wliq(j+1) = wliq(j+1) + wliq(j)
                  wice(j+1) = wice(j+1) + wice(j)
                  if (j > snl+1 .AND. snl < -1) then
                        do i = j, snl+2, -1
                              tss(i)  = tss(i-1)
                              wliq(i) = wliq(i-1)
                              wice(i) = wice(i-1)
                              dz(i)   = dz(i-1)
                        enddo
                  endif
                  snl = snl + 1
            endif
      enddo

      if (snl == 0) then
            scv    = 0.
            snowdp = 0.
            return
      else
            scv    = 0.
            snowdp = 0.
            zwice  = 0.
            zwliq  = 0.
            do j = snl+1, 0
                  scv    = scv    + wice(j)
                  snowdp = snowdp + dz(j)
                  zwice  = zwice  + wice(j)
                  zwliq  = zwliq  + wliq(j)
            enddo
      endif

      if (snowdp < 0.01) then  ! All snow gone
            snl = 0
            scv = zwice
            if (scv <= 0.) snowdp = 0.
            wliq(1) = wliq(1) + zwliq
            return
      else
            if (snl < -1) then
                  msn_old = snl
                  mssi    = 1
                  do i = msn_old+1, 0
                        if (dz(i) < dzmin(mssi)) then
                              if (i == snl+1) then
                                    neibor = i + 1
                              else if (i == 0) then
                                    neibor = i - 1
                              else
                                    neibor = i + 1
                                    if ((dz(i-1)+dz(i)) < (dz(i+1)+dz(i))) neibor = i-1
                              endif

                              if (neibor > i) then
                                    j = neibor;  l = i
                              else
                                    j = i;       l = neibor
                              endif
                              call combo_glc(dz(j), wliq(j), wice(j), tss(j), &
                                             dz(l), wliq(l), wice(l), tss(l))

                              if (j-1 > snl+1) then
                                    do k = j-1, snl+2, -1
                                          tss(k)  = tss(k-1)
                                          wice(k) = wice(k-1)
                                          wliq(k) = wliq(k-1)
                                          dz(k)   = dz(k-1)
                                    enddo
                              endif
                              snl = snl + 1
                              if (snl >= -1) EXIT
                        else
                              mssi = mssi + 1
                        endif
                  enddo
            end if

            ! Reset node depths and interface depths
            zi(0) = 0.
            do k = 0, snl+1, -1
                  z(k)    = zi(k) - 0.5*dz(k)
                  zi(k-1) = zi(k) - dz(k)
            enddo
      endif

      end subroutine snowlayerscombine_glc


      ! ======================================================================
      ! snowlayersdivide_glc
      ! Purpose: Subdivide snow layers when they exceed the maximum prescribed thickness.
      ! Authors: Shuwen Li and Xing Yuan
      ! ======================================================================
      subroutine snowlayersdivide_glc(lb, snl, z, dz, zi, wliq, wice, tss)

      USE shr_kind_mod, ONLY: r8 => shr_kind_rn
      implicit none

      integer,  INTENT(in)    :: lb
      integer,  INTENT(inout) :: snl
      real(r8), INTENT(inout) :: wice(lb:0), wliq(lb:0), tss(lb:0)
      real(r8), INTENT(inout) :: dz(lb:0), z(lb:0), zi(lb-1:0)

      real(r8) :: drr, dzsno(5), swice(5), swliq(5), tsno(5)
      integer  :: k, msno
      real(r8) :: zwice, zwliq, propor

      msno = abs(snl)
      do k = 1, msno
            dzsno(k) = dz  (k + snl)
            swice(k) = wice(k + snl)
            swliq(k) = wliq(k + snl)
            tsno(k)  = tss (k + snl)
      enddo

      if (msno == 1) then
            if (dzsno(1) > 0.03) then
                  msno = 2
                  dzsno(1) = dzsno(1)/2.;  swice(1) = swice(1)/2.;  swliq(1) = swliq(1)/2.
                  dzsno(2) = dzsno(1);     swice(2) = swice(1);     swliq(2) = swliq(1)
                  tsno(2)  = tsno(1)
            endif
      endif

      if (msno > 1) then
            if (dzsno(1) > 0.02) then
                  drr = dzsno(1) - 0.02
                  propor = drr/dzsno(1)
                  zwice  = propor*swice(1);  zwliq = propor*swliq(1)
                  propor = 0.02/dzsno(1)
                  swice(1) = propor*swice(1);  swliq(1) = propor*swliq(1);  dzsno(1) = 0.02
                  call combo_glc(dzsno(2), swliq(2), swice(2), tsno(2), drr, zwliq, zwice, tsno(1))
                  if (msno <= 2 .AND. dzsno(2) > 0.07) then
                        msno = 3
                        dzsno(2) = dzsno(2)/2.;  swice(2) = swice(2)/2.;  swliq(2) = swliq(2)/2.
                        dzsno(3) = dzsno(2);     swice(3) = swice(2);     swliq(3) = swliq(2)
                        tsno(3)  = tsno(2)
                  endif
            endif
      endif

      if (msno > 2) then
            if (dzsno(2) > 0.05) then
                  drr = dzsno(2) - 0.05
                  propor = drr/dzsno(2)
                  zwice  = propor*swice(2);  zwliq = propor*swliq(2)
                  propor = 0.05/dzsno(2)
                  swice(2) = propor*swice(2);  swliq(2) = propor*swliq(2);  dzsno(2) = 0.05
                  call combo_glc(dzsno(3), swliq(3), swice(3), tsno(3), drr, zwliq, zwice, tsno(2))
                  if (msno <= 3 .AND. dzsno(3) > 0.18) then
                        msno = 4
                        dzsno(3) = dzsno(3)/2.;  swice(3) = swice(3)/2.;  swliq(3) = swliq(3)/2.
                        dzsno(4) = dzsno(3);     swice(4) = swice(3);     swliq(4) = swliq(3)
                        tsno(4)  = tsno(3)
                  endif
            endif
      endif

      if (msno > 3) then
            if (dzsno(3) > 0.11) then
                  drr = dzsno(3) - 0.11
                  propor = drr/dzsno(3)
                  zwice  = propor*swice(3);  zwliq = propor*swliq(3)
                  propor = 0.11/dzsno(3)
                  swice(3) = propor*swice(3);  swliq(3) = propor*swliq(3);  dzsno(3) = 0.11
                  call combo_glc(dzsno(4), swliq(4), swice(4), tsno(4), drr, zwliq, zwice, tsno(3))
                  if (msno <= 4 .AND. dzsno(4) > 0.41) then
                        msno = 5
                        dzsno(4) = dzsno(4)/2.;  swice(4) = swice(4)/2.;  swliq(4) = swliq(4)/2.
                        dzsno(5) = dzsno(4);     swice(5) = swice(4);     swliq(5) = swliq(4)
                        tsno(5)  = tsno(4)
                  endif
            endif
      endif

      if (msno > 4) then
            if (dzsno(4) > 0.23) then
                  drr = dzsno(4) - 0.23
                  propor = drr/dzsno(4)
                  zwice  = propor*swice(4);  zwliq = propor*swliq(4)
                  propor = 0.23/dzsno(4)
                  swice(4) = propor*swice(4);  swliq(4) = propor*swliq(4);  dzsno(4) = 0.23
                  call combo_glc(dzsno(5), swliq(5), swice(5), tsno(5), drr, zwliq, zwice, tsno(4))
            endif
      endif

      snl = -msno
      do k = snl+1, 0
            dz(k)   = dzsno(k - snl)
            wice(k) = swice(k - snl)
            wliq(k) = swliq(k - snl)
            tss(k)  = tsno (k - snl)
      enddo

      zi(0) = 0.
      do k = 0, snl+1, -1
            z(k)    = zi(k) - 0.5*dz(k)
            zi(k-1) = zi(k) - dz(k)
      enddo

      end subroutine snowlayersdivide_glc


      ! ======================================================================
      ! combo_glc
      ! Purpose: Combine two snow/ice elements, conserving mass and enthalpy.
      ! Authors: Shuwen Li and Xing Yuan
      ! ======================================================================
      subroutine combo_glc(dz, wliq, wice, t, dz2, wliq2, wice2, t2)

      USE shr_kind_mod, ONLY: r8 => shr_kind_rn
      use module_sf_clm_constant, only : cpice, cpliq, hfus, tfrz
      implicit none

      real(r8), INTENT(in)    :: dz2, wliq2, wice2, t2
      real(r8), INTENT(inout) :: dz, wliq, wice, t

      real(r8) :: dzc, wliqc, wicec, tc, h, h2, hc

      dzc   = dz + dz2
      wicec = wice + wice2
      wliqc = wliq + wliq2
      h     = (cpice*wice  + cpliq*wliq )*(t  - tfrz) + hfus*wliq
      h2    = (cpice*wice2 + cpliq*wliq2)*(t2 - tfrz) + hfus*wliq2
      hc    = h + h2

      if (hc < 0.) then
            tc = tfrz + hc/(cpice*wicec + cpliq*wliqc)
      else if (hc <= hfus*wliqc) then
            tc = tfrz
      else
            tc = tfrz + (hc - hfus*wliqc)/(cpice*wicec + cpliq*wliqc)
      endif

      dz   = dzc
      wice = wicec
      wliq = wliqc
      t    = tc

      end subroutine combo_glc


      ! ======================================================================
      ! winddriftcompaction_glc
      ! Purpose: Compute wind drift compaction for a single snow layer.
      ! Must be called top-to-bottom order to correctly update zpseudo and mobile.
      ! Authors: Shuwen Li and Xing Yuan
      ! ======================================================================
      subroutine winddriftcompaction_glc(bi, forc_wind, dz, zpseudo, mobile, compaction_rate)

      USE shr_kind_mod, ONLY: r8 => shr_kind_rn

      real(r8), intent(in)    :: bi, forc_wind, dz
      real(r8), intent(inout) :: zpseudo
      logical,  intent(inout) :: mobile
      real(r8), intent(out)   :: compaction_rate

      real(r8) :: Frho, MO, SI, gamma_drift, tau_inverse

      real(r8), parameter :: rho_min    = 50._r8
      real(r8), parameter :: rho_max    = 350._r8
      real(r8), parameter :: drift_gs   = 0.35e-3_r8
      real(r8), parameter :: drift_sph  = 1.0_r8
      real(r8), parameter :: tau_ref    = 48._r8 * 3600._r8

      if (mobile) then
            Frho = 1.25_r8 - 0.0042_r8*(max(rho_min, bi) - rho_min)
            ! Non-dendritic snow: dendricity=0, sphericity=1, grain size=0.35 mm
            MO = 0.34_r8 * (-0.583_r8*drift_gs - 0.833_r8*drift_sph + 0.833_r8) + 0.66_r8*Frho
            SI = -2.868_r8 * exp(-0.085_r8*forc_wind) + 1._r8 + MO

            if (SI > 0.0_r8) then
                  SI = min(SI, 3.25_r8)
                  zpseudo = zpseudo + 0.5_r8 * dz * (3.25_r8 - SI)
                  gamma_drift   = SI*exp(-zpseudo/0.1_r8)
                  tau_inverse   = gamma_drift / tau_ref
                  compaction_rate = -max(0.0_r8, rho_max - bi) * tau_inverse
                  zpseudo = zpseudo + 0.5_r8 * dz * (3.25_r8 - SI)
            else
                  mobile          = .false.
                  compaction_rate = 0._r8
            end if
      else
            compaction_rate = 0._r8
      end if

      end subroutine winddriftcompaction_glc


      ! ======================================================================
      ! qsadv_glc
      ! Purpose: Compute saturation mixing ratio and its temperature derivative.
      ! Authors: Shuwen Li and Xing Yuan
      ! ======================================================================
      subroutine qsadv_glc(T, p, es, esdT, qs, qsdT)

      USE shr_kind_mod, ONLY: r8 => shr_kind_rn
      implicit none

      real(r8), INTENT(in)  :: T     ! Temperature [K]
      real(r8), INTENT(in)  :: p     ! Surface pressure [Pa]
      real(r8), INTENT(out) :: es    ! Vapor pressure [Pa]
      real(r8), INTENT(out) :: esdT  ! d(es)/d(T) [Pa/K]
      real(r8), INTENT(out) :: qs    ! Specific humidity [kg/kg]
      real(r8), INTENT(out) :: qsdT  ! d(qs)/d(T) [1/K]

      real(r8) :: td, vp, vp1, vp2
      real(r8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8
      real(r8) :: b0,b1,b2,b3,b4,b5,b6,b7,b8
      real(r8) :: c0,c1,c2,c3,c4,c5,c6,c7,c8
      real(r8) :: d0,d1,d2,d3,d4,d5,d6,d7,d8

      ! Water vapor coefficients (0°C to 100°C)
      data a0/6.11213476   /,a1/ 0.444007856    /,a2/0.143064234e-01/  &
          ,a3/0.264461437e-03/,a4/ 0.305903558e-05/,a5/0.196237241e-07/  &
          ,a6/0.892344772e-10/,a7/-0.373208410e-12/,a8/0.209339997e-15/

      ! Water vapor derivative coefficients
      data b0/0.444017302  /,b1/ 0.286064092e-01/,b2/ 0.794683137e-03/ &
          ,b3/ 0.121211669e-04/,b4/ 0.103354611e-06/,b5/ 0.404125005e-09/ &
          ,b6/-0.788037859e-12/,b7/-0.114596802e-13/,b8/ 0.381294516e-16/

      ! Ice coefficients (-75°C to 0°C)
      data c0/6.11123516     /,c1/0.503109514    /,c2/0.188369801e-01/ &
          ,c3/0.420547422e-03/,c4/0.614396778e-05/,c5/0.602780717e-07/ &
          ,c6/0.387940929e-09/,c7/0.149436277e-11/,c8/0.262655803e-14/

      ! Ice derivative coefficients
      data d0/0.503277922    /,d1/0.377289173e-01/,d2/0.126801703e-02/ &
          ,d3/0.249468427e-04/,d4/0.313703411e-06/,d5/0.257180651e-08/ &
          ,d6/0.133268878e-10/,d7/0.394116744e-13/,d8/0.498070196e-16/

      td = T - 273.16

      if (td < -75.0 .or. td > 100.0) print *, "qsadv: abnormal temperature", T
      if (td < -75.0) td = -75.0
      if (td > 100.0) td = 100.0

      if (td >= 0.0) then
            es   = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
                 + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
            esdT = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
                 + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
      else
            es   = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
                 + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
            esdT = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
                 + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
      end if

      es   = es   * 100.
      esdT = esdT * 100.
      vp   = 1.0  / (p - 0.378*es)
      vp1  = 0.622 * vp
      vp2  = vp1  * vp
      qs   = es   * vp1
      qsdT = esdT * vp2 * p

      end subroutine qsadv_glc


      ! ======================================================================
      ! moninobukini_glc
      ! Purpose: Initialize Monin-Obukhov length.
      ! Authors: Shuwen Li and Xing Yuan
      ! ======================================================================
      subroutine moninobukini_glc(ur, th, thm, thv, dth, dqh, dthv, zldis, z0m, um, obu)

      USE shr_kind_mod, ONLY: r8 => shr_kind_rn
      use module_sf_clm_constant, only : grav, vonkar
      implicit none

      real(r8), INTENT(in)  :: ur, thm, th, thv, dth, dthv, dqh, zldis, z0m
      real(r8), INTENT(out) :: um, obu

      real(r8) :: wc, rib, zeta

      wc = 0.5
      if (dthv >= 0.) then
            um = max(ur, 0.1)
      else
            um = sqrt(ur*ur + wc*wc)
      endif

      rib = grav*zldis*dthv/(thv*um*um)

      if (rib >= 0.) then
            zeta = rib*log(zldis/z0m)/(1. - 5.*min(rib, 0.19))
            zeta = min(2., max(zeta, 1.e-6))
      else
            zeta = rib*log(zldis/z0m)
            zeta = max(-100., min(zeta, -1.e-6))
      endif
      obu = zldis/zeta

      end subroutine moninobukini_glc


      ! ======================================================================
      ! moninobuk_glc
      ! Purpose: Compute friction velocity and profile relations for potential temperature
      ! and humidity in the surface boundary layer.
      ! Authors: Shuwen Li and Xing Yuan
      ! ======================================================================
      subroutine moninobuk_glc(hu, ht, hq, displa, z0m, z0h, z0q, obu, um, &
                               ustar, fh2m, fq2m, fm10m, fm, fh, fq)

      USE shr_kind_mod, ONLY: r8 => shr_kind_rn
      use module_sf_clm_constant, only : vonkar
      implicit none

      real(r8), INTENT(in)  :: hu, ht, hq, displa, z0m, z0h, z0q, obu, um
      real(r8), INTENT(out) :: ustar, fh2m, fq2m, fm10m, fm, fh, fq

      real(r8) :: zldis, zetam, zetat, zeta

      ! Momentum profile
      zldis = hu - displa
      zeta  = zldis/obu
      zetam = 1.574
      if (zeta < -zetam) then
            fm    = log(-zetam*obu/z0m) - psi(1,-zetam) &
                  + psi(1,z0m/obu) + 1.14*((-zeta)**0.333 - (zetam)**0.333)
            ustar = vonkar*um/fm
      else if (zeta < 0.) then
            fm    = log(zldis/z0m) - psi(1,zeta) + psi(1,z0m/obu)
            ustar = vonkar*um/fm
      else if (zeta <= 1.) then
            fm    = log(zldis/z0m) + 5.*zeta - 5.*z0m/obu
            ustar = vonkar*um/fm
      else
            fm    = log(obu/z0m) + 5. - 5.*z0m/obu + (5.*log(zeta)+zeta-1.)
            ustar = vonkar*um/fm
      endif

      ! 10m wind profile
      zldis = 10. + z0m
      zeta  = zldis/obu
      if (zeta < -zetam) then
            fm10m = log(-zetam*obu/z0m) - psi(1,-zetam) &
                  + psi(1,z0m/obu) + 1.14*((-zeta)**0.333 - (zetam)**0.333)
      else if (zeta < 0.) then
            fm10m = log(zldis/z0m) - psi(1,zeta) + psi(1,z0m/obu)
      else if (zeta <= 1.) then
            fm10m = log(zldis/z0m) + 5.*zeta - 5.*z0m/obu
      else
            fm10m = log(obu/z0m) + 5. - 5.*z0m/obu + (5.*log(zeta)+zeta-1.)
      endif

      ! Temperature profile
      zldis = ht - displa
      zeta  = zldis/obu
      zetat = 0.465
      if (zeta < -zetat) then
            fh = log(-zetat*obu/z0h) - psi(2,-zetat) &
               + psi(2,z0h/obu) + 0.8*((zetat)**(-0.333) - (-zeta)**(-0.333))
      else if (zeta < 0.) then
            fh = log(zldis/z0h) - psi(2,zeta) + psi(2,z0h/obu)
      else if (zeta <= 1.) then
            fh = log(zldis/z0h) + 5.*zeta - 5.*z0h/obu
      else
            fh = log(obu/z0h) + 5. - 5.*z0h/obu + (5.*log(zeta)+zeta-1.)
      endif

      ! 2m temperature profile
      zldis = 2. + z0h
      zeta  = zldis/obu
      if (zeta < -zetat) then
            fh2m = log(-zetat*obu/z0h) - psi(2,-zetat) &
                 + psi(2,z0h/obu) + 0.8*((zetat)**(-0.333) - (-zeta)**(-0.333))
      else if (zeta < 0.) then
            fh2m = log(zldis/z0h) - psi(2,zeta) + psi(2,z0h/obu)
      else if (zeta <= 1.) then
            fh2m = log(zldis/z0h) + 5.*zeta - 5.*z0h/obu
      else
            fh2m = log(obu/z0h) + 5. - 5.*z0h/obu + (5.*log(zeta)+zeta-1.)
      endif

      ! Humidity profile
      zldis = hq - displa
      zeta  = zldis/obu
      if (zeta < -zetat) then
            fq = log(-zetat*obu/z0q) - psi(2,-zetat) &
               + psi(2,z0q/obu) + 0.8*((zetat)**(-0.333) - (-zeta)**(-0.333))
      else if (zeta < 0.) then
            fq = log(zldis/z0q) - psi(2,zeta) + psi(2,z0q/obu)
      else if (zeta <= 1.) then
            fq = log(zldis/z0q) + 5.*zeta - 5.*z0q/obu
      else
            fq = log(obu/z0q) + 5. - 5.*z0q/obu + (5.*log(zeta)+zeta-1.)
      endif

      ! 2m humidity profile
      zldis = 2. + z0h
      zeta  = zldis/obu
      if (zeta < -zetat) then
            fq2m = log(-zetat*obu/z0q) - psi(2,-zetat) &
                 + psi(2,z0q/obu) + 0.8*((zetat)**(-0.333) - (-zeta)**(-0.333))
      else if (zeta < 0.) then
            fq2m = log(zldis/z0q) - psi(2,zeta) + psi(2,z0q/obu)
      else if (zeta <= 1.) then
            fq2m = log(zldis/z0q) + 5.*zeta - 5.*z0q/obu
      else
            fq2m = log(obu/z0q) + 5. - 5.*z0q/obu + (5.*log(zeta)+zeta-1.)
      endif

      end subroutine moninobuk_glc


      ! ======================================================================
      ! rain_snow_temp
      ! Purpose: Determine precipitation temperature and new snow bulk density.
      ! Authors: Shuwen Li and Xing Yuan
      ! ======================================================================
      SUBROUTINE rain_snow_temp(scwater, tm, qm, psrf, us, vs, tcrit, t_precip, bifall)

      USE shr_kind_mod, ONLY: r8 => shr_kind_rn
      use module_sf_clm_constant, only : tfrz
      IMPLICIT NONE

      integer,  INTENT(in)  :: scwater
      real(r8), INTENT(in)  :: tm, qm, psrf, us, vs, tcrit
      real(r8), INTENT(out) :: t_precip, bifall

      ! Compute wet-bulb temperature as precipitation temperature
      call wetbulb(tm, psrf, qm, t_precip)

      ! New snow bulk density
      CALL NewSnowBulkDensity(tm, us, vs, bifall)

      ! Constrain precipitation temperature
      if (tm > 275.65) then
            if (t_precip < tfrz) t_precip = tfrz
      else
            t_precip = min(tfrz, t_precip)
      endif

      END SUBROUTINE rain_snow_temp


      ! ======================================================================
      ! NewSnowBulkDensity
      ! Purpose: Compute bulk density of newly fallen dry snow, including wind packing.
      ! Authors: Shuwen Li and Xing Yuan
      ! ======================================================================
      SUBROUTINE NewSnowBulkDensity(tm, us, vs, bifall)

      USE shr_kind_mod, ONLY: r8 => shr_kind_rn
      use module_sf_clm_constant, only : tfrz
      implicit none

      real(r8), INTENT(in)  :: tm, us, vs
      real(r8), INTENT(out) :: bifall

      real(r8) :: t_for_bifall_degC, forc_wind

      if (tm > tfrz + 2.0) then
            bifall = 50.0 + 1.7*(17.0)**1.5
      else if (tm > tfrz - 15.0) then
            bifall = 50.0 + 1.7*(tm - tfrz + 15.0)**1.5
      else
            ! Below ~-15°C: smaller crystals yield higher density
            if (tm > tfrz - 57.55) then
                  t_for_bifall_degC = (tm - tfrz)
            else
                  t_for_bifall_degC = -57.55  ! Lower bound to prevent function inversion
            end if
            bifall = -(50.0/15.0 + 0.0333*15.0)*t_for_bifall_degC - 0.0333*t_for_bifall_degC**2
      end if

      forc_wind = sqrt(us**2 + vs**2)
      if (forc_wind > 0.1) then
            bifall = bifall + (266.861 * ((1.0 + TANH(forc_wind/5.0))/2.0)**8.8)
      end if

      END SUBROUTINE NewSnowBulkDensity


      ! ======================================================================
      ! wetbulb
      ! Purpose: Compute wet-bulb temperature by iterative psychrometric equation.
      ! Author: Xing Yuan
      ! ======================================================================
      subroutine wetbulb(t, p, q, twc)

      USE shr_kind_mod, ONLY: r8 => shr_kind_rn
      use module_sf_clm_constant, only : tfrz, hvap, cpair
      implicit none

      real(r8), intent(in)  :: t    ! Air temperature [K]
      real(r8), intent(in)  :: p    ! Pressure [Pa]
      real(r8), intent(in)  :: q    ! Specific humidity [kg/kg]
      real(r8), intent(out) :: twc  ! Wet-bulb temperature [K]

      integer  :: i
      real(r8) :: es, esdT, qs, qsdT, r, rws

      call qsadv_glc(t, p, es, esdT, qs, qsdT)
      r = q/(1.0 - q)
      if (q >= qs) r = qs/(1.0 - qs)
      twc = t
      do i = 1, 6
            call qsadv_glc(twc, p, es, esdT, qs, qsdT)
            rws = qs/(1.0 - qs)
            twc = (twc + t + hvap/cpair*(r - rws))/2.0
      enddo

      end subroutine wetbulb


      ! ======================================================================
      ! tridia
      ! Purpose: Solve tridiagonal system by Thomas algorithm.
      ! Falls back to mass-lumped method on singular diagonal.
      ! Author: Xing Yuan, 2009-10-22
      ! ======================================================================
      SUBROUTINE tridia(n, a, b, c, r, u)

      integer,  intent(in)    :: n
      real(r8), intent(in)    :: a(1:n), b(1:n), c(1:n), r(1:n)
      real(r8), intent(inout) :: u(1:n)

      real(r8) :: gam(1:n), bet, u1(1:n)
      integer  :: j

      bet = b(1)
      if (abs(bet) < 1.e-15) then
            call tridiaml(n, a, b, c, r, u)
            return
      endif
      u1(1) = r(1) / bet
      do j = 2, n
            gam(j) = c(j-1) / bet
            bet     = b(j) - a(j) * gam(j)
            if (abs(bet) < 1.e-15) then
                  call tridiaml(n, a, b, c, r, u)
                  return
            endif
            u1(j) = (r(j) - a(j)*u1(j-1)) / bet
      enddo
      do j = n-1, 1, -1
            u1(j) = u1(j) - gam(j+1) * u1(j+1)
      enddo
      u = u1

      END SUBROUTINE tridia


      ! ======================================================================
      ! tridiaml
      ! Solve tridiagonal system by mass-lumped method (fallback for singularity).
      ! Author: Xing Yuan, 2009-10-22
      ! ======================================================================
      SUBROUTINE tridiaml(n, a, b, c, r, u)

      integer,  intent(in)    :: n
      real(r8), intent(in)    :: a(1:n), b(1:n), c(1:n), r(1:n)
      real(r8), intent(inout) :: u(1:n)

      real(r8) :: b1(1:n)
      integer  :: j

      b1 = a + b + c
      do j = 1, n
            if (abs(b1(j)) > 0.) u(j) = r(j)/b1(j)
      enddo

      END SUBROUTINE tridiaml


      ! ======================================================================
      ! psi
      ! Stability function for the unstable case (rib < 0).
      !   k=1: momentum
      !   k=2: heat/moisture
      ! ======================================================================
      real(r8) function psi(k, zeta)

      USE shr_kind_mod, ONLY: r8 => shr_kind_rn
      implicit none

      integer  :: k
      real(r8) :: zeta, chik

      chik = (1. - 16.*zeta)**0.25
      if (k == 1) then
            psi = 2.*log((1.+chik)*0.5) + log((1.+chik*chik)*0.5) - 2.*atan(chik) + 2.*atan(1.)
      else
            psi = 2.*log((1.+chik*chik)*0.5)
      endif

      end function psi


      ! ======================================================================
      ! calc_dewpoint
      ! Purpose: Calculate dew-point temperature from specific humidity and pressure
      ! using the Tetens formula.
      ! Authors: Xing Yuan and Shuwen Li, 2024-06-01
      ! Refrence:
      ! Li, S., and X. Yuan (2026): Evaluation of a Land Surface-Glacier
!     ! Coupled Model over the Three-River Headwaters Region. Water, 18, 1030.
      ! ======================================================================
      subroutine calc_dewpoint(qs, p, T_dew)

      real(r8), INTENT(in)  :: qs     ! Specific humidity [kg/kg]
      real(r8), INTENT(in)  :: p      ! Atmospheric pressure [Pa]
      real(r8), INTENT(out) :: T_dew  ! Dew-point temperature [degrees Celsius]

      real(r8) :: e, a, b, c, ln_e_a

      a = 6.11;  b = 17.27;  c = 237.7

  
      e = qs * p

      ln_e_a = log(e / a)
      T_dew  = (b * ln_e_a) / (c - ln_e_a)

      end subroutine calc_dewpoint


      ! ======================================================================
      ! meltf_glc (legacy phase change subroutine, kept for reference)
      ! Original author: Nan Wei, 04/2023
      ! Revisions: Nan Wei, 04/2023 (supercooled water option)
      ! Transfered by Shuwen Li
      ! Reference:
      ! [1] Li, S., and X. Yuan (2026): Evaluation of a Land Surface-Glacier
!     !     Coupled Model over the Three-River Headwaters Region. Water, 18, 1030.
      ! [2] Dai, Y., Wei, N., Yuan, H., et al., 2019: Evaluation of soil thermal conductivity 
      !     schemes for use in land surface modeling. JAMES, 11. 
      ! ======================================================================
      subroutine meltf_glc(scwater, lb, nl_ice, nl_snow, dtime, tm, ht, &
                  fact, brr, hs, dhsdT, tg, &
                  tss_bef, tss, wliq, wice, imelt, &
                  scv, snowdp, sm, xmf, melt, smelt, refreeze, dz, ipr, jpr)

      USE shr_kind_mod, ONLY: r8 => shr_kind_rn
      use module_sf_clm_constant, only : tfrz, hfus, grav

      IMPLICIT NONE

      integer,  INTENT(in)    :: scwater, nl_ice, nl_snow, lb
      integer                 :: ipr, jpr
      real(r8), INTENT(in)    :: dtime, tss_bef(lb:nl_ice), brr(lb:nl_ice), fact(lb:nl_ice)
      real(r8), INTENT(inout) :: hs
      real(r8), INTENT(in)    :: dhsdT, dz(1:nl_ice), ht, tm, tg
      real(r8), INTENT(inout) :: tss(lb:nl_ice), wice(lb:nl_ice), wliq(lb:nl_ice)
      real(r8), INTENT(inout) :: scv, snowdp
      real(r8), INTENT(out)   :: sm, xmf, melt, smelt, refreeze
      integer,  INTENT(out)   :: imelt(nl_snow+1:nl_ice)

      ! Local variables
      real(r8) :: hm(lb:nl_ice), xm(lb:nl_ice)
      real(r8) :: heatr, temp1, propor, tinc
      real(r8) :: wmass0(lb:nl_ice), wice0(lb:nl_ice), wliq0(lb:nl_ice)
      real(r8) :: scvold, we, sm_heat
      integer  :: j
      real(r8) :: dt_glc, excess_heat

      sm       = 0.
      xmf      = 0.
      sm_heat  = 0.
      refreeze = 0.
      do j = lb, nl_ice
            imelt(j)  = 0
            hm(j)     = 0.
            xm(j)     = 0.
            wice0(j)  = wice(j)
            wliq0(j)  = wliq(j)
            wmass0(j) = wice(j) + wliq(j)
      enddo

      scvold = scv
      we     = 0.
      if (lb <= 0) we = sum(wice(lb:nl_ice) + wliq(lb:nl_ice))

      ! Phase change identification
      do j = lb, nl_ice
            if (tss(j) > tfrz) then
                
                  ! cpice [J/(kg·K)] requires mass [kg/m2] to yield energy [J/m2]
                  dt_glc      = tg - tfrz
                  tss(j)      = tfrz
                  excess_heat = dt_glc*cpice
                  hs          = hs + excess_heat
                  imelt(j)    = 1
            endif
            if (j <= 0) then
                  if (wliq(j) > 0. .and. tss(j) < tfrz) then
                        imelt(j) = 2;  tss(j) = tfrz
                  endif
            else
                  if (wliq(j) > 0. .and. tss(j) < tfrz) then
                        imelt(j) = 2;  tss(j) = tfrz
                  endif
            end if
      enddo

      if (lb == 1 .and. scv > 0.) then
            if (tss(lb) > tfrz) then
                  imelt(1) = 1;  tss(1) = tfrz
            endif
      endif

      ! Energy residual for phase change
      do j = lb, nl_ice
            if (imelt(j) > 0) then
                  tinc = tss(j) - tss_bef(j)
                  if (j > lb) then
                        hm(j) = brr(j) - tinc/fact(j)
                  else
                        hm(j) = hs + dhsdT*tinc + brr(j) - tinc/fact(j)
                  endif
            endif
      enddo

      do j = lb, nl_ice
            if (imelt(j) == 1 .and. hm(j) < 0.) then
                  hm(j) = 0.;  imelt(j) = 0
            endif
      enddo

      heatr = 0.
      smelt = 0.

      do j = lb, nl_ice
            if (imelt(j) > 0 .and. abs(hm(j)) > 0.) then
                  xm(j) = hm(j)*dtime/hfus

                  ! Surface snow melt
                  if (j == lb .and. scv > 0. .and. xm(j) > 0.) then
                        temp1  = scv
                        scv    = max(0., temp1 - xm(j))
                        propor = scv/temp1
                        snowdp = propor * snowdp
                        heatr  = hm(j) - hfus*(temp1 - scv)/dtime
                        if (heatr > 0.) then
                              xm(j) = heatr*dtime/hfus
                              hm(j) = heatr
                        else
                              xm(j) = 0.;  hm(j) = 0.
                        endif
                        sm      = max(0., (temp1 - scv))/dtime
                        sm_heat = sm*dtime*hfus
                        xmf     = hfus*sm
                  endif

                  if (xm(j) > 1e-6) then
                        wice(j) = max(0., wice0(j) - xm(j))
                        heatr   = hm(j) - hfus*xm(j)/dtime
                  else
                        wice(j) = min(wmass0(j), wice0(j) - xm(j))
                  endif

                  wliq(j) = max(0., wmass0(j) - wice(j))

                  if (abs(heatr) > 0.) then
                        if (j > lb) then
                              tss(j) = tss(j) + fact(j)*heatr
                        else
                              tss(j) = tss(j) + fact(j)*heatr/(1. - fact(j)*dhsdT)
                        endif
                  endif
                  xmf = xmf + hfus * (wice0(j) - wice(j))/dtime
            endif
      enddo

      if (tg > tfrz - 2. .and. hs > 0.) smelt = max(0., (hs - sm_heat))*dtime/hfus

      melt = xmf*dtime/hfus

      if (lb <= 0) then
            we = sum(wice(lb:nl_ice) + wliq(lb:nl_ice)) - we
            if (abs(we) > 1.e-4) then
                  print*, 'meltf mass conservation error:', we, ipr, jpr, &
                          sum(wice(lb:nl_ice) + wliq(lb:nl_ice))
            endif
      endif

      end subroutine meltf_glc
! ----------------------------------------------------------------------------

!----------------------------------------------------------------------------
      END MODULE MODULE_GLACIER