#include "configs.inc"

! start of irrigation module
#ifdef IrrDmd

module module_irrigation
!=======================================================================
!
!   Purpose: Simulate irrigation impact on hydrological cycle
!   Authors: Chenyuan Li, Yang Jiao, Xing Yuan, 2023-10-04
!   Reference:
!   [1]  X. Yuan*., C. Li, P. Ji, XZ. Xi, SW. Li, Y. Jiao, H. Yi, XZ. Liang, QC. Zeng, 2026:  
!        Development and evaluation of Conjunctive Surface-Subsurface Process Model version 3
!        (CSSPv3) at global scale. Journal of Advances in Modeling Earth Systems. Revised.
!   [2]  Li, C., X. Yuan*, Y. Jiao, P.Ji, Z.Huang,2024: High-resolution Land Surface Modeling 
!        of the Irrigation Effects on Evapotranspiration over the Yellow River Basin. 
!        Journal of Hydrology, 633, 130986. 
!=======================================================================
  use shr_kind_mod, ONLY: r8 => shr_kind_rn
  implicit none
  integer :: NRESSTEP  ! operation per N timesteps (=24h)
  integer :: resstep = 0
  !SW2GW = False. Assume the defcit surface water is satisfied by other origins
  !SW2GW = True. Satisfy the deficit of surface water by groundwater to keep water balance
  logical :: SW2GW = .True.
  public :: NRESSTEP   ! nlocres
!   Public member functions:
  public :: cssp_irrg
  public :: cssp_irrg_local_swpump
#ifdef Swpump
  public :: cssp_irrg_region_swpump
#endif
!=======================================================================
contains
  subroutine cssp_irrg(xnl_soil_br,SC_WATER,nl_soil,nl_snow,dtime,xarea, &
                      xdz,xhksm,xbsw,xphi0,xporsl,xwdry,xwliq,xwice, &
                      xirrgpadfrac,xirrgnopfrac,xirrgswfrac,xirrgthresh,&
                      xirrgwtrdmd,xswpumprate,xgwpumprate, &
                      xirrgwtrdmd_paddy,xirrgwtrdmd_nonpaddy,xpumprate_paddy, &
                      im,id,time_zone,&
#ifdef mask
                      calcumask,                 &
#endif
                      ids,ide ,jds,jde ,kds,kde ,&
                      ims,ime ,jms,jme ,kms,kme ,&
                      its,ite ,jts,jte ,kts,kte )
!----------------------- Dummy argument -------------------------------------
  integer, intent(in) :: &
      ids,ide,     &! start/end index for i in domain
      jds,jde,     &! start/end index for j in domain
      kds,kde,     &! start/end index for k in domain
      ims,ime,     &! start/end index for i in memory
      jms,jme,     &! start/end index for j in memory
      kms,kme,     &! start/end index for k in memory
      its,ite,     &! start/end index for i in tile
      jts,jte,     &! start/end index for j in tile
      kts,kte,     &! start/end index for k in tile
      nl_soil,     &! number of soil layers (>0)
      nl_snow       ! max number of snow layers (<0)
  real, intent(in) :: &
      dtime         ! CLM time step [seconds]
#ifdef LonLat
  real, intent(in), dimension(ims:ime,jms:jme) :: &
      xarea          ! grid cell area [m2]
#else
  real, intent(in) :: &
      xarea          ! grid cell area [m2]
#endif
!----------------------- Need Halo 8/24 points ------------------------------
  real, intent(in), dimension(ims:ime,jms:jme) :: &                         ! halo 8
      SC_WATER      ! surface category of water characteristics
                    ! 1=urban and built-up, 2=soil, 3=wetland, 4=land ice
                    ! 5=shallow lake, 6=deep lake,  7=sea ice, 8=ocean
   
  integer, intent(in), dimension(ims:ime,jms:jme) :: &                      ! halo 8
      xnl_soil_br   ! BEDROCK depth  !jpeng
  real, intent(in), dimension(ims:ime,jms:jme) ::    &                      ! halo 8
#ifdef mask
      calcumask     ! mask
#endif
  real, intent(in), dimension(ims:ime,1:nl_soil,jms:jme) :: &               ! halo 8
      xbsw,        &! Clapp-Hornberger (b)
      xphi0,       &! soil minimum suction [mm]
      xporsl,      &! soil porosity (fraction that is voids) [-]
      xwdry,       &! soil dry or minimum residual water content [-]
      xhksm         ! soil hydraulic conductivity at saturation averaged in a layer [mm h2o/s]

  real, intent(in), dimension(ims:ime,1:nl_soil-nl_snow,jms:jme) :: &       ! halo 8
      xdz           ! land layer thickness [m]

  real, intent(in), dimension(ims:ime,1:nl_soil-nl_snow,jms:jme) :: &       ! halo 8
      xwice         ! land layer ice lens [kg/m2]

  real, intent(inout), dimension(ims:ime,1:nl_soil-nl_snow,jms:jme) :: &    ! halo 8
      xwliq         ! land layer liquid water [kg/m2]

!----------------------- Irrigation Values ------------------------------
  real(r8), intent(in), dimension(ims:ime,jms:jme,1:12) :: &
      xirrgpadfrac, &! irrigated paddy area ratio
      xirrgnopfrac   ! irrigated non-paddy area ratio [-]
  real(r8), intent(in), dimension(ims:ime,jms:jme) :: &
      xirrgswfrac,  &! surface water irrigated area in the area equipped for irrigation ratio [-]
      xirrgthresh    ! the weighting factor which can control the irrigation target and amount 
  real(r8), intent(inout), dimension(ims:ime,jms:jme) :: &
      xirrgwtrdmd_paddy, & ! water demand for paddy irrigation [10^6 m3]
      xirrgwtrdmd_nonpaddy,  & ! water demand for non-paddy irrigaion [10^6m3]
      xirrgwtrdmd,  &! water demand for irrgation [10^6 m3]
      xswpumprate,  &! surface water pumping/irrigate rate [mm h2o/s]
      xgwpumprate,  &! ground water pumping/irrigate rate [mm h2o/s] (pump groundwater to fulfill water demand which was not satisfied)
      xpumprate_paddy! total water pumping/irrigate rate for paddy [mm h2o/s] 
  integer, intent(in) :: im, id   ! month and days in month
  integer, intent(in), dimension(ims:ime,jms:jme) :: time_zone ! natural time zone

! calculate irrigation amount
    call cssp_irrg_waterdemand(  &
         xnl_soil_br,SC_WATER,nl_soil,nl_snow,dtime,xarea, &
         xdz,xhksm,xbsw,xphi0,xporsl,xwdry,xwliq,xwice, &
         xirrgpadfrac,xirrgnopfrac,xirrgthresh,xirrgwtrdmd,im,id, &
         xirrgwtrdmd_paddy,xirrgwtrdmd_nonpaddy,time_zone,&  !paddy
#ifdef mask
         calcumask,                 &
#endif
         ids,ide ,jds,jde ,kds,kde ,&
         ims,ime ,jms,jme ,kms,kme ,&
         its,ite ,jts,jte ,kts,kte )

! Allocate water sources based on census ratio
    call cssp_irrg_watersupply(dtime,xarea,time_zone,&
            xirrgwtrdmd,xirrgswfrac,xswpumprate,xgwpumprate,&
            xirrgwtrdmd_paddy,xirrgwtrdmd_nonpaddy,xpumprate_paddy,&
            ims,ime,jms,jme,its,ite,jts,jte)
    ! update reservoir timestep counter
    resstep = resstep + 1
    ! reset counter, in case resstep is too large (>2^31-1)
    if (mod(resstep, NRESSTEP) == 0) resstep = 0
  end subroutine cssp_irrg

  subroutine cssp_irrg_waterdemand(  &
             xnl_soil_br,SC_WATER,nl_soil,nl_snow,dtime,xarea, &
             xdz,xhksm,xbsw,xphi0,xporsl,xwdry,xwliq,xwice, &
             xirrgpadfrac,xirrgnopfrac,xirrgthresh,xirrgwtrdmd,im,id, &
             xirrgwtrdmd_paddy,xirrgwtrdmd_nonpaddy,time_zone,&  !paddy
#ifdef mask
             calcumask,                 &
#endif
             ids,ide ,jds,jde ,kds,kde ,&
             ims,ime ,jms,jme ,kms,kme ,&
             its,ite ,jts,jte ,kts,kte )
!=======================================================================
!   Purpose: calculate irrigation amount
!   Authors:  Chenyuan Li, Yang Jiao
!   Revised for consideration the weight factor, 2022-10
!   [1]  Li, C., X. Yuan*, Y. Jiao, P.Ji, Z.Huang, 2024: High-resolution Land Surface Modeling 
!        of the Irrigation Effects on Evapotranspiration over the Yellow River Basin. 
!        Journal of Hydrology, 633, 130986. 
!   [2]  X. Yuan*., C. Li, P. Ji, XZ. Xi, SW. Li, Y. Jiao, H. Yi, XZ. Liang, QC. Zeng, 2026:  
!        Development and evaluation of Conjunctive Surface-Subsurface Process Model version 3
!        (CSSPv3) at global scale. Journal of Advances in Modeling Earth Systems. Revised.
!=======================================================================
    USE module_sf_clm_constant, ONLY: denliq, denice, smpmax
    implicit none
    integer, intent(in) :: &
        ids,ide,     &! start/end index for i in domain
        jds,jde,     &! start/end index for j in domain
        kds,kde,     &! start/end index for k in domain
        ims,ime,     &! start/end index for i in memory
        jms,jme,     &! start/end index for j in memory
        kms,kme,     &! start/end index for k in memory
        its,ite,     &! start/end index for i in tile
        jts,jte,     &! start/end index for j in tile
        kts,kte,     &! start/end index for k in tile
        nl_soil,     &! number of soil layers (>0)
        nl_snow       ! max number of snow layers (<0)

    real, intent(in) :: &
        dtime         ! CLM time step [seconds]
#ifdef LonLat
    real, intent(in), dimension(ims:ime,jms:jme) :: &
        xarea          ! grid cell area [m2]
#else
    real, intent(in) :: &
        xarea          ! grid cell area [m2]
#endif

    real, intent(in), dimension(ims:ime,jms:jme) :: &
        SC_WATER      ! surface category of water characteristics
                      ! 1=urban and built-up, 2=soil, 3=wetland, 4=land ice
                      ! 5=shallow lake, 6=deep lake,  7=sea ice, 8=ocean

    integer, intent(in), dimension(ims:ime,jms:jme) :: &
        xnl_soil_br   ! BEDROCK depth  !jpeng

#ifdef mask
    real, intent(in), dimension(ims:ime,jms:jme) :: &
        calcumask     ! mask
#endif
    real, intent(in), dimension(ims:ime,1:nl_soil,jms:jme) :: &
        xbsw,        &! Clapp-Hornberger (b)
        xphi0,       &! soil minimum suction [mm]
        xporsl,      &! soil porosity (fraction that is voids) [-]
        xwdry,       &! soil dry or minimum residual water content [-]
        xhksm         ! soil hydraulic conductivity at saturation averaged in a layer [mm h2o/s]

    real, intent(in), dimension(ims:ime,1:nl_soil-nl_snow,jms:jme) :: &
        xdz,         &! land layer thickness [m]
        xwice,       &! land layer ice lens [kg/m2]
        xwliq         ! land layer liquid water [kg/m2]

    real(r8), intent(in), dimension(ims:ime,jms:jme,1:12) :: &
        xirrgpadfrac, &! irrigated paddy area ratio
        xirrgnopfrac   ! irrigated non-paddy area ratio

    real(r8), intent(in), dimension(ims:ime,jms:jme) :: &
        xirrgthresh    ! the weighting factor which can control the irrigation target and amount

    integer, intent(in) :: im, id   ! month and days in month
    integer, intent(in), dimension(ims:ime,jms:jme) :: time_zone ! natural time zone
  ! ---------------------------------------------------------------
  ! dmd is inout, because if resstep != 0, then dmd will not be changed!
    real(r8), intent(inout), dimension(ims:ime,jms:jme) :: &
        xirrgwtrdmd_paddy,  & ! water demand for paddy irrigation [10^6 m3]
        xirrgwtrdmd_nonpaddy,  & ! water demand for non-paddy irrigaion [10^6m3]
        xirrgwtrdmd    ! water demand for irrgation [10^6 m3]

!   Local variables:
    integer :: i, j, k, l
    real(r8) :: water_demand_nonpaddy_mm, water_demand_paddy_mm   ! [mm]
    real(r8) :: water_demand_nonpaddy_mcm, water_demand_paddy_mcm ! [10^6 m3]
    real(r8) :: ddz, vol_liq, vol_ice
    real(r8) :: ef_porsl, ef_wdry
    real(r8) :: vol_fc, ef_wfld    ! vsm and effective vsm for field capacity
    real(r8) :: vol_wlt, ef_wwlt   ! vsm and effective vsm for wilting point
    real(r8) :: ef_thresh , f_thresh  ! irrigation moisture threshold
    real(r8) :: area
    logical :: isnotsoil
    integer :: localtime
!=======================================================================
        do j = jms, jme
          do i = ims, ime 
            localtime = modulo(resstep + time_zone(i,j), NRESSTEP)
            if (localtime == nint(NRESSTEP/3.)) then !local 8am
            ! test if not soil
            isnotsoil = nint(SC_WATER(i,j)) >= 4
#ifdef mask
            isnotsoil = isnotsoil .or. nint(calcumask(i,j)) == 0
#endif
            if (isnotsoil .or. (xirrgnopfrac(i,j,im) + xirrgpadfrac(i,j,im) == 0.)) then
              xirrgwtrdmd(i,j) = 0.0_r8
              xirrgwtrdmd_paddy(i,j) = 0.0_r8
              xirrgwtrdmd_nonpaddy(i,j) = 0.0_r8
              cycle
            end if
            ! start to calculate demand
            water_demand_nonpaddy_mm = 0.0_r8
            water_demand_paddy_mm = 0.0_r8
            ! loop over all soil layers
            do k = 1 , nl_soil
              l = k - nl_snow
              ! only calculate first 6 layers
              if (k > 6) exit
              ddz      = denliq * xdz(i,l,j)              ! [kg/m2 or mm h2o]
              vol_liq  = xwliq(i,l,j) / ddz
              vol_ice  = min(xporsl(i,k,j), xwice(i,l,j)/(denice*xdz(i,l,j)))
              ef_porsl = max(0., xporsl(i,k,j)- vol_ice)
              ef_wdry  = max(0., xwdry(i,k,j) - vol_ice)

              ! wilting point: suction = 15 bar = 1.5 * 10^5 mmH2O
              ! WARNING: here smpmax value is negative
              vol_wlt = xporsl(i,k,j) * (abs(smpmax) / xphi0(i,k,j))**(-1.0_r8 /xbsw(i,k,j))
              vol_wlt = min(xporsl(i,k,j), vol_wlt)
              vol_wlt = max(xwdry(i,k,j), vol_wlt)

              ! field capacity: suction = 100~ 400 mbar = 1k ~ 4k mmH2O
              ! here use 100 mbar = 1.02 mH2O
              vol_fc = xporsl(i,k,j) * (3400.0_r8 / xphi0(i,k,j))**(-1.0_r8 /xbsw(i,k,j))
              vol_fc = min(xporsl(i,k,j), vol_fc)
              vol_fc = max(xwdry(i,k,j), vol_fc)
              ef_wwlt = max(0., vol_wlt - vol_ice)
              ef_wfld = max(0., vol_fc - vol_ice)
              f_thresh = xirrgthresh(i,j)  ! a tuning parameter
              ef_thresh = (f_thresh*(ef_wfld-ef_wwlt)+ef_wwlt)

              ! 1. non-paddy crop demand
              if (vol_liq < ef_thresh) then
                water_demand_nonpaddy_mm = water_demand_nonpaddy_mm +  &
                                           (ef_thresh-vol_liq)*ddz   ! [kg/m2 or mm h2o]
              end if
              
              ! 2. paddy crop demand (here assume equal non-paddy crop)
              if (vol_liq < ef_thresh) then
                water_demand_paddy_mm =  water_demand_paddy_mm +  &
                                         (ef_thresh-vol_liq)*ddz   ! [kg/m2 or mm h2o]

              end if
            end do
            ! check nan
            if (water_demand_nonpaddy_mm /= water_demand_nonpaddy_mm) then
                water_demand_nonpaddy_mm = 0.
            endif
            
            if (water_demand_paddy_mm /= water_demand_paddy_mm) then
                water_demand_paddy_mm = 0.
            endif
#ifdef LonLat
            area = xarea(i,j)
#else
            area = xarea
#endif
            ! summing and converting unit, mm -> 10^6 m3
            water_demand_paddy_mcm =  &
                  xirrgpadfrac(i,j,im) * area * 1.0e-6_r8 *  &
                  water_demand_paddy_mm * 1.0e-3_r8

            water_demand_nonpaddy_mcm =  &
                  xirrgnopfrac(i,j,im) * area * 1.0e-6_r8 *  &
                  water_demand_nonpaddy_mm * 1.0e-3_r8

            ! output
            xirrgwtrdmd_nonpaddy(i,j) = water_demand_nonpaddy_mcm
            xirrgwtrdmd_paddy(i,j) = water_demand_paddy_mcm
            xirrgwtrdmd(i,j) = water_demand_nonpaddy_mcm + water_demand_paddy_mcm
  
            end if
          end do
        end do
  end subroutine cssp_irrg_waterdemand

  subroutine cssp_irrg_watersupply(dtime,xarea,time_zone,&
             xirrgwtrdmd,irrgswfrac,xswpumprate,xgwpumprate,&
             xirrgwtrdmd_paddy,xirrgwtrdmd_nonpaddy,xpumprate_paddy,&
             ims,ime,jms,jme,its,ite,jts,jte)
!=======================================================================
!   Purpose: Calculate water supply for irrigation
!   Authors: Chenyuan Li, Yang Jiao
!   [1]  Li, C., X. Yuan*, Y. Jiao, P.Ji, Z.Huang,2024: High-resolution Land Surface Modeling 
!        of the Irrigation Effects on Evapotranspiration over the Yellow River Basin. 
!        Journal of Hydrology, 633, 130986. 
!   [2]  X. Yuan*., C. Li, P. Ji, XZ. Xi, SW. Li, Y. Jiao, H. Yi, XZ. Liang, QC. Zeng, 2026:  
!        Development and evaluation of Conjunctive Surface-Subsurface Process Model version 3
!        (CSSPv3) at global scale. Journal of Advances in Modeling Earth Systems. Revised.
!=======================================================================
    implicit none
    integer, intent(in) :: &
        ims,ime,     &! start/end index for i in memory
        jms,jme,     &! start/end index for j in memory
        its,ite,     &! start/end index for i in tile
        jts,jte       ! start/end index for j in tile

    real, intent(in) :: &
        dtime         ! CLM time step [seconds]
#ifdef LonLat
    real, intent(in), dimension(ims:ime,jms:jme) :: &
        xarea          ! grid cell area [m2]
#else
    real, intent(in) :: &
        xarea          ! grid cell area [m2]
#endif
    integer(r8), intent(in), dimension(ims:ime,jms:jme) :: &
        time_zone
    real(r8), intent(in), dimension(ims:ime,jms:jme) :: &
        xirrgwtrdmd_nonpaddy,  & ! water demand for nonpaddy irrgation [10^6 m3]
        xirrgwtrdmd_paddy,  & ! water demand for paddy irrgation [10^6 m3]
        xirrgwtrdmd,   &! water demand for irrgation [10^6 m3]
        irrgswfrac    ! surface water irrigated area in the area equipped for irrigation ratio [-]

    real(r8), intent(inout), dimension(ims:ime,jms:jme) :: &
        xswpumprate,  &! surface water irrigate rate [mm h2o/s]
        xgwpumprate,  &! ground water irrigate rate [mm h2o/s]
        xpumprate_paddy ! total water irrigate rate for paddy [mm h2o/s] 
!
!   Local variables:
    integer :: i,j
    integer :: localtime
    real(r8) :: area
    real(r8) :: water_supply, water_supply_paddy
    real(r8), dimension(ims:ime,jms:jme) :: &
        xirrgindaywgt,        &! irrigation faction within a day [-]
        xirrgindaywgt_flood    ! flood irrigation faction within a day [-]
!=======================================================================
    xgwpumprate = 0.0_r8
    xswpumprate = 0.0_r8 
    xpumprate_paddy = 0.0_r8
    ! update xirrgrate for each timestep in a day

    do j = jms, jme
       do i = ims, ime
       ! Determine the irrigation timing
       localtime = modulo(resstep + time_zone(i,j), NRESSTEP)
       ! Normal irrigation happened from 8:00~12:00 (1/3day~1/2day)
       if (localtime >= nint(NRESSTEP/3.) .and. localtime < nint(NRESSTEP/2.)) then
          xirrgindaywgt(i,j) = 1 / real(nint(NRESSTEP / 2.) - nint(NRESSTEP / 3.))
       else
          xirrgindaywgt(i,j) = 0
       end if
       if (localtime == nint(NRESSTEP/3.)) then
          xirrgindaywgt_flood(i,j) = 1
       else
          xirrgindaywgt_flood(i,j) = 0
       end if

#ifdef LonLat
            area = xarea(i,j)
#else
            area = xarea
#endif
       water_supply = (xirrgwtrdmd_nonpaddy(i,j) * xirrgindaywgt(i,j) + xirrgwtrdmd_paddy(i,j) * xirrgindaywgt_flood(i,j)) &
                      / (area * 1.0e-6_r8) * 1.0e3_r8 / dtime
       water_supply_paddy = xirrgwtrdmd_paddy(i,j) * xirrgindaywgt_flood(i,j) &
                      / (area * 1.0e-6_r8) * 1.0e3_r8 / dtime

       ! Calucate Surface/Ground Water Withdrawal during a single timestep
       xswpumprate(i,j) = irrgswfrac(i,j) * water_supply
       xgwpumprate(i,j) = (1 - irrgswfrac(i,j)) * water_supply               
       ! Calucate Water Withdrawal for Paddy Flied during a single timestep
       xpumprate_paddy(i,j) = water_supply_paddy
      end do
    end do

  end subroutine cssp_irrg_watersupply

subroutine cssp_irrg_local_swpump(dtime,xarea,xswpumprate,xgwpumprate,xacswpumprate,xvolr,xdchrg,&
#ifdef Swpump
                            region_map,region_num,volr_region,swpump_region,xreswpump,&
#endif
#ifdef mask
                            calcumask,&
#endif
                            ims,ime,jms,jme,its,ite,jts,jte)
!=======================================================================
!   Purpose: surface water irrigation water pumping
!   Authors:  Chenyuan Li
!   Reference:
!   [1]  Li, C., X. Yuan*, Y. Jiao, P.Ji, Z.Huang,2024: High-resolution Land Surface Modeling 
!        of the Irrigation Effects on Evapotranspiration over the Yellow River Basin. 
!        Journal of Hydrology, 633, 130986. 
!   [2]  X. Yuan*., C. Li, P. Ji, XZ. Xi, SW. Li, Y. Jiao, H. Yi, XZ. Liang, QC. Zeng, 2026:  
!        Development and evaluation of Conjunctive Surface-Subsurface Process Model version 3
!        (CSSPv3) at global scale. Journal of Advances in Modeling Earth Systems. Revised.
! ===========================================================================
! Note:
!   CASE 1: Local pumping only (#undef SWPUMP)
!   - Both 'xswpumprate' and 'xgwpumprate' represent the source-specific 
!     irrigation amounts, which are equal to the local withdrawal rates.
!   - All water is extracted from and applied to the local grid cell.
!
!   CASE 2: Regional pumping considered (#define SWPUMP)
!   - 'xswpumprate': Represents ONLY the surface-water IRRIGATION amount applied 
!     to the local grid, not the actual withdrawal from it.
!   - 'xacswpumprate': Represents the ACTUAL surface-water WITHDRAWAL from the 
!     local grid cell, accounting for both local and regional water demands.
!   - 'xgwpumprate': Still represents both IRRIGATION and WITHDRAWAL, as 
!     groundwater is assumed to be extracted only from the local grid cell.

    implicit none
    integer, intent(in) :: &
        ims,ime,     &! start/end index for i in memory
        jms,jme,     &! start/end index for j in memory
        its,ite,     &! start/end index for i in tile
        jts,jte       ! start/end index for j in tile

    real, intent(in) :: &
        dtime          ! CLM time step [seconds]
#ifdef LonLat
    real, intent(in), dimension(ims:ime,jms:jme) :: &
        xarea          ! grid cell area [m2]
#else
    real, intent(in) :: &
        xarea          ! grid cell area [m2]
#endif
    real(r8), intent(inout), dimension(ims:ime,jms:jme) :: &
        xvolr,        &! cell h2o volume [m3]
        xdchrg,       &! total outflow [m3/s]
        xswpumprate,  &! surface water pumping/irrigate rate [mm h2o/s]
        xgwpumprate    ! ground water pumping/irrigate rate [mm h2o/s]
    real(r8), intent(out), dimension(ims:ime,jms:jme) :: &
        xacswpumprate  ! actual surface pumping rate [mm h2o/s]
#ifdef mask
    real, intent(in), dimension(ims:ime,jms:jme) :: &
        calcumask     ! mask
#endif
#ifdef Swpump
    integer,intent(in) :: region_num
    real, intent(in), dimension(ims:ime,jms:jme) :: &
        region_map
    real, intent(out), dimension(ims:ime,jms:jme) :: &
        xreswpump       ! The defcit of surface water pumping of the grid [m3]
    real, intent(out), dimension(region_num) :: &
        volr_region,   &! total cell h2o volume of the region [m3]
        swpump_region   ! The rest surface water pumping volume of the region [m3]
#endif
!   Local variables:
    integer :: i,j,k
    real(r8), dimension(ims:ime,jms:jme) :: &
        xreswpumprate  ! the deficit of surface water pumping [mm h2o/s]
    real :: area

! initialization 
xacswpumprate = 0.0_r8
xreswpumprate = 0.0_r8
#ifdef Swpump
    volr_region = 0.0_r8
    swpump_region = 0.0_r8
    xreswpump = 0.0_r8
#endif

   do j = jts,jte
      do i = its,ite
#ifdef LonLat
      area = xarea(i,j)
#else
      area = xarea
#endif
#ifdef mask
      if (nint(calcumask(i,j))==0) then
        cycle
      end if
#endif
      xacswpumprate(i,j) = min(xswpumprate(i,j),xvolr(i,j)/(1.e-3*area*dtime))
      xvolr(i,j) = max(xvolr(i,j)-xacswpumprate(i,j)*(1.e-3*area*dtime),0.)
      xreswpumprate(i,j) = xswpumprate(i,j) - xacswpumprate(i,j)
#ifndef Swpump
if (SW2GW) then
      xswpumprate(i,j) = xswpumprate(i,j) - xreswpumprate(i,j)
      xgwpumprate(i,j) = xgwpumprate(i,j) + xreswpumprate(i,j)
end if
#else
      k = region_map(i,j)
      if (k == 0) then !outside the basin
if (SW2GW) then
         xswpumprate(i,j) = xswpumprate(i,j) - xreswpumprate(i,j)
         xgwpumprate(i,j) = xgwpumprate(i,j) + xreswpumprate(i,j)   
end if
      else
! Using regions to satisfy the deficit of surface water over the irrigaiton areas
         xreswpump(i,j) = xreswpumprate(i,j) * (1.e-3*area*dtime)  !unit = m3
         if (xreswpump(i,j) /= xreswpump(i,j)) then !check nan 
            xreswpump(i,j) = 0.
         endif
         volr_region(k) = volr_region(k) + xvolr(i,j)
         swpump_region(k) = swpump_region(k) + xreswpump(i,j)
      end if
#endif
     xdchrg(i,j) = min(xvolr(i,j) / dtime, xdchrg(i,j)) ! revised the output, diagnostic variables
     end do
  end do

end subroutine cssp_irrg_local_swpump

#ifdef Swpump
subroutine cssp_irrg_region_swpump(dtime,xarea,xswpumprate,xgwpumprate,xacswpumprate,xvolr,xdchrg,&
                            region_map,region_num,volr_region,swpump_region,xreswpump,&
#ifdef mask
                            calcumask,&
#endif
                            ims,ime,jms,jme,its,ite,jts,jte)
!=======================================================================
!   Purpose: Dividing the sources of irrigation water regionally
!   Authors:  Chenyuan Li, Xing Yuan
!   Reference:
!   [1]  X. Yuan*., C. Li, P. Ji, XZ. Xi, SW. Li, Y. Jiao, H. Yi, XZ. Liang, QC. Zeng, 2026:  
!        Development and evaluation of Conjunctive Surface-Subsurface Process Model version 3
!        (CSSPv3) at global scale. Journal of Advances in Modeling Earth Systems. Revised.
! ===========================================================================
! Note:
!   CASE 1: Local pumping only (#undef SWPUMP)
!   - Both 'xswpumprate' and 'xgwpumprate' represent the source-specific 
!     irrigation amounts, which are equal to the local withdrawal rates.
!   - All water is extracted from and applied to the local grid cell.
!
!   CASE 2: Regional pumping considered (#define SWPUMP)
!   - 'xswpumprate': Represents ONLY the surface-water IRRIGATION amount applied 
!     to the local grid, not the actual withdrawal from it.
!   - 'xacswpumprate': Represents the ACTUAL surface-water WITHDRAWAL from the 
!     local grid cell, accounting for both local and regional water demands.
!   - 'xgwpumprate': Still represents both IRRIGATION and WITHDRAWAL, as 
!      groundwater is assumed to be extracted only from the local grid cell. 

    implicit none
    integer, intent(in) :: &
        ims,ime,     &! start/end index for i in memory
        jms,jme,     &! start/end index for j in memory
        its,ite,     &! start/end index for i in tile
        jts,jte       ! start/end index for j in tile

    real, intent(in) :: &
        dtime          ! CLM time step [seconds]
#ifdef LonLat
    real, intent(in), dimension(ims:ime,jms:jme) :: &
        xarea          ! grid cell area [m2]
#else
    real, intent(in) :: &
        xarea          ! grid cell area [m2]
#endif
    real(r8), intent(inout), dimension(ims:ime,jms:jme) :: &
        xvolr,        &! cell h2o volume [m3]
        xdchrg,       &! total outflow [m3/s]
        xswpumprate,  &! surface water irrigate rate [mm h2o/s]
        xgwpumprate    ! ground water irrigate rate [mm h2o/s]
    real(r8), intent(inout), dimension(ims:ime,jms:jme) :: &
        xacswpumprate  ! actual surface pumping rate [mm h2o/s]
    integer,intent(in) :: region_num
    real, intent(in), dimension(ims:ime,jms:jme) :: &
        region_map
    real, intent(inout), dimension(ims:ime,jms:jme) :: &
        xreswpump       ! The defcit of surface water pumping of the grid [m3]
    real, intent(in), dimension(region_num) :: &
        volr_region     ! total cell h2o volume of the region [m3]
    real, intent(in), dimension(region_num) :: &
        swpump_region   ! The rest surface water pumping volume of the region [m3]
#ifdef mask
    real, intent(in), dimension(ims:ime,jms:jme) :: &
        calcumask     ! mask
#endif
!   Local variables:
    integer :: i,j,k
    integer, dimension(region_num) :: &
        swpump_flag    ! wheather the volr_region can satisfy the swpumping demand
    real(r8), dimension(ims:ime,jms:jme) :: &
        xreswpumprate  ! the deficit of surface water pumping [mm h2o/s]
    real, dimension(region_num) :: &
        reswpump_region! the deficit of surface water in the region [m3]
    real :: area

! initialization
xreswpumprate = 0.

! whether the surface water withdrawal from the regional areas can meet the demand.
  do k = 1,region_num
     if ((volr_region(k) < swpump_region(k)))then
        swpump_flag(k) = 0 ! the demand can't be satisfied
        reswpump_region(k) = swpump_region(k)-volr_region(k)
     else
        swpump_flag(k) = 1 ! the demand can be satisfied
        reswpump_region(k) = 0
     end if
  end do

   do j = jts,jte
      do i = its,ite
#ifdef LonLat
      area = xarea(i,j)
#else
      area = xarea
#endif
#ifdef mask
      if (nint(calcumask(i,j))==0) then
        cycle
      end if
#endif
         k = region_map(i,j)
         if (k /= 0) then
           if (swpump_flag(k) == 0) then ! the demand can't be satisfied
             xacswpumprate(i,j) = xacswpumprate(i,j) + xvolr(i,j)/(1.e-3*area*dtime)
             xvolr(i,j) = 0.
             xreswpumprate(i,j) = reswpump_region(k)*(xreswpump(i,j)/swpump_region(k))/(1.e-3*area*dtime)
if (SW2GW) then
            xswpumprate(i,j) = xswpumprate(i,j) - xreswpumprate(i,j)
            xgwpumprate(i,j) = xgwpumprate(i,j) + xreswpumprate(i,j)
end if
           else
             if (volr_region(k) == 0.and.swpump_region(k) == 0) then ! to avoid both of region's volr and swpump = 0
                cycle
             end if
             xacswpumprate(i,j) = xacswpumprate(i,j) + (swpump_region(k)*(xvolr(i,j)/volr_region(k)))/(1.e-3*area*dtime)
             xvolr(i,j) = max((xvolr(i,j) - (swpump_region(k)*(xvolr(i,j)/volr_region(k)))), 0.)
           end if
         end if
      xreswpump(i,j) = xreswpumprate(i,j) * (1.e-3*area*dtime)  !unit = m3
      xdchrg(i,j) = min(xvolr(i,j) / dtime, xdchrg(i,j)) ! revised the output, diagnostic variables
      end do
   end do

end subroutine cssp_irrg_region_swpump
#endif

end module module_irrigation
#endif

! ===========================================================================
! Surface Water Withdrawal Coordination Logic
! NOTE: Conceptual reference only - Not executable
!   Authors: Chenyuan Li, Xing Yuan
!   Reference:
!   [1]  Li, C., X. Yuan*, Y. Jiao, P.Ji, Z.Huang,2024: High-resolution Land Surface Modeling 
!        of the Irrigation Effects on Evapotranspiration over the Yellow River Basin. 
!        Journal of Hydrology, 633, 130986. 
!   [2]  X. Yuan*., C. Li, P. Ji, XZ. Xi, SW. Li, Y. Jiao, H. Yi, XZ. Liang, QC. Zeng, 2026:  
!        Development and evaluation of Conjunctive Surface-Subsurface Process Model version 3
!        (CSSPv3) at global scale. Journal of Advances in Modeling Earth Systems. Revised.
! ===========================================================================
! Considering the issues of multi-core parallel interaction, a scheduling mechanism is implemented 
! when calculating the total regional volume. Please refer to the following code:
#ifdef IrrDmd
         call cssp_irrg_local_swpump(dtime,area,xswpumprate,xgwpumprate,xacswpumprate,xvolr,xdchrg,&
#ifdef Swpump 
                               region_map,region_num,volr_region,swpump_region,xreswpump,&
#endif
#ifdef mask
                               calcumask,&
#endif
                               ims,ime,jms,jme,its,ite,jts,jte)

#ifdef Swpump
allocate(fout(region_num))
call MPI_REDUCE(volr_region, fout, region_num, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
volr_region = fout
call MPI_BCAST(volr_region, region_num, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

call MPI_REDUCE(swpump_region, fout, region_num, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
swpump_region = fout
call MPI_BCAST(swpump_region, region_num, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
deallocate(fout)
        call cssp_irrg_region_swpump(dtime,area,xswpumprate,xgwpumprate,xacswpumprate,xvolr,xdchrg,&
                              region_map,region_num,volr_region,swpump_region,xreswpump,&
#ifdef mask
                              calcumask,&
#endif
                              ims,ime,jms,jme,its,ite,jts,jte)
#endif
#endif

! ===========================================================================
! CSSP-Irrigation Interaction Module Logic, which in the main physical module of CSSPv3
! NOTE: Conceptual reference only - Not executable
!   Authors: Chenyuan Li, Xing Yuan
!   Reference:
!   [1]  Li, C., X. Yuan*, Y. Jiao, P.Ji, Z.Huang,2024: High-resolution Land Surface Modeling 
!        of the Irrigation Effects on Evapotranspiration over the Yellow River Basin. 
!        Journal of Hydrology, 633, 130986. 
!   [2]  X. Yuan*., C. Li, P. Ji, XZ. Xi, SW. Li, Y. Jiao, H. Yi, XZ. Liang, QC. Zeng, 2026:  
!        Development and evaluation of Conjunctive Surface-Subsurface Process Model version 3
!        (CSSPv3) at global scale. Journal of Advances in Modeling Earth Systems. Revised.
! ===========================================================================

!============================================================================
! Irrigation Module: Ground Water Pumping
! Authors: Chenyuan Li, Xing Yuan
!============================================================================
!============================================================================
! Step 1: Water Table over Soil Layers
!============================================================================
   real(r8), intent(inout)
      tsw                  ! total soil water below the water table [mm]
   real(r8), intent(in)    :: 
      dtime              ,&! Time step [s]
      nl_soil            ,&! Number of soil layers
      dz(1:nl_soil)      ,&! Soil layer thickness [m]
      porsl(1:nl_soil)   ,&! Saturated volumetric water content [-]
      vol_liq(1:nl_soil) ,&! Current liquid water content [-]
      ef_wdry(1:nl_soil) ,&! Effective dry soil moisture limit [-]
      gwpumprate           ! ground water pumping rate [mm h2o/s]
   real(r8) :: &
      sz                 ,&! thickness of the soil layer below the water table
      jw                 ,&! index of the soil layer where the water table is located
      ws                 ,&! water used to fill soil air pores regardless of water content        
      total_pump         ,&! total amount of water remaining to be pumped from the soil layers
      ef_pump            ,&! the total water content can be pumped of k soil-layer
      qpump(1:nl_soil)     ! groundwater pumping from each soil layer(if water table < 5.67m) 

      sz = 0.         ! thickness of the soil layer below the water table
      jw = 0          ! index of the soil layer where the water table is located
      ws = 0.         ! water used to fill soil air pores regardless of water content        
      total_pump = 0. ! total amount of water remaining to be pumped from the soil layers
      ef_pump = 0.    ! the total water content can be pumped of k soil-layer
      qpump(:) = 0    ! groundwater pumping from each soil layer(if water table < 5.67m) 

     if (gwpumprate > 0) then
       if (tsw > 0.) then
          do k = nl_soil,1,-1
             ws = ws + porsl(k) * 1000 * dz(k) ![mm]
             if (tsw <= ws) then
               jw = k
               exit
             endif
          enddo
          sz = sum(dz(jw:nl_soil))                       ![m]
          total_pump = max(0.,min(tsw / dtime, gwpumprate)) ![mm/s]
        
          do k = jw,nl_soil,1
             ef_pump = (vol_liq(k) - ef_wdry(k)) * dz(k) * 1000 / dtime ![mm/s]
             qpump(k) = max(0., min(dz(k) / sz * total_pump, ef_pump))  ![mm/s]
             ! the rest layers continue to satisfy the water pumping demand
             total_pump = total_pump - qpump(k)
             sz = sz - dz(k)
           enddo
       endif
     end if
!============================================================================
! Step 2: qpump present groundwater pumping from each soil layer
!         which serves as a part of forcing term of tridiagonal matrix 
!         of the Richards equation.
!============================================================================
    rmx(1:nl_soil)      ,&! (r) forcing term of tridiagonal matrix
    
    do k = 2, nl_soil_br    ! inner nodes and bottom node
       rmx(k) = rmx(k) &
#ifdef IrrDmd
            - qpump(k)
#endif
    end 
!============================================================================
! Step 3: Renew Ground Water Storage
!============================================================================
tsw = tsw + qrchrg * dtime & ! qrchrg, groundwater recharge rate [mm h2o/s]
#ifdef IrrDmd
      - gwpumprate * dtime
#endif

!============================================================================
! Irrigation Module: Irrigation on Non-Paddy and Paddy Flieds
! NOTE: Conceptual reference only - Not executable
! Authors: Chenyuan Li, Xing Yuan
! Reference:
!   [1]  Li, C., X. Yuan*, Y. Jiao, P.Ji, Z.Huang,2024: High-resolution Land Surface Modeling 
!        of the Irrigation Effects on Evapotranspiration over the Yellow River Basin. 
!        Journal of Hydrology, 633, 130986. 
!   [2]  X. Yuan*., C. Li, P. Ji, XZ. Xi, SW. Li, Y. Jiao, H. Yi, XZ. Liang, QC. Zeng, 2026:  
!        Development and evaluation of Conjunctive Surface-Subsurface Process Model version 3
!        (CSSPv3) at global scale. Journal of Advances in Modeling Earth Systems. Revised.
!============================================================================
!------------------------- dummy argument -----------------------------------
  real(r8), intent(in) :: &
      dtime               ,&! time step [s]
      gwat                ,&! net water input from top [mm h2o/s]
      rainwater           ,&! rainwater rate [mm/s]
      swdepth               ! surface water depth [mm]

#ifdef IrrDmd
  real(r8), intent(in) :: &
      qpump(1:nl_soil)   ,&! groundwater pumping from each soil layer(if water table < 5.67m)
      swpumprate         ,&! surface water pumping rate [mm h2o/s]
      gwpumprate         ,&! ground water pumping rate [mm h2o/s]
      pumprate_paddy     ,&! total water pumping/irrigate rate for paddy [mmh2o/s]
      area               ,&! grid cell area [m2]
      qinmax             ,&! maximum possible infux on the top [mm/s]
      irrgpadfrac          ! irrigated paddy area ratio [-]
  real(r8), intent(inout) :: swdepth_paddy ! surface water depth in paddy field [mm]
#endif

!-------------------------- local variables ---------------------------------
  real(r8)
      rsur                ,&! hortonian (rainfall excess) surface runoff [mm/s]
      soil_inout          ,&! water inout the surface soil
      
      
#ifdef IrrDmd
  real(r8) :: & 
      swdepth_thresh      ,&! the ponding threshold of surface water
      ef_qinfl            ,&! the rest of possible infiltration rate for paddy field [mm/s]
      qinfl_paddy         ,&! the maximum possible infiltration rate for paddy field[mm/s]
      pumprate_nonpaddy     ! total water pumping/irrigate rate for non-paddy [mmh2o/s]
#endif

! Irrigation Processing
!============================================================================
! Step 1: Pre-partitioning of Water Input
!============================================================================
#ifdef IrrDmd
! Irrigation is added as precipitation
      pumprate_nonpaddy = swpumprate + gwpumprate - pumprate_paddy
      swdepth_paddy = swdepth_paddy + (gwat * irrgpadfrac + pumprate_paddy) * dtime
      rainwater = ((1 - irrgpadfrac) * gwat + pumprate_nonpaddy) + swdepth/dtime
#else
      rainwater = gwat + swdepth/dtime  ! [mm/s]
#endif

!============================================================================
! Step 2: Natural Hydrology Process (e.g., VIC Scheme)
! NOTE: This section simplifies the core VIC runoff and infiltration logic.
! Based on the Variable Infiltration Capacity (VIC) curve, the model estimates:
! rsur:  surface runoff [mm/s]
! The surface runoff generation and infiltration processes are driven by rainwater
! calculated in Step 1. In Step 1, rainwater has been partitioned into
! paddy field and other components.
! Only the rainwater over other areas contributes to runoff generation
! following the VIC scheme. The paddy field portion is assumed to be retained
! within the field (e.g., for ponding) and does not directly generate surface runoff.

!============================================================================
! Step 3: Paddy Modification on Surface Runoff and Infiltration
! When the water depth is below the threshold, the paddy field only undergoes 
! infiltration and evaporation; notably, the evaporation loss has already been
! incorporated into the net water input (rainwater/gwat).
!============================================================================
#ifdef IrrDmd
      ! calculate the infiltration of paddy field
      ! the water pond in paddy field < 100mm, keeping the surface water runoff
      ! to zero in paddy field
      swdepth_thresh = irrgpadfrac * 100.
      if (swdepth_paddy <= 0.) then
      ! swdepth_paddy have been evaporated
         qinfl_paddy = swdepth_paddy/dtime
         swdepth_paddy = 0.
      else
         ef_qinfl = irrgpadfrac * qinmax
         qinfl_paddy = min((swdepth_paddy/dtime),ef_qinfl)  
         swdepth_paddy = max(0., (swdepth_paddy - qinfl_paddy * dtime))
         ! renew surface runoff
         if (swdepth_paddy > swdepth_thresh) then
         rsur = rsur + max(0., (swdepth_paddy - swdepth_thresh)/dtime)
         swdepth_paddy = swdepth_thresh
         end if
      endif
#endif
!============================================================================
! Step 4: qinfl_paddy represents the infiltration from the paddy field,
!         which serves as a part of forcing term of tridiagonal matrix 
!         of the Richards equation.
!============================================================================
    rmx(1:nl_soil)        ! (r) forcing term of tridiagonal matrix
    k = 1                 ! surface node
    rmx(k) = rmx(k) &
#ifdef IrrDmd
            + qinfl_paddy
#endif