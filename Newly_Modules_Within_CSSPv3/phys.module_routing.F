#include "configs.inc"

! start of routing module

module module_routing
!=======================================================================
!
!   Purpose: Grid-based horizontal routing using linear reservoirs method
!   Linear reservoir routing: Xing Yuan, 2009-06-01
!   Revised for flow velocity and LonLat coordinate, Peng Ji, 2022-05-01
!   Revised for modularization and reservoir, Yang Jiao, 2023-04-30
!   Further updated: modularization and reservoir, Chenyuan Li, 2024-03-06
!   Reference:
!   [1]  X. Yuan*., C. Li, P. Ji, XZ. Xi, SW. Li, Y. Jiao, H. Yi, XZ. Liang, QC. Zeng, 2026:  
!        Development and evaluation of Conjunctive Surface-Subsurface Process Model version 3
!        (CSSPv3) at global scale. Journal of Advances in Modeling Earth Systems. Revised.
!=======================================================================
  use shr_kind_mod, ONLY: r8 => shr_kind_rn
  implicit none
  private
  public :: linear_transport           ! method: grid-based linear reservoirs
contains
  subroutine linear_transport(SC_WATER,dtime,deltx,delty,area,dirgds,xdirid, &
      xrnof,xvolr,xdchrg, &
#ifdef para_vic
      xwtslope,ROUT_VMIN,ROUT_PSLP,ROUT_VOUT, &
#endif
#ifdef ResOpn
      damcap,damstg,damyear,iy,im, &
#ifdef PRes
      outrate,ratio_min,damsupply,avgflw,&
#else
      outrate,Sflo,avgflw,VarH,&
#endif
      xinflow,xdelt_stg,&
#endif
#ifdef mask
      calcumask,                 &
#endif
#ifdef UpStrm
      upsnum, dirstrm, upxids, upyids, &
#endif
      ids,ide ,jds,jde ,kds,kde ,&
      ims,ime ,jms,jme ,kms,kme ,&
      its,ite ,jts,jte ,kts,kte )
!=======================================================================
!
!   Routing with effective water flow velocity from Dingman, 2002
!
!=======================================================================
!   Uses:
    USE module_sf_clm_constant, ONLY: idlist,jdlist,spv_water
#ifdef ResOpn
#ifdef PRes
    USE module_reservoir, ONLY: resvr_release_default
#else
    USE module_reservoir, ONLY: resvr_release_management
#endif
#endif
!
    implicit none
!
!   Input/Output variables:
    integer, intent(in) :: &
        ids,ide,     &! start/end index for i in domain
        jds,jde,     &! start/end index for j in domain
        kds,kde,     &! start/end index for k in domain
        ims,ime,     &! start/end index for i in memory
        jms,jme,     &! start/end index for j in memory
        kms,kme,     &! start/end index for k in memory
        its,ite,     &! start/end index for i in tile
        jts,jte,     &! start/end index for j in tile
        kts,kte       ! start/end index for k in tile

    real, intent(in), dimension(ims:ime,jms:jme) :: &                   ! halo 8
        SC_WATER      ! surface category of water characteristics
                      ! 1=urban and built-up, 2=soil, 3=wetland, 4=land ice
                      ! 5=shallow lake, 6=deep lake,  7=sea ice, 8=ocean

#ifdef mask
    real, intent(in), dimension(ims:ime,jms:jme) :: &
        calcumask     ! mask
#endif

    real, intent(in) :: &
        dtime         ! model time step [seconds]
#ifdef LonLat
    real, intent(in), dimension(ims:ime,jms:jme) :: &
        deltx,       &! 2D x-dir grid spacing  [m]
        delty,       &! 2D y-dir grid spacing  [m]
        area          ! grid cell area [m2]
#else
    real, intent(in) :: &
        deltx,       &! x-dir grid spacing [m]
        delty,       &! y-dir grid spacing [m]
        area          ! grid cell area [m2]
#endif

#ifdef LonLat
    real, intent(in), dimension(ims:ime,8,jms:jme) :: &
        dirgds        ! 3D directional grid spacing [m]
#else
    real, intent(in), dimension(8) :: &
        dirgds        ! directional grid spacing [m]
#endif
!----------------------- Need Halo 8/24 points ------------------------------
    integer, intent(in), dimension(ims:ime,8,jms:jme) :: &
        xdirid        ! 8-directional flow direction index              ! halo 8
                      !    +---+---+---+
                      !    | 6 | 7 | 8 |
                      !    +---+---+---+
                      !    | 5 | o | 1 |
                      !    +---+---+---+
                      !    | 4 | 3 | 2 |
                      !    +---+---+---+
                      ! note xdirid(i,l,j)>0 means inflow from
                      !    the neighbor grid(l) to the center grid(i,j)
                      !    & xdirid(i,l,j)<0 means outflow from
                      !    the center grid(i,j) to the neighbor grid(l)

    real, intent(in), dimension(ims:ime,jms:jme) :: &
        xrnof         ! total runoff [mm h2o/s]

#ifdef para_vic
    real, intent(in), dimension(ims:ime,jms:jme) :: &
        xwtslope,    &! water table slope, the same as surface slope
        ROUT_VMIN,   &! mininum effective water flow velocity [m/s]
        ROUT_PSLP     ! factor for slope
    real, intent(in), dimension(ims:ime,jms:jme) :: ROUT_VOUT !+cyli
#endif
    real, intent(inout), dimension(ims:ime,jms:jme) :: &
        xvolr,       &! cell h2o volume [m3]
        xdchrg        ! total outflow [m3/s]
#ifdef ResOpn
    real(r8), intent(in), dimension(ims:ime,jms:jme) :: &
        damcap        ! reservoir capacity [10^6 m3]
    real(r8), intent(inout), dimension(ims:ime,jms:jme) :: &
        damstg,      &! reservoir storage [10^6 m3]
        xdelt_stg      ! changes in reservoir stotage [m3]
    real(r8), intent(in), dimension(ims:ime,jms:jme) :: &
        damyear       ! the year when reservoir built
    integer, intent(in) :: iy,im ! model integration year and month
#ifdef PRes
    real(r8), intent(in), dimension(ims:ime,jms:jme) :: &
        avgflw,      &! long-term averaged inflow [m3/s]
        outrate,     &! the ratio of 'outlfow' based on 'avginflow' and 'inflow' (default = 0.5)
        ratio_min     ! the mininum release coefficient(default = 0.5)
    real(r8), intent(in), dimension(ims:ime,jms:jme,1:12) :: &
        damsupply   ! the water supplied by the reservoir to the downstream [m3/s]
#else
    real(r8), intent(in), dimension(ims:ime,jms:jme) :: &
        outrate,     &! the ratio of 'outlfow' based on 'avginflow' and 'inflow' (default = 0.5)
        Sflo          !
    real(r8), intent(in), dimension(ims:ime,jms:jme,1:12) :: & 
        avgflw,      &! long-term averaged inflow [m3/s]
        VarH
#endif
    real(r8),intent(out), dimension(ims:ime,jms:jme) :: xinflow
    real(r8) :: avginflow
#endif

#ifdef UpStrm
    integer, intent(in) :: &
        upsnum,      &! number of upstream watersheds (MAX 2)
        upxids(2),   &! upstream inflow at dim x
        upyids(2)     ! upstream inflow at dim y
    character(len=255) :: &
        dirstrm       ! inflow dir
#endif

!   Local variables:
    real(r8), dimension(ims:ime,jms:jme) :: &
        pvolr,       &! previous cell h2o volume [m3]
        pdchrg        ! previous total outflow [m3/s]

    real(r8), dimension(ims:ime,jms:jme) :: &
        vout          ! effective water flow velocity [m/s]

    real(r8) :: inflow, sumq
    real(r8) :: vout_default, vmin_default, pslp_default

#ifdef ResOpn
    real(r8) :: Qresout
#endif

#ifdef UpStrm
    integer :: iups
#endif

    integer :: i, j, k, l, id, ifr, jfr
    integer :: il, ir, jl, ju             ! local edges with len=BUFLEN-1

!=======================================================================
    vout_default = 0.8_r8
    vmin_default = 0.05_r8
    pslp_default = 1.0_r8

    ! memory-local-tile, determine local edges
    il  = max(ids  , its-2)
    ir  = min(ide-1, ite+2)
    jl  = max(jds  , jts-2)
    ju  = min(jde-1, jte+2)

    ! set default for initialization (needed at the domain edges)
    pvolr  = xvolr
    pdchrg = xdchrg
    ! set parameter
    do i = ims, ime
      do j = jms, jme
#ifdef para_vic
        vout(i,j) = ROUT_VOUT(i,j) !+cyli
#else
        vout(i,j) = vout_default
#endif
      end do
    end do

!
!---do the predictor step
!
    ! inflow: calculate cell h2o volume
    do i = its, ite
      do j = jts, jte
#ifdef mask
         if (nint(calcumask(i,j))==0) then 
           xvolr(i,j) = spv_water
           cycle
         end if
#endif
         inflow = 0.0_r8
#ifdef UpStrm
         if (upsnum > 0) then
           do iups = 1, upsnum
             if (i == upxids(iups) .and. j == upyids(iups)) then
             end if
           end do
         end if
#endif
        do l = 1,8
          id = xdirid(i,l,j)
          if (id <= 0) exit
          !-inflow part
          ifr = i + idlist(id)
          jfr = j + jdlist(id)
          !--no inflow from outside of the model domain
          if (jfr < jl .or. jfr > ju .or. ifr < il .or. ifr > ir) cycle
          if (nint(calcumask(ifr,jfr))==1) inflow = inflow + pdchrg(ifr,jfr)
        end do

#ifdef LonLat
          pvolr(i,j) = pvolr(i,j) + dtime * (inflow - pdchrg(i,j) + xrnof(i,j) * 1.0e-3_r8 * area(i,j))
#else
          pvolr(i,j) = pvolr(i,j) + dtime * (inflow - pdchrg(i,j) + xrnof(i,j) * 1.0e-3_r8 * area)
#endif
        pvolr(i,j) = max(0.0_r8, pvolr(i,j))

        ! update xvolr
        xvolr(i,j)  = pvolr(i,j)
      end do
    end do

    ! outflow: calculate total outflow
    do i = its, ite
      do j = jts, jte
#ifdef mask
          if (nint(calcumask(i,j))==0) then 
            xdchrg(i,j) = spv_water
            cycle
          end if
#endif
        sumq = 0.0_r8
        do l = 1,8
          id = xdirid(i,l,j)
          if (id < 0) then    !-outflow part
#ifdef LonLat
            sumq = sumq + vout(i,j) * pvolr(i,j) / dirgds(i,-id,j)
#else
            sumq = sumq + vout(i,j) * pvolr(i,j) / dirgds(-id)
#endif
            exit  ! only one outlet
          end if
        end do
        
        pdchrg(i,j) = min(pvolr(i,j) / dtime, sumq)
        xdchrg(i,j) = pdchrg(i,j)
#ifdef ResOpn
        xinflow(i,j)= 0
#ifdef PRes
        avginflow = avgflw(i,j)
#else
        avginflow = avgflw(i,j,im)
#endif
        if (damcap(i,j) > 0.0_r8 .and. REAL(iy) >= damyear(i,j) ) then
              xinflow(i,j)= pdchrg(i,j) 
          if (avginflow <= 0.) then
              avginflow = pdchrg(i,j) ! use current pdchrg as the avgflw
          endif
          call resvr_release(dtime,damcap(i,j),damstg(i,j),avginflow,xinflow(i,j),Qresout,&
                             outrate(i,j),ratio_min(i,j),damsupply(i,j,im))

          !(a) Qresout < pdchrg(i,j) means water moves from grid to reservoir
          !(b) Qresout > qdchrg(i,j) means water moves from reservoir to grid
          xdchrg(i,j) = Qresout
          xdelt_stg(i,j) = (pdchrg(i,j) - xdchrg(i,j)) * dtime
          xvolr(i,j) = xvolr(i,j) - xdelt_stg(i,j)
        end if
#endif
      end do
    end do

  end subroutine linear_transport

end module module_routing
