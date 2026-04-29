#include "configs.inc"

! start of reservoir module

#ifdef ResOpn

module module_reservoir
!=======================================================================
!   Purpose: Simulate reservoir influence on hydrological cycle with a 
!   prospective reservoir operation scheme.
!   Authors: Yang Jiao, Li Chenyuan, Xing Yuan, 2023-08-29
!   Revised for considering outrate, mininum outflow, Li Chenyuan, 2024-02-17
!   Revised for considering watersupply, Li Chenyuan, 2024-03-06
!   Reference:
!   [1]  X. Yuan*., C. Li, P. Ji, XZ. Xi, SW. Li, Y. Jiao, H. Yi, XZ. Liang, QC. Zeng, 2026:  
!        Development and evaluation of Conjunctive Surface-Subsurface Process Model version 3
!        (CSSPv3) at global scale. Journal of Advances in Modeling Earth Systems. Revised.
!=======================================================================
  use shr_kind_mod, ONLY: r8 => shr_kind_rn
  implicit none
  private
  public :: resvr_release
contains

  subroutine resvr_release(dtime,cap,stg,avginflow,inflow,outflow,&
                           outrate,ratio_min,watersupply)
!=======================================================================
!   Purpose: calculate reservoir outflow based on a simple model 
!=======================================================================    
    implicit none

!   Input/Output variables:
    real, intent(in) :: &
        dtime         ! model time step [seconds]
    real(r8), intent(in) :: &
        cap,         &! reservoir capacity [10^6 m3]
        avginflow,   &! averaged reservoir inflow [m3/s]
        inflow        ! reservoir inflow [m3/s]
    real(r8), intent(inout) :: &
        stg           ! reservoir storage [10^6 m3]
    real(r8), intent(out) :: &
        outflow       ! reservoir outflow [m3/s]

!   Extra input variables:
    real(r8), intent(in) :: &
        outrate,      &! the ratio of 'outlfow' based on 'avginflow' and 'inflow' (default = 0.5)
        ratio_min,    &! the mininum release coefficient(default = 0.15)
        watersupply    ! the water supplied by the reservoir to the downstream [m3/s]
  
!   Local variables:
    real(r8) :: Smax        ! total capacity (100%)
    real(r8) :: Smin        ! dead capacity (10%)
    real(r8) :: Shigh       ! flood (~flood control high level) (35+40%)
    real(r8) :: factor_Q    ! release factor if flooding
    real(r8) :: outmin      ! the mininum outflow (generate electricity and water supply)
    real(r8) :: pstg, ratio 
!=======================================================================

    Smax  = 1.00_r8 * cap
    Smin  = 0.20_r8 * cap
    Shigh = 0.80_r8 * cap
    factor_Q = 1.5_r8
    pstg = stg

    ! 1. default release: based on long-term (at least 5 years) averaged streamflow
    ratio = (stg - Smin) / (Shigh - Smin)
    ratio = max(0.0_r8, min(1.0_r8, ratio))

    ! outflow = ratio * avginflow, No used PCR-Scheme.
    outmin = ratio_min * avginflow
    outflow = outrate * (ratio * avginflow) + (1 - outrate) * inflow + watersupply
    outflow = max(outflow,outmin)
    stg = pstg + (inflow - outflow) * dtime / 1.0e6_r8

    ! 2. check storage
    if (stg > Shigh) then
      ! larger than capacity, release extra flow
      if (stg > Smax) then
         outflow = outflow + (stg - Smax) * 1.0e6_r8 / dtime
      end if
      
      ! larger than flood control high level, use flooding factor
      ratio = (stg - Shigh) / (Smax - Shigh)
      ratio = max(0.0_r8, min(1.0_r8, ratio))
      outflow = outflow + factor_Q * ratio * avginflow

      ! outflow should be larger than average inflow
      outflow = max(avginflow, outflow)
      stg = pstg + (inflow - outflow) * dtime / 1.0e6_r8
    end if

    ! lower than dead capacity, reduce outflow
    if (stg < Smin) then
       if ( (Smin - stg) * 1.0e6_r8 / dtime > outflow) then
          outflow = 0.
          stg = pstg + inflow * dtime / 1.0e6_r8
       else
          outflow = outflow - (Smin - stg) * 1.0e6_r8 / dtime
          stg = Smin
       endif
    end if

  end subroutine resvr_release

end module module_reservoir

#endif
