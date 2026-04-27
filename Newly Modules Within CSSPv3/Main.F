! ===========================================================================
! Main Controller of CSSPv3
! NOTE: Conceptual reference only - Not executable
! ===========================================================================
#ifdef IrrDmd
    time_zone = NINT(XLON/15 * NRESSTEP/24)  !for N time step
#endif
       call CSSP(...) 
       ! Main physical component of CSSP
       ! Includes soil water and groundwater equations, soil water–groundwater exchange,
       ! VIC-runoff, snow processes, vegetation-related processes and energy processes...
       ! as well as certain components of newly developed modules in CSSPv3 (interactively coupled parts):
       ! [1] Irrigation-related processes (e.g., groundwater withdrawal irrigation 
       !     over the irrigated field, paddy field ponding), directly coupled with soil water and 
       !     groundwater equations, VIC-runoff.
       ! [2] Carbon Cycle Processes, directly coupled with vegetation-related processes
       ! [3] Glacier Processes, includes glacier surface energy balance, mass balance, 
       !     and meltwater routing, coupled with energy processes, snow processes, and hydrology modules
       
#ifdef IrrDmd
       call cssp_irrg(xnl_soil_br,SC_WATER,nl_soil,nl_snow,dtime,area, &
            xdz,xhksm,xbsw,xphi0,xporsl,xwdry,xwliq,xwice, &
            xirrgpadfrac,xirrgnopfrac,xirrgswfrac,xirrgthresh, &
            xirrgwtrdmd,xswpumprate,xgwpumprate, &
            xirrgwtrdmd_paddy,xirrgwtrdmd_nonpaddy,xpumprate_paddy, &
            im,id,time_zone,&  ! Irrigation amount and water allocation
#ifdef mask
            calcumask,                 &
#endif
            ids,ide ,jds,jde ,kds,kde ,&
            ims,ime ,jms,jme ,kms,kme ,&
            its,ite ,jts,jte ,kts,kte )
#endif

       call linear_transport(xrnof, xvolr,xdchrg,...) ! Routing and Reservoirs
       ! xrnof,        &! runoff
       ! xvolr,        &! cell h2o volume [m3]
       ! xdchrg         ! total outflow [m3/s]

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