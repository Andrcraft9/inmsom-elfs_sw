!allocation of arrays
subroutine model_grid_allocate
    use mpi_parallel_tools
    use basin_grid
    implicit none

    call allocate_block2D_real8(lu, 0.0d0)  !mask of t-grid
    call allocate_block2D_real8(lu1, 0.0d0) !mask of t-grid (1 everywhere)
    call allocate_block2D_real8(luu, 0.0d0) !mask of h-grid (0 on boundary)
    call allocate_block2D_real8(luh, 0.0d0) !mask of h-grid (1 on boundary)
    call allocate_block2D_real8(lcu, 0.0d0) !mask of u-grid (0 on boundary)
    call allocate_block2D_real8(lcv, 0.0d0) !mask of v-grid (0 on boundary)
    call allocate_block2D_real8(llu, 0.0d0) !mask of u-grid (0 on boundary)
    call allocate_block2D_real8(llv, 0.0d0) !mask of v-grid (0 on boundary)

    call allocate_block2D_real8(hhh, 0.0d0)   !ocean depth on luh (h-points)
    call allocate_block2D_real8(hhhp, 0.0d0)  !ocean depth on luh (h-points) at previous step
    call allocate_block2D_real8(hhhn, 0.0d0)  !ocean depth on luh (h-points) at pre-previous step
    call allocate_block2D_real8(hhq_rest, 0.0d0)  !ocean depth on lu  (t-points) at rest state
    call allocate_block2D_real8(hhq, 0.0d0)       !ocean depth on lu  (t-points)
    call allocate_block2D_real8(hhqp, 0.0d0)  !ocean depth on lu  (t-points) at previous step
    call allocate_block2D_real8(hhqn, 0.0d0)  !ocean depth on lu  (t-points) at pre-previous step
    call allocate_block2D_real8(hhu, 0.0d0)   !ocean depth on lcu (u-points)
    call allocate_block2D_real8(hhup, 0.0d0)  !ocean depth on lcu (u-points) at previous step
    call allocate_block2D_real8(hhun, 0.0d0)  !ocean depth on lcu (u-points) at pre-previous step
    call allocate_block2D_real8(hhv, 0.0d0)   !ocean depth on lcv (v-points)
    call allocate_block2D_real8(hhvp, 0.0d0)  !ocean depth on lcv (v-points) at previous step
    call allocate_block2D_real8(hhvn, 0.0d0)  !ocean depth on lcv (v-points) at pre-previous step
    call allocate_block2D_real8(rlh_s, 0.0d0)  !coriolis-1 parameter on edge (t-centers) points
    call allocate_block2D_real8(rlh_c, 0.0d0)  !coriolis-2 parameter on edge (t-centers) points

    call allocate_block2D_real8(dxt, 0.0d0)
    call allocate_block2D_real8(dyt, 0.0d0)  !horizontal grid steps between   t-points (in radians or meters)
    call allocate_block2D_real8(dx, 0.0d0)
    call allocate_block2D_real8(dy, 0.0d0)   !horizontal grid steps between u,v-points (in radians or meters)
    call allocate_block2D_real8(dxh, 0.0d0)
    call allocate_block2D_real8(dyh, 0.0d0)  !horizontal grid steps between   h-points (in radians or meters)
    call allocate_block2D_real8(dxb, 0.0d0)
    call allocate_block2D_real8(dyb, 0.0d0)  !horizontal grid steps between v,u-points (in radians or meters)
    call allocate_block2D_real8(xt, 0.0d0)
    call allocate_block2D_real8(yt, 0.0d0)   !horizontal t-grid            x- and y-coordinates (in degrees)
    call allocate_block2D_real8(xu, 0.0d0)
    call allocate_block2D_real8(yv, 0.0d0)   !horizontal u-grid and v-grid x- and y-coordinates (in degrees)

    call allocate_block2D_real8(geo_lon_t, 0.0d0)    !geographical longitudes of T-points
    call allocate_block2D_real8(geo_lat_t, 0.0d0)    !geographical latitudes  of T-points
    call allocate_block2D_real8(geo_lon_u, 0.0d0)    !geographical longitudes of U-points
    call allocate_block2D_real8(geo_lat_u, 0.0d0)    !geographical latitudes  of U-points
    call allocate_block2D_real8(geo_lon_v, 0.0d0)    !geographical longitudes of V-points
    call allocate_block2D_real8(geo_lat_v, 0.0d0)    !geographical latitudes  of V-points
    call allocate_block2D_real8(geo_lon_h, 0.0d0)    !geographical longitudes of H-points
    call allocate_block2D_real8(geo_lat_h, 0.0d0)    !geographical latitudes  of H-points

endsubroutine model_grid_allocate
!-------------------------------------------------------------------------------

!deallocation of arrays
subroutine model_grid_deallocate
    use basin_grid
    use mpi_parallel_tools
    implicit none

    call deallocate_block2D_real8(lu)
    call deallocate_block2D_real8(lu1)
    call deallocate_block2D_real8(luu)
    call deallocate_block2D_real8(luh)
    call deallocate_block2D_real8(lcu)
    call deallocate_block2D_real8(lcv)
    call deallocate_block2D_real8(llu)
    call deallocate_block2D_real8(llv)

    call deallocate_block2D_real8(hhh)
    call deallocate_block2D_real8(hhhp)
    call deallocate_block2D_real8(hhhn)
    call deallocate_block2D_real8(hhq_rest)
    call deallocate_block2D_real8(hhq)
    call deallocate_block2D_real8(hhqp)
    call deallocate_block2D_real8(hhqn)
    call deallocate_block2D_real8(hhu)
    call deallocate_block2D_real8(hhup)
    call deallocate_block2D_real8(hhun)
    call deallocate_block2D_real8(hhv)
    call deallocate_block2D_real8(hhvp)
    call deallocate_block2D_real8(hhvn)
    call deallocate_block2D_real8(rlh_s)
    call deallocate_block2D_real8(rlh_c)

    call deallocate_block2D_real8(dxt)
    call deallocate_block2D_real8(dyt)
    call deallocate_block2D_real8(dx)
    call deallocate_block2D_real8(dy)
    call deallocate_block2D_real8(dxh)
    call deallocate_block2D_real8(dyh)
    call deallocate_block2D_real8(dxb)
    call deallocate_block2D_real8(dyb)
    call deallocate_block2D_real8(xt)
    call deallocate_block2D_real8(yt)
    call deallocate_block2D_real8(xu)
    call deallocate_block2D_real8(yv)

    call deallocate_block2D_real8(geo_lon_t)
    call deallocate_block2D_real8(geo_lat_t)
    call deallocate_block2D_real8(geo_lon_u)
    call deallocate_block2D_real8(geo_lat_u)
    call deallocate_block2D_real8(geo_lon_v)
    call deallocate_block2D_real8(geo_lat_v)
    call deallocate_block2D_real8(geo_lon_h)
    call deallocate_block2D_real8(geo_lat_h)

endsubroutine model_grid_deallocate
!-------------------------------------------------------------------------------

!allocation of arrays
subroutine ocean_variables_allocate
    use mpi_parallel_tools
    use ocean_variables
    implicit none

    call allocate_block2D_real8(ssh, 0.0d0)     !sea surface height (SSH) at current  time step [m] (internal mode)
    call allocate_block2D_real8(ubrtr, 0.0d0)   !barotropic velocity      zonal[m/s] at current  time step [m] (internal mode)
    call allocate_block2D_real8(vbrtr, 0.0d0)   !barotropic velocity meridional[m/s] at current  time step [m] (internal mode)
    call allocate_block2D_real8(RHSx2d, 0.0d0)  !x-component of external force(barotropic)
    call allocate_block2D_real8(RHSy2d, 0.0d0)  !y-component of external force(barotropic)

    call allocate_block2D_real8(sshn, 0.0d0)
    call allocate_block2D_real8(sshp, 0.0d0)
    call allocate_block2D_real8(ubrtrn, 0.0d0)
    call allocate_block2D_real8(ubrtrp, 0.0d0)
    call allocate_block2D_real8(vbrtrn, 0.0d0)
    call allocate_block2D_real8(vbrtrp, 0.0d0)

    call allocate_block2D_real8(wf_tot, 0.0d0) !total water flux

    call allocate_block2D_real8(BottomFriction, 0.0d0)
    call allocate_block2D_real8(r_diss, 0.0d0)

    call allocate_block2D_real8(amuv2d, 0.0d0)    !depth mean lateral viscosity
    call allocate_block2D_real8(amuv42d, 0.0d0)   !depth mean lateral viscosity
    call allocate_block2D_real8(r_vort2d, 0.0d0)  !relative vorticity of depth mean velocity
    call allocate_block2D_real8(stress_t2d, 0.0d0)          !Horizontal tension tensor component (barotropic)
    call allocate_block2D_real8(stress_s2d, 0.0d0)          !Horizontal shearing tensor component(barotropic)
    call allocate_block2D_real8(RHSx2d_tran_disp, 0.0d0)    !dispersion x-component of external force(barotropic)
    call allocate_block2D_real8(RHSy2d_tran_disp, 0.0d0)    !dispersion x-component of external force(barotropic)
    call allocate_block2D_real8(RHSx2d_diff_disp, 0.0d0)    !dispersion x-component of external force(barotropic)
    call allocate_block2D_real8(RHSy2d_diff_disp, 0.0d0)    !dispersion x-component of external force(barotropic)
    call allocate_block2D_real8(RHSx2d_bfc, 0.0d0)
    call allocate_block2D_real8(RHSy2d_bfc, 0.0d0)

endsubroutine ocean_variables_allocate
!-------------------------------------------------------------------------------

!deallocation of arrays
subroutine ocean_variables_deallocate
    use mpi_parallel_tools
    use ocean_variables
    implicit none

    call deallocate_block2D_real8(ssh)
    call deallocate_block2D_real8(ubrtr)
    call deallocate_block2D_real8(vbrtr)
    call deallocate_block2D_real8(RHSx2d)
    call deallocate_block2D_real8(RHSy2d)

    call deallocate_block2D_real8(sshn)
    call deallocate_block2D_real8(sshp)
    call deallocate_block2D_real8(ubrtrn)
    call deallocate_block2D_real8(ubrtrp)
    call deallocate_block2D_real8(vbrtrn)
    call deallocate_block2D_real8(vbrtrp)

    call deallocate_block2D_real8(wf_tot)

    call deallocate_block2D_real8(BottomFriction)
    call deallocate_block2D_real8(r_diss)

    call deallocate_block2D_real8(amuv2d)
    call deallocate_block2D_real8(amuv42d)
    call deallocate_block2D_real8(r_vort2d)
    call deallocate_block2D_real8(stress_t2d)
    call deallocate_block2D_real8(stress_s2d)
    call deallocate_block2D_real8(RHSx2d_tran_disp)
    call deallocate_block2D_real8(RHSy2d_tran_disp)
    call deallocate_block2D_real8(RHSx2d_diff_disp)
    call deallocate_block2D_real8(RHSy2d_diff_disp)
    call deallocate_block2D_real8(RHSx2d_bfc)
    call deallocate_block2D_real8(RHSy2d_bfc)

endsubroutine ocean_variables_deallocate
!-------------------------------------------------------------------------------
