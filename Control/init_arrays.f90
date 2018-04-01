!allocation of arrays
subroutine model_grid_allocate
    use mpi_parallel_tools
    use basin_grid
    implicit none

    call allocate_block2D(lu, 0.0d0)  !mask of t-grid
    call allocate_block2D(lu1, 0.0d0) !mask of t-grid (1 everywhere)
    call allocate_block2D(luu, 0.0d0) !mask of h-grid (0 on boundary)
    call allocate_block2D(luh, 0.0d0) !mask of h-grid (1 on boundary)
    call allocate_block2D(lcu, 0.0d0) !mask of u-grid (0 on boundary)
    call allocate_block2D(lcv, 0.0d0) !mask of v-grid (0 on boundary)
    call allocate_block2D(llu, 0.0d0) !mask of u-grid (0 on boundary)
    call allocate_block2D(llv, 0.0d0) !mask of v-grid (0 on boundary)

    call allocate_block2D(hhh, 0.0d0)   !ocean depth on luh (h-points)
    call allocate_block2D(hhhp, 0.0d0)  !ocean depth on luh (h-points) at previous step
    call allocate_block2D(hhhn, 0.0d0)  !ocean depth on luh (h-points) at pre-previous step
    call allocate_block2D(hhq_rest, 0.0d0)  !ocean depth on lu  (t-points) at rest state
    call allocate_block2D(hhq, 0.0d0)       !ocean depth on lu  (t-points)
    call allocate_block2D(hhqp, 0.0d0)  !ocean depth on lu  (t-points) at previous step
    call allocate_block2D(hhqn, 0.0d0)  !ocean depth on lu  (t-points) at pre-previous step
    call allocate_block2D(hhu, 0.0d0)   !ocean depth on lcu (u-points)
    call allocate_block2D(hhup, 0.0d0)  !ocean depth on lcu (u-points) at previous step
    call allocate_block2D(hhun, 0.0d0)  !ocean depth on lcu (u-points) at pre-previous step
    call allocate_block2D(hhv, 0.0d0)   !ocean depth on lcv (v-points)
    call allocate_block2D(hhvp, 0.0d0)  !ocean depth on lcv (v-points) at previous step
    call allocate_block2D(hhvn, 0.0d0)  !ocean depth on lcv (v-points) at pre-previous step
    call allocate_block2D(rlh_s, 0.0d0)  !coriolis-1 parameter on edge (t-centers) points
    call allocate_block2D(rlh_c, 0.0d0)  !coriolis-2 parameter on edge (t-centers) points

    call allocate_block2D(dxt, 0.0d0)
    call allocate_block2D(dyt, 0.0d0)  !horizontal grid steps between   t-points (in radians or meters)
    call allocate_block2D(dx, 0.0d0)
    call allocate_block2D(dy, 0.0d0)   !horizontal grid steps between u,v-points (in radians or meters)
    call allocate_block2D(dxh, 0.0d0)
    call allocate_block2D(dyh, 0.0d0)  !horizontal grid steps between   h-points (in radians or meters)
    call allocate_block2D(dxb, 0.0d0)
    call allocate_block2D(dyb, 0.0d0)  !horizontal grid steps between v,u-points (in radians or meters)
    call allocate_block2D(xt, 0.0d0)
    call allocate_block2D(yt, 0.0d0)   !horizontal t-grid            x- and y-coordinates (in degrees)
    call allocate_block2D(xu, 0.0d0)
    call allocate_block2D(yv, 0.0d0)   !horizontal u-grid and v-grid x- and y-coordinates (in degrees)

    call allocate_block2D(geo_lon_t, 0.0d0)    !geographical longitudes of T-points
    call allocate_block2D(geo_lat_t, 0.0d0)    !geographical latitudes  of T-points
    call allocate_block2D(geo_lon_u, 0.0d0)    !geographical longitudes of U-points
    call allocate_block2D(geo_lat_u, 0.0d0)    !geographical latitudes  of U-points
    call allocate_block2D(geo_lon_v, 0.0d0)    !geographical longitudes of V-points
    call allocate_block2D(geo_lat_v, 0.0d0)    !geographical latitudes  of V-points
    call allocate_block2D(geo_lon_h, 0.0d0)    !geographical longitudes of H-points
    call allocate_block2D(geo_lat_h, 0.0d0)    !geographical latitudes  of H-points

endsubroutine model_grid_allocate
!-------------------------------------------------------------------------------

!deallocation of arrays
subroutine model_grid_deallocate
    use basin_grid
    use mpi_parallel_tools
    implicit none

    call deallocate_blocks2D(lu)
    call deallocate_blocks2D(lu1)
    call deallocate_blocks2D(luu)
    call deallocate_blocks2D(luh)
    call deallocate_blocks2D(lcu)
    call deallocate_blocks2D(lcv)
    call deallocate_blocks2D(llu)
    call deallocate_blocks2D(llv)

    call deallocate_blocks2D(hhh)
    call deallocate_blocks2D(hhhp)
    call deallocate_blocks2D(hhhn)
    call deallocate_blocks2D(hhq_rest)
    call deallocate_blocks2D(hhq)
    call deallocate_blocks2D(hhqp)
    call deallocate_blocks2D(hhqn)
    call deallocate_blocks2D(hhu)
    call deallocate_blocks2D(hhup)
    call deallocate_blocks2D(hhun)
    call deallocate_blocks2D(hhv)
    call deallocate_blocks2D(hhvp)
    call deallocate_blocks2D(hhvn)
    call deallocate_blocks2D(rlh_s)
    call deallocate_blocks2D(rlh_c)

    call deallocate_blocks2D(dxt)
    call deallocate_blocks2D(dyt)
    call deallocate_blocks2D(dx)
    call deallocate_blocks2D(dy)
    call deallocate_blocks2D(dxh)
    call deallocate_blocks2D(dyh)
    call deallocate_blocks2D(dxb)
    call deallocate_blocks2D(dyb)
    call deallocate_blocks2D(xt)
    call deallocate_blocks2D(yt)
    call deallocate_blocks2D(xu)
    call deallocate_blocks2D(yv)

    call deallocate_blocks2D(geo_lon_t)
    call deallocate_blocks2D(geo_lat_t)
    call deallocate_blocks2D(geo_lon_u)
    call deallocate_blocks2D(geo_lat_u)
    call deallocate_blocks2D(geo_lon_v)
    call deallocate_blocks2D(geo_lat_v)
    call deallocate_blocks2D(geo_lon_h)
    call deallocate_blocks2D(geo_lat_h)

endsubroutine model_grid_deallocate
!-------------------------------------------------------------------------------

!allocation of arrays
subroutine ocean_variables_allocate
    use mpi_parallel_tools
    use ocean_variables
    implicit none

    call allocate_block2D(ssh, 0.0d0)     !sea surface height (SSH) at current  time step [m] (internal mode)
    call allocate_block2D(ubrtr, 0.0d0)   !barotropic velocity      zonal[m/s] at current  time step [m] (internal mode)
    call allocate_block2D(vbrtr, 0.0d0)   !barotropic velocity meridional[m/s] at current  time step [m] (internal mode)
    call allocate_block2D(RHSx2d, 0.0d0)  !x-component of external force(barotropic)
    call allocate_block2D(RHSy2d, 0.0d0)  !y-component of external force(barotropic)

    call allocate_block2D(sshn, 0.0d0)
    call allocate_block2D(sshp, 0.0d0)
    call allocate_block2D(ubrtrn, 0.0d0)
    call allocate_block2D(ubrtrp, 0.0d0)
    call allocate_block2D(vbrtrn, 0.0d0)
    call allocate_block2D(vbrtrp, 0.0d0)

    call allocate_block2D(wf_tot, 0.0d0) !total water flux

    call allocate_block2D(BottomFriction, 0.0d0)
    call allocate_block2D(r_diss, 0.0d0)

    call allocate_block2D(amuv2d, 0.0d0)    !depth mean lateral viscosity
    call allocate_block2D(amuv42d, 0.0d0)   !depth mean lateral viscosity
    call allocate_block2D(r_vort2d, 0.0d0)  !relative vorticity of depth mean velocity
    call allocate_block2D(stress_t2d, 0.0d0)          !Horizontal tension tensor component (barotropic)
    call allocate_block2D(stress_s2d, 0.0d0)          !Horizontal shearing tensor component(barotropic)
    call allocate_block2D(RHSx2d_tran_disp, 0.0d0)    !dispersion x-component of external force(barotropic)
    call allocate_block2D(RHSy2d_tran_disp, 0.0d0)    !dispersion x-component of external force(barotropic)
    call allocate_block2D(RHSx2d_diff_disp, 0.0d0)    !dispersion x-component of external force(barotropic)
    call allocate_block2D(RHSy2d_diff_disp, 0.0d0)    !dispersion x-component of external force(barotropic)
    call allocate_block2D(RHSx2d_bfc, 0.0d0)
    call allocate_block2D(RHSy2d_bfc, 0.0d0)

endsubroutine ocean_variables_allocate
!-------------------------------------------------------------------------------

!deallocation of arrays
subroutine ocean_variables_deallocate
    use mpi_parallel_tools
    use ocean_variables
    implicit none

    call deallocate_blocks2D(ssh)
    call deallocate_blocks2D(ubrtr)
    call deallocate_blocks2D(vbrtr)
    call deallocate_blocks2D(RHSx2d)
    call deallocate_blocks2D(RHSy2d)

    call deallocate_blocks2D(sshn)
    call deallocate_blocks2D(sshp)
    call deallocate_blocks2D(ubrtrn)
    call deallocate_blocks2D(ubrtrp)
    call deallocate_blocks2D(vbrtrn)
    call deallocate_blocks2D(vbrtrp)

    call deallocate_blocks2D(wf_tot)

    call deallocate_blocks2D(BottomFriction)
    call deallocate_blocks2D(r_diss)

    call deallocate_blocks2D(amuv2d)
    call deallocate_blocks2D(amuv42d)
    call deallocate_blocks2D(r_vort2d)
    call deallocate_blocks2D(stress_t2d)
    call deallocate_blocks2D(stress_s2d)
    call deallocate_blocks2D(RHSx2d_tran_disp)
    call deallocate_blocks2D(RHSy2d_tran_disp)
    call deallocate_blocks2D(RHSx2d_diff_disp)
    call deallocate_blocks2D(RHSy2d_diff_disp)
    call deallocate_blocks2D(RHSx2d_bfc)
    call deallocate_blocks2D(RHSy2d_bfc)

endsubroutine ocean_variables_deallocate
!-------------------------------------------------------------------------------
