module init_arrays_routes
    implicit none

contains

subroutine model_grid_allocate
    use mpi_parallel_tools
    use basin_grid
    implicit none

    ! Mask of entire basin area
    call allocate_block2D_real4(block_lu, 0.0)
    call allocate_block2D_real4(block_lu1, 0.0)
    call allocate_block2D_real4(block_luu, 0.0)
    call allocate_block2D_real4(block_luh, 0.0)
    call allocate_block2D_real4(block_lcu, 0.0)
    call allocate_block2D_real4(block_lcv, 0.0)
    call allocate_block2D_real4(block_llu, 0.0)
    call allocate_block2D_real4(block_llv, 0.0)

    ! Mask of local basin area (regional)
    call allocate_block2D_real4(block_lu_output, 0.0)
    call allocate_block2D_real4(block_lcu_output, 0.0)
    call allocate_block2D_real4(block_lcv_output, 0.0)
    call allocate_block2D_real4(block_llu_output, 0.0)
    call allocate_block2D_real4(block_llv_output, 0.0)

    ! Height of sea level
    call allocate_block2D_real8(block_hhh, 0.0d0)
    call allocate_block2D_real8(block_hhhp, 0.0d0)
    call allocate_block2D_real8(block_hhhn, 0.0d0)
    call allocate_block2D_real8(block_hhq_rest, 0.0d0)
    call allocate_block2D_real8(block_hhq, 0.0d0)
    call allocate_block2D_real8(block_hhqp, 0.0d0)
    call allocate_block2D_real8(block_hhqn, 0.0d0)
    call allocate_block2D_real8(block_hhu, 0.0d0)
    call allocate_block2D_real8(block_hhup, 0.0d0)
    call allocate_block2D_real8(block_hhun, 0.0d0)
    call allocate_block2D_real8(block_hhv, 0.0d0)
    call allocate_block2D_real8(block_hhvp, 0.0d0)
    call allocate_block2D_real8(block_hhvn, 0.0d0)

    ! Coriolis
    call allocate_block2D_real8(block_rlh_s, 0.0d0)
    call allocate_block2D_real8(block_rlh_c, 0.0d0)

    ! Grid steps
    call allocate_block2D_real8(block_dxt, 0.0d0)
    call allocate_block2D_real8(block_dyt, 0.0d0)
    call allocate_block2D_real8(block_dx, 0.0d0)
    call allocate_block2D_real8(block_dy, 0.0d0)
    call allocate_block2D_real8(block_dxh, 0.0d0)
    call allocate_block2D_real8(block_dyh, 0.0d0)
    call allocate_block2D_real8(block_dxb, 0.0d0)
    call allocate_block2D_real8(block_dyb, 0.0d0)

    ! Grid points
    ! Allocate 1D arrays with specific boundaries for each block
    allocate(block_xt(bcount))
    allocate(block_yt(bcount)) !horizontal t-grid            x- and y-coordinates (in degrees)
    allocate(block_xu(bcount))
    allocate(block_yv(bcount)) !horizontal u-grid and v-grid x- and y-coordinates (in degrees)
    do k = 1, bcount
       call set_block(k)
       allocate(block_xt(k)%vals(bnd_x1:bnd_x2))
       allocate(block_yt(k)%vals(bnd_y1:bnd_y2))
       allocate(block_xu(k)%vals(bnd_x1:bnd_x2))
       allocate(block_yv(k)%vals(bnd_y1:bnd_y2))
    enddo

    ! Geo coordinates
    call allocate_block2D_real8(block_geo_lon_t, 0.0d0)    !geographical longitudes of T-points
    call allocate_block2D_real8(block_geo_lat_t, 0.0d0)    !geographical latitudes  of T-points
    call allocate_block2D_real8(block_geo_lon_u, 0.0d0)    !geographical longitudes of U-points
    call allocate_block2D_real8(block_geo_lat_u, 0.0d0)    !geographical latitudes  of U-points
    call allocate_block2D_real8(block_geo_lon_v, 0.0d0)    !geographical longitudes of V-points
    call allocate_block2D_real8(block_geo_lat_v, 0.0d0)    !geographical latitudes  of V-points
    call allocate_block2D_real8(block_geo_lon_h, 0.0d0)    !geographical longitudes of H-points
    call allocate_block2D_real8(block_geo_lat_h, 0.0d0)    !geographical latitudes  of H-points

    ! Rotvec coeff
    call allocate_block3D_real8(block_rotvec_coeff, 4, 0.0d0)

    ! Set zero pointers to data
    do k = 1, bcount
       call set_block_zero(k)
    enddo

endsubroutine model_grid_allocate
!-------------------------------------------------------------------------------

subroutine model_grid_deallocate
    use basin_grid
    use mpi_parallel_tools
    implicit none

    ! Mask of entire basin area
    call deallocate_block2D_real4(block_lu)
    call deallocate_block2D_real4(block_lu1)
    call deallocate_block2D_real4(block_luu)
    call deallocate_block2D_real4(block_luh)
    call deallocate_block2D_real4(block_lcu)
    call deallocate_block2D_real4(block_lcv)
    call deallocate_block2D_real4(block_llu)
    call deallocate_block2D_real4(block_llv)

    ! Mask of local basin area (regional)
    call deallocate_block2D_real4(block_lu_output)
    call deallocate_block2D_real4(block_lcu_output)
    call deallocate_block2D_real4(block_lcv_output)
    call deallocate_block2D_real4(block_llu_output)
    call deallocate_block2D_real4(block_llv_output)

    ! Height of sea level
    call deallocate_block2D_real8(block_hhh)
    call deallocate_block2D_real8(block_hhhp)
    call deallocate_block2D_real8(block_hhhn)
    call deallocate_block2D_real8(block_hhq_rest)
    call deallocate_block2D_real8(block_hhq)
    call deallocate_block2D_real8(block_hhqp)
    call deallocate_block2D_real8(block_hhqn)
    call deallocate_block2D_real8(block_hhu)
    call deallocate_block2D_real8(block_hhup)
    call deallocate_block2D_real8(block_hhun)
    call deallocate_block2D_real8(block_hhv)
    call deallocate_block2D_real8(block_hhvp)
    call deallocate_block2D_real8(block_hhvn)

    ! Coriolis
    call deallocate_block2D_real8(block_rlh_s)
    call deallocate_block2D_real8(block_rlh_c)

    ! Grid steps
    call deallocate_block2D_real8(block_dxt)
    call deallocate_block2D_real8(block_dyt)
    call deallocate_block2D_real8(block_dx)
    call deallocate_block2D_real8(block_dy)
    call deallocate_block2D_real8(block_dxh)
    call deallocate_block2D_real8(block_dyh)
    call deallocate_block2D_real8(block_dxb)
    call deallocate_block2D_real8(block_dyb)
    ! Grid points
    call deallocate_block1D_real8(block_xt)
    call deallocate_block1D_real8(block_yt)
    call deallocate_block1D_real8(block_xu)
    call deallocate_block1D_real8(block_yv)

    ! Geo coordinates
    call deallocate_block2D_real8(block_geo_lon_t)
    call deallocate_block2D_real8(block_geo_lat_t)
    call deallocate_block2D_real8(block_geo_lon_u)
    call deallocate_block2D_real8(block_geo_lat_u)
    call deallocate_block2D_real8(block_geo_lon_v)
    call deallocate_block2D_real8(block_geo_lat_v)
    call deallocate_block2D_real8(block_geo_lon_h)
    call deallocate_block2D_real8(block_geo_lat_h)

endsubroutine model_grid_deallocate
!-------------------------------------------------------------------------------

subroutine ocean_variables_allocate
    use mpi_parallel_tools
    use ocean_variables
    implicit none

    ! Barotropic dynamics arrays
    call allocate_block2D_real8(ssh, 0.0d0)
    call allocate_block2D_real8(ubrtr, 0.0d0)
    call allocate_block2D_real8(vbrtr, 0.0d0)
    call allocate_block2D_real8(xxt, 0.0d0)
    call allocate_block2D_real8(yyt, 0.0d0)
    call allocate_block2D_real8(RHSx2d, 0.0d0)
    call allocate_block2D_real8(RHSy2d, 0.0d0)

    call allocate_block2D_real8(sshn, 0.0d0)
    call allocate_block2D_real8(sshp, 0.0d0)
    call allocate_block2D_real8(ubrtrn, 0.0d0)
    call allocate_block2D_real8(ubrtrp, 0.0d0)
    call allocate_block2D_real8(vbrtrn, 0.0d0)
    call allocate_block2D_real8(vbrtrp, 0.0d0)

    ! Sea surface boundary condition
    call allocate_block2D_real4(tflux_surf, 0.0d0)
    call allocate_block2D_real4(tflux_bot, 0.0d0)
    call allocate_block2D_real4(sflux_surf, 0.0d0)
    call allocate_block2D_real4(sflux_bot, 0.0d0)
    call allocate_block2D_real4(surf_stress_x, 0.0d0)
    call allocate_block2D_real4(surf_stress_y, 0.0d0)
    call allocate_block2D_real4(bot_stress_x, 0.0d0)
    call allocate_block2D_real4(bot_stress_y, 0.0d0)
    call allocate_block2D_real4(dkft, 0.0d0)
    call allocate_block2D_real4(dkfs, 0.0d0)
    call allocate_block2D_real4(sensheat, 0.0d0)
    call allocate_block2D_real4(latheat, 0.0d0)
    call allocate_block2D_real4(lw_bal, 0.0d0)
    call allocate_block2D_real4(sw_bal, 0.0d0)
    call allocate_block2D_real4(hf_tot, 0.0d0)
    call allocate_block2D_real4(wf_tot, 0.0d0)

    ! Atmospheric arrays for bulk-formulae
    call allocate_block2D_real4(tatm, 0.0)
    call allocate_block2D_real4(qatm, 0.0)
    call allocate_block2D_real4(rain, 0.0)
    call allocate_block2D_real4(snow, 0.0)
    call allocate_block2D_real4(evap, 0.0)
    call allocate_block2D_real4(wind, 0.0)
    call allocate_block2D_real4( lwr, 0.0)
    call allocate_block2D_real4( swr, 0.0)
    call allocate_block2D_real4(slpr, 0.0)
    call allocate_block2D_real4(uwnd, 0.0)
    call allocate_block2D_real4(vwnd, 0.0)
    call allocate_block2D_real4(taux, 0.0)
    call allocate_block2D_real4(tauy, 0.0)

    ! Rayleigh dissipation
    call allocate_block2D_real8(r_diss, 0.0d0)

    call allocate_block2D_real8(amuv2d, 0.0d0)
    call allocate_block2D_real8(amuv42d, 0.0d0)
    call allocate_block2D_real8(r_vort2d, 0.0d0)
    call allocate_block2D_real8(stress_t2d, 0.0d0)
    call allocate_block2D_real8(stress_s2d, 0.0d0)
    call allocate_block2D_real8(RHSx2d_tran_disp, 0.0d0)
    call allocate_block2D_real8(RHSy2d_tran_disp, 0.0d0)
    call allocate_block2D_real8(RHSx2d_diff_disp, 0.0d0)
    call allocate_block2D_real8(RHSy2d_diff_disp, 0.0d0)
    call allocate_block2D_real8(RHSx2d_bfc, 0.0d0)
    call allocate_block2D_real8(RHSy2d_bfc, 0.0d0)

    ! Maximum amplitutes
    call allocate_block2D_real8(ssh_max_amplitude, 0.0d0)
    call allocate_block2D_real8(ubrtr_max_amplitude, 0.0d0)
    call allocate_block2D_real8(vbrtr_max_amplitude, 0.0d0)

    ! Temporary array
    call allocate_block2D_real4(array4_2d, 0.0)

endsubroutine ocean_variables_allocate
!-------------------------------------------------------------------------------

subroutine ocean_variables_deallocate
    use mpi_parallel_tools
    use ocean_variables
    implicit none

    ! Barotropic dynamics arrays
    call deallocate_block2D_real8(ssh)
    call deallocate_block2D_real8(ubrtr)
    call deallocate_block2D_real8(vbrtr)
    call deallocate_block2D_real8(xxt)
    call deallocate_block2D_real8(yyt)
    call deallocate_block2D_real8(RHSx2d)
    call deallocate_block2D_real8(RHSy2d)

    call deallocate_block2D_real8(sshn)
    call deallocate_block2D_real8(sshp)
    call deallocate_block2D_real8(ubrtrn)
    call deallocate_block2D_real8(ubrtrp)
    call deallocate_block2D_real8(vbrtrn)
    call deallocate_block2D_real8(vbrtrp)

    ! Sea surface boundary condition
    call deallocate_block2D_real4(tflux_surf)
    call deallocate_block2D_real4(tflux_bot)
    call deallocate_block2D_real4(sflux_surf)
    call deallocate_block2D_real4(sflux_bot)
    call deallocate_block2D_real4(surf_stress_x)
    call deallocate_block2D_real4(surf_stress_y)
    call deallocate_block2D_real4(bot_stress_x)
    call deallocate_block2D_real4(bot_stress_y)
    call deallocate_block2D_real4(dkft)
    call deallocate_block2D_real4(dkfs)
    call deallocate_block2D_real4(sensheat)
    call deallocate_block2D_real4(latheat)
    call deallocate_block2D_real4(lw_bal)
    call deallocate_block2D_real4(sw_bal)
    call deallocate_block2D_real4(hf_tot)
    call deallocate_block2D_real4(wf_tot)

    ! Atmospheric arrays for bulk-formulae
    call deallocate_block2D_real4(tatm)
    call deallocate_block2D_real4(qatm)
    call deallocate_block2D_real4(rain)
    call deallocate_block2D_real4(snow)
    call deallocate_block2D_real4(evap)
    call deallocate_block2D_real4(wind)
    call deallocate_block2D_real4( lwr)
    call deallocate_block2D_real4( swr)
    call deallocate_block2D_real4(slpr)
    call deallocate_block2D_real4(uwnd)
    call deallocate_block2D_real4(vwnd)
    call deallocate_block2D_real4(taux)
    call deallocate_block2D_real4(tauy)

    ! Rayleigh dissipation
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

    ! Maximum amplitutes
    call deallocate_block2D_real8(ssh_max_amplitude)
    call deallocate_block2D_real8(ubrtr_max_amplitude)
    call deallocate_block2D_real8(vbrtr_max_amplitude)

    ! Temporary array
    call deallocate_block2D_real4(array4_2d)

endsubroutine ocean_variables_deallocate
!-------------------------------------------------------------------------------

end module
