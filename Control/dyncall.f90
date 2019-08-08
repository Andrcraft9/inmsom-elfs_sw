!--------------------- SUBROUTINE FOR: -----------------------------------------!
!----------------- Only shallow water solving ----------------------------------!
!-------------------------------------------------------------------------------!
subroutine shallow_water_model_step(tau)
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables
    use shallow_water
    use flux_routes
    use time_integration
    use key_switches

    implicit none
    integer :: m, n, k, ierr
    real*8 :: tau
    real*8 :: time_count

    !surf_stress_x = 0.0d0
    !surf_stress_y = 0.0d0

!---------------------- Shallow water equ solver -------------------------------
    !call start_timer(time_count)
    if (ksw_atmforc == 1) then
        !Computation of sea surface boundary conditions
        if(ksw_ssbc > 0) then
            call sea_surface_fluxes
            !call sea_surface_fluxes_simple
        endif
        !Computing bottom stresses
        !if(type_fric>0) then
        !    call sea_bottom_fluxes
        !endif
        do k = 1, bcount
            call set_block(k)
            call set_block_lu(k)
            call set_block_h(k)
            call set_block_dxdy(k)
            call compute_fluxes_rhs(surf_stress_x(k)%vals, surf_stress_y(k)%vals, slpr(k)%vals, RHSx2d(k)%vals, RHSy2d(k)%vals)
        enddo
    endif

    do k = 1, bcount
        call set_block(k)
        amuv2d(k)%vals  = lvisc_2
        amuv42d(k)%vals = lvisc_4
    enddo

    call expl_shallow_water(tau,     &
                            ubrtr,   &
                            ubrtrp,  &
                            ubrtrn,  &
                            vbrtr,   &
                            vbrtrp,  &
                            vbrtrn,  &
                            ssh,     &
                            sshp,    &
                            sshn,    &
                            RHSx2d,  &
                            RHSy2d,  &
                            wf_tot,  &
                            amuv2d,  &
                           amuv42d,  &
                          r_vort2d,  &
                        stress_t2d,  &
                        stress_s2d,  &
                               xxt,  &
                               yyt,  &
                            r_diss,  &
                  RHSx2d_tran_disp,  &
                  RHSy2d_tran_disp,  &
                  RHSx2d_diff_disp,  &
                  RHSy2d_diff_disp,  &
                  RHSx2d_bfc,        &
                  RHSy2d_bfc)
    !call end_timer(time_count)
    !time_barotrop = time_barotrop + time_count

    ! Compute maximum amplitude
    do k = 1, bcount
        call set_block(k)
        call set_block_lu(k)
        call compute_max_amplitude(ssh(k)%vals, ssh_max_amplitude(k)%vals, lu)
        call compute_max_amplitude(ubrtr(k)%vals, ubrtr_max_amplitude(k)%vals, lcu)
        call compute_max_amplitude(vbrtr(k)%vals, vbrtr_max_amplitude(k)%vals, lcv)
    enddo

    ! Check errors
    do k = 1, bcount
        call set_block(k)
        call set_block_lu(k)
        call check_ssh_err(ssh(k)%vals)
    enddo

endsubroutine shallow_water_model_step

subroutine compute_max_amplitude(var, var_max_amplitude, varlu)
    implicit none
    real*8 :: var(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
              var_max_amplitude(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
    real*4 :: varlu(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
    integer :: m, n

    do n = ny_start, ny_end
        do m = nx_start, nx_end
            if (varlu(m, n)>0.5) then
                if ( var_max_amplitude(m, n) < abs(var(m, n)) ) then
                    var_max_amplitude(m, n) = abs(var(m, n))
                endif
            endif
        enddo
    enddo
end subroutine


subroutine check_ssh_err(ssh)
    use mpi_parallel_tools
    implicit none
    real*8 :: ssh(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
    integer :: m, n, ierr

    do n = ny_start,ny_end
        do m = nx_start,nx_end
            if (lu(m,n)>0.5) then
                if (ssh(m,n)<10000.0d0 .and. ssh(m,n)>-10000.0d0) then
                  continue
                else
                    write(*,*) rank, 'ERROR!!! In the point m=', m, 'n=', n, 'ssh=', ssh(m,n)
                    !write(*,*) rank, 'ERR: Block k=', k, 'In the point m=', m, 'n=', n, 'ssh=', ssh(k)%vals(m,n),   &
                    !    'step: ', num_step, 'lon: ', geo_lon_t(k)%vals(m, n), 'lat: ', geo_lat_t(k)%vals(m, n)

                    call mpi_abort(cart_comm, 1, ierr)
                    stop
                endif
            endif
        enddo
    enddo
end subroutine
