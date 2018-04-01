!--------------------- SUBROUTINE FOR: -----------------------------------------!
!----------------- Only shallow water solving ----------------------------------!
!-------------------------------------------------------------------------------!
subroutine shallow_water_model_step(tau)
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables
    use shallow_water

    use time_integration

    implicit none
    integer :: m, n, k, ierr
    real*8 :: tau, diffslpr, surf_stress
    real*8 :: time_count

    diffslpr = 1.0d0
    surf_stress = 0.0d0

!---------------------- Shallow water equ solver -------------------------------
    !call start_timer(time_count)
    do k = 1, bcount
        call set_block_boundary(k)
        wf_tot(k)%vals = 0.0d0
        do n=ny_start,ny_end
            do m=nx_start,nx_end
                if (lu(k)%vals(m,n) > 0.5) then
                    RHSx2d(k)%vals(m, n) = (surf_stress)*dxt(k)%vals(m,n)*dyh(k)%vals(m,n) &
                        - (diffslpr)*hhu(k)%vals(m,n)*dyh(k)%vals(m,n)/RefDen
                    RHSy2d(k)%vals(m, n) = (surf_stress)*dyt(k)%vals(m,n)*dxh(k)%vals(m,n) &
                        - (diffslpr)*hhv(k)%vals(m,n)*dxh(k)%vals(m,n)/RefDen
                endif
            enddo
        enddo
    enddo

    do k = 1, bcount
        amuv2d(k)%vals  = lvisc_2
        amuv42d(k)%vals = lvisc_4
    enddo

    call expl_shallow_water(tau,     &
                          ksw_lat4,  &
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
                            r_diss,  &
                  RHSx2d_tran_disp,  &
                  RHSy2d_tran_disp,  &
                  RHSx2d_diff_disp,  &
                  RHSy2d_diff_disp,  &
                  RHSx2d_bfc,        &
                  RHSy2d_bfc)
    !call end_timer(time_count)
    !time_barotrop = time_barotrop + time_count

    ! Check errors
    do k = 1, bcount
        call set_block_boundary(k)
        do n=ny_start,ny_end
            do m=nx_start,nx_end
                if (lu(k)%vals(m,n) > 0.5) then
                    if (ssh(k)%vals(m,n) < 10000.0d0 .and. ssh(k)%vals(m,n) > -10000.0d0) then
                        continue
                    else
                        write(*,*) rank, 'ERR: Block k=', k, 'In the point m=', m, 'n=', n, 'ssh=', ssh(k)%vals(m,n),   &
                            'step: ', num_step, 'lon: ', geo_lon_t(k)%vals(m, n), 'lat: ', geo_lat_t(k)%vals(m, n)
                        call parallel_finalize; stop
                    endif
                endif
            enddo
        enddo
    enddo
endsubroutine shallow_water_model_step
