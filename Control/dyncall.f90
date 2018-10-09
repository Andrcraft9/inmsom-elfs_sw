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

    diffslpr = 1000.0d0
    surf_stress = 0.0d0

!---------------------- Shallow water equ solver -------------------------------
    !call start_timer(time_count)
    do k = 1, bcount
        call set_block_boundary(k)

        wf_tot(k)%vals = 0.0d0
        call compute_rhs(diffslpr, surf_stress, RHSx2d(k)%vals, RHSy2d(k)%vals,   &
                            hhu(k)%vals, hhv(k)%vals,                             &
                            dxt(k)%vals, dyt(k)%vals, dxh(k)%vals, dyh(k)%vals, lu(k)%vals)
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
        call check_ssh_err(ssh(k)%vals, lu(k)%vals)
    enddo

endsubroutine shallow_water_model_step


subroutine check_ssh_err(ssh, lu)
    use mpi_parallel_tools
    implicit none

    real*8 :: ssh(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
              lu(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

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
