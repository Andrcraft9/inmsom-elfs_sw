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
    real*8 :: diffslpr, tau
    real*8 :: time_count

    !diffslpr = 0.0d0
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

        do n=ny_start,ny_end
            do m=nx_start,nx_end
                if(lu(m,n)>0.5) then
                    !RHSx2d(m, n) = ( surf_stress_x(m,n)+bot_stress_x(m,n) )*dxt(m,n)*dyh(m,n)    &
                    RHSx2d(m, n) = (surf_stress_x(m,n))*dxt(m,n)*dyh(m,n)    &
                             -(slpr(m+1,n)-slpr(m,n))*hhu(m,n)*dyh(m,n)/RefDen

                    !RHSy2d(m, n) = ( surf_stress_y(m,n)+bot_stress_y(m,n) )*dyt(m,n)*dxh(m,n)    &
                    RHSy2d(m, n) = (surf_stress_y(m,n))*dyt(m,n)*dxh(m,n)    &
                             -(slpr(m,n+1)-slpr(m,n))*hhv(m,n)*dxh(m,n)/RefDen
                endif
            enddo
        enddo
    else
        !wf_tot = 0.0d0
        !do n=ny_start,ny_end
        !    do m=nx_start,nx_end
        !        if(lu(m,n)>0.5) then
        !            RHSx2d(m, n) = (surf_stress_x(m,n))*dxt(m,n)*dyh(m,n) -(diffslpr)*hhu(m,n)*dyh(m,n)/RefDen
        !            RHSy2d(m, n) = (surf_stress_y(m,n))*dyt(m,n)*dxh(m,n) -(diffslpr)*hhv(m,n)*dxh(m,n)/RefDen
        !        endif
        !    enddo
        !enddo
    endif

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
!    do n = ny_start, ny_end
!        do m = nx_start, nx_end
!            if (lu(m, n)>0.5) then
!                if ( ssh_max_amplitude(m, n) < abs(ssh(m, n)) ) then
!                    ssh_max_amplitude(m, n) = abs(ssh(m, n))
!                endif
!            endif
!
!            if (lcu(m, n)>0.5) then
!                if ( ubrtr_max_amplitude(m, n) < abs(ubrtr(m, n)) ) then
!                    ubrtr_max_amplitude(m, n) = abs(ubrtr(m, n))
!                endif
!            endif
!
!            if (lcv(m, n)>0.5) then
!                if ( vbrtr_max_amplitude(m, n) < abs(vbrtr(m, n)) ) then
!                    vbrtr_max_amplitude(m, n) = abs(vbrtr(m, n))
!                endif
!            endif
!        enddo
!    enddo

    ! Check errors
    do n=ny_start,ny_end
      do m=nx_start,nx_end
          if(lu(m,n)>0.5) then
              if(ssh(m,n)<10000.0d0 .and. ssh(m,n)>-10000.0d0) then
                  continue
              else
                  write(*,*) rank, 'ERROR!!! In the point m=', m, 'n=', n, 'ssh=', ssh(m,n),   &
                    'step: ', num_step, 'lon: ', geo_lon_t(m, n), 'lat: ', geo_lat_t(m, n)

                  call mpi_abort(cart_comm, 1, ierr)
                  stop
              endif
          endif
      enddo
    enddo

endsubroutine shallow_water_model_step
