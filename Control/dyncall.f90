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
    real(8) tau, diffslpr
    real*8 :: time_count

    diffslpr = 1.0d+3
    !diffslpr = 0.0d0

!---------------------- Shallow water equ solver -------------------------------
    !call start_timer(time_count)
    if (atm_forcing_on == 1) then
        !Computation of sea surface boundary conditions
        if(ksw_ssbc > 0) then
            call sea_surface_fluxes
            !call sea_surface_fluxes_simple

            call syncborder_extra_real8(surf_stress_x, 1, bnd_length)
            call syncborder_extra_real8(surf_stress_y, 1, bnd_length)
        endif
        !Computing bottom stresses
        !if(type_fric>0) then
        !    call sea_bottom_fluxes
        !endif

        ! Need correct slpr arrays!...
        ! TO DO: forc_atm.f90 -> to LCHC mode!...
        ! ...
        do n=bnd_y1+2, bnd_y2-2
            do m=bnd_x1+2, bnd_x2-2
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
        wf_tot = 0.0d0
        do n=bnd_y1+2, bnd_y2-2
            do m=bnd_x1+2, bnd_x2-2
                if(lu(m,n)>0.5) then
                    RHSx2d(m, n) = -(diffslpr)*hhu(m,n)*dyh(m,n)/RefDen
                    RHSy2d(m, n) = -(diffslpr)*hhv(m,n)*dxh(m,n)/RefDen
                endif
            enddo
        enddo
    endif

    amuv2d  = lvisc_2
    amuv42d = lvisc_4

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

    ! Check errors
    do n=ny_start,ny_end
      do m=nx_start,nx_end
          if(lu(m,n)>0.5) then
              if(ssh(m,n)<10000.0d0.and.ssh(m,n)>-10000.0d0) then
                  continue
              else
                  write(*,*) rank, 'ERROR!!! In the point m=', m, 'n=', n, 'ssh=', ssh(m,n),   &
                    'step: ', num_step, 'lon: ', geo_lon_t(m, n), 'lat: ', geo_lat_t(m, n),    &
                    'BFCx: ', RHSx2d_bfc(m, n), 'BFCy: ', RHSy2d_bfc(m, n),      &
                    'hhh: ', hhh(m, n), hhh(m-1, n), hhh(m, n-1),                &
                    'dxb: ', dxb(m, n), dxb(m-1, n), dxb(m, n-1),                &
                    'dyb: ', dyb(m, n), dyb(m-1, n), dyb(m, n-1)

                  call mpi_finalize(ierr)
                  call exit(0)
              endif
          endif
      enddo
    enddo

endsubroutine shallow_water_model_step
