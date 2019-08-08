module flux_routes
implicit none

contains

subroutine sea_surface_fluxes_simple
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables
    use atm_forcing
    implicit none

    integer :: m, n, k
    real*8 :: wnd, wnd_mod
    real*8 :: coeff_surf_fric

    coeff_surf_fric = 3.25d-6

    do k = 1, bcount
        call set_block(k)
        call set_block_lu(k)
        call set_block_dxdy(k)
        do n=ny_start,ny_end
            do m=nx_start,nx_end

            if(lcu(m,n)>0.5) then
                wnd = (uwnd(k)%vals(m  ,n)*dx(m  ,n)*dy(m  ,n)                              &
                     + uwnd(k)%vals(m+1,n)*dx(m+1,n)*dy(m+1,n))/2.0d0/dxt(m,n)/dyh(m,n)
                wnd_mod = (wind(k)%vals(m  ,n)*dx(m  ,n)*dy(m  ,n)                          &
                         + wind(k)%vals(m+1,n)*dx(m+1,n)*dy(m+1,n))/2.0d0/dxt(m,n)/dyh(m,n)
                surf_stress_x(m,n) = wnd * wnd_mod * coeff_surf_fric
            endif

            if(lcv(m,n)>0.5) then
                wnd = (vwnd(k)%vals(m,n  )*dx(m,n  )*dy(m,n  )                              &
                     + vwnd(k)%vals(m,n+1)*dx(m,n+1)*dy(m,n+1))/2.0d0/dxh(m,n)/dyt(m,n)
                wnd_mod = (wind(k)%vals(m,n  )*dx(m,n  )*dy(m,n  )                          &
                         + wind(k)%vals(m,n+1)*dx(m,n+1)*dy(m,n+1))/2.0d0/dxh(m,n)/dyt(m,n)
                surf_stress_y(m,n) = wnd * wnd_mod * coeff_surf_fric

            endif
            enddo
        enddo
    enddo

endsubroutine sea_surface_fluxes_simple

subroutine sea_surface_fluxes
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables
    use atm_forcing
    use ocalg_routes
    implicit none

    integer m,n,k,ierr

    real(8) evap_rate
    real(8) wf_ave, sf_ave, m_calc, wf, wnd
    real(8) tmp_real8

    if(ksw_ssbc==1 .or. ksw_ssbc==2) then
        do k = 1, bcount
            call set_block(k)
            call set_block_lu(k)
            call set_block_dxdy(k)
            do n=ny_start,ny_end
                do m=nx_start,nx_end

                if(lcu(m,n)>0.5) then
                    surf_stress_x(m,n)= (taux(k)%vals(m  ,n)*dx(m  ,n)*dy(m  ,n)   &
                                        +taux(k)%vals(m+1,n)*dx(m+1,n)*dy(m+1,n))/RefDen/2.0d0/dxt(m,n)/dyh(m,n)
                endif

                if(lcv(m,n)>0.5) then
                    surf_stress_y(m,n)= (tauy(k)%vals(m,n  )*dx(m,n  )*dy(m,n  )   &
                                        +tauy(k)%vals(m,n+1)*dx(m,n+1)*dy(m,n+1))/RefDen/2.0d0/dxh(m,n)/dyt(m,n)
                endif

                enddo
            enddo
        enddo
    endif

    if(ksw_ssbc==3) then
        do k = 1, bcount
            call set_block(k)
            call set_block_lu(k)
            call set_block_dxdy(k)
            do n=ny_start,ny_end
                do m=nx_start,nx_end
                    if(lu(m,n)>0.5) then

                        call air_sea_turbulent_fluxes(  wind(k)%vals(m,n),         &   ! wind modulo, m/s
                                                        slpr(k)%vals(m,n),         &   ! sea level pressure, Pa
                                                        tatm(k)%vals(m,n),         &   ! air temperature,  �C
                                                        tatm(k)%vals(m,n),         &   ! sea surface temp, �C
                                                        qatm(k)%vals(m,n),         &   ! air specific humidity, kg/kg
                                                         u_height,                 &   ! height of wind datasets, m
                                                         t_height,                 &   ! height of tair datasets, m
                                                         q_height,                 &   ! height of qair datasets, m
                                                        taux(k)%vals(m,n),         &   !      zonal wind stress, Pa
                                                        tauy(k)%vals(m,n)          )   ! meridional wind stress, Pa

                        wf_tot(k)%vals(m,n) = rain(k)%vals(m,n) + snow(k)%vals(m,n)

                    endif
                enddo
            enddo
        enddo

        call syncborder_block2D_real8(taux)
        call syncborder_block2D_real8(tauy)
        if(periodicity_x/=0) then
            call cyclize_x_block2D_real8(taux)
        endif
        if(periodicity_y/=0) then
            call cyclize_y_block2D_real8(tauy)
        endif

        wf_ave=0.0d0
        sf_ave=0.0d0
        m_calc=0.0d0
        if(ksw_wflux>0) then
            do k = 1, bcount
                call set_block(k)
                call set_block_lu(k)
                call set_block_dxdy(k)
                do n=ny_start,ny_end
                    do m=nx_start,nx_end
                        if(lu(m,n)>0.5) then
                            wf_ave = wf_ave + wf_tot(k)%vals(m,n)*dx(m,n)*dy(m,n)
                            m_calc = m_calc + dx(m,n)*dy(m,n)
                        endif
                    enddo
                enddo
            enddo
            tmp_real8 = wf_ave
            call mpi_allreduce(tmp_real8, wf_ave, 1, mpi_real8, mpi_sum, cart_comm, ierr)
            tmp_real8 = m_calc
            call mpi_allreduce(tmp_real8, m_calc, 1, mpi_real8, mpi_sum, cart_comm, ierr)

            wf_tot = wf_tot - wf_ave/m_calc
        endif

        call syncborder_block2D_real8(wf_tot)
        if(periodicity_x/=0) then
            call cyclize_x_block2D_real8(wf_tot)
        endif
        if(periodicity_y/=0) then
            call cyclize_y_block2D_real8(wf_tot)
        endif

        do k = 1, bcount
            call set_block(k)
            call set_block_lu(k)
            call set_block_dxdy(k)
            do n=ny_start,ny_end
                do m=nx_start,nx_end
                    if(lcu(m,n)>0.5) then
                        wnd = (uwnd(k)%vals(m  ,n)*dx(m  ,n)*dy(m  ,n)         &
                             + uwnd(k)%vals(m+1,n)*dx(m+1,n)*dy(m+1,n))/2.0d0/dxt(m,n)/dyh(m,n)

                        wf = (wf_tot(k)%vals(m  ,n)*uwnd(k)%vals(m  ,n)*dx(m  ,n)*dy(m  ,n)         &
                            + wf_tot(k)%vals(m+1,n)*uwnd(k)%vals(m+1,n)*dx(m+1,n)*dy(m+1,n))/2.0d0/dxt(m,n)/dyh(m,n)

                        surf_stress_x(k)%vals(m,n) = (taux(k)%vals(m,n) + taux(k)%vals(m+1,n))/RefDen/2.0d0     &
                                                    *( wnd ) + wf/RefDen
                    endif

                    if(lcv(m,n)>0.5) then
                        wnd = (vwnd(k)%vals(m,n  )*dx(m,n  )*dy(m,n  )         &
                             + vwnd(k)%vals(m,n+1)*dx(m,n+1)*dy(m,n+1))/2.0d0/dxh(m,n)/dyt(m,n)

                        wf = (wf_tot(k)%vals(m,n  )*vwnd(k)%vals(m,n  )*dx(m,n  )*dy(m,n  )         &
                            + wf_tot(k)%vals(m,n+1)*vwnd(k)%vals(m,n+1)*dx(m,n+1)*dy(m,n+1))/2.0d0/dxh(m,n)/dyt(m,n)

                        surf_stress_y(k)%vals(m,n) = (tauy(k)%vals(m,n) + tauy(k)%vals(m,n+1))/RefDen/2.0d0     &
                                                    *( wnd ) + wf/RefDen
                    endif

                enddo
            enddo
        enddo

    endif

endsubroutine sea_surface_fluxes

endmodule flux_routes
