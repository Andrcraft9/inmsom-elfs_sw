module vel_ssh
    use main_basin_pars
    use mpi_parallel_tools
    implicit none

contains

    ! RHS for implicit bfc scheme
    subroutine uv_bfc(u, v, hq, hu, hv, hh, RHSx, RHSy, dxb, dyb, lcu, lcv)
        implicit none

        real*8 :: u(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &  !Transporting zonal velocity
                  v(bnd_x1:bnd_x2,bnd_y1:bnd_y2)           !Transporting meridional velocity

        real*8 :: RHSx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &
                  RHSy(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

        real*8 :: hq(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  hu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  hv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  hh(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

        real*8 :: dxb(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  dyb(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  lcu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  lcv(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

        integer :: m, n

        real*8 :: k_bfc, s
        real*8 :: k1, k2

        !$omp parallel do private(m, n, k_bfc, s, k1, k2)
        do n=ny_start, ny_end
            do m=nx_start, nx_end
                if (lcu(m,n)>0.5) then
                    ! Discretization in h-points
                    k_bfc = FreeFallAcc * (nbfc**2) / (hh(m, n)**(1.0/3.0))
                    s = 0.5d0 * sqrt( (u(m, n) + u(m, n+1))**2 + (v(m, n) + v(m+1, n))**2 )
                    k1 = -dxb(m, n) * dyb(m, n) * k_bfc * s
                    !k1 = k1 * 0.5d0*(u(m, n) + u(m, n+1))

                    ! Discretization in h-points
                    k_bfc = FreeFallAcc * (nbfc**2) / (hh(m, n-1)**(1.0/3.0))
                    s = 0.5d0 * sqrt( (u(m, n) + u(m, n-1))**2 + (v(m, n-1) + v(m+1, n-1))**2 )
                    k2 = -dxb(m, n-1) * dyb(m, n-1) * k_bfc * s
                    !k2 = k2 * 0.5d0*(u(m, n) + u(m, n-1))

                    ! Discretization in u-points
                    RHSx(m, n) = 0.5d0 * (k1 + k2)
                endif

                if (lcv(m,n)>0.5) then
                    ! Discretization in h-points
                    k_bfc = FreeFallAcc * (nbfc**2) / (hh(m, n)**(1.0/3.0))
                    s = 0.5d0 * sqrt( (u(m, n) + u(m, n+1))**2 + (v(m, n) + v(m+1, n))**2 )
                    k1 = -dxb(m, n) * dyb(m, n) * k_bfc * s
                    !k1 = k1 * 0.5d0*(v(m, n) + v(m+1, n))

                    ! Discretization in h-points
                    k_bfc = FreeFallAcc * (nbfc**2) / (hh(m-1, n)**(1.0/3.0))
                    s = 0.5d0 * sqrt( (u(m-1, n) + u(m-1, n+1))**2 + (v(m, n) + v(m-1, n))**2 )
                    k2 = -dxb(m-1, n) * dyb(m-1, n) * k_bfc * s
                    !k2 = k2 * 0.5d0*(v(m, n) + v(m-1, n))

                    ! Discretization in v-points
                    RHSy(m, n) = 0.5d0 * (k1 + k2)
                endif
            enddo
        enddo
        !$omp end parallel do

    end subroutine uv_bfc

    subroutine uv_vort(u, v, vort, dxt, dyt, dxb, dyb, luu)
        implicit none

        real(8) u(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &   !Transporting zonal velocity
                v(bnd_x1:bnd_x2,bnd_y1:bnd_y2)      !Transporting meridional velocity

        real(8) vort(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

        real*8 :: dxt(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  dyt(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  dxb(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  dyb(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  luu(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

        integer m, n

        do n = ny_start, ny_end
            do m = nx_start, nx_end
                if (luu(m,n)>0.5) then
                    vort(m,n)= (v(m+1,n)*dyt(m+1,n)-v(m,n)*dyt(m,n))     &
                              -(u(m,n+1)*dxt(m,n+1)-u(m,n)*dxt(m,n))     &
                              -((v(m+1,n)-v(m,n))*dyb(m,n)-(u(m,n+1)-u(m,n))*dxb(m,n))
                endif
            enddo
        enddo

    end subroutine


    subroutine uv_trans(u, v, vort, hq, hu, hv, hh, RHSx, RHSy, dxh, dyh, lcu, lcv, luu)
        implicit none

        integer nlev

        real(8) u(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &    !Transporting zonal velocity
                v(bnd_x1:bnd_x2,bnd_y1:bnd_y2)        !Transporting meridional velocity

        real(8) RHSx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      & !Zonal source function
                RHSy(bnd_x1:bnd_x2,bnd_y1:bnd_y2)         !meridional source function


        real(8) hq(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                hu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                hv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                hh(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

        real(8) vort(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

        real*8 :: dxh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  dyh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  lcu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  lcv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  luu(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

        real(8) fx_p, fx_m, fy_p, fy_m   !fluxes through cell edges

        integer m, n

        do n = ny_start,ny_end
            do m = nx_start,nx_end
                ! Zonal velocity
                if (lcu(m,n)>0.5) then
                    fx_p = (u(m, n  )*dyh(m,n)*hu(m,n) + u(m+1,n  )*dyh(m+1,n)*hu(m+1,n))/2.0d0   &
                          *(u(m, n  ) + u(m+1,n  ))/2.0d0

                    fx_m = (u(m  ,n  )*dyh(m,n)*hu(m,n) + u(m-1,n  )*dyh(m-1,n)*hu(m-1,n))/2.0d0   &
                          *(u(m  ,n  ) + u(m-1,n  ))/2.0d0

                    fy_p = (v(m  ,n  )*dxh(m,n  )*hv(m,n  ) + v(m+1,n  )*dxh(m+1,n  )*hv(m+1,n  ))/2.0d0   &
                          *(u(m  ,n+1) + u(m  ,n  ))/2.0d0*dble(luu(m,n  ))

                    fy_m = (v(m  ,n-1)*dxh(m,n-1)*hv(m,n-1) + v(m+1,n-1)*dxh(m+1,n-1)*hv(m+1,n-1))/2.0d0   &
                          *(u(m  ,n-1) + u(m  ,n  ))/2.0d0*dble(luu(m,n-1))

                    RHSx(m,n) = - (fx_p - fx_m + fy_p - fy_m)           &
                        + ( vort(m,n  )*hh(m,n  )*(v(m+1,n  )+v(m,n  ))              &
                        +   vort(m,n-1)*hh(m,n-1)*(v(m+1,n-1)+v(m,n-1))  )/4.0d0
                end if

                ! Meridional velocity
                if (lcv(m,n)>0.5) then
                    fy_p = (v(m  ,n  )*dxh(m,n)*hv(m,n) + v(m  ,n+1)*dxh(m,n+1)*hv(m,n+1))/2.0d0    &
                          *(v(m  ,n  ) + v(m  ,n+1))/2.0d0

                    fy_m = (v(m  ,n  )*dxh(m,n)*hv(m,n) + v(m  ,n-1)*dxh(m,n-1)*hv(m,n-1))/2.0d0    &
                          *(v(m  ,n  ) + v(m  ,n-1))/2.0d0

                    fx_p = (u(m  ,n  )*dyh(m  ,n)*hu(m  ,n) + u(m  ,n+1)*dyh(m  ,n+1)*hu(m  ,n+1))/2.0d0    &
                          *(v(m+1,n  ) + v(m  ,n  ))/2.0d0

                    fx_m = (u(m-1,n  )*dyh(m-1,n)*hu(m-1,n) + u(m-1,n+1)*dyh(m-1,n+1)*hu(m-1,n+1))/2.0d0    &
                          *(v(m-1,n  ) + v(m  ,n  ))/2.0d0

                    RHSy(m,n) = - (fx_p - fx_m + fy_p - fy_m)          &
                        - ( vort(m  ,n)*hh(m  ,n)*(u(m  ,n+1)+u(m  ,n))               &
                        +   vort(m-1,n)*hh(m-1,n)*(u(m-1,n+1)+u(m-1,n))  )/4.0d0
                end if
            end do
        end do

    endsubroutine uv_trans

    !===========================================================================================
    subroutine uv_diff2(mu, str_t, str_s, hq, hu, hv, hh, RHSx, RHSy, &
                        dx, dy, dxt, dyt, dxh, dyh, dxb, dyb, lcu, lcv)
        implicit none

        real(8) muh_p, muh_m

        real(8) mu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      & !lateral viscosity coefficient
                RHSx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    & !Zonal source function
                RHSy(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    & !meridional source function
                str_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   & !Tension stress
                str_s(bnd_x1:bnd_x2,bnd_y1:bnd_y2)      !Shearing stress

        real(8) hq(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                hu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                hv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                hh(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

        real*8 :: dx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  dy(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  dxt(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  dyt(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  dxh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  dyh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  dxb(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  dyb(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  lcu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  lcv(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

        integer m, n

        do n = ny_start,ny_end
            do m = nx_start,nx_end
                ! Zonal velocity
                if (lcu(m,n)>0.5) then
                    muh_p = (mu(m,n)+mu(m+1,n)+mu(m,n+1)+mu(m+1,n+1))/4.0d0
                    muh_m = (mu(m,n)+mu(m+1,n)+mu(m,n-1)+mu(m+1,n-1))/4.0d0

                    RHSx(m,n) = ( dy(m+1,n)**2*mu(m+1,n)*hq(m+1,n)*str_t(m+1,n)             &
                                 -dy(m  ,n)**2*mu(m  ,n)*hq(m  ,n)*str_t(m  ,n) )/dyh(m,n)  &
                         + (dxb(m,n  )**2*muh_p*hh(m,n  )*str_s(m,n  )                   &
                           -dxb(m,n-1)**2*muh_m*hh(m,n-1)*str_s(m,n-1) )/dxt(m,n)
                end if

                ! Meridional velocity
                if (lcv(m,n)>0.5) then
                    muh_p = (mu(m,n)+mu(m+1,n)+mu(m,n+1)+mu(m+1,n+1))/4.0d0
                    muh_m = (mu(m,n)+mu(m-1,n)+mu(m,n+1)+mu(m-1,n+1))/4.0d0

                    RHSy(m,n) = -( dx(m,n+1)**2*mu(m,n+1)*hq(m,n+1)*str_t(m,n+1)              &
                                  -dx(m,n  )**2*mu(m,n  )*hq(m,n  )*str_t(m,n  ) ) /dxh(m,n)  &
                           + (dyb(m  ,n)**2*muh_p*hh(m  ,n)*str_s(m  ,n)                    &
                             -dyb(m-1,n)**2*muh_m*hh(m-1,n)*str_s(m-1,n) ) /dyt(m,n)
                end if
            end do
        end do

    endsubroutine uv_diff2

end module vel_ssh
