module depth
    implicit none

contains

    subroutine hh_init(hq, hqp, hqn,    &
                       hu, hup, hun,    &
                       hv, hvp, hvn,    &
                       hh, hhp, hhn,    &
                       sh, shp, h_r)
        use main_basin_pars
        use mpi_parallel_tools
        use basin_grid
        implicit none

        type(block2D), dimension(:), pointer :: hq, hqp, hqn,    &
                                                hu, hup, hun,    &
                                                hv, hvp, hvn,    &
                                                hh, hhp, hhn,    &
                                                sh, shp, h_r

        real(8) slu
        integer k, m, n

        do k = 1, bcount
            hq(k)%vals = h_r(k)%vals + sh(k)%vals * dfloat(full_free_surface)
            hqp(k)%vals = h_r(k)%vals + shp(k)%vals * dfloat(full_free_surface)
            hqn(k)%vals = h_r(k)%vals
        enddo

        !$omp parallel do
        do k = 1, bcount
            call set_block_boundary(k)
            !do n=ny_start-2, ny_end+1
            !do m=nx_start-2, nx_end+1
            do n=ny_start-1,ny_end
                do m=nx_start-1,nx_end
                    if (llu(k)%vals(m,n) > 0.5) then
                        ! interpolating hhq given on T-grid(lu) to hhu given on u-grid(lcu).
                        slu = dble(lu(k)%vals(m,n) + lu(k)%vals(m+1,n))

                        hu(k)%vals(m,n) = ( hq(k)%vals(m  ,n)*dx(k)%vals(m  ,n)*dy(k)%vals(m  ,n)*dble(lu(k)%vals(m  ,n))   &
                                          + hq(k)%vals(m+1,n)*dx(k)%vals(m+1,n)*dy(k)%vals(m+1,n)*dble(lu(k)%vals(m+1,n)) ) &
                                          /slu/dxt(k)%vals(m,n)/dyh(k)%vals(m,n)

                        hup(k)%vals(m,n) = ( hqp(k)%vals(m  ,n)*dx(k)%vals(m  ,n)*dy(k)%vals(m  ,n)*dble(lu(k)%vals(m  ,n))   &
                                           + hqp(k)%vals(m+1,n)*dx(k)%vals(m+1,n)*dy(k)%vals(m+1,n)*dble(lu(k)%vals(m+1,n)) ) &
                                           /slu/dxt(k)%vals(m,n)/dyh(k)%vals(m,n)

                        hun(k)%vals(m,n) = ( hqn(k)%vals(m  ,n)*dx(k)%vals(m  ,n)*dy(k)%vals(m  ,n)*dble(lu(k)%vals(m  ,n))   &
                                           + hqn(k)%vals(m+1,n)*dx(k)%vals(m+1,n)*dy(k)%vals(m+1,n)*dble(lu(k)%vals(m+1,n)) ) &
                                           /slu/dxt(k)%vals(m,n)/dyh(k)%vals(m,n)
                    endif

                    if (llv(k)%vals(m,n) > 0.5) then
                        ! interpolating hhq given on T-grid(lu) to hhv given on v-grid(lcv).
                        slu = dble(lu(k)%vals(m,n) + lu(k)%vals(m,n+1))

                        hv(k)%vals(m,n) = ( hq(k)%vals(m,n  )*dx(k)%vals(m,n  )*dy(k)%vals(m,n  )*dble(lu(k)%vals(m,n  ))       &
                                          + hq(k)%vals(m,n+1)*dx(k)%vals(m,n+1)*dy(k)%vals(m,n+1)*dble(lu(k)%vals(m,n+1)) ) &
                                          /slu/dxh(k)%vals(m,n)/dyt(k)%vals(m,n)

                        hvp(k)%vals(m,n) = ( hqp(k)%vals(m,n  )*dx(k)%vals(m,n  )*dy(k)%vals(m,n  )*dble(lu(k)%vals(m,n  ))       &
                                           + hqp(k)%vals(m,n+1)*dx(k)%vals(m,n+1)*dy(k)%vals(m,n+1)*dble(lu(k)%vals(m,n+1)) ) &
                                           /slu/dxh(k)%vals(m,n)/dyt(k)%vals(m,n)

                        hvn(k)%vals(m,n) = ( hqn(k)%vals(m,n  )*dx(k)%vals(m,n  )*dy(k)%vals(m,n  )*dble(lu(k)%vals(m,n  ))       &
                                           + hqn(k)%vals(m,n+1)*dx(k)%vals(m,n+1)*dy(k)%vals(m,n+1)*dble(lu(k)%vals(m,n+1)) ) &
                                           /slu/dxh(k)%vals(m,n)/dyt(k)%vals(m,n)
                    endif

                    if (luh(k)%vals(m,n) > 0.5) then
                        ! interpolating hhq given on T-grid(lu) to hhh given on h-grid(luu).
                        slu = dble(lu(k)%vals(m,n)+lu(k)%vals(m+1,n)+lu(k)%vals(m,n+1)+lu(k)%vals(m+1,n+1))

                        hh(k)%vals(m,n)=( hq(k)%vals(m  ,n  )*dx(k)%vals(m  ,n  )*dy(k)%vals(m  ,n  )*dble(lu(k)%vals(m  ,n  ))       &
                                        + hq(k)%vals(m+1,n  )*dx(k)%vals(m+1,n  )*dy(k)%vals(m+1,n  )*dble(lu(k)%vals(m+1,n  ))       &
                                        + hq(k)%vals(m  ,n+1)*dx(k)%vals(m  ,n+1)*dy(k)%vals(m  ,n+1)*dble(lu(k)%vals(m  ,n+1))       &
                                        + hq(k)%vals(m+1,n+1)*dx(k)%vals(m+1,n+1)*dy(k)%vals(m+1,n+1)*dble(lu(k)%vals(m+1,n+1)) ) &
                                        /slu/dxb(k)%vals(m,n)/dyb(k)%vals(m,n)

                        hhp(k)%vals(m,n) = ( hqp(k)%vals(m  ,n  )*dx(k)%vals(m  ,n  )*dy(k)%vals(m  ,n  )*dble(lu(k)%vals(m  ,n  ))       &
                                           + hqp(k)%vals(m+1,n  )*dx(k)%vals(m+1,n  )*dy(k)%vals(m+1,n  )*dble(lu(k)%vals(m+1,n  ))       &
                                           + hqp(k)%vals(m  ,n+1)*dx(k)%vals(m  ,n+1)*dy(k)%vals(m  ,n+1)*dble(lu(k)%vals(m  ,n+1))       &
                                           + hqp(k)%vals(m+1,n+1)*dx(k)%vals(m+1,n+1)*dy(k)%vals(m+1,n+1)*dble(lu(k)%vals(m+1,n+1)) ) &
                                           /slu/dxb(k)%vals(m,n)/dyb(k)%vals(m,n)

                        hhn(k)%vals(m,n) = ( hqn(k)%vals(m  ,n  )*dx(k)%vals(m  ,n  )*dy(k)%vals(m  ,n  )*dble(lu(k)%vals(m  ,n  ))       &
                                           + hqn(k)%vals(m+1,n  )*dx(k)%vals(m+1,n  )*dy(k)%vals(m+1,n  )*dble(lu(k)%vals(m+1,n  ))       &
                                           + hqn(k)%vals(m  ,n+1)*dx(k)%vals(m  ,n+1)*dy(k)%vals(m  ,n+1)*dble(lu(k)%vals(m  ,n+1))       &
                                           + hqn(k)%vals(m+1,n+1)*dx(k)%vals(m+1,n+1)*dy(k)%vals(m+1,n+1)*dble(lu(k)%vals(m+1,n+1)) ) &
                                           /slu/dxb(k)%vals(m,n)/dyb(k)%vals(m,n)
                    endif

                end do
            end do
        enddo
        !$omp end parallel do

        call syncborder_block2D(hu)
        call syncborder_block2D(hup)
        call syncborder_block2D(hun)
        call syncborder_block2D(hv)
        call syncborder_block2D(hvp)
        call syncborder_block2D(hvn)
        call syncborder_block2D(hh)
        call syncborder_block2D(hhp)
        call syncborder_block2D(hhn)

    endsubroutine hh_init

    !============================================================
    subroutine hh_update(hqn, hun, hvn, hhn, sh, h_r)
        use main_basin_pars
        use mpi_parallel_tools
        use basin_grid
        implicit none

        type(block2D), dimension(:), pointer :: hqn, hun,    &
                                                hvn, hhn,    &
                                                sh,  h_r

        integer k, m, n
        real(8) slu

        do k = 1, bcount
            hqn(k)%vals = h_r(k)%vals + sh(k)%vals
        enddo

        !$omp parallel do
        do k = 1, bcount
            call set_block_boundary(k)

            !do n=ny_start-2, ny_end+1
            !do m=nx_start-2, nx_end+1
            do n=ny_start-1,ny_end
                do m=nx_start-1,nx_end
                    if (llu(k)%vals(m, n) > 0.5) then
                        ! interpolating hhq given on T-grid(lu) to hhu given on u-grid(lcu).
                        slu = dble(lu(k)%vals(m,n)+lu(k)%vals(m+1,n))
                        hun(k)%vals(m,n) = ( hqn(k)%vals(m  ,n)*dx(k)%vals(m  ,n)*dy(k)%vals(m  ,n)*dble(lu(k)%vals(m  ,n))   &
                                           + hqn(k)%vals(m+1,n)*dx(k)%vals(m+1,n)*dy(k)%vals(m+1,n)*dble(lu(k)%vals(m+1,n)) ) &
                                           /slu/dxt(k)%vals(m,n)/dyh(k)%vals(m,n)
                    endif
                    if (llv(k)%vals(m,n) > 0.5) then
                        ! interpolating hhq given on T-grid(lu) to hhv given on v-grid(lcv).
                        slu = dble(lu(k)%vals(m,n)+lu(k)%vals(m,n+1))
                        hvn(k)%vals(m,n) = ( hqn(k)%vals(m,n  )*dx(k)%vals(m,n  )*dy(k)%vals(m,n  )*dble(lu(k)%vals(m,n  ))       &
                                           + hqn(k)%vals(m,n+1)*dx(k)%vals(m,n+1)*dy(k)%vals(m,n+1)*dble(lu(k)%vals(m,n+1)) ) &
                                           /slu/dxh(k)%vals(m,n)/dyt(k)%vals(m,n)
                    endif
                    if (luh(k)%vals(m,n) > 0.5) then
                        ! interpolating hhq given on T-grid(lu) to hhh given on h-grid(luu).
                        slu = dble(lu(k)%vals(m,n)+lu(k)%vals(m+1,n)+lu(k)%vals(m,n+1)+lu(k)%vals(m+1,n+1))
                        hhn(k)%vals(m,n) = ( hqn(k)%vals(m  ,n  )*dx(k)%vals(m  ,n  )*dy(k)%vals(m  ,n  )*dble(lu(k)%vals(m  ,n  ))       &
                                           + hqn(k)%vals(m+1,n  )*dx(k)%vals(m+1,n  )*dy(k)%vals(m+1,n  )*dble(lu(k)%vals(m+1,n  ))       &
                                           + hqn(k)%vals(m  ,n+1)*dx(k)%vals(m  ,n+1)*dy(k)%vals(m  ,n+1)*dble(lu(k)%vals(m  ,n+1))       &
                                           + hqn(k)%vals(m+1,n+1)*dx(k)%vals(m+1,n+1)*dy(k)%vals(m+1,n+1)*dble(lu(k)%vals(m+1,n+1)) ) &
                                           /slu/dxb(k)%vals(m,n)/dyb(k)%vals(m,n)
                    endif

                end do
            end do
        enddo
        !$omp end parallel do

        call syncborder_block2D(hun)
        call syncborder_block2D(hvn)
        call syncborder_block2D(hhn)

    endsubroutine hh_update

    subroutine hh_shift(hq, hqp, hqn,   &
                        hu, hup, hun,   &
                        hv, hvp, hvn,   &
                        hh, hhp, hhn)
        use main_basin_pars
        use mpi_parallel_tools
        use basin_grid
        implicit none

        type(block2D), dimension(:), pointer :: hq, hqp, hqn,    &
                                                hu, hup, hun,    &
                                                hv, hvp, hvn,    &
                                                hh, hhp, hhn

        integer k, m, n

        !$omp parallel do
        do k = 1, bcount
            call set_block_boundary(k)
            do n=ny_start-1,ny_end+1
                do m=nx_start-1,nx_end+1
                    if (llu(k)%vals(m,n) > 0.5) then
                        hup(k)%vals(m,n) = hu(k)%vals(m,n)  &
                            + time_smooth*(hun(k)%vals(m,n) - 2.0d0*hu(k)%vals(m,n) + hup(k)%vals(m,n))/2.0d0

                        hu(k)%vals(m,n)= hun(k)%vals(m,n)
                    endif
                    if (llv(k)%vals(m,n) > 0.5) then
                        hvp(k)%vals(m,n) = hv(k)%vals(m,n)  &
                            + time_smooth*(hvn(k)%vals(m,n) - 2.0d0*hv(k)%vals(m,n) + hvp(k)%vals(m,n))/2.0d0

                        hv(k)%vals(m,n)= hvn(k)%vals(m,n)
                    endif
                    if (lu(k)%vals(m,n) > 0.5) then
                        hqp(k)%vals(m,n) = hq(k)%vals(m,n)  &
                            + time_smooth*(hqn(k)%vals(m,n) - 2.0d0*hq(k)%vals(m,n) + hqp(k)%vals(m,n))/2.0d0

                        hq(k)%vals(m,n)= hqn(k)%vals(m,n)
                    endif
                    if (luh(k)%vals(m,n) > 0.5) then
                        hhp(k)%vals(m,n) = hh(k)%vals(m,n)  &
                            + time_smooth*(hhn(k)%vals(m,n) - 2.0d0*hh(k)%vals(m,n) + hhp(k)%vals(m,n))/2.0d0

                        hh(k)%vals(m,n)= hhn(k)%vals(m,n)
                    endif
                end do
            end do
        enddo
        !$omp end parallel do

    endsubroutine hh_shift

end module depth
