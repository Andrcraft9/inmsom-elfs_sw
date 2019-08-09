module mixing_routes
    use mpi_parallel_tools
    use basin_grid
    use ocalg_routes
    implicit none

contains

    subroutine stress_components(u, v, str_t, str_s)
        implicit none
        real(8) u(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                v(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                str_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &
                str_s(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

        integer m, n

        do n=ny_start, ny_end
            do m=nx_start, nx_end
                if(lu(m,n)>0.5) then
                        str_t(m,n) = dy(m,n)/dx(m,n)*(u(m,n)/dyh(m,n) - u(m-1,n)/dyh(m-1,n))     &
                                    -dx(m,n)/dy(m,n)*(v(m,n)/dxh(m,n) - v(m,n-1)/dxh(m,n-1))
                endif
            enddo
        enddo

        do n=ny_start, ny_end
            do m=nx_start, nx_end
                if(luu(m,n)>0.5) then
                        str_s(m,n) = dxb(m,n)/dyb(m,n)*(u(m,n+1)/dxt(m,n+1) - u(m,n)/dxt(m,n))     &
                                    +dyb(m,n)/dxb(m,n)*(v(m+1,n)/dyt(m+1,n) - v(m,n)/dyt(m,n))
                endif
            enddo
        enddo

    endsubroutine stress_components

endmodule mixing_routes
