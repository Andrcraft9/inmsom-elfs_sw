!====================================================================================
subroutine stress_components(u, v, str_t, str_s, nlev, bnd_step)
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    implicit none

    integer nlev
    real(8) u(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev),    &
            v(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev),    &
            str_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev),    &
            str_s(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev)

    integer m,n,k
    integer bnd_step

    !$omp parallel do private(m,n,k)
    do n = max(bnd_y1, ny_start - bnd_step), min(bnd_y2, ny_end + bnd_step)
        do m = max(bnd_x1, nx_start - bnd_step), min(bnd_x2, nx_end + bnd_step)
            if(lu(m,n)>0.5) then
                do k=1,nlev
                    str_t(m,n,k)=dy(m,n)/dx(m,n)*(u(m,n,k)/dyh(m,n)-u(m-1,n,k)/dyh(m-1,n))     &
                        -dx(m,n)/dy(m,n)*(v(m,n,k)/dxh(m,n)-v(m,n-1,k)/dxh(m,n-1))
                enddo
            endif
        enddo
    enddo
    !$omp end parallel do

    !$omp parallel do private(m,n,k)
    do n = max(bnd_y1, ny_start - bnd_step), min(bnd_y2, ny_end + bnd_step)
        do m = max(bnd_x1, nx_start - bnd_step), min(bnd_x2, nx_end + bnd_step)
            if(luu(m,n)>0.5) then
                do k=1,nlev
                    str_s(m,n,k)=dxb(m,n)/dyb(m,n)*(u(m,n+1,k)/dxt(m,n+1)-u(m,n,k)/dxt(m,n))     &
                        +dyb(m,n)/dxb(m,n)*(v(m+1,n,k)/dyt(m+1,n)-v(m,n,k)/dyt(m,n))
                enddo
            endif
        enddo
    enddo
    !$omp end parallel do

endsubroutine stress_components
