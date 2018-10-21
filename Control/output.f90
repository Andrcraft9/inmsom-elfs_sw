module output_routes
implicit none

contains

subroutine test_sync
    use mpi_parallel_tools
    implicit none

    integer :: k, m, n, ierr
    type(block2D), dimension(:), pointer :: testv
    character(len=1024) :: filename

    call allocate_block2D(testv, 0.0d0)
    do k = 1, bcount

        do m = bbnd_x1(k), bbnd_x2(k)
            do n = bbnd_y1(k), bbnd_y2(k)
                testv(k)%vals(m, n) = 9 !bindx(k, 1)
            enddo
        enddo

        do m = bnx_start(k), bnx_end(k)
            do n = bny_start(k), bny_end(k)
                testv(k)%vals(m, n) = lbasins(m, n)
            enddo
        enddo
    enddo

    call syncborder_block2D(testv)

    do k = 1, bcount
        write(filename, "(A5, I1, A1, I2, A9)") 'data/', rank, '_', k, '_test.txt'
        open(12, file=trim(filename), status='replace')
        !do n = bny_end(k), bny_start(k), -1
        !    do m = bnx_start(k), bnx_end(k)
        do n = bbnd_y2(k), bbnd_y1(k), -1
            do m = bbnd_x1(k), bbnd_x2(k)
                write(12, '(I1)', advance='no') int(testv(k)%vals(m, n))
            enddo
            write(12, *) ''
        enddo
        close(12)
    enddo

    call deallocate_blocks2D(testv)
end subroutine

subroutine print_grid
    use mpi_parallel_tools
    use basin_grid
    implicit none

    integer :: k, m, n
    character(len=1024) :: filename

    do k = 1, bcount
        write(filename, "(A5, I1, A1, I2, A8)") 'data/', rank, '_', k, '_lu.txt'
        open(12, file=trim(filename), status='replace')
        !do n = bny_end(k), bny_start(k), -1
        !    do m = bnx_start(k), bnx_end(k)
        do n = bbnd_y2(k), bbnd_y1(k), -1
            do m = bbnd_x1(k), bbnd_x2(k)
                write(12, '(I1)', advance='no') int(lu(k)%vals(m, n))
            enddo
            write(12, *) ''
        enddo
        close(12)
    enddo

end subroutine

subroutine parallel_check_point()
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables

    implicit none
    include 'locout.fi'

    integer :: m, n, k
    integer :: rb(2)
    real*8 :: lon, lat

    if (points_output > 0) then
        do k = 1, nloc_points
            lon = lon_loc_points(k)
            lat = lat_loc_points(k)

            m = floor((lon - rlon) / dxst) + mmm
            n = floor((lat - rlat) / dyst) + nnn

            !r = get_rank_by_point(m, n)
            call get_block_and_rank_by_point(m, n, rb)

            if (rank .eq. rb(1)) then
                print *, 'rank, k, lon, lat, geo_lon_t, geo_lat_t, hhq_rest'
                print *,  rank, rb(2), lon, lat, geo_lon_t(rb(2))%vals(m, n), geo_lat_t(rb(2))%vals(m, n), hhq_rest(rb(2))%vals(m, n)
            endif
        enddo
    endif

end subroutine


subroutine parallel_point_output(path2data, point_time)
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables
    use iodata_routes

    implicit none
    include 'locout.fi'

    character fname*256
    character*(*) path2data
    real*8 :: point_time
    integer :: k, m, n, ierr
    integer :: rb(2)
    real*8 :: lon, lat

    if (points_output > 0) then
        do k = 1, nloc_points
            lon = lon_loc_points(k)
            lat = lat_loc_points(k)

            m = floor((lon - rlon) / dxst) + mmm
            n = floor((lat - rlat) / dyst) + nnn
            !r = get_rank_by_point(m, n)
            call get_block_and_rank_by_point(m, n, rb)

            if (rank .eq. rb(1)) then
                call fulfname(fname, path2data, name_points(k), ierr)
                open(40, file=fname, status='unknown', position='append')
                write(40, *) point_time, ssh(rb(2))%vals(m, n)
                close(40)
            endif
        enddo
    endif

    return
end subroutine parallel_point_output

subroutine parallel_energy_output(path2data, point_time)
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables
    use iodata_routes

    implicit none
    include 'locout.fi'

    character fname*256
    character*(*) path2data
    real*8 :: point_time
    integer :: k, m, n, ierr
    real*8 :: kinetic_e, potential_e
    real*8 :: buf_k, buf_p

    if (energy_output > 0) then
        ! TotalEnergy = 0.5 rho dx dy Sum[h_u u**2 + h_v v**2 + g ssh**2]
        kinetic_e = 0.0d0; potential_e = 0.0d0
        do k = 1, bcount
            call set_block_boundary(k)
            do n=ny_start,ny_end
                do m=nx_start,nx_end
                    if (lcu(k)%vals(m,n)>0.5) then
                        kinetic_e = kinetic_e + 0.5d0*RefDen*dxt(k)%vals(m,n)*dyh(k)%vals(m,n)*hhu(k)%vals(m,n)*(ubrtr(k)%vals(m,n)**2)
                    endif
                    if (lcv(k)%vals(m,n)>0.5) then
                        kinetic_e = kinetic_e + 0.5d0*RefDen*dxh(k)%vals(m,n)*dyt(k)%vals(m,n)*hhv(k)%vals(m,n)*(vbrtr(k)%vals(m,n)**2)
                    endif

                    if (lu(k)%vals(m,n)>0.5) then
                        potential_e = potential_e + 0.5d0*RefDen*dx(k)%vals(m,n)*dy(k)%vals(m,n)*FreeFallAcc*(ssh(k)%vals(m,n)**2)
                    endif
                enddo
            enddo
        enddo
        buf_k = kinetic_e; buf_p = potential_e
        call mpi_allreduce(buf_k, kinetic_e, 1, mpi_real8, mpi_sum, cart_comm, ierr)
        call mpi_allreduce(buf_p, potential_e, 1, mpi_real8, mpi_sum, cart_comm, ierr)

        if (rank .eq. 0) then
            call fulfname(fname, path2data, 'kinetic_potential_energy', ierr)
            open(40, file=fname, status='unknown', position='append')
            write(40, *) point_time, kinetic_e, potential_e
            close(40)
        endif
    endif
    
    return
end subroutine

endmodule output_routes
