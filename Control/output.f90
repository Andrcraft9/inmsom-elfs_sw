module output_routes
    implicit none

contains

    subroutine init_model_output()
        use mpi_parallel_tools
        use basin_grid
        implicit none
        include 'locout.fi'

        integer :: ierr
        integer :: r, locr

        ! First point of entire basin area for output
        call get_tpoint(xtm1, ytn1, m1loc, n1loc)
        call get_uvpoint(xum1, yvn1, m1loc-1, n1loc-1)

    end subroutine

    subroutine parallel_check_point()
        use mpi_parallel_tools
        use basin_grid
        use ocean_variables

        implicit none
        include 'locout.fi'

        integer :: m, n, k
        real*8 :: lon, lat
        integer :: rb(2)

        if (points_output > 1) then
            do k = 1, nloc_points
                lon = lon_loc_points(k)
                lat = lat_loc_points(k)
                m = floor((lon - rlon) / dxst) + mmm
                n = floor((lat - rlat) / dyst) + nnn
                call get_block_and_rank_by_point(m, n, rb)
                if (rank .eq. rb(1)) then
                    call set_block(rb(2))
                    call set_block_lu(rb(2))
                    call set_block_h(rb(2))
                    call set_block_geo(rb(2))
                    print *, 'rank ', ' lon ', ' lat ', ' geo_lon_t ', ' geo_lat_t', 'hhq_rest'
                    print *,  rank,   lon,   lat,   geo_lon_t(m, n), geo_lat_t(m, n), hhq_rest(m, n)
                endif
            enddo
        endif

    end subroutine

    subroutine parallel_point_output(path2data, point_time)
        use mpi_parallel_tools
        use basin_grid
        use ocean_variables
        use iodata_routes
        use rw_ctl_routes

        implicit none
        include 'locout.fi'

        character fname*256
        character*(*) path2data
        real*8 :: point_time
        integer :: m, n, ierr, k
        real*8 :: lon, lat
        real*8 :: uuu, vvv
        integer :: rb(2)

        if (points_output > 0) then
            do k = 1, nloc_points
                lon = lon_loc_points(k)
                lat = lat_loc_points(k)
                m = floor((lon - rlon) / dxst) + mmm
                n = floor((lat - rlat) / dyst) + nnn
                call get_block_and_rank_by_point(m, n, rb)
                if (rank .eq. rb(1)) then
                    call set_block(rb(2))
                    call set_block_lu(rb(2))
                    call set_block_h(rb(2))
                    call set_block_dxdy(rb(2))

                    uuu = ( ubrtr(rb(2))%vals(m  ,n)*dyh(m  ,n)*hhu(m  ,n)    &
                           +ubrtr(rb(2))%vals(m-1,n)*dyh(m-1,n)*hhu(m-1,n) )/2.0/hhq(m,n)/dy(m,n)
                    vvv = ( vbrtr(rb(2))%vals(m,n  )*dxh(m,n  )*hhv(m,n  )    &
                           +vbrtr(rb(2))%vals(m,n-1)*dxh(m,n-1)*hhv(m,n-1) )/2.0/hhq(m,n)/dx(m,n)

                    call fulfname(fname, path2data, name_points(k), ierr)
                    open(40, file=fname, status='unknown', position='append')
                    write(40, *) point_time, ssh(rb(2))%vals(m, n), uuu, vvv
                    close(40)
                endif
            enddo
        endif

        return
    end subroutine parallel_point_output

    subroutine parallel_energy_output(path2data, point_time)
        use mpi_parallel_tools
        use basin_grid
        use ocean_variables
        use iodata_routes
        use rw_ctl_routes

        implicit none
        include 'locout.fi'

        character fname*256
        character*(*) path2data
        real*8 :: point_time
        integer :: m, n, k, ierr
        real*8 :: kinetic_e, potential_e
        real*8 :: buf_k, buf_p

        if (energy_output > 0) then
            ! TotalEnergy = 0.5 rho dx dy Sum[h_u u**2 + h_v v**2 + g ssh**2]
            kinetic_e = 0.0d0; potential_e = 0.0d0
            do k = 1, bcount
                call set_block(k)
                call set_block_lu_output(k)
                call set_block_h(k)
                call set_block_dxdy(k)
                do n = ny_start, ny_end
                    do m = nx_start, nx_end
                        if (lcu_output(m,n)>0.5) then
                            kinetic_e = kinetic_e + 0.5d0*RefDen*dxt(m,n)*dyh(m,n)*hhu(m,n)*(ubrtr(k)%vals(m,n)**2)
                        endif
                        if (lcv_output(m,n)>0.5) then
                            kinetic_e = kinetic_e + 0.5d0*RefDen*dxh(m,n)*dyt(m,n)*hhv(m,n)*(vbrtr(k)%vals(m,n)**2)
                        endif

                        if (lu_output(m,n)>0.5) then
                            potential_e = potential_e + 0.5d0*RefDen*dx(m,n)*dy(m,n)*FreeFallAcc*(ssh(k)%vals(m,n)**2)
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

    subroutine parallel_local_output(path2data,  &
                                     nrec,       &
                                     year,       &
                                    month,       &
                                      day,       &
                                     hour,       &
                                   minute,       &
                                   tstep,        &
                                   calendar  )
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables
    use rec_length
    use iodata_routes
    use rw_ctl_routes

    implicit none
    include 'locout.fi'

    integer nrec, year, month, day, hour, minute, calendar, ierr
    character fname*256
    character*(*) path2data
    real(4) tstep
    integer m,n,k
    real(8) z0(1), z1(1)

    z0 = 0.0d0; z1 = 1.0d0

    if (rank .eq. 0) write(*,*) 'Writing local output, record number ', nrec

    if(nrec==1) then
     if (rank == 0) then
         print *, "first x-value: ",  xtm1, "first y-value: ",  ytn1
     endif
     !writing HHQ
     ierr=0
     do k = 1, bcount
       call set_block(k)
       array4_2d(k)%vals = sngl(block_hhq_rest(k)%vals)
     enddo
     call pwdstd2D(path2data,'LOCAL/hhq.dat',nrec,array4_2d,block_lu_output,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,cart_comm,ierr)
     call fulfname(fname,path2data,'LOCAL/hhq.dat',ierr)
     if (rank .eq. 0) then
         call ctl_file_write(fname,    &     !file name
                         undef,    &     !value for undefined points
                         nx_loc,   &     !x-dimension
                         ny_loc,   &     !y-dimension
                              1,   &     !z-dimension
                           nrec,   &     !t-dimension
                       xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                           xtm1,   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                           ytn1,   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          0,       &     !z-grid type (0 - linear, 1 - levels)
                          z0,      &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                       calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                           year,   &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                           hour,   &     !hour   of the first field
                         minute,   &     !minute of the first field
                          tstep,   &     !time step (in seconds)
                       'HHQ, m',   &     !title of dataset
                          'hhq'   )      !variable name
     endif
    endif

    if (ssh_output>0) then
    ! writing SSH
     ierr=0
     do k = 1, bcount
       call set_block(k)
       array4_2d(k)%vals = sngl(ssh(k)%vals)
     enddo
    ! call wdstd(path2data,'LOCAL/ssh.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
     call pwdstd2D(path2data,'LOCAL/ssh.dat',nrec,array4_2d,block_lu_output,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,cart_comm,ierr)
     call fulfname(fname,path2data,'LOCAL/ssh.dat',ierr)
     if (rank .eq. 0) then
         call ctl_file_write(fname,    &     !file name
                         undef,    &     !value for undefined points
                         nx_loc,   &     !x-dimension
                         ny_loc,   &     !y-dimension
                              1,   &     !z-dimension
                           nrec,   &     !t-dimension
                       xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                           xtm1,   &     !first x-value (if linear) or x-array (if levels)
                          dxst,    &     !x-step (if linear)
                      ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                           ytn1,   &     !first y-value (if linear) or x-array (if levels)
                          dyst,    &     !y-step (if linear)
                          0,       &     !z-grid type (0 - linear, 1 - levels)
                          z0,      &     !first z-value (if linear) or x-array (if levels)
                          1.0d0,   &     !z-step (if linear)
                       calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                           year,   &     !year   of the first field
                          month,   &     !month  of the first field
                            day,   &     !day    of the first field
                           hour,   &     !hour   of the first field
                         minute,   &     !minute of the first field
                          tstep,   &     !time step (in seconds)
                       'SSH, m',   &     !title of dataset
                          'ssh'   )      !variable name
     endif
    endif

    if (uv_output>0) then
        if (grid_shift == 0) then
            ! writing on model grid
            ierr = 0
            do k = 1, bcount
              call set_block(k)
              array4_2d(k)%vals = sngl(ubrtr(k)%vals)
            enddo
            call pwdstd2D(path2data,'LOCAL/ubrtr.dat',nrec,array4_2d,block_llu_output,nx,ny,1,m1loc-1,m2loc,n1loc,n2loc,1,1,cart_comm,ierr)
            call fulfname(fname,path2data,'LOCAL/ubrtr.dat',ierr)
            if (rank .eq. 0) then
                call ctl_file_write(fname,    &     !file name
                                    undef,    &     !value for undefined points
                                  nx_loc+1,   &     !x-dimension
                                    ny_loc,   &     !y-dimension
                                         1,   &     !z-dimension
                                      nrec,   &     !t-dimension
                                  xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                                      xum1,   &     !first x-value (if linear) or x-array (if levels)
                                     dxst,    &     !x-step (if linear)
                                 ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                                      ytn1,   &     !first y-value (if linear) or x-array (if levels)
                                     dyst,    &     !y-step (if linear)
                                     0,       &     !z-grid type (0 - linear, 1 - levels)
                                     z0,      &     !first z-value (if linear) or x-array (if levels)
                                     1.0d0,   &     !z-step (if linear)
                                  calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                                      year,   &     !year   of the first field
                                     month,   &     !month  of the first field
                                       day,   &     !day    of the first field
                                      hour,   &     !hour   of the first field
                                    minute,   &     !minute of the first field
                                     tstep,   &     !time step (in seconds
                      'zonal velocity, m/s',  &     !title of dataset
                                       'u'   )      !variable name
            endif
        else
            ! writing on T-grid
            do k = 1, bcount
                call set_block(k)
                call set_block_lu(k)
                call set_block_dxdy(k)
                call set_block_h(k)
                do n=ny_start, ny_end
                  do m=nx_start, nx_end
                    if (lu(m,n)>0.5) then
                        array4_2d(k)%vals(m, n) = sngl( (ubrtr(k)%vals(m  ,n)*dyh(m  ,n)*hhu(m  ,n)   &
                                                        +ubrtr(k)%vals(m-1,n)*dyh(m-1,n)*hhu(m-1,n) )/2.0/hhq(m,n)/dy(m,n) )
                    endif
                  enddo
                enddo
            enddo

            call pwdstd2D(path2data,'LOCAL/ubrtr.dat',nrec,array4_2d,block_lu_output,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,cart_comm,ierr)
            call fulfname(fname,path2data,'LOCAL/ubrtr.dat',ierr)
            if (rank == 0) then
              call ctl_file_write(fname,    &     !file name
                                  undef,    &     !value for undefined points
                                  nx_loc,   &     !x-dimension
                                  ny_loc,   &     !y-dimension
                                       1,   &     !z-dimension
                                    nrec,   &     !t-dimension
                                xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                                   xtm1,    &     !first x-value (if linear) or x-array (if levels)
                                  dxst,     &     !x-step (if linear)
                              ygr_type,     &     !y-grid type (0 - linear, 1 - levels)
                                   ytn1,    &     !first y-value (if linear) or x-array (if levels)
                                  dyst,     &     !y-step (if linear)
                                      0,    &     !z-grid type (0 - linear, 1 - levels)
                                      z0,   &     !first z-value (if linear) or x-array (if levels)
                                  1.0d0,    &     !z-step (if linear)
                                calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                                    year,   &     !year   of the first field
                                  month,    &     !month  of the first field
                                    day,    &     !day    of the first field
                                    hour,   &     !hour   of the first field
                                  minute,   &     !minute of the first field
                                  tstep,    &     !time step (in seconds
                    'zonal velocity, m/s',  &     !title of dataset
                                    'u'    )      !variable name
            endif
        endif

        if (grid_shift == 0) then
            ierr=0
            do k = 1, bcount
              call set_block(k)
              array4_2d(k)%vals = sngl(vbrtr(k)%vals)
            enddo
            call pwdstd2D(path2data,'LOCAL/vbrtr.dat',nrec,array4_2d,block_llv_output,nx,ny,1,m1loc,m2loc,n1loc-1,n2loc,1,1,cart_comm,ierr)
            call fulfname(fname,path2data,'LOCAL/vbrtr.dat',ierr)
            if (rank .eq. 0) then
                call ctl_file_write(fname,    &     !file name
                                    undef,    &     !value for undefined points
                                    nx_loc,   &     !x-dimension
                                  ny_loc+1,   &     !y-dimension
                                         1,   &     !z-dimension
                                      nrec,   &     !t-dimension
                                  xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                                      xtm1,   &     !first x-value (if linear) or x-array (if levels)
                                     dxst,    &     !x-step (if linear)
                                 ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                                      yvn1,   &     !first y-value (if linear) or x-array (if levels)
                                     dyst,    &     !y-step (if linear)
                                     0,       &     !z-grid type (0 - linear, 1 - levels)
                                     z0,      &     !first z-value (if linear) or x-array (if levels)
                                     1.0d0,   &     !z-step (if linear)
                                  calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                                      year,   &     !year   of the first field
                                     month,   &     !month  of the first field
                                       day,   &     !day    of the first field
                                      hour,   &     !hour   of the first field
                                    minute,   &     !minute of the first field
                                     tstep,   &     !time step (in seconds
                 'meridional velocity, m/s',  &     !title of dataset
                                       'v'   )      !variable name
            endif
        else
            !writing on T-grid
            do k = 1, bcount
                call set_block(k)
                call set_block_lu(k)
                call set_block_dxdy(k)
                call set_block_h(k)
                do n=ny_start, ny_end
                  do m=nx_start, nx_end
                    if (lu(m,n)>0.5) then
                        array4_2d(k)%vals(m,n) = sngl( (vbrtr(k)%vals(m,n  )*dxh(m,n  )*hhv(m,n  )    &
                                                       +vbrtr(k)%vals(m,n-1)*dxh(m,n-1)*hhv(m,n-1))/2.0/hhq(m,n)/dx(m,n) )
                    endif
                  enddo
                enddo
            enddo

            call pwdstd2D(path2data,'LOCAL/vbrtr.dat',nrec,array4_2d,block_lu_output,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,cart_comm,ierr)
            call fulfname(fname,path2data,'LOCAL/vbrtr.dat',ierr)
            if (rank == 0) then
              call ctl_file_write(fname,     &     !file name
                                  undef,     &     !value for undefined points
                                  nx_loc,    &     !x-dimension
                                  ny_loc,    &     !y-dimension
                                      1,     &     !z-dimension
                                    nrec,    &     !t-dimension
                                xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                                     xtm1,   &     !first x-value (if linear) or x-array (if levels)
                                    dxst,    &     !x-step (if linear)
                                ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                                     ytn1,   &     !first y-value (if linear) or x-array (if levels)
                                    dyst,    &     !y-step (if linear)
                                        0,   &     !z-grid type (0 - linear, 1 - levels)
                                      z0,    &     !first z-value (if linear) or x-array (if levels)
                                    1.0d0,   &     !z-step (if linear)
                                calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                                    year,    &     !year   of the first field
                                    month,   &     !month  of the first field
                                      day,   &     !day    of the first field
                                    hour,    &     !hour   of the first field
                                  minute,    &     !minute of the first field
                                    tstep,   &     !time step (in seconds
                'meridional velocity, m/s',  &     !title of dataset
                                      'v'   )      !variable name
            endif
        endif
    endif

    endsubroutine parallel_local_output


    subroutine parallel_global_output(path2data,  &
                                             nrec,       &
                                             year,       &
                                            month,       &
                                              day,       &
                                             hour,       &
                                           minute,       &
                                           tstep,        &
                                           calendar  )
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables
    use rec_length
    use iodata_routes
    use rw_ctl_routes

    implicit none
    include 'locout.fi'

    integer nrec, year, month, day, hour, minute, calendar, ierr
    character fname*256
    character*(*) path2data
    real(4) array4_2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
    real(4) tstep
    integer m,n,k
    real(8) z0(1), z1(1)

    z0 = 0.0d0; z1 = 1.0d0

    if (rank .eq. 0) write(*,*) 'Writing global output, record number ', nrec

    if (ssh_max_amplitude_output>0) then
    ! writing SSH
     ierr=0
     do k = 1, bcount
       call set_block(k)
       array4_2d(k)%vals = sngl(ssh_max_amplitude(k)%vals)
     enddo
     call pwdstd2D(path2data,'LOCAL/ssh_max_amplitude.dat',nrec,array4_2d,lu_output,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,cart_comm,ierr)
     call fulfname(fname,path2data,'LOCAL/ssh_max_amplitude.dat',ierr)
     if (rank .eq. 0) then
         call ctl_file_write(fname,    &     !file name
                             undef,    &     !value for undefined points
                             nx_loc,   &     !x-dimension
                             ny_loc,   &     !y-dimension
                                  1,   &     !z-dimension
                               nrec,   &     !t-dimension
                           xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                               xtm1,   &     !first x-value (if linear) or x-array (if levels)
                              dxst,    &     !x-step (if linear)
                          ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                               ytn1,   &     !first y-value (if linear) or x-array (if levels)
                              dyst,    &     !y-step (if linear)
                              0,       &     !z-grid type (0 - linear, 1 - levels)
                              z0,      &     !first z-value (if linear) or x-array (if levels)
                              1.0d0,   &     !z-step (if linear)
                           calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                               year,   &     !year   of the first field
                              month,   &     !month  of the first field
                                day,   &     !day    of the first field
                               hour,   &     !hour   of the first field
                             minute,   &     !minute of the first field
                              tstep,   &     !time step (in seconds)
                           'SSH, m',   &     !title of dataset
                              'ssh'   )      !variable name
     endif
    endif

    if (uv_max_amplitude_output>0) then
        if (grid_shift == 0) then
            ! writing on model grid
            ierr = 0
            do k = 1, bcount
              call set_block(k)
              array4_2d(k)%vals = sngl(ubrtr_max_amplitude(k)%vals)
            enddo
            call pwdstd2D(path2data,'LOCAL/ubrtr_max_amplitude.dat',nrec,array4_2d,llu_output,nx,ny,1,m1loc-1,m2loc,n1loc,n2loc,1,1,cart_comm,ierr)
            call fulfname(fname,path2data,'LOCAL/ubrtr_max_amplitude.dat',ierr)
            if (rank .eq. 0) then
                call ctl_file_write(fname,    &     !file name
                                    undef,    &     !value for undefined points
                                  nx_loc+1,   &     !x-dimension
                                    ny_loc,   &     !y-dimension
                                         1,   &     !z-dimension
                                      nrec,   &     !t-dimension
                                  xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                                      xum1,   &     !first x-value (if linear) or x-array (if levels)
                                     dxst,    &     !x-step (if linear)
                                 ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                                      ytn1,   &     !first y-value (if linear) or x-array (if levels)
                                     dyst,    &     !y-step (if linear)
                                     0,       &     !z-grid type (0 - linear, 1 - levels)
                                     z0,      &     !first z-value (if linear) or x-array (if levels)
                                     1.0d0,   &     !z-step (if linear)
                                  calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                                      year,   &     !year   of the first field
                                     month,   &     !month  of the first field
                                       day,   &     !day    of the first field
                                      hour,   &     !hour   of the first field
                                    minute,   &     !minute of the first field
                                     tstep,   &     !time step (in seconds
                      'zonal velocity, m/s',  &     !title of dataset
                                       'u'   )      !variable name
            endif
        else
            ! writing on T-grid
            do k = 1, bcount
                call set_block(k)
                call set_block_lu(k)
                call set_block_dxdy(k)
                call set_block_h(k)
                do n=ny_start, ny_end
                  do m=nx_start, nx_end
                    if (lu(m,n)>0.5) then
                        array4_2d(k)%vals(m, n) = sngl( (ubrtr_max_amplitude(k)%vals(m  ,n)*dyh(m  ,n)*hhu(m  ,n)   &
                                                        +ubrtr_max_amplitude(k)%vals(m-1,n)*dyh(m-1,n)*hhu(m-1,n) )/2.0/hhq(m,n)/dy(m,n) )
                    endif
                  enddo
                enddo
            enddo

            call pwdstd2D(path2data,'LOCAL/ubrtr_max_amplitude.dat',nrec,array4_2d,lu_output,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,cart_comm,ierr)
            call fulfname(fname,path2data,'LOCAL/ubrtr_max_amplitude.dat',ierr)
            if (rank == 0) then
              call ctl_file_write(fname,    &     !file name
                                  undef,    &     !value for undefined points
                                  nx_loc,   &     !x-dimension
                                  ny_loc,   &     !y-dimension
                                       1,   &     !z-dimension
                                    nrec,   &     !t-dimension
                                xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                                   xtm1,    &     !first x-value (if linear) or x-array (if levels)
                                  dxst,     &     !x-step (if linear)
                              ygr_type,     &     !y-grid type (0 - linear, 1 - levels)
                                   ytn1,    &     !first y-value (if linear) or x-array (if levels)
                                  dyst,     &     !y-step (if linear)
                                      0,    &     !z-grid type (0 - linear, 1 - levels)
                                      z0,   &     !first z-value (if linear) or x-array (if levels)
                                  1.0d0,    &     !z-step (if linear)
                                calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                                    year,   &     !year   of the first field
                                  month,    &     !month  of the first field
                                    day,    &     !day    of the first field
                                    hour,   &     !hour   of the first field
                                  minute,   &     !minute of the first field
                                  tstep,    &     !time step (in seconds
                    'zonal velocity, m/s',  &     !title of dataset
                                    'u'    )      !variable name
            endif
        endif

        if (grid_shift == 0) then
            ierr=0
            do k = 1, bcount
              call set_block(k)
              array4_2d(k)%vals = sngl(vbrtr_max_amplitude(k)%vals)
            enddo
            call pwdstd2D(path2data,'LOCAL/vbrtr_max_amplitude.dat',nrec,array4_2d,llv_output,nx,ny,1,m1loc,m2loc,n1loc-1,n2loc,1,1,cart_comm,ierr)
            call fulfname(fname,path2data,'LOCAL/vbrtr_max_amplitude.dat',ierr)
            if (rank .eq. 0) then
                call ctl_file_write(fname,    &     !file name
                                    undef,    &     !value for undefined points
                                    nx_loc,   &     !x-dimension
                                  ny_loc+1,   &     !y-dimension
                                         1,   &     !z-dimension
                                      nrec,   &     !t-dimension
                                  xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                                      xtm1,   &     !first x-value (if linear) or x-array (if levels)
                                     dxst,    &     !x-step (if linear)
                                 ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                                      yvn1,   &     !first y-value (if linear) or x-array (if levels)
                                     dyst,    &     !y-step (if linear)
                                     0,       &     !z-grid type (0 - linear, 1 - levels)
                                     z0,      &     !first z-value (if linear) or x-array (if levels)
                                     1.0d0,   &     !z-step (if linear)
                                  calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                                      year,   &     !year   of the first field
                                     month,   &     !month  of the first field
                                       day,   &     !day    of the first field
                                      hour,   &     !hour   of the first field
                                    minute,   &     !minute of the first field
                                     tstep,   &     !time step (in seconds
                 'meridional velocity, m/s',  &     !title of dataset
                                       'v'   )      !variable name
            endif
        else
            !writing on T-grid
            do k = 1, bcount
                call set_block(k)
                call set_block_lu(k)
                call set_block_dxdy(k)
                call set_block_h(k)
                do n=ny_start, ny_end
                  do m=nx_start, nx_end
                    if (lu(m,n)>0.5) then
                        array4_2d(k)%vals(m,n) = sngl( (vbrtr_max_amplitude(k)%vals(m,n  )*dxh(m,n  )*hhv(m,n  )    &
                                                       +vbrtr_max_amplitude(k)%vals(m,n-1)*dxh(m,n-1)*hhv(m,n-1))/2.0/hhq(m,n)/dx(m,n) )
                    endif
                  enddo
                enddo
            enddo

            call pwdstd2D(path2data,'LOCAL/vbrtr_max_amplitude.dat',nrec,array4_2d,lu_output,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,cart_comm,ierr)
            call fulfname(fname,path2data,'LOCAL/vbrtr_max_amplitude.dat',ierr)
            if (rank == 0) then
              call ctl_file_write(fname,     &     !file name
                                  undef,     &     !value for undefined points
                                  nx_loc,    &     !x-dimension
                                  ny_loc,    &     !y-dimension
                                      1,     &     !z-dimension
                                    nrec,    &     !t-dimension
                                xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                                     xtm1,   &     !first x-value (if linear) or x-array (if levels)
                                    dxst,    &     !x-step (if linear)
                                ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                                     ytn1,   &     !first y-value (if linear) or x-array (if levels)
                                    dyst,    &     !y-step (if linear)
                                        0,   &     !z-grid type (0 - linear, 1 - levels)
                                      z0,    &     !first z-value (if linear) or x-array (if levels)
                                    1.0d0,   &     !z-step (if linear)
                                calendar,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                                    year,    &     !year   of the first field
                                    month,   &     !month  of the first field
                                      day,   &     !day    of the first field
                                    hour,    &     !hour   of the first field
                                  minute,    &     !minute of the first field
                                    tstep,   &     !time step (in seconds
                'meridional velocity, m/s',  &     !title of dataset
                                      'v'   )      !variable name
            endif
        endif
    endif
    endsubroutine parallel_global_output

endmodule output_routes
