module output_routes
    implicit none

    contains

    subroutine init_model_output()
        use main_basin_pars
        use mpi_parallel_tools
        use basin_grid
        implicit none
        include 'locout.fi'

        integer :: ierr
        integer :: r, locr

        ! First point of area for output
        r = get_rank_by_point(m1loc, n1loc)
        if (r < 0) then
            print *, "Parallel output error! Wrong values of m1loc or n1loc"
            call mpi_abort(cart_comm, 1, ierr)
            stop
        endif
        if (rank .eq. r) then
            xtm1 = xt(m1loc)
            ytn1 = yt(n1loc)
            xum1 = xu(m1loc-1)
            yvn1 = yv(n1loc-1)
        endif
        call mpi_bcast(xtm1, 1, mpi_real8, r, cart_comm, ierr)
        call mpi_bcast(ytn1, 1, mpi_real8, r, cart_comm, ierr)
        call mpi_bcast(xum1, 1, mpi_real8, r, cart_comm, ierr)
        call mpi_bcast(yvn1, 1, mpi_real8, r, cart_comm, ierr)

        ! First point of local area for output
        locr = get_rank_by_point(m1loc_local, n1loc_local)
        if (locr < 0) then
            print *, "Parallel output error! Wrong values of m1loc_local or n1loc_local"
            call mpi_abort(cart_comm, 1, ierr)
            stop
        endif
        if (rank .eq. locr) then
            xtm1loc = xt(m1loc_local)
            ytn1loc = yt(n1loc_local)
            xum1loc = xu(m1loc_local-1)
            yvn1loc = yv(n1loc_local-1)
        endif
        call mpi_bcast(xtm1loc, 1, mpi_real8, locr, cart_comm, ierr)
        call mpi_bcast(ytn1loc, 1, mpi_real8, locr, cart_comm, ierr)
        call mpi_bcast(xum1loc, 1, mpi_real8, locr, cart_comm, ierr)
        call mpi_bcast(yvn1loc, 1, mpi_real8, locr, cart_comm, ierr)

    end subroutine

    subroutine parallel_check_point()
        use main_basin_pars
        use mpi_parallel_tools
        use basin_grid
        use ocean_variables

        implicit none
        include 'locout.fi'

        integer :: m, n, r, k
        real*8 :: lon, lat

        if (points_output > 1) then
            do k = 1, nloc_points
                lon = lon_loc_points(k)
                lat = lat_loc_points(k)

                m = floor((lon - rlon) / dxst) + mmm
                n = floor((lat - rlat) / dyst) + nnn

                r = get_rank_by_point(m, n)

                if (rank .eq. r) then
                    print *, 'rank ', ' lon ', ' lat ', ' geo_lon_t ', ' geo_lat_t', 'hhq_rest'
                    print *,  rank,   lon,   lat,   geo_lon_t(m, n), geo_lat_t(m, n), hhq_rest(m, n)
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
        use rw_ctl_routes

        implicit none
        include 'locout.fi'

        character fname*256
        character*(*) path2data
        real*8 :: point_time
        integer :: m, n, r, ierr, k
        real*8 :: lon, lat
        real*8 :: uuu, vvv

        if (points_output > 0) then
            do k = 1, nloc_points
                lon = lon_loc_points(k)
                lat = lat_loc_points(k)

                m = floor((lon - rlon) / dxst) + mmm
                n = floor((lat - rlat) / dyst) + nnn

                r = get_rank_by_point(m, n)

                if (rank .eq. r) then
                    uuu = ( ubrtr(m  ,n)*dyh(m  ,n)*hhu(m  ,n)    &
                           +ubrtr(m-1,n)*dyh(m-1,n)*hhu(m-1,n) )/2.0/hhq(m,n)/dy(m,n)
                    vvv = ( vbrtr(m,n  )*dxh(m,n  )*hhv(m,n  )    &
                           +vbrtr(m,n-1)*dxh(m,n-1)*hhv(m,n-1) )/2.0/hhq(m,n)/dx(m,n)

                    call fulfname(fname, path2data, name_points(k), ierr)
                    open(40, file=fname, status='unknown', position='append')
                    write(40, *) point_time, ssh(m, n), uuu, vvv
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
        use rw_ctl_routes

        implicit none
        include 'locout.fi'

        character fname*256
        character*(*) path2data
        real*8 :: point_time
        integer :: m, n, ierr
        real*8 :: kinetic_e, potential_e
        real*8 :: buf_k, buf_p

        if (energy_output > 0) then
            ! TotalEnergy = 0.5 rho dx dy Sum[h_u u**2 + h_v v**2 + g ssh**2]
            kinetic_e = 0.0d0; potential_e = 0.0d0
            do n = ny_start, ny_end
                do m = nx_start, nx_end
                    if (lcu(m,n)>0.5) then
                        kinetic_e = kinetic_e + 0.5d0*RefDen*dxt(m,n)*dyh(m,n)*hhu(m,n)*(ubrtr(m,n)**2)
                    endif
                    if (lcv(m,n)>0.5) then
                        kinetic_e = kinetic_e + 0.5d0*RefDen*dxh(m,n)*dyt(m,n)*hhv(m,n)*(vbrtr(m,n)**2)
                    endif

                    if (lu(m,n)>0.5) then
                        potential_e = potential_e + 0.5d0*RefDen*dx(m,n)*dy(m,n)*FreeFallAcc*(ssh(m,n)**2)
                    endif
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

        if (energy_local_output > 0) then
            ! TotalEnergy = 0.5 rho dx dy Sum[h_u u**2 + h_v v**2 + g ssh**2]
            kinetic_e = 0.0d0; potential_e = 0.0d0
            do n = ny_start, ny_end
                do m = nx_start, nx_end
                    if (lcu_local(m,n)>0.5) then
                        kinetic_e = kinetic_e + 0.5d0*RefDen*dxt(m,n)*dyh(m,n)*hhu(m,n)*(ubrtr(m,n)**2)
                    endif
                    if (lcv_local(m,n)>0.5) then
                        kinetic_e = kinetic_e + 0.5d0*RefDen*dxh(m,n)*dyt(m,n)*hhv(m,n)*(vbrtr(m,n)**2)
                    endif

                    if (lu_local(m,n)>0.5) then
                        potential_e = potential_e + 0.5d0*RefDen*dx(m,n)*dy(m,n)*FreeFallAcc*(ssh(m,n)**2)
                    endif
                enddo
            enddo
            buf_k = kinetic_e; buf_p = potential_e
            call mpi_allreduce(buf_k, kinetic_e, 1, mpi_real8, mpi_sum, cart_comm, ierr)
            call mpi_allreduce(buf_p, potential_e, 1, mpi_real8, mpi_sum, cart_comm, ierr)

            if (rank .eq. 0) then
                call fulfname(fname, path2data, 'kinetic_potential_energy_local', ierr)
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
    use main_basin_pars
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

    if (rank .eq. 0) write(*,*) 'Writing local output, record number ', nrec

    if(nrec==1) then
     if (rank == 0) then
         print *, "first x-value: ",  xtm1, "first y-value: ",  ytn1
     endif
     !writing HHQ
     ierr=0
     array4_2d=sngl(hhq_rest)
     !call wdstd(path2data,'LOCAL/hhq.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
     call pwdstd(path2data,'LOCAL/hhq.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,cart_comm,ierr)
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
     array4_2d=sngl(ssh)
    ! call wdstd(path2data,'LOCAL/ssh.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
     call pwdstd(path2data,'LOCAL/ssh.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,cart_comm,ierr)
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
            array4_2d=sngl(ubrtr)
            call pwdstd(path2data,'LOCAL/ubrtr.dat',nrec,array4_2d,llu,nx,ny,1,m1loc-1,m2loc,n1loc,n2loc,1,1,cart_comm,ierr)
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
            !$omp parallel do 
            do n=ny_start, ny_end
              do m=nx_start, nx_end
                if (lu(m,n)>0.5) then
                    array4_2d(m, n) = sngl( (ubrtr(m  ,n)*dyh(m  ,n)*hhu(m  ,n)   &
                                            +ubrtr(m-1,n)*dyh(m-1,n)*hhu(m-1,n) )/2.0/hhq(m,n)/dy(m,n) ) 
                endif
              enddo
            enddo
            !$omp end parallel do
          
            call pwdstd(path2data,'LOCAL/ubrtr.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,cart_comm,ierr)
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
            array4_2d=sngl(vbrtr)
            call pwdstd(path2data,'LOCAL/vbrtr.dat',nrec,array4_2d,llv,nx,ny,1,m1loc,m2loc,n1loc-1,n2loc,1,1,cart_comm,ierr)
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
            !$omp parallel do 
            do n=ny_start, ny_end
              do m=nx_start, nx_end
                if (lu(m,n)>0.5) then
                    array4_2d(m,n) = sngl( (vbrtr(m,n  )*dxh(m,n  )*hhv(m,n  )    &
                                           +vbrtr(m,n-1)*dxh(m,n-1)*hhv(m,n-1))/2.0/hhq(m,n)/dx(m,n) )
                endif
              enddo
            enddo
            !$omp end parallel do
        
            call pwdstd(path2data,'LOCAL/vbrtr.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,cart_comm,ierr)
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

    ! ---------- LOCAL OUTPUT -----------------!
    !------------------------------------------!
    !------------------------------------------!
    if (ssh_local_output>0) then
        ! writing SSH
         ierr=0
         array4_2d=sngl(ssh)
         call pwdstd(path2data,'LOCAL/sshloc.dat',nrec,array4_2d,lu_local,nx,ny,1,  &
                     m1loc_local,m2loc_local,n1loc_local,n2loc_local,1,1, local_output_comm, ierr)
         call fulfname(fname,path2data,'LOCAL/sshloc.dat',ierr)
         if (rank .eq. 0) then
             call ctl_file_write(fname,    &     !file name
                             undef,    &     !value for undefined points
                             nx_loc_local,   &     !x-dimension
                             ny_loc_local,   &     !y-dimension
                                  1,   &     !z-dimension
                               nrec,   &     !t-dimension
                           xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                            xtm1loc,   &     !first x-value (if linear) or x-array (if levels)
                              dxst,    &     !x-step (if linear)
                          ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                            ytn1loc,   &     !first y-value (if linear) or x-array (if levels)
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

    if (uv_local_output>0) then
        if (grid_shift == 0) then
            ! writing on model grid
            ierr = 0
            array4_2d=sngl(ubrtr)
            call pwdstd(path2data,'LOCAL/ubrtrloc.dat',nrec,array4_2d,llu_local,nx,ny,1,  &
                        m1loc_local-1,m2loc_local,n1loc_local,n2loc_local,1,1, local_output_comm, ierr)
            call fulfname(fname,path2data,'LOCAL/ubrtrloc.dat',ierr)
            if (rank .eq. 0) then
                call ctl_file_write(fname,    &     !file name
                                    undef,    &     !value for undefined points
                                    nx_loc_local + 1,   &     !x-dimension
                                    ny_loc_local,   &     !y-dimension
                                            1,   &     !z-dimension
                                         nrec,   &     !t-dimension
                                     xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                                      xum1loc,   &     !first x-value (if linear) or x-array (if levels)
                                        dxst,    &     !x-step (if linear)
                                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                                      ytn1loc,   &     !first y-value (if linear) or x-array (if levels)
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
            !$omp parallel do 
            do n=ny_start, ny_end
                do m=nx_start, nx_end
                if (lu_local(m,n)>0.5) then
                    array4_2d(m, n) = sngl( (ubrtr(m  ,n)*dyh(m  ,n)*hhu(m  ,n)   &
                                            +ubrtr(m-1,n)*dyh(m-1,n)*hhu(m-1,n) )/2.0/hhq(m,n)/dy(m,n) ) 
                endif
                enddo
            enddo
            !$omp end parallel do
            
            call pwdstd(path2data,'LOCAL/ubrtrloc.dat',nrec,array4_2d,lu_local,nx,ny,1,  &
                        m1loc_local,m2loc_local,n1loc_local,n2loc_local,1,1, local_output_comm, ierr)
            call fulfname(fname,path2data,'LOCAL/ubrtrloc.dat',ierr)
            if (rank == 0) then
                call ctl_file_write(fname,    &     !file name
                                    undef,    &     !value for undefined points
                                    nx_loc_local,   &     !x-dimension
                                    ny_loc_local,   &     !y-dimension
                                         1,   &     !z-dimension
                                      nrec,   &     !t-dimension
                                  xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                                  xtm1loc,    &     !first x-value (if linear) or x-array (if levels)
                                    dxst,     &     !x-step (if linear)
                                ygr_type,     &     !y-grid type (0 - linear, 1 - levels)
                                  ytn1loc,    &     !first y-value (if linear) or x-array (if levels)
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
            array4_2d=sngl(vbrtr)
            call pwdstd(path2data,'LOCAL/vbrtrloc.dat',nrec,array4_2d,llv_local,nx,ny,1,  &
                        m1loc_local,m2loc_local,n1loc_local-1,n2loc_local,1,1, local_output_comm, ierr)
            call fulfname(fname,path2data,'LOCAL/vbrtrloc.dat',ierr)
            if (rank .eq. 0) then
                call ctl_file_write(fname,    &     !file name
                                    undef,    &     !value for undefined points
                                    nx_loc_local,   &     !x-dimension
                                    ny_loc_local + 1,   &     !y-dimension
                                            1,   &     !z-dimension
                                         nrec,   &     !t-dimension
                                     xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                                      xtm1loc,   &     !first x-value (if linear) or x-array (if levels)
                                        dxst,    &     !x-step (if linear)
                                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                                      yvn1loc,   &     !first y-value (if linear) or x-array (if levels)
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
            !$omp parallel do 
            do n=ny_start, ny_end
                do m=nx_start, nx_end
                if (lu_local(m,n)>0.5) then
                    array4_2d(m,n) = sngl( ( vbrtr(m,n  )*dxh(m,n  )*hhv(m,n  )    &
                                            +vbrtr(m,n-1)*dxh(m,n-1)*hhv(m,n-1))/2.0/hhq(m,n)/dx(m,n) )
                endif
                enddo
            enddo
            !$omp end parallel do
        
            call pwdstd(path2data,'LOCAL/vbrtrloc.dat',nrec,array4_2d,lu_local,nx,ny,1,  &
                        m1loc_local,m2loc_local,n1loc_local,n2loc_local,1,1, local_output_comm, ierr)
            call fulfname(fname,path2data,'LOCAL/vbrtrloc.dat',ierr)
            if (rank == 0) then
                call ctl_file_write(fname,     &     !file name
                                    undef,     &     !value for undefined points
                                    nx_loc_local,    &     !x-dimension
                                    ny_loc_local,    &     !y-dimension
                                      1,     &     !z-dimension
                                    nrec,    &     !t-dimension
                                xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                                  xtm1loc,   &     !first x-value (if linear) or x-array (if levels)
                                    dxst,    &     !x-step (if linear)
                                ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                                  ytn1loc,   &     !first y-value (if linear) or x-array (if levels)
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
    use main_basin_pars
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
     array4_2d=sngl(ssh_max_amplitude)
    ! call wdstd(path2data,'LOCAL/ssh.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
     call pwdstd(path2data,'LOCAL/ssh_max_amplitude.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,cart_comm,ierr)
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
            array4_2d=sngl(ubrtr_max_amplitude)
            call pwdstd(path2data,'LOCAL/ubrtr_max_amplitude.dat',nrec,array4_2d,llu,nx,ny,1,m1loc-1,m2loc,n1loc,n2loc,1,1,cart_comm,ierr)
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
            !$omp parallel do 
            do n=ny_start, ny_end
              do m=nx_start, nx_end
                if (lu(m,n)>0.5) then
                    array4_2d(m, n) = sngl( (ubrtr_max_amplitude(m  ,n)*dyh(m  ,n)*hhu(m  ,n)   &
                                            +ubrtr_max_amplitude(m-1,n)*dyh(m-1,n)*hhu(m-1,n) )/2.0/hhq(m,n)/dy(m,n) ) 
                endif
              enddo
            enddo
            !$omp end parallel do
          
            call pwdstd(path2data,'LOCAL/ubrtr_max_amplitude.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,cart_comm,ierr)
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
            array4_2d=sngl(vbrtr_max_amplitude)
            call pwdstd(path2data,'LOCAL/vbrtr_max_amplitude.dat',nrec,array4_2d,llv,nx,ny,1,m1loc,m2loc,n1loc-1,n2loc,1,1,cart_comm,ierr)
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
            !$omp parallel do 
            do n=ny_start, ny_end
              do m=nx_start, nx_end
                if (lu(m,n)>0.5) then
                    array4_2d(m,n) = sngl( (vbrtr_max_amplitude(m,n  )*dxh(m,n  )*hhv(m,n  )    &
                                           +vbrtr_max_amplitude(m,n-1)*dxh(m,n-1)*hhv(m,n-1))/2.0/hhq(m,n)/dx(m,n) )
                endif
              enddo
            enddo
            !$omp end parallel do
        
            call pwdstd(path2data,'LOCAL/vbrtr_max_amplitude.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,cart_comm,ierr)
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

    ! ---------- LOCAL OUTPUT -----------------!
    !------------------------------------------!
    !------------------------------------------!
    if (ssh_max_amplitude_local_output>0) then
        ! writing SSH
         ierr=0
         array4_2d=sngl(ssh_max_amplitude)
         call pwdstd(path2data,'LOCAL/sshloc_max_amplitude.dat',nrec,array4_2d,lu_local,nx,ny,1,  &
                     m1loc_local,m2loc_local,n1loc_local,n2loc_local,1,1, local_output_comm, ierr)
         call fulfname(fname,path2data,'LOCAL/sshloc_max_amplitude.dat',ierr)
         if (rank .eq. 0) then
             call ctl_file_write(fname,    &     !file name
                             undef,    &     !value for undefined points
                             nx_loc_local,   &     !x-dimension
                             ny_loc_local,   &     !y-dimension
                                  1,   &     !z-dimension
                               nrec,   &     !t-dimension
                           xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                            xtm1loc,   &     !first x-value (if linear) or x-array (if levels)
                              dxst,    &     !x-step (if linear)
                          ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                            ytn1loc,   &     !first y-value (if linear) or x-array (if levels)
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

    if (uv_max_amplitude_local_output>0) then
        if (grid_shift == 0) then
            ! writing on model grid
            ierr = 0
            array4_2d=sngl(ubrtr_max_amplitude)
            call pwdstd(path2data,'LOCAL/ubrtrloc_max_amplitude.dat',nrec,array4_2d,llu_local,nx,ny,1,  &
                        m1loc_local-1,m2loc_local,n1loc_local,n2loc_local,1,1, local_output_comm, ierr)
            call fulfname(fname,path2data,'LOCAL/ubrtrloc_max_amplitude.dat',ierr)
            if (rank .eq. 0) then
                call ctl_file_write(fname,    &     !file name
                                    undef,    &     !value for undefined points
                                    nx_loc_local + 1,   &     !x-dimension
                                    ny_loc_local,   &     !y-dimension
                                            1,   &     !z-dimension
                                         nrec,   &     !t-dimension
                                     xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                                      xum1loc,   &     !first x-value (if linear) or x-array (if levels)
                                        dxst,    &     !x-step (if linear)
                                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                                      ytn1loc,   &     !first y-value (if linear) or x-array (if levels)
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
            !$omp parallel do 
            do n=ny_start, ny_end
                do m=nx_start, nx_end
                if (lu_local(m,n)>0.5) then
                    array4_2d(m, n) = sngl( (ubrtr_max_amplitude(m  ,n)*dyh(m  ,n)*hhu(m  ,n)   &
                                            +ubrtr_max_amplitude(m-1,n)*dyh(m-1,n)*hhu(m-1,n) )/2.0/hhq(m,n)/dy(m,n) ) 
                endif
                enddo
            enddo
            !$omp end parallel do
            
            call pwdstd(path2data,'LOCAL/ubrtrloc_max_amplitude.dat',nrec,array4_2d,lu_local,nx,ny,1,  &
                        m1loc_local,m2loc_local,n1loc_local,n2loc_local,1,1, local_output_comm, ierr)
            call fulfname(fname,path2data,'LOCAL/ubrtrloc_max_amplitude.dat',ierr)
            if (rank == 0) then
                call ctl_file_write(fname,    &     !file name
                                    undef,    &     !value for undefined points
                                    nx_loc_local,   &     !x-dimension
                                    ny_loc_local,   &     !y-dimension
                                         1,   &     !z-dimension
                                      nrec,   &     !t-dimension
                                  xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                                  xtm1loc,    &     !first x-value (if linear) or x-array (if levels)
                                    dxst,     &     !x-step (if linear)
                                ygr_type,     &     !y-grid type (0 - linear, 1 - levels)
                                  ytn1loc,    &     !first y-value (if linear) or x-array (if levels)
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
            array4_2d=sngl(vbrtr_max_amplitude)
            call pwdstd(path2data,'LOCAL/vbrtrloc_max_amplitude.dat',nrec,array4_2d,llv_local,nx,ny,1,  &
                        m1loc_local,m2loc_local,n1loc_local-1,n2loc_local,1,1, local_output_comm, ierr)
            call fulfname(fname,path2data,'LOCAL/vbrtrloc_max_amplitude.dat',ierr)
            if (rank .eq. 0) then
                call ctl_file_write(fname,    &     !file name
                                    undef,    &     !value for undefined points
                                    nx_loc_local,   &     !x-dimension
                                    ny_loc_local + 1,   &     !y-dimension
                                            1,   &     !z-dimension
                                         nrec,   &     !t-dimension
                                     xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                                      xtm1loc,   &     !first x-value (if linear) or x-array (if levels)
                                        dxst,    &     !x-step (if linear)
                                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                                      yvn1loc,   &     !first y-value (if linear) or x-array (if levels)
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
            !$omp parallel do 
            do n=ny_start, ny_end
                do m=nx_start, nx_end
                if (lu_local(m,n)>0.5) then
                    array4_2d(m,n) = sngl( (vbrtr_max_amplitude(m,n  )*dxh(m,n  )*hhv(m,n  )    &
                                           +vbrtr_max_amplitude(m,n-1)*dxh(m,n-1)*hhv(m,n-1))/2.0/hhq(m,n)/dx(m,n) )
                endif
                enddo
            enddo
            !$omp end parallel do
        
            call pwdstd(path2data,'LOCAL/vbrtrloc_max_amplitude.dat',nrec,array4_2d,lu_local,nx,ny,1,  &
                        m1loc_local,m2loc_local,n1loc_local,n2loc_local,1,1, local_output_comm, ierr)
            call fulfname(fname,path2data,'LOCAL/vbrtrloc_max_amplitude.dat',ierr)
            if (rank == 0) then
                call ctl_file_write(fname,     &     !file name
                                    undef,     &     !value for undefined points
                                    nx_loc_local,    &     !x-dimension
                                    ny_loc_local,    &     !y-dimension
                                      1,     &     !z-dimension
                                    nrec,    &     !t-dimension
                                xgr_type,    &     !x-grid type (0 - linear, 1 - levels)
                                  xtm1loc,   &     !first x-value (if linear) or x-array (if levels)
                                    dxst,    &     !x-step (if linear)
                                ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                                  ytn1loc,   &     !first y-value (if linear) or x-array (if levels)
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
