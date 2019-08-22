module ocpar_routes
    implicit none

contains
    subroutine model_finalize()
        use mpi_parallel_tools
        use init_arrays_routes
        implicit none
        ! deallocating the arrays
        call ocean_variables_deallocate
        call model_grid_deallocate
        call atm_arrays_deallocate
    end subroutine

    subroutine ocean_model_parameters(tau)
        use mpi_parallel_tools
        use init_arrays_routes
        use grid_parameters
        use basin_grid
        use ocean_variables
        use key_switches
        use iodata_routes
        use rwpar_routes
        use gridcon_routes
        use ocalg_routes
        implicit none

        character(256) t_mask_file, &  !name of file with temperature point sea-land mask
                       bottom_topography_file, &  !name of file with bottom topography
                       t_mask_file_local, &  !name of file with temperature point sea-land mask (LOCAL area)
                       help_string
        integer m, n, k, ierr
        real(8) tau
        real(8) hx2, hy2

        ! define parameters of task
        ! description of parameters see in file with mame filepar
        if (rank .eq. 0) then
            open (90, file='phys_proc.par', status='old')
            read (90, *) ksw_atmforc       !Atmospheric forcing (0 - no, 1 - yes)
            read (90, *) ksw_bfc           !Bottom friction (0 - no, 1 - yes)
            read (90, *) ksw_lat           !Lateral 2nd order mix parametrization (0 - no, 1 - yes)
            read (90, *) ksw_lat4          !Lateral 4nd order momentum mix parametrization (0 - no, 1 - yes)
            read (90, *) ksw_ssbc          !Type of surface boundary conditions (1 - surface T&S and wind stress are set; 2 - T&S fluxes and wind stress are set; 3 - T&S fluxes and wind stress are computed
            read (90, *) ksw_wflux         !normalize global mean salt balance (0 - no, 1 - normalize water flux, 2 - normalize salinity flux)
            read (90, *) ksw_lbc_ts        !open boundary conditions for T&S using (0 - no, 1 - yes)
            read (90, *) ksw_lbc_uv        !open boundary conditions for U&V using (0 - no, 1 - yes)
            read (90, *) ksw_lbc_ssh       !open boundary conditions for SSH using (0 - no, 1 - yes)

            read (90, *) lvisc_2           !lateral  vicosity(2nd order)[m**2/s]
            read (90, *) lvisc_4           !lateral  vicosity(4th order) [undimensional]
            read (90, *) nbfc              !Bottom friction coeff (Manning's roughness)

            help_string = ' '
            read (90, '(a)') help_string   ! file with t-mask'
            call get_first_lexeme(help_string, t_mask_file)

            help_string = ' '
            read (90, '(a)') help_string  ! file with bottom topography'
            call get_first_lexeme(help_string, bottom_topography_file)

            help_string = ' '
            read (90, '(a)') help_string   ! file with t-mask'
            call get_first_lexeme(help_string, t_mask_file_local)

            close (90)
        endif

        call mpi_bcast(ksw_atmforc, 1, mpi_integer, 0, cart_comm, ierr)
        call mpi_bcast(ksw_bfc, 1, mpi_integer, 0, cart_comm, ierr)
        call mpi_bcast(ksw_lat, 1, mpi_integer, 0, cart_comm, ierr)
        call mpi_bcast(ksw_lat4, 1, mpi_integer, 0, cart_comm, ierr)
        call mpi_bcast(ksw_ssbc, 1, mpi_integer, 0, cart_comm, ierr)
        call mpi_bcast(ksw_wflux, 1, mpi_integer, 0, cart_comm, ierr)
        call mpi_bcast(ksw_lbc_ts, 1, mpi_integer, 0, cart_comm, ierr)
        call mpi_bcast(ksw_lbc_uv, 1, mpi_integer, 0, cart_comm, ierr)
        call mpi_bcast(ksw_lbc_ssh, 1, mpi_integer, 0, cart_comm, ierr)

        call mpi_bcast(lvisc_2, 1, mpi_real8, 0, cart_comm, ierr)
        call mpi_bcast(lvisc_4, 1, mpi_real8, 0, cart_comm, ierr)
        call mpi_bcast(nbfc, 1, mpi_real8, 0, cart_comm, ierr)

        call mpi_bcast(t_mask_file, 256, mpi_character, 0, cart_comm, ierr)
        call mpi_bcast(bottom_topography_file, 256, mpi_character, 0, cart_comm, ierr)
        call mpi_bcast(t_mask_file_local, 256, mpi_character, 0, cart_comm, ierr)

        if (rank .eq. 0) then
            write (*, '(i7,a)') ksw_atmforc, ' - Atmospheric forcing'
            write (*, '(i7,a)') ksw_bfc, ' - Bottom friction'
            write (*, '(i7,a)') ksw_lat, ' - Lateral 2nd order mix parametrization'
            write (*, '(i7,a)') ksw_lat4, ' - Lateral 4nd order momentum mix parametrization'
            write (*, '(i7,a)') ksw_ssbc, ' - Type of surface boundary conditions'
            write (*, '(i7,a)') ksw_wflux, ' - Normalize global mean salt balance'
            write (*, '(i7,a)') ksw_lbc_ts, ' - Open boundary conditions for T&S using'
            write (*, '(i7,a)') ksw_lbc_uv, ' - Open boundary conditions for U&V using'
            write (*, '(i7,a)') ksw_lbc_ssh, ' - Open boundary conditions for SSH using'
            write (*, '(e12.4,a)') lvisc_2, ' - Lateral  vicosity(2nd order)[m**2/s]'
            write (*, '(e12.4,a)') lvisc_4, ' - Lateral  vicosity(4th order)[undim]'
            write (*, '(e12.4,a)') nbfc, ' - Bottom friction coeff (Mannings roughness)'
            write (*, '(a,a)') ' file with T-point sea-land mask: ', t_mask_file(1:len_trim(t_mask_file))
            write (*, '(a,a)') '     file with bottom topography: ', bottom_topography_file(1:len_trim(bottom_topography_file))
            write (*, '(a,a)') ' file with T-point sea-land mask for LOCAL area: ', t_mask_file_local(1:len_trim(t_mask_file_local))
        endif

        ! Read sea-land mask, init blocks grid and distribute blocks to processors according to good load-balancing property
        call parallel_read_mask(t_mask_file)
        call parallel_blocks_distribution

        ! Allocate main data arrays
        call model_grid_allocate
        call ocean_variables_allocate
        call atm_arrays_allocate

        ! Initialization of varios masks of basin area
        call gridcon

        ! Read regional mask if necessary
        if (t_mask_file_local /= 'NONE') then
            call gridcon_output(t_mask_file_local)
        else
            do k = 1, bcount
                call set_block(k)
                call set_block_lu(k)
                call set_block_lu_output(k)
                lu_output = lu
                lcu_output = lcu
                lcv_output = lcv
                llu_output = llu
                llv_output = llv
            enddo
        endif

        ! Define grid geographical coordinates, steps and coriolis parameters
        call basinpar

        if (bottom_topography_file .eq. 'NONE') then
            if (rank .eq. 0) then
                print *, 'Topography is not given!'
                print *, "Set depth equals to 500m for entire basin area!"
            endif
            do k = 1, bcount
                call set_block(k)
                call set_block_lu(k)
                call set_block_h(k)
                hhq_rest = 500.0d0
            enddo
        else
            call prdstd2D(' ', bottom_topography_file, 1, array4_2d, block_lu, nx, ny, mmm, mm, nnn, nn, ierr)
            do k = 1, bcount
                call set_block(k)
                call set_block_lu(k)
                call set_block_h(k)
                hhq_rest = dble(array4_2d(k)%vals)
            enddo
        endif
        call syncborder_block2D_real8(block_hhq_rest)

        ! Rayleigh friction initialization
        do k = 1, bcount
            call set_block(k)
            call set_block_lu(k)
            call set_block_h(k)
            call set_block_dxdy(k)
            do n = ny_start,ny_end
                do m = nx_start,nx_end
                    if (lu(m,n)>0.5) then
                        hx2=  ((hhq_rest(m+1,n)-hhq_rest(m  ,n))/dxt(m  ,n))**2     &
                             +((hhq_rest(m  ,n)-hhq_rest(m-1,n))/dxt(m-1,n))**2
                        hy2=  ((hhq_rest(m,n+1)-hhq_rest(m,n  ))/dyt(m,n  ))**2     &
                             +((hhq_rest(m,n  )-hhq_rest(m,n-1))/dyt(m,n-1))**2
                        r_diss(k)%vals(m,n) = r_fric*dble( sqrt((hx2+hy2)/2.0) )
                    endif
                enddo
            enddo
        enddo
        call syncborder_block2D_real8(r_diss)

    endsubroutine ocean_model_parameters

    ! Init sea level, zero velocities.
    ! If flag_init > 0 then subroutine reads file slf.dat in results dir
    ! else subroutine sets sea level (ssh) equals to zero
    subroutine sw_only_inicond(flag_init, path2ocp)
        use mpi_parallel_tools
        use basin_grid
        use ocean_variables
        use depth_routes
        use iodata_routes
        implicit none

        integer :: flag_init
        character*(*) path2ocp
        integer :: k, ierr

        do k = 1, bcount
            call set_block(k)
            ubrtr(k)%vals = 0.0d0
            vbrtr(k)%vals = 0.0d0
        enddo

        ! Read init sea level
        if (flag_init > 0) then
            if (rank .eq. 0) print *, "Read init sea level"
            call prdstd2D(path2ocp, 'slf.dat', 1, array4_2d, block_lu, nx, ny, mmm, mm, nnn, nn, ierr)
            do k = 1, bcount
                call set_block(k)
                ssh(k)%vals = dble(array4_2d(k)%vals)
            enddo
            call syncborder_block2D_real8(ssh)

        else
            if (rank .eq. 0) print *, "Init sea level is zero"
            do k = 1, bcount
                ssh(k)%vals = 0.0d0
            enddo
        endif

        ! Set previous time step equls initial conditions
        do k = 1, bcount
            call set_block(k)
            ubrtrp(k)%vals = ubrtr(k)%vals
            vbrtrp(k)%vals = vbrtr(k)%vals
            sshp(k)%vals = ssh(k)%vals
        enddo

        do k = 1, bcount
            call set_block(k)
            call set_block_lu(k)
            call set_block_h(k)
            call set_block_dxdy(k)
            call hh_init(hhq, hhqp, hhqn, &
                         hhu, hhup, hhun, &
                         hhv, hhvp, hhvn, &
                         hhh, hhhp, hhhn, &
                         ssh(k)%vals, sshp(k)%vals, &
                         hhq_rest)
        enddo
        call syncborder_block2D_real8(block_hhu)
        call syncborder_block2D_real8(block_hhup)
        call syncborder_block2D_real8(block_hhun)
        call syncborder_block2D_real8(block_hhv)
        call syncborder_block2D_real8(block_hhvp)
        call syncborder_block2D_real8(block_hhvn)
        call syncborder_block2D_real8(block_hhh)
        call syncborder_block2D_real8(block_hhhp)
        call syncborder_block2D_real8(block_hhhn)

    endsubroutine sw_only_inicond

end module ocpar_routes
