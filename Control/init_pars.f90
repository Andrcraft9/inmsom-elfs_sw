subroutine ocean_model_parameters(tau)
use main_basin_pars
use mpi_parallel_tools
use basin_grid
use ocean_variables
use key_switches
use iodata_routes
use rwpar_routes
use gridcon_routes
use basinpar_routes
use ocalg_routes
implicit none

character(256)    t_mask_file,       &  !name of file with temperature point sea-land mask
       bottom_topography_file,       &  !name of file with bottom topography
            t_mask_file_local,       &  !name of file with temperature point sea-land mask (LOCAL area)
            help_string
real(4) array4(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
integer  m, n, k, ierr
real(8) tau
real(8) :: hx2, hy2

! define parameters of task
! description of parameters see in file with mame filepar
 if (rank .eq. 0) then
     open (90,file='phys_proc.par',status='old')
     read(90,*) ksw_atmforc       !Atmospheric forcing (0 - no, 1 - yes)
     read(90,*) ksw_bfc           !Bottom friction (0 - no, 1 - yes)
     read(90,*) ksw_lat           !Lateral 2nd order mix parametrization (0 - no, 1 - yes)
     read(90,*) ksw_lat4          !Lateral 4nd order momentum mix parametrization (0 - no, 1 - yes)
     read(90,*) ksw_ssbc          !Type of surface boundary conditions (1 - surface T&S and wind stress are set; 2 - T&S fluxes and wind stress are set; 3 - T&S fluxes and wind stress are computed
     read(90,*) ksw_wflux         !normalize global mean salt balance (0 - no, 1 - normalize water flux, 2 - normalize salinity flux)
     read(90,*) ksw_lbc_ts        !open boundary conditions for T&S using (0 - no, 1 - yes)
     read(90,*) ksw_lbc_uv        !open boundary conditions for U&V using (0 - no, 1 - yes)
     read(90,*) ksw_lbc_ssh       !open boundary conditions for SSH using (0 - no, 1 - yes)

     read(90,*) lvisc_2           !lateral  vicosity(2nd order)[m**2/s]
     read(90,*) lvisc_4           !lateral  vicosity(4th order) [undimensional]
     read(90,*) nbfc              !Bottom friction coeff (Manning's roughness)

     help_string =' '
     read (90,'(a)') help_string   ! file with t-mask'
     call get_first_lexeme(help_string ,t_mask_file   )

     help_string =' '
     read (90,'(a)') help_string  ! file with bottom topography'
     call get_first_lexeme(help_string ,bottom_topography_file  )

     help_string =' '
     read (90,'(a)') help_string   ! file with t-mask'
     call get_first_lexeme(help_string ,t_mask_file_local   )

     close(90)
 endif

 call mpi_bcast(ksw_atmforc   , 1, mpi_integer, 0, cart_comm, ierr)
 call mpi_bcast(ksw_bfc       , 1, mpi_integer, 0, cart_comm, ierr)
 call mpi_bcast(ksw_lat       , 1, mpi_integer, 0, cart_comm, ierr)
 call mpi_bcast(ksw_lat4      , 1, mpi_integer, 0, cart_comm, ierr)
 call mpi_bcast(ksw_ssbc      , 1, mpi_integer, 0, cart_comm, ierr)
 call mpi_bcast(ksw_wflux     , 1, mpi_integer, 0, cart_comm, ierr)
 call mpi_bcast(ksw_lbc_ts    , 1, mpi_integer, 0, cart_comm, ierr)
 call mpi_bcast(ksw_lbc_uv    , 1, mpi_integer, 0, cart_comm, ierr)
 call mpi_bcast(ksw_lbc_ssh   , 1, mpi_integer, 0, cart_comm, ierr)

 call mpi_bcast(lvisc_2, 1, mpi_real8, 0, cart_comm, ierr)
 call mpi_bcast(lvisc_4, 1, mpi_real8, 0, cart_comm, ierr)
 call mpi_bcast(nbfc   , 1, mpi_real8, 0, cart_comm, ierr)
 
 call mpi_bcast(t_mask_file, 256, mpi_character, 0, cart_comm, ierr)
 call mpi_bcast(bottom_topography_file, 256, mpi_character, 0, cart_comm, ierr)
 call mpi_bcast(t_mask_file_local, 256, mpi_character, 0, cart_comm, ierr)

 if (rank .eq. 0) then
     write(*,'(i7,a)') ksw_atmforc,  ' - Atmospheric forcing'
     write(*,'(i7,a)') ksw_bfc,      ' - Bottom friction'
     write(*,'(i7,a)') ksw_lat,      ' - Lateral 2nd order mix parametrization'
     write(*,'(i7,a)') ksw_lat4,     ' - Lateral 4nd order momentum mix parametrization'
     write(*,'(i7,a)') ksw_ssbc,     ' - Type of surface boundary conditions'
     write(*,'(i7,a)') ksw_wflux,    ' - Normalize global mean salt balance'
     write(*,'(i7,a)') ksw_lbc_ts,   ' - Open boundary conditions for T&S using'
     write(*,'(i7,a)') ksw_lbc_uv,   ' - Open boundary conditions for U&V using'
     write(*,'(i7,a)') ksw_lbc_ssh,  ' - Open boundary conditions for SSH using'
     write(*,'(e12.4,a)') lvisc_2,       ' - Lateral  vicosity(2nd order)[m**2/s]'
     write(*,'(e12.4,a)') lvisc_4,       ' - Lateral  vicosity(4th order)[undim]'
     write(*,'(e12.4,a)') nbfc,          ' - Bottom friction coeff (Mannings roughness)'
     write(*,'(a,a)')  ' file with T-point sea-land mask: ', t_mask_file(1:len_trim (t_mask_file))
     write(*,'(a,a)')  '     file with bottom topography: ', bottom_topography_file(1:len_trim (bottom_topography_file))
     write(*,'(a,a)')  ' file with T-point sea-land mask for LOCAL area: ', t_mask_file_local(1:len_trim (t_mask_file_local))
 endif

   ! Init block grid with sea-land mask
   call parallel_read_mask(t_mask_file)
   call parallel_blocks_distribution

   ! Allocating main arrays
   call model_grid_allocate
   call ocean_variables_allocate

   ! area mask initialization
   call gridcon
   if (t_mask_file_local /= 'NONE') then
    call gridcon_local(t_mask_file_local)
   else
        lu_local = lu
        lcu_local = lcu
        lcv_local = lcv
        llu_local = llu
        llv_local = llv
   endif

   ! define grid geographical coordinates, steps and coriolis parameters
   call basinpar

   if (bottom_topography_file .eq. 'NONE') then
       if (rank .eq. 0) then
           print *, 'none topography !'
           print *, "HHQ_REST = 500m, topo file was ingored !"
       endif
       do k = 1, bcount
           hhq_rest(k)%vals = 500.0d0
       enddo
   else
       call allocate_block2D_real4(array4, 0.0)
       call prdstd2D(' ', bottom_topography_file, 1, array4, lu, nx, ny, mmm, mm, nnn, nn, ierr)
       do k = 1, bcount
           hhq_rest(k)%vals = dble(array4(k)%vals)
       enddo
       call deallocate_block2D_real4(array4)
   endif

   call syncborder_block2D_real8(hhq_rest)
   if(periodicity_x/=0) then
       call cyclize8_x(hhq_rest,nx,ny,1,mmm,mm)
   end if
   if(periodicity_y/=0) then
       call cyclize8_y(hhq_rest,nx,ny,1,nnn,nn)
   end if

    !--------------Rayleigh friction initialization
    !$omp parallel do
    do k = 1, bcount
        call set_block_boundary(k)
        do n=ny_start,ny_end
            do m=nx_start,nx_end
                if (lu(k)%vals(m,n)*lu(k)%vals(m+1,n)>0.5 .and. lu(k)%vals(m,n)*lu(k)%vals(m,n+1)>0.5) then
                    hx2= ( ((hhq_rest(k)%vals(m+1,n) - hhq_rest(k)%vals(m  ,n))/dxt(k)%vals(m  ,n))**2 * dble(lcu(k)%vals(m  ,n))    &
                          +((hhq_rest(k)%vals(m  ,n) - hhq_rest(k)%vals(m-1,n))/dxt(k)%vals(m-1,n))**2 * dble(lcu(k)%vals(m-1,n)) )  &
                          /dble(lcu(k)%vals(m,n)+lcu(k)%vals(m-1,n))

                    hy2= ( ((hhq_rest(k)%vals(m,n+1) - hhq_rest(k)%vals(m,n  ))/dyt(k)%vals(m,n  ))**2 * dble(lcv(k)%vals(m,n  ))    &
                          +((hhq_rest(k)%vals(m,n  ) - hhq_rest(k)%vals(m,n-1))/dyt(k)%vals(m,n-1))**2 * dble(lcv(k)%vals(m,n-1)) )  &
                          /dble(lcv(k)%vals(m,n)+lcv(k)%vals(m,n-1))

                    r_diss(k)%vals(m,n) = r_fric*dsqrt(hx2 + hy2)
                endif
            enddo
        enddo
    enddo
    !$omp end parallel do
    call syncborder_block2D_real8(r_diss)

   if(periodicity_x/=0) then
       call cyclize8_x(r_diss, nx, ny, 1, mmm, mm)
   end if
   if(periodicity_y/=0) then
       call cyclize8_y(r_diss, nx, ny, 1, nnn, nn)
   end if

endsubroutine ocean_model_parameters

subroutine test_init
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables
    implicit none

    integer :: k, m, n

    do k = 1, bcount
        call set_block_boundary(k)
        ssh(k)%vals = -1
        do m = nx_start, nx_end
            do n = ny_start, ny_end
                ssh(k)%vals(m, n) = rank
            enddo
        enddo
    enddo

end subroutine

subroutine zero_sw_init
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables
    use depth
    implicit none

    integer :: k

    do k = 1, bcount
        ubrtr(k)%vals = 0.0
        ubrtrp(k)%vals = 0.0

        vbrtr(k)%vals = 0.0
        vbrtrp(k)%vals = 0.0

        ssh(k)%vals = 0.0
        sshp(k)%vals = 0.0
    enddo

    do k = 1, bcount
        call set_block_boundary(k)
        call hh_init(hhq(k)%vals, hhqp(k)%vals, hhqn(k)%vals,     &
                     hhu(k)%vals, hhup(k)%vals, hhun(k)%vals,     &
                     hhv(k)%vals, hhvp(k)%vals, hhvn(k)%vals,     &
                     hhh(k)%vals, hhhp(k)%vals, hhhn(k)%vals,     &
                     ssh(k)%vals, sshp(k)%vals, hhq_rest(k)%vals, &
                     dx(k)%vals, dy(k)%vals, dxt(k)%vals, dyt(k)%vals,   &
                     dxh(k)%vals, dyh(k)%vals, dxb(k)%vals, dyb(k)%vals, &
                     lu(k)%vals, llu(k)%vals, llv(k)%vals, luh(k)%vals)
    enddo
    call syncborder_block2D_real8(hhu)
    call syncborder_block2D_real8(hhup)
    call syncborder_block2D_real8(hhun)
    call syncborder_block2D_real8(hhv)
    call syncborder_block2D_real8(hhvp)
    call syncborder_block2D_real8(hhvn)
    call syncborder_block2D_real8(hhh)
    call syncborder_block2D_real8(hhhp)
    call syncborder_block2D_real8(hhhn)

end subroutine

subroutine sw_only_inicond(flag_init, path2ocp)
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables
    use depth
    use iodata_routes
    implicit none
    integer :: flag_init
    character*(*) path2ocp
    type(block2D_real4), dimension(:), pointer :: array4

    integer :: k, ierr
    !real(4) array4(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

    do k = 1, bcount
        ubrtr(k)%vals = 0.0d0
        vbrtr(k)%vals = 0.0d0
    enddo

    ! Read init sea level
    if (flag_init > 0) then
        if (rank .eq. 0) print *, "Read init sea level"
        call allocate_block2D_real4(array4, 0.0)
        call prdstd2D(path2ocp, 'slf.dat', 1, array4, lu, nx, ny, mmm, mm, nnn, nn, ierr)
        do k = 1, bcount
            ssh(k)%vals = dble(array4(k)%vals)
        enddo
        call deallocate_block2D_real4(array4)
        call syncborder_block2D_real8(ssh)
    else
        if (rank.eq.0) print *, "Init sea level is zero"
        do k = 1, bcount
            ssh(k)%vals = 0.0d0
        enddo
    endif

    do k = 1, bcount
        ubrtrp(k)%vals = ubrtr(k)%vals
        vbrtrp(k)%vals = vbrtr(k)%vals
        sshp(k)%vals = ssh(k)%vals
    enddo

    do k = 1, bcount
        call set_block_boundary(k)
        call hh_init(hhq(k)%vals, hhqp(k)%vals, hhqn(k)%vals,     &
                     hhu(k)%vals, hhup(k)%vals, hhun(k)%vals,     &
                     hhv(k)%vals, hhvp(k)%vals, hhvn(k)%vals,     &
                     hhh(k)%vals, hhhp(k)%vals, hhhn(k)%vals,     &
                     ssh(k)%vals, sshp(k)%vals, hhq_rest(k)%vals, &
                     dx(k)%vals, dy(k)%vals, dxt(k)%vals, dyt(k)%vals,   &
                     dxh(k)%vals, dyh(k)%vals, dxb(k)%vals, dyb(k)%vals, &
                     lu(k)%vals, llu(k)%vals, llv(k)%vals, luh(k)%vals)
    enddo
    call syncborder_block2D_real8(hhu)
    call syncborder_block2D_real8(hhup)
    call syncborder_block2D_real8(hhun)
    call syncborder_block2D_real8(hhv)
    call syncborder_block2D_real8(hhvp)
    call syncborder_block2D_real8(hhvn)
    call syncborder_block2D_real8(hhh)
    call syncborder_block2D_real8(hhhp)
    call syncborder_block2D_real8(hhhn)

endsubroutine sw_only_inicond
