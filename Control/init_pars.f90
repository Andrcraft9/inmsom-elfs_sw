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

! igrzts_surf  = min(IABS(ksw_ssbc),2) ! type of condition for T and S on sea surface
! igrzts_bot = 2                   ! type of condition for T and S on sea bottom

   ! area mask initialization
   call gridcon(t_mask_file)
   if (t_mask_file_local /= 'NONE') then
        call gridcon_local(t_mask_file_local)
   else
        lu_local = lu
        lcu_local = lcu
        lcv_local = lcv
        llu_local = llu
        llv_local = llv
   endif
   !if (rank .eq. 0) print *, "--------------------END OF GRIDCON----------------------"

   ! define grid geographical coordinates, steps and coriolis parameters
   call basinpar
   !if (rank .eq. 0) print *, "--------------------END OF BASINPAR----------------------"

   if (bottom_topography_file .eq. 'NONE') then
       if (rank .eq. 0) print *, 'NONE topography !'
       hhq_rest = 500.0d0
       if (rank .eq. 0) print *, "!!! HHQ_REST = 500m, topo file was ingored !!!"
   else
       array4=0.0
       call prdstd(' ',bottom_topography_file,1,array4,lu,nx,ny,1, mmm,mm,nnn,nn,1,1,ierr)
       hhq_rest=dble(array4)
   endif

   call syncborder_real8(hhq_rest, 1)
   if(periodicity_x/=0) then
       call cyclize8_x(hhq_rest,nx,ny,1,mmm,mm)
   end if
   if(periodicity_y/=0) then
       call cyclize8_y(hhq_rest,nx,ny,1,nnn,nn)
   end if

   call init_computational_domains(lu)

   ! Rayleigh friction initialization
   !$omp parallel do private(m, n, hx2, hy2)
   do n=ny_start,ny_end
       do m=nx_start,nx_end
           if (lu(m,n)*lu(m+1,n)>0.5 .and. lu(m,n)*lu(m,n+1)>0.5) then
               hx2= ( ((hhq_rest(m+1,n)-hhq_rest(m  ,n))/dxt(m  ,n))**2 * dble(lcu(m  ,n))    &
                    +((hhq_rest(m  ,n)-hhq_rest(m-1,n))/dxt(m-1,n))**2 * dble(lcu(m-1,n)) )/dble(lcu(m,n)+lcu(m-1,n))
               hy2= ( ((hhq_rest(m,n+1)-hhq_rest(m,n  ))/dyt(m,n  ))**2 * dble(lcv(m,n  ))    &
                    +((hhq_rest(m,n  )-hhq_rest(m,n-1))/dyt(m,n-1))**2 * dble(lcv(m,n-1)) )/dble(lcv(m,n)+lcv(m,n-1))

               r_diss(m,n)=r_fric*dsqrt(hx2+hy2)
           endif
       enddo
   enddo
   !$omp end parallel do

   call syncborder_real8(r_diss, 1)
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

    integer :: m, n

    ssh = -1
    do m = nx_start, nx_end
        do n = ny_start, ny_end
            ssh(m, n) = rank
        enddo
    enddo

end subroutine

subroutine sw_test2
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables
    use depth_routes
    use ocalg_routes
    use iodata_routes
    implicit none
    integer :: m, n, k, ierr
    real*8 :: hx2, hy2
    real*8 :: u0
    real   :: a
    real*8 :: d2r

    d2r= Pi / 180.0d0

    a = 0.0d0
    u0 = 2d0*Pi*RadEarth / (12.0d0*24*60*60)

    ssh = 0.0d0
    ubrtr = 0.0d0
    vbrtr = 0.0d0

    sshp = 0.0d0
    ubrtrp = 0.0d0
    vbrtrp = 0.0d0

!---------------------Test 2: --------------------------------------------------!
    do n=ny_start, ny_end
      do m=nx_start, nx_end
          if (lcu(m, n) > 0.5) then
              ubrtr(m,n) = u0 * (dcos(d2r*geo_lat_u(m, n))*dcos(d2r*a)            &
                  - dcos(d2r*geo_lon_u(m,n))*dsin(d2r*geo_lat_u(m,n))*dsin(d2r*a))
          endif

          if (lcv(m, n) > 0.5) then
              vbrtr(m,n) = u0 * dsin(d2r*geo_lon_v(m, n))*dsin(d2r*a)
          endif

          if (lu(m ,n) > 0.5) then
              ssh(m,n) = -(1.0d0 / FreeFallAcc)                                 &
                * (RadEarth*EarthAngVel*u0 + 0.5d0*u0*u0)                         &
                  * (( dcos(d2r*geo_lon_t(m,n))*dcos(d2r*geo_lat_t(m,n))*dsin(d2r*a) &
                      + dsin(d2r*geo_lat_t(m,n))*dcos(d2r*a) )**2)
          endif
      enddo
    enddo

    call syncborder_real8(ubrtr, 1)
    call syncborder_real8(vbrtr, 1)
    call syncborder_real8(ssh, 1)
    if(periodicity_x/=0) then
        call cyclize8_x(ssh, nx, ny, 1, mmm, mm)
        call cyclize8_x(ubrtr, nx, ny, 1, mmm, mm)
        call cyclize8_x(vbrtr, nx, ny, 1, mmm, mm)
    end if
    if(periodicity_y/=0) then
        call cyclize8_y(ssh, nx, ny, 1, nnn, nn)
        call cyclize8_y(ubrtr, nx, ny, 1, nnn, nn)
        call cyclize8_y(vbrtr, nx, ny, 1, nnn, nn)
    end if

    ubrtrp = ubrtr
    vbrtrp = vbrtr
    sshp = ssh

!initialize depth for internal mode
    call hh_init(hhq, hhqp, hhqn,    &
                 hhu, hhup, hhun,    &
                 hhv, hhvp, hhvn,    &
                 hhh, hhhp, hhhn,    &
                 ssh, sshp, hhq_rest)

endsubroutine sw_test2

subroutine sw_only_inicond(flag_init, path2ocp)
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables
    use depth_routes
    use ocalg_routes
    use iodata_routes
    implicit none
    integer :: flag_init
    character*(*) path2ocp

    integer :: ierr
    real(4) array4(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

    ubrtr = 0.0d0
    vbrtr = 0.0d0

! Read init sea level
    if (flag_init > 0) then
        if (rank.eq.0) print *, "Read init sea level"
        array4 = 0.0
        call prdstd(path2ocp, 'slf.dat', 1, array4, lu,nx,ny,1, mmm,mm,nnn,nn,1,1,ierr)
        ssh = dble(array4)
    else
        if (rank.eq.0) print *, "Init sea level is zero"
        ssh = 0.0d0
    endif

    call syncborder_real8(ssh, 1)
    if(periodicity_x/=0) then
        call cyclize8_x(ssh, nx, ny, 1, mmm, mm)
    end if
    if(periodicity_y/=0) then
        call cyclize8_y(ssh, nx, ny, 1, nnn, nn)
    end if

    ubrtrp = ubrtr
    vbrtrp = vbrtr
    sshp = ssh

    !initialize depth for internal mode
    call hh_init(hhq, hhqp, hhqn,    &
                 hhu, hhup, hhun,    &
                 hhv, hhvp, hhvn,    &
                 hhh, hhhp, hhhn,    &
                 ssh, sshp, hhq_rest)

endsubroutine sw_only_inicond
