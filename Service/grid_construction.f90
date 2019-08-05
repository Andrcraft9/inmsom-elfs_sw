module gridcon_routes
implicit none

contains

!======================================================================
! grid consruction module by temperature mask.
subroutine gridcon()
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use rec_length
    implicit none
    ! subroutin for construction pass boundary, velosity and bottom masks
    ! using temperature mask in diogin standart
    !--------------------------------------------------------------------
    ! temporary integer indexes
    integer k, m, n, ierr

    ! conversion integer diogin mask to real model mask
    do k = 1, bcount
        call set_block_boundary(k)

        do n = bnd_y1, bnd_y2
            do m = bnd_x1, bnd_x2
                if(lbasins(m,n) == 0) then
                    lu(k)%vals(m,n) = 1.0d0
                endif
            end do
        end do
        lu1(k)%vals = 1.0d0
    enddo

    !  forming mask for depth grid points
    !  forming luh from lu, which have land neibours in luh.
    !  constructing array luh for relief hh.
    if (rank .eq. 0) then
        write(*,*) 'Construction of H-grid masks: '
        write(*,*) 'LUH (includes boundary) and LUU (does not include boundary)'
    endif

    do k = 1, bcount
        call set_block_boundary(k)

        !do n=ny_start-1,ny_end
        !do m=nx_start-1,nx_end
        do n = bnd_y1, bnd_y2-1
            do m = bnd_x1, bnd_x2-1
                if (lu(k)%vals(m,n) + lu(k)%vals(m+1,n) + lu(k)%vals(m,n+1) + lu(k)%vals(m+1,n+1) > 0.5d0) then
                    luh(k)%vals(m,n) = 1.0d0
                endif
            enddo
        enddo

        !do n=ny_start-1,ny_end
        !do m=nx_start-1,nx_end
        do n = bnd_y1, bnd_y2-1
            do m = bnd_x1, bnd_x2-1
                if (lu(k)%vals(m,n)*lu(k)%vals(m+1,n)*lu(k)%vals(m,n+1)*lu(k)%vals(m+1,n+1) > 0.5d0) then
                    luu(k)%vals(m,n) = 1.0d0
                endif
            enddo
        enddo
    enddo

    if (rank .eq. 0) then
        write(*,*) 'Construction of U- and V-grid masks: '
        write(*,*) 'LCU and LCV (do not include boundary) and LLU and LLV (include boundary)'
    endif

    do k = 1, bcount
        call set_block_boundary(k)

        !do n=ny_start-1,ny_end
        !do m=nx_start-1,nx_end
        do n=bnd_y1, bnd_y2-1
            do m=bnd_x1, bnd_x2-1
                if (lu(k)%vals(m,n) + lu(k)%vals(m+1,n) > 0.5d0) then
                    llu(k)%vals(m,n) = 1.0d0
                endif

                if (lu(k)%vals(m,n) + lu(k)%vals(m,n+1) > 0.5d0) then
                    llv(k)%vals(m,n) = 1.0d0
                endif

                if (lu(k)%vals(m,n)*lu(k)%vals(m+1,n) > 0.5d0) then
                    lcu(k)%vals(m,n) = 1.0d0
                endif

                if (lu(k)%vals(m,n)*lu(k)%vals(m,n+1) > 0.5d0) then
                    lcv(k)%vals(m,n) = 1.0d0
                endif
            enddo
        enddo
    enddo

    return
endsubroutine gridcon

!======================================================================
! Local grid consruction module by temperature mask.
subroutine gridcon_local(ftemask)
  use main_basin_pars
  use mpi_parallel_tools
  use basin_grid
  use rec_length
  use ocalg_routes
  implicit none
  ! subroutin for construction pass boundary, velosity and bottom masks
  ! using temperature mask in diogin standart
  !--------------------------------------------------------------------
  character*(*) ftemask
  character frmt*16,comment*80
  integer m, n, ierr
  integer, allocatable:: lbasins_local(:,:)

  allocate(lbasins_local(nx, ny))
  lbasins_local = 0
  
  write(frmt,1000) nx
  1000  format('(',i9,'i1)')
  
  ! reading mask from:
  if (rank .eq. 0) then
      open (11,file=ftemask,status='old',recl=nx*lrecl)
       read (11,  '(a)') comment(1:min(80,nx))
       if (rank .eq. 0) write(*,'(1x,a)') comment
       do n=ny,1,-1
        read(11,frmt,end=99) (lbasins_local(m,n),m=1,nx)
       enddo
      close(11)
  endif
  call mpi_bcast(lbasins_local, nx*ny, mpi_integer, 0, cart_comm, ierr)
  
  ! conversion integer diogin mask to real model mask
  do n=bnd_y1, bnd_y2
     do m=bnd_x1, bnd_x2
      if(lbasins_local(m,n)==0) then
       lu_local(m,n)=1.0
      endif
     end do
  end do
  
  
  if(periodicity_x/=0) then
    call cyclize_x(lu_local,nx,ny,1,mmm,mm)
  endif
  
  if(periodicity_y/=0) then
    call cyclize_y(lu_local,nx,ny,1,nnn,nn)
  endif
  
  if (rank .eq. 0) then
      write(*,*) 'Construction of U- and V-grid masks for LOCAL area: '
      write(*,*) 'LCU and LCV (do not include boundary) and LLU and LLV (include boundary)'
  endif
  
  !do n=ny_start-1,ny_end
  !do m=nx_start-1,nx_end
  do n=bnd_y1, bnd_y2-1
     do m=bnd_x1, bnd_x2-1
      if(lu_local(m,n) + lu_local(m+1,n)>0.5)  then
       llu_local(m,n)=1.0
      endif
  
      if(lu_local(m,n) + lu_local(m,n+1)>0.5)  then
       llv_local(m,n)=1.0
      endif
  
      if(lu_local(m,n) * lu_local(m+1,n)>0.5)  then
       lcu_local(m,n)=1.0
      endif
  
      if(lu_local(m,n) * lu_local(m,n+1)>0.5)  then
       lcv_local(m,n)=1.0
      endif
  
     enddo
  enddo
  
  if (periodicity_x/=0) then
   if (rank .eq. 0) write(*,*)'  set periodicity to u-grid mask(lcu,llu).'
   call cyclize_x(lcu_local,nx,ny,1,mmm,mm)
   call cyclize_x(llu_local,nx,ny,1,mmm,mm)
   if (rank .eq. 0) write(*,*)'  set periodicity to v-grid mask(lcv,llv).'
   call cyclize_x(lcv_local,nx,ny,1,mmm,mm)
   call cyclize_x(llv_local,nx,ny,1,mmm,mm)
  endif
  
  if (periodicity_y/=0) then
   call cyclize_y(lcu_local,nx,ny,1,nnn,nn)
   call cyclize_y(llu_local,nx,ny,1,nnn,nn)
   call cyclize_y(lcv_local,nx,ny,1,nnn,nn)
   call cyclize_y(llv_local,nx,ny,1,nnn,nn)
  endif
  
  deallocate(lbasins_local)

  return

  99    write(*,*)'  error in reading file ',ftemask(1:len_trim(ftemask))
  call mpi_abort(cart_comm, 1, ierr)
  stop 1
  endsubroutine gridcon_local

endmodule gridcon_routes
