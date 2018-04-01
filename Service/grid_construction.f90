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
