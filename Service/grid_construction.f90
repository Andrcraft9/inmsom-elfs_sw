module gridcon_routes
    implicit none

contains

    !======================================================================
    ! grid consruction module by temperature mask.
    subroutine gridcon()
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
            call set_block(k)
            call set_block_lu(k)
            do n = bnd_y1, bnd_y2
                do m = bnd_x1, bnd_x2
                    if (lbasins(m,n)==0) then
                        lu(m,n)=1.0
                    endif
                end do
            end do
            lu1=1.0
        enddo
        call syncborder_block2D_real4(block_lu)

        !  forming mask for depth grid points
        !  forming luh from lu, which have land neibours in luh.
        ! constructing array luh for relief hh.

        if (rank == 0) then
            write(*,*) 'Construction of H-grid masks: '
            write(*,*) 'LUH (includes boundary) and LUU (does not include boundary)'
        endif

        do k = 1, bcount
            call set_block(k)
            call set_block_lu(k)
            do n = bnd_y1, bnd_y2-1
                do m = bnd_x1, bnd_x2-1
                    if (lu(m,n)+lu(m+1,n)+lu(m,n+1)+lu(m+1,n+1)>0.5)  then
                        luh(m,n)=1.0
                    endif
                enddo
            enddo
            do n = bnd_y1, bnd_y2-1
                do m = bnd_x1, bnd_x2-1
                    if(lu(m,n)*lu(m+1,n)*lu(m,n+1)*lu(m+1,n+1)>0.5)  then
                        luu(m,n)=1.0
                    endif
                enddo
            enddo
        enddo
        call syncborder_block2D_real4(block_luh)
        call syncborder_block2D_real4(block_luu)

        if (rank == 0) then
            write(*,*) 'Construction of U- and V-grid masks: '
            write(*,*) 'LCU and LCV (do not include boundary) and LLU and LLV (include boundary)'
        endif

        do k = 1, bcount
            call set_block(k)
            call set_block_lu(k)
            do n = bnd_y1, bnd_y2-1
                do m = bnd_x1, bnd_x2-1
                    if (lu(m,n)+lu(m+1,n)>0.5) then
                        llu(m,n) = 1.0
                    endif

                    if (lu(m,n)+lu(m,n+1)>0.5) then
                        llv(m,n) = 1.0
                    endif

                    if (lu(m,n)*lu(m+1,n)>0.5) then
                        lcu(m,n) = 1.0
                    endif

                    if (lu(m,n)*lu(m,n+1)>0.5) then
                        lcv(m,n) = 1.0
                    endif
                enddo
            enddo
        enddo
        call syncborder_block2D_real4(block_llu)
        call syncborder_block2D_real4(block_llv)
        call syncborder_block2D_real4(block_lcu)
        call syncborder_block2D_real4(block_lcv)

        return
    endsubroutine gridcon

    !======================================================================
    ! Construction u and v masks for output area from temperature mask
    subroutine gridcon_output(ftemask)
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
        integer m, n, k, ierr
        integer, allocatable:: lbasins_output(:,:)

        allocate(lbasins_output(nx, ny))
        lbasins_output = 0

        write(frmt,1000) nx
1000    format('(',i9,'i1)')

        ! reading mask from:
        if (rank .eq. 0) then
            open (11,file=ftemask,status='old',recl=nx*lrecl)
            read (11,  '(a)') comment(1:min(80,nx))
            if (rank .eq. 0) write(*,'(1x,a)') comment
            do n=ny,1,-1
                read(11,frmt,end=99) (lbasins_output(m,n),m=1,nx)
            enddo
            close(11)
        endif
        call mpi_bcast(lbasins_output, nx*ny, mpi_integer, 0, cart_comm, ierr)

        ! conversion integer diogin mask to real model mask
        do k = 1, bcount
            call set_block(k)
            call set_block_lu_output(k)
            do n = bnd_y1, bnd_y2
                do m = bnd_x1, bnd_x2
                    if(lbasins_output(m,n)==0) then
                        lu_output(m,n)=1.0
                    endif
                enddo
            enddo
        enddo
        call syncborder_block2D_real4(block_lu_output)

        if (rank .eq. 0) then
            write(*,*) 'Construction of U- and V-grid masks for output area: '
            write(*,*) 'LCU and LCV (do not include boundary) and LLU and LLV (include boundary)'
        endif

        do k = 1, bcount
            call set_block(k)
            call set_block_lu_output(k)
            do n = bnd_y1, bnd_y2-1
                do m = bnd_x1, bnd_x2-1
                    if (lu_output(m,n)+lu_output(m+1,n)>0.5) then
                        llu_output(m,n) = 1.0
                    endif

                    if (lu_output(m,n)+lu_output(m,n+1)>0.5) then
                        llv_output(m,n) = 1.0
                    endif

                    if (lu_output(m,n)*lu_output(m+1,n)>0.5) then
                        lcu_output(m,n) = 1.0
                    endif

                    if (lu_output(m,n)*lu_output(m,n+1)>0.5) then
                        lcv_output(m,n) = 1.0
                    endif
                enddo
            enddo
        enddo
        call syncborder_block2D_real4(block_llu_output)
        call syncborder_block2D_real4(block_llv_output)
        call syncborder_block2D_real4(block_lcu_output)
        call syncborder_block2D_real4(block_lcv_output)

        deallocate(lbasins_output)
        return

99      write(*,*)'  error in reading file ',ftemask(1:len_trim(ftemask))
        call mpi_abort(cart_comm, 1, ierr)
        stop 1
    endsubroutine gridcon_output

endmodule gridcon_routes
