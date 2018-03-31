!-----------module for definition of array dimensions and boundaries------------
module mpi_parallel_tools
    use mpi
    implicit none

    !include 'mpif.h'
    include "omp_lib.h"

    integer :: recommended_tot_blocks
    integer :: parallel_dbg

    integer :: nx_start, nx_end, ny_start, ny_end
    integer :: bnd_x1, bnd_x2, bnd_y1, bnd_y2

    type :: block2D
        real*8, dimension(:, :), pointer :: vals
    end type

    type :: block1D
        real*8, dimension(:), pointer :: vals
    end type

    integer, allocatable :: lbasins(:, :)

    ! There are two types of block numeration:
    !   Local:
    !       for each proc k = 1, ..., bcount. In custom order!
    !   Global:
    !       Cart coords of blocks (1, 1), ..., (bnx, bny)

    ! Cart grid of blocks
    integer :: bnx, bny
    ! Map: block coords to proc
    ! If bglob_proc(m, n) == -1 => (m, n) block is land-block!
    integer, allocatable :: bglob_proc(:, :)

    integer :: bcount
    ! Total blocks for computational
    integer :: total_blocks
    ! Map: local block number to block coords
    integer, allocatable :: bindx(:, :)
    ! significant point area in blocks, local block numeration
    integer, allocatable :: bnx_start(:), bnx_end(:),    &
                            bny_start(:), bny_end(:)
    ! array boundary in blocks, local block numeration
    integer, allocatable :: bbnd_x1(:), bbnd_x2(:),  &
                            bbnd_y1(:), bbnd_y2(:)

    ! Proc grid information
    integer :: rank, procs
    integer :: cart_comm
    integer, dimension(2) :: p_size, p_coord
    logical, dimension(2) :: p_period

    ! MPI buffers
    integer, allocatable :: reqsts(:), statuses(:, :)
    type(block1D), dimension(:), pointer :: sync_buf8_send_nyp, sync_buf8_recv_nyp
    type(block1D), dimension(:), pointer :: sync_buf8_send_nxp, sync_buf8_recv_nxp
    type(block1D), dimension(:), pointer :: sync_buf8_send_nym, sync_buf8_recv_nym
    type(block1D), dimension(:), pointer :: sync_buf8_send_nxm, sync_buf8_recv_nxm
    real*8 :: sync_edge_buf8_recv_nxp_nyp
    real*8 :: sync_edge_buf8_recv_nxp_nym
    real*8 :: sync_edge_buf8_recv_nxm_nyp
    real*8 :: sync_edge_buf8_recv_nxm_nym

    ! Timers
    real*8 :: time_barotrop
    real*8 :: time_model_step, time_output
    real*8 :: time_sync

contains

    subroutine parallel_check_err(err)
        implicit none

        integer :: err, totalerr
        integer :: ierr

        call mpi_allreduce(err, totalerr, 1, mpi_integer,  &
                           mpi_sum, cart_comm, ierr)
        if (totalerr >= 1) then
            if (rank .eq. 0) print *, 'Error! ABORT PROGRAM!'
            call parallel_finalize
            stop
        endif
    end subroutine

    subroutine parallel_init()
        implicit none

        integer :: count_threads, num_thread
        integer :: ierr, rank_cart

        call mpi_init(ierr)
        call mpi_comm_rank(mpi_comm_world, rank, ierr)
        if (rank .eq. 0) then
            print *, 'Read parallel.par...'

            open(90, file='parallel.par', status='old')
            read(90, *) recommended_tot_blocks
            print *, 'recommended_tot_blocks=', recommended_tot_blocks
            read(90, *) parallel_dbg
            print *, 'parallel_dbg=', parallel_dbg

            close(90)
        endif
        call mpi_bcast(recommended_tot_blocks, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(parallel_dbg, 1, mpi_integer, 0, mpi_comm_world, ierr)

        if (rank .eq. 0) then
            if (.not. MPI_SUBARRAYS_SUPPORTED) then
                print *, 'MPI_SUBARRAYS_SUPPORTED = FALSE'
            endif
        endif

        p_period = (/.true., .true./)
        p_size = (/0,0/)
        ierr = 0
        call mpi_comm_size(mpi_comm_world, procs, ierr)
        call mpi_dims_create(procs, 2, p_size, ierr)
        call mpi_cart_create(mpi_comm_world, 2, p_size, p_period, .false., cart_comm, ierr)
        call mpi_cart_coords(cart_comm, rank, 2, p_coord, ierr)

        call mpi_comm_rank(cart_comm, rank_cart, ierr)

        print *, 'rank_world, rank_cart and cord:', rank, rank_cart, p_coord
        call mpi_barrier(cart_comm, ierr)

        !$omp parallel
        count_threads = omp_get_num_threads()
        num_thread = omp_get_thread_num()
        if (num_thread .eq. 0) print *, "OMP Threads: ", count_threads
        !$omp end parallel
        call mpi_barrier(cart_comm, ierr)
    end subroutine

    subroutine parallel_finalize
        implicit none

        integer :: k, ierr

        if (allocated(bnx_start)) deallocate(bnx_start)
        if (allocated(bnx_end)) deallocate(bnx_end)
        if (allocated(bny_start)) deallocate(bny_start)
        if (allocated(bny_end)) deallocate(bny_end)

        if (allocated(bbnd_x1)) deallocate(bbnd_x1)
        if (allocated(bbnd_x2)) deallocate(bbnd_x2)
        if (allocated(bbnd_y1)) deallocate(bbnd_y1)
        if (allocated(bbnd_y2)) deallocate(bbnd_y2)

        if (allocated(bindx)) deallocate(bindx)
        if (allocated(bglob_proc)) deallocate(bglob_proc)
        if (allocated(lbasins)) deallocate(lbasins)

        if (allocated(reqsts)) deallocate(reqsts)
        if (allocated(statuses)) deallocate(statuses)

        !if (allocated(sync_buf8_send_nyp)) then
        do k = 1, bcount
            deallocate(sync_buf8_send_nyp(k)%vals, sync_buf8_recv_nyp(k)%vals)
            deallocate(sync_buf8_send_nxp(k)%vals, sync_buf8_recv_nxp(k)%vals)
            deallocate(sync_buf8_send_nym(k)%vals, sync_buf8_recv_nym(k)%vals)
            deallocate(sync_buf8_send_nxm(k)%vals, sync_buf8_recv_nxm(k)%vals)
        enddo
        deallocate(sync_buf8_send_nyp, sync_buf8_recv_nyp)
        deallocate(sync_buf8_send_nxp, sync_buf8_recv_nxp)
        deallocate(sync_buf8_send_nym, sync_buf8_recv_nym)
        deallocate(sync_buf8_send_nxm, sync_buf8_recv_nxm)
        !endif
        !if (allocated(sync_buf8_send_nyp)) deallocate(sync_buf8_send_nyp)
        !if (allocated(sync_buf8_recv_nyp)) deallocate(sync_buf8_recv_nyp)

        call mpi_finalize(ierr)
    end subroutine

    subroutine parallel_read_mask(ftemask)
        use main_basin_pars
        use rec_length
        implicit none

        character*(*) ftemask
        character frmt*16,comment*80
        integer :: m, n, ierr

        write(frmt,1000) nx
1000    format('(',i9,'i1)')

        allocate(lbasins(nx, ny))
        ! reading mask from:
        if (rank .eq. 0) then
            open(11, file=ftemask, status='old', recl=nx*lrecl)
            read(11, '(a)') comment(1:min(80,nx))
            if (rank .eq. 0) write(*,'(1x,a)') comment
            do n = ny, 1, -1
                read(11, frmt, end=99) (lbasins(m,n),m=1,nx)
            enddo
            close(11)
        endif
        call mpi_bcast(lbasins, nx*ny, mpi_integer, 0, cart_comm, ierr)

        return
99      write(*,*) 'error in reading file', ftemask(1:len_trim(ftemask))
        stop 1
    end subroutine

    subroutine parallel_uniform_distribution()
        use main_basin_pars
        implicit none

        integer, allocatable :: bweight(:)
        real*8 :: weight

        integer :: loc_bnx, loc_bny

        integer :: m, n, i, j, k, glob_k, locn
        integer :: xblock_start, yblock_start
        integer :: land_blocks
        integer :: ierr
        ! Buffers for MPI subroutines
        integer, allocatable :: buf_int(:, :)

        ! Set Cart grid of blocks
        bnx = floor(sqrt(recommended_tot_blocks * real(nx-4) / real(ny-4)))
        bny = floor(sqrt(recommended_tot_blocks * real(ny-4) / real(nx-4)))
        bnx = bnx - mod(bnx, p_size(1))
        bny = bny - mod(bny, p_size(2))
        if (parallel_dbg >= 1 .and. rank == 0) print *, 'bnx, bny and Total blocks:', bnx, bny, bnx*bny
        call mpi_barrier(cart_comm, ierr)

        ! Unifrom distribute blocks to procs
        ierr = 0
        if (mod(bnx, p_size(1)) /= 0) then
            print *, 'Uniform blocks distribution is impossible',     &
                        'Please check procs and total blocks numbers!'
            ierr = 1
        endif
        if (mod(bny, p_size(2)) /= 0) then
            print *, 'Uniform blocks distribution is impossible',     &
                        'Please check procs and total blocks numbers!'
            ierr = 1
        endif
        call parallel_check_err(ierr)
        loc_bnx = bnx / p_size(1)
        loc_bny = bny / p_size(2)
        bcount = loc_bnx*loc_bny
        if (parallel_dbg >= 1) print *, rank, 'loc_bnx, loc_bny and Blocks per proc: ', loc_bnx, loc_bny, bcount
        call mpi_barrier(cart_comm, ierr)

        xblock_start = 1 + p_coord(1)*loc_bnx
        yblock_start = 1 + p_coord(2)*loc_bny
        if (parallel_dbg >= 2) print *, rank, 'xb_start, yb_start', xblock_start, yblock_start
        call mpi_barrier(cart_comm, ierr)

        allocate(bglob_proc(bnx, bny))
        allocate(bindx(bcount, 2))
        allocate(bnx_start(bcount), bnx_end(bcount),     &
                 bny_start(bcount), bny_end(bcount))
        allocate(bbnd_x1(bcount), bbnd_x2(bcount),       &
                 bbnd_y1(bcount), bbnd_y2(bcount))
        allocate(bweight(bcount))

        bindx = 0
        bnx_start = 0; bnx_end = 0; bny_start = 0; bny_end = 0
        bbnd_x1 = 0; bbnd_x2 = 0; bbnd_y1 = 0; bbnd_y2 = 0
        bweight = 0
        bglob_proc = -1

        k = 1
        land_blocks = 0
        weight = 0.0d0
        do m = xblock_start, xblock_start + loc_bnx - 1
            do n = yblock_start, yblock_start + loc_bny - 1

                ! Map block coords to procs
                bglob_proc(m, n) = rank
                ! Map local block numeration to block coords
                bindx(k, 1) = m
                bindx(k, 2) = n

                locn = floor(real(nx - 4)/real(bnx))
                bnx_start(k) = locn*(m-1) + 1 + 2
                if (m .eq. bnx) then
                    locn = (nx - 2) - bnx_start(k) + 1
                endif
                bnx_end(k) = bnx_start(k) + locn - 1
                bnx_start(k) = bnx_start(k)
                ! border area
                bbnd_x1(k) = bnx_start(k) - 2
                bbnd_x2(k) = bnx_end(k) + 2

                locn = floor(real(ny - 4)/real(bny))
                bny_start(k) = locn*(n-1) + 1 + 2
                if (n .eq. bny) then
                    locn = (ny - 2) - bny_start(k) + 1
                endif
                bny_end(k) = bny_start(k) + locn - 1
                bny_start(k) = bny_start(k)
                ! border area
                bbnd_y1(k) = bny_start(k) - 2
                bbnd_y2(k) = bny_end(k) + 2

                ! Compute load-balance
                do i = bnx_start(k), bnx_end(k)
                    do j = bny_start(k), bny_end(k)
                        bweight(k) = bweight(k) + (1.0d0 - real(lbasins(i, j)))
                    enddo
                enddo
                weight = weight + real(bweight(k)) !/ real((bnx_end(k) - bnx_start(k) + 1)*(bny_end(k) - bny_start(k) + 1))
                ! Ignore only-land blocks
                if (bweight(k) == 0.0) then
                    land_blocks = land_blocks + 1
                    bglob_proc(m, n) = -1
                else
                    k = k + 1
                endif
            enddo
        enddo

        if (parallel_dbg >= 1) print *, rank, 'Total blocks with land:', land_blocks, 'Total Load-Balancing:', weight / ((nx-4)*(ny-4))
        call mpi_barrier(cart_comm, ierr)
        bcount = bcount - land_blocks
        if (parallel_dbg >= 1) print *, rank, 'Updated block count', bcount
        call mpi_barrier(cart_comm, ierr)
        ierr = 0
        if (bcount <= 0) then
            print *, rank, 'All blocks with land blocks!'
            ierr = 1
        endif
        call parallel_check_err(ierr)
        call mpi_allreduce(bcount, total_blocks, 1, mpi_integer,      &
                           mpi_sum, cart_comm, ierr)
        if (parallel_dbg >= 1 .and. rank == 0) print *, 'Updated Total block count', total_blocks
        call mpi_barrier(cart_comm, ierr)

        ! Sync bglob_proc array
        if (parallel_dbg >= 2) call parallel_int_output(bglob_proc, 1, bnx, 1, bny, 'bglob_proc before')
        allocate(buf_int(bnx, bny))
        buf_int = bglob_proc + 1
        call mpi_allreduce(buf_int, bglob_proc, bnx*bny, mpi_integer, mpi_sum, cart_comm, ierr)
        bglob_proc = bglob_proc - 1
        if (parallel_dbg >= 2) call parallel_int_output(bglob_proc, 1, bnx, 1, bny, 'bglob_proc after')

        !allocate(buf_int(bcount), recv_buf_int(procs), displs(procs))
        !call mpi_allgather(bcount, 1, mpi_integer, recv_buf_int, 1, mpi_integer, cart_comm, ierr)
        !displs(1) = 1
        !do k = 2, procs
        !    displs(k) = displs(k-1) + recv_buf_int(k)
        !enddo
        !buf_int = bglob_proc(1:bcount)
        !call mpi_allgatherv(buf_int, bcount, mpi_integer, bglob_proc, recv_buf_int, displs, mpi_integer, cart_comm, ierr)

        call allocate_mpi_buffers()

        deallocate(bweight)
        deallocate(buf_int)

    end subroutine

    subroutine allocate_mpi_buffers
        implicit none

        integer :: k, nx_dir_size, ny_dir_size

        allocate(sync_buf8_send_nyp(bcount), sync_buf8_recv_nyp(bcount))
        allocate(sync_buf8_send_nxp(bcount), sync_buf8_recv_nxp(bcount))
        allocate(sync_buf8_send_nym(bcount), sync_buf8_recv_nym(bcount))
        allocate(sync_buf8_send_nxm(bcount), sync_buf8_recv_nxm(bcount))
        do k = 1, bcount
            ny_dir_size = bnx_end(k) - bnx_start(k) + 1
            nx_dir_size = bny_end(k) - bny_start(k) + 1
            allocate(sync_buf8_send_nyp(k)%vals(ny_dir_size), sync_buf8_recv_nyp(k)%vals(ny_dir_size))
            allocate(sync_buf8_send_nxp(k)%vals(nx_dir_size), sync_buf8_recv_nxp(k)%vals(nx_dir_size))
            allocate(sync_buf8_send_nym(k)%vals(ny_dir_size), sync_buf8_recv_nym(k)%vals(ny_dir_size))
            allocate(sync_buf8_send_nxm(k)%vals(nx_dir_size), sync_buf8_recv_nxm(k)%vals(nx_dir_size))
        enddo

        allocate(reqsts(8*bcount), statuses(MPI_STATUS_SIZE, 8*bcount))
    end subroutine

    subroutine parallel_int_output(arr, x1, x2, y1, y2, msg)
        implicit none

        integer :: x1, x2, y1, y2
        integer :: arr(x1:x2, y1:y2)
        character*(*) msg
        integer :: k, m, n, ierr

        if (rank .eq. 0) print *, msg
        call mpi_barrier(cart_comm, ierr)

        do k = 0, procs-1
            if (rank .eq. k) then
                print *, rank
                do m = x1, x2
                    do n = y1, y2
                        print *, rank, 'm, n, a(m, n)', m, n, arr(m, n)
                    enddo
                enddo
            endif
            call mpi_barrier(cart_comm, ierr)
        enddo

    end subroutine

    subroutine allocate_block2D(blks, val)
        implicit none

        type(block2D), dimension(:), pointer :: blks
        real*8 :: val
        integer :: k

        allocate(blks(bcount))
        do k = 1, bcount
            allocate(blks(k)%vals(bbnd_x1(k):bbnd_x2(k), bbnd_y1(k):bbnd_y2(k)))
            blks(k)%vals = val
        enddo
    end subroutine

    subroutine deallocate_blocks2D(blks)
        implicit none

        type(block2D), dimension(:), pointer :: blks
        integer :: k

        do k = 1, bcount
            deallocate(blks(k)%vals)
        enddo
        deallocate(blks)
    end subroutine

    subroutine init_times
        implicit none
        time_model_step = 0.0d0
        time_barotrop = 0.0d0
        time_sync = 0.0d0
        time_output = 0.0d0
        return
    end subroutine

    subroutine print_times
        implicit none
        if (rank .eq. 0) then
            print *, "Time model step: ", time_model_step
            print *, "Time barotropic: ", time_barotrop
            print *, "Time sync: ", time_sync
            print *, "Time output: ", time_output
        endif
        return
    end subroutine

    subroutine start_timer(time)
        implicit none

        real*8, intent(inout) :: time

        time = mpi_wtime()
        return
    end subroutine

    subroutine end_timer(time)
        implicit none

        real*8, intent(inout) :: time
        real*8 :: outtime
        integer :: ierr

        time = mpi_wtime() - time
        call mpi_allreduce(time, outtime, 1, mpi_real8,      &
                           mpi_max, cart_comm, ierr)
        time = outtime
        return
    end subroutine

    integer function get_local_block_number(m, n)
        implicit none

        integer :: m, n
        integer :: k

        get_local_block_number = -1
        do k = 1, bcount
            if (m == bindx(k, 1) .and. n == bindx(k, 2)) then
                get_local_block_number = k
                exit
            endif
        enddo

        return
    end function


!    subroutine irecv_real8(k, src_block, blks, tag, dx1, dx2, dy1, dy2, reqst, stat)
!        implicit none
!        integer :: p_src, k, tag
!        integer, dimension(2) :: src_block
!        type(block2D), dimension(:), pointer :: blks
!        integer :: dx1, dx2, dy1, dy2
!        integer :: reqst, stat
!        integer :: ierr, flg_recv
!        integer, dimension(2) :: bgrid
                !call mpi_irecv(blks(k)%vals(dx1:dx2, dy1:dy2), (dx2 - dx1 + 1)*(dy2 - dy1 + 1), &
!                                             mpi_real8, p_src, tag, cart_comm, reqst, ierr)

    subroutine irecv_real8(k, src_block, src_proc, buff, buff_size, tag, reqst)
        implicit none
        integer, dimension(2) :: src_block
        integer :: k, src_proc, tag
        integer :: buff_size
        real*8 :: buff(:)
        integer :: reqst
        integer :: ierr

        if (parallel_dbg >= 2) print *, rank, 'IRECV. block: ', bindx(k, 1), bindx(k, 2), 'src_block:', src_block(1), src_block(2), 'src_p:', src_proc, 'tag', tag
        !call mpi_irecv(blks(k)%vals(dx1:dx2, dy1:dy2), (dx2 - dx1 + 1)*(dy2 - dy1 + 1), &
        !               mpi_real8, p_src, tag, cart_comm, reqst, ierr)
        call mpi_irecv(buff, buff_size, mpi_real8, src_proc, tag, cart_comm, reqst, ierr)

    end subroutine

    subroutine isend_real8(k, dist_block, dist_proc, buff, buff_size, tag, reqst)
        implicit none

        integer, dimension(2) :: dist_block
        integer :: k, dist_proc, tag
        integer :: buff_size
        real*8 :: buff(:)
        integer :: reqst
        integer :: ierr

        if (parallel_dbg >= 2) print *, rank, 'ISEND. block: ', bindx(k, 1), bindx(k, 2), 'dst_block:', dist_block(1), dist_block(2), 'dst_p:', dist_proc, 'tag', tag
        !call mpi_isend(blks(k)%vals(sx1:sx2, sy1:sy2), (sx2 - sx1 + 1)*(sy2 - sy1 + 1), &
        !               mpi_real8, p_dist, tag, cart_comm, reqst, ierr)
        call mpi_isend(buff, buff_size, mpi_real8, dist_proc, tag, cart_comm, reqst, ierr)

    end subroutine

    subroutine check_block_status(b_coords, p)
        implicit none

        integer, dimension(2), intent(in) :: b_coords
        integer, intent(out) :: p
        integer, dimension(2) :: bgrid

        bgrid(1) = bnx; bgrid(2) = bny
        if (check_cart_coord(b_coords - 1, bgrid) /= 1) then
            ! Out of range, block does not exist
            p = -2
            return
        endif

        p = bglob_proc(b_coords(1), b_coords(2))
        return

    end subroutine

    subroutine syncborder_block2D(blks)
        implicit none

        type(block2D), dimension(:), pointer :: blks
        integer :: k, lock, bm, bn, reqst
        integer, dimension(2) :: src_block, dist_block
        integer :: src_p, dist_p
        integer :: send_size, recv_size
        integer :: ierr
        integer :: tag, icount

        icount = 1

        ! Non-blocking Recv calls
        do k = 1, bcount
            bm = bindx(k, 1)
            bn = bindx(k, 2)

            ! irecv in ny+
            src_block(1) = bm; src_block(2) = bn - 1
            tag = src_block(2) + (src_block(1) - 1)*bnx
            call check_block_status(src_block, src_p)
            if (src_p >= 0) then
                if (src_p /= rank) then
                    recv_size = bnx_end(k) - bnx_start(k) + 1
                    call irecv_real8(k, src_block, src_p, sync_buf8_recv_nyp(k)%vals, recv_size, tag, reqst)
                    !call irecv_real8(k, src_block, src_p, blks(k)%vals(bnx_start(k):bnx_end(k), bbnd_y1(k) + 1), recv_size, bn, reqst)
                    reqsts(icount) = reqst
                    icount = icount + 1
                endif
            endif
            ! irecv in nx+
            src_block(1) = bm - 1; src_block(2) = bn
            tag = src_block(2) + (src_block(1) - 1)*bnx
            call check_block_status(src_block, src_p)
            if (src_p >= 0) then
                if (src_p /= rank) then
                    recv_size = bny_end(k) - bny_start(k) + 1
                    call irecv_real8(k, src_block, src_p, sync_buf8_recv_nxp(k)%vals, recv_size, tag, reqst)
                    !call irecv_real8(k, src_block, src_p, blks(k)%vals(bbnd_x1(k) + 1, bny_start(k):bny_end(k)), recv_size, bn, reqst)
                    reqsts(icount) = reqst
                    icount = icount + 1
                endif
            endif
            ! irecv in ny-
            src_block(1) = bm; src_block(2) = bn + 1
            tag = src_block(2) + (src_block(1) - 1)*bnx
            call check_block_status(src_block, src_p)
            if (src_p >= 0) then
                if (src_p /= rank) then
                    recv_size = bnx_end(k) - bnx_start(k) + 1
                    call irecv_real8(k, src_block, src_p, sync_buf8_recv_nym(k)%vals, recv_size, tag, reqst)
                    !call irecv_real8(k, src_block, src_p, blks(k)%vals(bnx_start(k):bnx_end(k), bbnd_y2(k) - 1), recv_size, bn, reqst)
                    reqsts(icount) = reqst
                    icount = icount + 1
                endif
            endif
            ! irecv in nx-
            src_block(1) = bm + 1; src_block(2) = bn
            tag = src_block(2) + (src_block(1) - 1)*bnx
            call check_block_status(src_block, src_p)
            if (src_p >= 0) then
                if (src_p /= rank) then
                    recv_size = bny_end(k) - bny_start(k) + 1
                    call irecv_real8(k, src_block, src_p, sync_buf8_recv_nxm(k)%vals, recv_size, tag, reqst)
                    !call irecv_real8(k, src_block, src_p, blks(k)%vals(bbnd_x2(k) - 1, bny_start(k):bny_end(k)), recv_size, bn, reqst)
                    reqsts(icount) = reqst
                    icount = icount + 1
                endif
            endif

            !----------------------------Edge points----------------------------!
            ! Edge irecv in nx+ ny+
            src_block(1) = bm - 1; src_block(2) = bn - 1
            tag = 10*(src_block(2) + (src_block(1) - 1)*bnx)
            call check_block_status(src_block, src_p)
            if (src_p >= 0) then
                if (src_p /= rank) then
                    call mpi_irecv(sync_edge_buf8_recv_nxp_nyp, 1, mpi_real8, src_p, tag, cart_comm, reqst, ierr)
                    reqsts(icount) = reqst
                    icount = icount + 1
                endif
            endif
            ! Edge irecv in nx+ ny-
            src_block(1) = bm - 1; src_block(2) = bn + 1
            tag = 10*(src_block(2) + (src_block(1) - 1)*bnx)
            call check_block_status(src_block, src_p)
            if (src_p >= 0) then
                if (src_p /= rank) then
                    call mpi_irecv(sync_edge_buf8_recv_nxp_nym, 1, mpi_real8, src_p, tag, cart_comm, reqst, ierr)
                    reqsts(icount) = reqst
                    icount = icount + 1
                endif
            endif
            ! Edge irecv in nx- ny+
            src_block(1) = bm + 1; src_block(2) = bn - 1
            tag = 10*(src_block(2) + (src_block(1) - 1)*bnx)
            call check_block_status(src_block, src_p)
            if (src_p >= 0) then
                if (src_p /= rank) then
                    call mpi_irecv(sync_edge_buf8_recv_nxm_nyp, 1, mpi_real8, src_p, tag, cart_comm, reqst, ierr)
                    reqsts(icount) = reqst
                    icount = icount + 1
                endif
            endif
            ! Edge irecv in nx- ny-
            src_block(1) = bm + 1; src_block(2) = bn + 1
            tag = 10*(src_block(2) + (src_block(1) - 1)*bnx)
            call check_block_status(src_block, src_p)
            if (src_p >= 0) then
                if (src_p /= rank) then
                    call mpi_irecv(sync_edge_buf8_recv_nxm_nym, 1, mpi_real8, src_p, tag, cart_comm, reqst, ierr)
                    reqsts(icount) = reqst
                    icount = icount + 1
                endif
            endif

        enddo
        if (parallel_dbg >= 2) print *, rank, 'icount recv:', icount-1

        ! Non-blocking Send calls
        do k = 1, bcount
            bm = bindx(k, 1)
            bn = bindx(k, 2)

            tag = bn + (bm - 1)*bnx
            ! isend in ny+
            dist_block(1) = bm; dist_block(2) = bn + 1
            call check_block_status(dist_block, dist_p)
            if (dist_p >= 0) then
                if (dist_p /= rank) then
                    send_size = bnx_end(k) - bnx_start(k) + 1
                    sync_buf8_send_nyp(k)%vals = blks(k)%vals(bnx_start(k):bnx_end(k), bny_end(k))
                    call isend_real8(k, dist_block, dist_p, sync_buf8_send_nyp(k)%vals, send_size, tag, reqst)
                    !call isend_real8(k, dist_block, dist_p, blks(k)%vals(bnx_start(k):bnx_end(k), bny_end(k)), send_size, dist_block(2), reqst)
                    reqsts(icount) = reqst
                    icount = icount + 1
                else
                    lock = get_local_block_number(dist_block(1), dist_block(2))
                    blks(lock)%vals(bnx_start(lock):bnx_end(lock), bbnd_y1(lock) + 1) = blks(k)%vals(bnx_start(k):bnx_end(k), bny_end(k))
                endif
            endif
            ! isend in nx+
            dist_block(1) = bm + 1; dist_block(2) = bn
            call check_block_status(dist_block, dist_p)
            if (dist_p >= 0) then
                if (dist_p /= rank) then
                    send_size = bny_end(k) - bny_start(k) + 1
                    sync_buf8_send_nxp(k)%vals = blks(k)%vals(bnx_end(k), bny_start(k):bny_end(k))
                    call isend_real8(k, dist_block, dist_p, sync_buf8_send_nxp(k)%vals, send_size, tag, reqst)
                    !call isend_real8(k, dist_block, dist_p, blks(k)%vals(bnx_end(k), bny_start(k):bny_end(k)), send_size, dist_block(2), reqst)
                    reqsts(icount) = reqst
                    icount = icount + 1
                else
                    lock = get_local_block_number(dist_block(1), dist_block(2))
                    blks(lock)%vals(bbnd_x1(lock) + 1, bny_start(lock):bny_end(lock)) = blks(k)%vals(bnx_end(k), bny_start(k):bny_end(k))
                endif
            endif
            ! isend in ny-
            dist_block(1) = bm; dist_block(2) = bn - 1
            call check_block_status(dist_block, dist_p)
            if (dist_p >= 0) then
                if (dist_p /= rank) then
                    send_size = bnx_end(k) - bnx_start(k) + 1
                    sync_buf8_send_nym(k)%vals = blks(k)%vals(bnx_start(k):bnx_end(k), bny_start(k))
                    call isend_real8(k, dist_block, dist_p, sync_buf8_send_nym(k)%vals, send_size, tag, reqst)
                    !call isend_real8(k, dist_block, dist_p, blks(k)%vals(bnx_start(k):bnx_end(k), bny_start(k)), send_size, dist_block(2), reqst)
                    reqsts(icount) = reqst
                    icount = icount + 1
                else
                    lock = get_local_block_number(dist_block(1), dist_block(2))
                    blks(lock)%vals(bnx_start(lock):bnx_end(lock), bbnd_y2(lock) - 1) = blks(k)%vals(bnx_start(k):bnx_end(k), bny_start(k))
                endif
            endif
            ! isend in nx-
            dist_block(1) = bm - 1; dist_block(2) = bn
            call check_block_status(dist_block, dist_p)
            if (dist_p >= 0) then
                if (dist_p /= rank) then
                    send_size = bny_end(k) - bny_start(k) + 1
                    sync_buf8_send_nxm(k)%vals = blks(k)%vals(bnx_start(k), bny_start(k):bny_end(k))
                    call isend_real8(k, dist_block, dist_p, sync_buf8_send_nxm(k)%vals, send_size, tag, reqst)
                    !call isend_real8(k, dist_block, dist_p, blks(k)%vals(bnx_start(k), bny_start(k):bny_end(k)), send_size, dist_block(2), reqst)
                    reqsts(icount) = reqst
                    icount = icount + 1
                else
                    lock = get_local_block_number(dist_block(1), dist_block(2))
                    blks(lock)%vals(bbnd_x2(lock) - 1, bny_start(lock):bny_end(lock)) = blks(k)%vals(bnx_start(k), bny_start(k):bny_end(k))
                endif
            endif

            !----------------------------Edge points----------------------------!
            tag = 10*(bn + (bm - 1)*bnx)
            ! Edge isend in nx+ ny+
            dist_block(1) = bm + 1; dist_block(2) = bn + 1
            call check_block_status(dist_block, dist_p)
            if (dist_p >= 0) then
                if (dist_p /= rank) then
                    call mpi_isend(blks(k)%vals(bnx_end(k), bny_end(k)), 1, mpi_real8, dist_p, tag, cart_comm, reqst, ierr)
                    reqsts(icount) = reqst
                    icount = icount + 1
                else
                    lock = get_local_block_number(dist_block(1), dist_block(2))
                    blks(lock)%vals(bbnd_x1(lock) + 1, bbnd_y1(lock) + 1) = blks(k)%vals(bnx_end(k), bny_end(k))
                endif
            endif
            ! Edge isend in nx+ ny-
            dist_block(1) = bm + 1; dist_block(2) = bn - 1
            call check_block_status(dist_block, dist_p)
            if (dist_p >= 0) then
                if (dist_p /= rank) then
                    call mpi_isend(blks(k)%vals(bnx_end(k), bny_start(k)), 1, mpi_real8, dist_p, tag, cart_comm, reqst, ierr)
                    reqsts(icount) = reqst
                    icount = icount + 1
                else
                    lock = get_local_block_number(dist_block(1), dist_block(2))
                    blks(lock)%vals(bbnd_x1(lock) + 1, bbnd_y2(lock) - 1) = blks(k)%vals(bnx_end(k), bny_start(k))
                endif
            endif
            ! Edge isend in nx- ny+
            dist_block(1) = bm - 1; dist_block(2) = bn + 1
            call check_block_status(dist_block, dist_p)
            if (dist_p >= 0) then
                if (dist_p /= rank) then
                    call mpi_isend(blks(k)%vals(bnx_start(k), bny_end(k)), 1, mpi_real8, dist_p, tag, cart_comm, reqst, ierr)
                    reqsts(icount) = reqst
                    icount = icount + 1
                else
                    lock = get_local_block_number(dist_block(1), dist_block(2))
                    blks(lock)%vals(bbnd_x2(lock) - 1, bbnd_y1(lock) + 1) = blks(k)%vals(bnx_start(k), bny_end(k))
                endif
            endif
            ! Edge isend in nx- ny-
            dist_block(1) = bm - 1; dist_block(2) = bn - 1
            call check_block_status(dist_block, dist_p)
            if (dist_p >= 0) then
                if (dist_p /= rank) then
                    call mpi_isend(blks(k)%vals(bnx_start(k), bny_start(k)), 1, mpi_real8, dist_p, tag, cart_comm, reqst, ierr)
                    reqsts(icount) = reqst
                    icount = icount + 1
                else
                    lock = get_local_block_number(dist_block(1), dist_block(2))
                    blks(lock)%vals(bbnd_x2(lock) - 1, bbnd_y2(lock) - 1) = blks(k)%vals(bnx_start(k), bny_start(k))
                endif
            endif

        enddo
        if (parallel_dbg >= 2) print *, rank, 'icount totl:', icount-1

        ! Wait all, sync point
        call mpi_waitall(icount-1, reqsts, statuses, ierr)

        do k = 1, bcount
            bm = bindx(k, 1)
            bn = bindx(k, 2)

            ! ny+
            src_block(1) = bm; src_block(2) = bn - 1
            call check_block_status(src_block, src_p)
            if (src_p >= 0) then
                if (src_p /= rank) then
                    blks(k)%vals(bnx_start(k):bnx_end(k), bbnd_y1(k) + 1) = sync_buf8_recv_nyp(k)%vals
                endif
            endif
            ! nx+
            src_block(1) = bm - 1; src_block(2) = bn
            call check_block_status(src_block, src_p)
            if (src_p >= 0) then
                if (src_p /= rank) then
                    blks(k)%vals(bbnd_x1(k) + 1, bny_start(k):bny_end(k)) = sync_buf8_recv_nxp(k)%vals
                endif
            endif
            ! ny-
            src_block(1) = bm; src_block(2) = bn + 1
            call check_block_status(src_block, src_p)
            if (src_p >= 0) then
                if (src_p /= rank) then
                    blks(k)%vals(bnx_start(k):bnx_end(k), bbnd_y2(k) - 1) = sync_buf8_recv_nym(k)%vals
                endif
            endif
            ! nx-
            src_block(1) = bm + 1; src_block(2) = bn
            call check_block_status(src_block, src_p)
            if (src_p >= 0) then
                if (src_p /= rank) then
                    blks(k)%vals(bbnd_x2(k) - 1, bny_start(k):bny_end(k)) = sync_buf8_recv_nxm(k)%vals
                endif
            endif

            !----------------------------Edge points----------------------------!
            ! Edge in nx+ ny+
            src_block(1) = bm - 1; src_block(2) = bn - 1
            call check_block_status(src_block, src_p)
            if (src_p >= 0) then
                if (src_p /= rank) then
                    blks(k)%vals(bbnd_x1(k) + 1, bbnd_y1(k) + 1) = sync_edge_buf8_recv_nxp_nyp
                endif
            endif
            ! Edge in nx+ ny-
            src_block(1) = bm - 1; src_block(2) = bn + 1
            call check_block_status(src_block, src_p)
            if (src_p >= 0) then
                if (src_p /= rank) then
                    blks(k)%vals(bbnd_x1(k) + 1, bbnd_y2(k) - 1) = sync_edge_buf8_recv_nxp_nym
                endif
            endif
            ! Edge irecv in nx- ny+
            src_block(1) = bm + 1; src_block(2) = bn - 1
            call check_block_status(src_block, src_p)
            if (src_p >= 0) then
                if (src_p /= rank) then
                    blks(k)%vals(bbnd_x2(k) - 1, bbnd_y1(k) + 1) = sync_edge_buf8_recv_nxm_nyp
                endif
            endif
            ! Edge irecv in nx- ny-
            src_block(1) = bm + 1; src_block(2) = bn + 1
            call check_block_status(src_block, src_p)
            if (src_p >= 0) then
                if (src_p /= rank) then
                    blks(k)%vals(bbnd_x2(k) - 1, bbnd_y2(k) - 1) = sync_edge_buf8_recv_nxm_nym
                endif
            endif
        enddo

    end subroutine

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
    integer function get_rank_by_point(m, n)
        implicit none

        integer :: m, n
        integer :: flag_r, r, ierr

        flag_r = -1
        if (m >= nx_start .and. m <= nx_end) then
            if (n >= ny_start .and. n <= ny_end) then
                flag_r = rank
            endif
        endif

        call mpi_allreduce(flag_r, r, 1, mpi_integer,      &
                           mpi_max, cart_comm, ierr)

        get_rank_by_point = r
        return
    end function


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
    integer function check_cart_coord(coord, grid_size)
        implicit none
        integer, dimension(2), intent(in) :: coord, grid_size

        check_cart_coord = 0
        !write(*,*) coord,all(coord.ge.0),all((p_size-coord).ge.1)
        !print *, coord, p_size - coord, all((p_size-coord).ge.1)
        if (all(coord .ge. 0) .and. all((grid_size - coord) .ge. 1)) then
            check_cart_coord = 1
        endif
        return
    end function

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
    subroutine directsync_real8(field, p_dist, xd1, xd2, yd1, yd2,             &
                                       p_src,  xs1, xs2, ys1, ys2, nz)
        implicit none
        integer :: nz
        real*8, intent(in out) :: field(bnd_x1:bnd_x2, bnd_y1:bnd_y2, nz)
        integer, dimension(2), intent(in) :: p_dist, p_src
        integer :: xd1, xd2, yd1, yd2 ! bound of array which sending to p_dist
        integer :: xs1, xs2, ys1, ys2 ! bound of array which recieving from p_src

        integer :: dist_rank, src_rank
        integer :: flag_dist, flag_src
        integer :: ierr, debg
        integer :: stat(mpi_status_size)

        debg = 0

        if ( ((xd1-xd2+1)*(yd1-yd2+1)) .ne. (xs1-xs2+1)*(ys1-ys2+1) ) then
            print *, "Error in sync arrays size!"
        endif

        flag_dist = check_cart_coord(p_dist, p_size)
        flag_src = check_cart_coord(p_src, p_size)

        if ( (flag_src .eq. 1) .and. (flag_dist .eq. 1) ) then
            call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)
            call mpi_cart_rank(cart_comm, p_src, src_rank, ierr)

            call mpi_sendrecv(field(xd1:xd2, yd1:yd2, 1:nz),                          &
                              (xd2 - xd1 + 1)*(yd2 - yd1 + 1)*nz,                 &
                              mpi_real8, dist_rank, 1,                         &
                              field(xs1:xs2, ys1:ys2, 1:nz),                          &
                              (xs2 - xs1 + 1)*(ys2 - ys1 + 1)*nz,                 &
                              mpi_real8, src_rank, 1,                          &
                              cart_comm, stat, ierr)
!            print *, rank, "rsendecv", ierr
        else
            if (flag_src .eq. 1) then
                call mpi_cart_rank(cart_comm,p_src,src_rank,ierr)

                call mpi_recv(field(xs1:xs2, ys1:ys2, 1:nz),                          &
                              (xs2 - xs1 + 1)*(ys2 - ys1 + 1)*nz,                 &
                              mpi_real8, src_rank, 1,                          &
                              cart_comm, stat, ierr)
!                print *, rank, src_rank, "recv", xs1, xs2, ys1, ys2, field(xs1:xs2, ys1:ys2)
            endif

            if (flag_dist .eq. 1) then
                call mpi_cart_rank(cart_comm,p_dist,dist_rank,ierr)

                call mpi_send(field(xd1:xd2, yd1:yd2, 1:nz),                          &
                             (xd2 - xd1 + 1)*(yd2 - yd1 + 1)*nz,                  &
                             mpi_real8, dist_rank, 1,                          &
                             cart_comm, ierr)
!                print *, rank, dist_rank, "send", xd1, xd2, yd1, yd2, field(xd1:xd2, yd1:yd2)
            endif
        endif

    end subroutine directsync_real8

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
    subroutine syncborder_real8(field, nz)
        implicit none
        integer :: nz
        real*8, intent(in out) :: field(bnd_x1:bnd_x2, bnd_y1:bnd_y2, nz)

        integer, dimension(2) :: p_dist, p_src
        real*8 :: time_count

        !call start_timer(time_count)
!------------------ send-recv in ny+ -------------------------------------------
        p_dist(1) = p_coord(1)
        p_dist(2) = p_coord(2) + 1
        p_src(1) = p_coord(1)
        p_src(2) = p_coord(2) - 1
        call directsync_real8(field, p_dist, nx_start, nx_end, ny_end, ny_end,       &
                                     p_src,  nx_start, nx_end, bnd_y1 + 1, bnd_y1 + 1, nz)
!------------------ send-recv in nx+ -------------------------------------------
        p_dist(1) = p_coord(1) + 1
        p_dist(2) = p_coord(2)
        p_src(1) = p_coord(1) - 1
        p_src(2) = p_coord(2)
        call directsync_real8(field, p_dist, nx_end, nx_end, ny_start, ny_end,       &
                                     p_src,  bnd_x1 + 1, bnd_x1 + 1, ny_start, ny_end, nz)
!------------------ send-recv in ny- -------------------------------------------
        p_dist(1) = p_coord(1)
        p_dist(2) = p_coord(2) - 1
        p_src(1) = p_coord(1)
        p_src(2) = p_coord(2) + 1
        call directsync_real8(field, p_dist, nx_start, nx_end, ny_start, ny_start,   &
                                     p_src,  nx_start, nx_end, bnd_y2 - 1, bnd_y2 - 1, nz)
!------------------ send-recv in nx- -------------------------------------------
        p_dist(1) = p_coord(1) - 1
        p_dist(2) = p_coord(2)
        p_src(1) = p_coord(1) + 1
        p_src(2) = p_coord(2)
        call directsync_real8(field, p_dist, nx_start, nx_start, ny_start, ny_end,   &
                                     p_src,  bnd_x2 - 1, bnd_x2 - 1, ny_start, ny_end, nz)


!------------------ Sync edge points (EP) --------------------------------------
!------------------ send-recv EP in nx+,ny+ ------------------------------------
         p_dist(1) = p_coord(1) + 1
         p_dist(2) = p_coord(2) + 1
         p_src(1) = p_coord(1) - 1
         p_src(2) = p_coord(2) - 1
         call directsync_real8(field, p_dist, nx_end, nx_end, ny_end, ny_end,   &
                                      p_src,  bnd_x1 + 1, bnd_x1 + 1, bnd_y1 + 1, bnd_y1 + 1, nz)
!------------------ send-recv EP in nx+,ny- ------------------------------------
         p_dist(1) = p_coord(1) + 1
         p_dist(2) = p_coord(2) - 1
         p_src(1) = p_coord(1) - 1
         p_src(2) = p_coord(2) + 1
         call directsync_real8(field, p_dist, nx_end, nx_end, ny_start, ny_start,   &
                                      p_src,  bnd_x1 + 1, bnd_x1 + 1, bnd_y2 - 1 , bnd_y2 - 1, nz)
!------------------ send-recv EP in nx-,ny- ------------------------------------
         p_dist(1) = p_coord(1) - 1
         p_dist(2) = p_coord(2) - 1
         p_src(1) = p_coord(1) + 1
         p_src(2) = p_coord(2) + 1
         call directsync_real8(field, p_dist, nx_start, nx_start, ny_start, ny_start,   &
                                      p_src,  bnd_x2 - 1, bnd_x2 - 1, bnd_y2 - 1, bnd_y2 - 1, nz)

!------------------ send-recv EP in nx-,ny+ ------------------------------------
         p_dist(1) = p_coord(1) - 1
         p_dist(2) = p_coord(2) + 1
         p_src(1) = p_coord(1) + 1
         p_src(2) = p_coord(2) - 1
         call directsync_real8(field, p_dist, nx_start, nx_start, ny_end, ny_end,  &
                                      p_src,  bnd_x2 - 1, bnd_x2 - 1, bnd_y1 + 1, bnd_y1 + 1, nz)

        !call end_timer(time_count)
        !time_sync = time_sync + time_count
        return
    end subroutine syncborder_real8

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
    subroutine directsync_real(field, p_dist, xd1, xd2, yd1, yd2,             &
                                       p_src,  xs1, xs2, ys1, ys2, nz)
        implicit none
        integer :: nz
        real*4, intent(in out) :: field(bnd_x1:bnd_x2, bnd_y1:bnd_y2, nz)
        integer, dimension(2), intent(in) :: p_dist, p_src
        integer :: xd1, xd2, yd1, yd2 ! bound of array which sending to p_dist
        integer :: xs1, xs2, ys1, ys2 ! bound of array which recieving from p_src

        integer :: dist_rank, src_rank
        integer :: flag_dist, flag_src
        integer :: ierr, debg
        integer :: stat(mpi_status_size)

        debg = 0

        if ( ((xd1-xd2+1)*(yd1-yd2+1)) .ne. (xs1-xs2+1)*(ys1-ys2+1) ) then
            print *, "Error in sync arrays size!"
        endif

        flag_dist = check_cart_coord(p_dist, p_size)
        flag_src = check_cart_coord(p_src, p_size)

        if ( (flag_src .eq. 1) .and. (flag_dist .eq. 1) ) then
            call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)
            call mpi_cart_rank(cart_comm, p_src, src_rank, ierr)

            call mpi_sendrecv(field(xd1:xd2, yd1:yd2, 1:nz),                          &
                              (xd2 - xd1 + 1)*(yd2 - yd1 + 1)*nz,                 &
                              mpi_real4, dist_rank, 1,                         &
                              field(xs1:xs2, ys1:ys2, 1:nz),                          &
                              (xs2 - xs1 + 1)*(ys2 - ys1 + 1)*nz,                 &
                              mpi_real4, src_rank, 1,                          &
                              cart_comm, stat, ierr)
!            print *, rank, "rsendecv", ierr
        else
            if (flag_src .eq. 1) then
                call mpi_cart_rank(cart_comm,p_src,src_rank,ierr)

                call mpi_recv(field(xs1:xs2, ys1:ys2, 1:nz),                          &
                              (xs2 - xs1 + 1)*(ys2 - ys1 + 1)*nz,                 &
                              mpi_real4, src_rank, 1,                          &
                              cart_comm, stat, ierr)
!                print *, rank, src_rank, "recv", xs1, xs2, ys1, ys2, field(xs1:xs2, ys1:ys2)
            endif

            if (flag_dist .eq. 1) then
                call mpi_cart_rank(cart_comm,p_dist,dist_rank,ierr)

                call mpi_send(field(xd1:xd2, yd1:yd2, 1:nz),                          &
                             (xd2 - xd1 + 1)*(yd2 - yd1 + 1)*nz,                  &
                             mpi_real4, dist_rank, 1,                          &
                             cart_comm, ierr)
!                print *, rank, dist_rank, "send", xd1, xd2, yd1, yd2, field(xd1:xd2, yd1:yd2)
            endif
        endif
        return
    end subroutine directsync_real

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
    subroutine syncborder_real(field, nz)
        implicit none
        integer :: nz
        real*4, intent(in out) :: field(bnd_x1:bnd_x2, bnd_y1:bnd_y2, nz)

        integer, dimension(2) :: p_dist, p_src

!------------------ send-recv in ny+ -------------------------------------------
        p_dist(1) = p_coord(1)
        p_dist(2) = p_coord(2) + 1
        p_src(1) = p_coord(1)
        p_src(2) = p_coord(2) - 1
        call directsync_real(field, p_dist, nx_start, nx_end, ny_end, ny_end,       &
                                     p_src,  nx_start, nx_end, bnd_y1 + 1, bnd_y1 + 1, nz)
!------------------ send-recv in nx+ -------------------------------------------
        p_dist(1) = p_coord(1) + 1
        p_dist(2) = p_coord(2)
        p_src(1) = p_coord(1) - 1
        p_src(2) = p_coord(2)
        call directsync_real(field, p_dist, nx_end, nx_end, ny_start, ny_end,       &
                                     p_src,  bnd_x1 + 1, bnd_x1 + 1, ny_start, ny_end, nz)
!------------------ send-recv in ny- -------------------------------------------
        p_dist(1) = p_coord(1)
        p_dist(2) = p_coord(2) - 1
        p_src(1) = p_coord(1)
        p_src(2) = p_coord(2) + 1
        call directsync_real(field, p_dist, nx_start, nx_end, ny_start, ny_start,   &
                                     p_src,  nx_start, nx_end, bnd_y2 - 1, bnd_y2 - 1, nz)
!------------------ send-recv in nx- -------------------------------------------
        p_dist(1) = p_coord(1) - 1
        p_dist(2) = p_coord(2)
        p_src(1) = p_coord(1) + 1
        p_src(2) = p_coord(2)
        call directsync_real(field, p_dist, nx_start, nx_start, ny_start, ny_end,   &
                                     p_src,  bnd_x2 - 1, bnd_x2 - 1, ny_start, ny_end, nz)


!------------------ Sync edge points (EP) --------------------------------------
!------------------ send-recv EP in nx+,ny+ ------------------------------------
         p_dist(1) = p_coord(1) + 1
         p_dist(2) = p_coord(2) + 1
         p_src(1) = p_coord(1) - 1
         p_src(2) = p_coord(2) - 1
         call directsync_real(field, p_dist, nx_end, nx_end, ny_end, ny_end,   &
                                      p_src,  bnd_x1 + 1, bnd_x1 + 1, bnd_y1 + 1, bnd_y1 + 1, nz)
!------------------ send-recv EP in nx+,ny- ------------------------------------
         p_dist(1) = p_coord(1) + 1
         p_dist(2) = p_coord(2) - 1
         p_src(1) = p_coord(1) - 1
         p_src(2) = p_coord(2) + 1
         call directsync_real(field, p_dist, nx_end, nx_end, ny_start, ny_start,   &
                                      p_src,  bnd_x1 + 1, bnd_x1 + 1, bnd_y2 - 1 , bnd_y2 - 1, nz)
!------------------ send-recv EP in nx-,ny- ------------------------------------
         p_dist(1) = p_coord(1) - 1
         p_dist(2) = p_coord(2) - 1
         p_src(1) = p_coord(1) + 1
         p_src(2) = p_coord(2) + 1
         call directsync_real(field, p_dist, nx_start, nx_start, ny_start, ny_start,   &
                                      p_src,  bnd_x2 - 1, bnd_x2 - 1, bnd_y2 - 1, bnd_y2 - 1, nz)

!------------------ send-recv EP in nx-,ny+ ------------------------------------
         p_dist(1) = p_coord(1) - 1
         p_dist(2) = p_coord(2) + 1
         p_src(1) = p_coord(1) + 1
         p_src(2) = p_coord(2) - 1
         call directsync_real(field, p_dist, nx_start, nx_start, ny_end, ny_end,  &
                                      p_src,  bnd_x2 - 1, bnd_x2 - 1, bnd_y1 + 1, bnd_y1 + 1, nz)

        return
    end subroutine syncborder_real


endmodule mpi_parallel_tools
!---------------------end module for definition of array dimensions and boundaries-----------------
