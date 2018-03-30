!-----------module for definition of array dimensions and boundaries------------
module mpi_parallel_tools
    use mpi
    implicit none

    !include 'mpif.h'
    include "omp_lib.h"

    integer :: recommended_tot_blocks

    integer :: nx_start, nx_end, ny_start, ny_end
    integer :: bnd_x1, bnd_x2, bnd_y1, bnd_y2

    type :: block2D
        real*8, dimension(:, :), pointer :: vals
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
    integer :: sync_buff_size
    real*8, allocatable :: sync_buf8_send_nyp(:, :), sync_buf8_recv_nyp(:, :)

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

            close(90)
        endif
        call mpi_bcast(recommended_tot_blocks, 1, mpi_integer, 0, mpi_comm_world, ierr)

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

        ! Allocate buffers
        allocate(reqsts(8*procs), statuses(MPI_STATUS_SIZE, 8*procs))

    end subroutine

    subroutine parallel_finalize
        implicit none

        integer :: ierr

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

        if (allocated(sync_buf8_send_nyp)) deallocate(sync_buf8_send_nyp)
        if (allocated(sync_buf8_recv_nyp)) deallocate(sync_buf8_recv_nyp)

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
        print *, rank, 'bnx, bny and Total blocks:', bnx, bny, bnx*bny
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
        print *, rank, 'loc_bnx, loc_bny and Blocks per proc: ', loc_bnx, loc_bny, bcount
        call mpi_barrier(cart_comm, ierr)

        xblock_start = 1 + p_coord(1)*loc_bnx
        yblock_start = 1 + p_coord(2)*loc_bny
        print *, rank, 'xb_start, yb_start', xblock_start, yblock_start
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
                !print *, rank, 'k, bnx_start, bnx_end', k, bnx_start(k), bnx_end(k)
                !print *, rank, 'k, bny_start, bny_end', k, bny_start(k), bny_end(k)
                do i = bnx_start(k), bnx_end(k)
                    do j = bny_start(k), bny_end(k)
                        bweight(k) = bweight(k) + (1.0d0 - real(lbasins(i, j)))
                    enddo
                enddo
                weight = weight + real(bweight(k)) / real((bnx_end(k) - bnx_start(k) + 1)*(bny_end(k) - bny_start(k) + 1))
                ! Ignore only-land blocks
                if (bweight(k) == 0.0) then
                    land_blocks = land_blocks + 1
                else
                    k = k + 1
                endif
            enddo
        enddo

        print *, rank, 'Total blocks with land:', land_blocks, 'Total Load-Balancing:', weight / bcount
        call mpi_barrier(cart_comm, ierr)
        bcount = bcount - land_blocks
        print *, rank, 'Updated block count', bcount
        call mpi_barrier(cart_comm, ierr)
        ierr = 0
        if (bcount <= 0) then
            print *, rank, 'All blocks with land blocks!'
            ierr = 1
        endif
        call parallel_check_err(ierr)
        call mpi_allreduce(bcount, total_blocks, 1, mpi_integer,      &
                           mpi_sum, cart_comm, ierr)
        print *, rank, 'Updated total block count', total_blocks
        call mpi_barrier(cart_comm, ierr)

        ! Sync bglob_proc array
        call parallel_int_output(bglob_proc, 1, bnx, 1, bny, 'bglob_proc before')
        allocate(buf_int(bnx, bny))
        buf_int = bglob_proc + 1
        call mpi_allreduce(buf_int, bglob_proc, bnx*bny, mpi_integer, mpi_sum, cart_comm, ierr)
        bglob_proc = bglob_proc - 1
        call parallel_int_output(bglob_proc, 1, bnx, 1, bny, 'bglob_proc after')

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

        integer :: k, maxnx, maxny

        maxnx = 0; maxny = 0
        do k = 1, bcount
            if ((bnx_end(k) - bnx_start(k) + 1) > maxnx) then
                maxnx = bnx_end(k) - bnx_start(k) + 1
            endif
            if ((bny_end(k) - bny_start(k) + 1) > maxny) then
                maxny = bny_end(k) - bny_start(k) + 1
            endif
        enddo

        sync_buff_size = max(maxnx, maxny)
        print *, rank, 'allocate buffers. maxnx, maxny, buffsize: ', maxnx, maxny, sync_buff_size
        allocate(sync_buf8_send_nyp(bcount, sync_buff_size), sync_buf8_recv_nyp(bcount, sync_buff_size))

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

    subroutine irecv_real8(k, src_block, buff, buff_size, tag, reqst, stat)
        implicit none
        integer :: k, tag
        integer, dimension(2) :: src_block
        integer :: buff_size
        real*8 :: buff(:)
        integer :: reqst, stat
        integer :: p_src, ierr, flg_recv
        integer, dimension(2) :: bgrid
        stat = 0

        bgrid(1) = bnx; bgrid(2) = bny
        flg_recv = check_cart_coord(src_block-1, bgrid)
        if (flg_recv /= 1) then
            print *, rank, 'NOIRECV - block dont exist:', src_block(1), src_block(2)
            stat = 1
            return
        endif

        p_src = bglob_proc(src_block(1), src_block(2))
        if (p_src < 0) then
            print *, rank, 'NOIRECV - land block: ', src_block(1), src_block(2)
            stat = 1
            return
        endif

        if (p_src /= rank) then
            print *, rank, 'IRECV. block: ', bindx(k, 1), bindx(k, 2), 'src_block:', src_block(1), src_block(2), 'src_p:', p_src
            !call mpi_irecv(blks(k)%vals(dx1:dx2, dy1:dy2), (dx2 - dx1 + 1)*(dy2 - dy1 + 1), &
            !               mpi_real8, p_src, tag, cart_comm, reqst, ierr)
            call mpi_irecv(buff, buff_size, &
                           mpi_real8, p_src, tag, cart_comm, reqst, ierr)
        else
            ! Just Copy
            print *, rank, 'IRECV JUST COPY. block: ', bindx(k, 1), bindx(k, 2), 'src_block:', src_block(1), src_block(2), 'src_p:', p_src
            stat = 1
        endif
    end subroutine

    subroutine isend_real8(k, dist_block, buff, buff_size, tag, reqst, stat)
        implicit none

        integer :: k, tag
        integer, dimension(2) :: dist_block
        integer :: buff_size
        real*8 :: buff(:)
        integer :: reqst, stat
        integer :: p_dist, ierr, flg_send
        integer, dimension(2) :: bgrid

        stat = 0

        bgrid(1) = bnx; bgrid(2) = bny
        flg_send = check_cart_coord(dist_block-1, bgrid)
        if (flg_send /= 1) then
            print *, rank, 'NOISEND - block dont exist:', dist_block(1), dist_block(2)
            stat = 1
            return
        endif

        p_dist = bglob_proc(dist_block(1), dist_block(2))
        if (p_dist < 0) then
            print *, rank, 'NOISEND - land block: ', dist_block(1), dist_block(2)
            stat = 1
            return
        endif

        if (p_dist /= rank) then
            print *, rank, 'ISEND. block: ', bindx(k, 1), bindx(k, 2), 'dst_block:', dist_block(1), dist_block(2), 'dst_p:', p_dist
            !call mpi_isend(blks(k)%vals(sx1:sx2, sy1:sy2), (sx2 - sx1 + 1)*(sy2 - sy1 + 1), &
            !               mpi_real8, p_dist, tag, cart_comm, reqst, ierr)
            call mpi_isend(buff, buff_size,  &
                           mpi_real8, p_dist, tag, cart_comm, reqst, ierr)
        else
            ! Just Copy
            print *, rank, 'ISEND JUST COPY. block: ', bindx(k, 1), bindx(k, 2), 'dst_block:', dist_block(1), dist_block(2), 'dst_p:', p_dist
            stat = 1

        endif

    end subroutine

    subroutine syncborder_block2D(blks)
        implicit none

        type(block2D), dimension(:), pointer :: blks
        integer :: k, bm, bn, stat, reqst
        integer :: ierr
        integer, dimension(2) :: src_block, dist_block
        integer :: send_size, recv_size
        integer :: icount

        icount = 1

        ! Non-blocking Recv calls
        do k = 1, bcount
            bm = bindx(k, 1)
            bn = bindx(k, 2)

            ! irecv in ny+
            src_block(1) = bm; src_block(2) = bn - 1
            recv_size = bnx_end(k) - bnx_start(k) + 1
            call irecv_real8(k, src_block, blks(k)%vals(bnx_start(k):bnx_end(k), bbnd_y1(k) + 1), recv_size, bn, reqst, stat)
            if (stat == 0) then
                reqsts(icount) = reqst
                icount = icount + 1
            endif

            ! irecv in nx+
            ! irecv in ny-
            ! irecv in nx-
        enddo
        print *, rank, 'icount recv:', icount-1

        ! Non-blocking Send calls
        do k = 1, bcount
            bm = bindx(k, 1)
            bn = bindx(k, 2)

            ! isend in ny+
            dist_block(1) = bm; dist_block(2) = bn + 1
            send_size = bnx_end(k) - bnx_start(k) + 1
            call isend_real8(k, dist_block, blks(k)%vals(bnx_start(k):bnx_end(k), bny_end(k)), send_size, dist_block(2), reqst, stat)
            if (stat == 0) then
                reqsts(icount) = reqst
                icount = icount + 1
            endif

            ! isend in nx+
            ! isend in ny-
            ! isend in nx-
        enddo
        print *, rank, 'icount total:', icount-1

        ! Wait all, sync point
        !MPI_WAITALL(COUNT, ARRAY_OF_REQUESTS, ARRAY_OF_STATUSES, IERROR)
        !INTEGER    COUNT, ARRAY_OF_REQUESTS(*)
        !INTEGER    ARRAY_OF_STATUSES(MPI_STATUS_SIZE,*), IERROR
        call mpi_waitall(icount-1, reqsts, statuses, ierr)
        !deallocate(reqsts, statuses)

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
