module mpi_parallel_tools
    use mpi
    implicit none

    !include 'mpif.h'
    include "omp_lib.h"

    integer :: nx_start,    &   !first significant point in x-direction
               nx_end,      &   !last  significant point in x-direction
               ny_start,    &   !first significant point in y-direction
               ny_end           !last  significant point in y-direction

    integer :: bnd_x1,      &   !left   array boundary in x-direction
               bnd_x2,      &   !right  array boundary in x-direction
               bnd_y1,      &   !bottom array boundary in y-direction
               bnd_y2           !top    array boundary in y-direction

    integer :: rank, procs
    integer :: cart_comm, local_output_comm
    integer, dimension(2) :: p_size, p_coord
    logical, dimension(2) :: period

    integer, allocatable :: comput_domain(:)

    real*8 :: time_barotrop
    real*8 :: time_model_step, time_output
    real*8 :: time_sync

    contains

    subroutine parallel_init()
        use main_basin_pars
        
        implicit none
        
        include 'locout.fi'

        integer :: ierr, i
        integer :: locn
        integer :: count_threads, num_thread
        integer :: colr

        call mpi_init(ierr)

        period = (/.true., .true./)
        p_size = (/0,0/)
        ierr = 0

        call mpi_comm_rank(mpi_comm_world, rank, ierr)
        call mpi_comm_size(mpi_comm_world, procs, ierr)
        call mpi_dims_create(procs, 2, p_size, ierr)
        call mpi_cart_create(mpi_comm_world, 2, p_size, period, .false., cart_comm, ierr)
        call mpi_cart_coords(cart_comm, rank, 2, p_coord, ierr)

        !-----------------------------------NX----------------------------------
        locn = floor(real(nx - 4)/real(p_size(1)))
        nx_start = locn*p_coord(1) + 1 + 2
        if ( p_coord(1) .EQ. p_size(1) - 1 ) then
            locn = (nx - 2) - nx_start + 1
        endif
        nx_end = nx_start + locn - 1
        nx_start = nx_start
        ! border area
        bnd_x1 = nx_start - 2
        !if (bnd_x1 < 1) bnd_x1 = 1
        bnd_x2 = nx_end + 2
        !if (bnd_x2 > nx) bnd_x2 = nx

        !-----------------------------------NY----------------------------------
        locn = floor(real(ny - 4)/real(p_size(2)))
        ny_start = locn*p_coord(2) + 1 + 2
        if ( p_coord(2) .EQ. p_size(2) - 1 ) then
            locn = (ny - 2) - ny_start + 1
        endif
        ny_end = ny_start + locn - 1
        ny_start = ny_start
        ! border area
        bnd_y1 = ny_start - 2
        !if (bnd_y1 < 1) bnd_y1 = 1
        bnd_y2 = ny_end + 2
        !if (bnd_y2 > ny) bnd_y2 = ny

        if (rank .eq. 0) print *, "MPI pocs: ", procs, " Domain decomposition:"
        do i = 0, procs-1
            if (rank .eq. i) then
                print *, "nx ", rank, p_coord, nx_start, nx_end, ny_start, ny_end
                print *, "bnd", rank, p_coord, bnd_x1, bnd_x2, bnd_y1, bnd_y2
            endif
            call mpi_barrier(cart_comm, ierr)
        enddo

        !$omp parallel
        count_threads = omp_get_num_threads()
        num_thread = omp_get_thread_num()
        if (num_thread .eq. 0) print *, "OMP Threads: ", count_threads
        !$omp end parallel
        call mpi_barrier(cart_comm, ierr)

        allocate(comput_domain(procs))
        comput_domain = 1

        ! LOCAL OUTPUT COMM
        if (rank == 0) print *, "Create output communicator for LOCAL area output only ..."
        colr = 1
        if (nx_end  < m1loc_local .or. nx_start > m2loc_local) colr = mpi_undefined
        if (ny_end  < n1loc_local .or. ny_start > n2loc_local) colr = mpi_undefined
        call mpi_comm_split(cart_comm, colr, rank, local_output_comm, ierr)
        !print *, rank, colr

    end subroutine

    subroutine parallel_finalize()
        implicit none

        integer :: ierr

        if (allocated(comput_domain)) deallocate(comput_domain)

        call mpi_finalize(ierr)

    end subroutine

    subroutine init_computational_domains(lu)
        implicit none

        real*4 :: lu(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
        real*4 :: weigth
        integer :: m, n, ierr
        integer, allocatable :: int_buf(:)

        if (rank == 0) print *, 'INIT COMPUTATIONAL DOMAINS ...'

        allocate(int_buf(procs))
        int_buf = 0

        weigth = 0
        do n = ny_start, ny_end
            do m = nx_start, nx_end
                weigth = weigth + lu(m, n)
            enddo
        enddo

        if (weigth == 0.0) then
            int_buf(rank+1) = 0
        else
            int_buf(rank+1) = 1
        endif

        call mpi_allreduce(int_buf, comput_domain, procs, mpi_integer, mpi_sum, cart_comm, ierr)

        if (rank == 0) print *, 'Total procs with computational domains:', sum(comput_domain)
        call mpi_barrier(cart_comm, ierr)

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
    integer function check_p_coord(coord)
        implicit none
        integer, dimension(2), intent(in) :: coord

        integer :: coord_rank, ierr

        check_p_coord = 0
        !write(*,*) coord,all(coord.ge.0),all((p_size-coord).ge.1)
        !print *, coord, p_size - coord, all((p_size-coord).ge.1)
        if (all(coord.ge.0) .and. all((p_size-coord).ge.1)) then
            call mpi_cart_rank(cart_comm, coord, coord_rank, ierr)
            if (comput_domain(coord_rank + 1) > 0) check_p_coord = 1
        endif

        return
    end function

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
        integer :: ierr
        integer :: stat(mpi_status_size)

        if ( ((xd1-xd2+1)*(yd1-yd2+1)) .ne. (xs1-xs2+1)*(ys1-ys2+1) ) then
            print *, rank, "Error in sync arrays size!"
        endif

        flag_dist = check_p_coord(p_dist)
        flag_src = check_p_coord(p_src)

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
        else
            if (flag_src .eq. 1) then
                call mpi_cart_rank(cart_comm,p_src,src_rank,ierr)

                call mpi_recv(field(xs1:xs2, ys1:ys2, 1:nz),                          &
                              (xs2 - xs1 + 1)*(ys2 - ys1 + 1)*nz,                 &
                              mpi_real8, src_rank, 1,                          &
                              cart_comm, stat, ierr)
            endif

            if (flag_dist .eq. 1) then
                call mpi_cart_rank(cart_comm,p_dist,dist_rank,ierr)

                call mpi_send(field(xd1:xd2, yd1:yd2, 1:nz),                          &
                             (xd2 - xd1 + 1)*(yd2 - yd1 + 1)*nz,                  &
                             mpi_real8, dist_rank, 1,                          &
                             cart_comm, ierr)
            endif
        endif

    end subroutine directsync_real8

!-------------------------------------------------------------------------------
    subroutine syncborder_real8(field, nz)
        implicit none
        integer :: nz
        real*8, intent(in out) :: field(bnd_x1:bnd_x2, bnd_y1:bnd_y2, nz)

        integer, dimension(2) :: p_dist, p_src
        real*8 :: time_count

        if (comput_domain(rank + 1) == 0.0) then
            return
        endif

        !call start_timer(time_count)
        !------------------ send-recv in ny+ -----------------------------------
        p_dist(1) = p_coord(1)
        p_dist(2) = p_coord(2) + 1
        p_src(1) = p_coord(1)
        p_src(2) = p_coord(2) - 1
        call directsync_real8(field, p_dist, nx_start, nx_end, ny_end, ny_end,       &
                                     p_src,  nx_start, nx_end, bnd_y1 + 1, bnd_y1 + 1, nz)
        !------------------ send-recv in nx+ -----------------------------------
        p_dist(1) = p_coord(1) + 1
        p_dist(2) = p_coord(2)
        p_src(1) = p_coord(1) - 1
        p_src(2) = p_coord(2)
        call directsync_real8(field, p_dist, nx_end, nx_end, ny_start, ny_end,       &
                                     p_src,  bnd_x1 + 1, bnd_x1 + 1, ny_start, ny_end, nz)
        !------------------ send-recv in ny- -----------------------------------
        p_dist(1) = p_coord(1)
        p_dist(2) = p_coord(2) - 1
        p_src(1) = p_coord(1)
        p_src(2) = p_coord(2) + 1
        call directsync_real8(field, p_dist, nx_start, nx_end, ny_start, ny_start,   &
                                     p_src,  nx_start, nx_end, bnd_y2 - 1, bnd_y2 - 1, nz)
        !------------------ send-recv in nx- -----------------------------------
        p_dist(1) = p_coord(1) - 1
        p_dist(2) = p_coord(2)
        p_src(1) = p_coord(1) + 1
        p_src(2) = p_coord(2)
        call directsync_real8(field, p_dist, nx_start, nx_start, ny_start, ny_end,   &
                                     p_src,  bnd_x2 - 1, bnd_x2 - 1, ny_start, ny_end, nz)


         !------------------ Sync edge points (EP) -----------------------------
         !------------------ send-recv EP in nx+,ny+ ---------------------------
         p_dist(1) = p_coord(1) + 1
         p_dist(2) = p_coord(2) + 1
         p_src(1) = p_coord(1) - 1
         p_src(2) = p_coord(2) - 1
         call directsync_real8(field, p_dist, nx_end, nx_end, ny_end, ny_end,   &
                                      p_src,  bnd_x1 + 1, bnd_x1 + 1, bnd_y1 + 1, bnd_y1 + 1, nz)
         !------------------ send-recv EP in nx+,ny- ---------------------------
         p_dist(1) = p_coord(1) + 1
         p_dist(2) = p_coord(2) - 1
         p_src(1) = p_coord(1) - 1
         p_src(2) = p_coord(2) + 1
         call directsync_real8(field, p_dist, nx_end, nx_end, ny_start, ny_start,   &
                                      p_src,  bnd_x1 + 1, bnd_x1 + 1, bnd_y2 - 1 , bnd_y2 - 1, nz)
         !------------------ send-recv EP in nx-,ny- ---------------------------
         p_dist(1) = p_coord(1) - 1
         p_dist(2) = p_coord(2) - 1
         p_src(1) = p_coord(1) + 1
         p_src(2) = p_coord(2) + 1
         call directsync_real8(field, p_dist, nx_start, nx_start, ny_start, ny_start,   &
                                      p_src,  bnd_x2 - 1, bnd_x2 - 1, bnd_y2 - 1, bnd_y2 - 1, nz)

         !------------------ send-recv EP in nx-,ny+ ---------------------------
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
        integer :: ierr
        integer :: stat(mpi_status_size)

        if ( ((xd1-xd2+1)*(yd1-yd2+1)) .ne. (xs1-xs2+1)*(ys1-ys2+1) ) then
            print *, rank, "Error in sync arrays size!"
        endif

        flag_dist = check_p_coord(p_dist)
        flag_src = check_p_coord(p_src)

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
        else
            if (flag_src .eq. 1) then
                call mpi_cart_rank(cart_comm,p_src,src_rank,ierr)

                call mpi_recv(field(xs1:xs2, ys1:ys2, 1:nz),                          &
                              (xs2 - xs1 + 1)*(ys2 - ys1 + 1)*nz,                 &
                              mpi_real4, src_rank, 1,                          &
                              cart_comm, stat, ierr)
            endif

            if (flag_dist .eq. 1) then
                call mpi_cart_rank(cart_comm,p_dist,dist_rank,ierr)

                call mpi_send(field(xd1:xd2, yd1:yd2, 1:nz),                          &
                             (xd2 - xd1 + 1)*(yd2 - yd1 + 1)*nz,                  &
                             mpi_real4, dist_rank, 1,                          &
                             cart_comm, ierr)
            endif
        endif
        return
    end subroutine directsync_real

!-------------------------------------------------------------------------------
    subroutine syncborder_real(field, nz)
        implicit none
        integer :: nz
        real*4, intent(in out) :: field(bnd_x1:bnd_x2, bnd_y1:bnd_y2, nz)

        integer, dimension(2) :: p_dist, p_src

        if (comput_domain(rank + 1) == 0.0) then
            return
        endif

        !------------------ send-recv in ny+ -----------------------------------
        p_dist(1) = p_coord(1)
        p_dist(2) = p_coord(2) + 1
        p_src(1) = p_coord(1)
        p_src(2) = p_coord(2) - 1
        call directsync_real(field, p_dist, nx_start, nx_end, ny_end, ny_end,       &
                                     p_src,  nx_start, nx_end, bnd_y1 + 1, bnd_y1 + 1, nz)
        !------------------ send-recv in nx+ -----------------------------------
        p_dist(1) = p_coord(1) + 1
        p_dist(2) = p_coord(2)
        p_src(1) = p_coord(1) - 1
        p_src(2) = p_coord(2)
        call directsync_real(field, p_dist, nx_end, nx_end, ny_start, ny_end,       &
                                     p_src,  bnd_x1 + 1, bnd_x1 + 1, ny_start, ny_end, nz)
        !------------------ send-recv in ny- -----------------------------------
        p_dist(1) = p_coord(1)
        p_dist(2) = p_coord(2) - 1
        p_src(1) = p_coord(1)
        p_src(2) = p_coord(2) + 1
        call directsync_real(field, p_dist, nx_start, nx_end, ny_start, ny_start,   &
                                     p_src,  nx_start, nx_end, bnd_y2 - 1, bnd_y2 - 1, nz)
        !------------------ send-recv in nx- -----------------------------------
        p_dist(1) = p_coord(1) - 1
        p_dist(2) = p_coord(2)
        p_src(1) = p_coord(1) + 1
        p_src(2) = p_coord(2)
        call directsync_real(field, p_dist, nx_start, nx_start, ny_start, ny_end,   &
                                     p_src,  bnd_x2 - 1, bnd_x2 - 1, ny_start, ny_end, nz)


         !------------------ Sync edge points (EP) -----------------------------
         !------------------ send-recv EP in nx+,ny+ ---------------------------
         p_dist(1) = p_coord(1) + 1
         p_dist(2) = p_coord(2) + 1
         p_src(1) = p_coord(1) - 1
         p_src(2) = p_coord(2) - 1
         call directsync_real(field, p_dist, nx_end, nx_end, ny_end, ny_end,   &
                                      p_src,  bnd_x1 + 1, bnd_x1 + 1, bnd_y1 + 1, bnd_y1 + 1, nz)
         !------------------ send-recv EP in nx+,ny- ---------------------------
         p_dist(1) = p_coord(1) + 1
         p_dist(2) = p_coord(2) - 1
         p_src(1) = p_coord(1) - 1
         p_src(2) = p_coord(2) + 1
         call directsync_real(field, p_dist, nx_end, nx_end, ny_start, ny_start,   &
                                      p_src,  bnd_x1 + 1, bnd_x1 + 1, bnd_y2 - 1 , bnd_y2 - 1, nz)
         !------------------ send-recv EP in nx-,ny- ---------------------------
         p_dist(1) = p_coord(1) - 1
         p_dist(2) = p_coord(2) - 1
         p_src(1) = p_coord(1) + 1
         p_src(2) = p_coord(2) + 1
         call directsync_real(field, p_dist, nx_start, nx_start, ny_start, ny_start,   &
                                      p_src,  bnd_x2 - 1, bnd_x2 - 1, bnd_y2 - 1, bnd_y2 - 1, nz)

         !------------------ send-recv EP in nx-,ny+ ---------------------------
         p_dist(1) = p_coord(1) - 1
         p_dist(2) = p_coord(2) + 1
         p_src(1) = p_coord(1) + 1
         p_src(2) = p_coord(2) - 1
         call directsync_real(field, p_dist, nx_start, nx_start, ny_end, ny_end,  &
                                      p_src,  bnd_x2 - 1, bnd_x2 - 1, bnd_y1 + 1, bnd_y1 + 1, nz)

        return
    end subroutine syncborder_real


endmodule mpi_parallel_tools
