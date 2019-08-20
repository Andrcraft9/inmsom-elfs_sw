!---------------------module for definition of basin grid arrays-----------------
module basin_grid
    use mpi_parallel_tools
    implicit none

    include 'basinpar.fi'

    ! Mask of entire basin area
    ! Temporary pointers
    real(4), dimension(:, :), pointer :: lu,        &  !mask of t-grid
                                         lu1,       &  !mask of t-grid (1 everywhere)
                                         luu,       &  !mask of h-grid (0 on boundary)
                                         luh,       &  !mask of h-grid (1 on boundary)
                                         lcu,       &  !mask of u-grid (0 on boundary)
                                         lcv,       &  !mask of v-grid (0 on boundary)
                                         llu,       &  !mask of u-grid (1 on boundary)
                                         llv           !mask of v-grid (1 on boundary)
    ! Main arrays
    type(block2D_real4), dimension(:), pointer :: block_lu,        &  !mask of t-grid
                                                  block_lu1,       &  !mask of t-grid (1 everywhere)
                                                  block_luu,       &  !mask of h-grid (0 on boundary)
                                                  block_luh,       &  !mask of h-grid (1 on boundary)
                                                  block_lcu,       &  !mask of u-grid (0 on boundary)
                                                  block_lcv,       &  !mask of v-grid (0 on boundary)
                                                  block_llu,       &  !mask of u-grid (1 on boundary)
                                                  block_llv           !mask of v-grid (1 on boundary)

    ! Mask of local basin area (regional).
    ! Used only for local (regional) output
    ! Temporary pointers
    real(4), dimension(:, :), pointer :: lu_output,        &  !mask of t-grid
                                         lcu_output,       &  !mask of u-grid (0 on boundary)
                                         lcv_output,       &  !mask of v-grid (0 on boundary)
                                         llu_output,       &  !mask of u-grid (1 on boundary)
                                         llv_output           !mask of v-grid (1 on boundary)
    ! Main arrays
    type(block2D_real4), dimension(:), pointer :: block_lu_output,        &  !mask of t-grid
                                                  block_lcu_output,       &  !mask of u-grid (0 on boundary)
                                                  block_lcv_output,       &  !mask of v-grid (0 on boundary)
                                                  block_llu_output,       &  !mask of u-grid (1 on boundary)
                                                  block_llv_output           !mask of v-grid (1 on boundary)

    ! Temporary pointers
    real(8), dimension(:, :), pointer :: hhh,      &  !ocean depth on luh (h-points)
                                         hhhp,     &  !ocean depth on luh (h-points) at previous step
                                         hhhn,     &  !ocean depth on luh (h-points) at next step
                                         hhq_rest,     &  !ocean depth on lu  (t-points) at rest state
                                         hhq,      &  !ocean depth on lu  (t-points)
                                         hhqp,     &  !ocean depth on lu  (t-points) at previous step
                                         hhqn,     &  !ocean depth on lu  (t-points) at next step
                                         hhu,      &  !ocean depth on lcu (u-points)
                                         hhup,     &  !ocean depth on lcu (u-points) at previous step
                                         hhun,     &  !ocean depth on lcu (u-points) at next step
                                         hhv,      &  !ocean depth on lcv (v-points)
                                         hhvp,     &  !ocean depth on lcv (v-points) at previous step
                                         hhvn,     &  !ocean depth on lcv (v-points) at next step
                                         rlh_s,    &  !main (sin) coriolis parameter on h-points
                                         rlh_c,    &  !2-nd (cos) coriolis parameter on h-points
                                         dxt, dyt,      &  !horizontal grid steps between   t-points (in meters)
                                         dx , dy ,      &  !horizontal grid steps between u,v-points (in meters)
                                         dxh, dyh,      &  !horizontal grid steps between   h-points (in meters)
                                         dxb, dyb          !horizontal grid steps between v,u-points (in meters)
    ! Main arrays
    type(block2D_real8), dimension(:), pointer :: block_hhh,      &  !ocean depth on luh (h-points)
                                                  block_hhhp,     &  !ocean depth on luh (h-points) at previous step
                                                  block_hhhn,     &  !ocean depth on luh (h-points) at next step
                                                  block_hhq_rest,     &  !ocean depth on lu  (t-points) at rest state
                                                  block_hhq,      &  !ocean depth on lu  (t-points)
                                                  block_hhqp,     &  !ocean depth on lu  (t-points) at previous step
                                                  block_hhqn,     &  !ocean depth on lu  (t-points) at next step
                                                  block_hhu,      &  !ocean depth on lcu (u-points)
                                                  block_hhup,     &  !ocean depth on lcu (u-points) at previous step
                                                  block_hhun,     &  !ocean depth on lcu (u-points) at next step
                                                  block_hhv,      &  !ocean depth on lcv (v-points)
                                                  block_hhvp,     &  !ocean depth on lcv (v-points) at previous step
                                                  block_hhvn,     &  !ocean depth on lcv (v-points) at next step
                                                  block_rlh_s,    &  !main (sin) coriolis parameter on h-points
                                                  block_rlh_c,    &  !2-nd (cos) coriolis parameter on h-points
                                                  block_dxt, block_dyt,      &  !horizontal grid steps between   t-points (in meters)
                                                  block_dx , block_dy ,      &  !horizontal grid steps between u,v-points (in meters)
                                                  block_dxh, block_dyh,      &  !horizontal grid steps between   h-points (in meters)
                                                  block_dxb, block_dyb          !horizontal grid steps between v,u-points (in meters)
    ! Temporary pointers
    real(8), dimension(:), pointer :: xt, yt,        &  !horizontal t-grid            x- and y-coordinates (in degrees)
                                      xu, yv            !horizontal u-grid and v-grid x- and y-coordinates (in degrees)
    ! Main arrays
    type(block1D_real8), dimension(:), pointer :: block_xt, block_yt,        &  !horizontal t-grid            x- and y-coordinates (in degrees)
                                                  block_xu, block_yv            !horizontal u-grid and v-grid x- and y-coordinates (in degrees)
    ! Temporary pointers
    real(8), dimension(:, :), pointer :: geo_lon_t,   &    !geographical longitudes of T-points
                                         geo_lat_t,   &    !geographical latitudes  of T-points
                                         geo_lon_u,   &    !geographical longitudes of U-points
                                         geo_lat_u,   &    !geographical latitudes  of U-points
                                         geo_lon_v,   &    !geographical longitudes of V-points
                                         geo_lat_v,   &    !geographical latitudes  of V-points
                                         geo_lon_h,   &    !geographical longitudes of H-points
                                         geo_lat_h         !geographical latitudes  of H-points
    ! Main arrays
    type(block2D_real8), dimension(:), pointer :: block_geo_lon_t,   &    !geographical longitudes of T-points
                                                  block_geo_lat_t,   &    !geographical latitudes  of T-points
                                                  block_geo_lon_u,   &    !geographical longitudes of U-points
                                                  block_geo_lat_u,   &    !geographical latitudes  of U-points
                                                  block_geo_lon_v,   &    !geographical longitudes of V-points
                                                  block_geo_lat_v,   &    !geographical latitudes  of V-points
                                                  block_geo_lon_h,   &    !geographical longitudes of H-points
                                                  block_geo_lat_h         !geographical latitudes  of H-points
    ! Temporary pointers
    real(8), dimension(:, :, :), pointer :: rotvec_coeff       !cos and sin of angles between coordinate lines
    ! Main arrays
    type(block3D_real8), dimension(:), pointer :: block_rotvec_coeff       !cos and sin of angles between coordinate lines

    ! First points of entire basin area
    real(8) :: xtm1(1), ytn1(1), xum1(1), yvn1(1)
    ! First points of local (regional) area for output
    real(8) :: xtm1loc(1), ytn1loc(1), xum1loc(1), yvn1loc(1)

contains

    ! Set all temporary pointers to null
    subroutine set_block_zero(k)
        implicit none
        integer :: k
        lu  =>null()
        lu1 => null()
        luu => null()
        luh => null()
        lcu => null()
        lcv => null()
        llu => null()
        llv => null()
        lu_output => null()
        lcu_output => null()
        lcv_output => null()
        llu_output => null()
        llv_output => null()
        hhh => null()
        hhhp => null()
        hhhn => null()
        hhq_rest => null()
        hhq => null()
        hhqp => null()
        hhqn => null()
        hhu => null()
        hhup => null()
        hhun => null()
        hhv => null()
        hhvp => null()
        hhvn => null()
        rlh_s => null()
        rlh_c => null()
        dxt => null()
        dyt => null()
        dx => null()
        dy => null()
        dxh => null()
        dyh => null()
        dxb => null()
        dyb => null()
        xt => null()
        yt => null()
        xu => null()
        yv => null()
        geo_lon_t => null()
        geo_lat_t => null()
        geo_lon_u => null()
        geo_lat_u => null()
        geo_lon_v => null()
        geo_lat_v => null()
        geo_lon_h => null()
        geo_lat_h => null()
        rotvec_coeff => null()
    end subroutine set_block_zero

    ! Set temporary pointers of masks to k-th block
    subroutine set_block_lu(k)
        implicit none
        integer :: k
        ! Set all null pointers, this is only for debugging
        call set_block_zero(k)
        lu  => block_lu(k)%vals
        lu1 => block_lu1(k)%vals
        luu => block_luu(k)%vals
        luh => block_luh(k)%vals
        lcu => block_lcu(k)%vals
        lcv => block_lcv(k)%vals
        llu => block_llu(k)%vals
        llv => block_llv(k)%vals
    end subroutine set_block_lu

    ! Set temporary pointers of regional masks to k-th block
    subroutine set_block_lu_output(k)
        implicit none
        integer :: k
        lu_output => block_lu_output(k)%vals
        lcu_output => block_lcu_output(k)%vals
        lcv_output => block_lcv_output(k)%vals
        llu_output => block_llu_output(k)%vals
        llv_output => block_llv_output(k)%vals
    end subroutine set_block_lu_output

    ! Set temporary pointers of heights to k-th block
    subroutine set_block_h(k)
        implicit none
        integer :: k
        hhh => block_hhh(k)%vals
        hhhp => block_hhhp(k)%vals
        hhhn => block_hhhn(k)%vals
        hhq_rest => block_hhq_rest(k)%vals
        hhq => block_hhq(k)%vals
        hhqp => block_hhqp(k)%vals
        hhqn => block_hhqn(k)%vals
        hhu => block_hhu(k)%vals
        hhup => block_hhup(k)%vals
        hhun => block_hhun(k)%vals
        hhv => block_hhv(k)%vals
        hhvp => block_hhvp(k)%vals
        hhvn => block_hhvn(k)%vals
    end subroutine set_block_h

    ! Set temporary pointers of coriolis to k-th block
    subroutine set_block_cor(k)
        implicit none
        integer :: k
        rlh_s => block_rlh_s(k)%vals
        rlh_c => block_rlh_c(k)%vals
    end subroutine set_block_cor

    ! Set temporary pointers of grid steps to k-th block
    subroutine set_block_dxdy(k)
        implicit none
        integer :: k
        dxt => block_dxt(k)%vals
        dyt => block_dyt(k)%vals
        dx => block_dx(k)%vals
        dy => block_dy(k)%vals
        dxh => block_dxh(k)%vals
        dyh => block_dyh(k)%vals
        dxb => block_dxb(k)%vals
        dyb => block_dyb(k)%vals
    end subroutine set_block_dxdy

    ! Set temporary pointers of grids to k-th block
    subroutine set_block_xy(k)
        implicit none
        integer :: k
        xt => block_xt(k)%vals
        yt => block_yt(k)%vals
        xu => block_xu(k)%vals
        yv => block_yv(k)%vals
    end subroutine set_block_xy

    ! Set temporary pointers of geo coordinates to k-th block
    subroutine set_block_geo(k)
        implicit none
        integer :: k
        geo_lon_t => block_geo_lon_t(k)%vals
        geo_lat_t => block_geo_lat_t(k)%vals
        geo_lon_u => block_geo_lon_u(k)%vals
        geo_lat_u => block_geo_lat_u(k)%vals
        geo_lon_v => block_geo_lon_v(k)%vals
        geo_lat_v => block_geo_lat_v(k)%vals
        geo_lon_h => block_geo_lon_h(k)%vals
        geo_lat_h => block_geo_lat_h(k)%vals
    end subroutine set_block_geo

    ! Set temporary pointers of rotvec to k-th block
    subroutine set_block_rotvec(k)
        implicit none
        integer :: k
        rotvec_coeff => block_rotvec_coeff(k)%vals
    end subroutine set_block_rotvec

    ! Get point in t-grid with number (m, n)
    subroutine get_tpoint(xtt, ytt, m, n)
        implicit none
        real*8, intent(out) :: xtt(1), ytt(1)
        integer, intent(in) :: m, n

        if (xgr_type==0) then
            xtt = rlon + dfloat(m-mmm)*dxst
        else !in case of irregular grid
            xtt = x_levels(m)
        endif
        if (ygr_type==0) then
            ytt = rlat + dfloat(n-nnn)*dyst
        else !in case of irregular grid
            ytt = y_levels(n)
        endif
    end subroutine

    ! Get point in uv-grid with number (m, n)
    subroutine get_uvpoint(xuu, yvv, m, n)
        implicit none
        real*8, intent(out) :: xuu(1), yvv(1)
        integer, intent(in) :: m, n
        real*8 :: xtt(1), ytt(1), xttp(1), yttp(1)

        call get_tpoint(xtt, ytt, m, n)
        call get_tpoint(xttp, yttp, m+1, n+1)
        xuu = (xtt + xttp)/2.0d0
        yvv = (ytt + yttp)/2.0d0
    end subroutine

endmodule basin_grid
!---------------------end module for definition of basin grid arrays-----------------
