module shallow_water
    use main_basin_pars
    use mpi_parallel_tools
    use mixing
    use vel_ssh
    use key_switches
    implicit none

contains

    subroutine compute_rhs(diffslpr, surf_stress, RHSx2d, RHSy2d, hhu, hhv, dxt, dyt, dxh, dyh, lu)
        implicit none

        real*8 :: diffslpr, surf_stress

        real*8 ::   RHSx2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                    RHSy2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                    hhu(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                    hhv(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                    dxt(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                    dyt(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                    dxh(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                    dyh(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                    lu(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

        integer :: m, n

        do n=ny_start,ny_end
            do m=nx_start,nx_end
                if (lu(m,n)>0.5) then
                    RHSx2d(m, n) = (surf_stress)*dxt(m,n)*dyh(m,n) -(diffslpr)*hhu(m,n)*dyh(m,n)/RefDen
                    RHSy2d(m, n) = (surf_stress)*dyt(m,n)*dxh(m,n) -(diffslpr)*hhv(m,n)*dxh(m,n)/RefDen
                endif
            enddo
        enddo
    end subroutine


    ! explicit shallow water equation sloving
    subroutine expl_shallow_water( tau,     &
                                  ksw4,     &
                                   ubrtr,   &
                                   ubrtrp,  &
                                   ubrtrn,  &
                                   vbrtr,   &
                                   vbrtrp,  &
                                   vbrtrn,  &
                                   ssh,     &
                                   sshp,    &
                                   sshn,    &
                                     RHSx,  &
                                     RHSy,  &
                                    wflux,  &
                                       mu,  &
                                      mu4,  &
                                     vort,  &
                                  str_t2d,  &
                                  str_s2d,  &
                                     rdis,  &
                                 RHSx_adv,  &
                                 RHSy_adv,  &
                                 RHSx_dif,  &
                                 RHSy_dif,  &
                                 RHSx_bfc,  &
                                 RHSy_bfc)

        !use vel_ssh
        use basin_grid
        use depth
        implicit none

        real(8) tau
        integer ksw4, k, m, n

        type(block2D), dimension(:), pointer :: ubrtr,      &
                                                ubrtrp,     &
                                                ubrtrn,     &
                                                 vbrtr,     &
                                                vbrtrp,     &
                                                vbrtrn,     &
                                                   ssh,     &
                                                  sshp,     &
                                                  sshn,     &
                                                   wflux,     &
                                                    RHSx,     &
                                                    RHSy,     &
                                                      mu,     &
                                                     mu4,     &
                                                    vort,     &
                                                 str_t2d,     &
                                                 str_s2d,     &
                                                    rdis,     &
                                                RHSx_adv,     &
                                                RHSy_adv,     &
                                                RHSx_dif,     &
                                                RHSy_dif,     &
                                                RHSx_bfc,     &
                                                RHSy_bfc

        real*8 bp, bp0, grx, gry, slx, sly, slxn, slyn

        real*8 time_count
        integer ierr

        ! Computing ssh
        !$omp parallel do
        do k = 1, bcount
            call set_block_boundary(k)
            call compute_ssh(tau, ubrtr(k)%vals, vbrtr(k)%vals, ssh(k)%vals, sshp(k)%vals, sshn(k)%vals ,  &
                             hhu(k)%vals, hhv(k)%vals, wflux(k)%vals,                                      &
                             dxh(k)%vals, dyh(k)%vals, dx(k)%vals, dy(k)%vals, lu(k)%vals)
        enddo
        !$omp end parallel do
        call syncborder_block2D(sshn)

        if (full_free_surface>0) then
            do k = 1, bcount
                call set_block_boundary(k)
                call hh_update(hhqn(k)%vals, hhun(k)%vals, hhvn(k)%vals, hhhn(k)%vals, sshn(k)%vals, hhq_rest(k)%vals, &
                               dx(k)%vals, dy(k)%vals, dxt(k)%vals, dyt(k)%vals,   &
                               dxh(k)%vals, dyh(k)%vals, dxb(k)%vals, dyb(k)%vals, &
                               lu(k)%vals, llu(k)%vals, llv(k)%vals, luh(k)%vals)
            enddo
            call syncborder_block2D(hhun)
            call syncborder_block2D(hhvn)
            call syncborder_block2D(hhhn)
        endif

        ! Computing advective and lateral-viscous terms for 2d-velocity
        if (trans_terms > 0) then
            !$omp parallel do
            do k = 1, bcount
                call set_block_boundary(k)
                call uv_vort(ubrtr(k)%vals, vbrtr(k)%vals, vort(k)%vals,        &
                             dxt(k)%vals, dyt(k)%vals, dxb(k)%vals, dyb(k)%vals, &
                             luu(k)%vals)
            enddo
            !$omp end parallel do
            call syncborder_block2D(vort)

            !$omp parallel do
            do k = 1, bcount
                call set_block_boundary(k)
                call uv_trans(ubrtr(k)%vals, vbrtr(k)%vals, vort(k)%vals,         &
                              hhq(k)%vals, hhu(k)%vals, hhv(k)%vals, hhh(k)%vals, &
                              RHSx_adv(k)%vals, RHSy_adv(k)%vals,                 &
                              dxh(k)%vals, dyh(k)%vals,                           &
                              lcu(k)%vals, lcv(k)%vals, luu(k)%vals)
            enddo
            !$omp end parallel do
        endif

        if (ksw_lat > 0) then
            !$omp parallel do
            do k = 1, bcount
                call set_block_boundary(k)
                call stress_components(ubrtrp(k)%vals, vbrtrp(k)%vals, str_t2d(k)%vals, str_s2d(k)%vals,  &
                                       dx(k)%vals, dy(k)%vals, dxt(k)%vals, dyt(k)%vals,   &
                                       dxh(k)%vals, dyh(k)%vals, dxb(k)%vals, dyb(k)%vals, &
                                       lu(k)%vals, luu(k)%vals)

            enddo
            !$omp end parallel do
            call syncborder_block2D(str_t2d)
            call syncborder_block2D(str_s2d)

            !$omp parallel do
            do k = 1, bcount
                call set_block_boundary(k)
                call uv_diff2(mu(k)%vals, str_t2d(k)%vals, str_s2d(k)%vals,        &
                              hhq(k)%vals, hhu(k)%vals, hhv(k)%vals, hhh(k)%vals,  &
                              RHSx_dif(k)%vals, RHSy_dif(k)%vals,                  &
                              dx(k)%vals, dy(k)%vals, dxt(k)%vals, dyt(k)%vals,    &
                              dxh(k)%vals, dyh(k)%vals, dxb(k)%vals, dyb(k)%vals,  &
                              lcu(k)%vals, lcv(k)%vals)

            enddo
            !$omp end parallel do
        endif

        ! compute BottomFriction (bfc)
        if (ksw_bfc > 0) then
            !$omp parallel do
            do k = 1, bcount
                call set_block_boundary(k)
                call uv_bfc(ubrtrp(k)%vals, vbrtrp(k)%vals, hhq(k)%vals, hhu(k)%vals, hhv(k)%vals, hhh(k)%vals,  &
                            RHSx_bfc(k)%vals, RHSy_bfc(k)%vals,                                                  &
                            dxb(k)%vals, dyb(k)%vals, lcu(k)%vals, lcv(k)%vals)
            enddo
            !$omp end parallel do
        endif

        ! Compute velocities
        !$omp parallel do
        do k = 1, bcount
            call set_block_boundary(k)
            call compute_vel(tau, ubrtr(k)%vals, ubrtrn(k)%vals, ubrtrp(k)%vals, vbrtr(k)%vals, vbrtrn(k)%vals, vbrtrp(k)%vals,  &
                                 ssh(k)%vals, hhu(k)%vals, hhun(k)%vals, hhup(k)%vals, hhv(k)%vals, hhvn(k)%vals, hhvp(k)%vals, hhh(k)%vals, &
                                 RHSx(k)%vals, RHSx_dif(k)%vals, RHSx_adv(k)%vals, RHSx_bfc(k)%vals, RHSy(k)%vals, RHSy_dif(k)%vals, RHSy_adv(k)%vals, RHSy_bfc(k)%vals, &
                                 dxt(k)%vals, dyt(k)%vals, dxh(k)%vals, dyh(k)%vals, dxb(k)%vals, dyb(k)%vals, rlh_s(k)%vals, &
                                 rdis(k)%vals, lcu(k)%vals, lcv(k)%vals)
        enddo
        !$omp end parallel do
        call syncborder_block2D(ubrtrn)
        call syncborder_block2D(vbrtrn)

        ! Shifting time indices
        !$omp parallel do
        do k = 1, bcount
            call set_block_boundary(k)
            call shift_var(ssh(k)%vals, sshp(k)%vals, sshn(k)%vals, lu(k)%vals)
            call shift_var(ubrtr(k)%vals, ubrtrp(k)%vals, ubrtrn(k)%vals, lcu(k)%vals)
            call shift_var(vbrtr(k)%vals, vbrtrp(k)%vals, vbrtrn(k)%vals, lcv(k)%vals)
        enddo
        !$omp end parallel do

        if (full_free_surface>0) then
            !$omp parallel do
            do k = 1, bcount
                call set_block_boundary(k)
                call hh_shift(hhq(k)%vals, hhqp(k)%vals, hhqn(k)%vals,   &
                              hhu(k)%vals, hhup(k)%vals, hhun(k)%vals,   &
                              hhv(k)%vals, hhvp(k)%vals, hhvn(k)%vals,   &
                              hhh(k)%vals, hhhp(k)%vals, hhhn(k)%vals,   &
                              lu(k)%vals, llu(k)%vals, llv(k)%vals, luh(k)%vals)

            enddo
            !$omp end parallel do
        endif

        if (full_free_surface>0) then
            ! Initialize depth for external mode
            !$omp parallel do
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
             !$omp end parallel do
             call syncborder_block2D(hhu)
             call syncborder_block2D(hhup)
             call syncborder_block2D(hhun)
             call syncborder_block2D(hhv)
             call syncborder_block2D(hhvp)
             call syncborder_block2D(hhvn)
             call syncborder_block2D(hhh)
             call syncborder_block2D(hhhp)
             call syncborder_block2D(hhhn)
        endif

    endsubroutine expl_shallow_water

    subroutine compute_vel(tau, ubrtr, ubrtrn, ubrtrp, vbrtr, vbrtrn, vbrtrp, ssh, hhu, hhun, hhup, hhv, hhvn, hhvp, hhh, &
                           RHSx, RHSx_dif, RHSx_adv, RHSx_bfc, RHSy, RHSy_dif, RHSy_adv, RHSy_bfc,  &
                           dxt, dyt, dxh, dyh, dxb, dyb, rlh_s, rdis, lcu, lcv)
        implicit none

        real(8)  ubrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                ubrtrp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                ubrtrn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                 vbrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                vbrtrp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                vbrtrn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                   ssh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                   hhu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                   hhun(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &
                   hhup(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &
                   hhv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                   hhvn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &
                   hhvp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &
                   hhh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                    RHSx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                    RHSy(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                    rdis(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                RHSx_adv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                RHSy_adv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                RHSx_dif(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                RHSy_dif(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                RHSx_bfc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                RHSy_bfc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                dxt(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                dyt(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                dxh(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                dyh(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                dxb(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                dyb(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                rlh_s(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                lcu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
                lcv(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

        real(8) tau, bp, bp0, grx, gry, slx, sly, slxn, slyn
        integer :: m, n

        do n=ny_start,ny_end
            do m=nx_start,nx_end
                !zonal flux
                if(lcu(m,n)>0.5) then
                    bp  = hhun(m,n)*dxt(m,n)*dyh(m,n)/2.0d0/tau
                    bp0 = hhup(m,n)*dxt(m,n)*dyh(m,n)/2.0d0/tau

                    slx = - FreeFallAcc * ( ssh(m+1,n) - ssh(m,n))*dyh(m,n)* hhu(m,n)
                    !slxn= - FreeFallAcc * (sshn(m+1,n) -sshn(m,n))*dyh(m,n)*hhun_e(m,n)
                    grx = RHSx(m,n) + slx + RHSx_dif(m,n) + RHSx_adv(m,n)      &
                        - (rdis(m,n)+rdis(m+1,n))/2.0d0 *ubrtrp(m,n)*dxt(m,n)*dyh(m,n)*hhu(m,n)        &
                        + ( rlh_s(m,n  )*hhh(m,n  )*dxb(m,n  )*dyb(m,n  )*(vbrtr(m+1,n  ) + vbrtr(m,n  ))          &
                        +   rlh_s(m,n-1)*hhh(m,n-1)*dxb(m,n-1)*dyb(m,n-1)*(vbrtr(m+1,n-1) + vbrtr(m,n-1)) )/4.0d0

                    ubrtrn(m,n) = (ubrtrp(m,n)*bp0 + grx )/(bp - RHSx_bfc(m, n))
                endif

                !meridional flux
                if(lcv(m,n)>0.5) then
                    bp  = hhvn(m,n)*dyt(m,n)*dxh(m,n)/2.0d0/tau
                    bp0 = hhvp(m,n)*dyt(m,n)*dxh(m,n)/2.0d0/tau

                    sly = - FreeFallAcc * ( ssh(m,n+1)- ssh(m,n))*dxh(m,n)* hhv(m,n)
                    !slyn= - FreeFallAcc * (sshn(m,n+1)-sshn(m,n))*dxh(m,n)*hhvn_e(m,n)
                    gry = RHSy(m,n) + sly  + RHSy_dif(m,n) + RHSy_adv(m,n)      &
                        - (rdis(m,n)+rdis(m,n+1))/2.0d0 *vbrtrp(m,n)*dxh(m,n)*dyt(m,n)*hhv(m,n)        &
                        - ( rlh_s(m  ,n)*hhh(m  ,n)*dxb(m  ,n)*dyb(m  ,n)*(ubrtr(m  ,n+1)+ubrtr(m  ,n))             &
                        +   rlh_s(m-1,n)*hhh(m-1,n)*dxb(m-1,n)*dyb(m-1,n)*(ubrtr(m-1,n+1)+ubrtr(m-1,n)) )/4.0d0

                    vbrtrn(m,n) = (vbrtrp(m,n)*bp0 + gry )/(bp - RHSy_bfc(m, n))
                endif
            enddo
        enddo
    end subroutine compute_vel

    subroutine compute_ssh(tau, ubrtr, vbrtr, ssh, sshp, sshn, hhu, hhv, wflux, dxh, dyh, dx, dy, lu)
        implicit none

        real(8)  tau

        real(8)  ubrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                 vbrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                   ssh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                   sshp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &
                   sshn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &
                   hhu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                   hhv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                   wflux(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                   dxh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                   dyh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                   dx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                   dy(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                   lu(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

        integer :: m, n

        do n=ny_start,ny_end
            do m=nx_start,nx_end
                if(lu(m,n)>0.5) then
                    sshn(m,n) = sshp(m,n) + 2.0d0*tau*( wflux(m,n)/RefDen*dfloat(full_free_surface)   &
                    - ( ubrtr(m,n)*hhu(m,n)*dyh(m,n) - ubrtr(m-1,n)*hhu(m-1,n)*dyh(m-1,n)             &
                      + vbrtr(m,n)*hhv(m,n)*dxh(m,n) - vbrtr(m,n-1)*hhv(m,n-1)*dxh(m,n-1) )/(dx(m,n)*dy(m,n))  )
                endif
            enddo
        enddo

    end subroutine

    subroutine shift_var(var, varp, varn, lvar)
        implicit none

        real*8 ::   var(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                    varp(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                    varn(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                    lvar(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

        integer :: m, n

        do n = ny_start-1,ny_end+1
            do m = nx_start-1,nx_end+1
                if (lvar(m,n)>0.5) then
                    varp(m,n) = var(m,n) + time_smooth*(varn(m,n) - 2.0d0*var(m,n) + varp(m,n))/2.0d0
                    var(m,n) = varn(m,n)
                endif
            enddo
        enddo

    end subroutine

endmodule
