module shallow_water
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid

    !use vel_ssh
    use depth

    implicit none

    contains

!===============================================================================
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

        !computing ssh
        !$omp parallel do
        do k = 1, bcount
            call set_block_boundary(k)
            do n=ny_start,ny_end
                do m=nx_start,nx_end
                    if (lu(k)%vals(m,n)>0.5) then
                        sshn(k)%vals(m,n) = sshp(k)%vals(m,n) + 2.0d0*tau * ( wflux(k)%vals(m,n)/RefDen*dfloat(full_free_surface)   &
                        -( ubrtr(k)%vals(m,n)*hhu(k)%vals(m,n)*dyh(k)%vals(m,n) - ubrtr(k)%vals(m-1,n)*hhu(k)%vals(m-1,n)*dyh(k)%vals(m-1,n)   &
                         + vbrtr(k)%vals(m,n)*hhv(k)%vals(m,n)*dxh(k)%vals(m,n) - vbrtr(k)%vals(m,n-1)*hhv(k)%vals(m,n-1)*dxh(k)%vals(m,n-1) ) &
                          /(dx(k)%vals(m,n)*dy(k)%vals(m,n)) )
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do

        call syncborder_block2D(sshn)

        if(full_free_surface>0) then
            call hh_update(hhqn, hhun, hhvn, hhhn, sshn, hhq_rest)
        endif

        !computing advective and lateral-viscous terms for 2d-velocity
        !call stress_components(ubrtrp, vbrtrp, str_t2d,str_s2d,1)

        !computing advective and lateral-viscous terms for 2d-velocity
        !call uv_trans(ubrtr, vbrtr, vort,     &
        !              hhq, hhu, hhv, hhh,     &
        !              RHSx_adv, RHSy_adv, 1)

        !call uv_diff2( mu, str_t2d, str_s2d,  &
        !               hhq, hhu, hhv, hhh,     &
        !               RHSx_dif, RHSy_dif, 1  )

        ! compute BottomFriction (bfc)
        !call uv_bfc(ubrtrp, vbrtrp, hhq, hhu, hhv, hhh, RHSx_bfc, RHSy_bfc)

        !$omp parallel do
        do k = 1, bcount
            call set_block_boundary(k)
            do n=ny_start,ny_end
                do m=nx_start,nx_end
                    !zonal flux
                    if (lcu(k)%vals(m,n) > 0.5) then
                        bp  = hhun(k)%vals(m,n)*dxt(k)%vals(m,n)*dyh(k)%vals(m,n)/2.0d0/tau
                        bp0 = hhup(k)%vals(m,n)*dxt(k)%vals(m,n)*dyh(k)%vals(m,n)/2.0d0/tau

                        slx = -FreeFallAcc * (ssh(k)%vals(m+1,n) - ssh(k)%vals(m,n))*dyh(k)%vals(m,n)* hhu(k)%vals(m,n)
                        !slxn= - FreeFallAcc * (sshn(m+1,n) -sshn(m,n))*dyh(m,n)*hhun_e(m,n)
                        grx = RHSx(k)%vals(m,n) + slx + RHSx_dif(k)%vals(m,n) + RHSx_adv(k)%vals(m,n) &
                            - (rdis(k)%vals(m,n)+rdis(k)%vals(m+1,n))/2.0d0 * ubrtrp(k)%vals(m,n)*dxt(k)%vals(m,n)*dyh(k)%vals(m,n)*hhu(k)%vals(m,n)  &
                            + ( rlh_s(k)%vals(m,n  )*hhh(k)%vals(m,n  )*dxb(k)%vals(m,n  )*dyb(k)%vals(m,n  )*(vbrtr(k)%vals(m+1,n  ) + vbrtr(k)%vals(m,n  ))   &
                            +   rlh_s(k)%vals(m,n-1)*hhh(k)%vals(m,n-1)*dxb(k)%vals(m,n-1)*dyb(k)%vals(m,n-1)*(vbrtr(k)%vals(m+1,n-1) + vbrtr(k)%vals(m,n-1)) )/4.0d0

                        ubrtrn(k)%vals(m,n) = (ubrtrp(k)%vals(m,n)*bp0 + grx )/(bp - RHSx_bfc(k)%vals(m, n))
                    endif

                    !meridional flux
                    if (lcv(k)%vals(m,n) > 0.5) then
                        bp  = hhvn(k)%vals(m,n)*dyt(k)%vals(m,n)*dxh(k)%vals(m,n)/2.0d0/tau
                        bp0 = hhvp(k)%vals(m,n)*dyt(k)%vals(m,n)*dxh(k)%vals(m,n)/2.0d0/tau

                        sly = -FreeFallAcc * (ssh(k)%vals(m,n+1) - ssh(k)%vals(m,n))*dxh(k)%vals(m,n)* hhv(k)%vals(m,n)
                        !slyn= - FreeFallAcc * (sshn(m,n+1)-sshn(m,n))*dxh(m,n)*hhvn_e(m,n)
                        gry = RHSy(k)%vals(m,n) + sly + RHSy_dif(k)%vals(m,n) + RHSy_adv(k)%vals(m,n)      &
                            - (rdis(k)%vals(m,n)+rdis(k)%vals(m,n+1))/2.0d0 * vbrtrp(k)%vals(m,n)*dxh(k)%vals(m,n)*dyt(k)%vals(m,n)*hhv(k)%vals(m,n)  &
                            - ( rlh_s(k)%vals(m  ,n)*hhh(k)%vals(m  ,n)*dxb(k)%vals(m  ,n)*dyb(k)%vals(m  ,n)*(ubrtr(k)%vals(m  ,n+1) + ubrtr(k)%vals(m  ,n))  &
                            +   rlh_s(k)%vals(m-1,n)*hhh(k)%vals(m-1,n)*dxb(k)%vals(m-1,n)*dyb(k)%vals(m-1,n)*(ubrtr(k)%vals(m-1,n+1)+ubrtr(k)%vals(m-1,n)) )/4.0d0

                        vbrtrn(k)%vals(m,n) = (vbrtrp(k)%vals(m,n)*bp0 + gry )/(bp - RHSy_bfc(k)%vals(m, n))
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do

        call syncborder_block2D(ubrtrn)
        call syncborder_block2D(vbrtrn)

        !shifting time indices
        !$omp parallel do
        do k = 1, bcount
            call set_block_boundary(k)
            do n=ny_start-1,ny_end+1
                do m=nx_start-1,nx_end+1
                    if (lu(k)%vals(m,n)>0.5) then
                        sshp(k)%vals(m,n) = ssh(k)%vals(m,n) &
                            +time_smooth*(sshn(k)%vals(m,n) - 2.0d0*ssh(k)%vals(m,n) + sshp(k)%vals(m,n))/2.0d0

                        ssh(k)%vals(m,n) = sshn(k)%vals(m,n)
                    endif
                    if (lcu(k)%vals(m,n)>0.5) then
                        !up(m,n) =  hhu_e(m,n)*u(m,n)+time_smooth*(hhun_e(m,n)*un(m,n)-2.0d0*hhu_e(m,n)*u(m,n)+hhup_e(m,n)*up(m,n))/2.0d0/dfloat(nstep)
                        ubrtrp(k)%vals(m,n) = ubrtr(k)%vals(m,n) &
                            + time_smooth*(ubrtrn(k)%vals(m,n) - 2.0d0*ubrtr(k)%vals(m,n) + ubrtrp(k)%vals(m,n))/2.0d0

                        ubrtr(k)%vals(m,n) = ubrtrn(k)%vals(m,n)
                    endif
                    if (lcv(k)%vals(m,n)>0.5) then
                        !vp(m,n) =  hhv_e(m,n)*v(m,n)+time_smooth*(hhvn_e(m,n)*vn(m,n)-2.0d0*hhv_e(m,n)*v(m,n)+hhvp_e(m,n)*vp(m,n))/2.0d0/dfloat(nstep)
                        vbrtrp(k)%vals(m,n) = vbrtr(k)%vals(m,n) &
                            + time_smooth*(vbrtrn(k)%vals(m,n) - 2.0d0*vbrtr(k)%vals(m,n) + vbrtrp(k)%vals(m,n))/2.0d0

                        vbrtr(k)%vals(m,n) = vbrtrn(k)%vals(m,n)
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do

        if(full_free_surface>0) then
            call hh_shift(hhq, hhqp, hhqn,   &
                          hhu, hhup, hhun,   &
                          hhv, hhvp, hhvn,   &
                          hhh, hhhp, hhhn)
        endif

        if(full_free_surface>0) then
            !initialize depth for external mode
            call hh_init(hhq, hhqp, hhqn,    &
                         hhu, hhup, hhun,    &
                         hhv, hhvp, hhvn,    &
                         hhh, hhhp, hhhn,    &
                         ssh, sshp, hhq_rest)
        endif

    endsubroutine expl_shallow_water

    subroutine compute_vel(tau, ubrtr, ubrtrn, ubrtrp, vbrtr, vbrtrn, vbrtrp, ssh, hhu, hhun, hhup, hhv, hhvn, hhvp, hhh, &
                           RHSx, RHSx_dif, RHSx_adv, RHSx_bfc, RHSy, RHSy_dif, RHSy_adv, RHSy_bfc,  &
                           dxt, dyt, dxh, dyh, dxb, dyb, rlh_s, rdis, lcu, lcv)
        use mpi_parallel_tools
        implicit none

        real(8)  ubrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                ubrtrp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                ubrtrn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                 vbrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                vbrtrp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                vbrtrn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                   ssh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                   hhu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &
                   hhun(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                   hhup(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                   hhv(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                   hhvn(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                   hhvp(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                   hhh(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
                    RHSx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                    RHSy(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                    rdis(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                RHSx_adv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                RHSy_adv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                RHSx_dif(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                RHSy_dif(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                RHSx_bfc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                RHSy_bfc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),  &
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

        !$omp parallel do private(bp, bp0, grx, gry, slx, sly, slxn, slyn)
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
!    call compute_vel(tau,   &
!                     ubrtr(k)%vals,  &
!                     ubrtrn(k)%vals,  &
!                     ubrtrp(k)%vals,  &
!                     vbrtr(k)%vals,  &
!                     vbrtrn(k)%vals,  &
!                     vbrtrp(k)%vals,  &
!                     ssh(k)%vals,  &
!                     hhu(k)%vals,  &
!                     hhun(k)%vals,     &
!                     hhup(k)%vals,     &
!                     hhv(k)%vals,         &
!                     hhvn(k)%vals,     &
!                     hhvp(k)%vals,     &
!                     hhh(k)%vals,     &
!                     RHSx(k)%vals,     &
!                     RHSx_dif(k)%vals,     &
!                     RHSx_adv(k)%vals,     &
!                     RHSx_bfc(k)%vals,     &
!                     RHSy(k)%vals,     &
!                     RHSy_dif(k)%vals,     &
!                     RHSy_adv(k)%vals,     &
!                     RHSy_bfc(k)%vals,     &
!                     dxt(k)%vals,     &
!                     dyt(k)%vals,     &
!                     dxh(k)%vals,     &
!                     dyh(k)%vals,     &
!                     dxb(k)%vals,     &
!                     dyb(k)%vals,     &
!                     rlh_s(k)%vals,     &
!                     rdis(k)%vals,     &
!                     lcu(k)%vals,     &
!                     lcv(k)%vals)

endmodule
