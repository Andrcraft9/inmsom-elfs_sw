module shallow_water
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid

    use key_switches

    use ocalg_routes
    use depth_routes
    use mixing_routes
    use flux_routes
    use velssh_routes
    implicit none

contains
    ! Compute RHS which includes wind forcing
    subroutine compute_fluxes_rhs(surf_stress_x, surf_stress_y, slpr, RHSx2d, RHSy2d)
        implicit none
        real(8) surf_stress_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                surf_stress_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                slpr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),              &
                RHSx2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),            &
                RHSy2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

        do n=ny_start,ny_end
            do m=nx_start,nx_end
                if(lu(m,n)>0.5) then
                    !RHSx2d(m, n) = ( surf_stress_x(m,n)+bot_stress_x(m,n) )*dxt(m,n)*dyh(m,n)    &
                    RHSx2d(m, n) = (surf_stress_x(m,n))*dxt(m,n)*dyh(m,n)    &
                             -(slpr(m+1,n)-slpr(m,n))*hhu(m,n)*dyh(m,n)/RefDen

                    !RHSy2d(m, n) = ( surf_stress_y(m,n)+bot_stress_y(m,n) )*dyt(m,n)*dxh(m,n)    &
                    RHSy2d(m, n) = (surf_stress_y(m,n))*dyt(m,n)*dxh(m,n)    &
                             -(slpr(m,n+1)-slpr(m,n))*hhv(m,n)*dxh(m,n)/RefDen
                endif
            enddo
        enddo
    end subroutine

    ! Explicit shallow water equation sloving, leap-frog scheme for time
    subroutine expl_shallow_water( tau,     &
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
                                       fx,  &
                                       fy,  &
                                     rdis,  &
                                 RHSx_adv,  &
                                 RHSy_adv,  &
                                 RHSx_dif,  &
                                 RHSy_dif,  &
                                 RHSx_bfc,  &
                                 RHSy_bfc)

        implicit none
        real(8) tau
        integer m, n
        type(block2D_real8), dimension(:), pointer :: ubrtr,        &
                                                      ubrtrp,       &
                                                      ubrtrn,       &
                                                      vbrtr,        &
                                                      vbrtrp,       &
                                                      vbrtrn,       &
                                                      ssh,          &
                                                      sshp,         &
                                                      sshn,         &
                                                      wflux,        &
                                                      RHSx,         &
                                                      RHSy,         &
                                                      mu,           &
                                                      mu4,          &
                                                      vort,         &
                                                      str_t2d,      &
                                                      str_s2d,      &
                                                      fx,           &
                                                      fy,           &
                                                      rdis,         &
                                                      RHSx_adv,     &
                                                      RHSy_adv,     &
                                                      RHSx_dif,     &
                                                      RHSy_dif,     &
                                                      RHSx_bfc,     &
                                                      RHSy_bfc
        real(8) bp, bp0, grx, gry, slx, sly, slxn, slyn
        integer ierr

        ! Computing sea level
        do k = 1, bcount
            call set_block(k)
            call set_block_lu(k)
            call set_block_h(k)
            call set_block_dxdy(k)
            call compute_ssh(tau, ubrtr(k)%vals, vbrtr(k)%vals, ssh(k)%vals, sshp(k)%vals, sshn(k)%vals, wflux(k)%vals)
        enddo
        call syncborder_block2D_real8(sshn)
        if(periodicity_x/=0) then
            call cyclize_x_block2D_real8(sshn)
        endif
        if(periodicity_y/=0) then
            call cyclize_y_block2D_real8(sshn)
        endif

        ! Update depths by new sea level
        if (full_free_surface>0) then
            do k = 1, bcount
                call set_block(k)
                call set_block_lu(k)
                call set_block_h(k)
                call set_block_dxdy(k)
                call hh_update(hhqn, hhun, hhvn, hhhn, sshn, hhq_rest)
            enddo
            call syncborder_block2D_real8(block_hhun)
            call syncborder_block2D_real8(block_hhvn)
            call syncborder_block2D_real8(block_hhhn)
            if(periodicity_x/=0) then
                call cyclize_x_block2D_real8(block_hhun)
                call cyclize_x_block2D_real8(block_hhvn)
                call cyclize_x_block2D_real8(block_hhhn)
            endif
            if(periodicity_y/=0) then
                call cyclize_y_block2D_real8(block_hhun)
                call cyclize_y_block2D_real8(block_hhvn)
                call cyclize_y_block2D_real8(block_hhhn)
            endif
        endif

        ! Computing advective and lateral-viscous terms for 2d-velocity
        if (trans_terms > 0) then
            do k = 1, bcount
                call set_block(k)
                call set_block_lu(k)
                call set_block_dxdy(k)
                call uv_vort(ubrtr(k)%vals, vbrtr(k)%vals, vort(k)%vals)
            enddo
            call syncborder_block2D_real8(vort)
            if(periodicity_x/=0) then
                call cyclize_x_block2D_real8(vort)
            endif
            if(periodicity_y/=0) then
                call cyclize_y_block2D_real8(vort)
            endif

            do k = 1, bcount
                call set_block(k)
                call set_block_lu(k)
                call set_block_h(k)
                call set_block_dxdy(k)
                call uv_trans(ubrtr(k)%vals, vbrtr(k)%vals, vort(k)%vals,    &
                              hhq, hhu, hhv, hhh,                            &
                              RHSx_adv(k)%vals, RHSy_adv(k)%vals)
            enddo
        endif

        ! Computing diffusion terms
        if (ksw_lat > 0) then
            do k = 1, bcount
                call set_block(k)
                call set_block_lu(k)
                call set_block_dxdy(k)
                call stress_components(ubrtrp(k)%vals, vbrtrp(k)%vals, str_t2d(k)%vals, str_s2d(k)%vals)

            enddo
            call syncborder_block2D_real8(str_t2d)
            call syncborder_block2D_real8(str_s2d)
            if(periodicity_x/=0) then
                call cyclize_x_block2D_real8(str_t2d)
                call cyclize_x_block2D_real8(str_s2d)
            endif
            if(periodicity_y/=0) then
                call cyclize_y_block2D_real8(str_t2d)
                call cyclize_y_block2D_real8(str_s2d)
            endif

            do k = 1, bcount
                call set_block(k)
                call set_block_lu(k)
                call set_block_h(k)
                call set_block_dxdy(k)
                call uv_diff2(mu(k)%vals, str_t2d(k)%vals, str_s2d(k)%vals,        &
                              hhq, hhu, hhv, hhh,  &
                              RHSx_dif(k)%vals, RHSy_dif(k)%vals)

            enddo

            if(ksw_lat4 > 0) then
                call uv_diff4( mu4, str_t2d, str_s2d,  &
                               fx, fy, hhq, hhu, hhv, hhh,    &
                               RHSx_dif, RHSy_dif, 1 )
            endif
        endif

        ! Compute bottom friction
        if (ksw_bfc > 0) then
            do k = 1, bcount
                call set_block(k)
                call set_block_lu(k)
                call set_block_h(k)
                call set_block_dxdy(k)
                call uv_bfc(ubrtrp(k)%vals, vbrtrp(k)%vals,   &
                            hhq, hhu, hhv, hhh,               &
                            RHSx_bfc(k)%vals, RHSy_bfc(k)%vals)
            enddo
        endif

        ! Compute velocities
        do k = 1, bcount
            call set_block_lu(k)
            call set_block_h(k)
            call set_block_dxdy(k)
            call set_block_cor(k)
            call compute_vel(tau, ubrtr(k)%vals, ubrtrn(k)%vals, ubrtrp(k)%vals,      &
                                  vbrtr(k)%vals, vbrtrn(k)%vals, vbrtrp(k)%vals,      &
                                  ssh(k)%vals, hhu, hhun, hhup, hhv, hhvn, hhvp, hhh, &
                                  RHSx(k)%vals, RHSx_dif(k)%vals, RHSx_adv(k)%vals, RHSx_bfc(k)%vals,  &
                                  RHSy(k)%vals, RHSy_dif(k)%vals, RHSy_adv(k)%vals, RHSy_bfc(k)%vals,  &
                                  rdis(k)%vals)
        enddo
        call syncborder_block2D_real8(ubrtrn)
        call syncborder_block2D_real8(vbrtrn)
        if(periodicity_x/=0) then
            call cyclize_x_block2D_real8(ubrtrn)
            call cyclize_x_block2D_real8(vbrtrn)
        endif
        if(periodicity_y/=0) then
            call cyclize_y_block2D_real8(ubrtrn)
            call cyclize_y_block2D_real8(vbrtrn)
        endif

        ! Transition to the nex time step for sea level and velocity variables.
        ! Leap-frog scheme with filtration
        do k = 1, bcount
            call set_block(k)
            call set_block_lu(k)
            call shift_var(ssh(k)%vals, sshp(k)%vals, sshn(k)%vals, lu)
            call shift_var(ubrtr(k)%vals, ubrtrp(k)%vals, ubrtrn(k)%vals, lcu)
            call shift_var(vbrtr(k)%vals, vbrtrp(k)%vals, vbrtrn(k)%vals, lcv)
        enddo

        ! Transition to the nex time step for depths variables
        if(full_free_surface>0) then
            do k = 1, bcount
                call set_block(k)
                call set_block_lu(k)
                call set_block_h(k)
                call hh_shift(hhq, hhqp, hhqn,   &
                              hhu, hhup, hhun,   &
                              hhv, hhvp, hhvn,   &
                              hhh, hhhp, hhhn)
          enddo
        endif

        ! Initialize depth for external mode
        if(full_free_surface>0) then
            do k = 1, bcount
                call set_block(k)
                call set_block_lu(k)
                call set_block_h(k)
                call set_block_dxdy(k)
                call hh_init(hhq, hhqp, hhqn,    &
                             hhu, hhup, hhun,    &
                             hhv, hhvp, hhvn,    &
                             hhh, hhhp, hhhn,    &
                             ssh, sshp, hhq_rest)
            enddo
            call syncborder_block2D_real8(block_hhu)
            call syncborder_block2D_real8(block_hhup)
            call syncborder_block2D_real8(block_hhun)
            call syncborder_block2D_real8(block_hhv)
            call syncborder_block2D_real8(block_hhvp)
            call syncborder_block2D_real8(block_hhvn)
            call syncborder_block2D_real8(block_hhh)
            call syncborder_block2D_real8(block_hhhp)
            call syncborder_block2D_real8(block_hhhn)
            if(periodicity_x/=0) then
                call cyclize_x_block2D_real8(block_hhu )
                call cyclize_x_block2D_real8(block_hhup)
                call cyclize_x_block2D_real8(block_hhun)
                call cyclize_x_block2D_real8(block_hhv )
                call cyclize_x_block2D_real8(block_hhvp)
                call cyclize_x_block2D_real8(block_hhvn)
                call cyclize_x_block2D_real8(block_hhh )
                call cyclize_x_block2D_real8(block_hhhp)
                call cyclize_x_block2D_real8(block_hhhn)
            end if
            if(periodicity_y/=0) then
                call cyclize_y_block2D_real8(block_hhu )
                call cyclize_y_block2D_real8(block_hhup)
                call cyclize_y_block2D_real8(block_hhun)
                call cyclize_y_block2D_real8(block_hhv )
                call cyclize_y_block2D_real8(block_hhvp)
                call cyclize_y_block2D_real8(block_hhvn)
                call cyclize_y_block2D_real8(block_hhh )
                call cyclize_y_block2D_real8(block_hhhp)
                call cyclize_y_block2D_real8(block_hhhn)
            end if
        endif
    endsubroutine expl_shallow_water

    subroutine compute_vel(tau, ubrtr, ubrtrn, ubrtrp, vbrtr, vbrtrn, vbrtrp,  &
                           ssh, hhu, hhun, hhup, hhv, hhvn, hhvp, hhh,         &
                           RHSx, RHSx_dif, RHSx_adv, RHSx_bfc, RHSy, RHSy_dif, RHSy_adv, RHSy_bfc,  &
                           rdis)
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
                RHSy_bfc(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
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

    subroutine compute_ssh(tau, ubrtr, vbrtr, ssh, sshp, sshn, hhu, hhv, wflux)
        implicit none
        real(8)  tau
        real(8)  ubrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                 vbrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                   ssh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                   sshp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &
                   sshn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &
                   hhu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                   hhv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
                   wflux(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
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
        real*8 :: var(bnd_x1:bnd_x2,bnd_y1:bnd_y2), &
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
