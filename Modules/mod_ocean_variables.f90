!input physical parameters of the run
module key_switches
implicit none

integer ksw_atmforc,    &     !Atmospheric forcing (0 - no, 1 - yes)
        ksw_bfc,        &     !Bottom friction (0 - no, 1 - yes)
        ksw_lat,        &     !Lateral 2nd order mix parametrization (0 - no, 1 - yes)
        ksw_lat4,       &     !Lateral 4nd order momentum mix parametrization (0 - no, 1 - yes)
        ksw_ssbc,       &     !Type of surface boundary conditions (1 - surface T&S and wind stress are set; 2 - T&S fluxes and wind stress are set; 3 - T&S fluxes and wind stress are computed
        ksw_wflux,      &     !normalize global mean salt balance (0 - no, 1 - normalize water flux, 2 - normalize salinity flux)
        ksw_lbc_ts,     &     !open boundary conditions for T&S using (0 - no, 1 - yes)
        ksw_lbc_uv,     &     !open boundary conditions for U&V using (0 - no, 1 - yes)
        ksw_lbc_ssh           !open boundary conditions for SSH using (0 - no, 1 - yes)

real(8) lvisc_2,        &     !lateral  vicosity(2nd order)[m**2/s]
        lvisc_4,        &     !lateral  vicosity(4th order) [undim]
        nbfc                  !Bottom friction coeff (Manning's roughness)

endmodule key_switches

!-------------module for description common ogcm variables and task control parameters---------
module ocean_variables
    use key_switches
    use mpi_parallel_tools
    implicit none
    !barotropic dynamics arrays
    type(block2D_real8), dimension(:), pointer ::   ssh,     &  !sea surface height (SSH) at current  time step [m] (internal mode)
                                                  ubrtr,     &  !barotropic velocity      zonal[m/s] at current time step (internal mode)
                                                  vbrtr,     &  !barotropic velocity meridional[m/s] at current time step (internal mode)
                                                  RHSx2d,    &  !x-component of external force(barotropic)
                                                  RHSy2d        !y-component of external force(barotropic)

    type(block2D_real8), dimension(:), pointer ::  sshn,    &
                                                   sshp,    &  !sea surface height (SSH) at previous time step [m] (external mode)
                                                  ubrtrn,   &
                                                  ubrtrp,   &  !barotropic velocity      zonal[m/s] at previous time step (external mode)
                                                  vbrtrn,   &
                                                  vbrtrp       !barotropic velocity meridional[m/s] at previous time step (external mode)

    ! sea surface boundary condition
    type(block2D_real8), dimension(:), pointer :: wf_tot !total water flux

    type(block2D_real8), dimension(:), pointer :: BottomFriction,  &  !Bottom friction rate (m/s)
                                                  r_diss              !Rayleigh friction scale (1/s)

    type(block2D_real8), dimension(:), pointer :: amuv2d,               &    !depth mean lateral viscosity
                                                  amuv42d,              &    !depth mean lateral viscosity
                                                  r_vort2d,             &    !relative vorticity of depth mean velocity
                                                  stress_t2d,           &    !Horizontal tension tensor component (barotropic)
                                                  stress_s2d,           &    !Horizontal shearing tensor component(barotropic)
                                                  RHSx2d_tran_disp,     &    !dispersion x-component of external force(barotropic)
                                                  RHSy2d_tran_disp,     &    !dispersion y-component of external force(barotropic)
                                                  RHSx2d_diff_disp,     &    !dispersion x-component of external force(barotropic)
                                                  RHSy2d_diff_disp,     &    !dispersion y-component of external force(barotropic)
                                                  RHSx2d_bfc,           &
                                                  RHSy2d_bfc

real(8), allocatable :: ssh_max_amplitude(:, :),    &
                        ubrtr_max_amplitude(:, :),  &
                        vbrtr_max_amplitude(:, :)
                        
endmodule ocean_variables
!-------------end module for description common ogcm variables and task control parameters---------
