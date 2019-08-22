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
    type(block2D_real8), dimension(:), pointer ::   ssh,     &  !sea surface height (SSH) at current  time step [m]
                                                  ubrtr,     &  !barotropic velocity      zonal[m/s] at current time step
                                                  vbrtr,     &  !barotropic velocity meridional[m/s] at current time step
                                                  xxt,       &  !temporary arrays
                                                  yyt,       &  !temporary arrays
                                                  RHSx2d,    &  !x-component of external force(barotropic)
                                                  RHSy2d        !y-component of external force(barotropic)

    type(block2D_real8), dimension(:), pointer ::  sshn,     &
                                                   sshp,     &  !sea surface height (SSH) at previous time step [m]
                                                   ubrtrn,   &
                                                   ubrtrp,   &  !barotropic velocity      zonal[m/s] at previous time step
                                                   vbrtrn,   &
                                                   vbrtrp       !barotropic velocity meridional[m/s] at previous time step

    ! sea surface boundary condition
    type(block2D_real4), dimension(:), pointer :: tflux_surf,      &       !total surface heat flux [�C*m/s]
                                                  tflux_bot,       &       !total bottom heat flux [�C*m/s]
                                                  sflux_surf,      &       !total surface salt flux [psu*m/s]
                                                  sflux_bot,       &       !total bottom salt flux [psu*m/s]
                                              surf_stress_x,       &       !wind      zonal stress per water density [m^2/s^2]
                                              surf_stress_y,       &       !wind meridional stress per water density [m^2/s^2]
                                               bot_stress_x,       &       !bottom    zonal stress per water density [m^2/s^2]
                                               bot_stress_y,       &       !bottom meridional stress per water density [m^2/s^2]
                                                       dkft,       &       !relaxation coefficient for SST, [m/s]
                                                       dkfs,       &       !relaxation coefficient for SSS, [m/s]
                                                   sensheat,       &       !sensible heat flux
                                                    latheat,       &       !latent heat flux
                                                     lw_bal,       &       !longwave radiation balance
                                                     sw_bal,       &       !shortwave radiation balance
                                                     hf_tot,       &       !total heat flux
                                                     wf_tot                !total water flux

    ! Atmospheric arrays for bulk-formulae
    type(block2D_real4), dimension(:), pointer :: tatm,   &    !Air temperature, [�C]
                                                  qatm,   &    !Air humidity, [kg/kg]
                                                  rain,   &    !rain, [kg/m^2/s]
                                                  snow,   &    !snow, [kg/m^2/s]
                                                  evap,   &    !evaporation, [kg/m^2/s]
                                                  wind,   &    !Wind speed module, [m/s]
                                                   lwr,   &    !Downward  longwave radiation, [W/m^2]
                                                   swr,   &    !Downward shortwave radiation, [W/m^2]
                                                  slpr,   &    !Sea level pressure, [Pa]
                                                  uwnd,   &    !Zonal      wind speed, [m/s]
                                                  vwnd,   &    !Meridional wind speed, [m/s]
                                                  taux,   &    !Zonal      wind stress, [Pa]
                                                  tauy         !Meridional wind stress, [Pa]

    type(block2D_real8), dimension(:), pointer :: r_diss       !Rayleigh friction scale (1/s)

    type(block2D_real8), dimension(:), pointer :: amuv2d,               &    !depth mean lateral viscosity
                                                  amuv42d,              &    !depth mean lateral viscosity
                                                  r_vort2d,             &    !relative vorticity of depth mean velocity
                                                  stress_t2d,           &    !Horizontal tension tensor component (barotropic)
                                                  stress_s2d,           &    !Horizontal shearing tensor component(barotropic)
                                                  RHSx2d_tran_disp,     &    !dispersion x-component of external force(barotropic)
                                                  RHSy2d_tran_disp,     &    !dispersion y-component of external force(barotropic)
                                                  RHSx2d_diff_disp,     &    !dispersion x-component of external force(barotropic)
                                                  RHSy2d_diff_disp,     &    !dispersion y-component of external force(barotropic)
                                                  RHSx2d_bfc,           &    !bottom firction x-component
                                                  RHSy2d_bfc                 !bottom firction y-component

    type(block2D_real8), dimension(:), pointer  :: ssh_max_amplitude,    &   ! Maximum amplitute of sea surface height
                                                   ubrtr_max_amplitude,  &   ! Maximum amplitute of zonal velocity
                                                   vbrtr_max_amplitude       ! Maximum amplitute of meridional velocity

    type(block2D_real4), dimension(:), pointer :: array4_2d ! Temporary real4 array
    !type(block3D_real4), dimension(:), pointer :: array4_3d ! Temporary real4 array

endmodule ocean_variables
!-------------end module for description common ogcm variables and task control parameters---------
