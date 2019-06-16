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
implicit none
!barotropic dynamics arrays
real(8),allocatable::   ssh(:,:),     &  !sea surface height (SSH) at current  time step [m] (internal mode)
                       pgrx(:,:),     &  !pressure gradient x-component for RHS
                       pgry(:,:),     &  !pressure gradient y-component for RHS
                      ubrtr(:,:),     &  !barotropic velocity      zonal[m/s] at current time step (internal mode)
                      vbrtr(:,:),     &  !barotropic velocity meridional[m/s] at current time step (internal mode)
                     RHSx2d(:,:),     &  !x-component of external force(barotropic)
                     RHSy2d(:,:)         !y-component of external force(barotropic)

real(8), allocatable::  sshn(:,:),   &
                        sshp(:,:),   &  !sea surface height (SSH) at previous time step [m] (external mode)
                      ubrtrn(:,:),   &
                      ubrtrp(:,:),   &  !barotropic velocity      zonal[m/s] at previous time step (external mode)
                      vbrtrn(:,:),   &
                      vbrtrp(:,:)       !barotropic velocity meridional[m/s] at previous time step (external mode)

!3d dynamics arrays
real(8),allocatable:: xxt(:,:,:),   &  !auxiliary array 1
                      yyt(:,:,:)       !auxiliary array 2

! sea surface boundary condition
real(8), allocatable:: tflux_surf(:,:),      &       !total surface heat flux [�C*m/s]
                       tflux_bot(:,:),       &       !total bottom heat flux [�C*m/s]
                       sflux_surf(:,:),      &       !total surface salt flux [psu*m/s]
                       sflux_bot(:,:),       &       !total bottom salt flux [psu*m/s]
                   surf_stress_x(:,:),       &       !wind      zonal stress per water density [m^2/s^2]
                   surf_stress_y(:,:),       &       !wind meridional stress per water density [m^2/s^2]
                    bot_stress_x(:,:),       &       !bottom    zonal stress per water density [m^2/s^2]
                    bot_stress_y(:,:),       &       !bottom meridional stress per water density [m^2/s^2]
                      divswrad(:,:,:),       &       !shortwave radiation divergence coefficients
                            dkft(:,:),       &       !relaxation coefficient for SST, [m/s]
                            dkfs(:,:),       &       !relaxation coefficient for SSS, [m/s]
                        sensheat(:,:),       &       !sensible heat flux
                         latheat(:,:),       &       !latent heat flux
                          lw_bal(:,:),       &       !longwave radiation balance
                          sw_bal(:,:),       &       !shortwave radiation balance
                          hf_tot(:,:),       &       !total heat flux
                         wf_tot(:,:)                 !total water flux

!Atmospheric arrays for bulk-formulae
real(8),allocatable:: tatm(:,:),   &    !Air temperature, [�C]
                      qatm(:,:),   &    !Air humidity, [kg/kg]
                      rain(:,:),   &    !rain, [kg/m^2/s]
                      snow(:,:),   &    !snow, [kg/m^2/s]
                      wind(:,:),   &    !Wind speed module, [m/s]
                       lwr(:,:),   &    !Downward  longwave radiation, [W/m^2]
                       swr(:,:),   &    !Downward shortwave radiation, [W/m^2]
                      slpr(:,:),   &    !Sea level pressure, [Pa]
                      uwnd(:,:),   &    !Zonal      wind speed, [m/s]
                      vwnd(:,:),   &    !Meridional wind speed, [m/s]
                      taux(:,:),   &    !Zonal      wind stress, [Pa]
                      tauy(:,:)         !Meridional wind stress, [Pa]

real(8), allocatable:: BottomFriction(:,:),    &    !Bottom friction rate (m/s)
                               r_diss(:,:)          !Rayleigh friction scale (1/s)

real(8), allocatable:: amuv2d(:,:),     &    !depth mean lateral viscosity
                      amuv42d(:,:),     &    !depth mean lateral viscosity
                     r_vort2d(:,:),     &    !relative vorticity of depth mean velocity
                   stress_t2d(:,:),     &    !Horizontal tension tensor component (barotropic)
                   stress_s2d(:,:),     &    !Horizontal shearing tensor component(barotropic)
             RHSx2d_tran_disp(:,:),     &    !dispersion x-component of external force(barotropic)
             RHSy2d_tran_disp(:,:),     &    !dispersion y-component of external force(barotropic)
             RHSx2d_diff_disp(:,:),     &    !dispersion x-component of external force(barotropic)
             RHSy2d_diff_disp(:,:),     &    !dispersion y-component of external force(barotropic)
             RHSx2d_bfc(:, :),          &
             RHSy2d_bfc(:, :)

real(8), allocatable :: ssh_max_amplitude(:, :),    &
                        ubrtr_max_amplitude(:, :),  &
                        vbrtr_max_amplitude(:, :)
                        
endmodule ocean_variables
!-------------end module for description common ogcm variables and task control parameters---------
