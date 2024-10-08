program main_HDS

    ! A modified version (v2) of the Hysteretic Depressional Storage (HDS) model to simulate prairie fill and spill mechanism.
    ! The original HDS (v1) was based on the equations listed in Ahmed et al. (2023), doi: https://doi.org/10.1016/j.envsoft.2023.105769.
    ! This version of HDS (v2) is revised based on the equations presented by Clark and Shook (2022), doi: https://doi.org/10.1029/2022WR032694,
    ! which is cleaner and more robust than the one proposed in HDS (v1)

    USE type_HDS
    USE ascii_util, only: read_csv

    implicit none

    ! filenames
    character(len=100)               :: fName_spatial        ! file name that contains depressional storage information
    character(len=100)               :: fName_forcing        ! file name that contains time series inputs
    ! information in the spatial csv file
    integer(i4b),       parameter    :: ixCatArea=1,ixDepArea=2,ixDepVol=3, ixDepCatFrac= 4, ixDep_b=5, ixDep_p=6    ! index of variables in data structures
    character(len=128), parameter    :: s1='total_catchment_m2', s2='depression_area_m2', s3='depression_volume_m3', s4='depressionCatchAreaFrac', s5='depression_b', s6='depression_p' ! variable name (header) in the file
    ! information in the forcing csv file
    integer(i4b),       parameter    :: ixP=1,ixET=2,ixQ=3                    ! index of variables in data structures
    character(len=128), parameter    :: f1='pRate', f2='etPond', f3='qSeas'   ! NOTE: exclude quotes
    ! data structures
    type(hStruct), allocatable       :: spatialData(:)       ! spatial information
    type(hStruct), allocatable       :: forcingData(:)       ! forcing information
    integer(i4b)                     :: nBasins              ! number of basins
    integer(i4b)                     :: nTime                ! number of time steps
    integer(i4b)                     :: iVar                 ! variable index
    ! error control
    integer(i4b)                     :: ierr                 ! error code
    character(len=256)               :: cmessage             ! error message for downwind routine
    ! program starts here

    ! ----- initialize --------------------------------------------------------------------------------------------

    ! define filenames
    fName_spatial = 'SCRB_1subbasin_depressions_info.csv'
    fName_forcing = 'syntheticForcing.csv'

    ! read the spatial csv file
    call read_csv(fName_spatial,                                      & ! input file name
                  [ixCatArea, ixDepArea, ixDepVol, ixDepCatFrac, ixDep_b, ixDep_p], & ! index of desired variable in data structure
                  [s1       , s2       , s3      , s4          , s5     , s6     ], & ! names of desired variables
                  spatialData, nBasins,                               & ! populated data structures
                  ierr, cMessage)                                       ! error control
    call handle_err(ierr, cMessage)

    ! read the forcing csv file
    call read_csv(fName_forcing,                    & ! input file name
                  [ixP, ixET, ixQ],                 & ! index of desired variable in data structure
                  [f1 , f2  , f3 ],                 & ! names of desired variables
                  forcingData, nTime,               & ! populated data structures
                  ierr, cMessage)                     ! error control
    call handle_err(ierr, cMessage)

    ! check
    print *, 'nBasins = ', nBasins, 'nTimesteps = ', nTime

    ! ----- run HDS -----------------------------------------------------------------------------------------------

    ! run the HDS model
    call run_HDS(nBasins, nTime,                                                                    & ! space/time dimensions
                 spatialData(ixCatArea)%dat, spatialData(ixDepArea)%dat, spatialData(ixDepVol)%dat, & ! spatial information
                 spatialData(ixDepCatFrac)%dat, spatialData(ixDep_b)%dat, spatialData(ixDep_p)%dat, & ! spatial information (cont.)
                 forcingData(ixP)%dat, forcingData(ixET)%dat, forcingData(ixQ)%dat,                 & ! forcing information
                 ierr, cMessage)                                                                      ! error control
    call handle_err(ierr, cMessage)

    ! ----- finalize ----------------------------------------------------------------------------------------------

    ! deallocate space
    do iVar=1,3
     deallocate(spatialData(iVar)%dat, forcingData(iVar)%dat, stat=ierr)
     call handle_err(ierr,'problem deallocating space for dat vector')
    end do
    deallocate(spatialData, forcingData, stat=ierr)
    call handle_err(ierr,'problem deallocating space for data structures')

end program

! =========================================================================
! * subroutine to run HDS
! =========================================================================

subroutine run_HDS(nBasins, nTimesteps,                                & ! space/time dimensions
                   catchmentArea, depressionArea, depressionVol,       & ! spatial information
                   depCatchAreaFrac,b, p,                              & ! spatial information (cont.)
                   pRate, etPond, qSeas,                               & ! forcing information
                   ierr, message)                                        ! error control

    USE type_HDS
    USE HDS

    implicit none
    ! input/output of the subroutine
    integer(i4b), intent(in)  :: nBasins, nTimesteps       ! number of sub-basins and timesteps included in the analysis
    real(rkind),  intent(in)  :: catchmentArea(nBasins)    ! catchment area of the depression in m^2
    real(rkind),  intent(in)  :: depressionArea(nBasins)   ! depression area in m^2
    real(rkind),  intent(in)  :: depressionVol(nBasins)    ! depression volume in m^3
    real(rkind),  intent(in)  :: depCatchAreaFrac(nBasins) ! fractional area (of the landArea= basinArea - depressionArea) that drains to the depressions
    real(rkind),  intent(in)  :: b(nBasins)                ! shape of the fractional contributing area curve [-]
    real(rkind),  intent(in)  :: p(nBasins)                ! shape of the slope profile [-]
    real(rkind),  intent(in)  :: pRate(nTimesteps)         ! forcing data = precipitation [mm/day]
    real(rkind),  intent(in)  :: etPond(nTimesteps)        ! forcing data = ET            [mm/day]
    real(rkind),  intent(in)  :: qSeas(nTimesteps)         ! forcing data = runoff        [mm/day]
    integer(i4b),intent(out)  :: ierr                      ! error code
    character(*),intent(out)  :: message                   ! error message
    ! local variables -- parameters
    real(rkind), parameter    :: tau = 0.01_rkind          ! time constant linear reservoir [days-1]
    real(rkind), parameter    :: dt  = 1.0_rkind           ! time step [days]
    ! local variables
    integer(i4b)              :: ibasin, itime             ! loop counter for basin and timeseries
    real(rkind)               :: landArea(nBasins)         ! land area of the basin in m^2 = catchmentArea - depressionArea
    real(rkind)               :: upslopeArea_dep(nBasins)  ! upslope (upland) area of the depression in m^2 = landArea * depCatchAreaFrac
    real(rkind)               :: upslopeArea_rvr(nBasins)  ! upslope (upland) area of the river in m^2 = landArea * (1-depCatchAreaFrac) part of land area contributing to the river directly
    real(rkind)               :: vMin(nBasins)             ! minimum pond volume below which contributing area is zero [m3]
    real(rkind)               :: depConAreaFrac(nBasins)   ! contributing area fraction of pothole dominated area [-]
    real(rkind)               :: basinConAreaFrac(nBasins) ! contributing area fraction per the entire subbasin from pothole and non-pothole areas [-]
    real(rkind)               :: volFrac(nBasins)          ! volume fraction per subbasin [-]
    real(rkind)               :: pondVol(nBasins)          ! pond volume [m3]
    real(rkind)               :: pondArea(nBasins)         ! pond area [m2] 
    real(rkind)               :: pondOutflow               ! pond outflow [m3]
    real(rkind)               :: basinOutflow              ! outflow from the entire basin (pond + non-pothole areas) [m3]
    real(rkind)               :: totEvap                   ! total evaporation for initialization of the pond [m]
    real(rkind)               :: Q_det_adj, Q_dix_adj      ! adjusted evapotranspiration & infiltration fluxes [L3 T-1] for mass balance closure (i.e., when losses > pondVol). Zero values mean no adjustment needed.
    ! initialize error control
    ierr=0; message='run_HDS/'

    !===============================
    ! initialization
    !===============================

    ! initial assignments
    landArea    = max(catchmentArea - depressionArea, zero)
    upslopeArea_dep  = landArea * depCatchAreaFrac
    upslopeArea_rvr  = landArea * (1._rkind - depCatchAreaFrac)
    vMin(:)     = zero       ! time varying model parameter, will be updated later
    totEvap     = zero       ! initialize pond volume using ET (set to zero to initialize using volume fraction)

    ! initialize variables (all variables will be updated by the model)
    volFrac  = 0.2_rkind     ! assume depressions are 20% full at time = 0 (initial condition)
    depConAreaFrac  = 0.2_rkind     ! assume contributing area is 20% at time = 0 (initial condition)
    !areaFrac = 0.2_rkind     ! assume areafrac 20% at time = 0 (initial condition)
    pondArea = zero          ! updated inside the initialization subroutine
    pondVol  = zero          ! updated inside the initialization subroutine

    ! loop though subbasins to initialize the variables
    do ibasin = 1, nBasins
        call init_pond_Area_Volume(depressionArea(ibasin), depressionVol(ibasin), totEvap, volFrac(ibasin), p(ibasin), pondVol(ibasin), pondArea(ibasin))
        vMin(ibasin) = pondVol(ibasin)
    enddo ! loop for subbasin initialization

    ! create output file and write header
    open(unit=10, file='HDS_output.csv', status='unknown', action='write')
    write(10,*) 'time, basinID, pondVol, volFrac, depConAreaFrac, basinConAreaFrac, vMin, pondArea, pondOutflow, basinOutflow,'

    !===============================
    ! start of time loop (timeseries simulation)
    !===============================
    do itime = 1, nTimesteps

        ! for each timestep, loop through subbasins
        do ibasin = 1, nBasins

            ! run the meta depression model for a single depression
            call runDepression(pondVol(ibasin), qSeas(itime), pRate(itime), etPond(itime), depressionArea(ibasin), depressionVol(ibasin), upslopeArea_dep(ibasin), &
                                p(ibasin), tau, b(ibasin), vMin(ibasin), dt, Q_det_adj, Q_dix_adj, volFrac(ibasin), depConAreaFrac(ibasin), pondArea(ibasin), pondOutflow)
            
            ! get integrated basin contributing area fraction (fraction from pothole area + non-pothole area)
            basinConAreaFrac(ibasin) = ((depConAreaFrac(ibasin) * (depressionArea(ibasin)+upslopeArea_dep(ibasin))) + upslopeArea_rvr(ibasin))/catchmentArea(ibasin)
            ! get basin outflow
            basinOutflow = (qSeas(itime) * upslopeArea_rvr(ibasin)) + pondOutflow
            ! uncomment and impelement the following in the LSM
            ! if(Q_det_adj>zero .or. Q_dix_adj>zero)then
                ! assign the new evaporation and infiltration volumes as the Q_det_adj, Q_dix_adj in the LSM
                ! to ensure mass balance closure
            ! end if
            
            ! save information (bookkeeping for next time step)
            write(*,*) 'time step = ', itime
            write(10,1110) itime, ibasin, pondVol(ibasin), volFrac(ibasin), depConAreaFrac(ibasin), basinConAreaFrac(ibasin), vMin(ibasin), pondArea(ibasin), pondOutflow, basinOutflow
            1110    format(9999(g15.7e2, ','))

        enddo ! loop for subbasin

    enddo ! loop for time

    close(10) ! close HDS_output file

end subroutine run_HDS


! **************************************************************************************************
! error handler
! **************************************************************************************************
subroutine handle_err(ierr,message)
 USE type_HDS
 implicit none
 integer(i4b),intent(in)            :: ierr            ! error code
 character(*),intent(in)            :: message         ! error message
 ! ---------------------------------------------------------------------------------------
 ! return if A-OK
 if(ierr==0) return

    ! process error messages
    if (ierr>0) then
     write(*,'(//a/)') 'FATAL ERROR: '//trim(message)
    else
     write(*,'(//a/)') 'WARNING: '//trim(message); print*,'(can keep going, but stopping anyway)'
    endif
    stop 1

end subroutine handle_err
