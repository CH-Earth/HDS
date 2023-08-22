program main_HDS
    ! A modified version (v2) of the Hysteretic Depressional Storage (HDS) model to simulate prairie fill and spill mechanism.
    ! The original HDS (v1) was based on the equations listed in Ahmed et al. (2023), doi: https://doi.org/10.1016/j.envsoft.2023.105769.
    ! This version of HDS (v2) is revised based on the equations presented by Clark and Shook (2022), doi: https://doi.org/10.1029/2022WR032694,
    ! which are more numerically effecient than the ones proposed in HDS (v1)

    use HDS

    implicit none
    character(len=100) :: fName ! file name that contains depressional storage information
    !character(len=100) :: trash ! reads extra (unneeded) variables
    real(rkind)  ::  dummy               ! dummy variable used to read extra data not included in the analysis
    integer :: GetNumberOfLines !function to get the number of lines of the input file
    integer :: nlines ! number of lines of the input file
    integer :: ibasin, itime ! loop counter for basin and timeseries
    integer :: nbasin, ntimesteps ! number of sub-basins and timesteps included in the analysis
    real(rkind),  allocatable :: depressionArea(:)  ! depression area in m^2
    real(rkind),  allocatable :: depressionVol(:)   ! depression volume in m^3
    real(rkind),  allocatable :: catchmentArea(:)   ! Catchment area of the depression in m^2
    real(rkind),  allocatable :: upslopeArea(:)   ! upslope (upland) area of the depression in m^2 = catchmentArea - depressionArea
    real(rkind),  allocatable :: qSeas(:), pRate(:), etPond(:)     ! forcing data = runoff, precipitation, ET [mm/day]
    real(rkind)  ::  dt ! time step [days]
    real(rkind)  ::  p ! shape of the slope profile [-]. Exponent for calculating the fractional wet area
    real(rkind)  ::  b ! shape of the fractional contributing area curve [-]
    real(rkind)  ::  tau ! time constant linear reservoir [days-1] for infiltration calculations
    !real(rkind)  ::  seed ! seed number used for generating stochastic simulations
    real(rkind)  ::  totEvap ! total evaporation for initialization of the pond [m]
    real(rkind),  allocatable :: conArea(:) ! contributing area fraction per subbasin [-]
    real(rkind),  allocatable :: volFrac(:) !volume fraction per subbasin [-]
    real(rkind),  allocatable :: areaFrac(:)  !area fraction per subbasin [-]
    real(rkind),  allocatable :: pondVol(:)  !pond volume [m3] [-]
    real(rkind),  allocatable :: vMin(:)  ! minimum pond volume below which contributing area is zero [m3]
    real(rkind)  ::  pondVol0 !, pondVol ! initial and final pond volume [m3]
    !real(rkind)  ::  vMin     ! minimum pond volume below which contributing area is zero [m3]

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !read inputs
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! read depressions spatial properties
    fName = 'SCRB_1subbasin_depressions_info.csv'
    nlines = GetNumberOfLines(fName) !get number of lines (number of basins in file)
    nbasin = nlines - 1
!    write(*,*) nlines
    allocate(depressionArea(nbasin), depressionVol(nbasin), catchmentArea(nbasin), upslopeArea(nbasin))
    !read file
    open(unit=1, file=trim(fName), status='old', action='read')
    read(1,*) !trash !skip header line
    do ibasin = 1, nbasin
!        write(*,*)i
        ! subbasinID	depression_area_m2	depression_volume_m3	total_catchment_m2
        read(1,*) dummy, depressionArea(ibasin), depressionVol(ibasin), catchmentArea(ibasin)
        !calculate the upslope (upland) area of each depression
        upslopeArea(ibasin) = max(catchmentArea(ibasin) - depressionArea (ibasin), zero )
    end do !loop for subbasin
    close(1)
!    write(*,*) depressionArea, depressionVol, catchmentArea, upslopeArea
    allocate(conArea(nbasin), volFrac(nbasin), areaFrac(nbasin), pondVol(nbasin), vMin(nbasin))

    !read synthetic forcings
    fName = 'syntheticForcing.csv'
    nlines = GetNumberOfLines(fName) !get number of lines (number of depressions in file)
    ntimesteps = nlines - 1
    allocate(qSeas(ntimesteps), pRate(ntimesteps), etPond(ntimesteps))
    !read file
    open(unit=1, file=trim(fName), status='old', action='read')
    read(1,*) ! trash !skip header line
    do itime = 1, ntimesteps
!        write(*,*)i
        ! t	qSeas	pRate	etPond
        read(1,*) dummy, qSeas(itime), pRate(itime), etPond(itime) !all in mm/day
    end do !loop for time
    dt = one  !time steps in days (based on the synthetic data)
    close(1)

    ! define parameters (for SCRB)
    p   = 1.72  ! shape of the slope profile [-]
    b   = 1.5   ! shape of the fractional contributing area curve [-]
    tau = 0.01_rkind  ! time constant linear reservoir [days-1]
    vMin(:) = zero  ! model parameter, will be updated later
    ! initialize pond volume using ET (set to zero to initialize using volume fraction)
    totEvap = zero  ! m
    !initialize variables (all variables will be updated by the model)
    volFrac = 0.2 ! assume depressions are 20% full at time = 0 (initial condition)
    conArea = 0.2 ! assume contributing area is 20% at time = 0 (initial condition)
    areaFrac = 0.2 ! assume areafrac 20% at time = 0 (initial condition)
    pondVol = zero  !updated inside the initialization subroutine

    !conArea = dblarr(nDepressions) at time = 0 (initial condition)
    !loop though subbasins to initialize the variables
    do ibasin = 1, nbasin
        call init_pondVolume(depressionArea(ibasin), depressionVol(ibasin), totEvap, volFrac(ibasin), p, pondVol(ibasin))
        pondVol0 = pondVol(ibasin)
        vMin(ibasin) = pondVol(ibasin)
!        write(*,*) pondVol0, vMin(ibasin)
    enddo ! loop for subbasin initialization
    ! initialize arrays
    ! volFrac     = dblarr(nDays)
    ! areaFrac    = dblarr(nDays)
    ! volFracPond = dblarr(nDepressions,nDays)
    !===============================
    ! create output file
    !===============================
    open(unit=10, file='HDS_output.csv', status='unknown', action='write')

    write(10,*) 'time, basinID, pondVol, volFrac, conArea, vMin, '

    !===============================
    ! start of time loop (timeseries simulation)
    !===============================
    do itime = 1, ntimesteps
        ! for each timestep, loop through subbasins
        write(*,*) 'time step = ', itime
        do ibasin = 1, nbasin
            ! run the meta depression model for a single depression

            call runDepression(pondVol(ibasin), qSeas(itime), pRate(itime), etPond(itime), depressionArea(ibasin), depressionVol(ibasin), upslopeArea(ibasin), &
                                p, tau, b, vMin(ibasin), dt, volFrac(ibasin), conArea(ibasin))

            ! save information (bookkeeping for next time step)
            !volFrac(ibasin)       = fVol
            !conArea(ibasin)      = fArea
            !volFracPond(ibasin) = pondVol/depVol
            write(10,1110) itime, ibasin, pondVol(ibasin), volFrac(ibasin), conArea(ibasin), vMin(ibasin)
            1110    format(9999(g15.7e2, ','))

        enddo !Loop for subbasin

    enddo !loop for time

    close(10) !close HDS_output file


end program
! ======================================================
! local functions used for pre/post processing
! ======================================================
! get the number of files
function GetNumberOfLines(filename) result(num_lines)
    implicit none
    character(len=100), intent(in) :: filename
    integer :: num_lines, i
    character(1000) :: line

    ! Open the file
    open(unit=1, file=trim(filename), status='old', action='read', iostat=i)

    ! Check if the file was opened successfully
    if (i /= 0) then
        write(*, *) "Error opening the file. Check if the file exists"
        stop
    end if

    ! Initialize the line counter
    num_lines = 0

    ! Loop through the file and count the lines
    do
        read(1, '(A)', iostat=i) line
        if (i == 0) then
            num_lines = num_lines + 1
        else if (i /= 0 .and. i /= -1) then
            write(*, *) "Error reading the file."
            exit
        else
            exit
        end if
    end do

    ! Close the file
    close(1)

end function GetNumberOfLines
