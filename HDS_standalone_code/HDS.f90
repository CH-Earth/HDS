module HDS
USE type_HDS

    contains

    subroutine init_pondVolume(depArea, depVol, totEvap, volFrac, p, pondVol)
        ! used to initialize pond volume
        ! initialize at capacity minus depth of evaporation
        implicit none
        !subroutine arguments
        real(rkind),  intent(in) :: depArea !depression area [m2]
        real(rkind),  intent(in) :: depVol !depression volume [m3]
        real(rkind),  intent(in) :: totEvap ! total evaporation [m] to be subtracted from the ponded water
        real(rkind),  intent(in) :: volFrac ! volume fraction [-]. Used to fill the depressions using percentage (i.e., 50% full)
        real(rkind),  intent(in) :: p ! shape of the slope profile [-]. Exponent for calculating the fractional wet area
        real(rkind),  intent(out) :: pondVol ! pond volume [m3]
        ! local variables
        real(rkind)  ::  depHeight ! depression height [m]
        real(rkind)  ::  pondHeight ! height of the water level in depression [m]
        ! estimate the maximum depth of the depression
        depHeight = depVol*(1.0_rkind  + 2.0_rkind /p)/depArea  ! depression height [m]
        !ixWater   = where(depHeight gt totEvap, nWater)
        !subtract ET from ponded water if ET > 0
        if(totEvap .gt. 0.0_rkind ) then
            ! estimate the height of the water level after evaporation
            pondHeight = max(depHeight - totEvap, 0.0_rkind ) !0.0_rkind  !replicate(0.d, nDepressions)
            !if(depHeight gt 0)then pondHeight[ixWater] = depHeight[ixWater] - totEvap

            ! get the volume of water in the depression
            pondVol = depVol*((pondHeight/depHeight)**(1.0_rkind  + 2.0_rkind /p))

        else !use volume frac to get initial pondVol
            pondVol = depVol*volFrac
        end if

    end subroutine init_pondVolume

    !=============================================================
    !=============================================================
    subroutine runDepression(pondVol,                  & ! input/output:  state variable = pond volume [m3]
                   qSeas, pRate, etPond,     & ! input:         forcing data = runoff, precipitation, ET [mm/day]
                   depArea, depVol, upsArea, & ! input:         spatial attributes = depression area [m2], depression volume [m3], upstream area [m2]
                   p, tau,                   & ! input:         model parameters = p [-] shape of the slope profile; tau [day-1] time constant linear reservoir
                   b, vmin,                  & ! input:         model parameters = b [-] shape of contributing fraction curve; vmin [m3] minimum volume
                   dt,                       & ! input:         model time step [days]
                   fVol, fArea)                ! output:        fractional volume and fractional contributing area [-]


        ! used to estimate the volumetric storage at the end of the time interval
        implicit none
        !subroutine arguments
        real(rkind),  intent(inout) :: pondVol, vmin               !state variable = pond volume [m3]; vmin [m3] minimum volume
        real(rkind),  intent(in)    :: qSeas, pRate, etPond        ! input:         forcing data = runoff, precipitation, ET [mm/day]
        real(rkind),  intent(in)    :: depArea, depVol, upsArea    ! input:         spatial attributes = depression area [m2], depression volume [m3], upstream area [m2]
        real(rkind),  intent(in)    :: p, tau                      ! input:         model parameters = p [-] shape of the slope profile; tau [day-1] time constant linear reservoir
        real(rkind),  intent(in)    :: b                           ! input:         model parameters = b [-] shape of contributing fraction curve;
        real(rkind),  intent(in)    :: dt                          ! input:         model time step [days]
        real(rkind),  intent(out)   :: fVol, fArea                ! output:        fractional volume and fractional contributing area [-]
        !local variables
        real(rkind)  ::  vMinNew, pondVolNew                        !new state variable = pond volume [m3]; vmin [m3] minimum volume
        ! return upon failure
        !on_error, 2

        ! define the desired depression
        !ixDesired = 8157

        ! get the number of depressions
        !nDepressions = n_elements(pondVol)

        ! initilaize the contributing area
        !conArea = 0.0_rkind  !replicate(0.d, nDepressions)

        ! loop through depressions
        !for iDep=0,nDepressions-1 do begin

        ! skip
        !print, iDep
        !if(iDep ne ixDesired)then continue

        ! run model for one time step and one depression

        ! ensemble depression
        !if keyword_set(isEnsemble) then $
        !runOnestep, pondVol[iDep], qSeas, pRate, etPond, depArea[iDep], depVol[iDep], upsArea[iDep], $
        !            p, tau, b, vmin[iDep], dt, pondVolNew, vMinNew, fArea, /ensemble

        ! meta depression
        !if keyword_set(isMeta) then begin
        call runOnestep(pondVol, qSeas, pRate, etPond, depArea, depVol, upsArea, &
                        p, tau, b, vmin, dt, pondVolNew, vMinNew, fArea)

        !endif

        ! save variables
        vMin = vMinNew
        pondVol = pondVolNew
        ! compute fractional volume and fractional area
        fVol  = pondVol/depVol!total(pondVol)/total(depVol)


        !print, 'fVol = ', pondVol[0:1]/depVol[0:1]
        !print, 'fVolVec = ', pondVol/depVol
        !print, 'forcing = ', qSeas, pRate, etPond

    end subroutine runDepression

    !=============================================================
    !=============================================================
    subroutine runOnestep(pondVol,                  & ! input:  state variable = pond volume [m3]
                qSeas, pRate, etPond,     & ! input:  forcing data = runoff, precipitation, ET [mm/day]
                depArea, depVol, upsArea, & ! input:  spatial attributes = depression area [m2], depression volume [m3], upstream area [m2]
                p, tau,                   & ! input:  model parameters = p [-] shape of the slope profile; tau [day-1] time constant linear reservoir
                b, vMinOld,               & ! input:  model parameters = b [-] shape of contributing fraction curve; vmin [m3] minimum volume
                dt,                       & ! input:  model time step [days]
                xVol, vMin, fArea)          ! output: pond volume at the end of the time step [m3], fractional contributing area [-]

        implicit none
        !subroutine arguments
        real(rkind),  intent(inout) :: pondVol               !state variable = pond volume [m3]
        real(rkind),  intent(in)    :: qSeas, pRate, etPond        ! input:         forcing data = runoff, precipitation, ET [mm/day]
        real(rkind),  intent(in)    :: depArea, depVol, upsArea    ! input:         spatial attributes = depression area [m2], depression volume [m3], upstream area [m2]
        real(rkind),  intent(in)    :: p, tau                      ! input:         model parameters = p [-] shape of the slope profile; tau [day-1] time constant linear reservoir
        real(rkind),  intent(in)    :: b, vMinOld                  ! input:         model parameters = b [-] shape of contributing fraction curve; vminold [m3] minimum volume
        real(rkind),  intent(in)    :: dt                          ! input:         model time step [days]
        real(rkind),  intent(out)   :: xVol, vMin, fArea           !output: pond volume at the end of the time step [m3], fractional contributing area [-]
        !local variables
        integer :: implicitEuler, shortSubsteps, solution ! flags to activate the required solution
        real(rkind)  ::  xConv ! convergence criteria (algorithm control parameter)
        integer :: nIter, iter ! number of iterations (algorithm control parameter), iteration counter
        integer :: nSub ! number of substeps (algorithm control parameter)
        real(rkind)  ::  xMin, xMax ! min and max storage
        real(rkind)  ::  Q_di, Q_det, Q_dix, Q_do !fluxes [L3 T-1]: sum of water inputs to the pond, evapotranspiration, infiltration, pond outflow
        real(rkind)  ::  cFrac, g, dgdv !contributing fraction, net fluxes and derivative
        real(rkind)  ::  xRes !residual
        integer :: failure !failure flag
        real(rkind)  ::  dtSub !dt for shortSubsteps solution
        integer :: iSub !counter for shortSubsteps solution
        ! define numerical options
        implicitEuler = 1001
        shortSubsteps = 1002
        solution      = implicitEuler
        !solution      = shortSubsteps

        ! define algorithmic control parameters
        xConv = 1.e-6 ! convergence criteria
        nIter = 100   ! number of iterations
        nSub  = 100   ! number of substeps

        ! initialize pond volume
        xVol = pondVol
        vMin = vMinOld

        ! ---------- option 1: implicit Euler ----------

        ! check if implicit Euler is desired
        if(solution .eq. implicitEuler)then

            ! initialize brackets
            xMin = 0.0_rkind
            xMax = depVol
            !if keyword_set(isMeta)     then xMax = depVol
            !if keyword_set(isEnsemble) then xMax = depVol*10.d ! the ensemble model computes an intermediate solution without spill

            ! iterate (can start with 1 because the index is not used)
            do iter=1,nIter

                ! compute fluxes and derivatives
                call computFlux(iter, xVol, qSeas, pRate, etPond, depArea, depVol, upsArea, p, tau, b, vMin, Q_di, Q_det, Q_dix, Q_do, cFrac, g, dgdv)
                !if keyword_set(isMeta)     then computFlux, iter, xVol, qSeas, pRate, etPond, depArea, depVol, upsArea, p, tau, b, vMin, Q_di, Q_det, Q_dix, Q_do, cFrac, g, dgdv, /meta
                !if keyword_set(isEnsemble) then computFlux, iter, xVol, qSeas, pRate, etPond, depArea, depVol, upsArea, p, tau, b, vMin, Q_di, Q_det, Q_dix, Q_do, cFrac, g, dgdv, /ensemble
                !print, g, dgdv

                ! compute residual
                xRes = (pondVol + g*dt) - xVol

                ! check convergence
                if(iter .gt. 1 .and. abs(xRes) .lt. xConv)then
                    failure=0
                    exit
                endif

                ! update constraints
                if(xRes .gt. 0.0_rkind ) xMin=xVol
                if(xRes .lt. 0.0_rkind ) xMax=xVol

                ! special case where xMax is too small
                if(xRes .gt. 0.0_rkind  .and. xVol .gt. 0.99*xMax) xMax = xMax*10.0_rkind

                ! update state (pondVol)
                xVol = xVol + xRes / (1.0_rkind  - dgdv*dt)

                ! use bi-section if violated constraints
                if(xVol .lt. xMin .or. xVol .gt. xMax) xVol=(xMin+xMax)/2.0_rkind
                !if(xVol /= xVol) xVol = 0.0_rkind  !NaN check MIA (does not work)
                !xVol = max(xVol, 0.0_rkind ) !prevent -ve values of xVol
                !print, iter, pondVol, xVol, depArea, depVol, Q_di, Q_det, Q_dix, Q_do, xRes, xMin, xMax

                ! assign failure
                if(iter .eq. nIter) failure=1

            enddo ! iterating

        endif  ! if implicit Euler

        ! ---------- option 2: short substeps ----------

        ! check if short substeps is desired
        if(solution .eq. shortSubsteps .or. failure .eq. 1)then

            ! define length of the substeps
            dtSub = 1.0_rkind  / real(nSub)

            ! loop through substeps
            do iSub=1,nSub
                call computFlux(iter, xVol, qSeas, pRate, etPond, depArea, depVol, upsArea, p, tau, b, vMin, Q_di, Q_det, Q_dix, Q_do, cFrac, g, dgdv)
                xVol = xVol + g*dtSub
                !print, iSub, qSeas, pRate, etPond, cFrac, xVol
            enddo  ! looping through substeps

        endif  ! if short substeps

        ! check
        !if(failure eq 1)then stop, 'check failure'

        ! ---------------------------------------------------------------------------------
        ! ---------------------------------------------------------------------------------
        ! ---------- fill and spill process for the meta depression model -----------------
        fArea = cFrac
        ! ---------------------------------------------------------------------------------


    end subroutine runOnestep

    !=============================================================
    !=============================================================
    subroutine computFlux(iter, pondVol,             & ! input:  iteration index, state variable = pond volume [m3]
                qSeas, pRate, etPond,      & ! input:  forcing data = runoff, precipitation, ET [mm/day]
                depArea, depVol, upsArea,  & ! input:  spatial attributes = depression area [m2], depression volume [m3], upstream area [m2]
                p, tau,                    & ! input:  model parameters = p [-] shape of the slope profile; tau [day-1] time constant linear reservoir
                b, vmin,                   & ! input:  model parameters = b [-] shape of contributing fraction curve; vmin [m3] minimum volume
                Q_di, Q_det, Q_dix, Q_do,  & ! output: individual model fluxes
                cFrac, g, dgdv)              ! output: contributing fraction, net fluxes and derivative

        implicit none
        ! subroutine arguments
        integer, intent(in) :: iter                     !iteration index
        real(rkind),  intent(in) :: pondVol                     !state variable = pond volume [m3]
        real(rkind),  intent(in) :: qSeas, pRate, etPond        ! input:  forcing data = runoff, precipitation, ET [mm/day]
        real(rkind),  intent(in) :: depArea, depVol, upsArea    ! input:  spatial attributes = depression area [m2], depression volume [m3], upstream area [m2]
        real(rkind),  intent(in) :: p, tau                      ! input:  model parameters = p [-] shape of the slope profile; tau [day-1] time constant linear reservoir
        real(rkind),  intent(in) :: b                           ! input:  model parameters = b [-] shape of contributing fraction curve; vmin [m3] minimum volume
        real(rkind),  intent(inout) :: vmin                     ! in/out: vmin [m3] minimum volume
        real(rkind),  intent(out) :: Q_di, Q_det, Q_dix, Q_do   ! output: individual model fluxes
        real(rkind),  intent(out) :: cFrac, g, dgdv             ! output: contributing fraction, net fluxes and derivative, pond area
        ! local variables
        real(rkind)  ::  ms ! smoothing parameter (algorithm control)
        real(rkind)  ::  pondArea ! calculated pond area
        real(rkind)  ::  rCoef    ! runoff coefficient [-]
        real(rkind)  ::  pInput  ! precipitation (or rain+melt) [m/day]
        real(rkind)  ::  qInput  ! surface runoff [m/day]
        real(rkind)  ::  eLosses  ! evaporation losses [m/day]
        real(rkind)  ::  vPrime !smoothed pondVol value for Euler solution
        real(rkind)  ::  vTry ! adjusted pondVol for derivatives calculations
        real(rkind)  ::  dadv, dpdv, dfdv, didv

        ! define algorithmic control parameters
        ms = 0.0001_rkind

        ! compute the pond area
        pondArea = depArea*((pondVol/depVol)**(2.0_rkind /(p + 2.0_rkind )))

        ! get the forcing
        rCoef   = 0.050_rkind ! (runoff coefficient)
        pInput  = pRate/1000.0_rkind  !mm -> m
        qInput  = qSeas/1000.0_rkind  + rCoef*pRate/1000.0_rkind  ! surface runoff !mm -> m
        eLosses = etPond/1000.0_rkind  ! evaporation losses mm -> m

        ! get volume fluxes from the host land model
        ! sum of water input to the depression (eq 11)
        Q_di  = upsArea*qInput + (depArea - pondArea)*qInput + pondArea*pInput
        ! evapotranspiration losses (eq 12)
        Q_det = pondArea*eLosses

        ! compute infiltration from the bottom of the pond (eq 13)
        Q_dix = tau*pondVol

        ! update the minimum parameter (force hysteresis)
        if(iter .eq. 1) then
            if(Q_di .lt. Q_det+Q_dix) vmin = pondVol
        endif

        ! compute the outflow from the meta depression
        ! smoothing pondVol calculation ! (eq 24)
        vPrime = 0.50*(vMin + pondVol + sqrt((pondVol-vMin)**2.0_rkind  + ms)) ! (eq 24)
        ! calculate contributing fraction !(eq 25)
        cFrac  = 1.0_rkind  - ((depVol - vPrime)/(depVol - vMin))**b !(eq 25)
        cFrac = max(cFrac, 0.0_rkind ) !prevent -ve values MIA
        ! calculate outlfow (eq 16)
        Q_do   = cFrac*Q_di

        ! compute the net flux (eq 10)
        g = Q_di - Q_det - Q_dix - Q_do

        ! compute derivate in pond area w.r.t. pond volume
        vTry = max(1.0e-12_rkind, pondVol)  ! avoid divide by zero
        dadv = ((2.0_rkind *depArea)/((p + 2.0_rkind )*depVol)) * ((vTry/depVol)**(-p/(p + 2.0_rkind )))

        ! compute derivatives in volume fluxes w.r.t. pond volume (meta depression)
        dpdv = 0.50*((pondVol-vMin)/sqrt((pondVol-vMin)**2.0_rkind  + ms) + 1.0_rkind )  ! max smoother
        dfdv = dpdv*(b*(depVol - pondVol)**(b - 1.0_rkind ))/((depVol - vMin)**b)  ! contributing fraction (note chain rule)
        didv = dadv*(pInput - qInput)
        dgdv = didv*(1.0_rkind  - cFrac) - dfdv*Q_di - dadv*eLosses - tau

        !write(*,*) g, pondVol, dadv, dpdv, dfdv, didv, dgdv

    end subroutine computFlux

end module HDS
