# HDS
A modified version (v2) of the Hysteretic Depressional Storage (HDS) model to simulate prairie fill and spill mechanism.

The [original HDS (v1)](https://github.com/UC-HAL/HYPE-HDS) was based on the equations listed in [Ahmed et al. (2023)](https://doi.org/10.1016/j.envsoft.2023.105769). This version of HDS (v2) is revised based on the equations presented by [Clark and Shook (2022)](https://doi.org/10.1029/2022WR032694), which is cleaner and more robust than the one proposed in HDS (v1).

# Notes on the implementation

HDS should be implemented as a module just before the routing module in the model (i.e., after all runoff fluxes have been calculated). By doing so, the modeller ensures that the runoff that is being routed (reaches the river network) is the net runoff (after subtracting the stored water in the depressions).

HDS mainly needs the following information:

* input fluxes/forcing: runoff, precipitation, ET (for pond/depression evaporation) in mm/timestep.
* Subbasin/grid properties: depressions area (m<sup>2</sup>), depressions volume (m<sup>3</sup>), and total catchment area (m<sup>2</sup>)
* Parameters (currently hard coded): `p` (shape of the slope profile [-]), `b` (shape of the fractional contributing area curve [-]), `tau` (time constant linear reservoir to simulate seepage losses from depressions [days<sup>-1</sup>]), `rCoef` (runoff coefficient [-]).

All these information are required at the sub-basin or grid scale. Please note that the two parameters (`tau` and `rCoef`) are conceptual, and it is better if not used when coupling HDS with a hydrological model. These parameters were used for comparison purposes against a synthetic test case from [Clark and Shook (2022)](https://doi.org/10.1029/2022WR032694). So, it is better to set their values to zero so that they do not affect the results. If the model has a way to handle seepage losses from water bodies, then that can be passed into the model as input flux, but it would require some modifications in the function input argument.

For the rest of the parameters (`p` and `b`), please leave their value unchanged as these are the typical values for the Smith Creek Basin.

Please be aware that the following variables are state variables within the code, and it is crucial to ensure their proper maintenance from one time step to another for proper bookkeeping:
| variable  | description  |
|---|---|
| `vMin`     |          ! minimum pond volume below which contributing area is zero [m<sup>3</sup>] |
| `conArea`  |          ! contributing area fraction per subbasin [-] |
| `pondVol`  |          ! pond volume [m<sup>3</sup>] |
| `pondArea` |          ! pond area [m<sup>2</sup>] |

# Change log:
Changes are sorted from newest to oldest:
## Oct 27, 2023
* Calculate adjusted ET and infiltation for mass balance closure when pond is dry (losses > pondVol)
## Oct 19, 2023
* add constrains to avoid `pondVol` being greater than `depVol` under extremely wet conditions
## Oct 18, 2023
* add if statements to prevent negative pondVol values under extremely dry conditions
## Sep 1, 2023
* produce the outflow as an output timeseries.

## Aug 30, 2023
1. simplify code to combine `runDepression` and `runOnestep` into one subroutine `runDepression` and remove `runOnestep`.
2. Produce the pond area as an output timeseries.
3. Remove extra internal variables that are not needed.

## (Aug 30, 2023
1. Modify the makefile to activate `release` and `debug` targets. The previous version only used `release` by default as the target and `debug` was not activated.
2. Fix indentation (space -> tab) in `ascii_util.f90` and `type_HDS.f90`.
