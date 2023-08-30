# HDS
A modified version (v2) of the Hysteretic Depressional Storage (HDS) model to simulate prairie fill and spill mechanism.

The [original HDS (v1)](https://github.com/UC-HAL/HYPE-HDS) was based on the equations listed in Ahmed et al. (2023), doi: https://doi.org/10.1016/j.envsoft.2023.105769. This version of HDS (v2) is revised based on the equations presented by Clark and Shook (2022), doi: https://doi.org/10.1029/2022WR032694, which is cleaner and more robust than the one proposed in HDS (v1).

# Changes in the repo
1. Modify the makefile to activate `release` and `debug` targets. The previous version only uses `release` by default as the target and `debug` was not activated.
2. Fix indentation (space -> tab) in `ascii_util.f90` and `type_HDS.f90`.
3. Add the old version of HDS (v1) and allow it to work with the same synthetic data used for HDS_v2.
