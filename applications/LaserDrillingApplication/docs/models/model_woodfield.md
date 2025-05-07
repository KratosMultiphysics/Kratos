
### Woodfield model [^1]
The energy of a pulse is deposited instantaneously within the material according to the Beer-Lambert law; time stepping is frozen.
Where the energy exceeds a threshold, the material is ablated away immediately, which in the simulation entails that the corresponding 
elements are deactivated instantaneously. Then, time integration of the heat equation resumes until the following pulse. 
The ablation (the removal of elements) leaves a crater. Since this crater is very shallow compared to the Rayleigh range of the laser, 
we neglect its depth when applying the next pulse.

[^1]:[Optical penetration models for practical prediction of femtosecond laser ablation of dental hard tissue](https://doi.org/10.1002/lsm.23784)
