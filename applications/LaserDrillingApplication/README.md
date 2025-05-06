 ## Laser Drilling Application  
The LaserDrillingApplication is a Kratos application built on top of the ConvectionDiffusionApplication to simulate the process of laser drilling with a femtosecond laser.

### Introduction
In essence, it works by distributing the energy of the laser pulse inside of the sample and letting the ConvectionDiffusionApplication perform the heat diffusion over time. Where the energy or the temperature exceed certain thresholds, the material is modified. For instance, if the temperature is high enough in a region, that part of the material suffers chemical changes and is carbonized. Or, where the energy is above a certain value, ablation happens and that part of the material is removed.


### Physical model
Currently (May 2025), it is based on the model proposed by Woodfield et al. in [Optical penetration models for practical prediction of femtosecond laser ablation of dental hard tissue](https://doi.org/10.1002/lsm.23784) using an axisymmetric geometry.

The energy of a pulse is deposited within the material according to the Beer-Lambert law. Where the energy exceeds a threshold, the material is ablated away immediately, which in the simulation entails that the corresponding elements are deactivated instantaneously. Then, the heat is propagated until the following pulse. The ablation (the removal of elements) leaves a crater. Since this crater is very shallow compared to the Rayleigh range of the laser, we neglect its depth when applying the next pulse.
