 ## Laser Drilling Application  
The LaserDrillingApplication is a Kratos application built on top of the ConvectionDiffusionApplication to simulate the process of laser drilling with a femtosecond laser.

### Introduction
In essence, it works by distributing the energy of the laser pulse inside of the sample and letting the ConvectionDiffusionApplication perform the heat diffusion over time. Where the energy or the temperature exceed certain thresholds, the material is modified. For instance, if the temperature is high enough in a region, that part of the material suffers chemical changes and is carbonized. Or, where the energy is above a certain value, ablation happens and that part of the material is removed.


### Physical models
The current model (May 2025), is based on the model proposed by Woodfield et al. using an axisymmetric geometry.

- [Woodfield model](./docs/models/model_woodfield.md)
