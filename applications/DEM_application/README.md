# DEM Application

This application focuses on the Discrete Element Method (DEM), a particles method for modeling the bulk behavior of granular materials and many geomaterials such as coal, ores, soil, rocks, aggregates, pellets, tablets and powders.

The [DEMpack Team](www.cimne.com/dem) at [CIMNE](www.cimne.com) is in charge of all developments related to the DEM.

For the coupling between DEM and Fluids, go to the [Swimming DEM Application](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/swimming_DEM_application).

## Getting Started

This application is part of the Kratos Multiphysics Platform. Instructions on how to get you a copy of the project up and running on your local machine for development and testing purposes are available for both [Linux](http://kratos-wiki.cimne.upc.edu/index.php/LinuxInstall) and [Windows](http://kratos-wiki.cimne.upc.edu/index.php/Windows_7_Download_and_Installation) systems.

### Prerequisites

Build [Kratos](https://github.com/KratosMultiphysics/Kratos/wiki) and make sure that you put

``` cmake
-DDEM_APPLICATION=ON
```

between the compilation options, so the DEM application is compiled.

No auxiliar external libraries are needed.

## Theory

The DEM is a numerical method that has been applied to simulate and analyze flow behavior in a wide range of disciplines including mechanical and process engineering, pharmaceutical, materials science, agricultural engineering and more.
Coupling with fluid is already available through the Swimming-DEM application also integrated in the Kratos Multiphysics Platform.

The fundamental theoretical background corresponding to the discontinuous (granular matter) part of the code can be found in the DEM literature easily.

### Contact Laws

The contact laws are implemented in [this folder](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/DEM_application/custom_constitutive). Note that a letter 'D' or 'd' in the file name means 'discontiuum'. It is, related to non cohesive or slightly cohesive contacts.

#### Linear Repulsive Force

The most simple representation of a repulsive contact force between a sphere and a wall is given by a linear law, where the force acting on the sphere when contacting the wall is a linear function of the indentation, which in turn would bring a quadratic dependence with the contact radius.

#### Non-Linear Repulsive Force

Hertz solved in 1882 the non-cohesive normal contact between a sphere and a plane. In 1971 Johnson, Kendall and Roberts presented the solution (JKR-Theory) for the same problem in this case adding cohesive behaviour. Not much later, Derjaguin, Müller and Toporov published similar results (DMT-Theory).
Both theories are very close and correct, while the JKR theory is adequate to the study of flexible, large spheres, the DMT theory is specially suited to represent the behaviour of rigid, small ones.

## Numerical approach

The application includes three types of DEM elements used for different purposes:

* Spheric Particle - Base element used to simulate granular materials (non cohesive or slightly cohesive)
* Spheric Continuum Particle - With specific build-in variables to simulate fracture in cohesive materials. It can also be understood as a discretization method of the continuum, with spheres.

Together with:

* FEM-based walls - Objects that can not be crossed by DEM spheres.
* Inlets - Generate new particles while running the simulation, with a certain velocity and linked to some material properties.

### DEM strategies

#### Non-cohesive materials Strategy

Once contact between two spheres occurs, the forces at the contact point are computed. The interaction between the two contacting spheres can be represented by two forces with the same module but opposite directions. This force F can be decomposed into its normal and shear components Fn and Fs, respectively.
The contact interface for the simplest formulation is characterized by the normal and tangential stiffness Kn and Ks, respectively, a frictional device obeying the Couloumb law with a frictional coefficient, and a dashpot defined by a contact damping coefficient.

In order to represent irregular particles with spheres, a numerical correction is used. the rolling friction imposes a virtual moment opposite to particle rotation and dependent on its size.

#### Continuum materials Strategy

For continuum materials simulations, the contact between particles can resist tractions up to a certain limit, when breakage occurs. Depending on the chosen constitutive law, the computation of the forces changes. In the basic versions, a bond strategy is used, but more advanced laws use a non-local stress-tensor based strategy.

### DEM integration schemes

The standard translational and rotational equations for the motion of rigid bodies are used to compute the dynamics of the spheres and clusters. The following schemes can be chosen separately for translation and rotation:

* Symplectic Euler Scheme
* Forward Euler Scheme
* Taylor scheme
* Verlet Velocity Scheme

### Contact search

The contact detection basically consists in determining, for every target object, what other objects are in contact with, and then, judge for the correspondent interaction. It is usually not needed to perform a search every time step, which is generally limited for the stability of the explicit integration of the equations of motion.
A bins based technique is used for this purpose.

## Available Interfaces

### G-DEMPack

G-DEMPack (formerly D-DEMpack) is the package that allows a user to create, run and analyze results of a DEM simulation for discontinuum / granular / little-cohesive materials. Requires [GiD](https://www.gidhome.com/) - Pre and Post Processing software.

You can read the G-DEMPack manual or follow the G-DEMPack Tutorials included in DEMPack Tutorials for fast learning on how to use the GUI.

### C-DEMPack

C-DEMPack combines the features of G-DEMPack with the simulation of continuum/cohesive materials. Refer to the C-DEMPack manual or follow the C-DEMPack Tutorials included in DEMPack Tutorials for more information.

### F-DEMPack

F-DEMPack allows you to simulate a wide spectrum of problems involving the interaction of granular DEM and fluid.

* **Miguel Angel Celigueta** - *Core Development* - [maceli@cimne.upc.edu](mailto:maceli@cimne.upc.edu)
* **Salva Latorre** - *Granular materials* - [latorre@cimne.upc.edu](mailto:latorre@cimne.upc.edu)
* **Ferran Arrufat** - *Cohesive materials* - [farrufat@cimne.upc.edu](mailto:farrufat@cimne.upc.edu)
* **Guillermo Casas** - *Fluid coupling* - [gcasas@cimne.upc.edu](mailto:gcasas@cimne.upc.edu)
* **Joaquín Irazabal** - *Particle clusters & DEM-Solid interaction* - [jirazabal@cimne.upc.edu](mailto:jirazabal@cimne.upc.edu)

## License

The DEM application is OPEN SOURCE. The main code and program structure is available and aimed to grow with the need of any user willing to expand it. The BSD (Berkeley Software Distribution) licence allows to use and distribute the existing code without any restriction, but with the possibility to develop new parts of the code on an open or close basis depending on the developers.

## FAQ

### What to do if particles behave strangely

* Check the Young Modulus. Materials with low elasticity may require smaller timesteps to ensure stability.
* Check the material density.
* Check the Time Step. If the time step is too big, the elements can fail to interact with each other. In worst case scenarios it may even crash the program.
* Check the frequency of neighbour search. If the search is not done frequently enough, new contacts may not be detected.
* Check the restitution coefficient of the material. Explicit integration schemes gain energy noticeably, unless you use a really small time step. If the time step is high (but stable), use the restitution coefficient to compensate the gain of energy and get more realistic results.
