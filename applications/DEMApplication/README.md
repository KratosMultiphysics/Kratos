# DEM Application

This application focuses on the Discrete Element Method (DEM), a particles method for modeling the bulk behavior of granular materials and many geomaterials such as coal, ores, soil, rocks, aggregates, pellets, tablets and powders.

The [DEMpack Team](www.cimne.com/dem) at [CIMNE](www.cimne.com) is in charge of all developments related to the DEM.

For the coupling between DEM and Fluids, go to the [Swimming DEM Application](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/SwimmingDEMApplication).

For the coupling between DEM and thermal effects, go to the [Thermal DEM Application](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/ThermalDEMApplication).

## Getting started

This application is part of the Kratos Multiphysics Platform. Instructions on how to get you a copy of the project up and running on your local machine for development and testing purposes are available for both [Linux](http://kratos-wiki.cimne.upc.edu/index.php/LinuxInstall) and [Windows](http://kratos-wiki.cimne.upc.edu/index.php/Windows_7_Download_and_Installation) systems.

### Prerequisites

Build [Kratos](https://github.com/KratosMultiphysics/Kratos/wiki) and, before that, make sure that you add

``` cmake
-DDEM_APPLICATION=ON
```

amongst the compilation options, so the DEM application is compiled.

No auxiliar external libraries are needed.

## Theory

The DEM is a numerical method that has been applied to simulate and analyze flow behavior in a wide range of disciplines including mechanical and process engineering, pharmaceutical, materials science, agricultural engineering and more.
Coupling with fluid is already available through the Swimming-DEM application, also integrated in the Kratos Multiphysics Platform.

The fundamental theoretical background corresponding to the discontinuous (granular matter) part of the code can be found in the DEM literature easily.

### Contact laws

The contact laws are implemented in [this folder](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/DEMApplication/custom_constitutive). Note that the letter 'D' or 'd' in the file name stands for 'discontiuum'. It is related to non cohesive or slightly cohesive contacts.

#### Linear repulsive force

The most simple representation of a repulsive contact force between a sphere and a wall is given by a linear law, where the force acting on the sphere when contacting the wall is a linear function of the indentation, which in turn would bring a quadratic dependence with the contact radius.

#### Non-Linear repulsive force

Hertz solved in 1882 the non-cohesive normal contact between a sphere and a plane. In 1971 Johnson, Kendall and Roberts presented the solution (JKR-Theory) for the same problem in this case adding cohesive behaviour. Not much later, Derjaguin, Müller and Toporov published similar results (DMT-Theory).
Both theories are very close and correct and, while the JKR theory is adequate to the study of flexible, large spheres, the DMT theory is specially suited to represent the behaviour of rigid, small ones.

## Numerical approach

The application includes two types of DEM elements used for different purposes:

* Spheric Particle - Base element used to simulate granular materials (non cohesive or slightly cohesive)
* Spheric Continuum Particle - With specific built-in variables to simulate fracture in cohesive materials. It can also be understood as a discretization method of the continuum by using spheres.

And has the following easy-to-use capabilities:

* Interaction with FEM-based walls - Objects that cannot be crossed by DEM spheres. The user can choose to impose Linear-periodic conditions or rigid body conditions.
* Inlets - Inject new particles while running the simulation linked to some material properties. With user defined granulometry, mass flow and particle type (single particle or clusters). Inlets are based on FEM-based walls and boundary conditions can also be applied to them.
* Initial conditions on particle elements.
* Boundary conditions on particle elements.

It also includes several predefined cluster formations to be used.

### DEM strategies

#### Non-cohesive materials strategy

Once contact between two spheres occurs, the forces at the contact point are computed. The interaction between the two contacting spheres can be represented by two forces with the same module but opposite directions. This force F can be decomposed into its normal and shear components Fn and Fs, respectively.
The contact interface for the simplest formulation is characterized by the normal and tangential stiffness Kn and Ks, respectively, a frictional device obeying the Coulomb law with a frictional coefficient, and a dashpot defined by a contact damping coefficient.

In order to represent irregular particles with spheres, a numerical correction is used. The rolling friction imposes a virtual moment opposite to particle rotation and dependent on its size.

#### Continuum materials strategy

For continuum materials simulations, the contact between particles can resist tractions up to a certain limit, when breakage occurs. Depending on the chosen constitutive law, the computation of the forces changes. In the basic versions, a bond strategy is used, but more advanced laws use a non-local stress-tensor based strategy.

### DEM integration schemes

The standard translational and rotational equations for the motion of rigid bodies are used to compute the dynamics of the spheres and clusters. The following schemes can be chosen separately for translation and rotation:

* Symplectic Euler
* Velocity Verlet
* Forward Euler
* Taylor

Also, two rotational specific integration schemes are available:
* [Runge-Kutta](https://link.springer.com/article/10.1007/s40571-019-00232-5)
* Quaternion based

### Contact search

The contact detection basically consists in determining, for every target object, which other objects are in contact with it, and then apply the corresponding interaction. It is usually not needed to perform a search at every time step, which is generally limited by the stability of the explicit integration of the equations of motion.
A bins based technique is used for this purpose.

## Available interfaces

### DEM

This is the package that allows a user to create, run and analyze results of a DEM simulation for discontinuum / granular / little-cohesive materials. Requires [GiD](https://www.gidhome.com/) - Pre and Post Processing software. It has both 2D and 3D versions. Check the manuals, follow the tutorials or play with the preloaded sample problems in order to learn how this application works.

### Cohesive-DEM

This package combines the features of the previous one also with the simulation of continuum/cohesive materials. It also offers the possibility of tackling both 2D and 3D problems. Check also the manuals or tutorials or load the test examples in the GUI in order to learn how this problem type works.

### Fluid-DEM

This package allows you to simulate a wide spectrum of problems involving the interaction of granular DEM and fluids. This application has only a 3D version. Check also for existing manuals or tutorials to get a feel of how to work with this application.


## Contact

* **Miguel Angel Celigueta** - *Core Development* - [maceli@cimne.upc.edu](mailto:maceli@cimne.upc.edu)
* **Salva Latorre** - *Granular materials* - [latorre@cimne.upc.edu](mailto:latorre@cimne.upc.edu)
* **Ferran Arrufat** - *Cohesive materials* - [farrufat@cimne.upc.edu](mailto:farrufat@cimne.upc.edu)
* **Guillermo Casas** - *Fluid coupling* - [gcasas@cimne.upc.edu](mailto:gcasas@cimne.upc.edu)
* **Joaquín Irazabal** - *Particle clusters & DEM-Solid interaction* - [jirazabal@cimne.upc.edu](mailto:jirazabal@cimne.upc.edu)
* **Joaquín González-Usúa** - *Fluid coupling* - [jgonzalez@cimne.upc.edu](mailto:jgonzalez@cimne.upc.edu)
* **Chengshun Shang** - *Bonded partcile models* - [cshang@cimne.upc.edu](mailto:cshang@cimne.upc.edu)


## License

The DEM application is OPEN SOURCE. The main code and program structure is available and aimed to grow with the need of any users willing to expand it. The BSD (Berkeley Software Distribution) licence allows to use and distribute the existing code without any restriction, but with the possibility to develop new parts of the code on an open or close basis depending on the developers.


## New GIDInterface for Kratos

The new GIDInterface currently under development can be found [here](https://github.com/KratosMultiphysics/GiDInterface). Based on the customLib, it includes the interfaces for most of the Kratos applications in addition to the new DEM interface.


## FAQ

### What to do if particles behave strangely

* Check the Young Modulus. Materials with high stiffness may require smaller time steps to ensure stability.
* Check the material density.
* Check the time step. If the time step is too large, the elements can fail to interact with each other. In the worst case scenarios, the simulation may even crash.
* Check the frequency of neighbours' search. If the search is not done frequently enough, new contacts may not be detected.
* Check the restitution coefficient of the material. Explicit integration schemes gain energy noticeably unless an enough small time step is used. If the time step is large (but stable), use the restitution coefficient to compensate for the gain of energy to obtain more realistic results.
