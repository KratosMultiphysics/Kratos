# Interface

This set of tests verifies the interface functionality.

## Interface properties' input

Currently, User Defined Soil Model (UDSM) is used to prescribe the interface properties. However, the standard properties like soil/water properties, porosity etc. must be prescribed as well when their values can be over-written by UDSM. In the following tests only two UMAT parameters are used, the first one is a stiffness and the third one is a cohesion.   

## Using interfaces between beam and soil

### Setup

This test models the beam movement in the soil. The model consists of 4230 3-noded elements for the soil (using the UPwSmallStrainElement2D3N class), 31 2-noded elements for the beam (GeoCrBeamElementLinear2D2N), and 62 4-noded elements for interfaces (UPwSmallStrainInterfaceElement2D4N) that connect the beam with the soil. A schematic representation can be found in the figure below:

![MeshStructure](interface_on_beam.svg)

All nodes on the vertical sides and the bottom are fixed in both the horizontal and the vertical directions. In the middle of the domain the vertical beam is placed. It is shown in red color and letters B and T are placed near it bottom and top, respectively. The beam bottom has a fixed support and a constant load in the positive x-direction with a magnitude of $q=10.0 \  \mathrm{[N/m]}$ is applied to the beam top. 


The soil is described with GeoLinearElasticPlaneStrain2DLaw:
-   A Young's modulus $E = 3.0e7 \  \mathrm{[Pa]}$ with Poisson's ratio $\nu = 0.20 \  \mathrm{[-]}$.
-   The soil and water density are $2000$ and $1000 \  \mathrm{[kg/m^3]}$ respectively. The porosity is $n=0.3$. 
-   The bulk modulus of solid $K = 1.0e12 \  \mathrm{[Pa]}$.
-   The dynamic viscosity of water $\mu = 10^{-3} \  \mathrm{[Pa \cdot s]}$ and the intrinsic permeability of the soil $\kappa = 4.5\cdot 10^{-30} \  \mathrm{[m^2]}$.


The beam is described with BeamConstitutiveLaw:
-   A Young's modulus $E = 2.07e13 \  \mathrm{[Pa]}$ with Poisson's ratio $\nu = 0.29 \  \mathrm{[-]}$.
-   Density of $7850.0 \  \mathrm{[kg/m^3]}$.
-   Cross area of $0.01 \  \mathrm{[m^2]}$,
-   Moment of inertia about Z axis $I33 = 8.33333\cdot 10^{-8} \  \mathrm{[kg m^2]}$

The interface is described with SmallStrainUDSM2DInterfaceLaw and it has a stiffness of $1e12 \  \mathrm{[Pa]}$.  

### Solution

Under the load the beam bends and acts on the soil. The following picture shows the bent beam and soil displacement.

![Displacement](interface_on_beam_deformation.svg)

### Assertions

The test asserts maximum values of displacement in X and Y directions. The values are compared with the correspondent solution obtained for a base case. The base case has the same settings, but it does not use the interface. 

## Weak Interface

### Setup

This test uses the same setting as the previous test except the interface stiffness is set to $1e2 \  \mathrm{[Pa]}$ (instead of $1e12 \  \mathrm{[Pa]}$). 

### Solution

This very small value of the stiffness leads to a large movement of the beam against the soil. The following picture shows that the beam is bent as in the previous test but the beam affects the soil displacement very little. As a result the interface elements are stretched. The brown lines depict the interfaces' horizontal lines that connect the interface soil part with interface beam part. 

![Displacement](weak_interface_on_beam_deformation.svg)

### Assertions

The test checks that the beam displacement shall be much larger than the soil displacement (a factor of $10^8$). 

## Interface cohesive side

This test models a vertical movement of the entire soil block which left hand side is connected to the interface. 

### Setup

The model consists of 275 3-noded elements for the soil (using the UPwSmallStrainElement2D3N class) and 21 4-noded elements for interface (UPwSmallStrainInterfaceElement2D4N). A schematic representation can be found in the figure below, the interface is shown with a pink color. 

![MeshStructure](box.svg)

The soil is dry, and it is described with GeoLinearElasticPlaneStrain2DLaw:
-   A Young's modulus $E = 1.0E9 \  \mathrm{[Pa]}$ with Poisson's ratio $\nu = 0.49 \  \mathrm{[-]}$.
-   The soil density is $2000 \  \mathrm{[kg/m^3]}$ and the porosity is $n=0.3$. 
-   The bulk modulus of solid $K = 1.0e12 \  \mathrm{[Pa]}$.

The interface is described with SmallStrainUDSM2DInterfaceLaw and it has a cohesion of $10 \  \mathrm{[kN]}$ and a stiffness of $1e10 \  \mathrm{[Pa]}$. A higher value of the stiffness leads to a solution divergence for the current setup.

The left side of the interface is fixed. The upper and bottom sides of the soil block move down with a prescribed vertical displacement of $-0.1 \  \mathrm{[m]}$  and a horizontal line load of $-1000 \  \mathrm{[kN]}$  is applied to the soil block right side. 

### Solution

A result of the prescribed displacement is shown below. The soil moved down, when the interface is squeezed because its left side is fixed.  

![Displacement](box-moved.svg)


### Assertions

The test checks the expected shear stress in the interface of $10 \  \mathrm{[kN]}$. 
