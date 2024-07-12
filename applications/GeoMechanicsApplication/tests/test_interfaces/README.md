# Interface

This set of tests verifies the interface functionality.


## Using ingerface between beam and soil

### Setup

This test models the beam movement in the soil. The beam is fixed on the bottom and loaded on the top. The model consists of 4230 3-noded elements for the soil (using the UPwSmallStrainElement2D3N class), 31 2-noded elements for the beam (GeoCrBeamElementLinear2D2N), and 62 4-noded elements for interfaces (UPwSmallStrainInterfaceElement2D4N) that connect the beam with the soil. A schematic representation can be found in the figure below:

![MeshStructure](interface_on_beam.svg)

All nodes on the sides have fixed horizontal displacements. The bottom nodes have fixed displacement in both the horizontal and the vertical
direction. The load on the bam top has a magnitude of $q=10.0 \mathrm{[N/m]}$. This load is kept constant during the whole analysis time. 

In the middel of the domain a vertical beam is placed. It is shown with red color.


Soil is decribed with GeoLinearElasticPlaneStrain2DLaw:
-   A Young's modulus $E = 3.0E7 \mathrm{[Pa]}$ with Poisson's ratio $\nu = 0.20 \mathrm{[-]}$.
-   The soil and water density are $2000$ and $1000 \mathrm{[kg/m^3]}$ respectively. The porosity is $n=0.3$. 
-   The bulk modulus of solid $K = 1.0e12 \mathrm{[Pa]}$.
-   The dynamic viscosity of water is given as $\mu = 1.0E-3 \mathrm{[Pa \cdot s]}$ and the intrinsic permeability of the soil as $\kappa = 4.5e-30 \mathrm{[m^2]}$.


The beam is described with BeamConstitutiveLaw:
-   A Young's modulus $E = 2.07E+13 \mathrm{[Pa]}$ with Poisson's ratio $\nu = 0.29 \mathrm{[-]}$.
-   Density of $7850.0 \mathrm{[kg/m^3]}$.
-   Cross area of $0.01 \mathrm{[m^2]}$,
-   Moment of interia about Z axis $I33 = 8.33333E-08 \mathrm{[kg m^2]}$

The interface is described with SmallStrainUDSM2DInterfaceLaw and it has a cohesion of $1000 \mathrm{[kN]}$ and a stiffness of $1e12 \mathrm{[Pa]}$. The stiffness is provided as the first UMAT parameter. Where is the cohesion defined? Inverse of MINIMUM_JOINT_WIDTH value? This value affects the week interface test.

### Solution

Under the load the beam bends and acts on the soil. The following picture shows the soil displacement.

![Displacement](interface_on_beam_deformation.svg)

### Assertions

The test asserts max displacement in X and Y directions between this case and the base case. The base case has the same settings but it does not use the interface. 

## Weak Interface

### Setup

This test uses the same setting as the previous test except the interface stiffness of $1e2 \mathrm{[Pa]}$. 

### Solution

The very stiffness value leads to a large movement of the beam against the soil. The following picture shows that the beam moves as before but the soil experiences very little effect. The brown lines show the interfaces' horizontal lines.

![Displacement](week-interface_on_beam_deformation.svg)

### Assertions

The test checks that the beam displacement shall be much larger than the soil displacement. 

## Interface cohesive side

This test models a vertical movement of the soil block which left hand side is connected to the interface. 

### Setup

The model consists of 275 3-noded elements for the soil (using the UPwSmallStrainElement2D3N class) and 21 4-noded elements for interface (UPwSmallStrainInterfaceElement2D4N). Al. A schematic representation can be found in the figure below, the interface is shown with a pink color. 

![MeshStructure](box.svg)

The soil is dry and it is decribed with GeoLinearElasticPlaneStrain2DLaw:
-   A Young's modulus $E = 1.0E9 \mathrm{[Pa]}$ with Poisson's ratio $\nu = 0.49 \mathrm{[-]}$.
-   The soil density is $2000 \mathrm{[kg/m^3]}$ and the porosity is $n=0.3$. 
-   The bulk modulus of solid $K = 1.0e12 \mathrm{[Pa]}$.

The interface is described with SmallStrainUDSM2DInterfaceLaw and it has a cohesion of $1000 \mathrm{[kN]}$ and a stiffness of $1e12 \mathrm{[Pa]}$. 

The left side of the interface is fixed. The soil moves down with a prescribed vertical displacement of $-0.1 \mathrm{[m]}$  and a horizontal line load of $-1 \mathrm{[kN]}$  is applied to the soil right side. 

### Solution

A result of the prescribed displacement is shown below. The soil moved down, when the interface is sqeezed because its left side is fixed.  

![Displacement](box-moved.svg)


### Assertions

The test checks the expected shear stress in the interface of $10 \mathrm{[kN]}$. 
