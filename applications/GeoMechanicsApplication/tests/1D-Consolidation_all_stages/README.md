# 1D Consolidation

This tests verifies that in a column of triangular plane strain elements, consolidation follows the solution given by Terzaghi in 1923.

## Setup

This test consists of 78 6-noded elements (using the SmallStrainUPwDiffOrderElement2D6N class) of b x h = 0.1 x 1.0 meter. A
schematic representation can be found in the figure below:

![MeshStructure](MeshStructure.svg)

All nodes on the sides have fixed horizontal displacements. The bottom nodes have fixed displacement in both the horizontal and the vertical
direction. At the top of the column a vertical compressive line load with a magnitude of 1.0 N/m is applied. This load is kept constant during the whole analysis time. Note that gravity is absent in this computation, the top load generates the excess pore pressure. In the first stage no Dirichlet boundaries for the water pressure D.O.F. are given, such that all water is contained within the column. In all later stages, the water pressure on the top of the column is specified to be 0. [Pa]. This creates an outflow boundary at the top.

The material is described using:
-   A linear elastic plane stress model (LinearElasticPlaneStrain2DLaw)
-   A Young's modulus $E = 1.0E3 \mathrm{[Pa]}$ with Poisson's ratio $\nu = 0.0 \mathrm{[-]}$.
-   The soil and water density are specified, but irrelevant due to the absence of gravity. The porosity is 0.3.
These material properties of the dry soil give a compression modulus $K = E / (3(1-2\nu)) = 3.33E2 \mathrm{[Pa]}$ and a shear modulus $G = E / (2( 1 + \nu )) = 5.0E2 \mathrm{[Pa]}$.
-   The dynamic viscosity of water is given as $\mu = 1.0E-6 \mathrm{[Pa \cdot s]}$ and the intrinsic permeability of the soil as $\kappa = 1.17982E-15 \mathrm{[m^2]}$
-   The bulk modulus of water $K_w = 2.0E9 \mathrm{[Pa]}$

The differential equation for the 1 D consolidation column under constant load as presented by Verruijt (2001) is :

$$ \frac{\partial p}{\partial t} = c_v \frac{\partial^2 p}{\partial y^2} $$

where the consolidation coefficient:

$$ c_v = \frac{kappa}{mu ( \frac{1}{K + 4G/3} + \frac{n}{K_w})}$$

The analytical solution given by Terzaghi reads:

$$ p = \frac{4 p_0}{\pi} \SIGMA \frac{(-1)^{j-1}}{2j-1} cos((2j-1) \frac{\pi y}{2 h}) exp(-(2j-1)^2 \frac{\pi^2 c_v t}{4 h^2}) $$

Analysis stages are choosen such that at the ends of the 2nd and further stages the dimensionless time parameter $\frac{c_v t}{h^2}$ follows the pattern 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0

## Assertions

The test asserts that the vertical distribution of pressure matches the given Terzaghi solution at all nodes at the end of every stage.

## Bibliography
Verruijt, A., 2001. [Soil Mechanics](https://ocw.tudelft.nl/wp-content/uploads/SoilMechBook.pdf). Delft University of Technology.