# Test Cases for Mohr-Coulomb with tension cutoff

**Author:** [Mohamed Nabi](https://github.com/mnabideltares)

**Source files:** [Thermal fixed temperature](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/GeoMechanicsApplication/tests/test_mohr_coulomb_with_tension_cutoff)

## Case Specification
In these test cases, a 1 x 1 m soil is considered. The pressure is fixed to be at 0 Pa and only the deformation is taken into account. The size of this block is changed by means of Dirichlet condition, and the stresses are then calculated. The purpose of these tests is to validate the response of different zones in the Mohr-Coulomb diagram (including the tension cutoff). 

The diagram includes one admissible and four inadmissible zones. Hence, five test cases are considered. Based on the imposed boundary conditions the trial stresses fall in one of those zones. These tests consider one quadrilateral element of 4 nodes, namely 2D4N.

The boundary conditions are shown below:

<img src="MeshStructure.svg" width="600">

All degrees of freedom have been prescribed. In horizontal and vertical directions, the element is compressed or extended  by imposing displacements. The elastic behavior of the elements is described by a Young's modulus of $`10o0 \mathrm{Pa}`$ and a Poisson's ratio of $`0.0`$. The material behavior is described using a cohesion of $`10.0 \mathrm{Pa}`$, a friction angle $`35^{\circ}`$, a dilatancy angle of $`20^{\circ}`$ and a tensile strength of $`10 \mathrm{Pa}`$.
