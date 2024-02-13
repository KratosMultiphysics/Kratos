# Lysmer absorbing boundary conditions on a column made of 3D elements

These 3D tests verify the P-wave propagation through 1D column.  

## Setup

The tests is performed for an 1D column. The column height and width/depth are 10 $m$ and 1 $m$, respectively. A
schematic representation can be found in the figure below:

![SetUp](SetUp.svg)

All Nodes have fixed displacement in the lateral directions. The column bottom is fixed in all directions. The absorbing boundary conditions are applied at the bottom and the unity absorbing factors are used there. An instant load $L$ of 10 $N/m^2$ is applied at the top of the column towards the bottom.

Two kinds of 3D elements are used to generate a mesh:

-   tetra elements: folders **tetra** and **tetra_horizontal**
-   hexa elements: folder **hexa**

The mesh with tetra elements is presented below for both cases. Note, the column is aligned along Z axis for **tetra** and Y axis for **tetra_horizontal**.
![MeshTetra](MeshTetra.svg) ![MeshTetraHorizontal](MeshTetraHorizontal.svg)

The mesh with hexa elements is given in the following figure.
![MeshHexa](MeshHexa.svg)

The material is described using:

-   A linear elastic 3D model (LinearElastic3DLaw),
-   Young's modulus is $E$=10000 $kN/m^2$
-   Poisson's ratio is $\nu$=0.2,
-   Density of solid is $\rho$=2.65 $ton/m^3$,
-   Porosity of $\Phi$=0.3.

P-wave velocity is $v=\sqrt{\frac{E \frac{1 - \nu}{(1 + \nu) (1 - 2 \nu)}}{\rho (1 - \Phi)}}$

Node velocity is $L / (v \rho (1 - \Phi))$

## Assertions

The instant load generates a P-wave. The test 

-   calculates time when the P-wave arrives at Nodes 31, 54, and 81 that are at 1/4, 1/2, 3/4 of the column length from the column top and
-   checks the velocity of these nodes at that time moments. 
