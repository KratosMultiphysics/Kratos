# Policies

## Introduction

In our element structure, the strategy-pattern (also known as policy-pattern), is used to gain control and
configurability for the element functionality. This is done by encapsulating the behavior in an interface and using a
pointer to the interface in the element. The behavior can then be changed by changing the pointer to a different
implementation of the interface. This is a powerful tool for creating flexible and reusable code. For more information
about the strategy-pattern, the reader is referred to the book written by 'the Gang of four' on the topic
of [design patterns](#bibliography).

Currently, the following policies are implemented (still under development):

1. The [StressStatePolicy](#stressstatepolicy) - this is used to calculate several stress-related properties depending
   on the configuration.
2. The [Integration coefficient calculator](#integration-coefficient-calculator)

## StressStatePolicy

The stress state policy is used to easily configure elements to have a specific stress state. The responsibility of the
stress state is to convert strain tensors to [strain vectors](#strain-vectors), define the strain-displacement relation
by filling the [B-matrix](#b-matrix).

Currently, there are five different stress state implementations:

1. `ThreeDimensionalStressState`
2. `PlaneStrainStressState`
3. `AxisymmetricStressState`
4. `Line2DInterfaceStressState`
5. `SurfaceInterfaceStressState`

They all derive from the same interface (`StressStatePolicy`) and can be used interchangeably to calculate the B-matrix and the Green Lagrange Strain vector. See the following schematic for the structure of the
stress state policies:

![stress_state_policies.svg](stress_state_policies.svg)

For simple code examples of the functionalities described in the next sections, we refer to the unit tests for [3D](../tests/cpp_tests/custom_elements/test_three_dimensional_stress_state.cpp), [plane strain](../tests/cpp_tests/custom_elements/test_plane_strain_stress_state.cpp), [axisymmetric](../tests/cpp_tests/custom_elements/test_axisymmetric_stress_state.cpp), and [interface](../tests/cpp_tests/custom_elements/test_interface_stress_state_policy.cpp) stress states.

### Strain vectors

For the different stress states, the strain vector is created by performing a conversion of the strain tensor, which is defined as:
```math
\boldsymbol{\epsilon} = \begin{bmatrix} \epsilon_{xx} & \epsilon_{xy} & \epsilon_{xz} \\
                        \epsilon_{yx} & \epsilon_{yy} & \epsilon_{yz} \\
                        \epsilon_{zx} & \epsilon_{zy} & \epsilon_{zz} \end{bmatrix}
```
Note that the strain tensor is symmetric, meaning that $\epsilon_{ij} = \epsilon_{ji}$. The strain tensor can be calculated using different methods, like Green Lagrange or Hencky Strain. This calculation is delegated to the `StressStrainUtilities`. However, the conversion to strain vectors depends on the type of stress state. Therefore, that conversion is defined in these policies and are described below. 

For the 3D stress state, the strain vector is defined as:
```math
\vec{\epsilon} = \begin{bmatrix} \epsilon_{xx} \\
                                 \epsilon_{yy} \\
                                 \epsilon_{zz} \\
                                 \epsilon_{xy} \\
                                 \epsilon_{yz} \\
                                 \epsilon_{xz} \end{bmatrix}
```

For the plane strain stress state, it is defined as:
```math
\vec{\epsilon} = \begin{bmatrix} \epsilon_{xx} \\
                                 \epsilon_{yy} \\
                                 0 \\
                                 \epsilon_{xy} \end{bmatrix}
```

For the axisymmetric stress state, it is defined as:
```math
\vec{\epsilon} = \begin{bmatrix} \epsilon_{xx} \\
                                 \epsilon_{yy} \\
                                 u_r/r \\
                                 \epsilon_{xy} \end{bmatrix}
```
Where $u_r$ is the radial displacement and $r$ is the radial coordinate (in our geomechanics code base, the radial coordinate is equal to $x$).

For the interface stress states, it is _not_ possible to calculate the Green-Lagrange strain based on the deformation gradient.  Interface elements use relative displacements as strain measure.  For two-dimensional line interface elements, the relative displacements are:
```math
\vec{\Delta u} = \begin{bmatrix}\Delta u_n \\
                                \Delta u_t \end{bmatrix}
```
and for surface interface elements, the relative displacements are:
```math
\vec{\Delta u} = \begin{bmatrix}\Delta u_n \\
                                \Delta u_t \\
                                \Delta u_s \end{bmatrix}
```
Here, $`\Delta u_n`$ is the relative displacement component in normal direction, and $`\Delta u_t`$ and $`\Delta u_s`$ are the relative displacement components in the tangential directions. 

### B-matrix
The B-matrix is used to relate strains and displacements. Therefore, its elements are filled with the spatial gradients of the shape functions ($N$).

For a 3D element consisting of 4 nodes, the gradient looks like:
```math
\nabla{N} =
\begin{bmatrix}
\partial N_1/\partial x & \partial N_1/\partial y & \partial N_1/\partial z\\
\partial N_2/\partial x & \partial N_2/\partial y & \partial N_2/\partial z\\
\partial N_3/\partial x & \partial N_3/\partial y & \partial N_3/\partial z\\
\partial N_4/\partial x & \partial N_4/\partial y & \partial N_4/\partial z
\end{bmatrix}
```
Resulting in the following B-matrix:
```math
B =
\begin{bmatrix}
\partial N_1/\partial x & 0 & 0 & && & \partial N_4/\partial x & 0 & 0 \\
0 & \partial N_1/\partial y & 0 & && & 0 & \partial N_4/\partial y & 0\\
0 & 0 & \partial N_1/\partial z & &\dots& & 0 & 0 & \partial N_4/\partial z\\
\partial N_1/\partial y & \partial N_1/\partial x & 0 & && & \partial N_4/\partial y & \partial N_4/\partial x & 0 \\
0 & \partial N_1/\partial z & \partial N_1/\partial y & && & 0 & \partial N_4/\partial z & \partial N_4/\partial y \\
\partial N_1/\partial z & 0 & \partial N_1/\partial x & && & \partial N_4/\partial z & 0 & \partial N_4/\partial x
\end{bmatrix}
```
Where the first three columns are for the first node, the next three columns for the second node and so on for the number of nodes in the element.

For a 2D element consisting of three nodes, the gradient looks like:
```math
\nabla{N} =
\begin{bmatrix}
\partial N_1/\partial x & \partial N_1/\partial y \\
\partial N_2/\partial x & \partial N_2/\partial y \\
\partial N_3/\partial x & \partial N_3/\partial y
\end{bmatrix}
```
Where the rows depict the different nodes in the element (i.e. $N_i$ is the shape function for node $i$). For the plane strain case, this would result in the following B-matrix:
```math
B =
\begin{bmatrix}
\partial N_1/\partial x & 0 & \partial N_2/\partial x & 0 & \partial N_3/\partial x & 0 \\
0 & \partial N_1/\partial y & 0 & \partial N_2/\partial y & 0 & \partial N_3/\partial y\\
0 & 0 & 0 & 0 & 0 & 0\\
\partial N_1/\partial y & \partial N_1/\partial x & \partial N_2/\partial y & \partial N_2/\partial x & \partial N_3/\partial y & \partial N_3/\partial x
\end{bmatrix}
```
For the axisymmetric stress state, this would result in the following B-matrix:
```math
B =
\begin{bmatrix}
\partial N_1/\partial x & 0 & \partial N_2/\partial x & 0 & \partial N_3/\partial x & 0 \\
0 & \partial N_1/\partial y & 0 & \partial N_2/\partial y & 0 & \partial N_3/\partial y\\
N_1 / r & 0 & N_2 / r & 0 & N_3 / r & 0\\
\partial N_1/\partial y & \partial N_1/\partial x & \partial N_2/\partial y & \partial N_2/\partial x & \partial N_3/\partial y & \partial N_3/\partial x
\end{bmatrix}
```
Note that in our geomechanics code base, the radial coordinate $r$ is equal to $x$.

For interface stress states, the $`B`$ matrix relates the nodal displacement vector $`u`$ (with length "number of nodes" times "number of displacement degrees of freedom per node") of an interface element to the relative displacement vector $`\Delta u`$ at an integration point (with length 2 or 3, depending on the exact interface configuration): $`\Delta u = B \cdot u`$.   In essence, to calculate a relative displacement component at an integration point we need to subtract the relevant displacement component at the first side from the corresponding one at the second side.  Given the shape function values $`N_i`$ at the integration point and the nodal displacement vector $`u`$ of the element, we can calculate these pairs of displacement components.  In addition, the $`B`$ matrix needs to take into account that the relative displacement(s) in normal direction (which precede(s) the one(s) in tangential direction) depends on degrees of freedom that are listed _after_ the degrees of freedom that are needed for the relative displacement(s) in tangential direction.  For instance, for a 2+2 two-dimensional line interface, the $`B`$ matrix looks as follows:
```math
B =
\begin{bmatrix}
 0   & -N_1 &  0   & -N_2 & 0   & N_1 & 0   & N_2 \\
-N_1 &  0   & -N_2 &  0   & N_1 & 0   & N_2 & 0
\end{bmatrix}
```
and for a 3+3 surface interface, the $`B`$ matrix looks as follows:
```math
B =
\begin{bmatrix}
 0   &  0   & -N_1 &  0   &  0   & -N_2 &  0   &  0   & -N_3 & 0   & 0   & N_1 & 0   & 0   & N_2 & 0   & 0   & N_3 \\
-N_1 &  0   & 0    & -N_2 &  0   &  0   & -N_3 &  0   &  0   & N_1 & 0   & 0   & N_2 & 0   & 0   & N_3 & 0   & 0   \\
 0   & -N_1 & 0    & 0    & -N_2 &  0   &  0   & -N_3 &  0   & 0   & N_1 & 0   & 0   & N_2 & 0   & 0   & N_3 & 0
\end{bmatrix}
```
The generalization here is how the components of the relative displacement vector (normal and tangential components) relate to the displacement degrees of freedom at the nodes.  This differs depending on the exact interface element configuration. 

## Integration coefficient calculator

Based on its name, the integration coefficient calculator calculates the integration coefficient. The integration coefficient (or integration volume) 
depends on the integration weight ($w_i$) and the determinant of the Jacobian ($J$). 

### Default calculation
There is a default version that works for plane and 3D elements. In case of the 3D element, the integration volume is calculated as:
$$V_i = w_i \det{(J)}$$

Since $w_i$ is the volume in isoparametric coordinates. For the plane element, $w_i$ is the area in isoparametric coordinates. 
However, since the thickness is by definition 1, the integration coefficient can be calculated using the same equation.

### Modifiers

For other elements the integration coefficient calculator has a modifier, which modifies the default integration coefficient value. Currently, there are two modifiers and they are described below. 

#### An axisymmetric element

This modifier takes $w_i$ as the isoparametric area and converts it to a volume by a multiplication with the axisymmetric circumference:
$$V_i = 2\pi r w_i \det{(J)}$$

#### A line element

There is a number of line elements and the line element modifier provides the default integration coefficient multiplied by the user provided CROSS_AREA.


## Bibliography
Gamma, E., Helm, R., Johnson, R., & Vlissides, J. (1994). _Design Patterns: Elements of Reusable Object-Oriented Software._
