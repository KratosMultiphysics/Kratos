# Policies

## Introduction

In our element structure, the strategy-pattern (also known as policy-pattern), is used to gain control and
configurability for the element functionality. This is done by encapsulating the behavior in an interface and using a
pointer to the interface in the element. The behavior can then be changed by changing the pointer to a different
implementation of the interface. This is a powerful tool for creating flexible and reusable code.

Currently, the following policies are implemented (still under development):

1. The [StressStatePolicy](#stressstatepolicy) - this is used to calculate several stress-related properties depending
   on the configuration. There are three flavours: `ThreeDimensionalStressState`, `PlaneStrainStressState` and `AxisymmetricStressState`.

## StressStatePolicy

The stress state policy is used to easily configure elements to have a specific stress state. The responsibility of the
stress state is to convert strain tensors to [strain vectors](#strain-vectors), define the strain-displacement relation
by filling the [B-matrix](#b-matrix) and calculate the [integration coefficient](#integration-coefficient) (which needs a correction for the axisymmetric stress state).

There are three different
stress state implementations:

1. `ThreeDimensionalStressState`
2. `PlaneStrainStressState`
3. `AxisymmetricStressState`

They all derive from the same interface (`StressStatePolicy`) and can be used interchangeably to calculate the B-matrix,
the integration coefficient and the Green Lagrange Strain vector. See the following schematic for the structure of the
stress state policies:

![stress_state_policies.svg](stress_state_policies.svg)

For simple code examples of the functionalities described in the next sections, we refer to the unit tests for [3d](../tests/cpp_tests/test_three_dimensional_stress_state.cpp), [plane strain](../tests/cpp_tests/test_plane_strain_stress_state.cpp) and [axisymmetric](../tests/cpp_tests/test_axisymmetric_stress_state.cpp) stress states.

### Strain vectors

For the different stress states, the strain vector is created by performing a conversion of the strain tensor, which is defined as:
```math
\epsilon = \begin{bmatrix} \epsilon_{xx} & \epsilon_{xy} & \epsilon_{xz} \\
                    \epsilon_{yx} & \epsilon_{yy} & \epsilon_{yz} \\
                    \epsilon_{zx} & \epsilon_{zy} & \epsilon_{zz} \end{bmatrix}
```
The strain tensor can be calculated using different methods, like Green Lagrange or Hencky Strain. This calculation is delegated to the `StressStrainUtilities`. However the conversion to strain vectors depends on the type of stress state. Therefore, that conversion is defined in these policies and are described below. 

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

Lastly, for the axisymmetric stress state, it is defined as:
```math
\vec{\epsilon} = \begin{bmatrix} \epsilon_{xx} \\
                                 \epsilon_{yy} \\
                                 u_r/r \\
                                 \epsilon_{xy} \end{bmatrix}
```
Where $u_r$ is the radial displacement and $r$ is the radial coordinate. 

### B-matrix
The B matrix is used to relate strains and displacements. Therefore, its elements are filled with the spatial gradients of the shape functions ($N$).

For a 3D element consisting of 4 nodes, the gradient looks like:
```math
\nabla{N} =
\begin{bmatrix}
\delta N_1/\delta x & \delta N_1/\delta y & \delta N_1/\delta z\\
\delta N_2/\delta x & \delta N_2/\delta y & \delta N_2/\delta z\\
\delta N_3/\delta x & \delta N_3/\delta y & \delta N_3/\delta z\\
\delta N_4/\delta x & \delta N_4/\delta y & \delta N_4/\delta z
\end{bmatrix}
```
Resulting in the following B-matrix:
```math
B =
\begin{bmatrix}
\delta N_1/\delta d_x & 0 & 0 & && & \delta N_4/\delta d_x & 0 & 0 \\
0 & \delta N_1/\delta d_y & 0 & && & 0 & \delta N_4/\delta d_y & 0\\
0 & 0 & \delta N_1/\delta d_z & &\dots& & 0 & 0 & \delta N_4/\delta d_z\\
\delta N_1/\delta d_y & \delta N_1/\delta d_x & 0 & && & \delta N_4/\delta d_y & \delta N_4/\delta d_x & 0 \\
0 & \delta N_1/\delta d_z & \delta N_1/\delta d_y & && & 0 & \delta N_4/\delta d_z & \delta N_4/\delta d_y \\
\delta N_1/\delta d_z & 0 & \delta N_1/\delta d_x & && & \delta N_4/\delta d_z & 0 & \delta N_4/\delta d_x
\end{bmatrix}
```
Where the first three columns are for the first node, the next three columns for the second node and so on for the number of nodes in the element.

For a 2D element consisting of three nodes, the gradient looks like:
```math
\nabla{N} =
\begin{bmatrix}
\delta N_1/\delta d_x & \delta N_1/\delta d_y \\
\delta N_2/\delta d_x & \delta N_2/\delta d_y \\
\delta N_3/\delta d_x & \delta N_3/\delta d_y
\end{bmatrix}
```
Where the rows depict the different nodes in the element (i.e. $N_i$ is the shape function for node $i$). For the plane strain case, this would result in the following B-matrix
```math
B =
\begin{bmatrix}
\delta N_1/\delta d_x & 0 & \delta N_2/\delta d_x & 0 & \delta N_3/\delta d_x & 0 \\
0 & \delta N_1/\delta d_y & 0 & \delta N_2/\delta d_y & 0 & \delta N_3/\delta d_y\\
0 & 0 & 0 & 0 & 0 & 0\\
\delta N_1/\delta d_y & \delta N_1/\delta d_x & \delta N_2/\delta d_y & \delta N_2/\delta d_x & \delta N_2/\delta d_y & \delta N_3/\delta d_x
\end{bmatrix}
```
For the axisymmetric stress state, this would result in the following B-matrix
```math
B =
\begin{bmatrix}
\delta N_1/\delta d_x & 0 & \delta N_2/\delta d_x & 0 & \delta N_3/\delta d_x & 0 \\
0 & \delta N_1/\delta d_y & 0 & \delta N_2/\delta d_y & 0 & \delta N_3/\delta d_y\\
N_1 / r & 0 & N_2 / r & 0 & N_3 / r & 0\\
\delta N_1/\delta d_y & \delta N_1/\delta d_x & \delta N_2/\delta d_y & \delta N_2/\delta d_x & \delta N_2/\delta d_y & \delta N_3/\delta d_x
\end{bmatrix}
```

### Integration coefficient

For the 3D and plane strain stress states, the integration coefficient (or integration volume) is the product of the weight of the
integration point and the determinant of the Jacobian:
$$V_i = w_i \det{(J)}$$
Where $w_i$ is the weight and $J$ is the Jacobian. For the axisymmetric stress state, this same calculation is done, but due to the symmetry, the integration volume
is corrected by multiplying it with the axisymmetric circumference:
$$V_i = 2\pi r w_i \det{(J)}$$
