# Policies

## Introduction

In our element structure, the strategy-pattern (also known as policy-pattern), is used to gain control and
configurability for the element functionality. This is done by encapsulating the behavior in an interface and using a
pointer to the interface in the element. The behavior can then be changed by changing the pointer to a different
implementation of the interface. This is a powerful tool for creating flexible and reusable code.

Currently, the following policies are implemented (still under development):

1. The [StressStatePolicy](#stressstatepolicy) - this is used to calculate several stress-related properties depending
   on the configuration. There are three flavours: the `AxisymmetricStressState`, the `PlaneStrainStressState` and
   the `ThreeDimensionalStressState`.

## StressStatePolicy

The stress state policy is used to easily configure elements to have a specific stress state. The responsibility of the
stress state is to convert strain tensors to [strain vectors](#strain-vectors), define the strain-displacement relation
by filling the [B-matrix](#b-matrix) and calculate the [integration coefficient](#integration-coefficient) (which needs a correction depending on the symmetry of the
stress state).

There are three different
stress state implementations:

1. `ThreeDimensionalStressState`
2. `PlaneStrainStressState`
3. `AxisymmetricStressState`

They all derive from the same interface (`StressStatePolicy`) and can be used interchangeably to calculate the B-matrix,
the integration coefficient and the Green Lagrange Strain vector. See the following schematic for the structure of the
stress state policies:

![stress_state_policies.svg](stress_state_policies.svg)

### Strain vectors

For the different stress states, the strain vector is created by performing a conversion of the strain tensor, which is defined as:
```math
E = \begin{bmatrix} \epsilon_{xx} & \epsilon_{xy} & \epsilon_{xz} \\
                    \epsilon_{yx} & \epsilon_{yy} & \epsilon_{yz} \\
                    \epsilon_{zx} & \epsilon_{zy} & \epsilon_{zz} \end{bmatrix}
```
The strain tensor can be calculated using different methods, like Green Lagrange or Hencky Strain. This calculation is delegated to the `StressStrainUtilities`. However the conversion to strain vectors depends on the type of stress state. Therefore that conversion is defined in these policies and are described below. 

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
                                 u/r \\
                                 \epsilon_{xy} \end{bmatrix}
```
Where $u$ is the displacement and $r$ is the radial coordinate. 

### B-matrix
The B matrix is used to relate strains and displacements. Therefore, its elements consist of numbers in the `rGradNpT` matrix (Gradient of the shape functions with respect to ???)

For a very simple 2D element consisting of three nodes, this gradient could look like:
```math
\nabla{N} =
\begin{bmatrix}
\delta N_1/\delta d_x & \delta N_1/\delta d_x \\
\delta N_2/\delta d_x & \delta N_2/\delta d_x \\
\delta N_3/\delta d_x & \delta N_3/\delta d_x
\end{bmatrix}
```
Where the rows depict the different nodes in the element.

### Integration coefficient

For the 3D and plane strain stress states, the integration coefficient is simply the product of the weight of the
integration point we're interested in and the determinant of the Jacobian.

For the axisymmetric stress state, this same calculation is done, but due to the symmetry, the integration coefficient
is corrected by multiplying it with the axisymmetric circumference ($2\pi r$, with $r$ being the value of the radial
coordinate).
