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
by filling the B-matrix and calculate the integration coefficient (needs a correction depending on the symmetry of the
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

For the different stress states, the strain tensor, defined as:

| $\epsilon_{xx}$ $\epsilon_{xy}$ $\epsilon_{xz}$ |\
| $\epsilon_{yx}$ $\epsilon_{yy}$ $\epsilon_{yz}$ | \
| $\epsilon_{zx}$ $\epsilon_{zy}$ $\epsilon_{zz}$ |

This is converted to a strain vector for the different cases. For 3D, it is defined as:\
| $\epsilon_{xx}$ |\
| $\epsilon_{yy}$ |\
| $\epsilon_{zz}$ |\
| $\epsilon_{xy}$ |\
| $\epsilon_{yz}$ |\
| $\epsilon_{xz}$ |

For plane strain, it is defined as:\
| $\epsilon_{xx}$ |\
| $\epsilon_{yy}$ |\
| $0$ |\
| $\epsilon_{xy}$ |

For axisymmetric, it is defined as:\
| $\epsilon_{xx}$ |\
| $\epsilon_{yy}$ |\
| $u/r$ |\
| $\epsilon_{xy}$ |

Where $u$ is the displacement and $r$ is the radial coordinate. 


### B-matrix



### Integration coefficient

For the 3D and plane strain stress states, the integration coefficient is simply the product of the weight of the
integration point we're interested in and the determinant of the Jacobian.

For the axisymmetric stress state, this same calculation is done, but due to the symmetry, the integration coefficient
is corrected by multiplying it with the axisymmetric circumference ($2\pi r$, with $r$ being the value of the radial
coordinate).
