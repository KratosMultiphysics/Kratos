# Policies

## Introduction
In our element structure, the strategy-pattern (also known as policy-pattern), is used to gain control and configurability for the element functionality. This is done by encapsulating the behavior in an interface and using a pointer to the interface in the element. The behavior can then be changed by changing the pointer to a different implementation of the interface. This is a powerful tool for creating flexible and reusable code.

Currently, the following policies are implemented (still under development):
1. The [StressStatePolicy](#stressstatepolicy) - this is used to calculate several stress-related properties depending on the configuration. There are three flavours: the `AxisymmetricStressState`, the `PlaneStrainStressState` and the `ThreeDimensionalStressState`.

## StressStatePolicy
The stress state policy is used to easily configure elements to have a specific stress state. There are three different stress state implementations:
1. `AxisymmetricStressState`
2. `PlaneStrainStressState`  
3. `ThreeDimensionalStressState` 

They all derive from the same interface (`StressStatePolicy`) and can be used interchangeably to calculate the B-matrix, the integration coefficient and the Green Lagrange Strain vector. See the following schematic for the structure of the stress state policies: 

![stress_state_policies.svg](stress_state_policies.svg)

