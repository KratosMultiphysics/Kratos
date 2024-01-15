# Geomechanics Time Integration Schemes

## Introduction
The `Scheme` class in Kratos is used to create time integration schemes. It is an abstract class, for which Geomechanics has implemented a number of flavors, subdivided in the Backward Euler and the Generalized Newmark families.


## Code Structure
The bulk of the functionality is in the `GeoMechanicsTimeIntegrationScheme` class. It contains two lists of variables that are used in the time integration schemes. A list of first order scalar variables (such as water pressure or temperature) and a list of second order vector variables (such as displacements or rotations). For the first order time derivatives, only the first time derivative is taken into account, while for the second order time derivatives, both the first and the second time derivatives are considered.

## Backward Euler


## Generalized Newmark


## Dynamic and damped schemes
