# Geomechanics Time Integration Schemes

## Introduction
The `Scheme` class in Kratos is used to create time integration schemes. It is an abstract class, for which Geomechanics has implemented a number of flavors, subdivided in the Backward Euler and the Generalized Newmark families.

![img.png](img.png)

The bulk of the functionality is in the `GeoMechanicsTimeIntegrationScheme` class. It contains two lists of variables that are used in the time integration schemes. 

