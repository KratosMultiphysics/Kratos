# Test K<sub>0</sub> procedure simple dike

**Author:** [Wijtze Pieter Kikstra](https://github.com/WPK4FEM)

## Case Specification
A test for the K<sub>0</sub> procedure process. In the first stage the underground vertical effective stresses in the column of soil follow from density, gravity constant, vertical distance to the surface and the presence of pore water pressure. For this a nonlinear ( Mohr Coulomb ) constitutive law is used. Horizontal effective stresses are derived from the vertical effective stresses by multiplication with K<sub>0,NC</sub> and checked at an integration point near the bottom of the column. In the second stage, a dam is built on top of the soil. Stresses far away from the dam location should not drastically change.
