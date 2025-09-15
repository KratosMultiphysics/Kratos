# Test K<sub>0</sub> procedure with umat

**Author:** [Wijtze Pieter Kikstra](https://github.com/WPK4FEM)

## Case Specification
A test for the K<sub>0</sub> procedure process. Vertical effective stresses in the column of soil follow from density, gravity constant, vertical distance to the surface and the presence of pore water pressure. Horizontal effective stresses are derived from the vertical effective stresses by multiplication with K<sub>0</sub> = 1 - sin(&Phi;) and checked at an integration point near the bottom of the column. Phi (&Phi;) is retrieved from the parameter list of the user supplied material.
