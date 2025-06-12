# Test K<sub>0</sub> procedure normal consolidation with skew layers and a dam on top

**Author:** [Wijtze Pieter Kikstra](https://github.com/WPK4FEM)

## Case Specification
A test for the K<sub>0</sub> procedure process. Vertical effective stresses in the column of soil layers and a dam follow from soil layer density, gravity constant, vertical distance to the surface and the presence of pore water pressure. The soil layers are positioned at an angle deviating from horizontal. Horizontal effective stresses are derived from the vertical effective stresses by multiplication with K<sub>0</sub><sup>nc</sup> and checked at an integration point near the bottoms of the soil layers. As the layers deviate from horizontal and there is a dam weighting those down, the resulting horizontal stress field is not in equilibrium.
