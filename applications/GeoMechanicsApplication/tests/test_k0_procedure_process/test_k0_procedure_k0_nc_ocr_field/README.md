# Test K<sub>0</sub> procedure normal consolidation with ocr field input.

**Author:** [Wijtze Pieter Kikstra](https://github.com/WPK4FEM)

## Case Specification
A test for the K<sub>0</sub> procedure process. Vertical effective stresses in the column of soil follow from density, gravity constant, vertical distance to the surface and the presence of pore water pressure. Horizontal effective stresses are derived from the vertical effective stresses by multiplication with K<sub>0</sub><sup>nc</sup> and OCR and checked at an integration point near the bottom of the column. The OCR values are supplied through a set_parameter_field_process.
