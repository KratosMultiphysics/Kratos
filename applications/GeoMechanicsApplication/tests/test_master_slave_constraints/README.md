# Test master-slave constraints

**Author:** [Anne van de Graaf](https://github.com/avdg81)

## Case specification

This test consists of two linear quadrilateral U-Pw elements that are connected using master-slave constraints.  Both elements have dimensions of 2.0 m by 2.0 m.  The model has been supported by a fixed node (at the bottom left corner) and a vertical slider (along the left edge).  The model is loaded by a uniform normal line load along the right edge, which leads to elongation of the elements.  Note that the two elements don't have any shared nodes.  They are tied together by applying master-slave constraints, where the master side is the right edge of the leftmost element, and the slave side is the left edge of the rightmost element.

To keep the model as simple as possible, a linear elastic material model has been applied.  Also, the groundwater pressure field has been fixed at a value of 0.0 throughout the entire domain.

The goal of this test is to show that master-slave constraints are correctly applied.  We verify that by comparing the actual displacements in $`x`$ direction along the loaded side with the analytical solution.  Furthermore, we compare the reaction forces in $`x`$ direction along the supported side with the analytical solution.