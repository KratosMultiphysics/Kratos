# Test master-slave constraints

**Author:** [Anne van de Graaf](https://github.com/avdg81)

## Case specification

This test consists of two linear quadrilateral U-Pw elements that are connected using master-slave constraints.  Both elements have dimensions of 2.0 m by 2.0 m.  The model has been supported by a fixed node (at the bottom left corner) and a vertical slider (along the left edge).  The model is loaded by a uniform normal line load along the right edge, which leads to elongation of the elements.  Note that the two elements don't have any shared nodes.  They are tied together by applying master-slave constraints, where the master side is the right edge of the leftmost element, and the slave side is the left edge of the rightmost element.  The following figure shows the adopted mesh.

<img src="mesh.svg" width="600">

To keep the model as simple as possible, a linear elastic material model has been applied.  Also, the groundwater pressure field has been fixed at a value of 0.0 throughout the entire domain.

The goal of this test is to show that master-slave constraints are correctly applied.  We verify that by comparing the actual displacements in $`x`$ direction along the loaded side with the analytical solution:
```math
u_{x} = \epsilon_{xx} \cdot L = \frac{\sigma_{xx}}{E} \cdot L
```
where $`u_{x}`$ equals the displacements in $`x`$ direction at the loaded side, $`\epsilon_{xx}`$ is the normal strain in $`x`$ direction, $`L`$ is the length of the model in $`x`$ direction, $`\sigma_{xx}`$ is the normal stress in $`x`$ direction (which equals the applied uniform normal edge load), and $`E`$ is the Young's modulus of the soil.  Note that the above formula implies that:
- We are using a linear strain measure.
- We assume a uniform deformation and stress state over the depth of the elements.  This allows us to reduce the problem formulation to an extensional bar.

Furthermore, we compare the reaction forces in $`x`$ direction along the supported side with the analytical solution:
```math
R_{x} = \int_{A} \sigma_{xx} dA
```
where $`R_{x}`$ is the reaction force component in $`x`$ direction, $`\sigma_{xx}`$ is the normal stress in $`x`$ direction (which equals the applied uniform normal edge load), and $`A`$ the area associated with a node for the integration of $`\sigma_{xx}`$.  With the above assumptions in place, we can simplify the above formula to: $`R_{x} = \sigma_{xx} \cdot \frac{h}{2}`$, where $`h`$ is the depth of an element (the thickness equals unity, for this is a plane strain problem).

Finally, we also check whether the "tied" nodes have equal displacements.