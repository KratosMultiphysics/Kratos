# Single element tests with Mohr-Coulomb material model

**Author:** [Anne van de Graaf](https://github.com/avdg81)

## Case Specification

Each test consists of a single quadrilateral element with sides of $`1 \mathrm{m}`$ by $`1 \mathrm{m}`$.  The tested element types include:

1. `UPwSmallStrainElement2D4N`
2. `SmallStrainUPwDiffOrderElement2D8N`

In both tests, all degrees of freedom have been prescribed.  In vertical direction, the elements are compressed by imposing a vertical displacement of the element's top edge by $`-0.015 \mathrm{m}`$.  In horizontal direction, the elements are extended by displacing the right edge by $`+0.015 \mathrm{m}`$.  The water pressure degrees of freedom have all been set to $`0.0 \mathrm{Pa}`$.

The elastic behavior of the elements is described by a Young's modulus of $`100 \mathrm{Pa}`$ and a Poisson's ratio of $`0.0`$.  The material behavior is described using a Mohr-Coulomb model with a cohesion $`c`$ of $`2.0 \mathrm{Pa}`$ and a friction angle $`\phi = 0.0^{\circ}`$.  The material model has been implemented by the shared library `MohrCoulomb64.dll`.

The prescribed displacement field results in the following stress field: $`\sigma_{xx} = +1.5 \mathrm{Pa}`$, $`\sigma_{yy} = -1.5 \mathrm{Pa}`$, and $`\sigma_{zz} = 0.0 \mathrm{Pa}`$.  Note that only elastic deformations occur. 

Both tests assert the prescribed displacement fields (to ensure that the correct Dirichlet boundary conditions have been applied).  In addition, these tests assert the shear capacity, which equals $`0.75`$ in each integration point.
