# Reset of nodal positions after a non-converged step

This test reproduces issue [#15007](https://github.com/KratosMultiphysics/Kratos/issues/15007): when a solution step does not converge, the state is reverted to the start of that step. When the mesh is not moved (`"move_mesh_flag": false`), this revert must not change the nodal positions. The test is run via both the Python route and the C++ route.

## Case Specification

The model consists of two `UPwSmallStrainElement2D4N` elements of $`1 \mathrm{m}`$ by $`1 \mathrm{m}`$ that are placed next to each other:

- The left element (`strong_soil`) is linear elastic.
- The right element (`weak_soil`) uses a Mohr-Coulomb model with a tension cut-off. The tensile strength is $`20 \mathrm{Pa}`$, while the cohesion is high enough to prevent shear failure.

Both elements have a Young's modulus of $`100 \mathrm{Pa}`$ and a Poisson's ratio of $`0.0`$. The water pressure is fixed to $`0.0 \mathrm{Pa}`$ at all nodes. The left edge is fixed horizontally, and the bottom edge is fixed vertically. The model is stretched by prescribing the horizontal displacement of the right edge:

| Time [s] | Displacement [m] |
|----------|------------------|
| 0.0      | 0.0              |
| 1.0      | 0.2              |
| 1.5      | 0.398            |
| 2.0      | 1.0              |

The time step is $`1.0 \mathrm{s}`$ and `max_iterations` is 17:

1. The first step (to $`t = 1.0 \mathrm{s}`$) is elastic and results in a horizontal displacement of $`0.1 \mathrm{m}`$ of the middle nodes.
2. The second step (to $`t = 2.0 \mathrm{s}`$) makes the weak element yield. Since the material model uses the elastic stiffness matrix, it takes 21 iterations to converge, so the first attempt fails and the state is reverted to the start of the step.
3. The step is halved. The step to $`t = 1.5 \mathrm{s}`$ is still elastic, and the step to $`t = 2.0 \mathrm{s}`$ converges in 14 iterations, since the elastic element is already close to the tensile strength.

## Expected results

At $`t = 2.0 \mathrm{s}`$, both elements carry a tensile stress equal to the tensile strength of $`20 \mathrm{Pa}`$. Therefore, the strain of the elastic element is $`20 / 100 = 0.2`$, and the middle nodes are displaced by $`0.2 \mathrm{m}`$. The weak element takes the remaining elongation, which results in a strain of $`0.8`$.

The test asserts the output times (to make sure the step was cut back), the displacements, the strains and the stresses. For the Python route, it also asserts that the nodal coordinates still equal the initial ones.

If the nodal positions were updated during the revert, the elements would have a length of $`1.1 \mathrm{m}`$ for the remainder of the analysis, which changes the computed strains and displacements.
