### Compression Cap Hardening

The standard Mohr–Coulomb (MC) yield surface characterizes shear failure in geomaterials by relating the shear stress $\tau$ on a potential failure plane to the corresponding normal stress $\sigma$. While suitable for frictional materials, the standard Mohr-Coulomb envelope lacks a mechanism to confine the admissible stress space under high compressive pressure. As a result, the pure Mohr-Coulomb model cannot represent the compaction, crushing, and plastic volumetric hardening that occur in soils and rocks under high confining stresses.

To address this limitation, a compression cap is introduced. The cap provides a smooth or piecewise-smooth closure of the yield surface in the high-compression regime and allows the plasticity model to incorporate hardening or crushing effects. 

This document describes the combined MC–cap yield surface, its formulation, and its role in hardening.

### Mohr–Coulomb yield surface

In the $`(\sigma, \tau)`$ stress space, the Mohr-Coulmb yiel surface is expressed as:

```math
    F_{MC}(\sigma, \tau) = \tau + \sigma \sin⁡{\phi} - c \cos⁡{\phi} = 0
```
where:

- $`\sigma`$ = normal stress component
- $`\tau`$ = shear stress component
- $`c`$ = cohesion of material
- $`\phi`$ = Internal friction angle

In stress-invariant form, the MC yield function is typically written as:

```math
    F_{MC}(p, q) = q + p \sin⁡{\phi} - c \cos⁡{\phi} = 0
```
where:

- $p = \frac{1}{3} tr(\sigma)$ is the mean effective stress
- $q = \sqrt{\frac{3}{2}s:s}$ is the diviotoric norm, with $𝑠$ the deviatoric stress tensor.

This defines a hexagonal pyramid in principal stress space, but is shown as a straight line in the $`(\sigma, \tau)`$ Mohr plane.

### Compression Cap Concept
At high confining pressures, real geomaterials exhibit compaction and crushing rather than unlimited shear strength. The Mohr-Coulomb envelope alone allows unbounded compressive stresses. A cap yield surface introduces a limit to admissible volumetric compression and establishes a mechanism for volumetric plastic deformation and hardening.

In stress-invariant space, the cap is usually defined as an ellipse (or a smooth rounded surface) closing the Morh-Coulomb shear wall in the compressive regime.

### Cap yield surface
A common choice is an elliptical cap:

```math
    F_{cap}(p, q) = \left( \frac{q}{X} \right)^2 + p^2 - p_c^2
```
where:

- $p_c$ = cap position (controls hardening/softening),
- $X$ = cap size parameter

The cap intersects the MC surface at a smooth transition point to ensure the overall yield surface is convex. The cap expands or contracts depending on the accumulated plastic volumetric strain:

```math
    p_c = p_{c0} + H \epsilon^p
```
where:
- $p_{c0}$ = the initial cap position
- $H$ = the hardening modulus
- $\epsilon^p$ = the plastic volumetric strain

This produces compaction hardening: as the material densifies, it supports higher compressive stresses.

### Combined Mohr–Coulomb + Cap Yield Surface








