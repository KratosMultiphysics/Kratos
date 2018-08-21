
## Structural Mechanics Application

The Structural Mechanics Application contains a series of structural elements, as well as solid elements, constitutive laws and the corresponding strategies and solvers within Kratos Multiphysics.

<p align="center">
  <img src="https://github.com/KratosMultiphysics/Examples/raw/master/structural_mechanics/validation/beam_roll_up/data/rollup.gif" alt="Solution" style="width: 600px;"/>
</p>

The application includes tests to check the proper functioning of the application

### Features:

- A set of *Neumann* conditions:
     * Point loads (loads applied directly on the nodes)
     * Point moment (a discret moment applied directly on the nodes)
     * Line load (a distributed load applied over a line)
     * Surface load (a distributed load applied over a face)
     * A simple point contact conditions based on the distance

- Solid elements:
    * Small displacement elements
    * Total Lagrangian elements
    * Updated Lagrangian elements
    * Total Lagrangian prismatic solid-shell element (SPrism)

- Structural elements:
    * Zero-dimensional elements :
        * Nodal concentrated element (both 2D/3D). Includes nodal damping, nodal mass and nodal stiffness
    * Uni-dimensional elements :
        * Spring-damper element (3D)
        * Cable element (3D)
        * Truss element (3D)
        * Corrotational beam element (both 2D/3D)
    * Two-dimensional elements :
        * Membrane (pre-stressed)
        * Isotropic shell element
        * Thin shell (Quadrilateral and triangular)
        * Thick shell (Quadrilateral and triangular)

- Constitutive laws:
    * Isotropic laws (Plane strain, plane stress and 3D)
    * Orthotropic law (Plane stress)
    * Hyperelastic laws:
        * Neo-Hookean
        * Kirchhoff
    * Small displacement isotropic plasticity laws (just 3D):
        * Combining:
            * Yield surfaces:
                * VonMises
                * ModifiedMohrCoulomb
                * Tresca
                * DruckerPrager
            * Plastic potential:
                * VonMises
                * ModifiedMohrCoulomb
                * Tresca
                * DruckerPrager
        * Complete list:
            * SmallStrainIsotropicPlasticity3DVonMisesVonMises
            * SmallStrainIsotropicPlasticity3DVonMisesModifiedMohrCoulomb
            * SmallStrainIsotropicPlasticity3DVonMisesDruckerPrager
            * SmallStrainIsotropicPlasticity3DVonMisesTresca
            * SmallStrainIsotropicPlasticity3DModifiedMohrCoulombVonMises
            * SmallStrainIsotropicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb
            * SmallStrainIsotropicPlasticity3DModifiedMohrCoulombDruckerPrager
            * SmallStrainIsotropicPlasticity3DModifiedMohrCoulombTresca
            * SmallStrainIsotropicPlasticity3DTrescaVonMises
            * SmallStrainIsotropicPlasticity3DTrescaModifiedMohrCoulomb
            * SmallStrainIsotropicPlasticity3DTrescaDruckerPrager
            * SmallStrainIsotropicPlasticity3DTrescaTresca
            * SmallStrainIsotropicPlasticity3DDruckerPragerVonMises
            * SmallStrainIsotropicPlasticity3DDruckerPragerModifiedMohrCoulomb
            * SmallStrainIsotropicPlasticity3DDruckerPragerDruckerPrager
            * SmallStrainIsotropicPlasticity3DDruckerPragerTresca
    * Small displacement isotropic damage laws (just 3D):
        * Combining:
            * Yield surfaces:
                * VonMises
                * ModifiedMohrCoulomb
                * Tresca
                * DruckerPrager
                * Rankine
                * SimoJu
            * Damage potential:
                * VonMises
                * ModifiedMohrCoulomb
                * Tresca
                * DruckerPrager
        * Complete list:
            * SmallStrainIsotropicDamage3DVonMisesVonMises
            * SmallStrainIsotropicDamage3DVonMisesModifiedMohrCoulomb
            * SmallStrainIsotropicDamage3DVonMisesDruckerPrager
            * SmallStrainIsotropicDamage3DVonMisesTresca
            * SmallStrainIsotropicDamage3DModifiedMohrCoulombVonMises
            * SmallStrainIsotropicDamage3DModifiedMohrCoulombModifiedMohrCoulomb
            * SmallStrainIsotropicDamage3DModifiedMohrCoulombDruckerPrager
            * SmallStrainIsotropicDamage3DModifiedMohrCoulombTresca
            * SmallStrainIsotropicDamage3DTrescaVonMises
            * SmallStrainIsotropicDamage3DTrescaModifiedMohrCoulomb
            * SmallStrainIsotropicDamage3DTrescaDruckerPrager
            * SmallStrainIsotropicDamage3DTrescaTresca
            * SmallStrainIsotropicDamage3DDruckerPragerVonMises
            * SmallStrainIsotropicDamage3DDruckerPragerModifiedMohrCoulomb
            * SmallStrainIsotropicDamage3DDruckerPragerDruckerPrager
            * SmallStrainIsotropicDamage3DDruckerPragerTresca
            * SmallStrainIsotropicDamage3DRankineVonMises
            * SmallStrainIsotropicDamage3DRankineModifiedMohrCoulomb
            * SmallStrainIsotropicDamage3DRankineDruckerPrager
            * SmallStrainIsotropicDamage3DRankineTresca
            * SmallStrainIsotropicDamage3DSimoJuVonMises
            * SmallStrainIsotropicDamage3DSimoJuModifiedMohrCoulomb
            * SmallStrainIsotropicDamage3DSimoJuDruckerPrager
            * SmallStrainIsotropicDamage3DSimoJuTresca

- Adjoint Sensitivity Analysis:
    * This feature provides the framework to compute sensitivities of structural responses (e.g. displacements, strain energy or stresses) with respect to different types of design variables (e.g. nodal coordinates, material or cross-sectional properties or load intensity) with the adjoint approach

- Strategies:
    * Formfinding strategies
    * Eigensolver strategy
    * Harmonic analysis strategies
    * Arc-length strategy

- Schemes:
    * Relaxation scheme
    * Eigen solver scheme

- Convergence criteria:
    * For displacement and other DoF
    * For displacement and rotation

- Utilities and processe:
    * A process to post-process eigenvalues
    * A GiDIO utility for eigen values
    * Process to compute the global mass of the system
    * Process to identify the neighbours in a prismatic mesh
    * Process to transform a pure shell mesh (local dimension equal to 2), to solid-shell mesh (pure 3D mesh)
