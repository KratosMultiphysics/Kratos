
# Constitutive Laws Application

 |             **Application**             |                                                                                    **Description**                                                                                    |                              **Status**                              | **Authors** |
|:---------------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:--------------------------------------------------------------------:|:-----------:|
| `ConstitutiveLawsApplication` | The *Constitutive Laws Application* contains a series of constitutive laws implementations within *Kratos Multiphysics*. | <img src="https://img.shields.io/badge/Status-%F0%9F%9A%80%20Actively%20developed-Green"  width="300px"> | Alejandro Cornejo Velázquez *(acornejo@cimne.upc.edu )* <br />  Sergio Jimenez Reyes *(sjimenez@cimne.upc.edu)* <br /> Riccardo Rossi *(rrossi@cimne.upc.edu)* <br /> Rubén Zorrilla Martínez *(rzorrilla@cimne.upc.edu)* <br /> Vicente Mataix Ferrándiz *(vmataix@altair.com)*    |


The application includes tests to check the proper functioning of the application.

## 😎 Features:

- **Constitutive laws**
    * *Orthotropic law (Plane stress)*
    * *Hyperelastic laws*
        * Neo-Hookean
        * Kirchhoff
    * *Small displacement isotropic plasticity laws (just 3D)*
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
    * *Small displacement isotropic damage laws (just 3D)*
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

- **Utilities**
    * *Generic constitutive laws utilities*
    * *Tangent operator AD* 

- **Processes**
    * *Automatic initial damage*
    * *Advance in time HCF*

- **Several python unittest, including Validation tests, and several cpp tests**

## ⚙️ Examples:

Examples can be found [in the same folder as the *Structural Mechanics Application*](https://github.com/KratosMultiphysics/Examples/tree/master/structural_mechanics).