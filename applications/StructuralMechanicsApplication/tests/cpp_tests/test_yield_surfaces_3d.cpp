// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "includes/process_info.h"

// Application includes
#include "custom_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/yield_surfaces/tresca_yield_surface.h"

namespace Kratos
{
namespace Testing
{

void GenerateTestVariables(
    Vector& rStressVector,
    Vector& rStrainVector,
    Properties& rMaterialProperties
)
{
    rStressVector = ZeroVector(6);
    rStressVector[0] = 1.0e6;
    rStressVector[1] = 1.0e6;
    rStressVector[2] = 2.0e6;
    rStressVector[3] = 0.5e6;
    rStressVector[4] = 0.5e6;
    rStressVector[5] = 0.0;

    rStrainVector = ZeroVector(6);
    rStrainVector[0] = 0.001;
    rStrainVector[1] = 0.001;
    rStrainVector[2] = 0.002;
    rStrainVector[3] = 0.00001;
    rStrainVector[4] = 0.0;
    rStrainVector[5] = 0.00001;

    rMaterialProperties.SetValue(YOUNG_MODULUS, 210e9);
    rMaterialProperties.SetValue(POISSON_RATIO, 0.22);
    rMaterialProperties.SetValue(YIELD_STRESS_COMPRESSION, 30.0e6);
    rMaterialProperties.SetValue(YIELD_STRESS_TENSION, 3.0e6);
    rMaterialProperties.SetValue(FRICTION_ANGLE, 32.0);
    rMaterialProperties.SetValue(DILATANCY_ANGLE, 16.0);
}


/** 
* Check the correct calculation of the uniaxial stress of the yield surfaces
*/
KRATOS_TEST_CASE_IN_SUITE(YieldSurfacesUniaxialStress, KratosStructuralMechanicsFastSuite)
{
    Vector Strain, Stress;
    Properties TestProperties;
    GenerateTestVariables(Stress, Strain, TestProperties);

    typedef ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential> MC;
    typedef VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential> VM;
    typedef DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential> DP;
    typedef RankineYieldSurface<ModifiedMohrCoulombPlasticPotential> R;
    typedef TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential> T;
    typedef SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential> SJ;

    // Analytical solutions of the yield surfaces
    double MCres, VMres, DPres, Rres, Tres, SJres;
    MCres = 2.1991e+07;
    VMres = 1.58114e+06;
    DPres = 1.17878e+07;
    Rres = 2.2406e+06;
    Tres = 1.82564e+06;
    SJres =  774.919;

    double TestMC, TestVM, TestDP, TestR, TestT, TestSJ;
    MC::CalculateEquivalentStress(rStressVector, rStrainVector, TestMC, rMaterialProperties);
    VM::CalculateEquivalentStress(rStressVector, rStrainVector, TestVM, rMaterialProperties);
    DP::CalculateEquivalentStress(rStressVector, rStrainVector, TestDP, rMaterialProperties);
    R::CalculateEquivalentStress(rStressVector, rStrainVector, TestR, rMaterialProperties);
    T::CalculateEquivalentStress(rStressVector, rStrainVector, TestT, rMaterialProperties);
    SJ::CalculateEquivalentStress(rStressVector, rStrainVector, TestSJ, rMaterialProperties);

    KRATOS_CHECK_NEAR(MCres, TestMC, 100.0);
    KRATOS_CHECK_NEAR(VMres, TestVM, 100.0);
    KRATOS_CHECK_NEAR(DPres, TestDP, 100.0);
    KRATOS_CHECK_NEAR(Rres, TestR, 100.0);
    KRATOS_CHECK_NEAR(Tres, TestT, 100.0);
    KRATOS_CHECK_NEAR(SJres, TestSJ, 1.0);
}















}
}