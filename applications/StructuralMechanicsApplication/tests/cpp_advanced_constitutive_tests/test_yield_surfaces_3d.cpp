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
#include "testing/testing.h"

// Application includes

// Yields
#include "custom_advanced_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/tresca_yield_surface.h"
// Plastic Potentials
#include "custom_advanced_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"

namespace Kratos
{
namespace Testing
{
typedef Node<3> NodeType;
typedef ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>> MC;
typedef VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>> VM;
typedef DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>> DP;
typedef RankineYieldSurface<ModifiedMohrCoulombPlasticPotential<6>> R;
typedef TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>> T;
typedef SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential<6>> SJ;

void GenerateTestVariables(
    array_1d<double, 6> &rStressVector,
    Vector &rStrainVector,
    Properties &rMaterialProperties)
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
    rMaterialProperties.SetValue(SOFTENING_TYPE, 0);
    rMaterialProperties.SetValue(FRACTURE_ENERGY, 1.0e3);
}

/**
* Check the correct calculation of the uniaxial stress of the yield surfaces
*/
KRATOS_TEST_CASE_IN_SUITE(YieldSurfacesUniaxialStress, KratosStructuralMechanicsFastSuite)
{
    Vector Strain;
    array_1d<double, 6> Stress;
    Properties material_properties;
    GenerateTestVariables(Stress, Strain, material_properties);

    // Analytical solutions of the yield surfaces
    double MCres, VMres, DPres, Rres, Tres, SJres;
    MCres = 2.1991e+07;
    VMres = 1.58114e+06;
    DPres = 5.77553e+06;
    Rres = 2.2406e+06;
    Tres = 1.82564e+06;
    SJres = 774.919;

    ProcessInfo dummy_process_info;
    Geometry<NodeType> dummy_geometry;
    ConstitutiveLaw::Parameters aux_param(dummy_geometry, material_properties, dummy_process_info);

    // Solutions to test...
    double TestMC = 0.0, TestVM = 0.0, TestDP = 0.0, TestR = 0.0, TestT = 0.0, TestSJ = 0.0;
    MC::CalculateEquivalentStress(Stress, Strain, TestMC, aux_param);
    VM::CalculateEquivalentStress(Stress, Strain, TestVM, aux_param);
    DP::CalculateEquivalentStress(Stress, Strain, TestDP, aux_param);
    R::CalculateEquivalentStress(Stress, Strain, TestR, aux_param);
    T::CalculateEquivalentStress(Stress, Strain, TestT, aux_param);
    SJ::CalculateEquivalentStress(Stress, Strain, TestSJ, aux_param);

    // Check the results!
    KRATOS_CHECK_NEAR(MCres, TestMC, 0.001e6);
    KRATOS_CHECK_NEAR(VMres, TestVM, 0.0001e6);
    KRATOS_CHECK_NEAR(DPres, TestDP, 0.001e6);
    KRATOS_CHECK_NEAR(Rres, TestR, 0.0001e6);
    KRATOS_CHECK_NEAR(Tres, TestT, 0.0001e6);
    KRATOS_CHECK_NEAR(SJres, TestSJ, 0.01);
}

/**
    * Check the correct calculation of the derivatives of the yield surfaces
    */
KRATOS_TEST_CASE_IN_SUITE(YieldSurfacesDerivatives, KratosStructuralMechanicsFastSuite)
{
    double I1, J2;
    Vector Strain;
    array_1d<double, 6> Stress;
    array_1d<double, 6> Deviator = ZeroVector(6);
    Properties material_properties;

    GenerateTestVariables(Stress, Strain, material_properties);
    ConstitutiveLawUtilities<6>::CalculateI1Invariant(Stress, I1);
    ConstitutiveLawUtilities<6>::CalculateJ2Invariant(Stress, I1, Deviator, J2);

    // Analytical solutions of the yield surfaces derivatives
    std::vector<double> MCres, VMres, DPres, Rres, Tres, SJres;
    MCres = {0.109261, 2.07822, 10.6714, 2.6863, 11.8748, 2.62528};
    VMres = {-0.316228, -0.316228, 0.632456, 0.948683, 0.948683, 0.0};
    DPres = {0.197647,0.197647,1.85929,1.66165,1.66165,0};
    Tres = {-0.369513, -0.364032, 0.733545, 1.08113, 1.10671, 0.0073077};

    array_1d<double, 6> TestMC = ZeroVector(6), TestVM = ZeroVector(6), TestDP = ZeroVector(6), TestT = ZeroVector(6);

    ProcessInfo dummy_process_info;
    Geometry<NodeType> dummy_geometry;
    ConstitutiveLaw::Parameters aux_param(dummy_geometry, material_properties, dummy_process_info);

    MC::CalculateYieldSurfaceDerivative(Stress, Deviator, J2, TestMC, aux_param);
    VM::CalculateYieldSurfaceDerivative(Stress, Deviator, J2, TestVM, aux_param);
    DP::CalculateYieldSurfaceDerivative(Stress, Deviator, J2, TestDP, aux_param);
    T::CalculateYieldSurfaceDerivative(Stress, Deviator, J2, TestT, aux_param);

    // Check the results!
    for (int comp = 0; comp < 6; comp++) {
        KRATOS_CHECK_NEAR(MCres[comp], TestMC[comp], 1.0e-3);
        KRATOS_CHECK_NEAR(VMres[comp], TestVM[comp], 1.0e-3);
        KRATOS_CHECK_NEAR(DPres[comp], TestDP[comp], 1.0e-3);
        KRATOS_CHECK_NEAR(Tres[comp], TestT[comp], 1.0e-3);
    }
}

/**
    * Check the correct calculation of the initial threshold of the yield surfaces
    */
KRATOS_TEST_CASE_IN_SUITE(YieldSurfacesInitialUniaxialThreshold, KratosStructuralMechanicsFastSuite)
{
    Properties material_properties;
    Vector Strain;
    array_1d<double, 6> Stress;
    GenerateTestVariables(Stress, Strain, material_properties);

    // Analytical solutions of the initial threshold
    double MCres, VMres, DPres, Rres, Tres, SJres;
    MCres = 30.0e6;
    VMres = 3.0e6;
    DPres = 7.50918e+06;
    Rres = 3.0e6;
    Tres = 3.0e6;
    SJres = 65.4654;

    ProcessInfo dummy_process_info;
    Geometry<NodeType> dummy_geometry;
    ConstitutiveLaw::Parameters aux_param(dummy_geometry, material_properties, dummy_process_info);

    double TestMC = 0.0, TestVM = 0.0, TestDP = 0.0, TestR = 0.0, TestT, TestSJ = 0.0;
    MC::GetInitialUniaxialThreshold(aux_param, TestMC);
    VM::GetInitialUniaxialThreshold(aux_param, TestVM);
    DP::GetInitialUniaxialThreshold(aux_param, TestDP);
    R::GetInitialUniaxialThreshold(aux_param, TestR);
    T::GetInitialUniaxialThreshold(aux_param, TestT);
    SJ::GetInitialUniaxialThreshold(aux_param, TestSJ);

    // Check the results!
    KRATOS_CHECK_NEAR(MCres, TestMC, 0.0001e6);
    KRATOS_CHECK_NEAR(VMres, TestVM, 0.0001e6);
    KRATOS_CHECK_NEAR(DPres, TestDP, 0.0001e6);
    KRATOS_CHECK_NEAR(Rres, TestR, 0.0001e6);
    KRATOS_CHECK_NEAR(Tres, TestT, 0.0001e6);
    KRATOS_CHECK_NEAR(SJres, TestSJ, 1.0e-3);
}

/**
    * Check the correct calculation of the Damage Parameter of the yield surfaces
    */
KRATOS_TEST_CASE_IN_SUITE(YieldSurfacesIDamageParameterLinear, KratosStructuralMechanicsFastSuite)
{
    Properties material_properties;
    Vector Strain;
    array_1d<double, 6> Stress;
    GenerateTestVariables(Stress, Strain, material_properties);
    double characteristic_length = 0.1;

    // Analytical solutions damage parameter Linear Softening
    double MCres, VMres, DPres, Rres, Tres, SJres;
    MCres = -0.00214286;
    VMres = -0.214286;
    DPres = -0.00214286;
    Rres = -0.214286;
    Tres = -0.00214286;
    SJres = -4.500000e+08;

    ProcessInfo dummy_process_info;
    Geometry<NodeType> dummy_geometry;
    ConstitutiveLaw::Parameters aux_param(dummy_geometry, material_properties, dummy_process_info);

    double TestMC = 0.0, TestVM = 0.0, TestDP = 0.0, TestR = 0.0, TestT = 0.0, TestSJ = 0.0;
    MC::CalculateDamageParameter(aux_param, TestMC, characteristic_length);
    VM::CalculateDamageParameter(aux_param, TestVM, characteristic_length);
    DP::CalculateDamageParameter(aux_param, TestDP, characteristic_length);
    R::CalculateDamageParameter(aux_param, TestR, characteristic_length);
    T::CalculateDamageParameter(aux_param, TestT, characteristic_length);
    SJ::CalculateDamageParameter(aux_param, TestSJ, characteristic_length);

    // Check the results!
    KRATOS_CHECK_NEAR(MCres, TestMC, 0.001);
    KRATOS_CHECK_NEAR(VMres, TestVM, 0.1);
    KRATOS_CHECK_NEAR(DPres, TestDP, 0.001);
    KRATOS_CHECK_NEAR(Rres, TestR, 0.001);
    KRATOS_CHECK_NEAR(Tres, TestT, 0.001);
    KRATOS_CHECK_NEAR(SJres, TestSJ, 0.1e8);
}

/**
	* Check the correct calculation of the Damage Parameter of the yield surfaces
	*/
KRATOS_TEST_CASE_IN_SUITE(YieldSurfacesIDamageParameterExponential, KratosStructuralMechanicsFastSuite)
{
    Properties material_properties;
    Vector Strain;
    array_1d<double, 6> Stress;
    GenerateTestVariables(Stress, Strain, material_properties);
    double characteristic_length = 0.1;

    material_properties.SetValue(SOFTENING_TYPE, 1);

    // Analytical solutions damage parameter Exponential Softening
    double MCres, VMres, DPres, Rres, Tres;//, SJres;
    MCres = 0.00429492;
    VMres = 0.545455;
    DPres = 0.00429492;
    Rres = 0.545455;
    Tres = 0.00429492;
//     SJres = -2.00000000;

    ProcessInfo dummy_process_info;
    Geometry<NodeType> dummy_geometry;
    ConstitutiveLaw::Parameters aux_param(dummy_geometry, material_properties, dummy_process_info);

    double TestMC = 0.0, TestVM = 0.0, TestDP = 0.0, TestR = 0.0, TestT = 0.0;//, TestSJ = 0.0;
    MC::CalculateDamageParameter(aux_param, TestMC, characteristic_length);
    VM::CalculateDamageParameter(aux_param, TestVM, characteristic_length);
    DP::CalculateDamageParameter(aux_param, TestDP, characteristic_length);
    R::CalculateDamageParameter(aux_param, TestR, characteristic_length);
    T::CalculateDamageParameter(aux_param, TestT, characteristic_length);

    // Check the results!
    KRATOS_CHECK_NEAR(MCres, TestMC, 0.0001);
    KRATOS_CHECK_NEAR(VMres, TestVM, 0.0001);
    KRATOS_CHECK_NEAR(DPres, TestDP, 0.0001);
    KRATOS_CHECK_NEAR(Rres, TestR, 0.0001);
    KRATOS_CHECK_NEAR(Tres, TestT, 0.0001);
}
} // namespace Testing
} // namespace Kratos
