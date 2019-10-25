// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Sergio Jimenez/Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"

// Application includes
// Integrator
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_damage.h"
// Yield surfaces
#include "custom_constitutive/yield_surfaces/von_mises_yield_surface.h"
// Plastic potentials
#include "custom_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
// Constitutive law
#include "custom_constitutive/generic_small_strain_high_cycle_fatigue_law.h"
#include "geometries/tetrahedra_3d_4.h"

namespace Kratos
{
namespace Testing
{
// We test the associated damage Constitutive laws...
typedef Node<3> NodeType;
typedef GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> FatigueLawType;

/**
    * Check the correct calculation of the integrated stress with the CL's
    */
KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawHighCycleFatigueExponential, KratosStructuralMechanicsFastSuite)
{
    ConstitutiveLaw::Parameters cl_parameters;
    Properties material_properties;
    Vector stress_vector, strain_vector;

    Model current_model;

    ModelPart& test_model_part = current_model.CreateModelPart("Main");

    NodeType::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    NodeType::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    NodeType::Pointer p_node_3 = test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    NodeType::Pointer p_node_4 = test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);

    Tetrahedra3D4<NodeType> Geom = Tetrahedra3D4<NodeType>(p_node_1, p_node_2, p_node_3, p_node_4);

    stress_vector = ZeroVector(6);

    strain_vector = ZeroVector(6);
    strain_vector[0] = 5.0e-03;
    strain_vector[1] = 0.0;
    strain_vector[2] = 0.0;
    strain_vector[3] = 0.0;
    strain_vector[4] = 0.0;
    strain_vector[5] = 0.0;

    material_properties.SetValue(DENSITY, 7850.0);
    material_properties.SetValue(YOUNG_MODULUS, 2.01e11);
    material_properties.SetValue(POISSON_RATIO, 0.3);
    material_properties.SetValue(FRACTURE_ENERGY, 20e5);
    material_properties.SetValue(YIELD_STRESS, 838.9e6);
    material_properties.SetValue(SOFTENING_TYPE, 1);
    Vector high_cycle_fatigue_coefficients(7);
    high_cycle_fatigue_coefficients[0]=0.5;
    high_cycle_fatigue_coefficients[1]=0.7;
    high_cycle_fatigue_coefficients[2]=0.5;
    high_cycle_fatigue_coefficients[3]=0.0068;
    high_cycle_fatigue_coefficients[4]=3.35;
    high_cycle_fatigue_coefficients[5]=0.0133;
    high_cycle_fatigue_coefficients[6]=0.0068;
    material_properties.SetValue(HIGH_CYCLE_FATIGUE_COEFFICIENTS, high_cycle_fatigue_coefficients);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions = cl_parameters.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    cl_parameters.SetElementGeometry(Geom);
    cl_parameters.SetProcessInfo(test_model_part.GetProcessInfo());
    cl_parameters.SetMaterialProperties(material_properties);
    cl_parameters.SetStrainVector(strain_vector);
    cl_parameters.SetStressVector(stress_vector);
    Matrix const_matrix;
    cl_parameters.SetConstitutiveMatrix(const_matrix);

    // Create the CL's
    FatigueLawType HCF_CL = FatigueLawType(8.8914e-01, 7.6380e+08, 6.9760e+08, 7.6380e+08, -2.2914e+08, 984, 5.3239E-07);

    std::vector<double> expected_value;
    expected_value = {1.07636e+09,4.61299e+08,4.61299e+08,0,0,0};

    Vector TestStress;
    Vector dummy;
    HCF_CL.InitializeMaterial(material_properties, Geom, dummy);
    HCF_CL.CalculateMaterialResponseCauchy(cl_parameters);
    HCF_CL.FinalizeMaterialResponseCauchy(cl_parameters);
    TestStress = cl_parameters.GetStressVector();

    //Check the results
    KRATOS_CHECK_VECTOR_NEAR(expected_value, TestStress, 0.0001e+08);
}

} // namespace Testing
} // namespace Kratos
