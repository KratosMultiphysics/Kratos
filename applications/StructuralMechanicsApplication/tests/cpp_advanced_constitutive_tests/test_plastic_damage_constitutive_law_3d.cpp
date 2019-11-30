// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Sergio Jimenez
//

// System includes

// External includes

// Project includes
#include "includes/process_info.h"
#include "testing/testing.h"
#include "containers/model.h"

// Application includes
#include "structural_mechanics_application_variables.h"

// Integrator
#include "custom_advanced_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity.h"
#include "custom_advanced_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_damage.h"

// Yield surfaces
#include "custom_advanced_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_advanced_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/tresca_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"

// Constitutive law
#include "custom_advanced_constitutive/generic_small_strain_plastic_damage_model.h"
#include "includes/model_part.h"
#include "geometries/tetrahedra_3d_4.h"

namespace Kratos
{
namespace Testing
{
// We test the associated damage Constitutive laws...
typedef Node<3> NodeType;

/**
    * Check the correct calculation of the integrated stress with the CL's
    */

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawPlasticDamageInternalVariables, KratosStructuralMechanicsFastSuite)
{
    //
    // Test: check correct behavior of internal and calculated variables
    //

    typedef GenericSmallStrainPlasticDamageModel<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> PD_CL;

    Model current_model;
    ModelPart& r_test_model_part = current_model.CreateModelPart("Main");
    PD_CL cl = PD_CL();

    KRATOS_CHECK(cl.Has(DAMAGE));  // = True
    KRATOS_CHECK(cl.Has(UNIAXIAL_STRESS));  // = True
    KRATOS_CHECK(cl.Has(PLASTIC_DISSIPATION));  // = True
    KRATOS_CHECK(cl.Has(PLASTIC_STRAIN_VECTOR));  // = True
    KRATOS_CHECK(cl.Has(INTERNAL_VARIABLES));  // = True

    Vector internal_variables_w(9);
    internal_variables_w[0] = 0.0;
    internal_variables_w[1] = 0.1;
    internal_variables_w[2] = 0.2;
    internal_variables_w[3] = 0.3;
    internal_variables_w[4] = 0.4;
    internal_variables_w[5] = 0.5;
    internal_variables_w[6] = 0.6;
    internal_variables_w[7] = 0.7;
    internal_variables_w[8] = 0.8;
    cl.SetValue(INTERNAL_VARIABLES, internal_variables_w, r_test_model_part.GetProcessInfo());
    Vector internal_variables_r;  // CL should internally resize it to 6
    cl.GetValue(INTERNAL_VARIABLES, internal_variables_r);

    KRATOS_CHECK_NEAR(internal_variables_r.size(), 9., 1.e-5);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[0], 0.0, 1.e-5);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[1], 0.1, 1.e-5);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[2], 0.2, 1.e-5);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[3], 0.3, 1.e-5);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[4], 0.4, 1.e-5);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[5], 0.5, 1.e-5);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[6], 0.6, 1.e-5);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[7], 0.7, 1.e-5);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[8], 0.8, 1.e-5);  // = True
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawPlasticDamageExponential, KratosStructuralMechanicsFastSuite)
{
    typedef GenericSmallStrainPlasticDamageModel<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> PD_CL;

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
    strain_vector[0] = 0.01;
    strain_vector[1] = 0.0;
    strain_vector[2] = 0.0;
    strain_vector[3] = 0.0;
    strain_vector[4] = 0.0;
    strain_vector[5] = 0.0;

    material_properties.SetValue(DENSITY, 7850.0);
    material_properties.SetValue(YOUNG_MODULUS, 2.01e11);
    material_properties.SetValue(POISSON_RATIO, 0.3);
    material_properties.SetValue(FRACTURE_ENERGY, 50e5);
    material_properties.SetValue(FRACTURE_ENERGY_DAMAGE_PROCESS, 50e5);
    material_properties.SetValue(YIELD_STRESS_COMPRESSION, 8.389e8);
    material_properties.SetValue(YIELD_STRESS_TENSION, 8.389e8);
    material_properties.SetValue(HARDENING_CURVE, 0);
    material_properties.SetValue(SOFTENING_TYPE, 0);

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
    PD_CL PlasticDamage_CL = PD_CL();

    std::vector<double> expected_value;
    expected_value = {1.60206e+09,8.81364e+08,8.81364e+08,0,0,0};

    Vector TestStress;
    Vector dummy;
    PlasticDamage_CL.InitializeMaterial(material_properties, Geom, dummy);
    PlasticDamage_CL.CalculateMaterialResponseCauchy(cl_parameters);
    PlasticDamage_CL.FinalizeMaterialResponseCauchy(cl_parameters);
    TestStress = cl_parameters.GetStressVector();

    //Check the results
    for (int comp = 0; comp < 6; comp++) {
        KRATOS_CHECK_NEAR(expected_value[comp], TestStress[comp], 0.0001e+08);
    }
}

} // namespace Testing
} // namespace Kratos
