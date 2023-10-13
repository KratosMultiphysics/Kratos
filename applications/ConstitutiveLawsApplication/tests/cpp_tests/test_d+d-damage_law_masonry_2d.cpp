// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Philip Kalkbrenner
//					 Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "includes/process_info.h"
#include "testing/testing.h"
#include "containers/model.h"

// Application includes
#include "constitutive_laws_application_variables.h"

// Constitutive law
#include "custom_constitutive/small_strains/damage/plane_stress_d_plus_d_minus_damage_masonry_2d.h"
#include "includes/model_part.h"
#include "geometries/triangle_2d_3.h"
#include "constitutive_laws_application_variables.h"

namespace Kratos
{
namespace Testing
{
// We test the d+ d- damage Constitutive Law (for masonry) 2D...
typedef Node NodeType;

KRATOS_TEST_CASE_IN_SUITE(DPlusDMinusMasonry2DInternalVariables, KratosConstitutiveLawsFastSuite)
{
    //
    // Test: check correct behavior of internal and calculated variables
    //

    Model current_model;
    ModelPart &test_model_part = current_model.CreateModelPart("Main");
    DamageDPlusDMinusMasonry2DLaw cl = DamageDPlusDMinusMasonry2DLaw();

    KRATOS_EXPECT_FALSE(cl.Has(INTEGRATED_STRESS_TENSOR));  // = False, in order to use CalculateValue())
    KRATOS_EXPECT_TRUE(cl.Has(DAMAGE_TENSION));  // = True
    KRATOS_EXPECT_TRUE(cl.Has(THRESHOLD_TENSION));  // = True
    KRATOS_EXPECT_TRUE(cl.Has(DAMAGE_COMPRESSION));  // = True
    KRATOS_EXPECT_TRUE(cl.Has(THRESHOLD_COMPRESSION));  // = True
    KRATOS_EXPECT_TRUE(cl.Has(UNIAXIAL_STRESS_COMPRESSION));  // = True
    KRATOS_EXPECT_TRUE(cl.Has(UNIAXIAL_STRESS_TENSION));  // = True
    KRATOS_EXPECT_TRUE(cl.Has(INTERNAL_VARIABLES));  // = True

    Vector internal_variables_w(6);
    internal_variables_w[0] = 0.0;
    internal_variables_w[1] = 0.1;
    internal_variables_w[2] = 0.2;
    internal_variables_w[3] = 0.3;
    internal_variables_w[4] = 0.4;
    internal_variables_w[5] = 0.5;
    cl.SetValue(INTERNAL_VARIABLES, internal_variables_w, test_model_part.GetProcessInfo());
}

KRATOS_TEST_CASE_IN_SUITE(DPlusDMinusMasonry2DPureTensionTest, KratosConstitutiveLawsFastSuite)
{
    ConstitutiveLaw::Parameters cl_parameters;
    Properties material_properties;
    ProcessInfo process_info;
    Vector stress_vector, strain_vector;

    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("Main");

    NodeType::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    NodeType::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.5, 0.0);
    NodeType::Pointer p_node_3 = test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);

    Triangle2D3<NodeType> Geom = Triangle2D3<NodeType>(p_node_1, p_node_2, p_node_3);

    stress_vector = ZeroVector(3);
    strain_vector = ZeroVector(3);
    strain_vector[0] = 0.0005;
    strain_vector[1] = 0.0;
    strain_vector[2] = 0.0;

    material_properties.SetValue(YOUNG_MODULUS, 3718.0e+6);
    material_properties.SetValue(POISSON_RATIO, 0.2);
    material_properties.SetValue(YIELD_STRESS_TENSION, 1.559e+6);
    material_properties.SetValue(FRACTURE_ENERGY_TENSION, 10.0e+2);
    material_properties.SetValue(DAMAGE_ONSET_STRESS_COMPRESSION, 10.0e+6);
    material_properties.SetValue(YIELD_STRESS_COMPRESSION, 17.99e+6);
    material_properties.SetValue(YIELD_STRAIN_COMPRESSION, 0.02);
    material_properties.SetValue(RESIDUAL_STRESS_COMPRESSION, 2.00e+6);
    material_properties.SetValue(BIAXIAL_COMPRESSION_MULTIPLIER, 1.20);
    material_properties.SetValue(FRACTURE_ENERGY_COMPRESSION, 8.0e+5);
    material_properties.SetValue(SHEAR_COMPRESSION_REDUCTOR, 0.16);
    material_properties.SetValue(BEZIER_CONTROLLER_C1, 0.65);
    material_properties.SetValue(BEZIER_CONTROLLER_C2,  0.45);
    material_properties.SetValue(BEZIER_CONTROLLER_C3, 1.5);
    material_properties.SetValue(INTEGRATION_IMPLEX, 0);
    material_properties.SetValue(TENSION_YIELD_MODEL, 0);
    material_properties.SetValue(DELAY_TIME, 1.0);

    Flags cl_options;
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    cl_parameters.SetElementGeometry(Geom);
    cl_parameters.SetMaterialProperties(material_properties);
    cl_parameters.SetStrainVector(strain_vector);
    cl_parameters.SetStressVector(stress_vector);
    cl_parameters.SetProcessInfo(process_info);
    cl_parameters.SetOptions(cl_options);
    Matrix const_matrix;
    cl_parameters.SetConstitutiveMatrix(const_matrix);

    // Create the CL
    DamageDPlusDMinusMasonry2DLaw masonry2d_cl = DamageDPlusDMinusMasonry2DLaw();

    //Initialize Material
    const Vector& r_dummy = ZeroVector(3);
    masonry2d_cl.InitializeMaterial(material_properties, Geom, r_dummy);

    std::vector<double> masonry2d_res;
    masonry2d_res = {1327960, 265592, 0};

    Vector test_masonry2d_stress;
    masonry2d_cl.CalculateMaterialResponseCauchy(cl_parameters);

    test_masonry2d_stress = cl_parameters.GetStressVector();

    // Check the results
    KRATOS_EXPECT_VECTOR_NEAR(test_masonry2d_stress, masonry2d_res, 0.01e6);
}

KRATOS_TEST_CASE_IN_SUITE(DPlusDMinusMasonry2DPureCompressionTest, KratosConstitutiveLawsFastSuite)
{
    ConstitutiveLaw::Parameters cl_parameters;
    Properties material_properties;
    ProcessInfo process_info;
    Vector stress_vector, strain_vector;

    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("Main");

    NodeType::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    NodeType::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.5, 0.0);
    NodeType::Pointer p_node_3 = test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);

    Triangle2D3<NodeType> Geom = Triangle2D3<NodeType>(p_node_1, p_node_2, p_node_3);

    stress_vector = ZeroVector(3);
    strain_vector = ZeroVector(3);
    strain_vector[0] = -0.0093;
    strain_vector[1] = 0.0;
    strain_vector[2] = 0.0;

    material_properties.SetValue(YOUNG_MODULUS, 3718.0e+6);
    material_properties.SetValue(POISSON_RATIO, 0.2);
    material_properties.SetValue(YIELD_STRESS_TENSION, 1.559e+6);
    material_properties.SetValue(FRACTURE_ENERGY_TENSION, 10.0e+2);
    material_properties.SetValue(DAMAGE_ONSET_STRESS_COMPRESSION, 10.0e+6);
    material_properties.SetValue(YIELD_STRESS_COMPRESSION, 17.99e+6);
    material_properties.SetValue(YIELD_STRAIN_COMPRESSION, 0.02);
    material_properties.SetValue(RESIDUAL_STRESS_COMPRESSION, 2.00e+6);
    material_properties.SetValue(BIAXIAL_COMPRESSION_MULTIPLIER, 1.20);
    material_properties.SetValue(FRACTURE_ENERGY_COMPRESSION, 8.0e+5);
    material_properties.SetValue(SHEAR_COMPRESSION_REDUCTOR, 0.16);
    material_properties.SetValue(BEZIER_CONTROLLER_C1, 0.65);
    material_properties.SetValue(BEZIER_CONTROLLER_C2,  0.45);
    material_properties.SetValue(BEZIER_CONTROLLER_C3, 1.5);
    material_properties.SetValue(INTEGRATION_IMPLEX, 0);
    material_properties.SetValue(TENSION_YIELD_MODEL, 0);
    material_properties.SetValue(DELAY_TIME, 1.0);

    Flags cl_options;
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    cl_parameters.SetElementGeometry(Geom);
    cl_parameters.SetMaterialProperties(material_properties);
    cl_parameters.SetStrainVector(strain_vector);
    cl_parameters.SetStressVector(stress_vector);
    cl_parameters.SetProcessInfo(process_info);
    cl_parameters.SetOptions(cl_options);
    Matrix const_matrix;
    cl_parameters.SetConstitutiveMatrix(const_matrix);

    // Create the CL
    DamageDPlusDMinusMasonry2DLaw masonry2d_cl = DamageDPlusDMinusMasonry2DLaw();

    //Initialize Material
    const Vector& r_dummy = ZeroVector(3);
    masonry2d_cl.InitializeMaterial(material_properties, Geom, r_dummy);

    std::vector<double> masonry2d_res;
    masonry2d_res = {-1.85669e+07,  -3.71337e+06, 0};

    Vector test_masonry2d_stress;
    masonry2d_cl.CalculateMaterialResponseCauchy(cl_parameters);
    test_masonry2d_stress = cl_parameters.GetStressVector();

    // Check the results
    KRATOS_EXPECT_VECTOR_NEAR(test_masonry2d_stress, masonry2d_res, 0.0001e6);
}

KRATOS_TEST_CASE_IN_SUITE(DPlusDMinusMasonry2DPMixedStateTest, KratosConstitutiveLawsFastSuite)
{
    ConstitutiveLaw::Parameters cl_parameters;
    Properties material_properties;
    ProcessInfo process_info;
    Vector stress_vector, strain_vector;

    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("Main");

    NodeType::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    NodeType::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.5, 0.0);
    NodeType::Pointer p_node_3 = test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);

    Triangle2D3<NodeType> Geom = Triangle2D3<NodeType>(p_node_1, p_node_2, p_node_3);

    stress_vector = ZeroVector(3);
    strain_vector = ZeroVector(3);
    strain_vector[0] = 0.0;
    strain_vector[1] = 0.0;
    strain_vector[2] = -0.0029;

    material_properties.SetValue(YOUNG_MODULUS, 3718.0e+6);
    material_properties.SetValue(POISSON_RATIO, 0.2);
    material_properties.SetValue(YIELD_STRESS_TENSION, 1.559e+6);
    material_properties.SetValue(FRACTURE_ENERGY_TENSION, 10.0e+2);
    material_properties.SetValue(DAMAGE_ONSET_STRESS_COMPRESSION, 10.0e+6);
    material_properties.SetValue(YIELD_STRESS_COMPRESSION, 17.99e+6);
    material_properties.SetValue(YIELD_STRAIN_COMPRESSION, 0.02);
    material_properties.SetValue(RESIDUAL_STRESS_COMPRESSION, 2.00e+6);
    material_properties.SetValue(BIAXIAL_COMPRESSION_MULTIPLIER, 1.20);
    material_properties.SetValue(FRACTURE_ENERGY_COMPRESSION, 8.0e+5);
    material_properties.SetValue(SHEAR_COMPRESSION_REDUCTOR, 0.16);
    material_properties.SetValue(BEZIER_CONTROLLER_C1, 0.65);
    material_properties.SetValue(BEZIER_CONTROLLER_C2,  0.45);
    material_properties.SetValue(BEZIER_CONTROLLER_C3, 1.5);
    material_properties.SetValue(INTEGRATION_IMPLEX, 0);
    material_properties.SetValue(TENSION_YIELD_MODEL, 0);
    material_properties.SetValue(DELAY_TIME, 1.0);

    Flags cl_options;
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    cl_parameters.SetElementGeometry(Geom);
    cl_parameters.SetMaterialProperties(material_properties);
    cl_parameters.SetStrainVector(strain_vector);
    cl_parameters.SetStressVector(stress_vector);
    cl_parameters.SetProcessInfo(process_info);
    cl_parameters.SetOptions(cl_options);
    Matrix const_matrix;
    cl_parameters.SetConstitutiveMatrix(const_matrix);

    // Create the CL
    DamageDPlusDMinusMasonry2DLaw masonry2d_cl = DamageDPlusDMinusMasonry2DLaw();

    //Initialize Material
    const Vector& r_dummy = ZeroVector(3);
    masonry2d_cl.InitializeMaterial(material_properties, Geom, r_dummy);

    std::vector<double> masonry2d_res;
    masonry2d_res = {-1.92013e+06, -1.92013e+06, -2.26567e+06};

    Vector test_masonry2d_stress;
    masonry2d_cl.CalculateMaterialResponseCauchy(cl_parameters);
    test_masonry2d_stress = cl_parameters.GetStressVector();

    // Check the results
    KRATOS_EXPECT_VECTOR_NEAR(test_masonry2d_stress, masonry2d_res, 0.0001e6);
}

} // namespace Testing
} // namespace Kratos
