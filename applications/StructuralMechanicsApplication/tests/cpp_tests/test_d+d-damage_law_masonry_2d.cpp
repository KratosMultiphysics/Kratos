// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
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

// Constitutive law
#include "custom_constitutive/plane_stress_d_plus_d_minus_damage_masonry_2d.h"
#include "includes/model_part.h"
#include "geometries/triangle_2d_3.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
namespace Testing
{
// We test the d+ d- damage Constitutive Law (for masonry) 2D...
typedef Node<3> NodeType;

/**
    * Check the correct calculation of the integrated stress with the CL's
    */
KRATOS_TEST_CASE_IN_SUITE(DPlusDMinusMasonry2DPureTensionTest, KratosStructuralMechanicsFastSuite)
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
    strain_vector[0] = 0.000787;
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
    material_properties.SetValue(DELAY_TIME, 1.0);

    Flags cl_options;
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    //process_info.SetValue(DELTA_TIME, 0.1);

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

    std::vector<double> masonry2d_res;
    masonry2d_res = {926316.0, 185263.0, 0};

    Vector test_masonry2d_stress;
    masonry2d_cl.CalculateMaterialResponseCauchy(cl_parameters);
    test_masonry2d_stress = cl_parameters.GetStressVector();

    // Check the results
    for (int comp = 0; comp < 3; comp++) {
        // KRATOS_CHECK_NEAR(test_masonry2d_stress[comp], masonry2d_res[comp], 0.0001e6);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DPlusDMinusMasonry2DPureCompressionTest, KratosStructuralMechanicsFastSuite)
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
    strain_vector[0] = -0.02227;
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
    material_properties.SetValue(DELAY_TIME, 1.0);

    Flags cl_options;
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    //process_info.SetValue(DELTA_TIME, 0.1);

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

    std::vector<double> masonry2d_res;
    masonry2d_res = {-2.06954e+07,  -4.13908e+06, 0};

    Vector test_masonry2d_stress;
    masonry2d_cl.CalculateMaterialResponseCauchy(cl_parameters);
    test_masonry2d_stress = cl_parameters.GetStressVector();

    // Check the results
    for (int comp = 0; comp < 3; comp++) {
        KRATOS_CHECK_NEAR(test_masonry2d_stress[comp], masonry2d_res[comp], 0.0001e6);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DPlusDMinusMasonry2DPMixedStateTest, KratosStructuralMechanicsFastSuite)
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
    strain_vector[2] = -0.005218;
 
    material_properties.SetValue(YOUNG_MODULUS, 3718.0e+6);
    material_properties.SetValue(POISSON_RATIO, 0.2);
    material_properties.SetValue(YIELD_STRESS_TENSION, 1.559e+6);
    material_properties.SetValue(FRACTURE_ENERGY_TENSION, 10.0e+3);
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
    material_properties.SetValue(DELAY_TIME, 1.0);

    Flags cl_options;
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    //process_info.SetValue(DELTA_TIME, 0.1);

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

    std::vector<double> masonry2d_res;
    masonry2d_res = {-1.57685e+06, -1.57685e+06, -5.93494e+06};

    Vector test_masonry2d_stress;
    masonry2d_cl.CalculateMaterialResponseCauchy(cl_parameters);
    test_masonry2d_stress = cl_parameters.GetStressVector();

    // Check the results
    for (int comp = 0; comp < 3; comp++) {
        KRATOS_CHECK_NEAR(test_masonry2d_stress[comp], masonry2d_res[comp], 0.0001e6);
    }
}


} // namespace Testing
} // namespace Kratos
