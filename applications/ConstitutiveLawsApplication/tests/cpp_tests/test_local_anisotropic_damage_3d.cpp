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
//  Main authors:    Alejandro Cornejo
//                   Marcelo Raschi
//

// System includes

// External includes

// Project includes
#include "includes/process_info.h"
#include "testing/testing.h"
#include "containers/model.h"

// Application includes
#include "structural_mechanics_application_variables.h"
#include "constitutive_laws_application_variables.h"

// Constitutive law
#include "custom_constitutive/elastic_local_anisotropic_damage_3d.h"
#include "includes/model_part.h"
#include "geometries/hexahedra_3d_8.h"

namespace Kratos
{
namespace Testing
{
typedef Node<3> NodeType;
KRATOS_TEST_CASE_IN_SUITE(LocalElasticAnisotropicVariables, KratosConstitutiveLawsFastSuite)
{
    // Create the CL
    ElasticAnisotropicDamage cl = ElasticAnisotropicDamage();

    //  TESTS
    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("Main");

    KRATOS_CHECK_IS_FALSE(cl.Has(STRAIN_ENERGY));  // = False, in order to use CalculateValue())
    KRATOS_CHECK_IS_FALSE(cl.Has(INTERNAL_VARIABLES));  // = False, in order to use CalculateValue())
    KRATOS_CHECK_IS_FALSE(cl.Has(DAMAGE_VARIABLE));  // = False
    KRATOS_CHECK(cl.Has(DAMAGE_VECTOR)); //True

    Vector damage_vector = ZeroVector(3);
    damage_vector[0] = 0.106572;
    damage_vector[0] = 0.21515445;
    damage_vector[0] = 0.25454542;
    cl.SetValue(DAMAGE_VECTOR, damage_vector, test_model_part.GetProcessInfo());
    Vector damage_vector_r;
    cl.GetValue(DAMAGE_VECTOR, damage_vector_r);
    for (std::size_t i = 0; i < 3; ++i)
        KRATOS_CHECK_NEAR(damage_vector_r[i], damage_vector[i], 1.e-15);

}
KRATOS_TEST_CASE_IN_SUITE(LocalElasticAnisotropicDamage, KratosConstitutiveLawsFastSuite)
{
    ConstitutiveLaw::Parameters cl_parameters;
    Properties material_properties;
    Vector stress_vector(6), strain_ref(6);
    Matrix const_matrix(6,6);
    // Create gauss point
    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("Main");
    NodeType::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    NodeType::Pointer p_node_2 = test_model_part.CreateNewNode(2, 0.2, 0.0, 0.0);
    NodeType::Pointer p_node_3 = test_model_part.CreateNewNode(3, 0.2, 0.2, 0.0);
    NodeType::Pointer p_node_4 = test_model_part.CreateNewNode(4, 0.0, 0.2, 0.0);
    NodeType::Pointer p_node_5 = test_model_part.CreateNewNode(5, 0.0, 0.0, 0.2);
    NodeType::Pointer p_node_6 = test_model_part.CreateNewNode(6, 0.2, 0.0, 0.2);
    NodeType::Pointer p_node_7 = test_model_part.CreateNewNode(7, 0.2, 0.2, 0.2);
    NodeType::Pointer p_node_8 = test_model_part.CreateNewNode(8, 0.0, 0.2, 0.2);
    Hexahedra3D8<NodeType> geometry = Hexahedra3D8<NodeType>(p_node_1, p_node_2, p_node_3, p_node_4, p_node_5, p_node_6, p_node_7, p_node_8);

    // Set material properties
    material_properties.SetValue(YOUNG_MODULUS, 31700.0);
    material_properties.SetValue(POISSON_RATIO, 0.2);
    material_properties.SetValue(YIELD_STRESS_TENSION, 3.48);
    material_properties.SetValue(YIELD_STRESS_COMPRESSION, 27.6);
    material_properties.SetValue(CHARACTERISTIC_INTERNAL_LENGTH, 0.2);
    material_properties.SetValue(DAMAGE_MODEL_PARAMETER_BETA1_TENSION, 0.85);
    material_properties.SetValue(DAMAGE_MODEL_PARAMETER_BETA2_TENSION, 0.18);
    material_properties.SetValue(DAMAGE_MODEL_PARAMETER_BETA1_COMPRESSION, 0.0);
    material_properties.SetValue(DAMAGE_MODEL_PARAMETER_BETA2_COMPRESSION, 0.095);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=cl_parameters.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Set required constitutive law parameters:
    cl_parameters.SetElementGeometry(geometry);
    cl_parameters.SetProcessInfo(test_model_part.GetProcessInfo());
    cl_parameters.SetMaterialProperties(material_properties);
    cl_parameters.SetStrainVector(strain_ref);
    cl_parameters.SetStressVector(stress_vector);
    cl_parameters.SetConstitutiveMatrix(const_matrix);

    // Create the CL
    ElasticAnisotropicDamage cl = ElasticAnisotropicDamage();

    // Check CalculateValue(STRAIN)
    strain_ref(0) = -0.00341;
    strain_ref(1) = -0.00341;
    strain_ref(2) = 0.0002;
    strain_ref(3) = 0.0;
    strain_ref(4) = 0.0;
    strain_ref(5) = 0.000131555;

    Vector strain_r;
    cl.CalculateValue(cl_parameters, STRAIN, strain_r);
    KRATOS_CHECK_VECTOR_NEAR(strain_r, strain_ref, 1e-15);

    //Check CalculateVector(STRESSES)
    Vector stress_r;
    cl.CalculateValue(cl_parameters, STRESSES, stress_r);
    KRATOS_CHECK_VECTOR_NEAR(stress_r, stress_vector, 1e-15);

    //
    // Test: damage evolution law
    //
    const double tolerance = 1e-5;
    Vector ref_damage_vector = ZeroVector(3);
    ref_damage_vector[0] = 0.554592;
    ref_damage_vector[1] = 0.554592;
    ref_damage_vector[2] = 0.207461;

    Vector ref_stress = ZeroVector(6);

    ref_stress(0) = -29.623288881;
    ref_stress(1) = -29.623288946;
    ref_stress(2) = -18.083304966;
    ref_stress(3) = 0.0;
    ref_stress(4) = 0.001087924;
    ref_stress(5) = 0.000000011;
    /* Matrix ref_C = ZeroMatrix(6, 6);
    ref_C(0,0) =  3.1672e+00; ref_C(0,1) =  6.2137e-01; ref_C(0,2) =  9.8019e-01; ref_C(0,3) = -1.7259e-01; ref_C(0,4) = -2.0169e-03; ref_C(0,5) =  1.0539e-01;
    ref_C(1,0) =  6.2137e-01; ref_C(1,1) =  1.3478e+00; ref_C(1,2) =  2.0793e-01; ref_C(1,3) = -4.2244e-01; ref_C(1,4) = -4.9367e-03; ref_C(1,5) =  2.5796e-01;
    ref_C(2,0) =  9.8019e-01; ref_C(2,1) =  2.0793e-01; ref_C(2,2) =  2.7511e+00; ref_C(2,3) = -2.5256e-01; ref_C(2,4) = -2.9514e-03; ref_C(2,5) =  1.5422e-01;
    ref_C(3,0) = -1.7259e-01; ref_C(3,1) = -4.2244e-01; ref_C(3,2) = -2.5256e-01; ref_C(3,3) =  9.2737e-01; ref_C(3,4) = -9.5492e-04; ref_C(3,5) =  4.9898e-02;
    ref_C(4,0) = -2.0169e-03; ref_C(4,1) = -4.9367e-03; ref_C(4,2) = -2.9514e-03; ref_C(4,3) = -9.5492e-04; ref_C(4,4) =  1.0090e+00; ref_C(4,5) =  5.8312e-04;
    ref_C(5,0) =  1.0539e-01; ref_C(5,1) =  2.5796e-01; ref_C(5,2) =  1.5422e-01; ref_C(5,3) =  4.9898e-02; ref_C(5,4) =  5.8312e-04; ref_C(5,5) =  9.7862e-01;
 */
    // Simulate the call sequence of the element
    Vector dummy_vector;
    cl.InitializeMaterial(material_properties, geometry, dummy_vector);

    if (cl.RequiresInitializeMaterialResponse()){
        cl.InitializeMaterialResponseCauchy(cl_parameters);
    }

    cl.CalculateMaterialResponsePK2(cl_parameters);

    if (cl.RequiresFinalizeMaterialResponse()){
        cl.FinalizeMaterialResponseCauchy(cl_parameters);
    }

    // Check damage variable and strain energy
    Vector values;
    cl.CalculateValue(cl_parameters, DAMAGE_VECTOR, values);
    for (std::size_t comp = 0; comp < 3; ++comp)
        KRATOS_CHECK_NEAR(ref_damage_vector(comp), values(comp), tolerance);

    // Check stress
    for (std::size_t comp = 0; comp < 6; ++comp)
        KRATOS_CHECK_NEAR(stress_vector(comp), ref_stress(comp), tolerance);
/*
    // Check constitutive tensor
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            KRATOS_CHECK_NEAR(const_matrix(i, j), ref_C(j, i), tolerance); */

}
} // namespace Testing
} // namespace Kratos
