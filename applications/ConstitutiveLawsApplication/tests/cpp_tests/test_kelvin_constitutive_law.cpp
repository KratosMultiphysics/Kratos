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
//

// System includes

// External includes

// Project includes
#include "includes/process_info.h"
#include "containers/model.h"

// Application includes
#include "tests/cpp_tests/constitutive_laws_fast_suite.h"

// Constitutive law
#include "custom_constitutive/small_strains/viscous/viscous_generalized_kelvin.h"
#include "includes/model_part.h"
#include "geometries/tetrahedra_3d_4.h"
#include "constitutive_laws_application_variables.h"

namespace Kratos::Testing
{

/**
* Check the correct calculation of the integrated stress with the CL's
*/

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawKelvinInternalVariables, KratosConstitutiveLawsFastSuite)
{
    //
    // Test: check correct behavior of internal and calculated variables
    //

    ViscousGeneralizedKelvin<ElasticIsotropic3D> cl = ViscousGeneralizedKelvin<ElasticIsotropic3D>();

    KRATOS_EXPECT_FALSE(cl.Has(INTEGRATED_STRESS_TENSOR));  // = False, in order to use CalculateValue())

    // This constitutive law does not use internal variables
    // TODO (marandra): check that this is compatible con API
    KRATOS_EXPECT_FALSE(cl.Has(INTERNAL_VARIABLES));  // = False
}


KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawKelvin, KratosConstitutiveLawsFastSuite)
{
    ConstitutiveLaw::Parameters cl_parameters;
    Properties material_properties;
    ProcessInfo process_info;
    Vector stress_vector, strain_vector;

    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("Main");

    Node::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    Node::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    Node::Pointer p_node_3 = test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    Node::Pointer p_node_4 = test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);

    Tetrahedra3D4<Node> Geom = Tetrahedra3D4<Node>(p_node_1, p_node_2, p_node_3, p_node_4);

    stress_vector = ZeroVector(6);
    strain_vector = ZeroVector(6);
    strain_vector[0] = 0.0;
    strain_vector[1] = 0.0;
    strain_vector[2] = -8.0e-5;
    strain_vector[3] = 0.0;
    strain_vector[4] = 0.0;
    strain_vector[5] = -1.6941e-21;

    material_properties.SetValue(YOUNG_MODULUS, 210e9);
    material_properties.SetValue(POISSON_RATIO, 0.22);
    material_properties.SetValue(DELAY_TIME, 1.0);

    Flags cl_options;
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    process_info.SetValue(DELTA_TIME, 0.1);

    cl_parameters.SetElementGeometry(Geom);
    cl_parameters.SetMaterialProperties(material_properties);
    cl_parameters.SetStrainVector(strain_vector);
    cl_parameters.SetStressVector(stress_vector);
    cl_parameters.SetProcessInfo(process_info);
    cl_parameters.SetOptions(cl_options);
    Matrix const_matrix;
    cl_parameters.SetConstitutiveMatrix(const_matrix);

    // Create the CL
    ViscousGeneralizedKelvin<ElasticIsotropic3D> kelvin_cl = ViscousGeneralizedKelvin<ElasticIsotropic3D>();

    std::vector<double> kelvin_res;
    kelvin_res = {-53828.8, -53828.8, -190847, 0, 0, -1.45077e-12};

    Vector test_kelvin_stress;
    kelvin_cl.CalculateMaterialResponseCauchy(cl_parameters);
    test_kelvin_stress = cl_parameters.GetStressVector();

    // Check the results
    KRATOS_EXPECT_VECTOR_NEAR(test_kelvin_stress, kelvin_res, 1.0);
}

} // namespace Kratos::Testing
