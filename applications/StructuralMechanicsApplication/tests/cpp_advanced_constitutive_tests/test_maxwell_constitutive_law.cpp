// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
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
#include "testing/testing.h"
#include "containers/model.h"
// Application includes

// Constitutive law
#include "custom_advanced_constitutive/viscous_generalized_maxwell.h"
#include "includes/model_part.h"
#include "geometries/tetrahedra_3d_4.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
namespace Testing
{
// We test the associated plasticity Constitutive laws...
typedef Node<3> NodeType;

/**
* Check the correct calculation of the integrated stress with the CL's
*/

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawMaxwellInternalVariables, KratosStructuralMechanicsFastSuite)
{
    //
    // Test: check correct behavior of internal and calculated variables
    //

    ViscousGeneralizedMaxwell<ElasticIsotropic3D> cl = ViscousGeneralizedMaxwell<ElasticIsotropic3D>();

    KRATOS_CHECK_IS_FALSE(cl.Has(INTEGRATED_STRESS_TENSOR));  // = False, in order to use CalculateValue())

    // This constitutive law does not use internal variables
    // TODO (marandra): check that this is compatible con API
    KRATOS_CHECK_IS_FALSE(cl.Has(INTERNAL_VARIABLES));  // = False
}


KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawMaxwell, KratosStructuralMechanicsFastSuite)
{
    ConstitutiveLaw::Parameters cl_parameters;
    Properties material_properties;
    ProcessInfo process_info;
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
    strain_vector[0] = 0.0;
    strain_vector[1] = 0.0;
    strain_vector[2] = -8.0e-5;
    strain_vector[3] = 0.0;
    strain_vector[4] = 0.0;
    strain_vector[5] = -1.6941e-21;

    material_properties.SetValue(YOUNG_MODULUS, 210e9);
    material_properties.SetValue(POISSON_RATIO, 0.22);
    material_properties.SetValue(DELAY_TIME, 1.0);
    material_properties.SetValue(VISCOUS_PARAMETER, 0.15);

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
    ViscousGeneralizedMaxwell<ElasticIsotropic3D> maxwell_cl = ViscousGeneralizedMaxwell<ElasticIsotropic3D>();

    std::vector<double> maxwell_res;
    maxwell_res = {-5.37455e+06, -5.37455e+06, -1.90552e+07, 0, 0, -1.44853e-10};

    Vector test_maxwell_stress;
    maxwell_cl.CalculateMaterialResponseCauchy(cl_parameters);
    test_maxwell_stress = cl_parameters.GetStressVector();

    // Check the results
    KRATOS_CHECK_VECTOR_NEAR(test_maxwell_stress, maxwell_res, 0.0001e6);
}
} // namespace Testing
} // namespace Kratos
