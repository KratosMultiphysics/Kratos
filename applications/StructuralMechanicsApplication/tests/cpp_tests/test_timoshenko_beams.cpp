// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Gennady Markelov
//

// Project includes
#include "containers/model.h"
#include "structural_mechanics_fast_suite.h"
#include "structural_mechanics_application_variables.h"

#include <utility>
#include <boost/numeric/ublas/assignment.hpp>

namespace Kratos::Testing
{

void CreateBeamModel2N(std::string TimoshenkoBeamElementName)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(ROTATION_Z);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        constexpr auto youngs_modulus = 2.0e+06;
        p_elem_prop->SetValue(YOUNG_MODULUS, youngs_modulus);
        constexpr auto cross_area = 1.0;
        p_elem_prop->SetValue(CROSS_AREA, cross_area);
        constexpr auto i33 = 1.0;
        p_elem_prop->SetValue(I33, i33);
        constexpr auto area_effective_y = 1.0;
        p_elem_prop->SetValue(AREA_EFFECTIVE_Y, area_effective_y);

        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("TimoshenkoBeamElasticConstitutiveLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        constexpr double directional_length = 2.0;
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, directional_length, directional_length, directional_length);

        for (auto& r_node : r_model_part.Nodes()){
            r_node.AddDof(DISPLACEMENT_X);
            r_node.AddDof(DISPLACEMENT_Y);
            r_node.AddDof(DISPLACEMENT_Z);
            r_node.AddDof(ROTATION_Z);
        }

        std::vector<ModelPart::IndexType> element_nodes {1,2};
        auto p_element = r_model_part.CreateNewElement(std::move(TimoshenkoBeamElementName), 1, element_nodes, p_elem_prop);
        const auto& r_process_info = r_model_part.GetProcessInfo();
        p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law

        constexpr auto induced_strain = 0.1;
        p_element->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT) += ScalarVector(3, induced_strain * directional_length);
        p_element->GetGeometry()[1].FastGetSolutionStepValue(ROTATION_Z) += 0.1;

        std::vector<Vector> stress_vectors;
        p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stress_vectors, r_process_info);

        constexpr auto expected_stress = induced_strain * youngs_modulus;
        constexpr auto expected_shear_stress = -37500.0;

        constexpr auto tolerance = 1.0e-5;
        Vector expected_stress_vector(3);
        expected_stress_vector <<= expected_stress, 29631.5,expected_shear_stress;
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress_vector, stress_vectors[0], tolerance);
        expected_stress_vector <<= expected_stress, 70710.7,expected_shear_stress;
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress_vector, stress_vectors[1], tolerance);
        expected_stress_vector <<= expected_stress, 111790,expected_shear_stress;
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress_vector, stress_vectors[2], tolerance);

        Vector pre_stress(3);
        pre_stress <<= 1.0e5, 1.0e4, 1.0e3;
        p_element->GetProperties().SetValue(TIMOSHENKO_BEAM_PRESTRESS_PK2, pre_stress);
        p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stress_vectors, r_process_info);
        expected_stress_vector += pre_stress;
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress_vector, stress_vectors[2], tolerance);
    }

    KRATOS_TEST_CASE_IN_SUITE(LinearTimoshenkoBeam2D2N_CalculatesPK2Stress, KratosStructuralMechanicsFastSuite)
    {
        CreateBeamModel2N("LinearTimoshenkoBeamElement2D2N");
    }
}
