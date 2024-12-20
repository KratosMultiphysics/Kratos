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

template<SizeType TNNodes>
void Create2DBeamModel_and_CheckPK2Stress(std::string TimoshenkoBeamElementName, const std::vector<double>& rExpectedShearStress, const std::vector<double>& rExpectedBendingMoment)
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
        p_elem_prop->SetValue(CROSS_AREA, 1.0);
        p_elem_prop->SetValue(I33, 1.0);
        p_elem_prop->SetValue(AREA_EFFECTIVE_Y, 5.0/6.0);

        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("TimoshenkoBeamElasticConstitutiveLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        constexpr double directional_length = 2.0;
        r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        r_model_part.CreateNewNode(2, directional_length, directional_length, 0.0);
        std::vector<ModelPart::IndexType> element_nodes{1,2};
        if(TNNodes==3) {
            r_model_part.CreateNewNode(3, directional_length/2, directional_length/2, 0.0);
            element_nodes.push_back(3);
        }

        for (auto& r_node : r_model_part.Nodes()){
            r_node.AddDof(DISPLACEMENT_X);
            r_node.AddDof(DISPLACEMENT_Y);
            r_node.AddDof(ROTATION_Z);
        }

        auto p_element = r_model_part.CreateNewElement(std::move(TimoshenkoBeamElementName), 1, element_nodes, p_elem_prop);
        const auto& r_process_info = r_model_part.GetProcessInfo();
        p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law

        constexpr auto induced_strain = 0.1;
        p_element->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT) += ScalarVector(2, induced_strain * directional_length);
        p_element->GetGeometry()[1].FastGetSolutionStepValue(ROTATION_Z) += 0.1;
        if(TNNodes==3) {
            p_element->GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT) += ScalarVector(2, induced_strain * directional_length/2.0);
            p_element->GetGeometry()[2].FastGetSolutionStepValue(ROTATION_Z) += 0.1/2.0;
        }

        std::vector<Vector> stress_vectors;
        p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stress_vectors, r_process_info);

        constexpr auto expected_stress = induced_strain * youngs_modulus;
        constexpr auto tolerance = 1.0e-5;
        Vector expected_stress_vector(3);
        expected_stress_vector <<= expected_stress, rExpectedBendingMoment[0],rExpectedShearStress[0];
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress_vector, stress_vectors[0], tolerance);
        expected_stress_vector <<= expected_stress, rExpectedBendingMoment[1],rExpectedShearStress[1];
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress_vector, stress_vectors[1], tolerance);
        if(stress_vectors.size()>2) {
            expected_stress_vector <<= expected_stress, rExpectedBendingMoment[2],rExpectedShearStress[2];
            KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress_vector, stress_vectors[2], tolerance);
        }

        expected_stress_vector <<= expected_stress, rExpectedBendingMoment[0],rExpectedShearStress[0];
        Vector pre_stress(3);
        pre_stress <<= 1.0e5, 1.0e4, 1.0e3;
        p_element->GetProperties().SetValue(BEAM_PRESTRESS_PK2, pre_stress);
        p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stress_vectors, r_process_info);
        expected_stress_vector += pre_stress;
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress_vector, stress_vectors[0], tolerance);
    }

    KRATOS_TEST_CASE_IN_SUITE(LinearTimoshenkoBeam2D2N_CalculatesPK2Stress, KratosStructuralMechanicsFastSuite)
    {
        const std::vector expected_shear_stress {-32608.7, -32608.7, -32608.7};
        const std::vector expected_bending_moment {34989.6, 70710.7, 106432.0};

        Create2DBeamModel_and_CheckPK2Stress<2>("LinearTimoshenkoBeamElement2D2N", expected_shear_stress, expected_bending_moment);
    }

    KRATOS_TEST_CASE_IN_SUITE(LinearTimoshenkodBeam2D3N_CalculatesPK2Stress, KratosStructuralMechanicsFastSuite)
    {
       const std::vector expected_shear_stress {2604.74, -21673.6, -48084.9};
       const std::vector expected_bending_moment {86720.4, 68516.6, 61527.5};
       Create2DBeamModel_and_CheckPK2Stress<3>("LinearTimoshenkoBeamElement2D3N", expected_shear_stress, expected_bending_moment);
    }

    KRATOS_TEST_CASE_IN_SUITE(LinearTimoshenkodCurvedBeam2D3N_CalculatesPK2Stress, KratosStructuralMechanicsFastSuite)
    {
        const std::vector expected_shear_stress {17610.4, 65722.9};
        const std::vector expected_bending_moment {70710.7, 70710.7};
        Create2DBeamModel_and_CheckPK2Stress<3>("LinearTimoshenkoCurvedBeamElement2D3N", expected_shear_stress, expected_bending_moment);
    }
}
