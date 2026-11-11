// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

// Project includes
#include "containers/model.h"
#include "structural_mechanics_fast_suite.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos::Testing {

namespace {

    /**
     * \brief Sets up a nodal concentrated element with nodal displacement, volume acceleration, mass and stiffness
     * \param  rModelPart Current model part
     * \param  rElementName Registered name of the element (NodalConcentratedElement2D1N or NodalConcentratedElement3D1N)
     * \return NodalConcentratedElement
     */
    Element::Pointer SetUpNodalConcentratedElement(ModelPart& rModelPart, const std::string& rElementName)
    {
        rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
        rModelPart.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

        // Set the element properties
        const auto p_elem_prop = rModelPart.CreateNewProperties(0);

        // Create the test element
        auto p_node = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
        p_node->FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{ 5.0, 6.0, 7.0 };
        p_node->FastGetSolutionStepValue(VOLUME_ACCELERATION) = array_1d<double, 3>{ 1.0, -9.81, 2.0 };

        const std::vector<ModelPart::IndexType> element_nodes{ 1 };
        auto p_element = rModelPart.CreateNewElement(rElementName, 1, element_nodes, p_elem_prop);

        p_element->SetValue(NODAL_MASS, 2.5);
        p_element->SetValue(NODAL_DISPLACEMENT_STIFFNESS, array_1d<double, 3>{ 2.0, 3.0, 4.0 });

        return p_element;
    }

}

    // Tests the internal forces of the NodalConcentratedElement2D1N
    KRATOS_TEST_CASE_IN_SUITE(NodalConcentratedElementInternalForces2D, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
        auto p_element = SetUpNodalConcentratedElement(r_model_part, "NodalConcentratedElement2D1N");

        Vector calculated_internal_forces;
        p_element->Calculate(INTERNAL_FORCES_VECTOR, calculated_internal_forces, r_model_part.GetProcessInfo());

        // Fint = K * u, the z-component is not part of a 2D element
        Vector expected_internal_forces = ZeroVector(2);
        expected_internal_forces(0) = 2.0 * 5.0;
        expected_internal_forces(1) = 3.0 * 6.0;

        KRATOS_EXPECT_VECTOR_NEAR(expected_internal_forces, calculated_internal_forces, 1e-10);
    }

    // Tests the internal forces of the NodalConcentratedElement3D1N
    KRATOS_TEST_CASE_IN_SUITE(NodalConcentratedElementInternalForces3D, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
        auto p_element = SetUpNodalConcentratedElement(r_model_part, "NodalConcentratedElement3D1N");

        Vector calculated_internal_forces;
        p_element->Calculate(INTERNAL_FORCES_VECTOR, calculated_internal_forces, r_model_part.GetProcessInfo());

        // Fint = K * u
        Vector expected_internal_forces = ZeroVector(3);
        expected_internal_forces(0) = 2.0 * 5.0;
        expected_internal_forces(1) = 3.0 * 6.0;
        expected_internal_forces(2) = 4.0 * 7.0;

        KRATOS_EXPECT_VECTOR_NEAR(expected_internal_forces, calculated_internal_forces, 1e-10);
    }

    // Tests the external forces of the NodalConcentratedElement2D1N
    KRATOS_TEST_CASE_IN_SUITE(NodalConcentratedElementExternalForces2D, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
        auto p_element = SetUpNodalConcentratedElement(r_model_part, "NodalConcentratedElement2D1N");

        Vector calculated_external_forces;
        p_element->Calculate(EXTERNAL_FORCES_VECTOR, calculated_external_forces, r_model_part.GetProcessInfo());

        // Fext = m * VOLUME_ACCELERATION, the z-component is not part of a 2D element
        Vector expected_external_forces = ZeroVector(2);
        expected_external_forces(0) = 2.5 * 1.0;
        expected_external_forces(1) = 2.5 * -9.81;

        KRATOS_EXPECT_VECTOR_NEAR(expected_external_forces, calculated_external_forces, 1e-10);
    }

    // Tests the external forces of the NodalConcentratedElement3D1N
    KRATOS_TEST_CASE_IN_SUITE(NodalConcentratedElementExternalForces3D, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
        auto p_element = SetUpNodalConcentratedElement(r_model_part, "NodalConcentratedElement3D1N");

        Vector calculated_external_forces;
        p_element->Calculate(EXTERNAL_FORCES_VECTOR, calculated_external_forces, r_model_part.GetProcessInfo());

        // Fext = m * VOLUME_ACCELERATION
        Vector expected_external_forces = ZeroVector(3);
        expected_external_forces(0) = 2.5 * 1.0;
        expected_external_forces(1) = 2.5 * -9.81;
        expected_external_forces(2) = 2.5 * 2.0;

        KRATOS_EXPECT_VECTOR_NEAR(expected_external_forces, calculated_external_forces, 1e-10);
    }

    // Tests that the rhs of the NodalConcentratedElement3D1N equals Fext - Fint
    KRATOS_TEST_CASE_IN_SUITE(NodalConcentratedElementRightHandSideEqualsExternalMinusInternalForces, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
        auto p_element = SetUpNodalConcentratedElement(r_model_part, "NodalConcentratedElement3D1N");

        const auto& r_process_info = r_model_part.GetProcessInfo();

        Vector internal_forces;
        p_element->Calculate(INTERNAL_FORCES_VECTOR, internal_forces, r_process_info);

        Vector external_forces;
        p_element->Calculate(EXTERNAL_FORCES_VECTOR, external_forces, r_process_info);

        Vector calculated_rhs;
        p_element->CalculateRightHandSide(calculated_rhs, r_process_info);

        const Vector expected_rhs = external_forces - internal_forces;

        KRATOS_EXPECT_VECTOR_NEAR(expected_rhs, calculated_rhs, 1e-10);
    }

    // Tests that calculating an unsupported vector variable throws an error
    KRATOS_TEST_CASE_IN_SUITE(NodalConcentratedElementCalculateUnsupportedVectorVariableThrows, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
        auto p_element = SetUpNodalConcentratedElement(r_model_part, "NodalConcentratedElement3D1N");

        Vector output;
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            p_element->Calculate(RESIDUAL_VECTOR, output, r_model_part.GetProcessInfo()),
            "Variable RESIDUAL_VECTOR not supported in element Element #1");
    }
}
