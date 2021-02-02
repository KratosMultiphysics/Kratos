//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "utilities/cpp_tests_utilities.h"

/* Processes */
#include "processes/set_initial_state_process.h"

namespace Kratos
{
    namespace Testing
    {

        /**
        * Checks the correct work of the nodal H compute
        * Test triangle
        */
        KRATOS_TEST_CASE_IN_SUITE(ImposingInitialStrain2D, KratosCoreFastSuite)
        {
            Model current_model;

            ModelPart& this_model_part = current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(2);

            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);

            Vector initial_E = ZeroVector(3);
            initial_E(0) = 0.01;
            initial_E(1) = 0.02;
            initial_E(2) = 0.03;

            Vector initial_S = ZeroVector(3);
            initial_S(0) = 1.0e6;
            initial_S(0) = 2.0e6;
            initial_S(0) = 3.0e6;

            Matrix initial_F = ZeroMatrix(2,2);
            initial_F(0,0) = 0.001;
            initial_F(0,1) = 0.0001;
            initial_F(1,0) = initial_F(0,1);
            initial_F(1,1) = 0.002;

            // Set the initial state
            auto process = SetInitialStateProcess<2>(this_model_part, initial_E, initial_S, initial_F);
            process.ExecuteInitializeSolutionStep();

            const double tolerance = 1.0e-4;
            std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector;
            for (auto i_element = this_model_part.ElementsBegin(); i_element != this_model_part.ElementsEnd(); i_element++) {
                const auto& r_integration_points = i_element->GetGeometry().IntegrationPoints(i_element->GetIntegrationMethod());
                i_element->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_law_vector, this_model_part.GetProcessInfo());
                
                for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
                    auto p_initial_state = constitutive_law_vector[point_number]->pGetpInitialState();
                    const auto& r_imposed_F = p_initial_state->GetInitialDeformationGradientMatrix();
                    const auto& r_imposed_E = p_initial_state->GetInitialStrainVector();
                    const auto& r_imposed_S = p_initial_state->GetInitialStressVector();

                    for (IndexType component = 0; component < 3; component++) {
                        KRATOS_CHECK_LESS_EQUAL(r_imposed_E(component) - initial_E(component), tolerance);
                        KRATOS_CHECK_LESS_EQUAL(r_imposed_S(component) - initial_S(component), tolerance);
                    }
                    KRATOS_CHECK_LESS_EQUAL(r_imposed_F(0,0) - initial_F(0,0), tolerance);
                    KRATOS_CHECK_LESS_EQUAL(r_imposed_F(0,1) - initial_F(0,1), tolerance);
                    KRATOS_CHECK_LESS_EQUAL(r_imposed_F(1,1) - initial_F(1,1), tolerance);
                }
            }
        }

        /**
        * Checks the correct work of the nodal H non-historical compute
        * Test triangle
        */
        KRATOS_TEST_CASE_IN_SUITE(ImposingInitialStress, KratosCoreFastSuite)
        {
            // Model current_model;

            // ModelPart& this_model_part = current_model.CreateModelPart("Main");
            // this_model_part.SetBufferSize(1);

            // CppTestsUtilities::Create2DGeometry(this_model_part);

            // // Compute NodalH
            // auto process = FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable>(this_model_part);
            // process.Execute();

            // // // DEBUG
            // // GiDIODebugNodalH(this_model_part);

            // const double tolerance = 1.0e-4;
            // KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(1)->GetValue(NODAL_H) - 1.0, tolerance);
            // KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(2)->GetValue(NODAL_H) - 1.0, tolerance);
            // KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(5)->GetValue(NODAL_H) - 1.0, tolerance);
            // KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(6)->GetValue(NODAL_H) - 1.0, tolerance);
        }

        /**
        * Checks the correct work of the nodal H compute
        * Test tetrahedra
        */
        KRATOS_TEST_CASE_IN_SUITE(ImposingInitialF, KratosCoreFastSuite)
        {
            // Model current_model;

            // ModelPart& this_model_part = current_model.CreateModelPart("Main");
            // this_model_part.SetBufferSize(2);

            // this_model_part.AddNodalSolutionStepVariable(NODAL_H);

            // auto& process_info = this_model_part.GetProcessInfo();
            // process_info[STEP] = 1;
            // process_info[NL_ITERATION_NUMBER] = 1;

            // CppTestsUtilities::Create3DGeometry(this_model_part);

            // // Compute NodalH
            // auto process = FindNodalHProcess<FindNodalHSettings::SaveAsHistoricalVariable>(this_model_part);
            // process.Execute();

            // const double tolerance = 1.0e-4;
            // KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(1)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            // KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(2)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            // KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(3)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            // KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(5)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            // KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(9)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            // KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(10)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            // KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(11)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            // KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(12)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
        }
    } // namespace Testing
}  // namespace Kratos.
