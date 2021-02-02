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
        KRATOS_TEST_CASE_IN_SUITE(ImposingInitialStrain, KratosCoreFastSuite)
        {
            // Model current_model;

            // ModelPart& this_model_part = current_model.CreateModelPart("Main");
            // this_model_part.SetBufferSize(2);

            // this_model_part.AddNodalSolutionStepVariable(NODAL_H);

            // auto& process_info = this_model_part.GetProcessInfo();
            // process_info[STEP] = 1;
            // process_info[NL_ITERATION_NUMBER] = 1;

            // CppTestsUtilities::Create2DGeometry(this_model_part);

            // // Compute NodalH
            // auto process = FindNodalHProcess<FindNodalHSettings::SaveAsHistoricalVariable>(this_model_part);
            // process.Execute();

            // const double tolerance = 1.0e-4;
            // KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(1)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            // KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(2)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            // KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(5)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            // KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(6)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
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
