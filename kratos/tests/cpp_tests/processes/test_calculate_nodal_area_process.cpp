//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "utilities/cpp_tests_utilities.h"
#include "testing/testing.h"
#include "includes/kratos_flags.h"
// #include "includes/gid_io.h"

/* Processes */
#include "processes/calculate_nodal_area_process.h"

namespace Kratos::Testing 
{
// void GiDIODebugNodalArea(ModelPart& ThisModelPart)
// {
//     GidIO<> gid_io("TEST_NODAL_AREA", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
//     const int nl_iter = ThisModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
//     const double label = static_cast<double>(nl_iter);

//     gid_io.InitializeMesh(label);
//     gid_io.WriteMesh(ThisModelPart.GetMesh());
//     gid_io.FinalizeMesh();
//     gid_io.InitializeResults(label, ThisModelPart.GetMesh());
//     gid_io.WriteNodalResults(NODAL_AREA, ThisModelPart.Nodes(), label, 0);
// }

/**
* Checks the correct work of the nodal gradient compute
* Test triangle
*/
KRATOS_TEST_CASE_IN_SUITE(NodalArea1, KratosCoreFastSuite)
{
    Model current_model;

    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.SetBufferSize(2);

    this_model_part.AddNodalSolutionStepVariable(NODAL_AREA);

    auto& process_info = this_model_part.GetProcessInfo();
    process_info[STEP] = 1;
    process_info[NL_ITERATION_NUMBER] = 1;

    CppTestsUtilities::Create2DGeometry(this_model_part);

    using HistNodalAreaProcess = CalculateNodalAreaProcess<CalculateNodalAreaSettings::SaveAsHistoricalVariable>; ;
    HistNodalAreaProcess hist_process(this_model_part, 2);
    hist_process.Execute();

    // // DEBUG
    // GiDIODebugNodalArea(this_model_part);

    const double tolerance = 1.0e-8;
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(1)->FastGetSolutionStepValue(NODAL_AREA)) - 1.0/3.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(2)->FastGetSolutionStepValue(NODAL_AREA)) - 1.0/3.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(5)->FastGetSolutionStepValue(NODAL_AREA)) - 1.0/3.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(6)->FastGetSolutionStepValue(NODAL_AREA)) - 1.0/6.0, tolerance);

    using NonHistNodalAreaProcess= CalculateNodalAreaProcess<CalculateNodalAreaSettings::SaveAsNonHistoricalVariable>; ;
    NonHistNodalAreaProcess nonhist_process(this_model_part, 2);
    nonhist_process.Execute();

    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(1)->GetValue(NODAL_AREA)) - 1.0/3.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(2)->GetValue(NODAL_AREA)) - 1.0/3.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(5)->GetValue(NODAL_AREA)) - 1.0/3.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(6)->GetValue(NODAL_AREA)) - 1.0/6.0, tolerance);
}

/**
* Checks the correct work of the nodal gradient compute
* Test tetrahedra
*/
KRATOS_TEST_CASE_IN_SUITE(NodalArea2, KratosCoreFastSuite)
{
    Model current_model;

    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.SetBufferSize(2);

    this_model_part.AddNodalSolutionStepVariable(DISTANCE);
    this_model_part.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);
    this_model_part.AddNodalSolutionStepVariable(NODAL_AREA);

    auto& process_info = this_model_part.GetProcessInfo();
    process_info[STEP] = 1;
    process_info[NL_ITERATION_NUMBER] = 1;

    CppTestsUtilities::Create3DGeometry(this_model_part);

    // Calculate nodal area
    using HistNodalAreaProcess = CalculateNodalAreaProcess<CalculateNodalAreaSettings::SaveAsHistoricalVariable>; ;
    HistNodalAreaProcess hist_process(this_model_part, 3);
    hist_process.Execute();

    // // DEBUG
    // GiDIODebugNodalArea(this_model_part);

    const double tolerance = 1.0e-8;
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(1)->FastGetSolutionStepValue(NODAL_AREA)) - 1.0/12.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(2)->FastGetSolutionStepValue(NODAL_AREA)) - 1.0/12.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(3)->FastGetSolutionStepValue(NODAL_AREA)) - 1.0/4.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(5)->FastGetSolutionStepValue(NODAL_AREA)) - 1.0/12.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(9)->FastGetSolutionStepValue(NODAL_AREA)) - 1.0/4.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(10)->FastGetSolutionStepValue(NODAL_AREA)) - 1.0/12.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(11)->FastGetSolutionStepValue(NODAL_AREA)) - 1.0/12.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(12)->FastGetSolutionStepValue(NODAL_AREA)) - 1.0/12.0, tolerance);

    using NonHistNodalAreaProcess= CalculateNodalAreaProcess<CalculateNodalAreaSettings::SaveAsNonHistoricalVariable>; ;
    NonHistNodalAreaProcess nonhist_process(this_model_part, 3);
    nonhist_process.Execute();

    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(1)->GetValue(NODAL_AREA)) - 1.0/12.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(2)->GetValue(NODAL_AREA)) - 1.0/12.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(3)->GetValue(NODAL_AREA)) - 1.0/4.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(5)->GetValue(NODAL_AREA)) - 1.0/12.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(9)->GetValue(NODAL_AREA)) - 1.0/4.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(10)->GetValue(NODAL_AREA)) - 1.0/12.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(11)->GetValue(NODAL_AREA)) - 1.0/12.0, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(12)->GetValue(NODAL_AREA)) - 1.0/12.0, tolerance);
}

/**
* Checks the correct work of the nodal gradient compute
* Test quadrilateral
*/
KRATOS_TEST_CASE_IN_SUITE(NodalArea3, KratosCoreFastSuite)
{
    Model current_model;

    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.SetBufferSize(2);

    this_model_part.AddNodalSolutionStepVariable(NODAL_AREA);

    auto& process_info = this_model_part.GetProcessInfo();
    process_info[STEP] = 1;
    process_info[NL_ITERATION_NUMBER] = 1;

    CppTestsUtilities::Create2DQuadrilateralsGeometry(this_model_part);

    using HistNodalAreaProcess = CalculateNodalAreaProcess<CalculateNodalAreaSettings::SaveAsHistoricalVariable>; ;
    HistNodalAreaProcess hist_process(this_model_part, 2);
    hist_process.Execute();

    // // DEBUG
    // GiDIODebugNodalArea(this_model_part);

    const double tolerance = 1.0e-8;
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(1)->FastGetSolutionStepValue(NODAL_AREA)) - 0.25, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(2)->FastGetSolutionStepValue(NODAL_AREA)) - 0.5, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(3)->FastGetSolutionStepValue(NODAL_AREA)) - 0.5, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(4)->FastGetSolutionStepValue(NODAL_AREA)) - 0.25, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(5)->FastGetSolutionStepValue(NODAL_AREA)) - 0.25, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(6)->FastGetSolutionStepValue(NODAL_AREA)) - 0.25, tolerance);

    using NonHistNodalAreaProcess= CalculateNodalAreaProcess<CalculateNodalAreaSettings::SaveAsNonHistoricalVariable>; ;
    NonHistNodalAreaProcess nonhist_process(this_model_part, 2);
    nonhist_process.Execute();

    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(1)->GetValue(NODAL_AREA)) - 0.25, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(2)->GetValue(NODAL_AREA)) - 0.5, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(3)->GetValue(NODAL_AREA)) - 0.5, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(4)->GetValue(NODAL_AREA)) - 0.25, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(5)->GetValue(NODAL_AREA)) - 0.25, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(6)->GetValue(NODAL_AREA)) - 0.25, tolerance);
}

/**
* Checks the correct work of the nodal gradient compute
* Test tetrahedra
*/
KRATOS_TEST_CASE_IN_SUITE(NodalArea4, KratosCoreFastSuite)
{
    Model current_model;

    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.SetBufferSize(2);

    this_model_part.AddNodalSolutionStepVariable(DISTANCE);
    this_model_part.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);
    this_model_part.AddNodalSolutionStepVariable(NODAL_AREA);

    auto& process_info = this_model_part.GetProcessInfo();
    process_info[STEP] = 1;
    process_info[NL_ITERATION_NUMBER] = 1;

    CppTestsUtilities::Create3DHexahedraGeometry(this_model_part);

    // Calculate nodal area
    using HistNodalAreaProcess = CalculateNodalAreaProcess<CalculateNodalAreaSettings::SaveAsHistoricalVariable>; ;
    HistNodalAreaProcess hist_process(this_model_part, 3);
    hist_process.Execute();

    // // DEBUG
    // GiDIODebugNodalArea(this_model_part);

    const double tolerance = 1.0e-8;
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(1)->FastGetSolutionStepValue(NODAL_AREA)) - 0.125, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(2)->FastGetSolutionStepValue(NODAL_AREA)) - 0.125, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(3)->FastGetSolutionStepValue(NODAL_AREA)) - 0.125, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(4)->FastGetSolutionStepValue(NODAL_AREA)) - 0.25, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(5)->FastGetSolutionStepValue(NODAL_AREA)) - 0.125, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(6)->FastGetSolutionStepValue(NODAL_AREA)) - 0.25, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(7)->FastGetSolutionStepValue(NODAL_AREA)) - 0.25, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(8)->FastGetSolutionStepValue(NODAL_AREA)) - 0.25, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(9)->FastGetSolutionStepValue(NODAL_AREA)) - 0.125, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(10)->FastGetSolutionStepValue(NODAL_AREA)) - 0.125, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(11)->FastGetSolutionStepValue(NODAL_AREA)) - 0125, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(12)->FastGetSolutionStepValue(NODAL_AREA)) - 0.125, tolerance);

    using NonHistNodalAreaProcess= CalculateNodalAreaProcess<CalculateNodalAreaSettings::SaveAsNonHistoricalVariable>; ;
    NonHistNodalAreaProcess nonhist_process(this_model_part, 3);
    nonhist_process.Execute();

    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(1)->GetValue(NODAL_AREA)) - 0.125, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(2)->GetValue(NODAL_AREA)) - 0.125, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(3)->GetValue(NODAL_AREA)) - 0.125, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(4)->GetValue(NODAL_AREA)) - 0.25, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(5)->GetValue(NODAL_AREA)) - 0.125, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(6)->GetValue(NODAL_AREA)) - 0.25, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(7)->GetValue(NODAL_AREA)) - 0.25, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(8)->GetValue(NODAL_AREA)) - 0.25, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(9)->GetValue(NODAL_AREA)) - 0.125, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(10)->GetValue(NODAL_AREA)) - 0.125, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(11)->GetValue(NODAL_AREA)) - 0125, tolerance);
    KRATOS_EXPECT_LE(std::abs(this_model_part.pGetNode(12)->GetValue(NODAL_AREA)) - 0.125, tolerance);
}
} // namespace Kratos::Testing
