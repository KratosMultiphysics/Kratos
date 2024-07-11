// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "structural_mechanics_fast_suite.h"
#include "containers/model.h"
#include "includes/kratos_flags.h"
// #include "includes/gid_io.h"
#include "tests/test_utilities/cpp_tests_utilities.h"

/* Processes */
#include "processes/integration_values_extrapolation_to_nodes_process.h"
#include "includes/mat_variables.h"

namespace Kratos::Testing
{
using GeometryType = Geometry<Node>;

// void GiDIODebugInternalExtrapolation(ModelPart& rModelPart, const std::string name = "")
// {
//     GidIO<> gid_io("TEST_INTERNAL_EXTRAPOLATION"+name, GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
//     const int nl_iter = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
//     const double label = static_cast<double>(nl_iter);

//     gid_io.InitializeMesh(label);
//     gid_io.WriteMesh(rModelPart.GetMesh());
//     gid_io.FinalizeMesh();
//     gid_io.InitializeResults(label, rModelPart.GetMesh());
//     auto r_variable = KratosComponents<Variable<double>>::Get("REFERENCE_DEFORMATION_GRADIENT_DETERMINANT");
//     gid_io.PrintOnGaussPoints(r_variable, rModelPart, label);
//     gid_io.WriteNodalResultsNonHistorical(r_variable, rModelPart.Nodes(), label);
//     gid_io.WriteNodalResultsNonHistorical(NODAL_AREA, rModelPart.Nodes(), label);
// }

void Create2DModelPartForExtrapolation(ModelPart& rModelPart)
{
    Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

    CppTestsUtilities::Create2DGeometry(rModelPart, "UpdatedLagrangianElement2D3N", false);

    ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw");
    auto p_this_law = r_clone_cl.Clone();
    p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);

    const auto& r_process_info = rModelPart.GetProcessInfo();

    // Initialize Elements
    for (auto& r_elem : rModelPart.Elements())
        r_elem.Initialize(r_process_info);
}

void Create3DModelPartForExtrapolation(ModelPart& rModelPart)
{
    Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

    CppTestsUtilities::Create3DGeometry(rModelPart, "UpdatedLagrangianElement3D4N", false);

    ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic3DLaw");
    auto p_this_law = r_clone_cl.Clone();
    p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);

    const auto& r_process_info = rModelPart.GetProcessInfo();

    // Initialize Elements
    for (auto& r_elem : rModelPart.Elements())
        r_elem.Initialize(r_process_info);
}

void CreateQuadratic3DModelPartForExtrapolation(ModelPart& rModelPart)
{
    Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

    CppTestsUtilities::Create3DQuadraticGeometry(rModelPart, "UpdatedLagrangianElement3D10N", false);

    ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic3DLaw");
    auto p_this_law = r_clone_cl.Clone();
    p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);

    const auto& r_process_info = rModelPart.GetProcessInfo();

    // Initialize Elements
    for (auto& r_elem : rModelPart.Elements())
        r_elem.Initialize(r_process_info);
}

/**
* Checks the correct work of the internal variable extrapolation process
* Test r_node
*/
KRATOS_TEST_CASE_IN_SUITE(TestIntegrationValuesExtrapolationToNodesProcessNode, KratosStructuralMechanicsFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
    ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
    r_current_process_info[DOMAIN_SIZE] = 2;

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);

    auto p_node = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    std::vector<std::size_t> nodes(1,1);
    auto p_element = r_model_part.CreateNewElement("NodalConcentratedElement3D1N", 1, nodes, p_elem_prop);

    const auto& r_process_info = r_model_part.GetProcessInfo();

    // Initialize Elements
    for (auto& r_elem : r_model_part.Elements())
        r_elem.Initialize(r_process_info);

    auto& process_info = r_model_part.GetProcessInfo();
    process_info[STEP] = 1;
    process_info[NL_ITERATION_NUMBER] = 1;

    // Compute extrapolation
    Parameters extrapolation_parameters = Parameters(R"(
    {
        "list_of_variables" : ["REFERENCE_DEFORMATION_GRADIENT_DETERMINANT"]
    })");
    IntegrationValuesExtrapolationToNodesProcess extrapolation_process(r_model_part, extrapolation_parameters);
    extrapolation_process.ExecuteBeforeSolutionLoop();
    extrapolation_process.ExecuteFinalizeSolutionStep();

    const auto& r_variable = KratosComponents<Variable<double>>::Get("REFERENCE_DEFORMATION_GRADIENT_DETERMINANT");
    for (auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_DOUBLE_EQ(r_node.GetValue(r_variable), 0.0);
    }

    extrapolation_process.ExecuteFinalize();
}

/**
* Checks the correct work of the internal variable extrapolation process
* Test triangle
*/
KRATOS_TEST_CASE_IN_SUITE(TestIntegrationValuesExtrapolationToNodesProcessTriangle, KratosStructuralMechanicsFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
    ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
    r_current_process_info[DOMAIN_SIZE] = 2;

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    Create2DModelPartForExtrapolation(r_model_part);

    auto& process_info = r_model_part.GetProcessInfo();
    process_info[STEP] = 1;
    process_info[NL_ITERATION_NUMBER] = 1;

    // Compute extrapolation
    Parameters extrapolation_parameters = Parameters(R"(
    {
        "list_of_variables"          : ["REFERENCE_DEFORMATION_GRADIENT_DETERMINANT"]
    })");
    IntegrationValuesExtrapolationToNodesProcess extrapolation_process(r_model_part, extrapolation_parameters);
    extrapolation_process.ExecuteBeforeSolutionLoop();
    extrapolation_process.ExecuteFinalizeSolutionStep();

    // // DEBUG
    // GiDIODebugInternalExtrapolation(r_model_part, "1");

    const double tolerance = 1.0e-8;
    const auto& r_variable = KratosComponents<Variable<double>>::Get("REFERENCE_DEFORMATION_GRADIENT_DETERMINANT");
    for (auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_LE(std::abs(r_node.GetValue(r_variable) - 1.0), tolerance);
    }

    extrapolation_process.ExecuteFinalize();
}

/**
* Checks the correct work of the internal variable extrapolation process
* Test tetrahedra
*/
KRATOS_TEST_CASE_IN_SUITE(TestIntegrationValuesExtrapolationToNodesProcessTetra, KratosStructuralMechanicsFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
    ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
    r_current_process_info[DOMAIN_SIZE] = 3;

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    Create3DModelPartForExtrapolation(r_model_part);

    auto& process_info = r_model_part.GetProcessInfo();
    process_info[STEP] = 1;
    process_info[NL_ITERATION_NUMBER] = 1;

    // Compute extrapolation
    Parameters extrapolation_parameters = Parameters(R"(
    {
        "list_of_variables"          : ["REFERENCE_DEFORMATION_GRADIENT_DETERMINANT"]
    })");
    IntegrationValuesExtrapolationToNodesProcess extrapolation_process(r_model_part, extrapolation_parameters);
    extrapolation_process.ExecuteBeforeSolutionLoop();
    extrapolation_process.ExecuteFinalizeSolutionStep();

    // // DEBUG
    // GiDIODebugInternalExtrapolation(r_model_part, "2");

    const double tolerance = 1.0e-8;
    const auto& r_variable = KratosComponents<Variable<double>>::Get("REFERENCE_DEFORMATION_GRADIENT_DETERMINANT");
    for (auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_LE(std::abs(r_node.GetValue(r_variable) - 1.0), tolerance);
    }

    extrapolation_process.ExecuteFinalize();
}

/**
* Checks the correct work of the internal variable extrapolation process
* Test quadratic tetrahedra
*/
KRATOS_TEST_CASE_IN_SUITE(TestIntegrationValuesExtrapolationToNodesProcessQuadTetra, KratosStructuralMechanicsFastSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
    ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
    r_current_process_info[DOMAIN_SIZE] = 3;

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    CreateQuadratic3DModelPartForExtrapolation(r_model_part);

    auto& process_info = r_model_part.GetProcessInfo();
    process_info[STEP] = 1;
    process_info[NL_ITERATION_NUMBER] = 1;

    // Compute extrapolation
    Parameters extrapolation_parameters = Parameters(R"(
    {
        "list_of_variables"          : ["REFERENCE_DEFORMATION_GRADIENT_DETERMINANT"]
    })");
    IntegrationValuesExtrapolationToNodesProcess extrapolation_process(r_model_part, extrapolation_parameters);
    extrapolation_process.ExecuteBeforeSolutionLoop();
    extrapolation_process.ExecuteFinalizeSolutionStep();

    // // DEBUG
    // GiDIODebugInternalExtrapolation(r_model_part, "3");

    const double tolerance = 1.0e-6;
    const auto& r_variable = KratosComponents<Variable<double>>::Get("REFERENCE_DEFORMATION_GRADIENT_DETERMINANT");
    for (auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_LE(std::abs(r_node.GetValue(r_variable) - 1.0), tolerance);
    }

    extrapolation_process.ExecuteFinalize();
}
}  // namespace Kratos::Testing.
