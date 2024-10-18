//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "processes/structured_mesh_generator_process.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "custom_utilities/derivatives_recovery_utility.h"
#include "shallow_water_application_variables.h"
#include "includes/gid_io.h"


namespace Kratos {

namespace Testing {

typedef ModelPart::NodeType NodeType;

void AddVariables(ModelPart& rModelPart)
{
    rModelPart.AddNodalSolutionStepVariable(DISTANCE);
    rModelPart.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY_LAPLACIAN);
    rModelPart.AddNodalSolutionStepVariable(FIRST_DERIVATIVE_WEIGHTS);
    rModelPart.AddNodalSolutionStepVariable(SECOND_DERIVATIVE_WEIGHTS);
}

void Generate2DMesh(
    Model& rModel,
    const int NumberOfDivisions)
{
    Node::Pointer p_point_1(new NodeType(1,  0.0,  0.0, 0.0));
    Node::Pointer p_point_2(new NodeType(2,  0.0, 10.0, 0.0));
    Node::Pointer p_point_3(new NodeType(3, 10.0, 10.0, 0.0));
    Node::Pointer p_point_4(new NodeType(4, 10.0,  0.0, 0.0));

    Quadrilateral2D4<NodeType> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions":  4,
        "element_name":         "Element2D3N",
        "condition_name":       "LineCondition"
    })");
    mesher_parameters["number_of_divisions"].SetInt(NumberOfDivisions);

    ModelPart& r_model_part = rModel.CreateModelPart("model_part");
    AddVariables(r_model_part);
    StructuredMeshGeneratorProcess(geometry, r_model_part, mesher_parameters).Execute();
}

void Generate3DMesh(
    Model& rModel,
    const int NumberOfDivisions)
{
    Node::Pointer p_point1(new NodeType(1,  0.0,  0.0,  0.0));
    Node::Pointer p_point2(new NodeType(2, 10.0,  0.0,  0.0));
    Node::Pointer p_point3(new NodeType(3, 10.0, 10.0,  0.0));
    Node::Pointer p_point4(new NodeType(4,  0.0, 10.0,  0.0));
    Node::Pointer p_point5(new NodeType(5,  0.0,  0.0, 10.0));
    Node::Pointer p_point6(new NodeType(6, 10.0,  0.0, 10.0));
    Node::Pointer p_point7(new NodeType(7, 10.0, 10.0, 10.0));
    Node::Pointer p_point8(new NodeType(8,  0.0, 10.0, 10.0));

    Hexahedra3D8<NodeType> geometry(p_point1, p_point2, p_point3, p_point4, p_point5, p_point6, p_point7, p_point8);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions":  4,
        "element_name":         "Element3D4N",
        "condition_name":       "SurfaceCondition"
    })");
    mesher_parameters["number_of_divisions"].SetInt(NumberOfDivisions);

    ModelPart& r_model_part = rModel.CreateModelPart("model_part");
    AddVariables(r_model_part);
    StructuredMeshGeneratorProcess(geometry, r_model_part, mesher_parameters).Execute();
}

void PrintDebugMesh(ModelPart& rModelPart)
{
    GidIO<> io("output_mesh",
        GiD_PostMode::GiD_PostAscii,
        MultiFileFlag::SingleFile,
        WriteDeformedMeshFlag::WriteUndeformed,
        WriteConditionsFlag::WriteElementsOnly
    );
    io.InitializeMesh(0.0);
    io.WriteMesh(rModelPart.GetMesh());
    io.FinalizeMesh();

    io.InitializeResults(0.0, rModelPart.GetMesh());
    io.WriteNodalResults(DISTANCE, rModelPart.Nodes(), 0.0, 0);
    io.WriteNodalResults(VELOCITY, rModelPart.Nodes(), 0.0, 0);
    io.WriteNodalResults(DISTANCE_GRADIENT, rModelPart.Nodes(), 0.0, 0);
    io.WriteNodalResults(VELOCITY_LAPLACIAN, rModelPart.Nodes(), 0.0, 0);
    io.FinalizeResults();
}

template<std::size_t TDim>
void GradientRecoveryTest(
    Model& rModel,
    const std::function<double(array_1d<double,3>)>& rPrimitiveFunction,
    const std::function<array_1d<double,3>(array_1d<double,3>)>& rDerivativeFunction,
    const double Tolerance)
{
    const bool print_debug_mesh = false;

    ModelPart& r_model_part = rModel.GetModelPart("model_part");
    block_for_each(r_model_part.Nodes(), [&](NodeType& rNode){
        rNode.FastGetSolutionStepValue(DISTANCE) = rPrimitiveFunction(rNode);
    });

    FindGlobalNodalNeighboursProcess(r_model_part).Execute();
    DerivativesRecoveryUtility<TDim>::CalculatePolynomialWeights(r_model_part);
    DerivativesRecoveryUtility<TDim>::RecoverGradient(r_model_part, DISTANCE, DISTANCE_GRADIENT);

    if (print_debug_mesh) {
        PrintDebugMesh(r_model_part);
    }

    block_for_each(r_model_part.Nodes(), [&](NodeType& rNode){
        array_1d<double,3> computed_gradient = rNode.FastGetSolutionStepValue(DISTANCE_GRADIENT);
        array_1d<double,3> exact_gradient = rDerivativeFunction(rNode);
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(computed_gradient, exact_gradient, Tolerance)
    });
}

template<std::size_t TDim>
void LaplacianRecoveryTest(
    Model& rModel,
    const std::function<array_1d<double,3>(array_1d<double,3>)>& rPrimitiveFunction,
    const std::function<array_1d<double,3>(array_1d<double,3>)>& rDerivativeFunction,
    const double Tolerance)
{
    const bool print_debug_mesh = false;

    ModelPart& r_model_part = rModel.GetModelPart("model_part");
    block_for_each(r_model_part.Nodes(), [&](NodeType& rNode){
        rNode.FastGetSolutionStepValue(VELOCITY) = rPrimitiveFunction(rNode);
    });

    FindGlobalNodalNeighboursProcess(r_model_part).Execute();
    DerivativesRecoveryUtility<TDim>::CalculatePolynomialWeights(r_model_part);
    DerivativesRecoveryUtility<TDim>::RecoverLaplacian(r_model_part, VELOCITY, VELOCITY_LAPLACIAN);

    if (print_debug_mesh) {
        PrintDebugMesh(r_model_part);
    }

    block_for_each(r_model_part.Nodes(), [&](NodeType& rNode){
        array_1d<double,3> computed_laplacian = rNode.FastGetSolutionStepValue(VELOCITY_LAPLACIAN);
        array_1d<double,3> exact_laplacian = rDerivativeFunction(rNode);
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(computed_laplacian, exact_laplacian, Tolerance)
    });
}

KRATOS_TEST_CASE_IN_SUITE(DerivativesRecoveryUtility2DGradientX, ShallowWaterApplicationFastSuite)
{
    Model model;
    const int number_of_divisions = 4;
    auto f = [] (array_1d<double,3> x) {return 0.1 * std::pow(x[0] - 4.0, 2);};
    auto df = [] (array_1d<double,3> x) {
        array_1d<double,3> derivative({0.2 * (x[0] - 4.0), 0.0, 0.0});
        return derivative;
    };
    const double tolerance = 1e-12;
    Generate2DMesh(model, number_of_divisions);
    GradientRecoveryTest<2>(model, f, df, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(DerivativesRecoveryUtility2DGradientXY, ShallowWaterApplicationFastSuite)
{
    Model model;
    const int number_of_divisions = 4;
    auto f = [] (array_1d<double,3> x) {return 0.1 * std::pow(x[0] - 6.0, 2) + 0.1 * std::pow(x[1] - 5.0, 2);};
    auto df = [] (array_1d<double,3> x) {
        array_1d<double,3> derivative({0.2 * (x[0] - 6.0), 0.2 * (x[1] - 5.0), 0.0});
        return derivative;
    };
    const double tolerance = 1e-12;
    Generate2DMesh(model, number_of_divisions);
    GradientRecoveryTest<2>(model, f, df, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(DerivativesRecoveryUtility2DGradientHyperbolicFunction, ShallowWaterApplicationFastSuite)
{
    // Note, the derivatives recovery cannot pass this test analytically, since it uses a polynomial approximation.
    Model model;
    const int number_of_divisions = 20;
    auto f = [] (array_1d<double,3> x) {
        const double k = 0.2;
        const double p = 1.0;
        return 0.5 / std::cosh(k*x[0] + p);};
    auto df = [] (array_1d<double,3> x) {
        const double k = 0.2;
        const double p = 1.0;
        array_1d<double,3> derivative({-0.5*k*std::tanh(k*x[0] + p)/std::cosh(k*x[0] + p), 0.0, 0.0});
        return derivative;
    };
    const double tolerance = 1e-2;
    Generate2DMesh(model, number_of_divisions);
    GradientRecoveryTest<2>(model, f, df, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(DerivativesRecoveryUtility3DGradientZ, ShallowWaterApplicationFastSuite)
{
    Model model;
    const int number_of_divisions = 4;
    auto f = [] (array_1d<double,3> x) {return 0.3 * std::pow(x[2] - 4.0, 2);};
    auto df = [] (array_1d<double,3> x) {
        array_1d<double,3> derivative({0.0, 0.0, 0.6 * (x[2] - 4.0)});
        return derivative;
    };
    const double tolerance = 1e-12;
    Generate3DMesh(model, number_of_divisions);
    GradientRecoveryTest<3>(model, f, df, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(DerivativesRecoveryUtility2DLaplacian, ShallowWaterApplicationFastSuite)
{
    Model model;
    const int number_of_divisions = 4;
    auto f = [] (array_1d<double,3> x) {
        array_1d<double,3> primitive ({0.1*std::pow(x[0]-6.0, 2), 0.2*std::pow(x[1]-5.0, 2), 0.0});
        return primitive;
    };
    auto df = [] (array_1d<double,3> x) {
        array_1d<double,3> derivative({0.2, 0.4, 0.0});
        return derivative;
    };
    const double tolerance = 1e-12;
    Generate2DMesh(model, number_of_divisions);
    LaplacianRecoveryTest<2>(model, f, df, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(DerivativesRecoveryUtility2DLaplacianXY, ShallowWaterApplicationFastSuite)
{
    Model model;
    const int number_of_divisions = 4;
    auto f = [] (array_1d<double,3> x) {
        array_1d<double,3> primitive ({0.1*x[0]*x[1], 0.0, 0.0});
        return primitive;
    };
    auto df = [] (array_1d<double,3> x) {
        array_1d<double,3> derivative({0.0, 0.1, 0.0});
        return derivative;
    };
    const double tolerance = 1e-12;
    Generate2DMesh(model, number_of_divisions);
    LaplacianRecoveryTest<2>(model, f, df, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(DerivativesRecoveryUtility3DLaplacian, ShallowWaterApplicationFastSuite)
{
    Model model;
    const int number_of_divisions = 4;
    auto f = [] (array_1d<double,3> x) {
        array_1d<double,3> primitive ({0.1*x[0]*x[2], 0.0, std::pow(x[2],2)});
        return primitive;
    };
    auto df = [] (array_1d<double,3> x) {
        array_1d<double,3> derivative({0.0, 0.0, 2.1});
        return derivative;
    };
    const double tolerance = 1e-12;
    Generate3DMesh(model, number_of_divisions);
    LaplacianRecoveryTest<3>(model, f, df, tolerance);
}

} // namespace Testing

} // namespace Kratos
