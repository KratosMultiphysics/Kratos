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

void Generate2DMesh(
    Model& rModel,
    const int NumberOfDivisions)
{
    Node<3>::Pointer p_point_1(new NodeType(1,  0.0,  0.0, 0.0));
    Node<3>::Pointer p_point_2(new NodeType(2,  0.0, 10.0, 0.0));
    Node<3>::Pointer p_point_3(new NodeType(3, 10.0, 10.0, 0.0));
    Node<3>::Pointer p_point_4(new NodeType(4, 10.0,  0.0, 0.0));

    Quadrilateral2D4<NodeType> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions":  4,
        "element_name":         "Element2D3N",
        "condition_name":       "LineCondition"
    })");
    mesher_parameters["number_of_divisions"].SetInt(NumberOfDivisions);

    ModelPart& r_model_part = rModel.CreateModelPart("model_part");
    r_model_part.AddNodalSolutionStepVariable(DISTANCE);
    r_model_part.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);
    r_model_part.AddNodalSolutionStepVariable(FIRST_DERIVATIVE_WEIGHTS);
    r_model_part.AddNodalSolutionStepVariable(SECOND_DERIVATIVE_WEIGHTS);
    StructuredMeshGeneratorProcess(geometry, r_model_part, mesher_parameters).Execute();
}

void Generate3DMesh(
    Model& rModel,
    const int NumberOfDivisions)
{
    Node<3>::Pointer p_point1(new NodeType(1,  0.0,  0.0,  0.0));
    Node<3>::Pointer p_point2(new NodeType(2, 10.0,  0.0,  0.0));
    Node<3>::Pointer p_point3(new NodeType(3, 10.0, 10.0,  0.0));
    Node<3>::Pointer p_point4(new NodeType(4,  0.0, 10.0,  0.0));
    Node<3>::Pointer p_point5(new NodeType(5,  0.0,  0.0, 10.0));
    Node<3>::Pointer p_point6(new NodeType(6, 10.0,  0.0, 10.0));
    Node<3>::Pointer p_point7(new NodeType(7, 10.0, 10.0, 10.0));
    Node<3>::Pointer p_point8(new NodeType(8,  0.0, 10.0, 10.0));

    Hexahedra3D8<NodeType> geometry(p_point1, p_point2, p_point3, p_point4, p_point5, p_point6, p_point7, p_point8);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions":  4,
        "element_name":         "Element3D4N",
        "condition_name":       "SurfaceCondition"
    })");
    mesher_parameters["number_of_divisions"].SetInt(NumberOfDivisions);

    ModelPart& r_model_part = rModel.CreateModelPart("model_part");
    r_model_part.AddNodalSolutionStepVariable(DISTANCE);
    r_model_part.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);
    r_model_part.AddNodalSolutionStepVariable(FIRST_DERIVATIVE_WEIGHTS);
    r_model_part.AddNodalSolutionStepVariable(SECOND_DERIVATIVE_WEIGHTS);
    StructuredMeshGeneratorProcess(geometry, r_model_part, mesher_parameters).Execute();
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
    DerivativesRecoveryUtility<TDim>::ExtendNeighborsPatch(r_model_part);
    DerivativesRecoveryUtility<TDim>::CalculatePolynomialWeights(r_model_part);
    DerivativesRecoveryUtility<TDim>::RecoverGradient(r_model_part, DISTANCE, DISTANCE_GRADIENT);

    if (print_debug_mesh) {
        GidIO<> io("output_mesh",
            GiD_PostMode::GiD_PostAscii,
            MultiFileFlag::SingleFile,
            WriteDeformedMeshFlag::WriteUndeformed,
            WriteConditionsFlag::WriteElementsOnly
        );
        io.InitializeMesh(0.0);
        io.WriteMesh(r_model_part.GetMesh());
        io.FinalizeMesh();

        io.InitializeResults(0.0, r_model_part.GetMesh());
        io.WriteNodalResults(DISTANCE, r_model_part.Nodes(), 0.0, 0);
        io.WriteNodalResults(DISTANCE_GRADIENT, r_model_part.Nodes(), 0.0, 0);
        io.FinalizeResults();
    }

    block_for_each(r_model_part.Nodes(), [&](NodeType& rNode){
        array_1d<double,3> computed_gradient = rNode.FastGetSolutionStepValue(DISTANCE_GRADIENT);
        array_1d<double,3> exact_gradient = rDerivativeFunction(rNode);
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(computed_gradient, exact_gradient, Tolerance)
    });
}

KRATOS_TEST_CASE_IN_SUITE(DerivativesRecoveryUtility2DParabolaX, ShallowWaterApplicationFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(DerivativesRecoveryUtility2DParabolaXY, ShallowWaterApplicationFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(DerivativesRecoveryUtility2DHyperbolicFunction, ShallowWaterApplicationFastSuite)
{
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

KRATOS_TEST_CASE_IN_SUITE(DerivativesRecoveryUtility3DParabolaZ, ShallowWaterApplicationFastSuite)
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

} // namespace Testing

} // namespace Kratos
