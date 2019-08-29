//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "custom_processes/auxiliary_processes/rans_logarithmic_y_plus_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_logarithmic_y_plus_velocity_sensitivities_process.h"
#include "includes/model_part.h"
#include "rans_modelling_application_variables.h"
#include "testing/testing.h"

// Application includes

namespace Kratos
{
namespace Testing
{
typedef ModelPart::NodeType NodeType;

typedef ModelPart::ElementType ElementType;

typedef Geometry<NodeType> GeometryType;

typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

/**
 * Auxiliar function to generate a triangular element to be tested.
 */
void GenerateTestModelPart(ModelPart& rModelPart)
{
    // Set buffer size
    rModelPart.SetBufferSize(2);

    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(DISTANCE);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(RANS_Y_PLUS);

    // Process info creation
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);

    // Set the element properties
    Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);
    p_elem_prop->SetValue(KINEMATIC_VISCOSITY, 3.0e-02);

    // Element creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elem_nodes{1, 2, 3};
    rModelPart.CreateNewElement("RANSEVMK2D3N", 1, elem_nodes, p_elem_prop);

    // Set the VELOCITY and PRESSURE nodal values
    array_1d<double, 3> v_1 = ZeroVector(3);
    array_1d<double, 3> v_2 = ZeroVector(3);
    array_1d<double, 3> v_3 = ZeroVector(3);
    v_1[0] = 10.0;
    v_1[1] = 20.0;
    v_2[0] = 1000.0;
    v_2[1] = 500.0;
    v_3[0] = 30.0;
    v_3[1] = 20.0;
    (rModelPart.GetNode(1)).GetSolutionStepValue(VELOCITY) = v_1;
    (rModelPart.GetNode(2)).GetSolutionStepValue(VELOCITY) = v_2;
    (rModelPart.GetNode(3)).GetSolutionStepValue(VELOCITY) = v_3;

    (rModelPart.GetNode(1)).GetSolutionStepValue(DISTANCE) = 0.01;
    (rModelPart.GetNode(2)).GetSolutionStepValue(DISTANCE) = 0.02;
    (rModelPart.GetNode(3)).GetSolutionStepValue(DISTANCE) = 0.05;

    // Set the DENSITY and DYNAMIC_VISCOSITY nodal values
    for (ModelPart::NodeIterator it_node = rModelPart.NodesBegin();
         it_node < rModelPart.NodesEnd(); ++it_node)
    {
        it_node->FastGetSolutionStepValue(KINEMATIC_VISCOSITY) =
            p_elem_prop->GetValue(KINEMATIC_VISCOSITY);
    }
}

/**
 * Checks the RANSYPlusModelLogarithmic process.
 */
KRATOS_TEST_CASE_IN_SUITE(RansLogarithmicYPlusVelocitySensitivitiesProcess, RANSYPlusModelSensitivities)
{
    Model model;
    ModelPart& r_model_part =
        model.CreateModelPart("RANSYPlusModelLogarithmic");
    GenerateTestModelPart(r_model_part);

    Parameters empty_parameters = Parameters(R"({
        "model_part_name" : "RANSYPlusModelLogarithmic"
    })");

    RansLogarithmicYPlusVelocitySensitivitiesProcess adjoint_process(model, empty_parameters);
    RansLogarithmicYPlusCalculationProcess primal_process(model, empty_parameters);

    auto& r_element = *r_model_part.ElementsBegin();
    auto& r_geometry = r_element.GetGeometry();

    // Calculate finite difference values

    // Calculate initial y_plus values
    primal_process.Check();
    primal_process.Execute();

    // Calculate adjoint values
    adjoint_process.Check();
    adjoint_process.Execute();
    Matrix& r_adjoint_values = r_element.GetValue(RANS_Y_PLUS_VELOCITY_DERIVATIVES);

    const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];
    const int number_of_nodes = r_model_part.NumberOfNodes();

    const double epsilon = 1e-8;
    primal_process.Execute();
    std::vector<double> y_plus_0(number_of_nodes);

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        y_plus_0[i_node] = r_geometry[i_node].FastGetSolutionStepValue(RANS_Y_PLUS);

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        for (int i_dim = 0; i_dim < domain_size; ++i_dim)
        {
            array_1d<double, 3>& r_velocity =
                r_geometry[i_node].FastGetSolutionStepValue(VELOCITY);
            r_velocity[i_dim] += epsilon;
            primal_process.Execute();
            const double y_plus = r_geometry[i_node].FastGetSolutionStepValue(RANS_Y_PLUS);
            const double y_plus_sensitivity = (y_plus - y_plus_0[i_node]) / epsilon;
            r_velocity[i_dim] -= epsilon;

            KRATOS_CHECK_NEAR(r_adjoint_values(i_node, i_dim), y_plus_sensitivity, 1e-6);
            KRATOS_CHECK_NOT_EQUAL(y_plus_sensitivity, 0.0);
        }
    }
}
} // namespace Testing
} // namespace Kratos.
