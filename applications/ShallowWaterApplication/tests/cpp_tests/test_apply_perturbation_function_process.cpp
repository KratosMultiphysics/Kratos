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
#include "includes/properties.h"
#include "includes/kratos_parameters.h"
#include "shallow_water_application_variables.h"
#include "custom_processes/apply_perturbation_function_process.h"

namespace Kratos {

namespace Testing {

typedef ModelPart::IndexType               IndexType;
typedef ModelPart::NodeIterator     NodeIteratorType;

double periodic_function(double& x)
{
    double distance = 1.0;
    double perturbation = 1.0;
    double half_wave_number = std::acos(-1) / distance;
    return 0.5 * perturbation * (1 + std::cos(half_wave_number * x));
}

KRATOS_TEST_CASE_IN_SUITE(ApplyPerturbationFunctionProcess, ShallowWaterApplicationFastSuite)
{
    Model model;
    ModelPart& model_part = model.CreateModelPart("main", 2);
    ModelPart& sub_model_part = model_part.CreateSubModelPart("sub");

    // Variables addition
    model_part.AddNodalSolutionStepVariable(HEIGHT);
    model_part.AddNodalSolutionStepVariable(FREE_SURFACE_ELEVATION);
    model_part.AddNodalSolutionStepVariable(WATER_SURFACE);
    model_part.AddNodalSolutionStepVariable(BATHYMETRY);

    // Process info creation
    const double gravity = 9.81;
    model_part.GetProcessInfo().SetValue(GRAVITY_Z, gravity);

    // Geometry creation
    model_part.CreateNewNode( 1, 0.0, 0.0, 0.0);
    model_part.CreateNewNode( 2, 0.3, 0.0, 0.0);
    model_part.CreateNewNode( 3, 1.1, 0.1, 0.0);
    model_part.CreateNewNode( 4, 1.6, 0.2, 0.0);
    model_part.CreateNewNode( 5, 2.0, 0.0, 0.0);
    model_part.CreateNewNode( 6, 0.0, 1.0, 0.0);
    model_part.CreateNewNode( 7, 0.8, 0.9, 0.0);
    model_part.CreateNewNode( 8, 1.2, 1.0, 0.0);
    model_part.CreateNewNode( 9, 1.7, 0.9, 0.0);
    model_part.CreateNewNode(10, 2.0, 1.0, 0.0);
    model_part.CreateNewNode(11, 2.0, 0.85, 0.0);
    model_part.CreateNewNode(12, 2.0, 0.15, 0.0);

    sub_model_part.AddNodes({5, 10, 11, 12});

    // Creation and execution of the process for scalar
    Parameters parameters = Parameters(R"({})");
    ApplyPerturbationFunctionProcess<Variable<double>>(
        model_part,
        sub_model_part.Nodes(),
        FREE_SURFACE_ELEVATION,
        parameters)();

    double tolerance = 0.01; // The tolerance depends on the nodes positions
    for (auto& node : model_part.Nodes())
    {
        double x = node.X();
        double exact_value = ((x < 1.0) ? (0.0) : (periodic_function(x)));
        double value = node.FastGetSolutionStepValue(FREE_SURFACE_ELEVATION);
        KRATOS_CHECK_NEAR(value, exact_value, tolerance);
    }

    // Creation and execution of the process for component
    ApplyPerturbationFunctionProcess<Variable<double>>(
        model_part,
        sub_model_part.Nodes(),
        WATER_SURFACE_Z,
        parameters)();

    for (auto& node : model_part.Nodes())
    {
        double x = node.X();
        double exact_value = ((x < 1.0) ? (0.0) : (periodic_function(x)));
        double value = node.FastGetSolutionStepValue(WATER_SURFACE_Z);
        KRATOS_CHECK_NEAR(value, exact_value, tolerance);
    }
}

} // namespace Testing

} // namespace Kratos
