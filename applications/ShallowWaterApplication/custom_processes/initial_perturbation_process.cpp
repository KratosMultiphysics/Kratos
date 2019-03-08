//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


// System includes


// External includes


// Project includes
#include "includes/checks.h"
#include "initial_perturbation_process.h"


namespace Kratos
{

InitialPerturbationProcess::InitialPerturbationProcess(
    ModelPart& rThisModelPart,
    Parameters& rThisParameters
) : Process(), mrModelPart(rThisModelPart)
{
    // default parameters
    Parameters default_parameters = Parameters(R"(
    {
        "variable_name"              : "FREE_SURFACE_ELEVATION",
        "default_value"              : 0.0,
        "source_point_coordinates"   : [0.0, 0.0, 0.0],
        "distance_of_influence"      : 1.0,
        "maximum_perturbation_value" : 1.0
    })");
    rThisParameters.ValidateAndAssignDefaults(default_parameters);

    mVariable = KratosComponents< Variable<double> >::Get(rThisParameters["error_variable"].GetString());
    mDefaultValue = rThisParameters["default_value"].GetDouble();
    mInfluenceDistance = rThisParameters["distance_of_influence"].GetDouble();
    mPerturbation = rThisParameters["maximum_perturbation_value"].GetDouble();

    Check();
}

int InitialPerturbationProcess::Check()
{
    KRATOS_TRY
    if (mrModelPart.Nodes().size() != 0) {
        const auto& r_node = *mrModelPart.NodesBegin();
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mVariable, r_node);
    }
    KRATOS_CHECK(mInfluenceDistance <= std::numeric_limits<double>::epsilon());
    return 0;
    KRATOS_CATCH("")
}

void InitialPerturbationProcess::Execute()
{
    ExecuteBeforeSolutionLoop();
}

void InitialPerturbationProcess::ExecuteBeforeSolutionLoop()
{
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++)
    {
        auto i_node = mrModelPart.NodesBegin() + i;
        double& r_value = i_node->FastGetSolutionStepValue(mVariable);
        double value = 0;
        // Do something with the value
        r_value = value;
    }    
}

}  // namespace Kratos.


