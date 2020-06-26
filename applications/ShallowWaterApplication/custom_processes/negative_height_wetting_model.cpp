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
#include "negative_height_wetting_model.h"
#include "shallow_water_application_variables.h"
#include "custom_utilities/shallow_water_utilities.h"


namespace Kratos
{

NegativeHeightWettingModel::NegativeHeightWettingModel(ModelPart& rModelPart, Parameters ThisParameters)
 : mrModelPart(rModelPart)
{
    Parameters default_parameters = Parameters(R"(
    {
        "model_name"     : "negative_height",
        "beta"           : 1e4
    })");
    ThisParameters.ValidateAndAssignDefaults(default_parameters);
    mBeta = ThisParameters["beta"].GetDouble();
    mDryHeight = mrModelPart.GetProcessInfo()[DRY_HEIGHT];
}

NegativeHeightWettingModel::NegativeHeightWettingModel(ModelPart& rModelPart, double Beta)
 : mrModelPart(rModelPart)
{
    mBeta = Beta;
}

void NegativeHeightWettingModel::ExecuteInitializeSolutionStep()
{
    // Definition of the first node iterator
    auto it_node_begin = mrModelPart.NodesBegin();

    // Execution of the implicit negative height method
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.NumberOfNodes()); ++i)
    {
        auto it_node = it_node_begin + i;
        const double height = it_node->FastGetSolutionStepValue(HEIGHT);
        const double manning = it_node->FastGetSolutionStepValue(MANNING);
        double& equivalent_manning = it_node->FastGetSolutionStepValue(EQUIVALENT_MANNING);
        double& porosity = it_node->FastGetSolutionStepValue(POROSITY);

        if (height > mDryHeight) {
            equivalent_manning = manning;
            porosity = 1.0;
        }
        else {
            equivalent_manning = manning * (1 - mBeta * (height - mDryHeight));
            porosity = 0.0;
        }
    }
}

void NegativeHeightWettingModel::ExecuteFinalizeSolutionStep()
{
    ShallowWaterUtilities().IdentifyWetDomain(mrModelPart, FLUID, 0.0);
}

}  // namespace Kratos
