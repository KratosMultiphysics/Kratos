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
#include <iomanip> // For std::boolalpha

// External includes

// Project includes
#include "containers/model.h"
#include "utilities/parallel_utilities.h"
#include "processes/calculate_nodal_distance_to_skin_process.h"
#include "spatial_containers/geometrical_objects_bins.h"

namespace Kratos
{

CalculateNodalDistanceToSkinProcess::CalculateNodalDistanceToSkinProcess(
    ModelPart& rVolumeModelPart,
    ModelPart& rSkinModelPart,
    const bool HistoricalVariable,
    const std::string& rDistanceVariableName
    ) : mrVolumeModelPart(rVolumeModelPart),
        mrSkinModelPart(rSkinModelPart),
        mHistoricalVariable(HistoricalVariable)
{
    // Assign distance variable
    if (rDistanceVariableName != "") {
        mpDistanceVariable = &KratosComponents<Variable<double>>::Get(rDistanceVariableName);
    }

    // Check it is serial
    KRATOS_ERROR_IF(mrVolumeModelPart.IsDistributed()) << "Distributed computation still not supported. Please update implementation as soon as MPI search is merged. See https://github.com/KratosMultiphysics/Kratos/pull/11719" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

CalculateNodalDistanceToSkinProcess::CalculateNodalDistanceToSkinProcess(
    Model& rModel,
    Parameters ThisParameters
    ) : mrVolumeModelPart(rModel.GetModelPart(ThisParameters["volume_model_part"].GetString())),
        mrSkinModelPart(rModel.GetModelPart(ThisParameters["skin_model_part"].GetString()))
{
    // Validate and assign defaults
    ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    // If historical or non-historical variables
    const std::string database = ThisParameters["distance_database"].GetString();
    if (database == "nodal_historical") {
        mHistoricalVariable = true;
    } else if (database == "nodal_non_historical") {
        mHistoricalVariable = false;
    } else {
        KRATOS_ERROR << "Only options are nodal_historical and nodal_non_historical, provided is: " << database << std::endl;
    }

    // Assign distance variable
    mpDistanceVariable = &KratosComponents<Variable<double>>::Get(ThisParameters["distance_variable"].GetString());

    // Check it is serial
    KRATOS_ERROR_IF(mrVolumeModelPart.IsDistributed()) << "Distributed computation still not supported. Please update implementation as soon as MPI search is merged. See https://github.com/KratosMultiphysics/Kratos/pull/11719" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void CalculateNodalDistanceToSkinProcess::Execute()
{
    // Define distance lambda
    const std::function<void(ModelPart::NodesContainerType& rNodes, GeometricalObjectsBins& rBins)> distance_lambda_historical = [this](ModelPart::NodesContainerType& rNodes, GeometricalObjectsBins& rBins) {
        const auto& r_variable = *mpDistanceVariable;
        block_for_each(rNodes, [&rBins, &r_variable](Node& rNode) {
            auto result = rBins.SearchNearest(rNode);
            rNode.FastGetSolutionStepValue(r_variable) = result.GetDistance();
        });
    };
    const std::function<void(ModelPart::NodesContainerType& rNodes, GeometricalObjectsBins& rBins)> distance_lambda_non_historical = [this](ModelPart::NodesContainerType& rNodes, GeometricalObjectsBins& rBins) {
        const auto& r_variable = *mpDistanceVariable;
        block_for_each(rNodes, [&rBins, &r_variable](Node& rNode) {
            auto result = rBins.SearchNearest(rNode);
            rNode.SetValue(r_variable, result.GetDistance());
        });
    };

    const auto* p_distance_lambda = mHistoricalVariable ? &distance_lambda_historical : &distance_lambda_non_historical;

    // First try to detect if elements or conditions
    const std::size_t number_of_elements_skin = mrSkinModelPart.NumberOfElements();
    const std::size_t number_of_conditions_skin = mrSkinModelPart.NumberOfConditions();
    if (number_of_elements_skin > 0) {
        KRATOS_WARNING_IF("CalculateNodalDistanceToSkinProcess", number_of_conditions_skin > 0) << "Skin model part has elements and conditions. Considering elements. " << std::endl;
        GeometricalObjectsBins bins(mrSkinModelPart.ElementsBegin(), mrSkinModelPart.ElementsEnd());

        // Call lambda
        (*p_distance_lambda)(mrVolumeModelPart.Nodes(), bins);
    } else if (number_of_conditions_skin > 0) {
        GeometricalObjectsBins bins(mrSkinModelPart.ConditionsBegin(), mrSkinModelPart.ConditionsEnd());

        // Call lambda
        (*p_distance_lambda)(mrVolumeModelPart.Nodes(), bins);
    }
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters CalculateNodalDistanceToSkinProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"({
        "volume_model_part" : "",
        "skin_model_part"   : "",
        "distance_database" : "nodal_historical",
        "distance_variable" : "DISTANCE"
    })");
    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

std::string CalculateNodalDistanceToSkinProcess::Info() const
{
    return "CalculateNodalDistanceToSkinProcess";
}

/***********************************************************************************/
/***********************************************************************************/


void CalculateNodalDistanceToSkinProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "CalculateNodalDistanceToSkinProcess" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/


void CalculateNodalDistanceToSkinProcess::PrintData(std::ostream& rOStream) const
{
    rOStream << "Volume model part:\n" << mrVolumeModelPart << "\nSkin model part:\n" << mrSkinModelPart << "\nHistorical variable: " << std::boolalpha << mHistoricalVariable << "\nDistance variable: " << mpDistanceVariable->Name() << std::endl;
}

} // namespace Kratos.