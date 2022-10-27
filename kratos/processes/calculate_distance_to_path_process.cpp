//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/geometrical_projection_utilities.h"
#include "processes/calculate_distance_to_path_process.h"

namespace Kratos
{

template<bool THistorical>
CalculateDistanceToPathProcess<THistorical>::CalculateDistanceToPathProcess(
    Model& rModel, 
    Parameters ThisParameters
    ) : mrModel(rModel),
        mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
void CalculateDistanceToPathProcess<THistorical>::Execute()
{
    /// TODO
    const std::string& r_distance_model_part_name = mThisParameters["distance_model_part_name"].GetString();
    auto& r_distance_model_part = mrModel.GetModelPart(r_distance_model_part_name);
    const std::string& r_path_model_part_name = mThisParameters["path_model_part_name"].GetString();
    auto& r_path_model_part = mrModel.GetModelPart(r_path_model_part_name);
    const bool brute_force_calculation = mThisParameters["brute_force_calculation"].GetBool();
    if (brute_force_calculation) {
        this->CalculateDistanceByBruteForce(r_distance_model_part);
    } else {
        this->CalculateDistance(r_distance_model_part);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
void CalculateDistanceToPathProcess<THistorical>::CalculateDistance(ModelPart& rModelPart)
{
    /// TODO
    const double radius_path = mThisParameters["radius_path"].GetDouble();
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
void CalculateDistanceToPathProcess<THistorical>::CalculateDistanceByBruteForce(ModelPart& rModelPart)
{
    /// TODO
    const double radius_path = mThisParameters["radius_path"].GetDouble();
}

/***********************************************************************************/
/***********************************************************************************/

template<bool THistorical>
const Parameters CalculateDistanceToPathProcess<THistorical>::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"({
        "distance_model_part_name" :  "",
        "path_model_part_name"     :  "",
        "brute_force_calculation"  : false,
        "radius_path"              : 0.0
    })" );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class CalculateDistanceToPathProcess<true>;
template class CalculateDistanceToPathProcess<false>;

} /// namespace Kratos