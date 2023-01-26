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
//

// System includes
#include <limits>

// External includes

// Project includes
#include "processes/check_same_modelpart_using_skin_distance.h"
// #include "processes/calculate_distance_to_skin_process.h" // Continuous version
#include "processes/calculate_discontinuous_distance_to_skin_process.h" // Discontinuous version

namespace Kratos
{
template<std::size_t TDim>
void CheckSameModelPartUsingSkinDistance<TDim>::Execute()
{
    KRATOS_TRY



    KRATOS_CATCH("")
}


/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
const Parameters CheckSameModelPartUsingSkinDistance<TDim>::GetDefaultParameters() const
{
    KRATOS_TRY

    const Parameters default_parameters = Parameters(R"(
    {
        "model_part_1_name"   : "PLEASE_SPECIFY_MODEL_PART_1_NAME",
        "model_part_2_name"   : "PLEASE_SPECIFY_MODEL_PART_2_NAME",
        "tolerance"           : 1.0e-3,
        "continuous_distance" : false
    })" );

    return default_parameters;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template class Kratos::CheckSameModelPartUsingSkinDistance<2>;
template class Kratos::CheckSameModelPartUsingSkinDistance<3>;

} // namespace Kratos
