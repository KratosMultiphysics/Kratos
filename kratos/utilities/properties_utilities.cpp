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
#include "utilities/properties_utilities.h"

namespace Kratos
{
namespace PropertiesUtilities
{
void CopyPropertiesValues(
    const Properties& rOriginProperties,
    Properties& rDestinatiionProperties
    )
{
    KRATOS_TRY

    const auto& r_origin_data = rOriginProperties.GetData();
    auto& r_destination_data = rDestinatiionProperties.GetData();
    if (!r_destination_data.IsEmpty()) r_destination_data.Clear();

    // Copy data
    r_destination_data = DataValueContainer(r_origin_data);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

} // namespace PropertiesUtilities
} // namespace Kratos
