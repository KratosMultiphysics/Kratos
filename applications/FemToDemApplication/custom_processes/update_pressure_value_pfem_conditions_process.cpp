//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:             BSD License
//                               Kratos default license:
//kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#include "custom_processes/update_pressure_value_pfem_conditions_process.h"

namespace Kratos {

template <SizeType TDim>
UpdatePressureValuePfemConditionsProcess<TDim>::UpdatePressureValuePfemConditionsProcess(
    ModelPart& rModelPart)
    : mrModelPart(rModelPart)
{
}

/***********************************************************************************/
/***********************************************************************************/
template <SizeType TDim>
void UpdatePressureValuePfemConditionsProcess<TDim>::Execute()
{

}

/***********************************************************************************/
/***********************************************************************************/

template class UpdatePressureValuePfemConditionsProcess<2>;
template class UpdatePressureValuePfemConditionsProcess<3>;

}  // namespace Kratos