//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//                   Ashish Darekar
//

// System includes
#include <type_traits>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"
#include "utilities/atomic_utilities.h"
#include "conversion_utilities.h"

namespace Kratos {

template<class TDataType>
void ConversionUtilities::ConvertElementalDataToNodalData(
    ModelPart& rModelPart,
    const Variable<TDataType>& rElementalVariable,
    const Variable<TDataType>& rNodalVariable )
{
    // prepare nodal variable
    VariableUtils().SetHistoricalVariableToZero(rNodalVariable, rModelPart.Nodes());

    block_for_each(rModelPart.Elements(), [&](Element& rElement){
        const auto& elem_rVariable =  rElement.GetValue(rElementalVariable);

        const std::size_t num_nodes = rElement.GetGeometry().PointsNumber();

        for (auto& r_node : rElement.GetGeometry().Points()){
            if constexpr(std::is_same_v<TDataType, double>) {
                AtomicAdd( r_node.FastGetSolutionStepValue(rNodalVariable), (elem_rVariable / static_cast<double>(num_nodes)) );
            } else if constexpr(std::is_same_v<TDataType, array_1d<double, 3>>) {
                AtomicAddVector( r_node.FastGetSolutionStepValue(rNodalVariable), (elem_rVariable / static_cast<double>(num_nodes)) );
            } else {
                static_assert(!std::is_same_v<TDataType, TDataType>, "Unsupported data type.");
            }
        }
    });

    rModelPart.GetCommunicator().AssembleCurrentData(rNodalVariable);
}

template<class TDataType>
void ConversionUtilities::ConvertNodalDataToElementalData(
    ModelPart& rModelPart,
    const Variable<TDataType>& rElementalVariable,
    const Variable<TDataType>& rNodalVariable )
{
    // prepare nodal variable
    VariableUtils().SetHistoricalVariableToZero(rNodalVariable, rModelPart.Nodes());

    block_for_each(rModelPart.Elements(), [&](Element& rElement){
        const std::size_t num_nodes = rElement.GetGeometry().PointsNumber();

        if constexpr(std::is_same_v<TDataType, double>) {
                double temp = 0.0;
                for (auto& r_node : rElement.GetGeometry().Points()){
                    temp += r_node.GetValue(rNodalVariable) / num_nodes;
                } 
                rElement.SetValue(rElementalVariable, temp);
            }
        else if constexpr(std::is_same_v<TDataType, array_1d<double, 3>>) {
            array_1d<double, 3> temp = ZeroVector(3);
            for (auto& r_node : rElement.GetGeometry().Points()){
                temp += r_node.GetValue(rNodalVariable) / num_nodes;
            }
            rElement.SetValue(rElementalVariable, temp);
        }
        
        else {
                static_assert(!std::is_same_v<TDataType, TDataType>, "Unsupported data type.");
            }
    });

    rModelPart.GetCommunicator().AssembleCurrentData(rElementalVariable);
}

// template instantiations
template void ConversionUtilities::ConvertElementalDataToNodalData<double>(ModelPart&, const Variable<double>&,  const Variable<double>&);
template void ConversionUtilities::ConvertElementalDataToNodalData<array_1d<double, 3>>(ModelPart&, const Variable<array_1d<double, 3>>&,  const Variable<array_1d<double, 3>>&);

template void ConversionUtilities::ConvertNodalDataToElementalData<double>(ModelPart&, const Variable<double>&,  const Variable<double>&);
template void ConversionUtilities::ConvertNodalDataToElementalData<array_1d<double, 3>>(ModelPart&, const Variable<array_1d<double, 3>>&,  const Variable<array_1d<double, 3>>&);


}  // namespace Kratos.