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
//                   Ihar Antonau

// System includes
#include <type_traits>

// External includes
#include <map>

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"
#include "utilities/atomic_utilities.h"
#include "conversion_utilities.h"

namespace Kratos {

template<class TDataType>
void ConversionUtilities::ConvertElementalDataToNodalDataTranspose(
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
void ConversionUtilities::ConvertElementalDataToNodalDataDirect(
    ModelPart& rModelPart,
    const Variable<TDataType>& rElementalVariable,
    const Variable<TDataType>& rNodalVariable )
{
    // prepare nodal variable
    VariableUtils().SetHistoricalVariableToZero(rNodalVariable, rModelPart.Nodes());
    std::map<int, int> node_element_count;

    block_for_each(rModelPart.Elements(), [&](Element& rElement){
        const auto& elem_rVariable =  rElement.GetValue(rElementalVariable);

        for (auto& r_node : rElement.GetGeometry().Points()){
            if constexpr(std::is_same_v<TDataType, double>) {
                AtomicAdd( r_node.FastGetSolutionStepValue(rNodalVariable), (elem_rVariable) );
                node_element_count[r_node.Id()] += 1;

            } else if constexpr(std::is_same_v<TDataType, array_1d<double, 3>>) {
                AtomicAddVector( r_node.FastGetSolutionStepValue(rNodalVariable), (elem_rVariable) );
                node_element_count[r_node.Id()] += 1;
            } else {
                static_assert(!std::is_same_v<TDataType, TDataType>, "Unsupported data type.");
            }
        }
    });

    block_for_each(rModelPart.Nodes(), [&](Node& rNode){
        if (node_element_count.find(rNode.Id()) != node_element_count.end()){
            rNode.FastGetSolutionStepValue(rNodalVariable) /= node_element_count[rNode.Id()];
        }
        else{
            KRATOS_ERROR << "Node " << rNode.Id() << " has no associated elements." << std::endl;
        }
    });

    rModelPart.GetCommunicator().AssembleCurrentData(rNodalVariable);
}

template<class TDataType>
void ConversionUtilities::ConvertNodalDataToElementalDataDirect(
    ModelPart& rModelPart,
    const Variable<TDataType>& rElementalVariable,
    const Variable<TDataType>& rNodalVariable )
{
    // prepare elemental variable
    VariableUtils().SetNonHistoricalVariableToZero(rElementalVariable, rModelPart.Elements());

    block_for_each(rModelPart.Elements(), [&](Element& rElement){
        const std::size_t num_nodes = rElement.GetGeometry().PointsNumber();

        if constexpr(std::is_same_v<TDataType, double>) {
                double temp = 0.0;
                for (auto& r_node : rElement.GetGeometry().Points()){
                    temp += r_node.FastGetSolutionStepValue(rNodalVariable) / num_nodes;
                } 
                rElement.SetValue(rElementalVariable, temp);
            }
        else if constexpr(std::is_same_v<TDataType, array_1d<double, 3>>) {
            array_1d<double, 3> temp = ZeroVector(3);
            for (auto& r_node : rElement.GetGeometry().Points()){
                temp += r_node.FastGetSolutionStepValue(rNodalVariable) / num_nodes;
            }
        
        else {
                static_assert(!std::is_same_v<TDataType, TDataType>, "Unsupported data type.");
            }
    });

    rModelPart.GetCommunicator().AssembleCurrentData(rElementalVariable);
}

template<class TDataType>
void ConversionUtilities::ConvertNodalDataToElementalDataTranspose(
    ModelPart& rModelPart,
    const Variable<TDataType>& rElementalVariable,
    const Variable<TDataType>& rNodalVariable )
{
    // prepare elemental variable
    VariableUtils().SetNonHistoricalVariableToZero(rElementalVariable, rModelPart.Elements());

    // count number of shared elements for each node
    std::map<int, int> node_element_count;
    block_for_each(rModelPart.Elements(), [&](Element& rElement){
        for (auto& r_node : rElement.GetGeometry().Points()){
            node_element_count[r_node.Id()] += 1;
        }
    });

    // sum nodal value contributions to each element
    block_for_each(rModelPart.Elements(), [&](Element& rElement){

        if constexpr(std::is_same_v<TDataType, double>) {
                double temp = 0.0;
                for (auto& r_node : rElement.GetGeometry().Points()){
                    temp += r_node.FastGetSolutionStepValue(rNodalVariable) / node_element_count[r_node.Id()];
                } 
                rElement.SetValue(rElementalVariable, temp);
            }
        else if constexpr(std::is_same_v<TDataType, array_1d<double, 3>>) {
            array_1d<double, 3> temp = ZeroVector(3);
            for (auto& r_node : rElement.GetGeometry().Points()){
                temp += r_node.FastGetSolutionStepValue(rNodalVariable) / node_element_count[r_node.Id()];
            }
            rElement.SetValue(rElementalVariable, temp);
        
        else {
                static_assert(!std::is_same_v<TDataType, TDataType>, "Unsupported data type.");
            }
    });

    rModelPart.GetCommunicator().AssembleCurrentData(rElementalVariable);
}

// template instantiations
template void ConversionUtilities::ConvertElementalDataToNodalDataTranspose<double>(ModelPart&, const Variable<double>&,  const Variable<double>&);
template void ConversionUtilities::ConvertElementalDataToNodalDataTranspose<array_1d<double, 3>>(ModelPart&, const Variable<array_1d<double, 3>>&,  const Variable<array_1d<double, 3>>&);

template void ConversionUtilities::ConvertElementalDataToNodalDataDirect<double>(ModelPart&, const Variable<double>&,  const Variable<double>&);
template void ConversionUtilities::ConvertElementalDataToNodalDataDirect<array_1d<double, 3>>(ModelPart&, const Variable<array_1d<double, 3>>&,  const Variable<array_1d<double, 3>>&);

template void ConversionUtilities::ConvertNodalDataToElementalDataDirect<double>(ModelPart&, const Variable<double>&,  const Variable<double>&);
template void ConversionUtilities::ConvertNodalDataToElementalDataDirect<array_1d<double, 3>>(ModelPart&, const Variable<array_1d<double, 3>>&,  const Variable<array_1d<double, 3>>&);

template void ConversionUtilities::ConvertNodalDataToElementalDataTranspose<double>(ModelPart&, const Variable<double>&,  const Variable<double>&);
template void ConversionUtilities::ConvertNodalDataToElementalDataTranspose<array_1d<double, 3>>(ModelPart&, const Variable<array_1d<double, 3>>&,  const Variable<array_1d<double, 3>>&);


}  // namespace Kratos.
