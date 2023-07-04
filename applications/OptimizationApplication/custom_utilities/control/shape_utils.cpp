//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl
//                   Suneth Warnakulasuriya
//

// System includes
#include <algorithm>

// Project includes
#include "expression/literal_flat_expression.h"
#include "includes/define.h"
#include "utilities/parallel_utilities.h"

// Application includes

// Include base h
#include "shape_utils.h"

namespace Kratos
{

void ShapeUtils::GetNodalCoordinates(
            ContainerExpression<ModelPart::NodesContainerType>& rInputExpression)
{
    const IndexType local_size = rInputExpression.GetItemComponentCount();
    const auto& r_container = rInputExpression.GetContainer();
    const IndexType number_of_entities = rInputExpression.GetContainer().size();

    auto p_flat_data_expression = LiteralFlatExpression<double>::Create(number_of_entities, rInputExpression.GetItemShape());
    rInputExpression.SetExpression(p_flat_data_expression);
    auto& r_new_input_expression = *p_flat_data_expression;

    IndexPartition<IndexType>(r_container.size()).for_each([&](const IndexType Index){
        const auto p_node = (r_container.begin() + Index);
        const IndexType current_data_begin = Index * local_size;
        for (IndexType j = 0; j < local_size; ++j) {
            double& current_index_value = *(r_new_input_expression.begin() + current_data_begin + j);
            if (j == 0)
                current_index_value = p_node->X0();
            else if(j == 1)
                current_index_value = p_node->Y0();
            else
                current_index_value = p_node->Z0();
        }
    });

}

void ShapeUtils::SetNodalCoordinates(
            ContainerExpression<ModelPart::NodesContainerType>& rInputExpression)
{
    const IndexType local_size = rInputExpression.GetItemComponentCount();
    const auto& r_container = rInputExpression.GetContainer();
    const IndexType number_of_entities = rInputExpression.GetContainer().size();
    const auto& expression = rInputExpression.GetExpression();

    IndexPartition<IndexType>(r_container.size()).for_each([&](const IndexType Index){
        const auto p_node = (r_container.begin() + Index);
        const IndexType current_data_begin = Index * local_size;
        for (IndexType j = 0; j < local_size; ++j) {
            const double value = expression.Evaluate(Index, Index * local_size, j);
            if (j == 0){
                p_node->X0() = value;
                p_node->X() = value;
            }
            else if(j == 1){
                p_node->Y0() = value;
                p_node->Y() = value;
            }
            else{
                p_node->Z0() = value;
                p_node->Z() = value;
            }
        }
    });
}

void ShapeUtils::UpdateNodalCoordinates(
            ContainerExpression<ModelPart::NodesContainerType>& rInputExpression)
{
    const IndexType local_size = rInputExpression.GetItemComponentCount();
    const auto& r_container = rInputExpression.GetContainer();
    const IndexType number_of_entities = rInputExpression.GetContainer().size();
    const auto& expression = rInputExpression.GetExpression();

    IndexPartition<IndexType>(r_container.size()).for_each([&](const IndexType Index){
        const auto p_node = (r_container.begin() + Index);
        const IndexType current_data_begin = Index * local_size;
        for (IndexType j = 0; j < local_size; ++j) {
            const double value = expression.Evaluate(Index, Index * local_size, j);
            if (j == 0){
                p_node->X0() += value;
                p_node->X() += value;
            }
            else if(j == 1){
                p_node->Y0() += value;
                p_node->Y() += value;
            }
            else{
                p_node->Z0() += value;
                p_node->Z() += value;
            }
        }
    });
}

///@}
}