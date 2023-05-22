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
//

// System includes
#include <cmath>

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "optimization_application.h"

// Include base h
#include "implicit_filter_utils.h"

namespace Kratos
{

void ImplicitFilterUtils::CalculateNodeNeighbourCount(
    ModelPart& rModelPart)
{
    // Calculate number of neighbour elements for each node.
    VariableUtils().SetNonHistoricalVariableToZero(NUMBER_OF_NEIGHBOUR_ELEMENTS, rModelPart.Nodes());
    block_for_each(rModelPart.Elements(), [&](ModelPart::ElementType& rElement) {
        auto& r_geometry = rElement.GetGeometry();
        for (unsigned j = 0; j < r_geometry.PointsNumber(); ++j) {
            double& r_num_neighbour =
                r_geometry[j].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
            AtomicAdd(r_num_neighbour, 1.0);
        }
    });
    rModelPart.GetCommunicator().AssembleNonHistoricalData(NUMBER_OF_NEIGHBOUR_ELEMENTS);
}

void ImplicitFilterUtils::SetBulkRadiusForShapeFiltering(
    ModelPart& rModelPart)
{
    ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
    rCurrentProcessInfo.SetValue(HELMHOLTZ_BULK_RADIUS_SHAPE,1.0);

    const double local_volume_strain_energy = block_for_each<SumReduction<double>>(rModelPart.Elements(), [&](auto& rElement) {
        double elem_val;
        rElement.Calculate(ELEMENT_STRAIN_ENERGY,elem_val,rCurrentProcessInfo);
        return elem_val;
    });

    const double local_surface_strain_energy = block_for_each<SumReduction<double>>(rModelPart.Conditions(), [&](auto& rCondition) {
        double cond_val;
        rCondition.Calculate(ELEMENT_STRAIN_ENERGY,cond_val,rCurrentProcessInfo);
        return cond_val;
    });

    double bulk_filter_size = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_surface_strain_energy);
    bulk_filter_size /= rModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_volume_strain_energy);

    rCurrentProcessInfo.SetValue(HELMHOLTZ_BULK_RADIUS_SHAPE,bulk_filter_size);

}

void ImplicitFilterUtils::AssignVectorNodalExpressionToScalarVariable(const Variable<double>& rVariable, const SpecializedContainerExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>& rNodalContainer, int cIndex)
{
    KRATOS_TRY

    auto& r_container = rNodalContainer.GetContainer();
    const IndexType number_of_entities = r_container.size();

    if (number_of_entities > 0) {
        const auto& r_expression = rNodalContainer.GetExpression();
        auto shape = r_expression.GetShape();

        KRATOS_ERROR_IF(shape.size()!=1 || shape[0] !=3)
        << "AssignVectorNodalExpressionToScalarVariable: the given container is not vector container\n";

        VariableExpressionDataIO<Vector> variable_flatten_data_io(shape);

        //initialize the container variables first
        VariableUtils().SetNonHistoricalVariableToZero(rVariable, rNodalContainer.GetModelPart().GetCommunicator().GhostMesh().Nodes());

        IndexPartition<IndexType>(number_of_entities).for_each(Vector{}, [&r_container, &rVariable, &r_expression, &variable_flatten_data_io, &cIndex](const IndexType Index, Vector& rValue){
            variable_flatten_data_io.Assign(rValue, r_expression, Index);
            ContainerDataIO<ContainerDataIOTags::NonHistorical>::SetValue(*(r_container.begin() + Index), rVariable, rValue[cIndex]);
        });

        auto& r_communicator = rNodalContainer.GetModelPart().GetCommunicator();
        // r_communicator.SynchronizeNonHistoricalVariable(rVariable);
    }

    KRATOS_CATCH("");
}

void ImplicitFilterUtils::AssignScalarVariableToVectorNodalExpression(const Variable<double>& rVariable,  SpecializedContainerExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>& rNodalContainer, int cIndex)
{
    KRATOS_TRY

    auto& r_container = rNodalContainer.GetContainer();
    const IndexType number_of_entities = r_container.size();

    if (number_of_entities > 0) {
        auto& r_expression = rNodalContainer.pGetExpression();
        auto shape = r_expression->GetShape();

        KRATOS_ERROR_IF(shape.size()!=1 || shape[0] !=3)
        << "AssignScalarVariableToVectorNodalExpression: the given container is not vector container\n";

        VariableExpressionDataIO<Vector> variable_flatten_data_io(shape);

        IndexPartition<IndexType>(number_of_entities).for_each(Vector{}, [&r_container, &rNodalContainer, &rVariable, &r_expression, &variable_flatten_data_io, &cIndex](const IndexType Index, Vector& rValue){
            variable_flatten_data_io.Assign(rValue, *r_expression, Index);
            rValue[cIndex] = ContainerDataIO<ContainerDataIOTags::NonHistorical>::GetValue(*(r_container.begin() + Index), rVariable);
            //assign the value back to the ????????
        });

        auto& r_communicator = rNodalContainer.GetModelPart().GetCommunicator();
        // r_communicator.SynchronizeNonHistoricalVariable(rVariable);
    }

    KRATOS_CATCH("");
}

}