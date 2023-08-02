//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System incldues
#include <variant>

// Project includes
#include "includes/communicator.h"
#include "includes/mesh.h"
#include "expression/container_data_io.h"
#include "expression/variable_expression_data_io.h"
#include "expression/literal_flat_expression.h"
#include "expression/traits.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

namespace Kratos {

class ExpressionIOUtils
{
public:
    ///@name Public static operations
    ///@{

    template <class TContainerType, class TContainerDataIO, class TVariableType>
    static Expression::Pointer ReadToExpression(
        const TContainerType& rContainer,
        TVariableType pVariable)
    {
        KRATOS_TRY

        const IndexType number_of_entities = rContainer.size();

        if (number_of_entities != 0) {
            // initialize the shape with the first entity value
            return std::visit([&rContainer, number_of_entities](auto pVariable) {
                using data_type = typename std::remove_const_t<std::remove_pointer_t<decltype(pVariable)>>::Type;

                using raw_data_type = std::conditional_t<std::is_same_v<data_type, int>, int, double>;

                VariableExpressionDataIO<data_type> variable_flatten_data_io(TContainerDataIO::GetValue(*rContainer.begin(), *pVariable));

                auto p_expression = LiteralFlatExpression<raw_data_type>::Create(number_of_entities, variable_flatten_data_io.GetItemShape());
                auto& r_expression = *p_expression;

                IndexPartition<IndexType>(number_of_entities).for_each([&rContainer, &pVariable, &variable_flatten_data_io, &r_expression](const IndexType Index) {
                        const auto& values = TContainerDataIO::GetValue(*(rContainer.begin() + Index), *pVariable);
                        variable_flatten_data_io.Read(r_expression, Index, values);
                    });

                return Kratos::intrusive_ptr<Expression>(&r_expression);

            }, pVariable);
        }

        return nullptr;

        KRATOS_CATCH("");
    }

    template <class TContainerType, class TContainerDataIO, class TVariableType>
    static void WriteFromExpression(
        TContainerType& rContainer,
        Communicator& rCommunicator,
        const Expression& rExpression,
        TVariableType pVariable)
    {
        KRATOS_TRY

        const IndexType number_of_entities = rContainer.size();

        if (number_of_entities > 0) {
            KRATOS_TRY
            std::visit([&, number_of_entities](auto pVariable) {
                using data_type = typename std::remove_const_t<std::remove_pointer_t<decltype(pVariable)>>::Type;

                VariableExpressionDataIO<data_type> variable_flatten_data_io(rExpression.GetItemShape());

                // initialize the container variables first
                if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
                    // initializes ghost nodes as for the later synchronization
                    // only, the nodal non historical values needs to be set unless
                    // they are properly initialized. Otherwise, in synchronization, the variables will
                    // not be there in the ghost nodes hence seg faults.

                    // the vectors and matrices needs to be initialized in historical and non-historical
                    // data containers because they need to be initialized with the correct size for synchronization
                    if constexpr(std::is_same_v<data_type, Vector> || std::is_same_v<data_type, Matrix>) {
                        data_type dummy_value{};
                        variable_flatten_data_io.Assign(dummy_value, rExpression, 0);
                        if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::Historical>>) {
                            VariableUtils().SetVariable(*pVariable, dummy_value, rCommunicator.GhostMesh().Nodes());
                        } else if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::NonHistorical>>) {
                            VariableUtils().SetNonHistoricalVariable(*pVariable, dummy_value, rCommunicator.GhostMesh().Nodes());
                        }
                    } else {
                        // if it is a static type, then it only needs to be initialized in the non-historical container with zeros.
                        // historical container should be initialized with the default values when the container is created.
                        if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::NonHistorical>>) {
                            VariableUtils().SetNonHistoricalVariableToZero(*pVariable, rCommunicator.GhostMesh().Nodes());
                        }
                    }
                }

                IndexPartition<IndexType>(number_of_entities).for_each(data_type{}, [&rContainer, &pVariable, &rExpression, &variable_flatten_data_io](const IndexType Index, data_type& rValue){
                    variable_flatten_data_io.Assign(rValue, rExpression, Index);
                    TContainerDataIO::SetValue(*(rContainer.begin() + Index), *pVariable, rValue);
                });

                if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
                    // synchronize nodal values
                    if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::Historical>>) {
                        rCommunicator.SynchronizeVariable(*pVariable);
                    } else if constexpr(std::is_same_v<TContainerDataIO, ContainerDataIO<ContainerDataIOTags::NonHistorical>>) {
                        rCommunicator.SynchronizeNonHistoricalVariable(*pVariable);
                    }
                }
            }, pVariable);
            KRATOS_CATCH(" Variable: " + std::visit([](const auto& rpVariable) -> std::string {return rpVariable->Name();}, pVariable))
        }

        KRATOS_CATCH("");
    }

    static ModelPart::MeshType& GetMesh(
        Communicator& rCommunicator,
        MeshType  rMeshType)
    {
        switch (rMeshType) {
            case MeshType::Local: {
                    return rCommunicator.LocalMesh();
                    break;
                }
            case MeshType::Interface: {
                    return rCommunicator.InterfaceMesh();
                    break;
                }
            case MeshType::Ghost: {
                    return rCommunicator.GhostMesh();
                    break;
                }
            default: {
                KRATOS_ERROR << "Invalid mesh type";
                break;
            }
        }
    }

    static const ModelPart::MeshType& GetMesh(
        const Communicator& rCommunicator,
        MeshType  rMeshType)
    {
        switch (rMeshType) {
            case MeshType::Local: {
                    return rCommunicator.LocalMesh();
                    break;
                }
            case MeshType::Interface: {
                    return rCommunicator.InterfaceMesh();
                    break;
                }
            case MeshType::Ghost: {
                    return rCommunicator.GhostMesh();
                    break;
                }
            default: {
                KRATOS_ERROR << "Invalid mesh type";
                break;
            }
        }
    }

    ///@}

};

} // namespace VariableExpressionIOHelperUtilities
