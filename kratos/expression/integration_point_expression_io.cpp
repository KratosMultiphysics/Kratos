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

// System incldues
#include <type_traits>

// Project includes
#include "expression/expression_io_utils.h"
#include "expression/variable_expression_data_io.h"
#include "includes/process_info.h"

// Include base h
#include "integration_point_expression_io.h"

namespace Kratos {

namespace IntegrationPointExpressionIOUtils
{

template <class TContainerType, class TVariableType>
static Expression::Pointer ReadToExpression(
    TContainerType& rContainer,
    TVariableType pVariable,
    const DataCommunicator& rDataCommunicator,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    const IndexType number_of_entities = rContainer.size();

    return std::visit([&rContainer, &rDataCommunicator, &rProcessInfo, number_of_entities](auto pVariable) {
        using data_type = typename std::remove_const_t<std::remove_pointer_t<decltype(pVariable)>>::Type;
        using raw_data_type = std::conditional_t<std::is_same_v<data_type, int>, int, double>;

        // first get the shape correctly
        std::vector<int> shape_info;
        if (number_of_entities != 0) {
            // initialize the shape with the first entity value
            std::vector<data_type> output;
            rContainer.begin()->CalculateOnIntegrationPoints(*pVariable, output, rProcessInfo);
            VariableExpressionDataIO<std::vector<data_type>> variable_flatten_data_io(output);
            const auto& shape = variable_flatten_data_io.GetItemShape();
            shape_info.resize(shape.size() + 1, number_of_entities);
            std::transform(shape.begin(), shape.end(), shape_info.begin() + 1, [](const auto v) -> int { return v; });
        } else {
            shape_info.resize(1, number_of_entities);
        }

        // now communicate the shape
        const auto& r_shapes_info_in_ranks = rDataCommunicator.AllGatherv(shape_info);

        // get the shape
        std::vector<IndexType> shape;
        for (const auto& r_shape_info_in_rank : r_shapes_info_in_ranks) {
            KRATOS_ERROR_IF(r_shape_info_in_rank.empty()) << "Shape info needs to have minimum one dimension.";
            if (r_shape_info_in_rank[0] != 0) {
                shape.resize(r_shape_info_in_rank.size() - 1);
                std::transform(r_shape_info_in_rank.begin() + 1, r_shape_info_in_rank.end(), shape.begin(), [](const auto v) -> IndexType { return v;});
                break;
            }
        }

        // first entry in shape_info contains the number of entities in the rContainer.
        // therefore, remove number of entities for check
        shape_info.erase(shape_info.begin());

        // cross check between all ranks the shape is the same
        IndexType local_index = 0;
        KRATOS_ERROR_IF(number_of_entities > 0 && !std::all_of(shape.begin(), shape.end(), [&local_index, &shape_info](const auto v) { return static_cast<int>(v) == shape_info[local_index++]; }))
            << "All the ranks should have values with the same shape.\n";

        auto p_expression = LiteralFlatExpression<raw_data_type>::Create(number_of_entities, shape);
        auto& r_expression = *p_expression;

        if (number_of_entities != 0) {
            // initialize the shape with the first entity value
            std::vector<data_type> values;
            rContainer.begin()->CalculateOnIntegrationPoints(*pVariable, values, rProcessInfo);
            VariableExpressionDataIO<std::vector<data_type>> variable_flatten_data_io(values);

            IndexPartition<IndexType>(number_of_entities).for_each(values, [&rContainer, &pVariable, &variable_flatten_data_io, &r_expression, &rProcessInfo](const IndexType Index, auto& rTLS) {
                    (rContainer.begin() + Index)->CalculateOnIntegrationPoints(*pVariable, rTLS, rProcessInfo);
                    variable_flatten_data_io.Read(r_expression, Index, rTLS);
                });
        }

        return Kratos::intrusive_ptr<Expression>(&r_expression);

    }, pVariable);

    KRATOS_CATCH("");
}

template <class TContainerType, class TVariableType>
static void WriteFromExpression(
    TContainerType& rContainer,
    Communicator& rCommunicator,
    const Expression& rExpression,
    const ProcessInfo& rProcessInfo,
    TVariableType pVariable)
{
    KRATOS_TRY

    const IndexType number_of_entities = rContainer.size();

    std::visit([&, number_of_entities](auto pVariable) {
        KRATOS_TRY

        using data_type = typename std::remove_const_t<std::remove_pointer_t<decltype(pVariable)>>::Type;

        VariableExpressionDataIO<std::vector<data_type>> variable_flatten_data_io(rExpression.GetItemShape());

        // initialize the container variables first
        IndexPartition<IndexType>(number_of_entities).for_each(std::vector<data_type>{}, [&rContainer, &pVariable, &rExpression, &variable_flatten_data_io, &rProcessInfo](const IndexType Index, std::vector<data_type>& rValue){
            variable_flatten_data_io.Assign(rValue, rExpression, Index);
            (rContainer.begin() + Index)->SetValuesOnIntegrationPoints(*pVariable, rValue, rProcessInfo);
        });

        KRATOS_CATCH(" Variable: " + pVariable->Name())

    }, pVariable);

    KRATOS_CATCH("");
}
} // IntegrationPointExpressionIOUtils

IntegrationPointExpressionIO::Input::Input(
    ModelPart& rModelPart,
    const VariableType& rVariable,
    Globals::DataLocation rCurrentLocation,
    MeshType CurrentMeshType)
    : mpModelPart(&rModelPart),
      mpVariable(rVariable),
      mDataLocation(rCurrentLocation),
      mMeshType(CurrentMeshType)
{
}

Expression::Pointer IntegrationPointExpressionIO::Input::Execute() const
{
    auto& r_mesh = ExpressionIOUtils::GetMesh(mpModelPart->GetCommunicator(), mMeshType);
    const auto& r_data_communicator = mpModelPart->GetCommunicator().GetDataCommunicator();

    switch (mDataLocation) {
        case Globals::DataLocation::Condition:
            return IntegrationPointExpressionIOUtils::ReadToExpression<ModelPart::ConditionsContainerType, const VariableType>(r_mesh.Conditions(), mpVariable, r_data_communicator, mpModelPart->GetProcessInfo());
        case Globals::DataLocation::Element:
            return IntegrationPointExpressionIOUtils::ReadToExpression<ModelPart::ElementsContainerType, const VariableType>(r_mesh.Elements(), mpVariable, r_data_communicator, mpModelPart->GetProcessInfo());
        default:
            KRATOS_ERROR << "Invalid container type. Only supports Condition and Element.";
    }

    return nullptr;
}

IntegrationPointExpressionIO::Output::Output(
    ModelPart& rModelPart,
    const VariableType& rVariable,
    Globals::DataLocation CurrentLocation,
    MeshType  CurrentMeshType)
    : mpModelPart(&rModelPart),
      mpVariable(rVariable),
      mDataLocation(CurrentLocation),
      mMeshType(CurrentMeshType)
{
}

void IntegrationPointExpressionIO::Output::Execute(const Expression& rExpression)
{
    KRATOS_TRY
    auto& r_communicator = mpModelPart->GetCommunicator();
    auto& r_mesh = ExpressionIOUtils::GetMesh(r_communicator, mMeshType);

    switch (mDataLocation) {
        case Globals::DataLocation::Condition:
            IntegrationPointExpressionIOUtils::WriteFromExpression<ModelPart::ConditionsContainerType, const VariableType>(r_mesh.Conditions(), r_communicator, rExpression, mpModelPart->GetProcessInfo(), mpVariable);
            break;
        case Globals::DataLocation::Element:
            IntegrationPointExpressionIOUtils::WriteFromExpression<ModelPart::ElementsContainerType, const VariableType>(r_mesh.Elements(), r_communicator, rExpression, mpModelPart->GetProcessInfo(), mpVariable);
            break;
        default:
            KRATOS_ERROR << "Invalid container type. Only supports Condition and Element";
    }
    KRATOS_CATCH("");
}

template<class TContainerType, MeshType TMeshType>
void IntegrationPointExpressionIO::Read(
    ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
    const VariableType& rVariable)
{
    auto p_expression =
        Input(rContainerExpression.GetModelPart(), rVariable,
                                std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>
                                    ? Globals::DataLocation::Condition
                                    : Globals::DataLocation::Element,
                                TMeshType)
            .Execute();

    rContainerExpression.SetExpression(p_expression);
}

template<class TContainerType, MeshType TMeshType>
void IntegrationPointExpressionIO::Write(
    const ContainerExpression<TContainerType, TMeshType>& rContainerExpression,
    const VariableType& rVariable)
{
    Output(*rContainerExpression.pGetModelPart(), rVariable,
                             std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>
                                 ? Globals::DataLocation::Condition
                                 : Globals::DataLocation::Element,
                             TMeshType)
        .Execute(rContainerExpression.GetExpression());
}

#define KRATOS_INSTANTIATE_ENTITY_CONTAINER_IO_METHODS(CONTAINER_TYPE, MESH_TYPE)                                                    \
    template KRATOS_API(KRATOS_CORE) void IntegrationPointExpressionIO::Read(ContainerExpression<CONTAINER_TYPE, MESH_TYPE>&, const IntegrationPointExpressionIO::VariableType&); \
    template KRATOS_API(KRATOS_CORE) void IntegrationPointExpressionIO::Write(const ContainerExpression<CONTAINER_TYPE, MESH_TYPE>&, const IntegrationPointExpressionIO::VariableType&);\

#define KRATOS_INSTANTIATE_CONTAINER_VARIABLE_EXPRESSION_IO(MESH_TYPE)                              \
    KRATOS_INSTANTIATE_ENTITY_CONTAINER_IO_METHODS(ModelPart::ConditionsContainerType, MESH_TYPE)   \
    KRATOS_INSTANTIATE_ENTITY_CONTAINER_IO_METHODS(ModelPart::ElementsContainerType, MESH_TYPE)

KRATOS_INSTANTIATE_CONTAINER_VARIABLE_EXPRESSION_IO(MeshType::Local)
KRATOS_INSTANTIATE_CONTAINER_VARIABLE_EXPRESSION_IO(MeshType::Interface)
KRATOS_INSTANTIATE_CONTAINER_VARIABLE_EXPRESSION_IO(MeshType::Ghost)

#undef KRATOS_INSTANTIATE_CONTAINER_VARIABLE_EXPRESSION_IO
#undef KRATOS_INSTANTIATE_ENTITY_CONTAINER_IO_METHODS

} // namespace Kratos
