//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes
#include <cmath>
#include <numeric>

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"
#include "utilities/atomic_utilities.h"
#include "utilities/data_type_traits.h"

// Application includes

// Include base h
#include "optimization_utils.h"

namespace Kratos
{

namespace OptimizationHelperUtils {

Properties::Pointer CloneProperties(
    IndexType& rPropertiesId,
    const Properties& rOriginProperties,
    ModelPart& rModelPart,
    const bool IsRecursive)
{
    KRATOS_TRY

    // first create the base properties
    auto p_cloned_properties = rModelPart.CreateNewProperties(rPropertiesId);
    // assign the base property values from the origin properties
    *p_cloned_properties = rOriginProperties;
    p_cloned_properties->SetId(rPropertiesId);
    ++rPropertiesId;

    if (IsRecursive) {
        const auto& r_origin_sub_properties_list = rOriginProperties.GetSubProperties();
        for (IndexType i_prop = 0; i_prop < r_origin_sub_properties_list.size(); ++i_prop) {
            auto& p_cloned_sub_property = *((p_cloned_properties->GetSubProperties().begin() + i_prop).base());
            p_cloned_sub_property = CloneProperties(
                rPropertiesId, *(r_origin_sub_properties_list.begin() + i_prop),
                rModelPart, IsRecursive);
        }
    }

    return p_cloned_properties;

    KRATOS_CATCH("");
}

template<class TDataType>
void UpdatePropertiesVariableRecursively(
    Properties& rProperties,
    const Variable<TDataType>& rVariable,
    const TDataType& rValue)
{
    KRATOS_TRY

    rProperties.SetValue(rVariable, rValue);

    for (auto& r_sub_properties : rProperties.GetSubProperties()) {
        UpdatePropertiesVariableRecursively(r_sub_properties, rVariable, rValue);
    }

    KRATOS_CATCH("");
}

} // namespace OptimizationHelperUtils

template<class TContainerType>
GeometryData::KratosGeometryType OptimizationUtils::GetContainerEntityGeometryType(
    const TContainerType& rContainer,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    if (rContainer.size() > 0) {
        const auto first_geometry_type = rContainer.begin()->GetGeometry().GetGeometryType();
        const bool local_value = block_for_each<MinReduction<bool>>(rContainer, [&](const auto& rEntity) -> IndexType {
            return first_geometry_type == rEntity.GetGeometry().GetGeometryType();
        });

        if (rDataCommunicator.AndReduceAll(local_value)) {
            return first_geometry_type;
        } else {
            return GeometryData::KratosGeometryType::Kratos_generic_type;
        }
    } else {
        return GeometryData::KratosGeometryType::Kratos_generic_type;
    }

    KRATOS_CATCH("")
}

template<class TContainerType, class TDataType>
bool OptimizationUtils::IsVariableExistsInAllContainerProperties(
    const TContainerType& rContainer,
    const Variable<TDataType>& rVariable,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    const bool local_value = block_for_each<MinReduction<bool>>(rContainer, [&](const auto& rEntity){
        return rEntity.GetProperties().Has(rVariable);
    });

    return rDataCommunicator.AndReduceAll(local_value);

    KRATOS_CATCH("");
}

template<class TContainerType, class TDataType>
bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(
    const TContainerType& rContainer,
    const Variable<TDataType>& rVariable,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    const bool local_value = block_for_each<MaxReduction<bool>>(rContainer, [&](const auto& rEntity){
        return rEntity.GetProperties().Has(rVariable);
    });

    return rDataCommunicator.OrReduceAll(local_value);

    KRATOS_CATCH("");
}

template<class TContainerType>
void OptimizationUtils::CreateEntitySpecificPropertiesForContainer(
    ModelPart& rModelPart,
    TContainerType& rContainer,
    const bool IsRecursive)
{
    KRATOS_TRY

    auto element_properties_id = block_for_each<MaxReduction<IndexType>>(rContainer, [](const auto& rEntity) {
        return rEntity.GetProperties().Id();
    });

    auto properties_id = block_for_each<MaxReduction<IndexType>>(rModelPart.GetRootModelPart().PropertiesArray(), [](const auto pProperties) {
        return pProperties->Id();
    });

    properties_id = std::max(element_properties_id, properties_id) + 1;

    for (auto& r_entity : rContainer) {
        r_entity.SetProperties(OptimizationHelperUtils::CloneProperties(
            properties_id, r_entity.GetProperties(), rModelPart, IsRecursive));
    }

    KRATOS_CATCH("");
}

template<class TContainerType, class TDataType>
void OptimizationUtils::UpdatePropertiesVariableWithRootValueRecursively(
    TContainerType& rContainer,
    const Variable<TDataType>& rVariable)
{
    KRATOS_TRY

    block_for_each(rContainer, [&rVariable](auto& rEntity) {
        auto& r_properties = rEntity.GetProperties();
        const auto& r_root_value = r_properties[rVariable];
        OptimizationHelperUtils::UpdatePropertiesVariableRecursively(r_properties, rVariable, r_root_value);
    });

    KRATOS_CATCH("");
}

template<>
KRATOS_API(OPTIMIZATION_APPLICATION) IndexType OptimizationUtils::GetVariableDimension(
    const Variable<double>& rVariable,
    const IndexType DomainSize)
{
    return 1;
}

template<>
KRATOS_API(OPTIMIZATION_APPLICATION) IndexType OptimizationUtils::GetVariableDimension(
    const Variable<array_1d<double, 3>>& rVariable,
    const IndexType DomainSize)
{
    return DomainSize;
}

void OptimizationUtils::SetSolutionStepVariablesList(
    ModelPart& rDestinationModelPart,
    const ModelPart& rOriginModelPart)
{
    rDestinationModelPart.SetNodalSolutionStepVariablesList(rOriginModelPart.pGetNodalSolutionStepVariablesList());
}

bool OptimizationUtils::IsSolutionStepVariablesListASubSet(
    const ModelPart& rMainSetModelPart,
    const ModelPart& rSubSetModelPart)
{
    for (const auto& r_sub_variable : rSubSetModelPart.GetNodalSolutionStepVariablesList()) {
        if (!rMainSetModelPart.GetNodalSolutionStepVariablesList().Has(r_sub_variable)) {
            return false;
        }
    }
    return true;
}

std::vector<std::string> OptimizationUtils::GetSolutionStepVariableNamesList(const ModelPart& rModelPart)
{
    KRATOS_TRY

    std::vector<std::string> result;
    for (const auto& r_var : rModelPart.GetNodalSolutionStepVariablesList()) {
        result.push_back(r_var.Name());
    }
    return result;

    KRATOS_CATCH("");
}

std::vector<std::vector<ModelPart*>> OptimizationUtils::GetComponentWiseModelParts(
    Model& rModel,
    Parameters Settings)
{
    KRATOS_TRY

    std::vector<std::vector<ModelPart*>> result;

    int number_of_components = -1;

    for (auto it = Settings.begin(); it != Settings.end(); ++it) {
        // first set the number of components by checking the first entry.
        if (number_of_components == -1) {
            number_of_components = it->size();
            result.resize(number_of_components);
        }

        KRATOS_ERROR_IF_NOT(static_cast<int>(it->size()) == number_of_components)
            << "Number of component mismatch [ Number of components required = "
            << number_of_components << ", number of components specified for model part "
            << it.name() << " = " << it->size() << " ]. Settings = \n" << Settings;

        IndexType i_component = 0;
        for (const auto& r_value : *it) {
            if (r_value.GetBool()) {
                result[i_component].push_back(&rModel.GetModelPart(it.name()));
            }
            ++i_component;
        }
    }

    return result;

    KRATOS_CATCH("");
}

TensorAdaptor<double>::Pointer OptimizationUtils::MapContainerDataToNodalData(
    const TensorAdaptor<double>& rInput,
    ModelPart::NodesContainerType::Pointer pNodes)
{
    KRATOS_TRY

    const auto& shape = rInput.Shape();
    const auto input_data = rInput.ViewData();

    DenseVector<unsigned int> output_shape(shape);
    output_shape[0] = pNodes->size();

    auto p_nd_data = Kratos::make_shared<NDData<double>>(output_shape, 0.0);
    auto p_result = Kratos::make_shared<TensorAdaptor<double>>(pNodes, p_nd_data, false);
    auto output_data = p_nd_data->ViewData();

    const auto number_of_components = std::accumulate(shape.begin() + 1, shape.end(), 1U, std::multiplies<IndexType>{});

    std::visit([&rInput, &input_data, &pNodes, &output_data, number_of_components](const auto pContainer) {
        using container_type = BareType<decltype(*pContainer)>;

        if constexpr(IsInList<container_type, ModelPart::ConditionsContainerType, ModelPart::ElementsContainerType>) {
            IndexPartition<IndexType>(pContainer->size()).for_each([&pContainer,&pNodes, &rInput, &output_data, &input_data, number_of_components](const auto EntityIndex) {
                const auto& r_geometry = (pContainer->begin() + EntityIndex)->GetGeometry();

                for (const auto& r_node : r_geometry) {
                    auto p_node_itr = pNodes->find(r_node.Id());
                    KRATOS_ERROR_IF(p_node_itr == pNodes->end())
                        << "The node with id = " << r_node.Id()
                        << " in entity with id  = " << (pContainer->begin() + EntityIndex)->Id()
                        << " not found in the provided neighbour entities container [ "
                        << " input tensor adaptor = " << rInput << " ].\n";

                    const auto node_index = std::distance(pNodes->begin(), pNodes->find(r_node.Id()));

                    for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
                        AtomicAdd(
                            output_data[node_index * number_of_components + i_comp],
                            input_data[EntityIndex * number_of_components + i_comp] / r_geometry.size());
                    }
                }
            });

        } else {
            KRATOS_ERROR << "The input tensor adaptor should have data based on conditions / elements [ input tensor adaptor = " << rInput << " ].\n";
        }
    }, rInput.GetContainer());

    return p_result;

    KRATOS_CATCH("");
}

template<class TContainerPointerType>
TensorAdaptor<double>::Pointer OptimizationUtils::MapNodalDataToContainerData(
    const TensorAdaptor<double>& rInput,
    TContainerPointerType pEntities,
    const TensorAdaptor<int>& rNeighbourCount)
{
    KRATOS_TRY

    const auto& shape = rInput.Shape();
    const auto input_data = rInput.ViewData();
    const auto neighbour_data = rNeighbourCount.ViewData();

    DenseVector<unsigned int> output_shape(shape);
    output_shape[0] = pEntities->size();

    auto p_nd_data = Kratos::make_shared<NDData<double>>(output_shape, 0.0);
    auto p_result = Kratos::make_shared<TensorAdaptor<double>>(pEntities, p_nd_data, false);
    auto output_data = p_nd_data->ViewData();

    const auto number_of_components = std::accumulate(shape.begin() + 1, shape.end(), 1U, std::multiplies<IndexType>{});

    if (!std::holds_alternative<ModelPart::NodesContainerType::Pointer>(rInput.GetContainer()) ||
        !std::holds_alternative<ModelPart::NodesContainerType::Pointer>(rNeighbourCount.GetContainer())) {
        KRATOS_ERROR
            << "Input container and neighbour container count tensor adaptors should be having nodal containers [ input tensor adaptor = "
            << rInput << ", neighbour count tensor adaptor = " << rNeighbourCount << " ].\n";
    }

    const auto& r_nodes = *std::get<ModelPart::NodesContainerType::Pointer>(rInput.GetContainer());

    if (&r_nodes != &*(std::get<ModelPart::NodesContainerType::Pointer>(rNeighbourCount.GetContainer()))) {
        KRATOS_ERROR
            << "The nodal containers in input and neighbour count tensor adaptors mismatch [ "
            << "input tensor adaptor = " << rInput << ", neighbour count tensor adaptor = " << rNeighbourCount << " ].\n";
    }

    IndexPartition<IndexType>(pEntities->size()).for_each([&input_data, &output_data, &neighbour_data, &r_nodes, &pEntities, &rInput, &rNeighbourCount, number_of_components](const auto EntityIndex) {
        const auto& r_geometry = (pEntities->begin() + EntityIndex)->GetGeometry();

        for (const auto& r_node : r_geometry) {
            auto p_node_itr = r_nodes.find(r_node.Id());

            KRATOS_ERROR_IF(p_node_itr == r_nodes.end())
                << "The node with id = " << r_node.Id() << " in entity with id = "
                << (pEntities->begin() + EntityIndex)->Id() << " is not found in the nodes of the input tensor adaptor [ "
                << "input tensor adaptor = " << rInput << ", neighbour count tensor adaptor = " << rNeighbourCount << " ].\n";

            const IndexType node_index = std::distance(r_nodes.begin(), p_node_itr);

            for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
                AtomicAdd(output_data[EntityIndex * number_of_components + i_comp],
                          input_data[node_index * number_of_components + i_comp] / neighbour_data[node_index]);
            }
        }
    });

    return p_result;

    KRATOS_CATCH("");
}

// template instantiations
#ifndef KRATOS_OPTIMIZATION_UTILS_UPDATE_PROPERTIES_VARIABLE_RECURSIVELY
#define KRATOS_OPTIMIZATION_UTILS_UPDATE_PROPERTIES_VARIABLE_RECURSIVELY(CONTAINER_TYPE)                                                                                                                                \
template KRATOS_API(OPTIMIZATION_APPLICATION) void OptimizationUtils::UpdatePropertiesVariableWithRootValueRecursively<CONTAINER_TYPE, bool>(CONTAINER_TYPE&, const Variable<bool>&);                                   \
template KRATOS_API(OPTIMIZATION_APPLICATION) void OptimizationUtils::UpdatePropertiesVariableWithRootValueRecursively<CONTAINER_TYPE, int>(CONTAINER_TYPE&, const Variable<int>&);                                     \
template KRATOS_API(OPTIMIZATION_APPLICATION) void OptimizationUtils::UpdatePropertiesVariableWithRootValueRecursively<CONTAINER_TYPE, double>(CONTAINER_TYPE&, const Variable<double>&);                               \
template KRATOS_API(OPTIMIZATION_APPLICATION) void OptimizationUtils::UpdatePropertiesVariableWithRootValueRecursively<CONTAINER_TYPE, array_1d<double, 3>>(CONTAINER_TYPE&, const Variable<array_1d<double, 3>>&);     \
template KRATOS_API(OPTIMIZATION_APPLICATION) void OptimizationUtils::UpdatePropertiesVariableWithRootValueRecursively<CONTAINER_TYPE, array_1d<double, 4>>(CONTAINER_TYPE&, const Variable<array_1d<double, 4>>&);     \
template KRATOS_API(OPTIMIZATION_APPLICATION) void OptimizationUtils::UpdatePropertiesVariableWithRootValueRecursively<CONTAINER_TYPE, array_1d<double, 6>>(CONTAINER_TYPE&, const Variable<array_1d<double, 6>>&);     \
template KRATOS_API(OPTIMIZATION_APPLICATION) void OptimizationUtils::UpdatePropertiesVariableWithRootValueRecursively<CONTAINER_TYPE, array_1d<double, 9>>(CONTAINER_TYPE&, const Variable<array_1d<double, 9>>&);     \
template KRATOS_API(OPTIMIZATION_APPLICATION) void OptimizationUtils::UpdatePropertiesVariableWithRootValueRecursively<CONTAINER_TYPE, Matrix>(CONTAINER_TYPE&, const Variable<Matrix>&);                               \
template KRATOS_API(OPTIMIZATION_APPLICATION) void OptimizationUtils::UpdatePropertiesVariableWithRootValueRecursively<CONTAINER_TYPE, Vector>(CONTAINER_TYPE&, const Variable<Vector>&);                               \

#endif

KRATOS_OPTIMIZATION_UTILS_UPDATE_PROPERTIES_VARIABLE_RECURSIVELY(ModelPart::ConditionsContainerType);
KRATOS_OPTIMIZATION_UTILS_UPDATE_PROPERTIES_VARIABLE_RECURSIVELY(ModelPart::ElementsContainerType);

#undef KRATOS_OPTIMIZATION_UTILS_UPDATE_PROPERTIES_VARIABLE_RECURSIVELY

template KRATOS_API(OPTIMIZATION_APPLICATION) GeometryData::KratosGeometryType OptimizationUtils::GetContainerEntityGeometryType(const ModelPart::ConditionsContainerType&, const DataCommunicator&);
template KRATOS_API(OPTIMIZATION_APPLICATION) GeometryData::KratosGeometryType OptimizationUtils::GetContainerEntityGeometryType(const ModelPart::ElementsContainerType&, const DataCommunicator&);

template KRATOS_API(OPTIMIZATION_APPLICATION) void OptimizationUtils::CreateEntitySpecificPropertiesForContainer(ModelPart&, ModelPart::ConditionsContainerType&, const bool);
template KRATOS_API(OPTIMIZATION_APPLICATION) void OptimizationUtils::CreateEntitySpecificPropertiesForContainer(ModelPart&, ModelPart::ElementsContainerType&, const bool);

template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ElementsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);
template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);

template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ElementsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);
template KRATOS_API(OPTIMIZATION_APPLICATION) bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);

template KRATOS_API(OPTIMIZATION_APPLICATION) TensorAdaptor<double>::Pointer OptimizationUtils::MapNodalDataToContainerData(const TensorAdaptor<double>&, ModelPart::ConditionsContainerType::Pointer, const TensorAdaptor<int>&);
template KRATOS_API(OPTIMIZATION_APPLICATION) TensorAdaptor<double>::Pointer OptimizationUtils::MapNodalDataToContainerData(const TensorAdaptor<double>&, ModelPart::ElementsContainerType::Pointer, const TensorAdaptor<int>&);

}