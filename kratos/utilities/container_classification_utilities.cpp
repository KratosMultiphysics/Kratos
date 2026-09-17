#include "utilities/container_classification_utilities.h"
#include <map>
#include <typeindex>
#include <type_traits>

namespace Kratos {

namespace {
    // Helper function to extract GeometryType enum across Elements, Conditions, and Geometries
    template <typename TContainerType>
    GeometryData::KratosGeometryType ExtractGeometryType(const typename TContainerType::data_type& rEntity) {
        if constexpr (std::is_same_v<TContainerType, ModelPart::GeometryContainerType>) {
            return rEntity.GetGeometryType();
        } else {
            return rEntity.GetGeometry().GetGeometryType();
        }
    }
}

template <typename TContainerType>
std::vector<TContainerType> ContainerClassificationUtilities::Classify(TContainerType& rContainer) 
{
    using EntityType = typename TContainerType::data_type;

    std::vector<TContainerType> result;

    if constexpr (std::is_same_v<TContainerType, ModelPart::GeometryContainerType>) {
        // Geometries are grouped directly by GeometryType enum
        std::map<GeometryData::KratosGeometryType, PointerVector<EntityType>> groups;

        for (auto it = rContainer.ptr_begin(); it != rContainer.ptr_end(); ++it) {
            groups[(*it)->GetGeometryType()].push_back(*it);
        }

        result.reserve(groups.size());
        for (auto& [key, ptr_vec] : groups) {
            TContainerType sub_container;
            sub_container.insert(ptr_vec.ptr_begin(), ptr_vec.ptr_end());
            result.push_back(std::move(sub_container));
        }
    } else {
        // Elements/Conditions grouped by (Entity dynamic type, GeometryType enum)
        using KeyType = std::pair<std::type_index, GeometryData::KratosGeometryType>;
        std::map<KeyType, PointerVector<EntityType>> groups;

        for (auto it = rContainer.ptr_begin(); it != rContainer.ptr_end(); ++it) {
            KeyType key{
                std::type_index(typeid(**it)),
                (*it)->GetGeometry().GetGeometryType()
            };
            groups[key].push_back(*it);
        }

        result.reserve(groups.size());
        for (auto& [key, ptr_vec] : groups) {
            TContainerType sub_container;
            sub_container.insert(ptr_vec.ptr_begin(), ptr_vec.ptr_end());
            result.push_back(std::move(sub_container));
        }
    }

    return result;
}

template <typename TContainerType>
std::vector<TContainerType> ContainerClassificationUtilities::ClassifyByGeometryType(TContainerType& rContainer) 
{
    using EntityType = typename TContainerType::data_type;
    constexpr std::size_t num_geom_types = static_cast<std::size_t>(
        GeometryData::KratosGeometryType::NumberOfGeometryTypes);

    // Direct indexed lookup vector acting as a perfect hash table
    std::vector<PointerVector<EntityType>> indexed_groups(num_geom_types);

    for (auto it = rContainer.ptr_begin(); it != rContainer.ptr_end(); ++it) {
        const auto geom_type = ExtractGeometryType<TContainerType>(**it);
        const std::size_t idx = static_cast<std::size_t>(geom_type);
        indexed_groups[idx].push_back(*it);
    }

    // Vector pre-sized to match total geometry types enum count
    std::vector<TContainerType> result(num_geom_types);

    IndexPartition<std::size_t>(num_geom_types).for_each([&](std::size_t i) {
        if (!indexed_groups[i].empty()) {
            result[i].insert(indexed_groups[i].ptr_begin(), indexed_groups[i].ptr_end());
        }
    });

    return result;
}

// Explicit template instantiations
template std::vector<ModelPart::ElementsContainerType> 
ContainerClassificationUtilities::Classify(ModelPart::ElementsContainerType&);

template std::vector<ModelPart::ConditionsContainerType> 
ContainerClassificationUtilities::Classify(ModelPart::ConditionsContainerType&);

template std::vector<ModelPart::GeometryContainerType> 
ContainerClassificationUtilities::Classify(ModelPart::GeometryContainerType&);

template std::vector<ModelPart::ElementsContainerType> 
ContainerClassificationUtilities::ClassifyByGeometryType(ModelPart::ElementsContainerType&);

template std::vector<ModelPart::ConditionsContainerType> 
ContainerClassificationUtilities::ClassifyByGeometryType(ModelPart::ConditionsContainerType&);

template std::vector<ModelPart::GeometryContainerType> 
ContainerClassificationUtilities::ClassifyByGeometryType(ModelPart::GeometryContainerType&);

} // namespace Kratos