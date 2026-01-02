//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes
#include <vector>
#include <type_traits>

// External includes

// Project includes
#include "tensor_adaptors/index_tensor_adaptor.h"
#include "utilities/parallel_utilities.h"
#include "tensor_adaptors/tensor_adaptor_utils.h"

namespace Kratos {

namespace {

template <class TContainerType>
constexpr bool IsSupportedContainer()
{
    return std::is_same_v<TContainerType, ModelPart::GeometryContainerType> ||
           std::is_same_v<TContainerType, ModelPart::ElementsContainerType> ||
           std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>;
}

template <class TEntity>
const auto& GetGeometry(const TEntity& rEntity)
{
    if constexpr (std::is_same_v<TEntity, Geometry<Node>>) {
        return rEntity;
    } else {
        return rEntity.GetGeometry();
    }
}

template <class TContainerType>
DenseVector<unsigned int> GetShape(const TContainerType& rContainer)
{
    if (rContainer.empty()) {
        return DenseVector<unsigned int>(0);
    }

    const auto& r_first = *rContainer.begin();
    const std::size_t num_nodes = GetGeometry(r_first).size();

    DenseVector<unsigned int> shape(2);
    shape[0] = static_cast<unsigned int>(rContainer.size());
    shape[1] = static_cast<unsigned int>(num_nodes);
    return shape;
}

template <class TContainerType>
void CheckContainer(const TContainerType& rContainer)
{
    if (rContainer.empty()) {
        return;
    }
    const auto& r_first = *rContainer.begin();
    const std::size_t num_nodes = GetGeometry(r_first).size();

    IndexPartition<std::size_t>(rContainer.size()).for_each([&](std::size_t i) {
        const auto& r_entity = *(rContainer.begin() + i);
        const auto& r_geometry = GetGeometry(r_entity);
        KRATOS_ERROR_IF(r_geometry.size() != num_nodes) 
             << "Entity #" << r_entity.Id() << " (Index: " << i << ") has " << r_geometry.size() 
             << " nodes, but expected " << num_nodes << "." << std::endl;
    });
}

template<class TContainerType>
void CollectIds(Kratos::span<int> Span, const TContainerType& rContainer, const std::size_t NumNodes)
{
    IndexPartition<std::size_t>(rContainer.size()).for_each([&](std::size_t i) {
        auto it = rContainer.begin() + i;
        const auto& r_geometry = GetGeometry(*it);

        // Check compatibility
        KRATOS_DEBUG_ERROR_IF(r_geometry.size() != NumNodes) 
            << "Entity at index " << i << " has " << r_geometry.size() 
            << " nodes, but expected " << NumNodes << "." << std::endl;

        for (std::size_t n = 0; n < NumNodes; ++n) {
            Span[i * NumNodes + n] = static_cast<int>(r_geometry[n].Id());
        }
    });
}



} // namespace

IndexTensorAdaptor::IndexTensorAdaptor(ContainerPointerType pContainer)
    : TensorAdaptor<int>()
{
    mpContainer = pContainer;
    std::visit([this](auto&& p_container) {
        using ContainerType = std::decay_t<decltype(*p_container)>;
        if constexpr (IsSupportedContainer<ContainerType>()) {
            mpStorage = Kratos::make_shared<NDData<int>>(GetShape(*p_container));
        } else {
             KRATOS_ERROR << "Unsupported container type for IndexTensorAdaptor. Only Geometries, Elements, and Conditions are supported." << std::endl;
        }
    }, *mpContainer);
}

IndexTensorAdaptor::IndexTensorAdaptor(
    const TensorAdaptor& rOther,
    const bool Copy)
    : TensorAdaptor<int>(rOther, Copy)
{
}

void IndexTensorAdaptor::Check() const
{
    std::visit([](auto&& p_container) {
        using ContainerType = std::decay_t<decltype(*p_container)>;
        if constexpr (IsSupportedContainer<ContainerType>()) {
            CheckContainer(*p_container);
        }
    }, *mpContainer);
}

void IndexTensorAdaptor::CollectData()
{
    if (!mpStorage) {
        // Should not happen if constructed correctly
        return;
    }
    
    const auto& shape = mpStorage->Shape();
    if(shape.size() < 2) return; // Should be [N, M]
    const std::size_t num_nodes = shape[1];

    std::visit([&](auto&& p_container) {
        using ContainerType = std::decay_t<decltype(*p_container)>;
        if constexpr (IsSupportedContainer<ContainerType>()) {
            CollectIds(this->ViewData(), *p_container, num_nodes);
        }
    }, *mpContainer);
}

void IndexTensorAdaptor::StoreData()
{
    KRATOS_ERROR << "StoreData is not implemented for IndexTensorAdaptor. Renumbering via TensorAdaptor is disabled." << std::endl;
}

std::string IndexTensorAdaptor::Info() const
{
    return "IndexTensorAdaptor";
}

} // namespace Kratos
