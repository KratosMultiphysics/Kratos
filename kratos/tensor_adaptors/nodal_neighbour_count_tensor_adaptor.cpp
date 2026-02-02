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

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "utilities/atomic_utilities.h"
#include "utilities/data_type_traits.h"

// Include base h
#include "nodal_neighbour_count_tensor_adaptor.h"

namespace Kratos {

template<class TContainerPointerType>
NodalNeighbourCountTensorAdaptor::NodalNeighbourCountTensorAdaptor(
    ModelPart::NodesContainerType::Pointer pNodes,
    TContainerPointerType pEntityContainer)
    : mpEntityContainer(pEntityContainer)
{
    this->mpContainer = pNodes;
    this->mpStorage = Kratos::make_shared<Storage>(DenseVector<unsigned int>(1, pNodes->size()));
}

template<class TContainerPointerType>
NodalNeighbourCountTensorAdaptor::NodalNeighbourCountTensorAdaptor(
    const TensorAdaptor& rOther,
    TContainerPointerType pEntityContainer,
    const bool Copy)
    : BaseType(rOther, Copy),
      mpEntityContainer(pEntityContainer)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(std::holds_alternative<ModelPart::NodesContainerType::Pointer>(rOther.GetContainer()))
        << "NodalNeighbourCountTensorAdaptor can only be used with tensor data having nodal containers. "
        << "[ tensor adaptor = " << rOther << " ].\n";

    const auto& shape = rOther.Shape();
    KRATOS_ERROR_IF(shape.size() != 1 || shape[0] != this->Size())
        << "NodalNeighbourCountTensorAdaptor needs to have tensor data with shape [n, 1] where n being number of nodes."
        << "[ tensor adaptor shape = " << rOther.Shape() << ", tensor adaptor = " << rOther << " ].\n";

    KRATOS_CATCH("");
}

void NodalNeighbourCountTensorAdaptor::CollectData()
{
    KRATOS_TRY

    auto& r_nodes = *std::get<ModelPart::NodesContainerType::Pointer>(this->GetContainer());
    auto data_view = this->ViewData();

    std::unordered_map<IndexType, int> id_index_map;
    for (IndexType i = 0; i < r_nodes.size(); ++i) {
        data_view[i] = 0;
        id_index_map[(r_nodes.begin() + i)->Id()] = i;
    }

    std::visit([&data_view, &id_index_map](auto pContainer){
        using container_type = BareType<decltype(*pContainer)>;
        block_for_each(*pContainer, [&data_view, &id_index_map](const auto& rEntity) {
            for (const auto& r_node : rEntity.GetGeometry()) {
                auto itr = id_index_map.find(r_node.Id());

                KRATOS_ERROR_IF(itr == id_index_map.end())
                    << "The node with id " << r_node.Id() << " in the "
                    << ModelPart::Container<container_type>::GetEntityName() << " with id " << rEntity.Id()
                    << " is not found in the nodes list.";

                AtomicAdd<int>(data_view[itr->second], 1);
            }
        });
    }, mpEntityContainer);

    KRATOS_CATCH("");
}

void NodalNeighbourCountTensorAdaptor::StoreData()
{
    KRATOS_ERROR << "StoreData method is not implemented.";
}

std::string NodalNeighbourCountTensorAdaptor::Info() const
{
    std::stringstream info;
    info << "NodalNeighbourCountTensorAdaptor:";
    std::visit([&info](auto pContainer) {
        info << " number of  " << ModelPart::Container<BareType<decltype(*pContainer)>>::GetEntityName() << " = " << pContainer->size();
    }, this->mpEntityContainer);
    info << ", " << BaseType::Info();
    return info.str();
}

template KRATOS_API(KRATOS_CORE) NodalNeighbourCountTensorAdaptor::NodalNeighbourCountTensorAdaptor(ModelPart::NodesContainerType::Pointer, ModelPart::ConditionsContainerType::Pointer);
template KRATOS_API(KRATOS_CORE) NodalNeighbourCountTensorAdaptor::NodalNeighbourCountTensorAdaptor(ModelPart::NodesContainerType::Pointer, ModelPart::ElementsContainerType::Pointer);
template KRATOS_API(KRATOS_CORE) NodalNeighbourCountTensorAdaptor::NodalNeighbourCountTensorAdaptor(const TensorAdaptor&, ModelPart::ConditionsContainerType::Pointer, const bool);
template KRATOS_API(KRATOS_CORE) NodalNeighbourCountTensorAdaptor::NodalNeighbourCountTensorAdaptor(const TensorAdaptor&, ModelPart::ElementsContainerType::Pointer, const bool);

} // namespace Kratos