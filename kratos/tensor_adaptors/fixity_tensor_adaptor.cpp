//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes
#include <type_traits>
#include <vector>

// External includes

// Project includes
#include "tensor_adaptors/fixity_tensor_adaptor.h"
#include "tensor_adaptors/tensor_adaptor_utils.h"
#include "utilities/parallel_utilities.h"

namespace Kratos {

DenseVector<unsigned int> FixityTensorAdaptor::GetShape(const ModelPart::NodesContainerType& rContainer)
{
    if (rContainer.empty()) {
        return zero_vector<unsigned int>(mDofsVarPointerList.size());
    }

    DenseVector<unsigned int> shape(2);
    shape[0] = static_cast<unsigned int>(rContainer.size());
    shape[1] = static_cast<unsigned int>(mDofsVarPointerList.size());
    return shape;
}

FixityTensorAdaptor::FixityTensorAdaptor(
    ContainerPointerType pContainer,
    const std::vector<const Variable<double>*>& rDofsVarPointerList)
    : TensorAdaptor<bool>()
    , mDofsVarPointerList(rDofsVarPointerList)
{
    // Set the container and the storage
    // Note that only nodes container is supported for fixity tensor adaptor as fixity only applies to nodes
    mpContainer = pContainer;
    std::visit(
        [this](auto &&p_container) {
            using ContainerType = BareType<decltype(*p_container)>;
            if constexpr (std::is_same<ContainerType, ModelPart::NodesContainerType>::value) {
                mpStorage = Kratos::make_shared<NDData<bool>>(GetShape(*p_container));
            } else {
                KRATOS_ERROR << "Unsupported container type for FixityTensorAdaptor. Only Nodes are supported." << std::endl;
            }
        },
        *mpContainer);
}

FixityTensorAdaptor::FixityTensorAdaptor(
    const BaseType &rOther,
    const bool Copy)
    : TensorAdaptor<bool>(rOther, Copy)
{
}

void FixityTensorAdaptor::Check() const
{
    std::visit(
        [&](auto &&p_container) {
            using ContainerType = BareType<decltype(*p_container)>;
            if constexpr (std::is_same<ContainerType, ModelPart::NodesContainerType>::value) {
                CheckContainer(this->mDofsVarPointerList, *p_container);
            }
        },
        *mpContainer);
}

TensorAdaptor<bool>::Pointer FixityTensorAdaptor::Clone() const
{
    return Kratos::make_shared<FixityTensorAdaptor>(*this);
}

void FixityTensorAdaptor::CollectData()
{
    std::visit(
        [&](auto &&p_container) {
            using ContainerType = BareType<decltype(*p_container)>;
            if constexpr (std::is_same<ContainerType, ModelPart::NodesContainerType>::value) {
                CollectFixityData(this->ViewData(), this->mDofsVarPointerList, *p_container);
            }
        },
        *mpContainer);
}

void FixityTensorAdaptor::StoreData()
{
    std::visit(
        [&](auto &&p_container) {
            using ContainerType = BareType<decltype(*p_container)>;
            if constexpr (std::is_same<ContainerType, ModelPart::NodesContainerType>::value) {
                StoreFixityData(this->ViewData(), this->mDofsVarPointerList, *p_container);
            }
        },
        *mpContainer);
}

std::string FixityTensorAdaptor::Info() const
{
  return "FixityTensorAdaptor";
}

void FixityTensorAdaptor::CollectFixityData(
    Kratos::span<bool> Span,
    const std::vector<const Variable<double>*>& rDofsVarPointerList,
    const ModelPart::NodesContainerType& rNodes)
{
    for (std::size_t i_var = 0; i_var < rDofsVarPointerList.size(); ++i_var) {
        auto p_var = rDofsVarPointerList[i_var];
        IndexPartition<std::size_t>(rNodes.size()).for_each([&](std::size_t i_node) {
            auto it = rNodes.begin() + i_node;
            Span[i_node * rDofsVarPointerList.size() + i_var] = it->IsFixed(*p_var);
        });
    }
}

void FixityTensorAdaptor::StoreFixityData(
    Kratos::span<bool> Span,
    const std::vector<const Variable<double>*>& rDofsVarPointerList,
    const ModelPart::NodesContainerType& rNodes)
{
    for (std::size_t i_var = 0; i_var < rDofsVarPointerList.size(); ++i_var) {
        auto p_var = rDofsVarPointerList[i_var];
        IndexPartition<std::size_t>(rNodes.size()).for_each([&](std::size_t i_node) {
            auto it = rNodes.begin() + i_node;
            if (Span[i_node * rDofsVarPointerList.size() + i_var]) {
                it->Fix(*p_var);
            } else {
                it->Free(*p_var);
            }
        });
    }
}

void FixityTensorAdaptor::CheckContainer(
    const std::vector<const Variable<double>*>& rDofsVarPointerList,
    const ModelPart::NodesContainerType& rNodes)
{
    if (rNodes.empty()) {
        return;
    }

    for (std::size_t i_var = 0; i_var < rDofsVarPointerList.size(); ++i_var) {
        auto p_var = rDofsVarPointerList[i_var];
        IndexPartition<std::size_t>(rNodes.size()).for_each([&](std::size_t i_node) {
            auto it = rNodes.begin() + i_node;
            KRATOS_ERROR_IF_NOT(it->HasDofFor(*p_var)) << "Node " << it->Id() << " does not have DOF for variable " << p_var->Name() << "." << std::endl;
        });
    }
}

} // namespace Kratos
