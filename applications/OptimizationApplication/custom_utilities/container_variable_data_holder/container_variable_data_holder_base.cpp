//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes
#include <sstream>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

// Application includes

// Include base h
#include "container_variable_data_holder_base.h"

namespace Kratos {

template <class TContainerType>
ContainerVariableDataHolderBase<TContainerType>::ContainerVariableDataHolderBase(ModelPart& rModelPart)
    : mDataDimension(0),
      mpModelPart(&rModelPart)
{
}

template <class TContainerType>
ContainerVariableDataHolderBase<TContainerType>::ContainerVariableDataHolderBase(
    const ContainerVariableDataHolderBase& rOther)
    : mDataDimension(rOther.mDataDimension),
      mpModelPart(rOther.mpModelPart)
{
    this->CopyDataFrom(rOther);
}
template <class TContainerType>
void ContainerVariableDataHolderBase<TContainerType>::CopyDataFrom(
    const ContainerVariableDataHolderBase<TContainerType>& rOther)
{
    KRATOS_ERROR_IF(this->GetModelPart() != rOther.GetModelPart())
        << "Modelpart mismatch between origin and destination for copy.\n"
        << "   Destination = " << this->GetModelPart().FullName() << "\n"
        << "   Origin = " << rOther.GetModelPart().FullName() << "\n.";

    if (this->mData.size() != rOther.mData.size()) {
        this->mData.resize(rOther.mData.size(), false);
    }

    this->mDataDimension = rOther.mDataDimension;

    IndexPartition<IndexType>(this->mData.size()).for_each([&](const IndexType Index) {
        this->mData[Index] = rOther.mData[Index];
    });
}

template <class TContainerType>
IndexType ContainerVariableDataHolderBase<TContainerType>::GetDataDimension() const
{
    return mDataDimension;
}

template <class TContainerType>
Vector& ContainerVariableDataHolderBase<TContainerType>::GetData()
{
    return mData;
}

template <class TContainerType>
const Vector& ContainerVariableDataHolderBase<TContainerType>::GetData() const
{
    return mData;
}

template <class TContainerType>
ModelPart& ContainerVariableDataHolderBase<TContainerType>::GetModelPart()
{
    return *mpModelPart;
}

template <class TContainerType>
const ModelPart& ContainerVariableDataHolderBase<TContainerType>::GetModelPart() const
{
    return *mpModelPart;
}

template <>
ModelPart::NodesContainerType& ContainerVariableDataHolderBase<ModelPart::NodesContainerType>::GetContainer()
{
    return mpModelPart->GetCommunicator().LocalMesh().Nodes();
}

template <>
ModelPart::ConditionsContainerType& ContainerVariableDataHolderBase<ModelPart::ConditionsContainerType>::GetContainer()
{
    return mpModelPart->GetCommunicator().LocalMesh().Conditions();
}

template <>
ModelPart::ElementsContainerType& ContainerVariableDataHolderBase<ModelPart::ElementsContainerType>::GetContainer()
{
    return mpModelPart->GetCommunicator().LocalMesh().Elements();
}

template <>
const ModelPart::NodesContainerType& ContainerVariableDataHolderBase<ModelPart::NodesContainerType>::GetContainer() const
{
    return mpModelPart->GetCommunicator().LocalMesh().Nodes();
}

template <>
const ModelPart::ConditionsContainerType& ContainerVariableDataHolderBase<ModelPart::ConditionsContainerType>::GetContainer() const
{
    return mpModelPart->GetCommunicator().LocalMesh().Conditions();
}

template <>
const ModelPart::ElementsContainerType& ContainerVariableDataHolderBase<ModelPart::ElementsContainerType>::GetContainer() const
{
    return mpModelPart->GetCommunicator().LocalMesh().Elements();
}

template <class TContainerType>
std::string ContainerVariableDataHolderBase<TContainerType>::Info() const
{
    std::stringstream msg;

    msg << "ContainerVariableDataHolderInfo: [ ModelPartName: "
        << this->GetModelPart().FullName()
        << ", Number of entities: " << this->GetContainer().size()
        << ", Data size: " << mData.size()
        << ", Data dimension: " << mDataDimension << " ].\n";

    return msg.str();
}

// template instantiations
template class ContainerVariableDataHolderBase<ModelPart::NodesContainerType>;
template class ContainerVariableDataHolderBase<ModelPart::ConditionsContainerType>;
template class ContainerVariableDataHolderBase<ModelPart::ElementsContainerType>;

} // namespace Kratos