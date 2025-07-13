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

// External includes

// Project includes
#include "tensor_adaptor_utils.h"
#include "utilities/container_io_utils.h"

// Include base h
#include "node_position_tensor_adaptor.h"

namespace Kratos {

NodePositionTensorAdaptor::NodePositionTensorAdaptor(
    ModelPart::NodesContainerType::Pointer pContainer,
    Globals::Configuration Configuration,
    const std::vector<unsigned int>& rDataShape)
    : mpContainer(pContainer),
      mConfiguration(Configuration)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rDataShape.size() == 1)
        << "Invalid data shape. Only allowed to following data shapes:\n\t[1]\n\t[2]\n\t[3]";

    DenseVector<unsigned int> tensor_shape(2);
    tensor_shape[0] = pContainer->size();
    tensor_shape[1] = rDataShape[0];

    KRATOS_ERROR_IF(tensor_shape[1] > 3)
        << "Invalid data shape. Only allowed to following data shapes:\n\t[1]\n\t[2]\n\t[3]";

    this->SetShape(tensor_shape);

    KRATOS_CATCH("");
}

NodePositionTensorAdaptor::NodePositionTensorAdaptor(
    ModelPart::NodesContainerType::Pointer pContainer,
    Globals::Configuration Configuration)
    : NodePositionTensorAdaptor(pContainer, Configuration, {3})
{
}

void NodePositionTensorAdaptor::CollectData()
{
    const auto& r_shape = this->Shape();
    switch (mConfiguration) {
        case Globals::Configuration::Current:
            ContainerIOUtils::CopyToContiguousArray<array_1d<double, 3>>(
                *mpContainer, this->ViewData(), r_shape.data().begin(),
                r_shape.data().begin() + r_shape.size(),
                [](auto& rCoordinates, const Node& rNode) {
                    noalias(rCoordinates) = rNode.Coordinates();
                });
            break;
        case Globals::Configuration::Initial:
            ContainerIOUtils::CopyToContiguousArray<array_1d<double, 3>>(
                *mpContainer, this->ViewData(), r_shape.data().begin(),
                r_shape.data().begin() + r_shape.size(),
                [](auto& rCoordinates, const Node& rNode) {
                    noalias(rCoordinates) = rNode.GetInitialPosition();
                });
    }
}

void NodePositionTensorAdaptor::StoreData()
{
    const auto& r_shape = this->Shape();
    switch (mConfiguration) {
        case Globals::Configuration::Current:
            ContainerIOUtils::CopyFromContiguousDataArray<array_1d<double, 3>>(
                *mpContainer, this->ViewData(), r_shape.data().begin(),
                r_shape.data().begin() + r_shape.size(),
                [&r_shape](Node& rNode) -> auto& {
                    // get the current coordinates
                    return rNode.Coordinates();
                });
            break;
        case Globals::Configuration::Initial:
            ContainerIOUtils::CopyFromContiguousDataArray<array_1d<double, 3>>(
                *mpContainer, this->ViewData(), r_shape.data().begin(),
                r_shape.data().begin() + r_shape.size(),
                [&r_shape](Node& rNode) -> auto& {
                    // get the current coordinates
                    return rNode.GetInitialPosition();
                });
    }
}

NodePositionTensorAdaptor::ContainerPointerType NodePositionTensorAdaptor::GetContainer() const
{
    return mpContainer;
}

std::string NodePositionTensorAdaptor::Info() const
{
    switch (mConfiguration) {
        case Globals::Configuration::Initial:
            return TensorAdaptorUtils::Info("NodePositionTensorAdaptor configuration = Initial, ", this->Shape(), *mpContainer);
        case Globals::Configuration::Current:
            return TensorAdaptorUtils::Info("NodePositionTensorAdaptor configuration = Current, ", this->Shape(), *mpContainer);
    }
}

} // namespace Kratos