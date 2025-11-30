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
    : mConfiguration(Configuration)
{
    KRATOS_TRY

    this->mpContainer = pContainer;

    KRATOS_ERROR_IF_NOT(rDataShape.size() == 1)
        << "Invalid data shape. Following data shapes are allowed:\n\t[1]\n\t[2]\n\t[3]";

    DenseVector<unsigned int> tensor_shape(2);
    tensor_shape[0] = pContainer->size();
    tensor_shape[1] = rDataShape[0];

    KRATOS_ERROR_IF(tensor_shape[1] > 3)
        << "Invalid data shape. Following data shapes are allowed:\n\t[1]\n\t[2]\n\t[3]";

    this->mpStorage = Kratos::make_shared<Storage>(tensor_shape);

    KRATOS_CATCH("");
}

NodePositionTensorAdaptor::NodePositionTensorAdaptor(
    ModelPart::NodesContainerType::Pointer pContainer,
    Globals::Configuration Configuration)
    : NodePositionTensorAdaptor(pContainer, Configuration, {3})
{
}

NodePositionTensorAdaptor::NodePositionTensorAdaptor(
    const TensorAdaptor& rOther,
    Globals::Configuration Configuration,
    const bool Copy)
    : BaseType(rOther, Copy),
      mConfiguration(Configuration)
{
    KRATOS_TRY

    const auto& r_data_shape = this->DataShape();

    KRATOS_ERROR_IF_NOT(std::holds_alternative<ModelPart::NodesContainerType::Pointer>(this->GetContainer()))
        << "NodePositionTensorAdaptor can only be used with tensor data having nodal containers "
        << "[ tensor adaptor = " << rOther << " ].\n";

    KRATOS_ERROR_IF_NOT(r_data_shape.size() == 1)
        << "Incompatible number of dimensions in the provided tensor adaptor. NodePositionTensorAdaptor"
        << " requires a tensor data having 2 dimensions [ tensor data = " << this->Info() << " ].\n";

    KRATOS_ERROR_IF(r_data_shape[0] > 3)
        << "Invalid data shape. Only allowed to following data shapes:\n\t[1]\n\t[2]\n\t[3]";

    KRATOS_CATCH("");
}

void NodePositionTensorAdaptor::CollectData()
{
    KRATOS_TRY

    auto p_container = std::get<ModelPart::NodesContainerType::Pointer>(this->GetContainer());
    const auto& r_tensor_shape = this->Shape();

    KRATOS_ERROR_IF_NOT(r_tensor_shape[0] == p_container->size())
        << "Underlying container of the tensor data has changed size [ tensor data = "
        << this->Info() << ", container size = " << p_container->size() << " ].\n";

    switch (this->mConfiguration) {
        case Globals::Configuration::Current:
            ContainerIOUtils::CopyToContiguousArray<array_1d<double, 3>>(
                *p_container, this->ViewData(), r_tensor_shape.data().begin(),
                r_tensor_shape.data().begin() + r_tensor_shape.size(),
                [](auto& rCoordinates, const Node& rNode) {
                    rCoordinates = rNode.Coordinates();
                });
            break;
        case Globals::Configuration::Initial:
            ContainerIOUtils::CopyToContiguousArray<array_1d<double, 3>>(
                *p_container, this->ViewData(), r_tensor_shape.data().begin(),
                r_tensor_shape.data().begin() + r_tensor_shape.size(),
                [](auto& rCoordinates, const Node& rNode) {
                    rCoordinates = rNode.GetInitialPosition();
                });
    }

    KRATOS_CATCH("");
}

void NodePositionTensorAdaptor::StoreData()
{
    KRATOS_TRY

    auto p_container = std::get<ModelPart::NodesContainerType::Pointer>(this->GetContainer());
    const auto& r_tensor_shape = this->Shape();

    KRATOS_ERROR_IF_NOT(r_tensor_shape[0] == p_container->size())
        << "Underlying container of the tensor data has changed size [ tensor data = "
        << this->Info() << ", container size = " << p_container->size() << " ].\n";

    switch (mConfiguration) {
        case Globals::Configuration::Current:
            ContainerIOUtils::CopyFromContiguousDataArray<array_1d<double, 3>>(
                *p_container, this->ViewData(), r_tensor_shape.data().begin(),
                r_tensor_shape.data().begin() + r_tensor_shape.size(),
                [](Node& rNode) -> auto& {
                    // get the current coordinates
                    return rNode.Coordinates();
                });
            break;
        case Globals::Configuration::Initial:
            ContainerIOUtils::CopyFromContiguousDataArray<array_1d<double, 3>>(
                *p_container, this->ViewData(), r_tensor_shape.data().begin(),
                r_tensor_shape.data().begin() + r_tensor_shape.size(),
                [](Node& rNode) -> auto& {
                    // get the current coordinates
                    return rNode.GetInitialPosition();
                });
    }

    KRATOS_CATCH("");
}

std::string NodePositionTensorAdaptor::Info() const
{
    std::stringstream info;
    info << "NodePositionTensorAdaptor:";
    switch (mConfiguration) {
        case Globals::Configuration::Initial:
            info << " Configuration = Initial, ";
            break;
        case Globals::Configuration::Current:
            info << " Configuration = Current, ";
            break;
    }
    info << BaseType::Info();
    return info.str();
}

} // namespace Kratos