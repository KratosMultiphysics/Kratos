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
#include <string>
#include <variant>
#include <type_traits>

// External includes

// Project includes
#include "tensor_adaptor_utils.h"
#include "utilities/data_type_traits.h"
#include "utilities/container_io_utils.h"

// Include base class
#include "variable_tensor_adaptor.h"

namespace Kratos {

namespace Detail
{
template<template<class> class TContainerIOType, class TDataType, class TContainerPointerType, class... TArgs>
TContainerIOType<TDataType> const * InitializeAndGetContainerIO(
    DenseVector<int>& rContainerShape,
    const std::vector<int>& rUserGivenShape,
    const Variable<TDataType> * pVariable,
    TContainerPointerType pContainer,
    TArgs&&... rArgs)
{
    using return_type = typename TContainerIOType<TDataType>::ReturnType;

    // creating the io
    auto p_container_io = new TContainerIOType<TDataType>(*pVariable, rArgs...);

    if (rUserGivenShape.size() != 1 || rUserGivenShape[0] != -1) {
        KRATOS_ERROR_IF_NOT(DataTypeTraits<return_type>::Dimension + 1 == rUserGivenShape.size())
            << "Dimensions mismatch for " << pVariable->Name()
            << " [ Required dimensions by variable = " << DataTypeTraits<return_type>::Dimension
            << ", shape = " << rUserGivenShape << " ].\n";
    }

    TensorAdaptorUtils::GetShape(rContainerShape, *pContainer, *p_container_io);
    if constexpr (!DataTypeTraits<TDataType>::IsDynamic) {
        // now we know the shape exactly because, rContainerShape should have the
        // exact sizes, hence we can check whether user has provided the shape
        // correctly. These types are double, array3, array4, ...

        if (rUserGivenShape.size() != 1 || rUserGivenShape[0] != -1) {
            for (IndexType i = 0; i < rUserGivenShape.size(); ++i) {
                KRATOS_ERROR_IF(rUserGivenShape[i] != -1 && rUserGivenShape[i] != rContainerShape[i])
                    << "Shape mismatch for " << pVariable->Name()
                    << " [ Required variable shape = " << rContainerShape
                    << ", shape = " << rUserGivenShape << " ].\n";
            }
        }
    }
    else {
        // here we may not know the exact sizes. so we use the user given sizes
        for (IndexType i = 0; i < rUserGivenShape.size(); ++i) {
            if (rUserGivenShape[i] != -1) {
                rContainerShape[i] = rUserGivenShape[i];
            }
        }
    }

    return p_container_io;
}
} // namespace Details

VariableTensorAdaptor::VariableTensorAdaptor(
    ModelPart::NodesContainerType::Pointer pContainer,
    VariableType pVariable,
    const int StepIndex,
    const std::vector<int>& rShape)
    : BaseType()
{
    mpContainerIO = std::visit([this, pContainer, StepIndex, &rShape](const auto pVariable) -> ContainerIOType {
        using data_type = typename std::remove_cv_t<std::decay_t<decltype(*pVariable)>>::Type;
        using primitive_data_type = typename DataTypeTraits<data_type>::PrimitiveType;

        auto p_container_io = Detail::InitializeAndGetContainerIO<HistoricalIO>(this->mShape, rShape, pVariable, pContainer, StepIndex);

        DenseVector<primitive_data_type> data;
        data.resize(TensorAdaptorUtils::GetFlatLength(this->mShape.data().begin(), this->mShape.data().end()));
        this->mData = std::move(data);

        return p_container_io;
    }, pVariable);

    mpContainer = pContainer;
}

VariableTensorAdaptor::VariableTensorAdaptor(
    ModelPart::NodesContainerType::Pointer pContainer,
    VariableType pVariable,
    const int StepIndex)
    : VariableTensorAdaptor(pContainer, pVariable, StepIndex, {-1})
{
}

template<class TContainerPointerType>
VariableTensorAdaptor::VariableTensorAdaptor(
    TContainerPointerType pContainer,
    VariableType pVariable,
    const std::vector<int>& rShape)
{
    mpContainerIO = std::visit([this, pContainer, &rShape](const auto pVariable) -> ContainerIOType {
        using data_type = typename std::remove_cv_t<std::decay_t<decltype(*pVariable)>>::Type;
        using primitive_data_type = typename DataTypeTraits<data_type>::PrimitiveType;

        auto p_container_io = Detail::InitializeAndGetContainerIO<NonHistoricalIO>(this->mShape, rShape, pVariable, pContainer);

        DenseVector<primitive_data_type> data;
        data.resize(TensorAdaptorUtils::GetFlatLength(this->mShape.data().begin(), this->mShape.data().end()));
        this->mData = std::move(data);

        return p_container_io;
    }, pVariable);

    mpContainer = pContainer;
}

template<class TContainerPointerType>
VariableTensorAdaptor::VariableTensorAdaptor(
    TContainerPointerType pContainer,
    VariableType pVariable)
    : VariableTensorAdaptor(pContainer, pVariable, {-1})
{
}

template<class TContainerPointerType>
VariableTensorAdaptor::VariableTensorAdaptor(
    TContainerPointerType pContainer,
    VariableType pVariable,
    const ProcessInfo& rProcessInfo,
    const std::vector<int>& rShape)
{
    mpContainerIO = std::visit([this, pContainer, &rProcessInfo, &rShape](const auto pVariable) -> ContainerIOType {
        using data_type = typename std::remove_cv_t<std::decay_t<decltype(*pVariable)>>::Type;
        using primitive_data_type = typename DataTypeTraits<data_type>::PrimitiveType;

        auto p_container_io = Detail::InitializeAndGetContainerIO<GaussPointIO>(this->mShape, rShape, pVariable, pContainer, rProcessInfo);

        DenseVector<primitive_data_type> data;
        data.resize(TensorAdaptorUtils::GetFlatLength(this->mShape.data().begin(), this->mShape.data().end()));
        this->mData = std::move(data);

        return p_container_io;
    }, pVariable);

    mpContainer = pContainer;
}

template<class TContainerPointerType>
VariableTensorAdaptor::VariableTensorAdaptor(
    TContainerPointerType pContainer,
    VariableType pVariable,
    const ProcessInfo& rProcessInfo)
    : VariableTensorAdaptor(pContainer, pVariable, rProcessInfo, {-1})
{
}

VariableTensorAdaptor::~VariableTensorAdaptor()
{
    std::visit([](const auto pContainerIO){ delete pContainerIO; }, this->mpContainerIO);
}

void VariableTensorAdaptor::CollectData()
{
    KRATOS_TRY

    std::visit([this](auto pContainer, auto pContainerIO){
        // sanity checks
        KRATOS_ERROR_IF_NOT(this->mShape[0] == static_cast<int>(pContainer->size()))
            << "First dimension of the initialized tensor adaptor mismatch with the container size [ "
            << "Tensor adapter shape = " << this->mShape << ", container size = " << pContainer->size()
            << ", TensorAdaptor = " << *this << " ].\n";

        using container_type = typename std::remove_cv_t<std::decay_t<decltype(*pContainer)>>;
        using container_io_type = std::remove_cv_t<std::decay_t<decltype(*pContainerIO)>>;
        using data_type = typename container_io_type::DataType;
        using primitive_type = typename DataTypeTraits<data_type>::PrimitiveType;

        auto& r_data = std::get<DenseVector<primitive_type>>(this->mData);

        if constexpr(std::is_same_v<container_io_type, HistoricalIO<data_type>>) {
            if constexpr(std::is_same_v<container_type, ModelPart::NodesContainerType>) {
                CopyToContiguousArray(*pContainer, *pContainerIO, r_data.data().begin(), r_data.size());
            }
        } else if constexpr(std::is_same_v<container_io_type, GaussPointIO<data_type>>) {
            if constexpr(std::is_same_v<container_type, ModelPart::ConditionsContainerType> || std::is_same_v<container_type, ModelPart::ElementsContainerType>) {
                CopyToContiguousArray(*pContainer, *pContainerIO, r_data.data().begin(), r_data.size());
            }
        } else {
            CopyToContiguousArray(*pContainer, *pContainerIO, r_data.data().begin(), r_data.size());
        }

    }, mpContainer, mpContainerIO);

    KRATOS_CATCH("");
}

void VariableTensorAdaptor::StoreData()
{
    KRATOS_TRY

    std::visit([this](auto pContainer, auto pContainerIO){
        // sanity checks
        KRATOS_ERROR_IF_NOT(this->mShape[0] == static_cast<int>(pContainer->size()))
            << "First dimension of the initialized tensor adaptor mismatch with the container size [ "
            << "Tensor adapter shape = " << this->mShape << ", container size = " << pContainer->size()
            << ", TensorAdaptor = " << *this << " ].\n";

        using container_type = typename std::remove_cv_t<std::decay_t<decltype(*pContainer)>>;
        using container_io_type = std::remove_cv_t<std::decay_t<decltype(*pContainerIO)>>;
        using data_type = typename container_io_type::DataType;
        using primitive_type = typename DataTypeTraits<data_type>::PrimitiveType;

        auto& r_data = std::get<DenseVector<primitive_type>>(this->mData);

        std::vector<unsigned int> shape;
        shape.resize(this->mShape.size() - 1);
        std::copy(this->mShape.begin() + 1, this->mShape.end(), shape.begin());


        if constexpr(std::is_same_v<container_io_type, HistoricalIO<data_type>>) {
            if constexpr(std::is_same_v<container_type, ModelPart::NodesContainerType>) {
                CopyFromContiguousDataArray(*pContainer, *pContainerIO, r_data.data().begin(), shape);
            }
        } else if constexpr(std::is_same_v<container_io_type, GaussPointIO<data_type>>) {
            if constexpr(std::is_same_v<container_type, ModelPart::ConditionsContainerType> || std::is_same_v<container_type, ModelPart::ElementsContainerType>) {
                CopyFromContiguousDataArray(*pContainer, *pContainerIO, r_data.data().begin(), shape);
            }
        } else {
            CopyFromContiguousDataArray(*pContainer, *pContainerIO, r_data.data().begin(), shape);
        }

    }, mpContainer, mpContainerIO);

    KRATOS_CATCH("");
}

VariableTensorAdaptor::ContainerType VariableTensorAdaptor::GetContainer() const
{
    return mpContainer;
}

std::string VariableTensorAdaptor::Info() const
{
    // return std::visit([this](auto pContainer) {
    //     std::stringstream msg;
    //     msg << "VariableTensorAdaptor: " << pContainer->Info() << " with " << this->mpContainer->size()
    //         << " " << ModelPart::Container<TContainerType>::GetEntityName() << "(s).";
    //     return msg.str();
    // }, mpContainerIO);
    return "";
}

// template instantiations
// non historical io
template VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::NodesContainerType::Pointer, VariableTensorAdaptor::VariableType);
template VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::NodesContainerType::Pointer, VariableTensorAdaptor::VariableType, const std::vector<int>&);
template VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, VariableTensorAdaptor::VariableType);
template VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, VariableTensorAdaptor::VariableType, const std::vector<int>&);
template VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::ElementsContainerType::Pointer, VariableTensorAdaptor::VariableType);
template VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::ElementsContainerType::Pointer, VariableTensorAdaptor::VariableType, const std::vector<int>&);
template VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::PropertiesContainerType::Pointer, VariableTensorAdaptor::VariableType);
template VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::PropertiesContainerType::Pointer, VariableTensorAdaptor::VariableType, const std::vector<int>&);
template VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::GeometriesMapType::Pointer, VariableTensorAdaptor::VariableType);
template VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::GeometriesMapType::Pointer, VariableTensorAdaptor::VariableType, const std::vector<int>&);

// gauss point io
template VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, VariableTensorAdaptor::VariableType, const ProcessInfo&);
template VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::ConditionsContainerType::Pointer, VariableTensorAdaptor::VariableType, const ProcessInfo&, const std::vector<int>&);
template VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::ElementsContainerType::Pointer, VariableTensorAdaptor::VariableType, const ProcessInfo&);
template VariableTensorAdaptor::VariableTensorAdaptor(ModelPart::ElementsContainerType::Pointer, VariableTensorAdaptor::VariableType, const ProcessInfo&, const std::vector<int>&);

} // namespace Kratos