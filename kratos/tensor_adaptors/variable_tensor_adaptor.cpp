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

        auto p_container_io = TensorAdaptorUtils::InitializeAndGetContainerIO<HistoricalIO>(this->mShape, rShape, pVariable, pContainer, StepIndex);

        DenseVector<primitive_data_type> data;
        data.resize(this->Size());
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

        auto p_container_io = TensorAdaptorUtils::InitializeAndGetContainerIO<NonHistoricalIO>(this->mShape, rShape, pVariable, pContainer);

        DenseVector<primitive_data_type> data;
        data.resize(this->Size());
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

        auto p_container_io = TensorAdaptorUtils::InitializeAndGetContainerIO<GaussPointIO>(this->mShape, rShape, pVariable, pContainer, rProcessInfo);

        DenseVector<primitive_data_type> data;
        data.resize(this->Size());
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

        if constexpr(container_io_type::template IsAllowedContainer<container_type>::value) {
            CopyToContiguousArray(*pContainer, *pContainerIO, r_data.data().begin(), r_data.size());
        } else {
            KRATOS_ERROR << "It is prohibited to use " << ModelPart::Container<container_type>::GetEntityName()
                         << " containers with " <<  pContainerIO->Info () << ".";
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


        if constexpr(container_io_type::template IsAllowedContainer<container_type>::value) {
            CopyFromContiguousDataArray(*pContainer, *pContainerIO, r_data.data().begin(), shape);
        } else {
            KRATOS_ERROR << "It is prohibited to use " << ModelPart::Container<container_type>::GetEntityName()
                         << " containers with " <<  pContainerIO->Info () << ".";
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
    return std::visit([this](auto pContainer, auto pContainerIO) {
        std::stringstream msg;
        msg << "VariableTensorAdaptor: " << pContainerIO->Info() << " with " << pContainer->size()
            << " " << ModelPart::Container<std::remove_cv_t<std::decay_t<decltype(*pContainer)>>>::GetEntityName() << "(s) having shape "
            << this->Shape() << " ].";
        return msg.str();
    }, mpContainer, mpContainerIO);
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