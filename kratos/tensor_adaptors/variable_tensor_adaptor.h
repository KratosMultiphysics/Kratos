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

#pragma once

// System includes
#include <string>
#include <variant>

// External includes

// Project includes
#include "tensor_adaptor.h"
#include "utilities/container_io_utils.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief This class is used to read and write Variables of double
 *        type from any type of container and any type of \p TIOType.
 *        Examples are:
 *              Reading and writing from nodal historical data
 *              Reading and writing from data value containers for nodes, conditions, elements, properties container, ...
 *              Reading and writing from gauss points in conditions and elements containers.
 *              ...
 *
 * @tparam TIOType      Type of the IO to be used.
 * @tparam TIOArgs      Types of arguments required to construct an object of TIOType.
 */
template<template<class> class TIOType, class... TIOArgs>
class KRATOS_API(KRATOS_CORE) VariableTensorAdaptor: public TensorAdaptor<double> {
private:
    ///@name Private class definitions
    ///@{

    template<class... TDataType>
    struct VariantTypeHelper
    {
        using VariableVariantType = std::variant<Variable<TDataType> const *...>;
        using IOVariantType = std::variant<typename TIOType<TDataType>::Pointer...>;
    };

    ///@}
public:

    ///@name Type definitions
    ///@{

    using Pointer = Kratos::intrusive_ptr<VariableTensorAdaptor>;

    using ConstPointer = Kratos::intrusive_ptr<const VariableTensorAdaptor>;

    using BaseType = TensorAdaptor<double>;

    using ContainerType = typename BaseType::ContainerType;

    using VariantHelper = VariantTypeHelper<double,
                                            array_1d<double, 3>,
                                            array_1d<double, 4>,
                                            array_1d<double, 6>,
                                            array_1d<double, 9>,
                                            Vector,
                                            Matrix>;


    using VariableType = typename VariantHelper::VariableVariantType;

    using IOType = typename VariantHelper::IOVariantType;

    ///@}
    ///@name Life cycle
    ///@{

    /**
     * @brief Construct a new Variable Tensor Adaptor object
     * @details These constructors are only enabled if the given @p TContainerPointerType
     *          is a valid container type to be used with the TIOType.
     *
     * @tparam TContainerPointerType        Pointer type of the container.
     * @param pContainer                    Pointer to the container.
     * @param pVariable                     Variable to be used in the tensor adaptor for collecting and/or storing data.
     * @param rDataShape                    Required data shape. This should be a valid shape.
     * @param rArgs                         Arguments required to construct an object of @p TIOType.
     */
    template<
        class TContainerPointerType,
        typename = std::enable_if_t<TIOType<double>::template IsAllowedContainer<typename TContainerPointerType::element_type>>>
    VariableTensorAdaptor(
        TContainerPointerType pContainer,
        VariableType pVariable,
        const std::vector<unsigned int>& rDataShape,
        TIOArgs&&... rArgs)
    {
        std::visit([&](auto pVariable) {
            auto p_container_io = this->CreateIO(*pVariable, rArgs...);
            this->InitImpl(rDataShape, pContainer, p_container_io);
        }, pVariable);
    }
    /**
     * @brief Construct a new Variable Tensor Adaptor object
     * @details These constructors are only enabled if the given @p TContainerPointerType
     *          is a valid container type to be used with the TIOType.
     *
     * @tparam TContainerPointerType        Pointer type of the container.
     * @param pContainer                    Pointer to the container.
     * @param pVariable                     Variable to be used in the tensor adaptor for collecting and/or storing data.
     * @param rArgs                         Arguments required to construct an object of @p TIOType.
     */
    template<
        class TContainerPointerType,
        typename = std::enable_if_t<TIOType<double>::template IsAllowedContainer<typename TContainerPointerType::element_type>>>
    VariableTensorAdaptor(
        TContainerPointerType pContainer,
        VariableType pVariable,
        TIOArgs&&... rArgs)
    {
        std::visit([&](auto pVariable) {
            auto p_container_io = this->CreateIO(*pVariable, rArgs...);
            this->InitImpl(this->GetDataShape(*p_container_io, *pContainer), pContainer, p_container_io);
        }, pVariable);
    }

    // Destructor
    ~VariableTensorAdaptor() override = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Fill the internal data from Kratos data structures
     */
    void CollectData() override
    {
        std::visit([this](auto pContainer, auto pContainerIO){
            this->CollectDataImpl(*pContainer, *pContainerIO);
        }, mpContainer, mpIO);
    }

    /**
     * @brief Store internal data to the given TContainerType container.
     *
     */
    void StoreData() override
    {
        std::visit([this](auto pContainer, auto pContainerIO){
            this->StoreDataImpl(*pContainer, *pContainerIO);
        }, mpContainer, mpIO);
    }

    ContainerType GetContainer() const override
    {
        return mpContainer;
    }

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override
    {
        return std::visit([this](auto pContainer, auto pContainerIO) {
            return this->InfoImpl(*pContainer, *pContainerIO);
        }, mpContainer, mpIO);
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    IOType mpIO;

    ContainerType mpContainer;

    ///@}
    ///@name Private operations
    ///@{

    template<class TDataType, class... TArgs>
    typename TIOType<TDataType>::Pointer CreateIO(
        const Variable<TDataType>& rVariable,
        TArgs&&... rArgs)
    {
        return Kratos::make_shared<TIOType<TDataType>>(rVariable, rArgs...);
    }

    template<class TCurrentIOType, class TContainerType>
    std::vector<unsigned int> GetDataShape(
        const TCurrentIOType& rIO,
        const TContainerType& rContainer)
    {
        using return_type = typename TCurrentIOType::ReturnType;

        return_type dummy{};

        if (!rContainer.empty()) {
            rIO.GetValue(dummy, rContainer.front());
        }

        return DataTypeTraits<return_type>::template Shape<unsigned int>(dummy);
    }

    template<class TContainerPointerType, class TCurrentIOPointerType>
    void InitImpl(
        const std::vector<unsigned int>& rDataShape,
        TContainerPointerType pContainer,
        TCurrentIOPointerType pIO)
    {
        using return_type = typename TCurrentIOPointerType::element_type::ReturnType;

        KRATOS_ERROR_IF_NOT(DataTypeTraits<return_type>::IsValidShape(rDataShape.data(), rDataShape.data() + rDataShape.size()))
            << "Invalid data shape provided. [ data shape provided = " << rDataShape
            << ", max possible sizes in each dimension  = "
            << DataTypeTraits<return_type>::Shape(return_type{}) << " ].\n";

        // setting the tensor shape
        this->mShape.resize(rDataShape.size() + 1);
        std::copy(rDataShape.begin(), rDataShape.end(), this->mShape.begin() + 1);
        this->mShape[0] = pContainer->size();

        this->mpIO = pIO;
        this->mpContainer = pContainer;
        this->mData.resize(this->Size(), false);
    }

    template<class TContainerType, class TContainerIOType>
    void CollectDataImpl(
        const TContainerType& rContainer,
        const TContainerIOType& rContainerIO)
    {
        KRATOS_TRY

        // sanity checks
        KRATOS_ERROR_IF_NOT(this->mShape[0] == rContainer.size())
            << "First dimension of the initialized tensor adaptor mismatch with the container size [ "
            << "Tensor adapter shape = " << this->mShape << ", container size = " << rContainer.size()
            << ", TensorAdaptor = " << *this << " ].\n";

        KRATOS_ERROR_IF_NOT(this->mData.size() == this->Size())
            << "Data container size mismatch [ mData.size() = "
            << this->mData.size() << ", shape size = " << this->Size()
            << ", shape = " << this->mShape << " ].\n";

        // The case with following constexpr becoming false will never be reached
        // at runtime because, this class's constructors are enabled only
        // for the TContainerTypes which the TIOType is allowed.
        if constexpr(TContainerIOType::template IsAllowedContainer<TContainerType>) {
            CopyToContiguousArray(rContainer, rContainerIO, this->mData.data().begin(), this->mShape.begin(), this->mShape.end());
        }

        KRATOS_CATCH("");
    }

    template<class TContainerType, class TContainerIOType>
    void StoreDataImpl(
        TContainerType& rContainer,
        const TContainerIOType& rContainerIO)
    {
        KRATOS_TRY

        // sanity checks
        KRATOS_ERROR_IF_NOT(this->mShape[0] == rContainer.size())
            << "First dimension of the initialized tensor adaptor mismatch with the container size [ "
            << "Tensor adapter shape = " << this->mShape << ", container size = " << rContainer.size()
            << ", TensorAdaptor = " << *this << " ].\n";

        KRATOS_ERROR_IF_NOT(this->mData.size() == this->Size())
            << "Data container size mismatch [ mData.size() = "
            << this->mData.size() << ", shape size = " << this->Size()
            << ", shape = " << this->mShape << " ].\n";

        const std::vector<unsigned int> shape(this->mShape.begin() + 1, this->mShape.end());

        // The case with following constexpr becoming false will never be reached
        // at runtime because, this class's constructors are enabled only
        // for the TContainerTypes which the TIOType is allowed.
        if constexpr(TContainerIOType::template IsAllowedContainer<TContainerType>) {
            CopyFromContiguousDataArray(rContainer, rContainerIO, this->mData.data().begin(), this->mShape.begin(), this->mShape.end());
        }

        KRATOS_CATCH("");
    }

    template<class TContainerType, class TContainerIOType>
    std::string InfoImpl(
        const TContainerType& rContainer,
        const TContainerIOType& rIO) const
    {
        std::stringstream msg;
        msg << "VariableTensorAdaptor: " << rIO.Info() << " with " << rContainer.size()
            << " " << ModelPart::Container<TContainerType>::GetEntityName() << "(s) having shape "
            << this->Shape() << " ].";
        return msg.str();
    }

    ///@}
};

/// @}
} // namespace Kratos