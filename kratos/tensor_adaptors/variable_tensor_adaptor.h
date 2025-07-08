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
#include "tensor_adaptor_utils.h"
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
        : mpContainer(pContainer)
    {
        std::visit([&](auto pVariable) {
            auto p_io = this->CreateIO(*pVariable, rArgs...);

            // set the shape
            this->mShape.resize(rDataShape.size() + 1);
            this->mShape[0] = pContainer->size();
            std::copy(rDataShape.begin(), rDataShape.end(), this->mShape.begin() + 1);

            TensorAdaptorUtils::InitializeData(*pContainer, *p_io, this->Shape(), this->mData);
            mpIO = p_io;
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
        : mpContainer(pContainer)
    {
        std::visit([&](auto pVariable) {
            auto p_io = this->CreateIO(*pVariable, rArgs...);

            // set the shape
            const auto& r_data_shape = this->GetDataShape(*p_io, *pContainer);
            this->mShape.resize(r_data_shape.size() + 1);
            this->mShape[0] = pContainer->size();
            std::copy(r_data_shape.begin(), r_data_shape.end(), this->mShape.begin() + 1);

            TensorAdaptorUtils::InitializeData(*pContainer, *p_io, this->Shape(), this->mData);
            mpIO = p_io;
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
        std::visit([this](auto pContainer, auto pIO){
            TensorAdaptorUtils::CollectData(*pContainer, *pIO, this->Shape(), this->mData);
        }, mpContainer, mpIO);
    }

    /**
     * @brief Store internal data to the given TContainerType container.
     *
     */
    void StoreData() override
    {
        std::visit([this](auto pContainer, auto pIO){
            TensorAdaptorUtils::StoreData(*pContainer, *pIO, this->Shape(), this->mData);
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
        return std::visit([this](auto pContainer, auto pIO) {
            return TensorAdaptorUtils::Info(*pContainer, *pIO, this->Shape());
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

    ///@}
};

/// @}
} // namespace Kratos