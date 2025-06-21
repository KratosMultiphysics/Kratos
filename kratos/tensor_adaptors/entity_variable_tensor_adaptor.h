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

// External includes

// Project includes
#include "tensor_adaptor.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Base class or all the tensor adaptor types.
 */
template<class TIOType, class... TArgs>
class EntityVariableTensorAdaptor: public TensorAdaptor<typename TIOType::ContainerType, typename TIOType::PrimitiveDataType> {
public:

    ///@name Type definitions
    ///@{

    using Pointer = Kratos::intrusive_ptr<EntityVariableTensorAdaptor>;

    using ConstPointer = Kratos::intrusive_ptr<const EntityVariableTensorAdaptor>;

    using BaseType = TensorAdaptor<typename TIOType::ContainerType, typename TIOType::PrimitiveDataType>;

    using PrimitiveDataType = typename BaseType::PrimitiveDataType;

    using ContainerType = typename BaseType::ContainerType;

    ///@}
    ///@name Life cycle
    ///@{

    EntityVariableTensorAdaptor(
        typename ContainerType::Pointer pContainer,
        typename TIOType::VariableType pVariable,
        TArgs&&... rArgs);

    ~EntityVariableTensorAdaptor() = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Fill the internal data from Kratos data structures
     */
    void CollectData() override;

    /**
     * @brief Store internal data to the given TContainerType container.
     *
     */
    void StoreData() override;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    typename TIOType::VariableType mpVariable;

    typename TIOType::ContainerIOType mContainerIO;

    ///@}
};

/// @}
} // namespace Kratos