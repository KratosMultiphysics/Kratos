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
#include "tensor_adaptors/tensor_adaptor.h"
#include "tensor_adaptors/tensor_adaptor_utils.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @class PropertiesVariableTensorAdaptor
 * @ingroup TensorAdaptors
 * @brief Adaptor class for handling values represented as variables in Kratos entities' properties.
 *
 * @details This class provides an interface to collect and store tensor data associated with Kratos entities' property variables
 *          specified by @p pVariable in the Properties. It extends TensorAdaptor<double> and allows for flexible data management,
 *          including initialization, data collection, and storage operations.
 *
 * @section PropertiesVariableTensorAdaptor_supported_container Supported container types
 * - @ref ModelPart::ConditionsContainerType
 * - @ref ModelPart::ElementsContainerType
 *
 * @section PropertiesVariableTensorAdaptor_usage Usage
 * - Use @ref Check to verify that the variable exists and they are holding unique memory locations in the entities' properties before collecting/storing data.
 * - Use @ref CollectData to fill internal data from Kratos data structures entities' properties.
 * - Use @ref StoreData to write internal data back to the container's entities' properties, adding the variable if necessary.
 *
 * @author Suneth Warnakulasuriya
 * @see @ref TensorAdaptor                  Base class.
 * @see @ref DataValueContainer::GetValue   Variable value retrieval/update method.
 * @see @ref DataValueContainer::SetValue   Variable value set method.
 */
class KRATOS_API(OPTIMIZATION_APPLICATION) PropertiesVariableTensorAdaptor: public TensorAdaptor<double> {
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(PropertiesVariableTensorAdaptor);

    using BaseType = TensorAdaptor<double>;

    using VariablePointerType = TensorAdaptorUtils::VariablePointerType;

    ///@}
    ///@name Life cycle
    ///@{

    template<class TContainerPointerType>
    PropertiesVariableTensorAdaptor(
        TContainerPointerType pContainer,
        VariablePointerType pVariable);

    template<class TContainerPointerType>
    PropertiesVariableTensorAdaptor(
        TContainerPointerType pContainer,
        VariablePointerType pVariable,
        const std::vector<unsigned int>& rDataShape);

    PropertiesVariableTensorAdaptor(
        const TensorAdaptor& rOther,
        VariablePointerType pVariable,
        const bool Copy = true);

    // Destructor
    ~PropertiesVariableTensorAdaptor() override = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Execution of this check is only required if this PropertiesVariableTensorAdaptor
     *        will be used to call the @ref CollectData method. If it is only required to be called
     *        with the @ref StoreData, then please do not execute this @ref Check. Because
     *        then this Check will disallow storing data by throwing an error, when the entities do not have the variables
     *        already defined.
     * @see StoreData
     */
    void Check() const override;

    /**
     * @brief Fill the internal data from Kratos data structures's properties
     * @details This will fill the internal data from Kratos data structures's properties. It is advised to call
     *          at least once the Check method to ensure there won't be any errors if the
     *          variable is not present in the entities' properties. This will return Variable::Zero()
     *          values for all the entities when collecting if the variable is not set before.
     */
    void CollectData() override;

    /**
     * @brief Store internal data to the given container's properties.
     * @details This method is designed to store data even if the variable is not already available in the
     *          entities' properties. If it is not present in the entities' properties, then a correctly shaped zero valued value
     *          will be set and then their relevant components will be overwritten by this method.
     *
     * @warning If the entities' properties are not unique for each entity, this will result in a race condition.
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

    VariablePointerType mpVariable;

    ///@}
};

/// @}
} // namespace Kratos