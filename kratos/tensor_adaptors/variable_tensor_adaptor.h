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
#include "tensor_adaptor_utils.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @class VariableTensorAdaptor
 * @ingroup TensorAdaptors
 * @brief Adaptor class for handling values represented as variables in Kratos entities' DataValueContainers.
 *
 * @details This class provides an interface to collect and store tensor data associated with Kratos entities' variables
 *          specified by @p pVariable in the DataValueContainer (i.e. non-historical container). It extends TensorAdaptor<double> and allows for flexible data management,
 *          including initialization, data collection, and storage operations.
 *
 * @warning When using StoreData() to add variables to entities (if they are not present already in the entities), avoid calling Check() to prevent raising an error which
 *          will say the variable is not found.
 *
 * @section VariableTensorAdaptor_supported_container Supported container types
 * - @ref ModelPart::NodesContainerType
 * - @ref ModelPart::ConditionsContainerType
 * - @ref ModelPart::ElementsContainerType
 * - @ref ModelPart::PropertiesContainerType
 * - @ref ModelPart::GeometryContainerType
 * - @ref ModelPart::MasterSlaveConstraintContainerType
 *
 * @section VariableTensorAdaptor_usage Usage
 * - Use @ref Check to verify that the variable exists in the entities before collecting data.
 * - Use @ref CollectData to fill internal data from Kratos data structures.
 * - Use @ref StoreData to write internal data back to the container, adding the variable if necessary.
 *
 * @author Suneth Warnakulasuriya
 * @see @ref TensorAdaptor                  Base class.
 * @see @ref DataValueContainer::GetValue   Variable value retrieval/update method.
 * @see @ref DataValueContainer::SetValue   Variable value set method.
 */
class KRATOS_API(KRATOS_CORE) VariableTensorAdaptor: public TensorAdaptor<double> {
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(VariableTensorAdaptor);

    using BaseType = TensorAdaptor<double>;

    using VariablePointerType = TensorAdaptorUtils::VariablePointerType;

    ///@}
    ///@name Life cycle
    ///@{

    template<class TContainerPointerType>
    VariableTensorAdaptor(
        TContainerPointerType pContainer,
        VariablePointerType pVariable);

    template<class TContainerPointerType>
    VariableTensorAdaptor(
        TContainerPointerType pContainer,
        VariablePointerType pVariable,
        const std::vector<unsigned int>& rDataShape);

    VariableTensorAdaptor(
        const TensorAdaptor& rOther,
        VariablePointerType pVariable,
        const bool Copy = true);

    // Destructor
    ~VariableTensorAdaptor() override = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Execution of this check is only required if this VariableTensorAdaptor
     *        will be used to call the @ref CollectData method. If it is only required to be called
     *        with the @ref StoreData, then please do not execute this @ref Check. Because
     *        then this Check will disallow storing data by throwing an error, when the entities do not have the variables
     *        already defined.
     * @see StoreData
     */
    void Check() const override;

    /**
     * @brief Fill the internal data from Kratos data structures
     * @details This will fill the internal data from Kratos data structures. It is advised to call
     *          at least once the Check method to ensure there won't be any errors if the
     *          variable is not present in the entities. This will return Variable::Zero()
     *          values for all the entities when collecting if the variable is not set before.
     */
    void CollectData() override;

    /**
     * @brief Store internal data to the given container.
     * @details This method is designed to store data even if the variable is not already available in the
     *          entities. If it is not present in the entities, then a correctly shaped zero valued value
     *          will be set and then their relevant components will be overwritten by this method.
     * @warning Do not call @ref Check if you intend to use VariableTensorAdaptor to add the variable
     *          to the entities, and avoid unnecessarily initializing dummy values.
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