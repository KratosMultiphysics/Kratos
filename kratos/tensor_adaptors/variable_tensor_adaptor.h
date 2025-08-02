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

/**
 * @class VariableTensorAdaptor
 * @brief Adaptor class for handling values represented as variables in Kratos entities' DataValueContainers.
 *
 * @details This class provides an interface to collect and store tensor data associated with Kratos entities,
 *          using a specified variable from the DataValueContainer (i.e. non-historical container). It extends TensorAdaptor<double> and allows for flexible data management,
 *          including initialization, data collection, and storage operations.
 *
 * @tparam double The underlying data type for tensor values.
 *
 * @section Usage
 * - Use Check() to verify that the variable exists in the entities before collecting data.
 * - Use CollectData() to fill internal data from Kratos data structures.
 * - Use StoreData() to write internal data back to the container, adding the variable if necessary.
 *
 * @warning When using StoreData() to add variables to entities, avoid calling Check() to prevent unnecessary initialization.
 */
namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) VariableTensorAdaptor: public TensorAdaptor<double> {
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(VariableTensorAdaptor);

    using BaseType = TensorAdaptor<double>;

    using ContainerPointerType = Storage::ContainerPointerType;

    using VariablePointerType = TensorAdaptorUtils::VariablePointerType;

    ///@}
    ///@name Life cycle
    ///@

    VariableTensorAdaptor(
        ContainerPointerType pContainer,
        VariablePointerType pVariable);

    VariableTensorAdaptor(
        ContainerPointerType pContainer,
        VariablePointerType pVariable,
        const std::vector<unsigned int>& rDataShape);

    VariableTensorAdaptor(
        const TensorAdaptor& rOther,
        VariablePointerType pVariable,
        const bool Copy = false);

    // Destructor
    ~VariableTensorAdaptor() override = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief This check only suppose to be done if this VariableTensorAdaptor
     *        is suppose to use with the CollectData method. If it is supposed to be used
     *        with the StoreData, then this Check will disallow storing data, when the entities
     *        does not have the variables.
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
     * @brief Store internal data to the given TContainerType container.
     * @details This method is designed to store data even if the variable is not already available in the
     *          entities. If it is not present in the entities, then a correctly shaped zero valued values
     *          will be set and then their relevant components will be overwritten by this method by the
     *          values from the underlying data from the flat vector.
     * @warning Please don't call the VariableTensorAdaptor::Check method if you intend to use VariableTensorAdaptor to add the variable
     *          to the entities, and avoid unnecessary dummy initialization of values.
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