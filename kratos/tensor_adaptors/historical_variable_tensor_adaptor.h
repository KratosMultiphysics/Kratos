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

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @ingroup TensorAdaptors
 * @brief Adaptor for accessing and manipulating historical variable data of nodes as tensors.
 *
 * @details This class provides an interface to access and modify historical variable data specified in @p pVariable
 *          stored in the nodes of a @ref ModelPart for a given @p StepIndex,
 *          exposing them as tensor data structures. It inherits from TensorAdaptor<double> and allows for efficient data
 *          collection and storage operations, ensuring compatibility with Kratos' historical data containers.
 *
 * @warning @ref CollectData and @ref StoreData may cause segmentation faults if the variable (i.e. @p pVariable) is not found in the historical
 *          data value container of a node or the buffer size is not adequate for the specified @p StepIndex.
 *          Always use @ref Check before performing these operations.
 *
 * @section HistoricalVariableTensorAdaptor_supported_container Supported container types
 * - @ref ModelPart::NodesContainerType
 *
 * @section HistoricalVariableTensorAdaptor_usage Usage
 * - Use @ref Check to verify variable availability and buffer size before data operations.
 * - Use @ref CollectData to fill internal tensor data from Kratos historical nodal data containers.
 * - Use @ref StoreData to write internal tensor data back to Kratos historical nodal data containers.
 *
 * @author Suneth Warnakulasuriya
 * @see @ref TensorAdaptor                  Base class.
 * @see @ref Node::FastGetSolutionStepValue Method used to retrieve data from nodal historical data value container.
 */
class KRATOS_API(KRATOS_CORE) HistoricalVariableTensorAdaptor: public TensorAdaptor<double> {
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(HistoricalVariableTensorAdaptor);

    using BaseType = TensorAdaptor<double>;

    using VariablePointerType = TensorAdaptorUtils::VariablePointerType;

    ///@}
    ///@name Life cycle
    ///@{

    HistoricalVariableTensorAdaptor(
        ModelPart::NodesContainerType::Pointer pContainer,
        VariablePointerType pVariable,
        const int StepIndex = 0);

    HistoricalVariableTensorAdaptor(
        ModelPart::NodesContainerType::Pointer pContainer,
        VariablePointerType pVariable,
        const std::vector<unsigned int>& rDataShape,
        const int StepIndex = 0);

    HistoricalVariableTensorAdaptor(
        const TensorAdaptor& rOther,
        VariablePointerType pVariable,
        const int StepIndex = 0,
        const bool Copy = true);

    // Destructor
    ~HistoricalVariableTensorAdaptor() override = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Checks if the given variable is available and the buffer size is enough in each node.
     */
    void Check() const override;

    /**
     * @brief Fill the internal data from Kratos data structures
     * @details This method will fill the internal data from the historical data value container of each node.
     * @warning This may SEGFAULT if the variable is not found in the historical data value container of a node.
     *          Therefore, it is suggested to use HistoricalVariableTensorAdaptor::Check method before.
     */
    void CollectData() override;

    /**
     * @details This method will fill the historical data value container data of each node from the internal data.
     * @warning This may SEGFAULT if the variable is not found in the historical data value container of a node.
     *          Therefore, it is suggested to use HistoricalVariableTensorAdaptor::Check method before.
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

    const int mStepIndex;

    ///@}
};

/// @}
} // namespace Kratos