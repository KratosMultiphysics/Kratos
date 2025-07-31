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

class KRATOS_API(KRATOS_CORE) HistoricalVariableTensorAdaptor: public TensorAdaptor<double> {
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(HistoricalVariableTensorAdaptor);

    using BaseType = TensorAdaptor<double>;

    using VariablePointerType = TensorAdaptorUtils::VariablePointerType;

    ///@}
    ///@name Life cycle
    ///@

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
        const bool Copy = false);


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