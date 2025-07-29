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
        TensorData<double>::Pointer pTensorData,
        VariablePointerType pVariable,
        const int StepIndex = 0);


    // Destructor
    ~HistoricalVariableTensorAdaptor() override = default;

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

    VariablePointerType mpVariable;

    const int mStepIndex;

    ///@}
};

/// @}
} // namespace Kratos