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
#include "includes/process_info.h"

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) GaussPointVariableTensorAdaptor: public TensorAdaptor<double> {
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GaussPointVariableTensorAdaptor);

    using BaseType = TensorAdaptor<double>;

    using VariablePointerType = TensorAdaptorUtils::VariablePointerType;

    ///@}
    ///@name Life cycle
    ///@{

    template<class TContainerPointerType>
    GaussPointVariableTensorAdaptor(
        TContainerPointerType pContainer,
        VariablePointerType pVariable,
        ProcessInfo::Pointer pProcessInfo);

    // Destructor
    ~GaussPointVariableTensorAdaptor() override = default;

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

    ProcessInfo::Pointer mpProcessInfo;

    ///@}
};

/// @}
} // namespace Kratos