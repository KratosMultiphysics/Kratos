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

/**
 * @class GaussPointVariableTensorAdaptor
 * @brief Adapts Kratos variables at Gauss points to tensor operations.
 *
 * This class provides an adaptor to retrieve Kratos variables computed at Gauss points
 * as tensors. It inherits from TensorAdaptor<double> and enables tensor operations on variables
 * associated with finite element integration points.
 *
 * @warning GaussPointVariableTensorAdaptor::StoreData method is throwing an unimplemented error
 *          since setting gauss points is not allowed.
 *
 * @section Usage
 * - Use GaussPointVariableTensorAdaptor::Check() to verify the tensor data shape compatibility.
 * - Use CollectData() to fill internal tensor data from gauss point computed values at each element or condition.
 * - Use StoreData()    Not used.
 */
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