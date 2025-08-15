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

/**
 * @ingroup TensorAdaptors
 * @brief Adaptor class for handling tensor data associated with gauss point values.
 *
 * @details This class provides an adaptor to retrieve Kratos variables computed at Gauss points
 *          as tensors. It inherits from TensorAdaptor<double> and enables tensor operations on variables
 *          associated with finite element integration points.
 *
 * @warning StoreData method is throwing an unimplemented error since setting gauss points is not allowed.
 *
 * @section GaussPointVariableTensorAdaptor_supported_container Supported container types
 * - @ref ModelPart::ConditionsContainerType
 * - @ref ModelPart::ElementsContainerType
 *
 * @section GaussPointVariableTensorAdaptor_usage Usage
 * - Use @ref Check to verify the tensor data shape compatibility.
 * - Use @ref CollectData to fill internal tensor data from gauss point computed values at each element or condition.
 * - Use @ref StoreData    Not possible.
 *
 * @author Suneth Warnakulasuriya
 * @see @ref TensorAdaptor                                Base class.
 * @see @ref Element::CalculateOnIntegrationPoints        Used to retrieve gauss point values from elements.
 * @see @ref Condition::CalculateOnIntegrationPoints      Used to retrieve gauss point values from elements.
 */
class KRATOS_API(KRATOS_CORE) GaussPointVariableTensorAdaptor: public TensorAdaptor<double> {
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(GaussPointVariableTensorAdaptor);

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

    GaussPointVariableTensorAdaptor(
        const TensorAdaptor& rOther,
        VariablePointerType pVariable,
        ProcessInfo::Pointer pProcessInfo,
        const bool Copy = true);

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