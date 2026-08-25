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
#include "includes/process_info.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @ingroup TensorAdaptors
 * @brief Adaptor class for handling equation IDs as tensor data within Kratos.
 *
 * @details This class inherits from TensorAdaptor<int> and provides specialized functionality
 *          for collecting and storing equation ID data from Kratos data conditions or elements. It is
 *          designed to interface with containers and process information relevant to equation IDs.
 *
 * @warning StoreData method is throwing an unimplemented error since setting equation ids is not allowed.
 *
 * @section EquationIdsTensorAdaptor_supported_container Supported container types
 * - @ref ModelPart::ConditionsContainerType
 * - @ref ModelPart::ElementsContainerType
 *
 * @section EquationIdsTensorAdaptor_usage Usage
 * - Use @ref Check to verify the tensor data shape compatibility.
 * - Use @ref CollectData to fill internal tensor data from equation ids read from elements or conditions.
 * - Use @ref StoreData Not possible.
 *
 * @author Suneth Warnakulasuriya
 * @see @ref TensorAdaptor                 Base class.
 * @see @ref Element::EquationIdVector     Equation id retrieval method from elements.
 * @see @ref Condition::EquationIdVector   Equation id retrieval method from conditions.
 */
class KRATOS_API(KRATOS_CORE) EquationIdsTensorAdaptor: public TensorAdaptor<int> {
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(EquationIdsTensorAdaptor);

    using BaseType = TensorAdaptor<int>;

    ///@}
    ///@name Life cycle
    ///@{

    template<class TContainerPointerType>
    EquationIdsTensorAdaptor(
        TContainerPointerType pContainer,
        ProcessInfo::Pointer pProcessInfo);

    EquationIdsTensorAdaptor(
        const TensorAdaptor& rOther,
        ProcessInfo::Pointer pProcessInfo,
        const bool Copy = true);

    EquationIdsTensorAdaptor(const EquationIdsTensorAdaptor& rOther) = default;

    // Destructor
    ~EquationIdsTensorAdaptor() override = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Clones the existing tensor adaptor.
     */
    TensorAdaptor::Pointer Clone() const override;

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

    ProcessInfo::Pointer mpProcessInfo;

    ///@}
};

/// @}
} // namespace Kratos