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
#include "containers/flags.h"
#include "tensor_adaptor.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @ingroup TensorAdaptors
 * @brief Adaptor class for handling tensor data associated with Kratos Flags.
 *
 * @details This class extends TensorAdaptor<int> to provide specialized operations for
 *          tensors whose elements are boolean values, typically representing flags in Kratos.
 *          It manages the association between tensor data and Kratos Flags, allowing for
 *          efficient data collection and storage from/to Kratos data structures.
 *
 *          The state of the flags are stored as int having following definitions for each int value:
 *              - @p -1 - The flag is not defined.
 *              - @p 0 - The flag is defined and the value is false.
 *              - @p 1 - The flag is defined and the value is true.
 *
 * @section FlagsTensorAdaptor_supported_container Supported container types
 * - @ref ModelPart::NodesContainerType
 * - @ref ModelPart::ConditionsContainerType
 * - @ref ModelPart::ElementsContainerType
 *
 * @section FlagsTensorAdaptor_usage Usage
 * - Use @ref Check to verify the tensor data shape compatibility.
 * - Use @ref CollectData to fill internal tensor data from flag values from nodes, conditions and elements.
 * - Use @ref StoreData to store back the internal tensor data to flag values of nodes, conditions and elements.
 *
 * @author Suneth Warnakulasuriya
 * @see @ref TensorAdaptor                 Base class.
 * @see @ref Flags::Is                     Method used to retrieve flag status from entities.
 * @see @ref Flags::IsDefined              Method used to retrieve flag defined status from entities.
 * @see @ref Flags::Set                    Method used to set the flag status in the entities.
 * @see @ref Flags::Reset                  Method used to rest the flag in entities.
 */
class KRATOS_API(KRATOS_CORE) FlagsTensorAdaptor: public TensorAdaptor<int> {
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(FlagsTensorAdaptor);

    using BaseType = TensorAdaptor<int>;

    ///@}
    ///@name Life cycle
    ///@{

    template<class TContainerPointerType>
    FlagsTensorAdaptor(
        TContainerPointerType pContainer,
        const Flags& rFlags);

    FlagsTensorAdaptor(
        const TensorAdaptor& rOther,
        const Flags& rFlags,
        const bool Copy = true);

    FlagsTensorAdaptor(const FlagsTensorAdaptor& rOther) = default;

    // Destructor
    ~FlagsTensorAdaptor() override = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Clones the existing tensor adaptor.
     */
    TensorAdaptor::Pointer Clone() const override;

    void Check() const override;

    /**
     * @brief Fill the internal data from Kratos data structures
     */
    void CollectData() override;

    /**
     * @brief Store internal data to the given TContainerType container.
     * @throws std::runtime_error if the values in the Tensor adaptor's internal data does not correspond to any of the @p -1,  @p 0 or @p 1 states
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

    Flags mFlags;

    ///@}
};

/// @}
} // namespace Kratos