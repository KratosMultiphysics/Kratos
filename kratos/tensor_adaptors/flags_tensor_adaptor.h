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
 * @class FlagsTensorAdaptor
 * @ingroup TensorAdaptors
 * @brief Adaptor class for handling tensor data associated with Kratos Flags.
 *
 * @details This class extends TensorAdaptor<bool> to provide specialized operations for
 *          tensors whose elements are boolean values, typically representing flags in Kratos.
 *          It manages the association between tensor data and Kratos Flags, allowing for
 *          efficient data collection and storage from/to Kratos data structures.
 *
 * @section supported_container Supported container types
 * - @ref ModelPart::NodesContainerType
 * - @ref ModelPart::ConditionsContainerType
 * - @ref ModelPart::ElementsContainerType
 *
 * @section Usage
 * - Use Check() to verify the tensor data shape compatibility.
 * - Use CollectData() to fill internal tensor data from flag values from nodes, conditions and elements.
 * - Use StoreData() to store back the internal tensor data to flag values of nodes, conditions and elements.
 *
 * @author Suneth Warnakulasuriya
 * @see @ref TensorAdaptor                 Base class.
 * @see @ref Flags::Is                     Method used to retrieve flag status from entities.
 * @see @ref Flags::Set                    Method used to set the flag status in the entities.
 */
class KRATOS_API(KRATOS_CORE) FlagsTensorAdaptor: public TensorAdaptor<bool> {
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(FlagsTensorAdaptor);

    using BaseType = TensorAdaptor<bool>;

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

    // Destructor
    ~FlagsTensorAdaptor() override = default;

    ///@}
    ///@name Public operations
    ///@{

    void Check() const override;

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

    Flags mFlags;

    ///@}
};

/// @}
} // namespace Kratos