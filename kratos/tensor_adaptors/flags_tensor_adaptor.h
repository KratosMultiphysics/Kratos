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
        const bool Copy = false);

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