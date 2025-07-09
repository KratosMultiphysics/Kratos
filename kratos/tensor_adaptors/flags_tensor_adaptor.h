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
#include "utilities/container_io_utils.h"

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) FlagsTensorAdaptor: public TensorAdaptor<bool> {
public:

    ///@name Type definitions
    ///@{

    using Pointer = Kratos::intrusive_ptr<FlagsTensorAdaptor>;

    using ConstPointer = Kratos::intrusive_ptr<const FlagsTensorAdaptor>;

    using BaseType = TensorAdaptor<bool>;

    ///@}
    ///@name Life cycle
    ///@{

    template<class TContainerPointerType>
    FlagsTensorAdaptor(
        TContainerPointerType pContainer,
        const Flags& rFlags);

    // Destructor
    ~FlagsTensorAdaptor() override = default;

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

    ContainerType GetContainer() const override;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    FlagsIO::Pointer mpIO;

    ContainerType mpContainer;

    ///@}
};

/// @}
} // namespace Kratos