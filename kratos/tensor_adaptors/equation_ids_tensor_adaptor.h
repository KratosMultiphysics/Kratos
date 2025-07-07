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

class KRATOS_API(KRATOS_CORE) EquationIdsTensorAdaptor: public TensorAdaptor<int> {
public:

    ///@name Type definitions
    ///@{

    using Pointer = Kratos::intrusive_ptr<EquationIdsTensorAdaptor>;

    using ConstPointer = Kratos::intrusive_ptr<const EquationIdsTensorAdaptor>;

    using BaseType = TensorAdaptor<int>;

    ///@}
    ///@name Life cycle
    ///@{

    template<class TContainerPointerType>
    EquationIdsTensorAdaptor(
        TContainerPointerType pContainer,
        ProcessInfo::Pointer pProcessInfo);

    // Destructor
    ~EquationIdsTensorAdaptor() override = default;

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

    ProcessInfo::Pointer mpProcessInfo;

    ContainerType mpContainer;

    ///@}
};

/// @}
} // namespace Kratos