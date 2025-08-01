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

/**
 * @class EquationIdsTensorAdaptor
 * @brief Adaptor class for handling equation IDs as tensor data within Kratos.
 *
 * This class inherits from TensorAdaptor<int> and provides specialized functionality
 * for collecting and storing equation ID data from Kratos data conditions or elements. It is
 * designed to interface with containers and process information relevant to equation IDs.
 *
 * @tparam TContainerPointerType Type of the pointer to the container holding equation IDs.
 *
 * @author Kratos Team
 */
namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) EquationIdsTensorAdaptor: public TensorAdaptor<int> {
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(EquationIdsTensorAdaptor);

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
        const bool Copy = false);

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