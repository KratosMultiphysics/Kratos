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
#include "includes/model_part.h"
#include "tensor_adaptor.h"
#include "includes/global_variables.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @ingroup TensorAdaptors
 * @brief Adapts node position data from a Kratos ModelPart to a tensor representation.
 *
 * @details This class provides an adaptor for accessing and manipulating node position data
 *          as a tensor, enabling integration with tensor-based algorithms and operations.
 *
 *          The @ref Globals::Configuration is used to indicate whether this tensor adaptor will work
 *          with the initial (@ref Globals::Configuration::Initial) nodal positions or the current (
 *          @ref Globals::Configuration::Current) nodal positions.
 *
 * @section NodePositionTensorAdaptor_supported_container Supported container types
 * - @ref ModelPart::NodesContainerType
 *
 * @section NodePositionTensorAdaptor_usage Usage
 * - Construct with a pointer to an @ref ModelPart::NodesContainerType "array of nodes" and a configuration.
 * - Use @ref CollectData to fill internal tensor data from Kratos nodes' positions.
 * - Use @ref StoreData to write tensor data back to a Kratos nodes' positions.
 *
 * @author Suneth Warnakulasuriya
 * @see @ref TensorAdaptor                  Base class.
 * @see @ref Node::Coordinates              Method used to retrieve/update nodes' current coordinates.
 * @see @ref Node::GetInitialPosition       Method used to retrieve/update nodes' initial coordinates.
 * @see @ref Globals::Configuration         Enum used to specify whether to work on current or initial coordinates.
 */
class KRATOS_API(KRATOS_CORE) NodePositionTensorAdaptor: public TensorAdaptor<double> {
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(NodePositionTensorAdaptor);

    using BaseType = TensorAdaptor<double>;

    ///@}
    ///@name Life cycle
    ///@{

    NodePositionTensorAdaptor(
        ModelPart::NodesContainerType::Pointer pContainer,
        Globals::Configuration Configuration);

    NodePositionTensorAdaptor(
        ModelPart::NodesContainerType::Pointer pContainer,
        Globals::Configuration Configuration,
        const std::vector<unsigned int>& rDataShape);

    NodePositionTensorAdaptor(
        const TensorAdaptor& rOther,
        Globals::Configuration Configuration,
        const bool Copy = true);

    // Destructor
    ~NodePositionTensorAdaptor() override = default;

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

    Globals::Configuration mConfiguration;

    ///@}
};

/// @}
} // namespace Kratos