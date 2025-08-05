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

/**
 * @class NodePositionTensorAdaptor
 * @brief Adapts node position data from a Kratos ModelPart to a tensor representation.
 *
 * This class provides an adaptor for accessing and manipulating node position data
 * as a tensor, enabling integration with tensor-based algorithms and operations.
 * It inherits from TensorAdaptor<double> and manages the mapping between Kratos
 * node containers and tensor data structures.
 *
 * @tparam double The data type for tensor elements.
 *
 * @section Usage
 * - Construct with a pointer to a ModelPart::NodesContainerType and a configuration.
 * - Optionally specify the shape of the data or copy from another TensorAdaptor.
 * - Use CollectData() to fill internal tensor data from Kratos nodes' positions.
 * - Use StoreData() to write tensor data back to a Kratos nodes' positions.
 */
namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) NodePositionTensorAdaptor: public TensorAdaptor<double> {
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(NodePositionTensorAdaptor);

    using BaseType = TensorAdaptor<double>;

    ///@}
    ///@name Life cycle
    ///@

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