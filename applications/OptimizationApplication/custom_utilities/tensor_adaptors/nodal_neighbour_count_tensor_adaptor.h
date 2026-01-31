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
#include "tensor_adaptors/tensor_adaptor.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @class NodalNeighbourCountTensorAdaptor
 * @ingroup TensorAdaptors
 * @brief Adaptor class for calculating neighbour entities for nodes.
 *
 * @details This class provides an interface to calculate neighbour entities ( @p pEntityContainer i.e. conditions / elements) for a given specific
 *          @p pNodes. tensor data associated with Kratos entities' property variables. This @ref TensorAdaptor only implements
 *          the @ref CollectData method.
 *
 * @section NodalNeighbourCountTensorAdaptor_supported_container Supported entity container types
 * - @ref ModelPart::ConditionsContainerType
 * - @ref ModelPart::ElementsContainerType
 *
 * @section NodalNeighbourCountTensorAdaptor_usage Usage
 * - Use @ref Check to verify that the variable exists and they are holding unique memory locations in the entities' properties before collecting/storing data.
 * - Use @ref CollectData to fill internal data with the number of neighbouring entities.
 *
 * @author Suneth Warnakulasuriya
 * @see @ref TensorAdaptor                  Base class.
 */
class KRATOS_API(OPTIMIZATION_APPLICATION) NodalNeighbourCountTensorAdaptor: public TensorAdaptor<int> {
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(NodalNeighbourCountTensorAdaptor);

    using BaseType = TensorAdaptor<int>;

    ///@}
    ///@name Life cycle
    ///@{

    template<class TContainerPointerType>
    NodalNeighbourCountTensorAdaptor(
        ModelPart::NodesContainerType::Pointer pNodes,
        TContainerPointerType pEntityContainer);

    template<class TContainerPointerType>
    NodalNeighbourCountTensorAdaptor(
        const TensorAdaptor& rOther,
        TContainerPointerType pEntityContainer,
        const bool Copy = true);

    // Destructor
    ~NodalNeighbourCountTensorAdaptor() override = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Fill the internal data from Kratos data structures's properties
     * @details This will fill the internal data from Kratos data structures's properties. It is advised to call
     *          at least once the Check method to ensure there won't be any errors if the
     *          variable is not present in the entities' properties. This will return Variable::Zero()
     *          values for all the entities when collecting if the variable is not set before.
     */
    void CollectData() override;

    /**
     * @brief Store internal data to the given container's properties.
     * @details This method is designed to store data even if the variable is not already available in the
     *          entities' properties. If it is not present in the entities' properties, then a correctly shaped zero valued value
     *          will be set and then their relevant components will be overwritten by this method.
     *
     * @warning If the entities' properties are not unique for each entity, this will result in a race condition.
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

    std::variant<ModelPart::ConditionsContainerType::Pointer, ModelPart::ElementsContainerType::Pointer> mpEntityContainer;

    ///@}
};

/// @}
} // namespace Kratos