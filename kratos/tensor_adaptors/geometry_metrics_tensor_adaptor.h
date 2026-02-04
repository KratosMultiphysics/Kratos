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
#include "tensor_adaptors/tensor_adaptor.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @class GeometryMetricsTensorAdaptor
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
class KRATOS_API(KRATOS_CORE) GeometryMetricsTensorAdaptor: public TensorAdaptor<double> {
public:

    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(GeometryMetricsTensorAdaptor);

    using BaseType = TensorAdaptor<double>;

    ///@}
    ///@name Enums
    ///@{

    enum DatumType
    {
        DomainSize
    };

    ///@}
    ///@name Life cycle
    ///@{

    template<class TContainerPointerType>
    GeometryMetricsTensorAdaptor(
        TContainerPointerType pContainer,
        const DatumType Datum);

    GeometryMetricsTensorAdaptor(
        const TensorAdaptor& rOther,
        const DatumType Datum,
        const bool Copy = true);

    // Destructor
    ~GeometryMetricsTensorAdaptor() override = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Fill the internal data with number of neightour entities of nodes
     * @details This will fill the internal data for each node with its available number of entities (i.e. conditions / elements).
     * @throws std::runtime_error if there is an entity having a node with id which is not present in the provided @p pNodes when this object is constructed.
     */
    void CollectData() override;

    /**
     * @brief Does not do anything
     * @throws std::runtime_error always.
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

    DatumType mDatum;

    ///@}
};

/// @}
} // namespace Kratos