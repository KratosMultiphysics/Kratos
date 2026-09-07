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
 * @ingroup TensorAdaptors
 * @brief Adaptor class for calculating geometry metrics.
 *
 * @details This class provides an interface to calculate geometry metrics for given entity container pointer @p pContainer .
 *          This @ref TensorAdaptor only implements the @ref CollectData method. Following geometry metrics are supported.
 *                  - @ref GeometryMetricsTensorAdaptor::DomainSize
 *
 * @section GeometryMetricsTensorAdaptor_supported_container Supported entity container types
 * - @ref ModelPart::GeometryContainerType
 * - @ref ModelPart::ConditionsContainerType
 * - @ref ModelPart::ElementsContainerType
 *
 * @section GeometryMetricsTensorAdaptor_usage Usage
 * - Use @ref CollectData to fill internal data with the metric from each entity given by the container pointer @p pContainer .
 *
 * @author Suneth Warnakulasuriya
 * @see @ref TensorAdaptor                  Base class.
 */
class KRATOS_API(KRATOS_CORE) GeometryMetricsTensorAdaptor: public TensorAdaptor<double> {
public:
    ///@name Enums
    ///@{

    enum class Metric
    {
        DomainSize
    };

    ///@}
    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(GeometryMetricsTensorAdaptor);

    using BaseType = TensorAdaptor<double>;

    using enum Metric;

    ///@}
    ///@name Life cycle
    ///@{

    template<class TContainerPointerType>
    GeometryMetricsTensorAdaptor(
        TContainerPointerType pContainer,
        const Metric Datum);

    GeometryMetricsTensorAdaptor(
        const TensorAdaptor& rOther,
        const Metric Datum,
        const bool Copy = true);

    // Destructor
    ~GeometryMetricsTensorAdaptor() override = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Clones the existing tensor adaptor.
     */
    TensorAdaptor::Pointer Clone() const override;

    /**
     * @brief Fill the internal data with metric of the entities in the container.
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

    Metric mMetric;

    ///@}
};

/// @}
} // namespace Kratos