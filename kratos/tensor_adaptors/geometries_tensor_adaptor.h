//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#pragma once

// System includes
#include <string>

// External includes

// Project includes
#include "geometries/geometry_data.h"
#include "tensor_adaptor.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class GeometriesTensorAdaptor
 * @ingroup TensorAdaptors
 * @brief Adaptor class for handling Geometry functions (ShapeFunctions,
 * Derivatives, Jacobians).
 *
 * @details This class provides an interface to evaluate and expose geometry
 * functions on Gauss points for various Kratos containers (Geometries,
 * Elements, Conditions). It extends TensorAdaptor<double>. Possible DatumTypes:
 *          - ShapeFunctions: Returns shape function values. Shape: [n_elem,
 * n_gauss, n_node].
 *          - ShapeFunctionDerivatives: Returns global derivatives (DN/DX).
 * Shape: [n_elem, n_gauss, n_node, working_dim]. This means that
 * DNDx(iel,igauss,inode,0) is the dN(inode)/dx
 *          - Jacobians: Returns Jacobians at gauss points. Shape: [n_elem,
 * n_gauss, working_dim, local_dim].
 *          - IntegrationWeights: Returns integration weights. Shape: [n_elem,
 * n_gauss].
 *
 * @author Riccardo Rossi and Antigravity AI
 */
class KRATOS_API(KRATOS_CORE) GeometriesTensorAdaptor
    : public TensorAdaptor<double>
{
public:
    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(GeometriesTensorAdaptor);

    using BaseType = TensorAdaptor<double>;

    enum class DatumType
    {
        ShapeFunctions,
        ShapeFunctionDerivatives,
        Jacobians,
        IntegrationWeights
    };

    ///@}
    ///@name Life cycle
    ///@{

    GeometriesTensorAdaptor(
        ContainerPointerType pContainer, DatumType Datum,
        GeometryData::IntegrationMethod IntegrationMethod =
            GeometryData::IntegrationMethod::NumberOfIntegrationMethods);

    GeometriesTensorAdaptor(const TensorAdaptor &rOther,
                            ContainerPointerType pContainer, DatumType Datum,
                            GeometryData::IntegrationMethod IntegrationMethod,
                            const bool Copy = true);

    // Destructor
    ~GeometriesTensorAdaptor() override = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Check if the container is valid.
     */
    void Check() const override;

    /**
     * @brief Fill the internal data from Kratos data structures.
     */
    void CollectData() override;

    /**
     * @brief Store internal data to the given container.
     */
    void StoreData() override;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override;

    ///@}

private:
    // Helper templates scoped to the class to avoid unity-build collisions
    template <class TContainerType>
    static constexpr bool IsSupportedContainer();

    template <class TEntity>
    static const auto& GetGeometry(const TEntity& rEntity);

    DatumType mDatum;
    GeometryData::IntegrationMethod mIntegrationMethod;
};

/// @}
} // namespace Kratos
