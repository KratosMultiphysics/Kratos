//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mario Caballero
//

#pragma once

// System includes

// External includes

// Project includes
#include "custom_conditions/gsbm_laplacian_condition_dirichlet.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(IGA_APPLICATION) GsbmLaplacianConditionDirichletStabilized
    : public GsbmLaplacianConditionDirichlet
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GsbmLaplacianConditionDirichletStabilized);

    using BaseType = GsbmLaplacianConditionDirichlet;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Id and geometry.
    GsbmLaplacianConditionDirichletStabilized(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /// Constructor with Id, geometry and properties.
    GsbmLaplacianConditionDirichletStabilized(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /// Default constructor.
    GsbmLaplacianConditionDirichletStabilized()
        : BaseType()
    {
    }

    /// Destructor.
    ~GsbmLaplacianConditionDirichletStabilized() override = default;

    /// Create with Id, pointer to geometry and pointer to properties
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GsbmLaplacianConditionDirichletStabilized>(
            NewId, pGeom, pProperties
        );
    }

    /// Create with Id, nodes and pointer to properties
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GsbmLaplacianConditionDirichletStabilized>(
            NewId, GetGeometry().Create(ThisNodes), pProperties
        );
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "GsbmLaplacianConditionDirichletStabilized #" + std::to_string(Id());
    }

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    /// Shifted shape functions: series of N from x_\Gamma, plus \rho^m(N)
    Vector CalculateShiftedShapeFunctions(
        const GeometryType& rSide,
        const Node& rProjectionNode,
        const array_1d<double, 3>& rDistance,
        const SizeType Order) const override;

    /// G = g - \sum_k (t^k / k!) D_k f, from VOLUME_SOURCE_DERIVATIVES
    double CalculateShiftedBoundaryValue(
        const Node& rProjectionNode,
        const array_1d<double, 3>& rDistance,
        const SizeType Order,
        const Variable<double>& rUnknownVariable) const override;

    ///@}

private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
    }

    ///@}
};

///@}

} // namespace Kratos
