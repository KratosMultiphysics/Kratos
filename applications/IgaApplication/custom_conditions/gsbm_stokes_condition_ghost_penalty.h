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
#include "includes/condition.h"
#include "includes/define.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(IGA_APPLICATION) GsbmStokesConditionGhostPenalty
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GsbmStokesConditionGhostPenalty);

    using SizeType = std::size_t;
    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Id and geometry.
    GsbmStokesConditionGhostPenalty(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {
    }

    /// Constructor with Id, geometry and properties.
    GsbmStokesConditionGhostPenalty(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    /// Default constructor.
    GsbmStokesConditionGhostPenalty()
        : Condition()
    {
    }

    /// Destructor.
    ~GsbmStokesConditionGhostPenalty() override = default;

    /// Create with Id, pointer to geometry and pointer to properties
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GsbmStokesConditionGhostPenalty>(
            NewId, pGeom, pProperties
        );
    }

    /// Create with Id, nodes and pointer to properties
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GsbmStokesConditionGhostPenalty>(
            NewId, GetGeometry().Create(ThisNodes), pProperties
        );
    }

    ///@}
    ///@name Operations
    ///@{

    /// Ghost penalty of the GSBM Stokes problem on a facet, in residual form:
    /// of each velocity component, and of the pressure s(q, p), or the mixed
    /// +s(q, div u) and -s(div v, p) with PRESSURE_GHOST_PENALTY_ON_DIVERGENCE.
    /// Each term penalizes the jump of D_n^k, with k the degree of the space
    /// across the facet, and is normalized by (k!)^2 (2k + 1)
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /// Left hand side of CalculateLocalSystem
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /// Right hand side of CalculateLocalSystem
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /// Equation Ids of u_x, u_y and p of the - side, then of the + side
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /// Dofs of u_x, u_y and p of the - side, then of the + side
    void GetDofList(
        DofsVectorType& rConditionDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Check
    ///@{

    /// Checks the geometry parts and the properties, and that the spaces are
    /// divergence-conforming with PRESSURE_GHOST_PENALTY_ON_DIVERGENCE
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "GsbmStokesConditionGhostPenalty #" + std::to_string(Id());
    }

    ///@}

private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
    }

    ///@}
};

///@}

} // namespace Kratos
