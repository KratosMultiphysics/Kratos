//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ricky Aristio
//

#pragma once

// System includes
#include "includes/define.h"

// External includes

// Project includes
#include "custom_conditions/coupling_nitsche_6p_condition.h"

namespace Kratos
{

/// Surrogate boundary consistency condition for domain-level
/// SBM coupling of the degenerated-solid Reissner-Mindlin shell_6p_element.
class SbmSurrogateConsistency6pCondition
    : public CouplingNitsche6pCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of SbmSurrogateConsistency6pCondition
    KRATOS_CLASS_POINTER_DEFINITION(SbmSurrogateConsistency6pCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Id and geometry
    SbmSurrogateConsistency6pCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : CouplingNitsche6pCondition(NewId, pGeometry)
    {};

    /// Constructor with Id, geometry and property
    SbmSurrogateConsistency6pCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : CouplingNitsche6pCondition(NewId, pGeometry, pProperties)
    {};

    /// Default constructor
    SbmSurrogateConsistency6pCondition()
        : CouplingNitsche6pCondition()
    {};

    /// Destructor.
    virtual ~SbmSurrogateConsistency6pCondition() = default;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Create with Id, pointer to geometry and pointer to property
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive<SbmSurrogateConsistency6pCondition>(
            NewId, pGeom, pProperties);
    };

    /// Create with Id, pointer to geometry and pointer to property
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive< SbmSurrogateConsistency6pCondition >(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    };

    ///@}
    ///@name Operations
    ///@{

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        MatrixType left_hand_side_matrix = Matrix(0, 0);

        CalculateAll(left_hand_side_matrix, rRightHandSideVector,
            rCurrentProcessInfo, false, true);
    }

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        VectorType right_hand_side_vector = Vector(0);

        CalculateAll(rLeftHandSideMatrix, right_hand_side_vector,
            rCurrentProcessInfo, true, false);
    }

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
            rCurrentProcessInfo, true, true);
    }

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    /// Calculates left (K) and right (u) hand sides, according to the flags
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    );

    ///@}
    ///@name Check
    ///@{

    /// Performs check if THICKNESS/YOUNG_MODULUS/POISSON_RATIO are provided.
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "\"SbmSurrogateConsistency6pCondition\" #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"SbmSurrogateConsistency6pCondition\" #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}

private:

    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, CouplingNitsche6pCondition);
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, CouplingNitsche6pCondition);
    }

    ///@}

}; // Class SbmSurrogateConsistency6pCondition

}  // namespace Kratos.
