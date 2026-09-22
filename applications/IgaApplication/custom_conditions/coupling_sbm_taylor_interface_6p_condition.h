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

/// Domain-level SBM mutual coupling, built directly on the TRUE shared
/// interface 
class CouplingSbmTaylorInterface6pCondition
    : public CouplingNitsche6pCondition
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(CouplingSbmTaylorInterface6pCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    CouplingSbmTaylorInterface6pCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : CouplingNitsche6pCondition(NewId, pGeometry)
    {};

    CouplingSbmTaylorInterface6pCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : CouplingNitsche6pCondition(NewId, pGeometry, pProperties)
    {};

    CouplingSbmTaylorInterface6pCondition()
        : CouplingNitsche6pCondition()
    {};

    virtual ~CouplingSbmTaylorInterface6pCondition() = default;

    ///@}
    ///@name Life Cycle
    ///@{

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive<CouplingSbmTaylorInterface6pCondition>(
            NewId, pGeom, pProperties);
    };

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive<CouplingSbmTaylorInterface6pCondition>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    };

    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief Assigns the two shift sources (master's and slave's own
    *        Gamma_tilde_h conditions) this condition's Taylor shift reads
    *        from.
    */
    void SetShiftSources(Condition::Pointer pMasterSource, Condition::Pointer pSlaveSource)
    {
        mpMasterShiftSource = pMasterSource;
        mpSlaveShiftSource = pSlaveSource;
    }

    /**
    * @brief Computes, once, the Taylor distance vector and order from each
    *        shift source's own evaluation point to this condition's own
    *        target point.
    */
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

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "\"CouplingSbmTaylorInterface6pCondition\" #" << Id();
        return buffer.str();
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"CouplingSbmTaylorInterface6pCondition\" #" << Id();
    }

    void PrintData(std::ostream& rOStream) const override {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}

protected:

    ///@name Protected member variables
    ///@{

    Condition::Pointer mpMasterShiftSource;
    Condition::Pointer mpSlaveShiftSource;
    Vector mDistanceVectorMaster;
    Vector mDistanceVectorSlave;
    SizeType mTaylorOrderMaster = 0;
    SizeType mTaylorOrderSlave = 0;

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

}; // Class CouplingSbmTaylorInterface6pCondition

}  // namespace Kratos.
