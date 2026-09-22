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

/// Shifted Boundary Method (SBM) coupling condition for the degenerated-solid
/// Reissner-Mindlin shell_6p_element, using a Moving-Least-Squares (MLS)
class CouplingSbmExtensionOperator6pCondition
    : public CouplingNitsche6pCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// pointer of CouplingSbmExtensionOperator6pCondition
    KRATOS_CLASS_POINTER_DEFINITION(CouplingSbmExtensionOperator6pCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Id and geometry
    CouplingSbmExtensionOperator6pCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : CouplingNitsche6pCondition(NewId, pGeometry)
    {};

    /// Constructor with Id, geometry and property
    CouplingSbmExtensionOperator6pCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : CouplingNitsche6pCondition(NewId, pGeometry, pProperties)
    {};

    /// Default constructor
    CouplingSbmExtensionOperator6pCondition()
        : CouplingNitsche6pCondition()
    {};

    /// Destructor.
    virtual ~CouplingSbmExtensionOperator6pCondition() = default;

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
        return Kratos::make_intrusive<CouplingSbmExtensionOperator6pCondition>(
            NewId, pGeom, pProperties);
    };

    /// Create with Id, pointer to geometry and pointer to property
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive< CouplingSbmExtensionOperator6pCondition >(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    };

    ///@}
    ///@name Operations
    ///@{

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

    /**
    * @brief Sets on rResult the ID's of the condition degrees of freedom.
    */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    /**
    * @brief Sets on rConditionDofList the degrees of freedom of the FULL
    *        (own + MLS cloud) DOF set
    */
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

    /**
    * @brief Assigns the two shift sources (master's and slave's own
    *        Gamma_tilde_h conditions) the true boundary adjoint term reads 
    */
    void SetShiftSources(Condition::Pointer pMasterSource, Condition::Pointer pSlaveSource)
    {
        mpMasterShiftSource = pMasterSource;
        mpSlaveShiftSource = pSlaveSource;
    }

    ///@}
    ///@name Check
    ///@{

    /// Performs check if NITSCHE_STABILIZATION_FACTOR is provided and if
    /// PrecomputeAndStoreCouplingExtensionOperators has already been run.
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "\"CouplingSbmExtensionOperator6pCondition\" #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"CouplingSbmExtensionOperator6pCondition\" #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}

protected:

    ///@name Protected member variables
    ///@{

    Condition::Pointer mpMasterShiftSource;
    Condition::Pointer mpSlaveShiftSource;

    ///@}

private:

    ///@name Private operations
    ///@{

    /// Number of DOF-block nodes in the stored SBM_MLS_ALL_DOF_NODES list.
    SizeType GetFullDofNodeCount() const;

    ///@}
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

}; // Class CouplingSbmExtensionOperator6pCondition

}  // namespace Kratos.
