//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Crescenzio
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "includes/serializer.h"
#include "includes/process_info.h"

// Application includes
#include "particle_mechanics_application_variables.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Base class for all Conditions.

/**
 * @brief Implements the Navier Slip boundary condition
 * @tparam TDim Number of dimensions
 * @tparam TNumNodes Number of nodes
 */
template<unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(PARTICLE_MECHANICS_APPLICATION) MPMGridNavierSlipFrictionCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MPMGridNavierSlipFrictionCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(MPMGridNavierSlipFrictionCondition);

    static constexpr std::size_t BlockSize = TDim;
    static constexpr std::size_t LocalSize = TNumNodes*BlockSize;

    using SizeType = std::size_t;

    using IndexType = std::size_t;

    using VectorType = Vector;

    using MatrixType = Matrix;

    using PropertiesType = Properties;

    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;

    using NodesArrayType = Geometry<NodeType>::PointsArrayType;

    using EquationIdVectorType = std::vector<std::size_t>;

    using DofsVectorType = std::vector<Dof<double>::Pointer>;

    using IntegrationMethod = GeometryData::IntegrationMethod;

    using IntegrationPointsArrayType = GeometryType::IntegrationPointsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /**
     @param NewId Index of the new condition
     */
    MPMGridNavierSlipFrictionCondition(
        IndexType NewId = 0
        ):Condition(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    MPMGridNavierSlipFrictionCondition(
        IndexType NewId,
        const NodesArrayType& ThisNodes
        ):Condition(NewId,ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    MPMGridNavierSlipFrictionCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        ):Condition(NewId,pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    MPMGridNavierSlipFrictionCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        ):Condition(NewId,pGeometry,pProperties)
    {
    }

    /// Copy constructor.
    MPMGridNavierSlipFrictionCondition(
        MPMGridNavierSlipFrictionCondition const& rOther
        ):Condition(rOther)
    {
    }

    /// Destructor.
    ~MPMGridNavierSlipFrictionCondition() override {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    MPMGridNavierSlipFrictionCondition & operator=(MPMGridNavierSlipFrictionCondition const& rOther)
    {
        Condition::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create a new MPMGridNavierSlipFrictionCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override
    {
        return Kratos::make_intrusive<MPMGridNavierSlipFrictionCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    /// Create a new MPMGridNavierSlipFrictionCondition object.
    /**
      @param NewId Index of the new condition
      @param pGeom A pointer to the condition's geometry
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override
    {
        return Kratos::make_intrusive<MPMGridNavierSlipFrictionCondition>(NewId, pGeom, pProperties);
    }

    /**
     * Clones the selected element variables, creating a new one
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& rThisNodes
        ) const override
    {
        Condition::Pointer pNewCondition = Create(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

        pNewCondition->SetData(this->GetData());
        pNewCondition->SetFlags(this->GetFlags());

        return pNewCondition;
    }

    /**
     * @brief Provides the global indices for each one of this condition's local rows
     * This determines the elemental equation ID vector for all elemental DOFs
     * @param rResult Reference to the vector containing the global Id of each row
     * @param rCurrentProcessInfo Reference to the current ProcessInfo container
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * @brief Returns a list of the element's Dofs
     *
     * @param rConditionDofList Reference to the DOF pointers list
     * @param CurrentProcessInfo Reference to the current ProcessInfo container
     */
    void GetDofList(
        DofsVectorType& rConditionDofList,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * Getting method to obtain the time derivative of variable which defines the degrees of freedom
     * @param rValues The values of velocities
     * @param Step The step to be computed
     */
    void GetFirstDerivativesVector(
        Vector& values,
        int Step = 0
        ) const override;

    /**
     * this is called during the assembling process in order
     * to calculate the condition damping matrix
     * @param rDampingMatrix the condition damping matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * this is called during the assembling process in order
     * to calculate the condition left hand side matrix only
     * @param rLeftHandSideMatrix the condition left hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * this is called during the assembling process in order
     * to calculate the condition right hand side vector only
     * @param rRightHandSideVector the condition right hand side vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * this is called during the assembling process in order
     * to calculate all condition contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix the condition left hand side matrix
     * @param rRightHandSideVector the condition right hand side
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix: The LHS
     * @param rRightHandSideVector: The RHS
     * @param rCurrentProcessInfo: The current process info instance
     * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
     */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag
        );

    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the contitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * Called at the beginning of each solution step
     * @param rCurrentProcessInfo: the current process info instance
     */
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}

private:

    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
    }

    ///@}

}; // Class Condition

} // namespace Kratos.
