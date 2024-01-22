//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


// System includes
#if !defined(KRATOS_MPM_PARTICLE_PENALTY_DIRICHLET_CONDITION_H_INCLUDED )
#define      KRATOS_MPM_PARTICLE_PENALTY_DIRICHLET_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_base_dirichlet_condition.h"
#include "mpm_application_variables.h"

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

/// Short class definition.
/** Detail class definition.
*/

class MPMParticlePenaltyDirichletCondition
    : public MPMParticleBaseDirichletCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MPMParticlePenaltyDirichletCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( MPMParticlePenaltyDirichletCondition );

    using MPMParticleBaseDirichletCondition::CalculateOnIntegrationPoints;
    using MPMParticleBaseDirichletCondition::SetValuesOnIntegrationPoints;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MPMParticlePenaltyDirichletCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        );

    MPMParticlePenaltyDirichletCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        );

    /// Destructor.
    ~MPMParticlePenaltyDirichletCondition() override;

    ///@}
    ///@name Operators
    ///@{

    /**
     * Called at the beginning of each solution step
     * @param rCurrentProcessInfo: the current process info instance
     */
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Called at the begining at each nonlinear iteration
     * @param rCurrentProcessInfo the current process info instance
     */
    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called for non-linear analysis at the end of the iteration process
     */
    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Called at the end of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        ) const override;

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;


    ///@}
    ///@name Access Get Values
    ///@{

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access Set Values
    ///@{

    void SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
        const std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    double m_penalty = 0.0;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

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
        ) override;

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{
    virtual void CalculateInterfaceContactForce(const ProcessInfo& rCurrentProcessInfo );


    ///@}
    ///@name Protected LifeCycle
    ///@{

    // A protected default constructor necessary for serialization
    MPMParticlePenaltyDirichletCondition() {};

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMParticleBaseDirichletCondition );
        rSerializer.save("penalty", m_penalty);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMParticleBaseDirichletCondition );
        rSerializer.load("penalty", m_penalty);
    }


}; // Class MPMParticlePenaltyDirichletCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}  // namespace Kratos.

#endif // KRATOS_MPM_PARTICLE_PENALTY_DIRICHLET_CONDITION_H_INCLUDED  defined


