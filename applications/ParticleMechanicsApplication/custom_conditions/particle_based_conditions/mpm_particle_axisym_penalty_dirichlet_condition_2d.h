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


#if !defined(KRATOS_MPM_PARTICLE_AXISYM_PENALTY_DIRICHLET_CONDITION_2D_H_INCLUDED )
#define      KRATOS_MPM_PARTICLE_AXISYM_PENALTY_DIRICHLET_CONDITION_2D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/particle_based_conditions/mpm_particle_penalty_dirichlet_condition.h"

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

/// Axisymmetric penalty dirichlet condition

class MPMParticleAxisymPenaltyDirichletCondition2D
    : public MPMParticlePenaltyDirichletCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MPMParticleAxisymPenaltyDirichletCondition2D
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(MPMParticleAxisymPenaltyDirichletCondition2D);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MPMParticleAxisymPenaltyDirichletCondition2D(IndexType NewId, GeometryType::Pointer pGeometry);
    MPMParticleAxisymPenaltyDirichletCondition2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    ~MPMParticleAxisymPenaltyDirichletCondition2D() override;

    ///@}
    ///@name Operators
    ///@{
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

    //std::string Info() const;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    //      virtual String Info() const;

    /// Print information about this object.
    //      virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    //      virtual void PrintData(std::ostream& rOStream) const;
    ///@}
    ///@name Friends
    ///@{
    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    MPMParticleAxisymPenaltyDirichletCondition2D() : MPMParticlePenaltyDirichletCondition()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{
    ///@}
    ///@name Protected  Access
    ///@{
    ///@}
    ///@name Protected Inquiry
    ///@{
    ///@}
    ///@name Protected LifeCycle
    ///@{
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

    /**
     * This functions computes the integration weight to consider
     * @param IntegrationPoints: The array containing the integration points
     * @param PointNumber: The id of the integration point considered
     * @param detJ: The determinant of the jacobian of the element
     */
    double GetIntegrationWeight() override;

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    /// Assignment operator.
    //MPMParticleAxisymPenaltyDirichletCondition2D& operator=(const MPMParticleAxisymPenaltyDirichletCondition2D& rOther);
    /// Copy constructor.
    //MPMParticleAxisymPenaltyDirichletCondition2D(const MPMParticleAxisymPenaltyDirichletCondition2D& rOther);
    ///@}

}; // Class MPMParticleAxisymPenaltyDirichletCondition2D

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_MPM_PARTICLE_AXISYM_PENALTY_DIRICHLET_CONDITION_2D_H_INCLUDED  defined
