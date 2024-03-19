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


#if !defined(KRATOS_MPM_PARTICLE_BASE_DIRICHLET_CONDITION_3D_H_INCLUDED )
#define      KRATOS_MPM_PARTICLE_BASE_DIRICHLET_CONDITION_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_base_condition.h"
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

class MPMParticleBaseDirichletCondition
    : public MPMParticleBaseCondition
{

public:

    ///@name Type Definitions
    typedef std::size_t SizeType;
    ///@{

    // Counted pointer of MPMParticleBaseDirichletCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( MPMParticleBaseDirichletCondition );

    using MPMParticleBaseCondition::CalculateOnIntegrationPoints;
    using MPMParticleBaseCondition::SetValuesOnIntegrationPoints;

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    MPMParticleBaseDirichletCondition()
    {};

    // Constructor using an array of nodes
    MPMParticleBaseDirichletCondition( IndexType NewId, GeometryType::Pointer pGeometry )
        : MPMParticleBaseCondition(NewId,pGeometry)
    {};

    // Constructor using an array of nodes with properties
    MPMParticleBaseDirichletCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : MPMParticleBaseCondition(NewId,pGeometry,pProperties)
    {};

    // Destructor
    ~MPMParticleBaseDirichletCondition() override
    {};

    ///@}
    ///@name Operators
    ///@{


    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check( const ProcessInfo& rCurrentProcessInfo ) const override;

    /**
     * Called at the beginning of each solution step
     * @param rCurrentProcessInfo: the current process info instance
     */
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Called at the end of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access Get Values
    ///@{

    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3 > >& rVariable,
        std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access Set Values
    ///@{

    void SetValuesOnIntegrationPoints(
        const Variable<array_1d<double, 3 > >& rVariable,
        const std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}

protected:
    ///@name Protected member Variables
    ///@{


    array_1d<double, 3> m_imposed_displacement;
    array_1d<double, 3> m_imposed_velocity;
    array_1d<double, 3> m_imposed_acceleration;
    array_1d<double, 3> m_contact_force;

    ///@}
    ///@name Protected Operations
    ///@{

    /// Calculate Shape Function Values as a vector
    virtual void MPMShapeFunctionPointValues(Vector& rResult) const override;

    ///@}

private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMParticleBaseCondition );
        rSerializer.save("imposed_displacement",m_imposed_displacement);
        rSerializer.save("imposed_velocity",m_imposed_velocity);
        rSerializer.save("imposed_acceleration",m_imposed_acceleration);
        rSerializer.save("contact_force",m_contact_force);

    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMParticleBaseCondition );
        rSerializer.load("imposed_displacement",m_imposed_displacement);
        rSerializer.load("imposed_velocity",m_imposed_velocity);
        rSerializer.load("imposed_acceleration",m_imposed_acceleration);
        rSerializer.load("contact_force",m_contact_force);
    }

    ///@}

}; // class MPMParticleBaseDirichletCondition.

///@}

} // namespace Kratos.

#endif // KRATOS_MPM_PARTICLE_BASE_DIRICHLET_CONDITION_3D_H_INCLUDED  defined
