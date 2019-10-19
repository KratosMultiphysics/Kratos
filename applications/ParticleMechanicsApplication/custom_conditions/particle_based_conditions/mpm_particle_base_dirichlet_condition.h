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

class MPMParticleBaseDirichletCondition
    : public MPMParticleBaseCondition
{

public:

    ///@name Type Definitions
    typedef std::size_t SizeType;
    ///@{

    // Counted pointer of MPMParticleBaseDirichletCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( MPMParticleBaseDirichletCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    MPMParticleBaseDirichletCondition()
    {};

    // Constructor using an array of nodes
    MPMParticleBaseDirichletCondition( IndexType NewId, GeometryType::Pointer pGeometry ):MPMParticleBaseCondition(NewId,pGeometry)
    {};

    // Constructor using an array of nodes with properties
    MPMParticleBaseDirichletCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ):MPMParticleBaseCondition(NewId,pGeometry,pProperties)
    {};

    // Destructor
    ~MPMParticleBaseDirichletCondition() override
    {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Called at the beginning of each solution step
     * @param rCurrentProcessInfo: the current process info instance
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Called at the end of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Calculate Shape Function Values in a given point
     */
    Vector& MPMShapeFunctionPointValues(Vector& rResult, const array_1d<double,3>& rPoint) override;

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMParticleBaseCondition );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMParticleBaseCondition );
    }

}; // class MPMParticleBaseDirichletCondition.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_MPM_PARTICLE_BASE_DIRICHLET_CONDITION_3D_H_INCLUDED  defined
