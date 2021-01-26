//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Veronika Singer
//


// System includes
#if !defined(KRATOS_MPM_PARTICLE_BASE_LOAD_CONDITION_H_INCLUDED )
#define      KRATOS_MPM_PARTICLE_BASE_LOAD_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_base_condition.h"
#include "includes/variables.h"

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

class MPMParticleBaseLoadCondition
    : public MPMParticleBaseCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MPMParticleBaseLoadCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( MPMParticleBaseLoadCondition );

    // Constructor void
    MPMParticleBaseLoadCondition()
    {};

    // Constructor using an array of nodes
    MPMParticleBaseLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry ):MPMParticleBaseCondition(NewId,pGeometry)
    {};

    // Constructor using an array of nodes with properties
    MPMParticleBaseLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ):MPMParticleBaseCondition(NewId,pGeometry,pProperties)
    {};

    // Destructor
    ~MPMParticleBaseLoadCondition() override
    {};

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


    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * It calcules the integration load for the point load
     */
    virtual double GetPointLoadIntegrationWeight();

    /**
     * Calculate Shape Function Values as a vector
     */

    virtual void MPMShapeFunctionPointValues(Vector& rResult) const override;

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


}; // Class MPMParticleBaseLoadCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}  // namespace Kratos.

#endif // KRATOS_MPM_PARTICLE_BASE_LOAD_CONDITION_H_INCLUDED  defined