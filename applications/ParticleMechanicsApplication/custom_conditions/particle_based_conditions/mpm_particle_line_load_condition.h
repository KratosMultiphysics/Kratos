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
#if !defined(KRATOS_MPM_PARTICLE_POINT_LINE_CONDITION_H_INCLUDED )
#define      KRATOS_MPM_PARTICLE_POINT_LINE_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_base_load_condition.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"

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

class MPMParticleLineLoadCondition
     : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MPMParticleLineLoadCondition
    KRATOS_CLASS_POINTER_DEFINITION( MPMParticleLineLoadCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MPMParticleLineLoadCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        );

    MPMParticleLineLoadCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        );

    /// Destructor.
    ~MPMParticleLineLoadCondition() override;

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

    ///@}
    ///@name Access
    ///@{

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
        std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    
    void SetValuesOnIntegrationPoints(
        const Variable<double>& rVariable,
        const std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void SetValuesOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
        const std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    
   


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

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{

    // A protected default constructor necessary for serialization
    MPMParticleLineLoadCondition() {};

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    array_1d<double, 3> m_line_load;
    array_1d<double, 3> m_delta_xg = {ZeroVector{3}};
    array_1d<double, 3> m_velocity_a;
    array_1d<double, 3> m_velocity_b;
    array_1d<double, 3> m_xg;
    array_1d<double, 3> m_acceleration_a;
    array_1d<double, 3> m_acceleration_b;
    array_1d<double, 3> m_normal;
    array_1d<double, 3> m_id_list;
    array_1d<double, 3> m_area_list;
    double m_area;
    array_1d<double, 3> m_point_load={ ZeroVector(3) };

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
        rSerializer.save("xg",m_xg);
        rSerializer.save("acceleration_a",m_acceleration_a);
        rSerializer.save("velocity_a",m_velocity_a);
        rSerializer.save("acceleration_b",m_acceleration_b);
        rSerializer.save("velocity_b",m_velocity_b);
        rSerializer.save("normal",m_normal);
        rSerializer.save("m_line_load",m_line_load);

    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
        rSerializer.load("xg",m_xg);
        rSerializer.load("acceleration_a",m_acceleration_a);
        rSerializer.load("acceleration_b",m_acceleration_b);
        rSerializer.load("velocity_a",m_velocity_a);
        rSerializer.load("velocity_b",m_velocity_b);
        rSerializer.load("normal",m_normal);
        rSerializer.load("m_line_load",m_line_load);
    }


}; // Class MPMParticleLineLoadCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}  // namespace Kratos.

#endif // KRATOS_MPM_PARTICLE_POINT_LOAD_CONDITION_H_INCLUDED  defined