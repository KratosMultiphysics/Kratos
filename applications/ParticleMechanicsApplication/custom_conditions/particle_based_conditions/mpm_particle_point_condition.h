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
#if !defined(KRATOS_MPM_PARTICLE_POINT_CONDITION_H_INCLUDED )
#define      KRATOS_MPM_PARTICLE_POINT_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"
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

/// Short class definition.
/** Detail class definition.
*/

class MPMParticlePointCondition
     : public Condition
{
    struct GeneralVariables
    {
    public:
        Vector  N;
        // Variables including all integration points
        Matrix  CurrentDisp;
    };
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MPMParticlePointCondition
    KRATOS_CLASS_POINTER_DEFINITION( MPMParticlePointCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MPMParticlePointCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        );

    MPMParticlePointCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        );

    /// Destructor.
    ~MPMParticlePointCondition() override;

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

    void CalculateOnIntegrationPoints(const Variable<int>& rVariable,
        std::vector<int>& rValues,
        const ProcessInfo& rCurrentProcessInfo)override;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
        std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void SetValuesOnIntegrationPoints(
        const Variable<int>& rVariable,
        const std::vector<int>& rValues,
        const ProcessInfo& rCurrentProcessInfo)override;

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
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called for non-linear analysis at the end of the iteration process
     */
    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Called at the end of eahc solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void FinalizeSolutionStep(const ProcessInfo& CurrentProcessInfo) override;

    /**
     * Calculation of the Current Displacement
     */
    Matrix& CalculateCurrentDisp(Matrix & rCurrentDisp, const ProcessInfo& rCurrentProcessInfo);


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
    MPMParticlePointCondition() {};

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    array_1d<double, 3> m_displacement;
    array_1d<double, 3> m_velocity;
    array_1d<double, 3> m_xg;
    array_1d<double, 3> m_acceleration;
    array_1d<double, 3> m_normal;
    array_1d<double, 3> m_contact_force;
    array_1d<double, 3> m_imposed_velocity;
    array_1d<double, 3> m_imposed_acceleration;
    array_1d<double, 3> m_imposed_displacement;
    double m_area;

    int m_corresponding_condition_id = 0;
    int m_corresponding_node_id = 0;

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
        rSerializer.save("acceleration",m_acceleration);
        rSerializer.save("velocity",m_velocity);
        rSerializer.save("normal",m_normal);
        rSerializer.save("area",m_area);
        rSerializer.save("displacement",m_displacement);
        rSerializer.save("imposed_displacement",m_imposed_displacement);
        rSerializer.save("imposed_velocity",m_imposed_velocity);
        rSerializer.save("imposed_acceleration",m_imposed_acceleration);
    

    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
        rSerializer.load("xg",m_xg);
        rSerializer.load("acceleration",m_acceleration);
        rSerializer.load("velocity",m_velocity);
        rSerializer.load("normal",m_normal);
        rSerializer.load("area",m_area);
        rSerializer.load("displacement",m_displacement);
        rSerializer.load("imposed_displacement",m_imposed_displacement);
        rSerializer.load("imposed_velocity",m_imposed_velocity);
        rSerializer.load("imposed_acceleration",m_imposed_acceleration);
        
    }


}; // Class MPMParticlePointCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}  // namespace Kratos.

#endif // KRATOS_MPM_PARTICLE_POINT_CONDITION_H_INCLUDED  defined