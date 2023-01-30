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
#if !defined(KRATOS_MPM_PARTICLE_POINT_LOAD_CONDITION_H_INCLUDED )
#define      KRATOS_MPM_PARTICLE_POINT_LOAD_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_base_load_condition.h"
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

class MPMParticlePointLoadCondition
    : public MPMParticleBaseLoadCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MPMParticlePointLoadCondition
    KRATOS_CLASS_POINTER_DEFINITION( MPMParticlePointLoadCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MPMParticlePointLoadCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        );

    MPMParticlePointLoadCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        );

    /// Destructor.
    ~MPMParticlePointLoadCondition() override;

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

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
        std::vector<array_1d<double, 3 > >& rValues,
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

    /**
     * It calcules the integration load for the point load
     */
    double GetPointLoadIntegrationWeight() override;

    virtual void MPMShapeFunctionPointValuesKinematic(Vector& rResult) const;

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
     * Calculation of the Nodal Force
     */
    Matrix& CalculateNodalForce(Matrix & rNodalForce, const ProcessInfo& rCurrentProcessInfo);

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
    MPMParticlePointLoadCondition() {};

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    array_1d<double, 3> m_point_load;
    array_1d<double, 3> m_delta_xg;

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMParticleBaseLoadCondition );
        rSerializer.save("point_load", m_point_load);
        rSerializer.save("delta_xg", m_delta_xg);

    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMParticleBaseLoadCondition );
        rSerializer.load("point_load", m_point_load);
        rSerializer.load("delta_xg", m_delta_xg);
    }


}; // Class MPMParticlePointLoadCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}  // namespace Kratos.

#endif // KRATOS_MPM_PARTICLE_POINT_LOAD_CONDITION_H_INCLUDED  defined