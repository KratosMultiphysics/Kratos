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
#if !defined(KRATOS_MPM_LINE_LOAD_CONDITION_2D_H_INCLUDED )
#define  KRATOS_MPM_LINE_LOAD_CONDITION_2D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/mpm_base_load_condition.h"
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

class MPMLineLoadCondition2D
    : public MPMBaseLoadCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MPMLineLoadCondition2D
    KRATOS_CLASS_POINTER_DEFINITION( MPMLineLoadCondition2D );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MPMLineLoadCondition2D( IndexType NewId, GeometryType::Pointer pGeometry );
    MPMLineLoadCondition2D( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

    /// Destructor.
    ~MPMLineLoadCondition2D() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
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
        ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag 
        ) override;

    void CalculateAndSubKp(
        Matrix& K,
        const Matrix& DN_De,
        const Vector& N,
        const double Pressure,
        const double IntegrationWeight
        );

    void CalculateAndAddPressureForce(
        VectorType& rRightHandSideVector,
        const Vector& N,
        const array_1d<double, 3>& Normal,
        const double Pressure,
        const double IntegrationWeight 
        );

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
    MPMLineLoadCondition2D() {};

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMBaseLoadCondition );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMBaseLoadCondition );
    }


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //MPMLineLoadCondition2D& operator=(const MPMLineLoadCondition2D& rOther);

    /// Copy constructor.
    //MPMLineLoadCondition2D(const MPMLineLoadCondition2D& rOther);


    ///@}

}; // Class MPMLineLoadCondition2D

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
        MPMLineLoadCondition2D& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
        const MPMLineLoadCondition2D& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_MPM_LINE_LOAD_CONDITION_2D_H_INCLUDED  defined 


