// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes
#if !defined(KRATOS_LINE_LOAD_CONDITION_2D_H_INCLUDED )
#define  KRATOS_LINE_LOAD_CONDITION_2D_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "custom_conditions/base_load_condition.h"
#include "includes/ublas_interface.h"
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

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION)  LineLoadCondition2D
    : public BaseLoadCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of LineLoadCondition2D
    KRATOS_CLASS_POINTER_DEFINITION( LineLoadCondition2D );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LineLoadCondition2D( IndexType NewId, GeometryType::Pointer pGeometry );
    LineLoadCondition2D( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

    /// Destructor.
    virtual ~LineLoadCondition2D();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    Condition::Pointer Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const;
    Condition::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const;

    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo );

    void CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo );

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
    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                       ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag,
                       bool CalculateResidualVectorFlag );

    void CalculateAndSubKp(
        Matrix& K,
        const Matrix& DN_De,
        const Vector& N,
        double pressure,
        double weight
    );


    void CalculateAndAdd_PressureForce(
        Vector& residualvector,
        const Vector& N,
        Vector& v3,
        double pressure,
        double weight );
    
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

    // A private default constructor necessary for serialization
    LineLoadCondition2D() {};

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseLoadCondition );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseLoadCondition );
    }


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //LineLoadCondition2D& operator=(const LineLoadCondition2D& rOther);

    /// Copy constructor.
    //LineLoadCondition2D(const LineLoadCondition2D& rOther);


    ///@}

}; // Class LineLoadCondition2D

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
        LineLoadCondition2D& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
        const LineLoadCondition2D& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_LINE_LOAD_CONDITION_2D_H_INCLUDED  defined 


