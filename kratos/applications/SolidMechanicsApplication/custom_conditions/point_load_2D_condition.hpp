//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_POINT_LOAD_2D_CONDITION_H_INCLUDED )
#define  KRATOS_POINT_LOAD_2D_CONDITION_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
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
class PointLoad2DCondition
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of PointLoad2DCondition
    KRATOS_CLASS_POINTER_DEFINITION(PointLoad2DCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PointLoad2DCondition(IndexType NewId, GeometryType::Pointer pGeometry);
    PointLoad2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Copy constructor
    PointLoad2DCondition( PointLoad2DCondition const& rOther);


    /// Destructor.
    virtual ~PointLoad2DCondition();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(IndexType NewId, NodesArrayType const&
                              ThisNodes,  PropertiesType::Pointer pProperties) const;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType&
                              rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo&
                                rCurrentProcessInfo);
    //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo&
                          rCurrentProcessInfo);

    void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo&
                    CurrentProcessInfo);

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

    // A protected default constructor necessary for serialization
    PointLoad2DCondition() {};

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


    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
    }

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
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //PointLoad2DCondition& operator=(const PointLoad2DCondition& rOther);

    /// Copy constructor.
    //PointLoad2DCondition(const PointLoad2DCondition& rOther);


    ///@}

}; // Class PointLoad2DCondition

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    PointLoad2DCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const PointLoad2DCondition& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_POINT_LOAD_2D_CONDITION_H_INCLUDED  defined 



