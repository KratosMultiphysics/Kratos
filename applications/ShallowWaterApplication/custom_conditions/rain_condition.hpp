//   
//   Project Name:        Kratos       
//   Last modified by:    Miguel Masó Sotomayor
//   Date:                July 31st 2017
//   Revision:            1.3
//
//

#if !defined(KRATOS_RAIN_CONDITION_H_INCLUDED)
#define  KRATOS_RAIN_CONDITION_H_INCLUDED 

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

  template< unsigned int TNumNodes >
  class RainCondition : public Condition
  {
  public:
     
    /// Counted pointer of RainCondition
    KRATOS_CLASS_POINTER_DEFINITION( RainCondition );


    /// Default constructor

    /// Constructor using Geometry
    RainCondition(IndexType NewId, GeometryType::Pointer pGeometry):
        Condition(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    RainCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
        Condition(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor
    RainCondition(RainCondition const& rOther):
        Condition(rOther)
    {
    }

    /// Destructor.
    virtual ~ RainCondition() {};

//----------------------------------------------------------------------

    /// Create a new RainCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
        return Condition::Pointer(new RainCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    void GetDofList(DofsVectorType& rConditionDofList,ProcessInfo& rCurrentProcessInfo);

    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);

    double mHeightUnitConvert;

//----------------------------------------------------------------------

  protected:

    void CalculateConsistentMassMatrix(boost::numeric::ublas::bounded_matrix<double,TNumNodes,TNumNodes>& rMassMatrix);

    void CalculateLumpedMassMatrix(boost::numeric::ublas::bounded_matrix<double,TNumNodes,TNumNodes>& rMassMatrix);

//----------------------------------------------------------------------

  private:

    friend class Serializer;

    // A private default constructor necessary for serialization
    RainCondition() : Condition()
    {
    }

       
  }; // Class RainCondition
}  // namespace Kratos.

#endif // KRATOS_RAIN_CONDITION_H_INCLUDED  defined
