//
//   Project Name:        Kratos       
//   Last modified by:    Miguel Masó Sotomayor
//   Date:                September 27th 2017
//   Revision:            1.0
//
//

#if !defined(KRATOS_WALL_CONDITION_H_INCLUDED)
#define  KRATOS_WALL_CONDITION_H_INCLUDED 

// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h" 
#include "includes/serializer.h"

namespace Kratos
{

  template< unsigned int TNumNodes >
  class WallCondition : public Condition
  {
  public:
     
    /// Counted pointer of WallCondition
    KRATOS_CLASS_POINTER_DEFINITION( WallCondition );

//----------------------------------------------------------------------

    /// Default constructor
    WallCondition():
        Condition()
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    WallCondition(IndexType NewId, GeometryType::Pointer pGeometry):
        Condition(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    WallCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
        Condition(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor
    WallCondition(WallCondition const& rOther):
        Condition(rOther)
    {
    }

    /// Destructor.
    virtual ~ WallCondition() {};

//----------------------------------------------------------------------

    /// Create a new WallCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
        return Condition::Pointer(new WallCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    void GetDofList(DofsVectorType& rConditionDofList,ProcessInfo& rCurrentProcessInfo);

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

//----------------------------------------------------------------------

  protected:



//----------------------------------------------------------------------

  private:

    friend class Serializer;

  }; // class WallCondition

}  // namespace Kratos.

#endif // KRATOS_WALL_CONDITION_H_INCLUDED  defined
