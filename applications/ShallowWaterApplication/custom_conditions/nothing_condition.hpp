//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: Kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#if !defined(KRATOS_NOTHING_CONDITION_H_INCLUDED)
#define  KRATOS_NOTHING_CONDITION_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

  template< unsigned int TNumNodes >
  class NothingCondition : public Condition
  {
  public:

    /// Counted pointer of NothingCondition
    KRATOS_CLASS_POINTER_DEFINITION( NothingCondition );

//----------------------------------------------------------------------

    /// Default constructor
    NothingCondition():
        Condition()
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    NothingCondition(IndexType NewId, GeometryType::Pointer pGeometry):
        Condition(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    NothingCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
        Condition(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor
    NothingCondition(NothingCondition const& rOther):
        Condition(rOther)
    {
    }

    /// Destructor.
    virtual ~ NothingCondition() {};

//----------------------------------------------------------------------

    /// Create a new NothingCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Condition::Pointer(new NothingCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& rConditionDofList,ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

//----------------------------------------------------------------------

  protected:



//----------------------------------------------------------------------

  private:

    friend class Serializer;

  }; // class NothingCondition

}  // namespace Kratos.

#endif // KRATOS_NOTHING_CONDITION_H_INCLUDED  defined
