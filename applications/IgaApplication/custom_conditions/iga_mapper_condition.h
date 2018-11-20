//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    TT, PB
//

#if !defined(KRATOS_IGA_MAPPER_CONDITION_H_INCLUDED )
#define KRATOS_IGA_MAPPER_CONDITION_H_INCLUDED


// System includes

// External includes


// Project includes
#include "includes/condition.h"
#include "iga_application_variables.h"
#include "custom_utilities/nurbs_brep_modeler.h"


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

class IGAMapperCondition : public Condition {
public:


  ///@name Type Definitions
  ///@{

  typedef Condition BaseType;

  ///@}
  ///@name Pointer Definitions
  /// Pointer definition ofIGAMapperCondition
  KRATOS_CLASS_POINTER_DEFINITION(IGAMapperCondition);

  ///@}
  ///@name Life Cycle
  ///@{

  /**
   * Constructor.
   */
  IGAMapperCondition(IndexType NewId, NurbsBrepModeler* pModeler)
      : Condition(NewId), mpModeler(pModeler) {}

  /**
   * Copy Constructor
   */
 IGAMapperCondition(IGAMapperCondition const& rOther) = delete;

  /**
   * Destructor
   */
 virtual ~IGAMapperCondition() {};

  ///@}
  ///@name Operators
  ///@{

  /// Assignment operator.
  IGAMapperCondition & operator=(IGAMapperCondition const& rOther) = delete;

  ///@}
  ///@name Operations
  ///@{

  void Calculate(const Variable<Matrix >& rVariable,
                 Matrix& Output,
                 const ProcessInfo& rCurrentProcessInfo) override;

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
  virtual std::string Info() const;

  /// Print information about this object.
  virtual void PrintInfo(std::ostream& rOStream) const;

  /// Print object's data.
  virtual void PrintData(std::ostream& rOStream) const;

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

    NurbsBrepModeler* mpModeler;


  ///@}
  ///@name Private Operators
  ///@{

  ///@}
  ///@name Private Operations
  ///@{

  ///@}
  ///@name Serialization
  ///@{

  friend class Serializer;

  //virtual void save(Serializer& rSerializer) const;
  //virtual void load(Serializer& rSerializer);

  ///@}
  ///@name Private  Access
  ///@{

  ///@}
  ///@name Private Inquiry
  ///@{

  ///@}
  ///@name Un accessible methods
  ///@{

  ///@}

}; // ClassIGAMapperCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_IGA_MAPPER_CONDITION_H_INCLUDED  defined
