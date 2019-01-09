//
//   Project Name:        KratosSolversApplication    $
//   Created by:          $Author:        JMCarbonell $
//   Last modified by:    $Co-Author:                 $
//   Date:                $Date:         January 2019 $
//   Revision:            $Revision:              0.0 $
//
//

#if !defined(KRATOS_SOLVERS_APPLICATION_H_INCLUDED )
#define  KRATOS_SOLVERS_APPLICATION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/kratos_application.h"

namespace Kratos {

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
class KratosSolversApplication : public KratosApplication {
 public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of KratosSolversApplication
  KRATOS_CLASS_POINTER_DEFINITION(KratosSolversApplication);

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  KratosSolversApplication();

  /// Destructor.
  ~KratosSolversApplication() override {}

  ///@}
  ///@name Operators
  ///@{
  ///@}
  ///@name Operations
  ///@{

  void Register() override;

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
  std::string Info() const override
  {
    return "KratosSolversApplication";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << Info();
    PrintData(rOStream);
  }

  ///// Print object's data.
  void PrintData(std::ostream& rOStream) const override
  {
    KRATOS_WATCH("in Solvers application");
    KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
    rOStream << "Variables:" << std::endl;
    KratosComponents<VariableData>().PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "Elements:" << std::endl;
    KratosComponents<Element>().PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "Conditions:" << std::endl;
    KratosComponents<Condition>().PrintData(rOStream);
  }

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
  KratosSolversApplication& operator=(KratosSolversApplication const& rOther);

  /// Copy constructor.
  KratosSolversApplication(KratosSolversApplication const& rOther);


  ///@}

}; // Class KratosSolversApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_SOLVERS_APPLICATION_H_INCLUDED  defined
