//--------------------------------------------------------
//          ___  __                                      .
//  KRATOS | _ \/ _|___ _ __                             .
//         |  _/  _/ -_) '  \                            .
//         |_| |_| \___|_|_|_| APPLICATION               .
//                                                       .
//  License:(BSD)         PfemApplication/license.txt    .
//  Main authors:         Josep Maria Carbonell          .
//                        ..                             .
//--------------------------------------------------------
//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_PFEM_APPLICATION_H_INCLUDED )
#define  KRATOS_PFEM_APPLICATION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_application.h"

//elements
#include "custom_elements/fluid_elements/updated_lagrangian_segregated_fluid_element.hpp"

#include "pfem_application_variables.h"

namespace Kratos
{
  ///@name Type	Definitions
  ///@{

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
  class KRATOS_API(PFEM_APPLICATION) KratosPfemApplication : public KratosApplication
  {
  public:


    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosPfemApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosPfemApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosPfemApplication    ();

    /// Destructor.
    virtual ~KratosPfemApplication    (){}

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
      return "KratosPfemApplication";
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
      KRATOS_WATCH( "in KratosPfemApplication" )
      KRATOS_WATCH( KratosComponents<VariableData>::GetComponents().size() )
      rOStream << "Variables:" << std::endl;
      KratosComponents<VariableData>().PrintData(rOStream);
      rOStream << std::endl;
      rOStream << "Elements:" << std::endl;
      KratosComponents<Element>().PrintData(rOStream);
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

    const UpdatedLagrangianSegregatedFluidElement  mUpdatedLagrangianSegregatedFluidElement2D3N;
    const UpdatedLagrangianSegregatedFluidElement  mUpdatedLagrangianSegregatedFluidElement3D4N;

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
    KratosPfemApplication& operator=(KratosPfemApplication const& rOther);

    /// Copy constructor.
    KratosPfemApplication(KratosPfemApplication const& rOther);

    ///@}

  }; // Class KratosPfemApplication

  ///@}


  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}


}  // namespace Kratos.

#endif // KRATOS_PFEM_APPLICATION_H_INCLUDED  defined


