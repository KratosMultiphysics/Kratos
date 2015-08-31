//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Last modified by:    $Author:               JMCarbonell $
//   Date:                $Date:              September 2015 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_PFEM_FLUID_DYNAMICS_APPLICATION_H_INCLUDED )
#define  KRATOS_PFEM_FLUID_DYNAMICS_APPLICATION_H_INCLUDED


// System includes

// External includes 

// Project includes

// Core applications
#include "solid_mechanics_application.h"

//conditions

//elements

//constitutive laws
#include "containers/flags.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"

// yield Criteria

//flow rule

//hardening laws

//constitutive laws

namespace Kratos
{
  ///@name Type	Definitions
  ///@{

  ///@name Kratos Globals
  ///@{ 

  //Define Variables

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
  //class KratosPfemFluidDynamicsApplication : public KratosPfemSolidMechanicsApplication
  class KratosPfemFluidDynamicsApplication : public KratosApplication
  {
  public:


    ///@name Type Definitions
    ///@{
		

    /// Pointer definition of KratosPfemFluidDynamicsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosPfemFluidDynamicsApplication);


    ///@}
    ///@name Life Cycle 
    ///@{ 

    /// Default constructor.
    KratosPfemFluidDynamicsApplication();

    /// Destructor.
    virtual ~KratosPfemFluidDynamicsApplication(){}


    ///@}
    ///@name Operators 
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual void Register();



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
    virtual std::string Info() const
      {
	return "KratosPfemFluidDynamicsApplication";
      }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
      rOStream << Info();
      PrintData(rOStream);
    }

    ///// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
      KRATOS_WATCH( "in KratosPfemFluidDynamicsApplication" ) 
      KRATOS_WATCH( KratosComponents<VariableData>::GetComponents().size() )
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
    KratosPfemFluidDynamicsApplication& operator=(KratosPfemFluidDynamicsApplication const& rOther);

    /// Copy constructor.
    KratosPfemFluidDynamicsApplication(KratosPfemFluidDynamicsApplication const& rOther);


    ///@}    

  }; // Class KratosPfemFluidDynamicsApplication 

  ///@} 


  ///@name Type Definitions       
  ///@{ 


  ///@} 
  ///@name Input and output 
  ///@{ 

  ///@} 


}  // namespace Kratos.

#endif // KRATOS_PFEM_FLUID_DYNAMICS_APPLICATION_H_INCLUDED  defined 


