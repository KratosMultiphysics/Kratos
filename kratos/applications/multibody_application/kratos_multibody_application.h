//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-12-10 18:45:27 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_KRATOS_APPLICATIONNAME_APPLICATION_H_INCLUDED )
#define  KRATOS_KRATOS_APPLICATIONNAME_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"


namespace Kratos
{

  ///@name Kratos Globals
  ///@{ 
  
  // Variables definition 
  KRATOS_DEFINE_VARIABLE(double, IS_BODY );
  

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
  class KratosMultibodyApplication : public KratosApplication
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of KratosMultibodyApplication
      KRATOS_CLASS_POINTER_DEFINITION(KratosMultibodyApplication);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      KratosMultibodyApplication(){}

      /// Destructor.
      virtual ~KratosMultibodyApplication(){}
      

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
	return "KratosMultibodyApplication";
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
      {
	rOStream << Info();
      }

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
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
      

      // Static references for elements and conditions
//       static const ApplicationElement  msApplicationElement; 
//       static const ApplicationCondition  msApplicationCondition; 
        
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
      KratosMultibodyApplication& operator=(KratosMultibodyApplication const& rOther);

      /// Copy constructor.
      KratosMultibodyApplication(KratosMultibodyApplication const& rOther);

        
      ///@}    
        
    }; // Class KratosMultibodyApplication 

  ///@} 
  

  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_KRATOS_APPLICATIONNAME_APPLICATION_H_INCLUDED  defined 


