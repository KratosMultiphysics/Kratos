//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_SOFTENING_HARDENING_CRITERIA)
#define  KRATOS_SOFTENING_HARDENING_CRITERIA



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "structural_application.h"


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
  class SofteningHardeningCriteria
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of SofteningHardeningCriteria
      KRATOS_CLASS_POINTER_DEFINITION(SofteningHardeningCriteria);
        
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      SofteningHardeningCriteria(){}

      /// Destructor.
      virtual ~SofteningHardeningCriteria(){}
      
      
      virtual boost::shared_ptr<SofteningHardeningCriteria> Clone() const
      {
      boost::shared_ptr<SofteningHardeningCriteria> p_clone(new SofteningHardeningCriteria());
      return p_clone;
      }
     
      /*
      virtual double FunctionSofteningHardeningBehavior(const double& A, const double& r_o, const double& r) 
	{ 
	      return 0.00;
	}
      */	
    
      // Calcula la funcion de abalandamiento o endurecimiento junto con su derivada respecto kp_punto (0 <= kp_punto <= 1.00) 
      virtual void FunctionSofteningHardeningBehavior(const double& capap, const double& sigma, double& Result, double& der_Result) 
	{ 
	      return;
	}

      virtual void Linear_Strain_Softening(Vector& Principal_Stress)
	{
	  return;
	} 
 
      virtual void   InitializeMaterial(const Properties& props){mprops = &props;}
      virtual double Calculate(Vector& Imput_Parameters) {return 0;}
	
  
      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      
      
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
	return "The Base Class of Softening-Hardening Behavior";
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream)const
      {
	 rOStream << Info();
      }

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
	 return;
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
      const Properties *mprops;
      
        
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
      SofteningHardeningCriteria& operator=(SofteningHardeningCriteria const& rOther)
      {
          mprops= rOther.mprops;
	  return *this;
      }
      

      /// Copy constructor.
      SofteningHardeningCriteria(SofteningHardeningCriteria const& rOther):
          mprops(rOther.mprops)
      {
      }

        
      ///@}    
        
    }; // Class SofteningHardeningCriteria 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
   inline std::istream& operator >> (std::istream& rIStream, 
 				    SofteningHardeningCriteria& rThis);
 
   /// output stream function
   inline std::ostream& operator << (std::ostream& rOStream, 
 				    const SofteningHardeningCriteria& rThis)
     {
       rThis.PrintInfo(rOStream);
       rOStream << std::endl;
       rThis.PrintData(rOStream);
       return rOStream;     
     }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_SOFTENING_HARDENING_CRITERIA defined 