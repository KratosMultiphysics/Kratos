//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_EPETRA_DEFAULT_UTILITY_H_INCLUDED )
#define  KRATOS_EPETRA_DEFAULT_UTILITY_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "AztecOO.h"
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"


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
  class EpetraDefaultSetter
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of EpetraDefaultSetter
      KRATOS_CLASS_POINTER_DEFINITION(EpetraDefaultSetter);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      EpetraDefaultSetter(){}

      /// Destructor.
      virtual ~EpetraDefaultSetter(){}
      
      void SetDefaults(Teuchos::ParameterList& rlist, std::string settings_name)
      {
//	if(settings_name == std::string("SA") )
	  ML_Epetra::SetDefaults(settings_name.c_str(),rlist);
//	else
//	  std::cout << "WARNING: no defaults were set!!" << std::endl;
	
      }
      

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
	std::stringstream buffer;
        buffer << "EpetraDefaultSetter" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "EpetraDefaultSetter";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const {}
      
            
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
      EpetraDefaultSetter& operator=(EpetraDefaultSetter const& rOther){return *this;}

      /// Copy constructor.
      EpetraDefaultSetter(EpetraDefaultSetter const& rOther){}

        
      ///@}    
        
    }; // Class EpetraDefaultSetter 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    EpetraDefaultSetter& rThis){return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const EpetraDefaultSetter& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_EPETRA_DEFAULT_UTILITY_H_INCLUDED  defined 


