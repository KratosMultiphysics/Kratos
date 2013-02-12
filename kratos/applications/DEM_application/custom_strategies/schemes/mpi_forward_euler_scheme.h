//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelson Lafontaine $
//   Date:                $Date: 2012-04-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//



#if !defined(KRATOS_MPI_FORWARD_EULER_SCHEME_H_INCLUDED )
#define  KRATOS_MPI_FORWARD_EULER_SCHEME_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "integration_scheme.h"
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"

#include "DEM_application.h"

#include "custom_strategies/schemes/forward_euler_scheme.h"

namespace Kratos
{
  
  
  class MpiForwardEulerScheme :  public ForwardEulerScheme
    {
    public:
      ///@name Type Definitions
      ///@{
      
      typedef ModelPart::NodesContainerType NodesArrayType;
    
    
      /// Pointer definition of ForwardEulerScheme
      KRATOS_CLASS_POINTER_DEFINITION(MpiForwardEulerScheme);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      MpiForwardEulerScheme(){}

      /// Destructor.
      virtual ~MpiForwardEulerScheme(){}
      
      virtual NodesArrayType& GetNodes(ModelPart& model_part)
      {
          return model_part.GetCommunicator().LocalMesh().Nodes(); 
      }
 
      /// Turn back information as a string.
      virtual std::string Info() const
      {
          std::stringstream buffer;
          buffer << "ForwardEulerScheme" ;
          return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "ForwardEulerScheme";}

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
     MpiForwardEulerScheme& operator=(MpiForwardEulerScheme const& rOther)
     {
       return *this;
     }

      /// Copy constructor.
     MpiForwardEulerScheme(MpiForwardEulerScheme const& rOther)
     {
       *this = rOther;
     }

        
      ///@}    
        
    }; // Class ForwardEulerScheme 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
                    MpiForwardEulerScheme& rThis){return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
                    const MpiForwardEulerScheme& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_MPI_FORWARD_EULER_SCHEME_H_INCLUDED  defined 

