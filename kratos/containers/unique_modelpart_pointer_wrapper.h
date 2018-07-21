//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Micheal Andre 
//

#if !defined(KRATOS_UNIQUE_MODELPART_POINTER_WRAPPER_H_INCLUDED )
#define  KRATOS_UNIQUE_MODELPART_POINTER_WRAPPER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "containers/model.h"
#include "includes/model_part.h"

namespace Kratos
{
  ///@addtogroup KratosCore
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
  
  /// This is a holder for a modelpart pointer which guarantees removal from Model
  /** This class holds internally a reference to a modelpart. On destruction of the
   * class the modelpart is removed from the OwnerModel
  */
  class UniqueModelPartPointerWrapper
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of UniqueModelPartPointerWrapper
      KRATOS_CLASS_POINTER_DEFINITION(UniqueModelPartPointerWrapper);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      UniqueModelPartPointerWrapper(Model& rModel,const std::string& Name)
      {
          if(rModel.HasModelPart(Name))
            rModel.DeleteModelPart(Name); //TODO: maybe we should throw an error here... to be decided
          mpModelPart = &rModel.CreateModelPart(Name);
      }

      UniqueModelPartPointerWrapper(Model& rModel,
            const std::string& Name,
            std::size_t buffer_size)
      {
          if(rModel.HasModelPart(Name))
            rModel.DeleteModelPart(Name); //TODO: maybe we should throw an error here... to be decided
          mpModelPart = &rModel.CreateModelPart(Name, buffer_size);
      }


      /// Destructor.
      virtual ~UniqueModelPartPointerWrapper()
      {
          auto& owner_model = mpModelPart->GetOwnerModel();
          const std::string name = mpModelPart->Name();
          owner_model.DeleteModelPart(name);
      }
      
      ModelPart& GetModelPart()
      {
          return *mpModelPart;
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
        buffer << "UniqueModelPartPointerWrapper containing model " << mpModelPart->Name() << std::endl;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "UniqueModelPartPointerWrapper";}

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
      ModelPart* mpModelPart;
        
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
      UniqueModelPartPointerWrapper& operator=(UniqueModelPartPointerWrapper const& rOther) = delete;

      /// Copy constructor.
      UniqueModelPartPointerWrapper(UniqueModelPartPointerWrapper const& rOther) = delete;

        
      ///@}    
        
    }; // Class UniqueModelPartPointerWrapper 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
//   inline std::istream& operator >> (std::istream& rIStream, 
// 				    UniqueModelPartPointerWrapper& rThis)
//                     {
//                         return rIStream;
//                     }

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const UniqueModelPartPointerWrapper& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_UNIQUE_MODELPART_POINTER_WRAPPER_H_INCLUDED  defined 


