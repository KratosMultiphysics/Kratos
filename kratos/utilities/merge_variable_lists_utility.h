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
//
	           
#if !defined(KRATOS_MERGE_VARIABLE_LISTS_H_INCLUDED )
#define  KRATOS_MERGE_VARIABLE_LISTS_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
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
  
  /// Merges the variable lists of the input modelparts.
  /** This function ensures that all the variables added to the variable list of modelpart1 
   * are also added to the variable list of modelpart2. Symmetrically every variable in modelpart2 
   * will also be added to modelpart1
  */
  class MergeVariableListsUtility
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of MergeVariableListsUtility
      KRATOS_CLASS_POINTER_DEFINITION(MergeVariableListsUtility);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      MergeVariableListsUtility(){}

      /// Destructor.
      virtual ~MergeVariableListsUtility(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      void Merge(ModelPart& model_part1, ModelPart& model_part2)
    {
        KRATOS_TRY
        
        auto& list1 = model_part1.GetNodalSolutionStepVariablesList();
        auto& list2 = model_part2.GetNodalSolutionStepVariablesList();

        for(const auto& var : list1)
        {
            model_part2.GetNodalSolutionStepVariablesList().Add(var);
        }

        for(const auto& var : list2)
        { 
            model_part1.GetNodalSolutionStepVariablesList().Add(var);
        }

        KRATOS_CATCH("")
    }
      
      
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
        buffer << "MergeVariableListsUtility" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "MergeVariableListsUtility";}

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
      MergeVariableListsUtility& operator=(MergeVariableListsUtility const& rOther) = delete;

      /// Copy constructor.
      MergeVariableListsUtility(MergeVariableListsUtility const& rOther) = delete;

        
      ///@}    
        
    }; // Class MergeVariableListsUtility 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    MergeVariableListsUtility& rThis)
    {
      return rIStream;
    }

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const MergeVariableListsUtility& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_MERGE_VARIABLE_LISTS_H_INCLUDED  defined 


