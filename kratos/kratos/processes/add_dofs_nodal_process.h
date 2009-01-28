/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:33 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_ADD_DOFS_NODAL_PROCESS_H_INCLUDED )
#define  KRATOS_ADD_DOFS_NODAL_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <sstream>


// External includes 


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "containers/pointer_vector_set.h"

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
  template<class TVariableType>
    class AddDofsNodalProcess : public Process
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of AddDofsNodalProcess
      KRATOS_CLASS_POINTER_DEFINITION(AddDofsNodalProcess);

      typedef Node<3> NodeType;

      typedef Dof<typename TVariableType::Type> DofType;

      typedef PointerVectorSet<NodeType, IndexedObject> NodesContainerType;

      typedef PointerVectorSet<DofType, IndexedObject> DofsContainerType;

      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Constructor.
      AddDofsNodalProcess(TVariableType const& NewVariable, NodesContainerType & rThisNodes, DofsContainerType& rThisDofs) 
	: mVariable(NewVariable), mNodes(rThisNodes), mDofs(rThisDofs){}

      /// Destructor.
      virtual ~AddDofsNodalProcess(){}
      

      ///@}
      ///@name Operators 
      ///@{

      
      ///@}
      ///@name Operations
      ///@{

      virtual void Execute()
	{
	  for(typename NodesContainerType::iterator i_node = mNodes.begin() ; i_node != mNodes.end() ; i_node++)
	    mDofs.push_back(i_node->pAddDof(mVariable));
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
	  buffer << "add "
		 << mVariable.Name()
		 << " dofs nodal process";
	  return buffer.str();
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
      TVariableType const& mVariable;

      NodesContainerType& mNodes;

      DofsContainerType& mDofs; 
        
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
      AddDofsNodalProcess& operator=(AddDofsNodalProcess const& rOther);

      /// Copy constructor.
      AddDofsNodalProcess(AddDofsNodalProcess const& rOther);

        
      ///@}    
        
    }; // Class AddDofsNodalProcess 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_ADD_DOFS_NODAL_PROCESS_H_INCLUDED  defined 


