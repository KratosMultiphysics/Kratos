//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_BODY_RECOGNITION_PROCESS_H_INCLUDED )
#define  KRATOS_BODY_RECOGNITION_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "kratos_multibody_application.h"
#include "processes/process.h"
#include "processes/find_nodal_neighbours_process.h"
#include "includes/model_part.h"


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
  
  /// This Process takes a mesh and extract out the individual bodies by analysing the connectivity.
  /** .
  */
  class BodyRecognitionProcess : public Process
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of BodyRecognitionProcess
      KRATOS_CLASS_POINTER_DEFINITION(BodyRecognitionProcess);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      BodyRecognitionProcess(ModelPart& rModelPart, int Dimension) : mrModelPart(rModelPart), mDimension(Dimension)
      {
	std::cout << "Is there any body out there!?!" << std::endl;
      }

      /// Copy constructor.
      BodyRecognitionProcess(BodyRecognitionProcess const& rOther) : mrModelPart(rOther.mrModelPart), mDimension(rOther.mDimension)
      {
      }

      /// Destructor.
      virtual ~BodyRecognitionProcess(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      void operator()()
	{
	  Execute();
	}
      
      
      ///@}
      ///@name Operations
      ///@{
      
      virtual void Execute()
      {
	FindNodalNeighboursProcess find_nodal_neighbours_process(mrModelPart);
	find_nodal_neighbours_process.Execute();
	
	ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
	
	int current_body_id = 1;
	for(ModelPart::NodesContainerType::iterator i_node = r_nodes.begin(); i_node!=r_nodes.end(); i_node++)
	{
	WeakPointerVector<Element >& r_neighbour_elements = i_node->GetValue(NEIGHBOUR_ELEMENTS);
	  if(i_node->GetSolutionStepValue(IS_BODY) == 0)
	    AssignBody(i_node, current_body_id++);
	}
      }
      
      void AssignBody(ModelPart::NodesContainerType::iterator iNode, int BodyId)
      {
	ModelPart::NodesContainerType front_nodes;
	WeakPointerVector<Element >& r_neighbour_elements = iNode->GetValue(NEIGHBOUR_ELEMENTS);
	for(WeakPointerVector<Element >::iterator i_neighbour_element = r_neighbour_elements.begin() ; i_neighbour_element != r_neighbour_elements.end() ; i_neighbour_element++)
	{
	  if(i_neighbour_element->GetValue(IS_BODY) == 0)
	  {
	    i_neighbour_element->SetValue(IS_BODY, BodyId);
		
	    Element::GeometryType& p_geometry = i_neighbour_element->GetGeometry();
				
	    for(unsigned int i = 0; i < p_geometry.size(); i++)
	      if(p_geometry[i].GetSolutionStepValue(IS_BODY) == 0)
	      {
		p_geometry[i].GetSolutionStepValue(IS_BODY) = BodyId;
		front_nodes.push_back(p_geometry(i));
	      }
	  }
	}
	while(!front_nodes.empty())
	{
	  ModelPart::NodesContainerType new_front_nodes;
	  for(ModelPart::NodesContainerType::iterator i_node = front_nodes.begin() ; i_node != front_nodes.end() ; i_node++)
	  {
	    WeakPointerVector<Element >& r_neighbour_elements = i_node->GetValue(NEIGHBOUR_ELEMENTS);
	    for(WeakPointerVector<Element >::iterator i_neighbour_element = r_neighbour_elements.begin() ; i_neighbour_element != r_neighbour_elements.end() ; i_neighbour_element++)
	    {
	      if(i_neighbour_element->GetValue(IS_BODY) == 0)
	      {
		i_neighbour_element->SetValue(IS_BODY, BodyId);
		
		Element::GeometryType& p_geometry = i_neighbour_element->GetGeometry();
				
		for(unsigned int i = 0; i < p_geometry.size(); i++)
		  if(p_geometry[i].GetSolutionStepValue(IS_BODY) == 0)
		  {
		    p_geometry[i].GetSolutionStepValue(IS_BODY) = BodyId;
		    new_front_nodes.push_back(p_geometry(i));
		  }
	      }
	    }
	  }
	  front_nodes = new_front_nodes;
	}
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
	return "BodyRecognitionProcess";
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const{}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const{}
      
            
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
        
	ModelPart& mrModelPart;
	int mDimension;
        
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
      BodyRecognitionProcess& operator=(BodyRecognitionProcess const& rOther)
      {
	return *this;
      }

        
      ///@}    
        
    }; // Class BodyRecognitionProcess 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    BodyRecognitionProcess& rThis)
				    {
				      return rIStream;
				    }

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const BodyRecognitionProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_BODY_RECOGNITION_PROCESS_H_INCLUDED  defined 


