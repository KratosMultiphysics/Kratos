//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-10-02 10:47:21 $
//   Revision:            $Revision: 1.8 $ 
//
//  this process save structural elements in a separate list 

#if !defined(KRATOS_SAVE_FLAG_MODEL_PART_PROCESS_INCLUDED )
#define  KRATOS_SAVE_FLAG_MODEL_PART_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"


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
		Update the PRESSURE_FORCE on the nodes

		
	*/

	class SaveFlagModelPartProcess 
		: public Process
	{
	   public:
	      SaveFlagModelPartProcess(ModelPart& fluid_model_part,ModelPart& changing_fluid_model_part, int domain_size, Kratos::Variable<double>& flag, double value )
			:Process(), fluid_model(fluid_model_part), changing_model(changing_fluid_model_part), domain_size(domain_size), f_value(value), m_flag(flag) 
		{
		}

		///@name Type Definitions
		///@{


		/// Destructor.
		virtual ~SaveFlagModelPartProcess()
		{
		}
		
	      void operator()()
		{
		  Execute();
		}
      

	      ///@}
	      ///@name Operators 
	      ///@{



		///@}
		///@name Operators 
		///@{

//		void operator()()
//		{
//			SaveStructure();
//		}


		///@}
		///@name Operations
		///@{
	   virtual void Execute()
		 {
// 		void  SaveStructure(ModelPart& fluid_model_part, ModelPart& changing_fluid_model_part,  int domain_size, Kratos::Variable<int>& flag, int value)
// 		{
			KRATOS_TRY
			
			changing_model.Elements().clear();
			changing_model.Nodes().clear();			
			
                        ModelPart::ElementsContainerType fix_model_element = changing_model.Elements() ;
                        ModelPart::NodesContainerType fix_model_nodes = changing_model.Nodes() ;			
			fix_model_element.clear();	
			fix_model_nodes.clear();			
/*			fix_model_element.resreve(fluid_model_part.Elements().size())*/
			//number of structure nodes
// 			KRATOS_WATCH("SAVING STRUCTURE**********")
//                        KRATOS_WATCH(fluid_model.Elements().size());
			unsigned int ele_id=1;
			int cnt = 0;
		    for(ModelPart::ElementsContainerType::iterator im = fluid_model.ElementsBegin() ; 
					im != fluid_model.ElementsEnd() ; ++im)
			{		
				//PointerVector<Element> struct_elements_list;
				//check number of structure nodes
				unsigned int elem_size = im->GetGeometry().size();
				int n_flag = 0;
				double nd_val = 0.0;
				
				if( elem_size == static_cast<unsigned int>(domain_size +1) )
				{
				  for(unsigned int i=0;i<elem_size;i++)
				  {
					  nd_val =  im->GetGeometry()[i].FastGetSolutionStepValue(m_flag) ;
					  if(nd_val == f_value || nd_val == f_value + 1.0 )
					     n_flag += 1;					  
				  }
				  
				  if( n_flag == (domain_size +1) ){
				           im->GetValue(IS_WATER_ELEMENT) = 1.0;
					   fix_model_element.push_back(*(im.base()));	}			
				  else{
				            im->SetId(ele_id);
					    ele_id++;
				            changing_model.Elements().push_back(*(im.base()));}
				  
				}
				else{
				            cnt++;
					    im->SetId(ele_id);
					    ele_id++;
					    changing_model.Elements().push_back(*(im.base()));}
				 
	
			}
/*			KRATOS_WATCH(cnt);
			KRATOS_WATCH(changing_model.Elements().size());
			KRATOS_WATCH(fix_model_element.size());	*/		
			changing_model.Elements().Sort();

			ele_id += 1000;
//                     fluid_model.Elements() = fix_model_element;
			fluid_model.Elements().clear();
		    for(ModelPart::ElementsContainerType::iterator im = fix_model_element.begin() ; 
					im != fix_model_element.end() ; ++im)
			{
			  im->SetId(ele_id);
			  ele_id++;			  
			  fluid_model.Elements().push_back(*(im.base()));
			  
			}
			fluid_model.Elements().Sort();
// 			KRATOS_WATCH(fluid_model.Elements().size());				

		    unsigned int nd_id=1;		    
		    for(ModelPart::NodesContainerType::iterator nd = fluid_model.NodesBegin() ; 
					nd != fluid_model.NodesEnd() ; ++nd)
			{
			     double nd_val = 0.0;

			     nd_val =  nd->FastGetSolutionStepValue(m_flag) ;	
			     if( nd_val == f_value)
			       fix_model_nodes.push_back(*(nd.base()));
			     else{
				   nd->SetId(nd_id);
			           nd_id++;			       
			           changing_model.Nodes().push_back(*(nd.base()));}			       
			  
			}
			changing_model.Nodes().Sort();

			nd_id += 1000;			
//                       fluid_model.Nodes() = fix_model_nodes;		
			fluid_model.Nodes().clear();
		    for(ModelPart::NodesContainerType::iterator nd = fix_model_nodes.begin() ; 
					nd != fix_model_nodes.end() ; ++nd)
			{
		          nd->SetId(nd_id);
			  nd_id++;			  
			  fluid_model.Nodes().push_back(*(nd.base()));
			  
			}		      
			fluid_model.Nodes().Sort();
		    
		    	fix_model_element.clear();	
			fix_model_nodes.clear();

		KRATOS_CATCH("")
		}
		

		private:
			ModelPart& fluid_model;
			ModelPart& changing_model;
		         int domain_size;
			 double f_value;
			Kratos::Variable<double>& m_flag;			

	};
}  // namespace Kratos.

#endif // KRATOS_SAVE_STRUCTURE_MODEL_PART_PROCESS_INCLUDED  defined 


