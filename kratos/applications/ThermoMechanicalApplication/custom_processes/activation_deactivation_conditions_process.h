//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2007-11-06 12:34:26 $
//   Revision:            $Revision: 1.4 $ 
//
//  this process save structural elements in a separate list 

#if !defined(KRATOS_ACTIVATION_DEACTIVATION_CONDITIONS_PROCESS )
#define  KRATOS_ACTIVATION_DEACTIVATION_CONDITIONS_PROCESS



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
#include "thermo_mechanical_application.h"
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

	class ActivationDeactivationConditionsProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of PushStructureProcess
		KRATOS_CLASS_POINTER_DEFINITION(ActivationDeactivationConditionsProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
// 		DuplicateInterfaceNodesCreateConditionsProcess()
// 		{
// 		}
	       ActivationDeactivationConditionsProcess(ModelPart& ThisModelPart,unsigned int max_id, const Matrix activation_table )
		:Process(), mr_model_part(ThisModelPart), mr_Nmax(max_id), mr_active_deactive_table(activation_table)
		{
		}		

		/// Destructor.
		virtual ~ActivationDeactivationConditionsProcess()
		{
		}


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
	            KRATOS_WATCH("INSIDE ActivationDeactivationConditionsProcess");
		   
		   //create Hash Matrix
		   Matrix hash_list;
		   unsigned int max_possible_combinations = mr_Nmax * (mr_Nmax + 1) ;
                   hash_list.resize(max_possible_combinations,2,false);
		   MakeHashList(hash_list);
/*		   KRATOS_WATCH(hash_list);*/
		    for (ModelPart::ConditionsContainerType::iterator cnd = mr_model_part.ConditionsBegin() ; 
					cnd != mr_model_part.ConditionsEnd() ; ++cnd)
		      {		   
			    unsigned int cond_ref_id = cnd->GetValue(REF_ID);
			    
			    if( hash_list(cond_ref_id - 1, 0) == 1.0)
			      cnd->SetValue(IS_INACTIVE ,  ! hash_list(cond_ref_id-1,1) );
			    
			   //free or fix tenperature just for teh exterior environment_contact conditions 
			    if(cond_ref_id <= mr_Nmax){
			      if( hash_list(cond_ref_id-1,1) == 1)
				(cnd->GetGeometry()[0]).Free(TEMPERATURE);
			      else
				(cnd->GetGeometry()[0]).Fix(TEMPERATURE);
			    }
		      }

		}//end of execute



			
		private:		
			  ModelPart& mr_model_part;
			  unsigned int mr_Nmax;
			  const Matrix mr_active_deactive_table;
			  
		//functions
		
		inline void MakeHashList(Matrix& hash_list)
		{
		   for(int jj = 0; jj<hash_list.size1(); ++jj ){
		     hash_list(jj,0) = 0.0;
		     hash_list(jj,1) = 0.0;
		   }
		   
		   int active_num = mr_active_deactive_table.size1();
		   
		   for(int ii = 0; ii< active_num; ++ii)
		   {
		     int first_nd = mr_active_deactive_table(ii,0);
		     int second_nd = mr_active_deactive_table(ii,1);
		     
		    if(first_nd == second_nd)
			  KRATOS_ERROR(std::logic_error,"inside activate_deactivate _ WRONG Activation TABLE","");			     
		     
		     if(first_nd != 0 && second_nd != 0)
		      {
			//mark for the 2-node conditions in the list
			int act_deact = first_nd*mr_Nmax + second_nd - 1;
			hash_list(act_deact,0) = 1.0;
			hash_list(act_deact,1) = mr_active_deactive_table(ii,2);	
			
			act_deact = second_nd*mr_Nmax + first_nd - 1;
			hash_list(act_deact,0) = 1.0;
			hash_list(act_deact,1) = mr_active_deactive_table(ii,2);		
			
			//mark the contrary for 1_node conditions in teh list
			act_deact = first_nd*mr_Nmax + first_nd - 1;
			hash_list(act_deact,0) = 1.0;
			hash_list(act_deact,1) = ! mr_active_deactive_table(ii,2);
			
			act_deact = second_nd*mr_Nmax + second_nd - 1;
			hash_list(act_deact,0) = 1.0;
			hash_list(act_deact,1) = ! mr_active_deactive_table(ii,2);		      
		      }
		      else if(first_nd != 0)
		      {
			int act_deact = first_nd - 1;
			hash_list(act_deact,0) = 1.0;
			hash_list(act_deact,1) = mr_active_deactive_table(ii,2);			
		      }
		      else
		      {
			int act_deact = second_nd - 1;
			hash_list(act_deact,0) = 1.0;
			hash_list(act_deact,1) = mr_active_deactive_table(ii,2);				
		      }
		     		     		     
		   }
		}
		
		 };//end of class
		  

}//end of namespace Kratos

#endif // KRATOS_DUPLICATE_INTERFACE_NODES_CREATE_CONDITIONS_PROCESS  defined 


