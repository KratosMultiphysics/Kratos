//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2007-11-06 12:34:26 $
//   Revision:            $Revision: 1.4 $ 
//
//  this process save structural elements in a separate list 

#if !defined(KRATOS_DUPLICATE_INTERFACE_NODES_CREATE_CONDITIONS_PROCESS )
#define  KRATOS_DUPLICATE_INTERFACE_NODES_CREATE_CONDITIONS_PROCESS



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

	class DuplicateInterfaceNodesCreateConditionsProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of PushStructureProcess
		KRATOS_CLASS_POINTER_DEFINITION(DuplicateInterfaceNodesCreateConditionsProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
// 		DuplicateInterfaceNodesCreateConditionsProcess()
// 		{
// 		}
	       DuplicateInterfaceNodesCreateConditionsProcess(ModelPart& ThisModelPart, char* heat_contact_condition, unsigned int max_id, const Matrix contact_table )
			:Process(), mr_model_part(ThisModelPart), mrCndHeat(KratosComponents<Condition>::Get(heat_contact_condition)), mr_Nmax(max_id), mr_contact_table(contact_table)
		{
		}		

		/// Destructor.
		virtual ~DuplicateInterfaceNodesCreateConditionsProcess()
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
	            KRATOS_WATCH("INSIDE DuplicateInterfaceNodesCreateConditionsProcess");
		   
		   //create Hash vector
		   Vector hash_list;
		   unsigned int max_possible_contact = mr_Nmax * (mr_Nmax + 1) ;
                   hash_list.resize(max_possible_contact);
		   MakeHashList(hash_list);
	   
		   //save Id to AUX_INDEX
		    for(ModelPart::NodesContainerType::iterator nd = mr_model_part.NodesBegin() ; 
					nd != mr_model_part.NodesEnd() ; ++nd)
			{
                                nd->GetValue(AUX_INDEX) = nd->Id();
			}		   
		   
		   //Reset IS_VISITED flag for elements
			for(ModelPart::ElementsContainerType::iterator im = mr_model_part.ElementsBegin() ; 
				im != mr_model_part.ElementsEnd() ; ++im)
			{
                                im->GetValue(IS_VISITED) = 0.0;
				
			}	
	
		   //List of created nodes 
		   typedef Node<3> PointType;
		   typedef std::vector<PointType::Pointer>           PointVector;

		   
 		   PointVector duplicated_nodes_list;
		   duplicated_nodes_list.reserve(mr_model_part.Nodes().size());
		   int nd_id = mr_model_part.Nodes().size();
		   int cnd_id = mr_model_part.Conditions().size();
		   const int step_data_size = mr_model_part.GetNodalSolutionStepDataSize();	
		   //Duplicate nodes
		    for(ModelPart::NodesContainerType::iterator nd = mr_model_part.NodesBegin() ; 
					nd != mr_model_part.NodesEnd() ; ++nd)
			{

		           WeakPointerVector< Element >& neighbor_els = nd->GetValue(NEIGHBOUR_ELEMENTS);
			   int num_material = 0;

			   //loop to find # of different material
			  for(WeakPointerVector< Element >::iterator ele_base = neighbor_els.begin(); ele_base!=neighbor_els.end(); ele_base++)
			  {	
			   double counted = ele_base->GetValue(IS_VISITED);
			   if(counted == 0.0)
			   {
			    num_material += 1; 
			    ele_base->GetValue(IS_VISITED) = 1.0;			    
			    int ele_base_mat = (ele_base->pGetProperties())->Id();

			     //loop over neighbors to find similar flag
			      for(WeakPointerVector< Element >::iterator ngh_ele = neighbor_els.begin(); ngh_ele!=neighbor_els.end(); ngh_ele++)
				if(ngh_ele->Id() != ele_base->Id())
				{
				  int ngh_ele_mat = (ngh_ele->pGetProperties())->Id();	
				  int is_contact = 1.0;
				  CheckHashList(ele_base_mat, ngh_ele_mat, hash_list, is_contact);
				  
				  if(ele_base_mat == ngh_ele_mat || is_contact == 0.0 )
				    ngh_ele->GetValue(IS_VISITED) = 1.0;
				}			    
			    }			   
			  } 
				
			    //reset is_visited
			    for(WeakPointerVector< Element >::iterator ele = neighbor_els.begin(); ele!=neighbor_els.end(); ele++)	
			      ele->GetValue(IS_VISITED) = 0.0;

			    //check if at least one material is created
			    if(num_material == 0)
				  KRATOS_ERROR(std::logic_error,"NO MATERIAL IS ASSIGNED","");			    
			    
			    
			    //Create nodes and replace 
			    if(num_material > 1)
			    {
				double x = nd->X();
				double y = nd->Y();
				double z = nd->Z();	
				
				//Exclude one material type because (num_material -1) copy is going to be created
				neighbor_els.begin()->GetValue(IS_VISITED) = 1.0;
			        int begin_elem_mat = (neighbor_els.begin()->pGetProperties())->Id();
			        for(WeakPointerVector< Element >::iterator ngh_ele = neighbor_els.begin()+1; ngh_ele!=neighbor_els.end(); ngh_ele++)
				  {
			            int ngh_elem_mat = (ngh_ele->pGetProperties())->Id();
				    int is_contact = 1.0;
				    CheckHashList(begin_elem_mat, ngh_elem_mat, hash_list, is_contact);
				    
				    if(begin_elem_mat == ngh_elem_mat || is_contact == 0.0)
				       ngh_ele->GetValue(IS_VISITED) = 1.0;

				  }

				//Duplicate for the rest of material types
				int kk = 0;
				for(WeakPointerVector< Element >::iterator ele_base = neighbor_els.begin(); ele_base!=neighbor_els.end(); ele_base++)
				{
				  kk+=1;
				  if(ele_base->GetValue(IS_VISITED) == 0.0)
				  {
				    //Create and interpolate
				    nd_id++;
			    
/*				    Node<3>::Pointer pnode = mr_model_part.CreateNewNode(nd_id,x,y,z);*/
                                    PointType::Pointer pnode = PointType::Pointer(new PointType(nd_id, x, y,z));				    
                                    AssignDataContainer( *(nd.base()), pnode, step_data_size);				    
				    
				    //Replace
				    ele_base->GetValue(IS_VISITED) = 1.0;					
				    int ele_base_mat = (ele_base->pGetProperties())->Id();
				    Geometry< Node<3> >& geom = ele_base->GetGeometry();
				    for(unsigned int iii = 0; iii<geom.size(); ++iii)
				      if( geom[iii].Id() == nd->GetValue(AUX_INDEX) )
					geom(iii) = pnode;						
			      				      
				      
				    //Look for similar Material   
				    for(WeakPointerVector< Element >::iterator ngh_ele = neighbor_els.begin()+kk; ngh_ele!=neighbor_els.end(); ngh_ele++)
				      if( ngh_ele->Id() != ele_base->Id() )
				      {	
					int ngh_ele_mat = (ngh_ele->pGetProperties())->Id();	
					int is_contact = 1.0;
					CheckHashList(ele_base_mat, ngh_ele_mat, hash_list, is_contact);				      				      
					
					if(ele_base_mat == ngh_ele_mat || is_contact == 0.0)
					{
					  ngh_ele->GetValue(IS_VISITED) = 1.0;	
					  Geometry< Node<3> >& geom = ngh_ele->GetGeometry();
					  for(unsigned int iii = 0; iii<geom.size(); ++iii)
					   if( geom[iii].Id() == nd->GetValue(AUX_INDEX) )
					     geom(iii) = pnode;						    
					}					
				      }	  //end of similar
				      
				      //add to the list
				      duplicated_nodes_list.push_back(pnode);

				   }
				 }//end of duplicate
	  
				}//End of Adding			
			         //reset is_visited
				for(WeakPointerVector< Element >::iterator ele = neighbor_els.begin(); ele!=neighbor_els.end(); ele++)	
				   ele->GetValue(IS_VISITED) = 0.0;
			      //add to model part reser is_visited for ne created node 
			    }//End Of loop over nodes
		    
			    //Add new nodes to the model part
		    	    for( PointVector::iterator it =  duplicated_nodes_list.begin(); it!=duplicated_nodes_list.end(); it++)
			      {
				mr_model_part.Nodes().push_back(*(it.base()));
			      }

			    //REset IS_VISITED flag and assign Material Poperty flag of element to its nodes
			    for(ModelPart::ElementsContainerType::iterator im = mr_model_part.ElementsBegin() ; 
				    im != mr_model_part.ElementsEnd() ; ++im)
			      {
				      im->GetValue(IS_VISITED) = 0.0;
				      
				      int ele_mat = (im->pGetProperties())->Id();	
				      Geometry< Node<3> >& geom = im->GetGeometry();
				      for(unsigned int iii = 0; iii<geom.size(); ++iii)	
					geom[iii].GetValue(NODE_PROPERTY_ID) = ele_mat;
			      }	
			      
			    for(ModelPart::NodesContainerType::iterator nd = mr_model_part.NodesBegin() ; 
						nd != mr_model_part.NodesEnd() ; ++nd)
				{
					nd->GetValue(IS_VISITED) = 0;
				}			      
			      
			 //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//
			 //   //Loop over elements to create the condition  //
			 //_________________________________________________//   
// 			Properties::Pointer cond_properties = mr_model_part.GetMesh().pGetProperties(1);
// 			Properties& temp_prop = mr_model_part.GetMesh().GetProperties(1);
// 			Properties::Pointer pCond_prop(new Properties(temp_prop));
			
	
			for(ModelPart::NodesContainerType::iterator nd = mr_model_part.NodesBegin() ; 
					    nd != mr_model_part.NodesEnd() ; ++nd)
			    {
				nd->GetValue(IS_VISITED) = 1.0;
			        int base_aux = nd->GetValue(AUX_INDEX);
				for(ModelPart::NodesContainerType::iterator sc_nd = mr_model_part.NodesBegin() ; 
						    sc_nd != mr_model_part.NodesEnd() ; ++sc_nd)
				    {	
				      int cp_aux = sc_nd->GetValue(AUX_INDEX);
				      if( base_aux==cp_aux && sc_nd->GetValue(IS_VISITED) == 0.0 )
				      {
				        //generate a face condition
					Condition::NodesArrayType temp;
				        temp.reserve(2);	
					
					//fill points
					temp.push_back(*(nd.base())); 
					temp.push_back(*(sc_nd.base()));	

					//Create new condition  and assign prop_id
					cnd_id++;
					unsigned int condition_ref_id;
					
					PairToId(nd->GetValue(NODE_PROPERTY_ID), sc_nd->GetValue(NODE_PROPERTY_ID), mr_Nmax, condition_ref_id);
					
/*					pCond_prop->SetId(condtion_prop_id);*/
					Condition::Pointer p_cond = mrCndHeat.Create(cnd_id, temp, mr_model_part.GetMesh().pGetProperties(1));	
				        p_cond->SetValue( IS_INACTIVE,0 );
				        p_cond->SetValue( REF_ID,condition_ref_id );					
					mr_model_part.Conditions().push_back(p_cond);	
					
					//create Ambient Temprerature conditions on FIRST duplicated nodes
					cnd_id++;
					Condition::NodesArrayType env_temp;
					env_temp.reserve(1);
					env_temp.push_back(*(nd.base())); 
					
					int prop_id = nd->GetValue(NODE_PROPERTY_ID);
					
					Condition::Pointer p_env_cond_first = KratosComponents<Condition>::Get("EnvironmentContact").Create(cnd_id, env_temp, mr_model_part.GetMesh().pGetProperties(prop_id));
					p_env_cond_first->SetValue( IS_INACTIVE,1 );
					p_env_cond_first->SetValue( REF_ID, (prop_id*mr_Nmax + prop_id) );
				        mr_model_part.Conditions().push_back(p_env_cond_first);
					
					
					//create Ambient Temprerature conditions on SECOND duplicated nodes					
					cnd_id++;
					env_temp.clear();
					env_temp.push_back(*(sc_nd.base())); 
					
					prop_id = sc_nd->GetValue(NODE_PROPERTY_ID);
					
					Condition::Pointer p_env_cond_second = KratosComponents<Condition>::Get("EnvironmentContact").Create(cnd_id, env_temp, mr_model_part.GetMesh().pGetProperties(prop_id));
					p_env_cond_second->SetValue( IS_INACTIVE,1 );
					p_env_cond_second->SetValue( REF_ID, (prop_id*mr_Nmax + prop_id) );						
					mr_model_part.Conditions().push_back(p_env_cond_second);
					
				      }
				    }

			    }//end of loop over nodes


			    KRATOS_WATCH(mr_model_part.Conditions().size());

			  //add nodal boundary condition for the Ambient Temperature
			for(ModelPart::NodesContainerType::iterator nd = mr_model_part.NodesBegin() ; 
					    nd != mr_model_part.NodesEnd() ; ++nd)
			    {
			      if(nd->FastGetSolutionStepValue(IS_BOUNDARY) == 1.0){
				//free DOF
// 				nd->Free(TEMPERATURE);
				
				cnd_id++;
				//generate a face condition
				Condition::NodesArrayType temp;
				temp.reserve(1);	
				
				//fill points
				temp.push_back(*(nd.base())); 	

// 				pCond_prop->SetId(nd->GetValue(NODE_PROPERTY_ID));
				int prop_id = nd->GetValue(NODE_PROPERTY_ID);
				Condition::Pointer p_cond = KratosComponents<Condition>::Get("EnvironmentContact").Create(cnd_id, temp, mr_model_part.GetMesh().pGetProperties(prop_id));
				p_cond->SetValue( IS_INACTIVE,1 );
				p_cond->SetValue( REF_ID, prop_id );					
				
				mr_model_part.Conditions().push_back(p_cond);
			      }
			      
			    }
			    KRATOS_WATCH(mr_model_part.Conditions().size());
			    
// 			for(ModelPart::ElementsContainerType::iterator Belem = mr_model_part.ElementsBegin(); Belem != mr_model_part.ElementsEnd(); ++Belem)
// 			{
// 		           WeakPointerVector< Element >& Belem_ngh = Belem->GetValue(NEIGHBOUR_ELEMENTS);
// 		           int Belem_mat = (Belem->pGetProperties())->Id();
// 			   Belem->GetValue(IS_VISITED) = 1.0;
// 		           Geometry< Node<3> >& Belem_geom = Belem->GetGeometry();			   
// 
// 			   //generate a face condition
// 			   Condition::NodesArrayType temp;
// 			   temp.reserve(2);			   
// 			   
// 			   //Check if the neighbors have different material
// 		           for(int iii=0; iii<3; ++iii)
// 			   {
// 			     if(Belem_ngh[iii].GetValue(IS_VISITED) == 0.0)
// 			       {
// 		                double ele_ngh_mat = (Belem_ngh[iii].GetProperties()).Id();
// 				if(ele_ngh_mat != Belem_mat)
// 				{
// 				  if(iii == 0) // 0 -> 1-2
// 				  {
// 
// 				    temp.push_back(Belem_geom(1)); 
// 				    temp.push_back(Belem_geom(2));	
// 				    
// 				    int first_nd_aux = Belem_geom[1].GetValue(AUX_INDEX);
// 				    int second_nd_aux = Belem_geom[2].GetValue(AUX_INDEX);
// 				    int third = 4;
// 				    int forth = 4;
// 
// 		                    Geometry< Node<3> >& Belem_ngh_geom = Belem_ngh[iii].GetGeometry();	
// 				    for(int kk=0;kk<3;++kk)
// 				    {
// 				      int ngh_nd_aux = Belem_ngh_geom[kk].GetValue(AUX_INDEX);
// 				      if(ngh_nd_aux == first_nd_aux)
// 				        forth = kk;
// 				      if(ngh_nd_aux == second_nd_aux)  
// 					third = kk;				      
// 				    }
// 				    if(third == 4 || forth == 4)
// 					  KRATOS_ERROR(std::logic_error,"one of the condition nodes are misssing","");					    
// 				    
// 				    temp.push_back(Belem_ngh_geom(forth));
// 				    temp.push_back(Belem_ngh_geom(third));				    
// 				  }
// 				  
// 				  if(iii == 1) // 1 -> 2-0
// 				  {
// 				    temp.push_back(Belem_geom(2)); 
// 				    temp.push_back(Belem_geom(0));	
// 				    
// 				    int first_nd_aux = Belem_geom[2].GetValue(AUX_INDEX);
// 				    int second_nd_aux = Belem_geom[0].GetValue(AUX_INDEX);
// 				    int third = 4;
// 				    int forth = 4;
// 
// 		                    Geometry< Node<3> >& Belem_ngh_geom = Belem_ngh[iii].GetGeometry();	
// 				    for(int kk=0;kk<3;++kk)
// 				    {
// 				      int ngh_nd_aux = Belem_ngh_geom[kk].GetValue(AUX_INDEX);
// 				      if(ngh_nd_aux == first_nd_aux)
// 				        forth = kk;
// 				      if(ngh_nd_aux == second_nd_aux)  
// 					third = kk;				      
// 				    }
// 				    if(third == 4 || forth == 4)
// 					  KRATOS_ERROR(std::logic_error,"one of the condition nodes are misssing","");					    
// 				    
// 				    temp.push_back(Belem_ngh_geom(forth));
// 				    temp.push_back(Belem_ngh_geom(third));				    
// 				  }
// 				  
// 				  if(iii == 2) // 2 -> 0-1
// 				  {
// 				    temp.push_back(Belem_geom(0)); 
// 				    temp.push_back(Belem_geom(1));	
// 				    
// 				    int first_nd_aux = Belem_geom[0].GetValue(AUX_INDEX);
// 				    int second_nd_aux = Belem_geom[1].GetValue(AUX_INDEX);
// 				    int third = 4;
// 				    int forth = 4;
// 
// 		                    Geometry< Node<3> >& Belem_ngh_geom = Belem_ngh[iii].GetGeometry();	
// 				    for(int kk=0;kk<3;++kk)
// 				    {
// 				      int ngh_nd_aux = Belem_ngh_geom[kk].GetValue(AUX_INDEX);
// 				      if(ngh_nd_aux == first_nd_aux)
// 				        forth = kk;
// 				      if(ngh_nd_aux == second_nd_aux)  
// 					third = kk;				      
// 				    }
// 				    if(third == 4 || forth == 4)
// 					  KRATOS_ERROR(std::logic_error,"one of the condition nodes are misssing","");					    
// 				    
// 				    temp.push_back(Belem_ngh_geom(forth));
// 				    temp.push_back(Belem_ngh_geom(third));				    
// 				  }	
// 				  
// 				//Create new condition  
// 				cnd_id++;
//                                 Condition::Pointer p_cond = mrCndHeat.Create(cnd_id, temp, properties);				
// 			        mr_model_part.Conditions().push_back(p_cond);
// 				
// 				//Assign MATERIAL PROPERTy of teh condition
// 				  
// 				}//Different material
// 			       }//Not visited
// 			     
// 			   }//Elemental neighbors
// 			   
// 			   
// 			  
// 			}// end of loop over elemenets
				
			}//end of execute





		void PairToId(unsigned int ii, unsigned int jj,unsigned int N_max, unsigned int& cond_prop_id)
		{
		  if( ii > N_max or jj > N_max )
		    	KRATOS_ERROR(std::logic_error," Beginning or end id is bigger than Max_Id","");	
		  if( ii== jj)
		    KRATOS_WATCH("Nodes of the created condition have the same NODE_PROPERTY_ID");
		  
		  cond_prop_id = ii*N_max + jj; 
		}
		
		void IdToPair(unsigned int cond_prop_id , unsigned int N_max, int& init_prop_id, int& end_prop_id)
		{
		  if( N_max == 0 or cond_prop_id < N_max)
		    	KRATOS_ERROR(std::logic_error,"Max_Id is zero or Condition_id is less than N_max","");	
		    
		  init_prop_id = cond_prop_id / N_max ;
		  end_prop_id = cond_prop_id % N_max;
		}
			
		private:		
			  ModelPart& mr_model_part;
			  Condition const& mrCndHeat;	
			  unsigned int mr_Nmax;
                          const Matrix mr_contact_table;			  
			  
		//functions
		inline void AssignDataContainer(Node<3>::Pointer pMaster_nd, Node<3>::Pointer pCreated_nd, const int data_size)
		{
			unsigned int buffer_size = pMaster_nd->GetBufferSize();		  
			pCreated_nd->SetSolutionStepVariablesList(pMaster_nd->pGetVariablesList());
			pCreated_nd->SetBufferSize(buffer_size);		  
		  
			 //Get and copy Dof 
			Node<3>::DofsContainerType& master_dofs = pMaster_nd->GetDofs();
			for(Node<3>::DofsContainerType::iterator iii = master_dofs.begin();    iii != master_dofs.end(); iii++)
					{
						Node<3>::DofType& rDof = *iii;
						Node<3>::DofType::Pointer p_new_dof = pCreated_nd->pAddDof( rDof );
						
// 						(p_new_dof)->FreeDof();
					}
			 
			for(unsigned int step = 0; step<buffer_size; step++)
			{					
						//getting the data of the solution step
				double* created_step_data = (pCreated_nd)->SolutionStepData().Data(step);
				double* master_step_data = (pMaster_nd)->SolutionStepData().Data(step);				
											
				
				//copying this data in the position of the vector we are interested in
				for( int j= 0; j<data_size; j++)
				{ 
					created_step_data[j] = master_step_data[j];								
				}
								
			}

				    pCreated_nd->GetValue(AUX_INDEX) = pMaster_nd->GetValue(AUX_INDEX);	
		} //end of inline
		
		inline void MakeHashList(Vector& hash_list)
		{
		   for(int jj = 0; jj>hash_list.size(); ++jj )
		     hash_list[jj] = 0.0;
		   
		   int contact_num = mr_contact_table.size1();
		   
		   for(int ii = 0; ii< contact_num; ++ii)
		   {
		    if(mr_contact_table(ii,0) == 0 || mr_contact_table(ii,1) == 0)
			  KRATOS_ERROR(std::logic_error,"inside duplicate _ WRONG CONTACT TABLE","");				     
		     
		     int contact_int = mr_contact_table(ii,0)*mr_Nmax + mr_contact_table(ii,1) - 1;
		     hash_list[contact_int] = 1.0;
		     contact_int = mr_contact_table(ii,1)*mr_Nmax + mr_contact_table(ii,0) - 1;
		     hash_list[contact_int] = 1.0;		     
		   }
		}
		
		inline void CheckHashList(int begin_elem_mat,int ngh_elem_mat,Vector hash_list, int& is_contact)
		{
		  int contact_int = begin_elem_mat*mr_Nmax + ngh_elem_mat - 1;
		  
		  if( hash_list[contact_int] ==  0.0)
		    is_contact = 0.0;
		}
		 };//end of class
		  

}//end of namespace Kratos

#endif // KRATOS_DUPLICATE_INTERFACE_NODES_CREATE_CONDITIONS_PROCESS  defined 


