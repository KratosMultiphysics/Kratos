//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_BOUNDARY_SKIN_BUILD_PROCESS_H_INCLUDED )
#define  KRATOS_BOUNDARY_SKIN_BUILD_PROCESS_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes
#include <boost/timer.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "custom_utilities/boundary_normals_calculation_utilities.hpp"
#include "pfem_solid_mechanics_application.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
	typedef  ModelPart::NodesContainerType NodesContainerType;
	typedef  ModelPart::ElementsContainerType ElementsContainerType;
	typedef  ModelPart::ConditionsContainerType ConditionsContainerType;

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
	class BoundarySkinBuildProcess
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of BoundarySkinBuildProcess
		KRATOS_CLASS_POINTER_DEFINITION( BoundarySkinBuildProcess );

		///@}
		///@name Life Cycle
		///@{

		/// Default constructor.

		// BoundarySkinBuildProcess(ModelPart& model_part,
		// 		     char * ReferenceConditionName,
		// 		     unsigned int dim,
		// 		     unsigned int preserve)
		// 	: mr_model_part(model_part), mr_reference_condition(KratosComponents<Condition>::Get("ReferenceConditionName"))
		// { 
		// 	std::cout<<" Reference Condition "<<mr_reference_condition<<std::endl;
		// 	m_preserve = preserve;
		// 	mdim=dim;
		// }
  
		//second constructor
		BoundarySkinBuildProcess(ModelPart& model_part,
				     unsigned int dim = 2,
				     unsigned int preserve = 1)
			: mr_model_part(model_part)
		{ 
			m_preserve = preserve;
			mdim=dim;
		}

		/// Destructor.
		virtual ~BoundarySkinBuildProcess()
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
			bool success=false;

			boost::timer auxiliary;
			
			unsigned int start=0;
			unsigned int NumberOfMeshes=mr_model_part.NumberOfMeshes();
			//unsigned int NumberOfMeshes=1;
			if(NumberOfMeshes>1) 
			  start=1;
			  
			for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
			  {
			    
			    //success=SkinSearch(MeshId);
			    success=UniqueSkinConditions(MeshId);
			    
			    if(!success)
			      {
				std::cout<<" ERROR:  Skin Search FAILED on mesh : ["<<MeshId<<"] "<<std::endl;
			      }
			    else
			      {
				std::cout<<"Skin Search on mesh ["<<MeshId<<"] performed in time = "<<auxiliary.elapsed()<<std::endl;
				//PrintSkin(MeshId);
			      }
			  }

			if(NumberOfMeshes>1) 
			  SetGlobalConditions();

			//ComputeBoundaryNormals BoundUtils;
			BoundaryNormalsCalculationUtilities BoundaryComputation;
			BoundaryComputation.CalculateBoundaryNormals(mr_model_part,2);
			
			std::cout<<" Boundary Normals Computed "<<std::endl;

		};


		void SetGlobalConditions()
		{

		  
		  std::cout<<" [ OLD TOTAL CONDITIONS: "<<mr_model_part.NumberOfConditions()<<"] "<<std::endl;

		  //contact conditions are located on Mesh_0
		  ModelPart::ConditionsContainerType KeepConditions;

		  unsigned int condId=1;
		  unsigned int start=0;
		  
		  unsigned int NumberOfMeshes=mr_model_part.NumberOfMeshes();
		  if(NumberOfMeshes>1) 
		    start=1;

		  for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
		    {
		      for(ModelPart::ConditionsContainerType::iterator i_cond = mr_model_part.ConditionsBegin(MeshId) ; i_cond != mr_model_part.ConditionsEnd(MeshId) ; i_cond++)
			{
			  // i_cond->PrintInfo(std::cout);
			  // std::cout<<" -- "<<std::endl;
			  KeepConditions.push_back(*(i_cond.base()));
			  KeepConditions.back().SetId(condId);
			  condId+=1;

			  // KeepConditions.back().PrintInfo(std::cout);
			  // std::cout<<std::endl;

			}
		    }


		  for(ModelPart::ConditionsContainerType::iterator i_cond = mr_model_part.ConditionsBegin(); i_cond!= mr_model_part.ConditionsEnd(); i_cond++)
		    {
		      if(i_cond->Is(CONTACT)){
			KeepConditions.push_back(*(i_cond.base()));
			KeepConditions.back().SetId(condId);
			condId+=1;

			//std::cout<<" -- "<<std::endl;
			//KeepConditions.back().PrintInfo(std::cout);
			//std::cout<<std::endl;
			  
		      }
		      
		    }
      
		  mr_model_part.Conditions().swap(KeepConditions);

		  std::cout<<" [ NEW TOTAL CONDITIONS: "<<mr_model_part.NumberOfConditions()<<"] "<<std::endl;

		}

		void ClearConditions(ModelPart::IndexType MeshId=0)
		{
		  
		  std::cout<<" [ PREVIOUS CONDITONS: "<<mr_model_part.NumberOfConditions(MeshId)<<"] "<<std::endl;
		  //clone previous conditions
		  //m_conditions = mr_model_part.Conditions(MeshId);
		  
		  //clear conditions
		  //mr_model_part.Conditions(MeshId).clear();

		  //swap conditions
		  m_conditions.clear();
		  m_conditions.reserve(mr_model_part.Conditions(MeshId).size());
		  m_conditions.swap(mr_model_part.Conditions(MeshId));
		}


		bool SearchConditionMasters(ModelPart::IndexType MeshId=0)
		{

		        unsigned int counter = 0;
			bool found=false;

			for(ModelPart::ConditionsContainerType::iterator ic = mr_model_part.ConditionsBegin(MeshId); ic != mr_model_part.ConditionsEnd(MeshId); ic++)
			  {
	  
			    //std::cout<<" Condition ("<<ic->Id()<<") : ME="<<ic->GetValue(MASTER_ELEMENTS)[0].Id()<<", MN= "<<ic->GetValue(MASTER_NODES)[0].Id()<<std::endl;

			    //********************************************************************

			    boost::numeric::ublas::matrix<unsigned int> lpofa; //points that define the faces

			    Geometry< Node<3> >& rConditionGeom = ic->GetGeometry();
			    unsigned int size=rConditionGeom.size();
			    
			    bool perform_search = true;
			    for(unsigned int i=0; i<size; i++)
			      if( rConditionGeom[i].SolutionStepsDataHas(RIGID_WALL) )
				perform_search = false;			   		     

			    if( size != 2 ) 
			      perform_search = false;

			    //********************************************************************
			    found=false;

			    if( perform_search )
			      {

				WeakPointerVector<Element >& rE1 = rConditionGeom[0].GetValue(NEIGHBOUR_ELEMENTS);				    
				WeakPointerVector<Element >& rE2 = rConditionGeom[1].GetValue(NEIGHBOUR_ELEMENTS);

				for(WeakPointerVector< Element >::iterator ie = rE1.begin(); ie!=rE1.end(); ie++)
				  {
				    for(WeakPointerVector< Element >::iterator ne = rE2.begin(); ne!=rE2.end(); ne++)
				      {

					if (ne->Id() == ie->Id() && !found)
					  {
					    WeakPointerVector< Element > MasterElements;
					    MasterElements.push_back(Element::WeakPointer( *(ie.base()) ) );
					    ic->SetValue(MASTER_ELEMENTS,MasterElements);
					 
					    Geometry< Node<3> >& rElementGeom = ie->GetGeometry();

					    rElementGeom.NodesInFaces(lpofa);

					    int node = 0;
					    for (unsigned int i=0; i<rElementGeom.size(); i++)
					      {
						if( (   rConditionGeom[0].Id() == rElementGeom[lpofa(1,i)].Id() 
							&& rConditionGeom[1].Id() == rElementGeom[lpofa(2,i)].Id() ) || 
						    (   rConditionGeom[0].Id() == rElementGeom[lpofa(2,i)].Id() 
							&& rConditionGeom[1].Id() == rElementGeom[lpofa(1,i)].Id() ) )
						  {
						    node=i;
						    found = true;
						    break;
						  }
					      }
						
					    if(found){
					      WeakPointerVector< Node<3> > MasterNodes;
					      MasterNodes.push_back( Node<3>::WeakPointer( rElementGeom(lpofa(0,node)) ) );
					      ic->SetValue(MASTER_NODES,MasterNodes);
					    }
					    else{						 
					      std::cout<<" MASTER_NODE not FOUND : something is wrong "<<std::endl;			  
					    }

					  }
				      }
				  }
														  
			      }

			    //********************************************************************

			    //std::cout<<" After Condition ("<<ic->Id()<<") : ME="<<ic->GetValue(MASTER_ELEMENTS)[0].Id()<<", MN= "<<ic->GetValue(MASTER_NODES)[0].Id()<<std::endl;

			    if(found)
			      counter++;
		    
			  }

			double totalcond=0;
			if(mr_model_part.Conditions(MeshId).size()>0)
			  totalcond = mr_model_part.Conditions(MeshId).size();
			  

			if(counter == totalcond){
			  std::cout<<"   Condition Masters (mesh "<<MeshId<<"): LOCATED ["<<counter<<"]"<<std::endl;
			  found=true;
			}
			else{
			  std::cout<<"   Condition Masters (mesh "<<MeshId<<"): not LOCATED ["<<counter-totalcond<<"]"<<std::endl;
			  found=false;
			}
			
			return found;

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
			return "BoundarySkinBuildProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "BoundarySkinBuildProcess";
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


		///@}
		///@name Member Variables
		///@{
		ModelPart& mr_model_part;
		//Condition& mr_reference_condition;
		ConditionsContainerType m_conditions;
		unsigned int m_preserve;
		unsigned int mdim;


		///@}
		///@name Private Operators
		///@{


		///@}
		///@name Private Operations
		///@{


		void PrintSkin (ModelPart::IndexType MeshId=0)
		{
			//PRINT SKIN:
			std::cout<<" CONDITIONS: geometry nodes ("<<mr_model_part.Conditions(MeshId).size()<<")"<<std::endl;

			ConditionsContainerType& rCond = mr_model_part.Conditions(MeshId);
			for(ConditionsContainerType::iterator ic = rCond.begin(); ic!= rCond.end(); ic++)
			{
			
				Geometry< Node<3> >& rConditionGeom = ic->GetGeometry();
				std::cout<<"["<<ic->Id()<<"]:"<<std::endl;
				//ic->PrintInfo(std::cout);
				std::cout<<"( ";
				for(unsigned int i = 0; i < rConditionGeom.size(); i++)
				{
					std::cout<< rConditionGeom[i].Id()<<", ";
				}
				std::cout<<" ): ";

				ic->GetValue(MASTER_ELEMENTS)[0].PrintInfo(std::cout);
							
				std::cout<<std::endl;

			}
			std::cout<<std::endl;

		}

		bool FindCondition(Geometry< Node<3> >& rConditionGeom,Geometry< Node<3> >& rGeom , boost::numeric::ublas::matrix<unsigned int>& lpofa, unsigned int& i)
		{
			//2D edges:
			if(mdim==2)
			{
				if( (   rConditionGeom[0].Id() == rGeom[lpofa(1,i)].Id() 
					&& rConditionGeom[1].Id() == rGeom[lpofa(2,i)].Id() ) || 
				    (   rConditionGeom[0].Id() == rGeom[lpofa(2,i)].Id() 
					&& rConditionGeom[1].Id() == rGeom[lpofa(1,i)].Id() ) )
				{
					return true;
				}
			}
      
			//3D faces:
			if(mdim==3)
			{
				if( (   rConditionGeom[0].Id() == rGeom[lpofa(1,i)].Id() 
					&& rConditionGeom[1].Id() == rGeom[lpofa(2,i)].Id()
					&& rConditionGeom[2].Id() == rGeom[lpofa(3,i)].Id() ) || 
				    (   rConditionGeom[0].Id() == rGeom[lpofa(3,i)].Id() 
					&& rConditionGeom[1].Id() == rGeom[lpofa(1,i)].Id()
					&& rConditionGeom[2].Id() == rGeom[lpofa(2,i)].Id() ) ||
				    (   rConditionGeom[0].Id() == rGeom[lpofa(2,i)].Id() 
					&& rConditionGeom[1].Id() == rGeom[lpofa(3,i)].Id()
					&& rConditionGeom[2].Id() == rGeom[lpofa(1,i)].Id() ) )
				{
					return true;
				}
	  
			}

			return false;
  
		}


	        bool FindConditionDomain(Geometry< Node<3> >& rConditionGeom, unsigned int& MeshId)
		{

		        //check if the conditions belongs to the MeshId checking the nodes Id
		        for(unsigned int i=0; i<rConditionGeom.size(); i++)
			  {
			    for(ModelPart::NodesContainerType::const_iterator in = mr_model_part.NodesBegin(MeshId); in!=mr_model_part.NodesEnd(MeshId); in++)
			      {			
				if( rConditionGeom[i].Id() == in->Id() )
				  return true;
			      }
			  }

			return false;
  		}


		bool SkinSearch(ModelPart::IndexType MeshId=0)
		{
			//properties to be used in the generation
  		        int number_properties = mr_model_part.NumberOfProperties();
			Properties::Pointer properties = mr_model_part.GetMesh().pGetProperties(number_properties-1);
			
			//Clear Previous Skin (preserve the conditions)
			ClearConditions(MeshId);

			//set consecutive ids in the mesh conditions
			unsigned int condId=1;
			for(ModelPart::ConditionsContainerType::iterator ic = m_conditions.begin(); ic!= m_conditions.end(); ic++)
			  {
			    ic->SetId(condId);
			    condId++;
			  }

			//control the previous mesh conditions
			std::vector<int> PreservedConditions( m_conditions.size() );
			std::fill( PreservedConditions.begin(), PreservedConditions.end(), 0 );
		
			//reset the boundary flag
			for(ModelPart::NodesContainerType::const_iterator in = mr_model_part.NodesBegin(MeshId); in!=mr_model_part.NodesEnd(MeshId); in++)
			{
				in->Reset(BOUNDARY);
			}

			//filling the elemental neighbours list (from now on the elements list can not change)
			//ModelPart::ElementsContainerType::iterator elements_end = mr_model_part.Elements().end();

			//mr_model_part.Elements().Unique();

			//now the boundary faces
			std::cout<<" mr_model_part.Conditions("<<MeshId<<").size(): "<<mr_model_part.Conditions(MeshId).size()<<std::endl;
			std::cout<<" m_conditions.size(): "<< m_conditions.size()<<std::endl;
			std::cout<<" Preserved Conditions "<<PreservedConditions.size()<<std::endl;
			int id=0;
			for(ModelPart::ElementsContainerType::iterator ie = mr_model_part.ElementsBegin(MeshId); ie != mr_model_part.ElementsEnd(MeshId); ie++)
			{
	  
				/*each face is opposite to the corresponding node number so in 2D
				  0 ----- 1 2
				  1 ----- 2 0
				  2 ----- 0 1
				*/

				//finding boundaries and creating the "skin"
				//
				//********************************************************************

				boost::numeric::ublas::matrix<unsigned int> lpofa; //points that define the faces

				Geometry< Node<3> >& rGeom = ie->GetGeometry();
				WeakPointerVector<Element >& rE = ie->GetValue(NEIGHBOUR_ELEMENTS);
				rGeom.NodesInFaces(lpofa);   
	    
				//unsigned int size=rGeom.size();

				//loop on neighbour elements of an element
				unsigned int i=0;
				for(WeakPointerVector< Element >::iterator ne = rE.begin(); ne!=rE.end(); ne++)
				{
					if (ne->Id() == ie->Id())
					{
						//if no neighnour is present => the face is free surface
						int sizei=lpofa.size2();
						for(int j=1; j<sizei; j++)
						{
							rGeom[lpofa(j,i)].Set(BOUNDARY);
						}
	
						//Get the correct ReferenceCondition
						Condition::Pointer pBoundaryCondition;
						bool condition_found = false;
						bool point_condition = false;

						
						for(ModelPart::ConditionsContainerType::iterator ic = m_conditions.begin(); ic!= m_conditions.end(); ic++)
						{
							Geometry< Node<3> >& rConditionGeom = ic->GetGeometry();
							
							if( FindCondition(rConditionGeom,rGeom,lpofa,i) ){
							    pBoundaryCondition = (*(ic.base())); //accessing boost::shared_ptr  get() to obtain the raw pointer
							    PreservedConditions[ic->Id()-1] += 1;
							    condition_found = true;
							    if( rConditionGeom.PointsNumber() == 1 )
  								point_condition = true;						
							    break;
							}
																	
						}

						if( !point_condition ){

							//Generate condition
							Condition::NodesArrayType face;
							face.reserve(sizei);
							for(int j=1; j<sizei; j++)
							{
								face.push_back(rGeom(lpofa(j,i)));
							}
				    
							id +=1;

							Condition::GeometryType::Pointer vertices = Condition::GeometryType::Pointer(new Geometry< Node<3> >(face) );
							Condition::Pointer p_cond;
							if(condition_found)
							{
								p_cond = pBoundaryCondition->Clone(id, face);
							}
							else
							{
								//							pBoundaryCondition = new Condition(id,vertices);
								p_cond =  Condition::Pointer(new Condition(id,vertices,properties) ); //pBoundaryCondition->Create(id, vertices, properties);
								//	delete pBoundaryCondition;

								//p_cond = mr_reference_condition.Create(id, face, properties);
							}
 
							//usually one MasterElement and one MasterNode in 2D in 3D can be more than one

							//p_cond->GetValue(MASTER_ELEMENTS).push_back( Element::WeakPointer( *(ie.base()) ) );
							WeakPointerVector< Element >& MasterElements = p_cond->GetValue(MASTER_ELEMENTS);
							MasterElements.push_back( Element::WeakPointer( *(ie.base()) ) );
							p_cond->SetValue(MASTER_ELEMENTS,MasterElements);

							//p_cond->GetValue(MASTER_NODES).push_back( Node<3>::WeakPointer( rGeom(lpofa(0,i)) ) );						
							WeakPointerVector< Node<3> >& MasterNodes = p_cond->GetValue(MASTER_NODES);
							MasterNodes.push_back( Node<3>::WeakPointer( rGeom(lpofa(0,i)) ) );
							p_cond->SetValue(MASTER_NODES,MasterNodes);


							//std::cout<<" Set preserved condition found "<<id<<std::endl;
							mr_model_part.Conditions(MeshId).push_back(p_cond);

						}

					}
					
					i+=1;
				}

			}



			//all previous conditions have to be added
			for(ModelPart::ConditionsContainerType::iterator ic = m_conditions.begin(); ic!= m_conditions.end(); ic++)
			{		    
			    if( PreservedConditions[ic->Id()-1] == 0 ){

				Geometry< Node<3> >& rGeom = ic->GetGeometry();
				
				if( FindConditionDomain(rGeom, MeshId) ){

				  Condition::NodesArrayType face;

				  face.reserve(rGeom.size() );

				  for(unsigned int j=0; j<rGeom.size(); j++)
				    {
				      face.push_back(rGeom(j));
				    }

				  PreservedConditions[ic->Id()-1] += 1;

				  id +=1;

				  mr_model_part.Conditions(MeshId).push_back(ic->Clone(id,face));

				}
				//std::cout<<" Set preserved condition not found "<<ic->Id()<<std::endl;

			    }
			}


			//control if previous conditions have been assigned
			bool all_assigned = true;
			for(unsigned int i=0; i<PreservedConditions.size(); i++)
			{
			    if( PreservedConditions[i] == 0 )
				all_assigned = false;
			}

			std::cout<<" [ NEW_CONDITIONS: "<<mr_model_part.NumberOfConditions(MeshId)<<"] "<<std::endl;

			if(all_assigned == true)
			    std::cout<<"Boundary Conditions RELOCATED "<<std::endl;
			else
			    std::cout<<"Boundary Conditions NOT relocated "<<std::endl;

	
			return true;
		}



		bool UniqueSkinConditions(ModelPart::IndexType MeshId=0)
		{
			//properties to be used in the generation
  		        int number_properties = mr_model_part.NumberOfProperties();
			Properties::Pointer properties = mr_model_part.GetMesh().pGetProperties(number_properties-1);
			
			//Clear Previous Skin (preserve the conditions)
			ClearConditions(MeshId);

			//set consecutive ids in the mesh conditions
			unsigned int condId=1;
			for(ModelPart::ConditionsContainerType::iterator ic = m_conditions.begin(); ic!= m_conditions.end(); ic++)
			  {
			    ic->SetId(condId);
			    condId++;
			  }

			//control the previous mesh conditions
			std::vector<int> PreservedConditions( m_conditions.size() );
			std::fill( PreservedConditions.begin(), PreservedConditions.end(), 0 );
		
			//reset the boundary flag
			for(ModelPart::NodesContainerType::const_iterator in = mr_model_part.NodesBegin(MeshId); in!=mr_model_part.NodesEnd(MeshId); in++)
			  {
			    in->Reset(BOUNDARY);
			  }

			//filling the elemental neighbours list (from now on the elements list can not change)
			//ModelPart::ElementsContainerType::iterator elements_end = mr_model_part.Elements().end();

			//mr_model_part.Elements().Unique();

			//now the boundary faces
			std::cout<<" mr_model_part.Conditions("<<MeshId<<").size(): "<<mr_model_part.Conditions(MeshId).size()<<std::endl;
			std::cout<<" m_conditions.size(): "<< m_conditions.size()<<std::endl;
			std::cout<<" Preserved Conditions "<<PreservedConditions.size()<<std::endl;
			int id=0;
			for(ModelPart::ElementsContainerType::iterator ie = mr_model_part.ElementsBegin(MeshId); ie != mr_model_part.ElementsEnd(MeshId); ie++)
			{
	  
				/*each face is opposite to the corresponding node number so in 2D (triangles):
				  0 ----- 1 2
				  1 ----- 2 0
				  2 ----- 0 1
				*/ 

				//finding boundaries and creating the "skin"
				//
				//********************************************************************

				boost::numeric::ublas::matrix<unsigned int> lpofa; //points that define the faces

				Geometry< Node<3> >& rGeom = ie->GetGeometry();
				WeakPointerVector<Element >& rE = ie->GetValue(NEIGHBOUR_ELEMENTS);
				rGeom.NodesInFaces(lpofa);   
	    
				//unsigned int size=rGeom.size();

				//loop on neighbour elements of an element
				unsigned int i=0;
				for(WeakPointerVector< Element >::iterator ne = rE.begin(); ne!=rE.end(); ne++)
				{
 					if (ne->Id() == ie->Id())
					{
					        //if no neighnour is present => the face is free surface
						int sizei=lpofa.size2();
						for(int j=1; j<sizei; j++)
						{
							rGeom[lpofa(j,i)].Set(BOUNDARY);
						}
	

						//Generate condition
						Condition::NodesArrayType face;
						face.reserve(sizei);
						for(int j=1; j<sizei; j++)
						{
							face.push_back(rGeom(lpofa(j,i)));
						}
				    
						id +=1;

						Condition::GeometryType::Pointer vertices = Condition::GeometryType::Pointer(new Geometry< Node<3> >(face) );
						CompositeCondition::Pointer composite_cond =  CompositeCondition::Pointer(new CompositeCondition(id,vertices,properties) ); 
						
						//Get the correct ReferenceCondition
						Condition::Pointer pBoundaryCondition;
						bool point_condition = false;
						
						for(ModelPart::ConditionsContainerType::iterator ic = m_conditions.begin(); ic!= m_conditions.end(); ic++)
						  {
						    Geometry< Node<3> >& rConditionGeom = ic->GetGeometry();
			
						    if( FindCondition(rConditionGeom,rGeom,lpofa,i) ){

						      composite_cond->AddChild(*(ic.base()));

						      PreservedConditions[ic->Id()-1] += 1;

						      if( rConditionGeom.PointsNumber() == 1 )
							      point_condition = true;						

						    }

						  }

						if( !point_condition ){
						  //usually one MasterElement and one MasterNode in 2D; in 3D can be more than one -> it has to be extended to 3D					  
						  //composite_cond->GetValue(MASTER_ELEMENTS).push_back( Element::WeakPointer( *(ie.base()) ) );
						  WeakPointerVector< Element >& MasterElements = composite_cond->GetValue(MASTER_ELEMENTS);
						  MasterElements.push_back( Element::WeakPointer( *(ie.base()) ));
						  composite_cond->SetValue(MASTER_ELEMENTS,MasterElements);
									    
						  //composite_cond->GetValue(MASTER_NODES).push_back( Node<3>::WeakPointer( rGeom(lpofa(0,i)) ) );						
						  WeakPointerVector< Node<3> >& MasterNodes = composite_cond->GetValue(MASTER_NODES);
						  MasterNodes.push_back( Node<3>::WeakPointer( rGeom(lpofa(0,i)) ) );
						  composite_cond->SetValue(MASTER_NODES,MasterNodes);

						}


						//std::cout<<" Set preserved condition found "<<id<<std::endl;
						Condition::Pointer p_cond = composite_cond;
						mr_model_part.Conditions(MeshId).push_back(p_cond);
						

					}
					
					i+=1;
				}

			}



			//all previous conditions have to be added
			for(ModelPart::ConditionsContainerType::iterator ic = m_conditions.begin(); ic!= m_conditions.end(); ic++)
			{		    
			    if( PreservedConditions[ic->Id()-1] == 0 ){

				Geometry< Node<3> >& rGeom = ic->GetGeometry();

				if( FindConditionDomain(rGeom, MeshId) ){
								  
				  Condition::NodesArrayType face;

				  face.reserve(rGeom.size() );

				  for(unsigned int j=0; j<rGeom.size(); j++)
				    {
				      face.push_back(rGeom(j));
				    }

				  PreservedConditions[ic->Id()-1] += 1;

				  id +=1;

				  mr_model_part.Conditions(MeshId).push_back(ic->Clone(id,face));

				  //std::cout<<" Set preserved condition not found "<<ic->Id()<<std::endl;

				}
			    }
			}

			
			//control if previous conditions have been assigned
			bool all_assigned = true;
			for(unsigned int i=0; i<PreservedConditions.size(); i++)
			{
			    if( PreservedConditions[i] == 0 )
				all_assigned = false;
			}

			std::cout<<" [ NEW_CONDITIONS: "<<mr_model_part.NumberOfConditions(MeshId)<<"] "<<std::endl;

			if(all_assigned == true)
			    std::cout<<"Boundary Conditions RELOCATED "<<std::endl;
			else
			    std::cout<<"Boundary Conditions NOT relocated "<<std::endl;

	
			return true;
		};


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
		BoundarySkinBuildProcess& operator=(BoundarySkinBuildProcess const& rOther);

		/// Copy constructor.
		//BoundarySkinBuildProcess(BoundarySkinBuildProcess const& rOther);


		///@}

	}; // Class BoundarySkinBuildProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
	inline std::istream& operator >> (std::istream& rIStream,
					  BoundarySkinBuildProcess& rThis);

/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream,
					  const BoundarySkinBuildProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
///@}


}  // namespace Kratos.

#endif // KRATOS_BOUNDARY_SKIN_BUILD_PROCESS_H_INCLUDED  defined 
