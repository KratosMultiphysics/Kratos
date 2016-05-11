//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_RECONSTRUCT_MESH_BOUNDARY_PROCESS_H_INCLUDED )
#define  KRATOS_RECONSTRUCT_MESH_BOUNDARY_PROCESS_H_INCLUDED


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

#include "geometries/line_2d_2.h"
#include "geometries/triangle_3d_3.h"

#include "custom_conditions/composite_condition.hpp"
#include "custom_utilities/boundary_normals_calculation_utilities.hpp"
#include "custom_utilities/modeler_utilities.hpp"

#include "pfem_base_application_variables.h"

///VARIABLES used:
//Data:     MASTER_ELEMENTS(set), MASTER_NODES(set)
//StepData: 
//Flags:    (checked) TO_ERASE, TO_REFINE, CONTACT, NEW_ENTITY
//          (set)     BOUNDARY(set),  [TO_REFINE(nodes), TO_ERASE(condition)]->locally to not preserve condition
//          (modified)  
//          (reset)   
// (set):=(set in this process)

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
  class ReconstructMeshBoundaryProcess
    : public Process
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ReconstructMeshBoundaryProcess
    KRATOS_CLASS_POINTER_DEFINITION( ReconstructMeshBoundaryProcess );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ReconstructMeshBoundaryProcess(ModelPart& rModelPart,
				   ModelerUtilities::MeshingParameters& rRemeshingParameters,
				   ModelPart::IndexType MeshId,
				   int EchoLevel)
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    { 
      mMeshId = MeshId;
      mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~ReconstructMeshBoundaryProcess()
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
	
      if( mEchoLevel >= 1 )
	std::cout<<" [ Skin Search on Mesh["<<mMeshId<<"] ]"<<std::endl;

      success=UniqueSkinSearch(mMeshId);
			    
      if(!success)
	{
          std::cout<<"  ERROR:  Skin Search FAILED on mesh : ["<<mMeshId<<"] "<<std::endl;
        }
      else
	{
          if( mEchoLevel >= 1 )
	    std::cout<<" [ Search performed in Time = "<<auxiliary.elapsed()<<" ]"<<std::endl;
            //PrintSkin(mMeshId);
        }
  	
    };


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
      return "ReconstructMeshBoundaryProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
      rOStream << "ReconstructMeshBoundaryProcess";
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
    ModelPart& mrModelPart;

    ModelerUtilities::MeshingParameters& mrRemesh;

    int mMeshId;
    int mEchoLevel;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    bool SearchConditionMasters(int MeshId = 0)
    {

      unsigned int counter = 0;
      bool found=false;

      for(ModelPart::ConditionsContainerType::iterator ic = mrModelPart.ConditionsBegin(MeshId); ic != mrModelPart.ConditionsEnd(MeshId); ic++)
	{
	  
	  //std::cout<<" Condition ("<<ic->Id()<<") : ME="<<ic->GetValue(MASTER_ELEMENTS)[0].Id()<<", MN= "<<ic->GetValue(MASTER_NODES)[0].Id()<<std::endl;

	  //********************************************************************

	  boost::numeric::ublas::matrix<unsigned int> lpofa; //points that define the faces

	  Geometry< Node<3> >& rConditionGeometry = ic->GetGeometry();
	  unsigned int size=rConditionGeometry.size();
			    
	  bool perform_search = true;
	  for(unsigned int i=0; i<size; i++)
	    {
	      if( rConditionGeometry[i].SolutionStepsDataHas(RIGID_WALL) ){
		if( rConditionGeometry[i].FastGetSolutionStepValue(RIGID_WALL) ) //if is a rigid wall do not search else do search
		  perform_search = false;
	      }
	    }		   		     

	  if( size != 2 ) 
	    perform_search = false;

	  //********************************************************************
	  found=false;

	  if( perform_search )
	    {

	      WeakPointerVector<Element >& rE1 = rConditionGeometry[0].GetValue(NEIGHBOUR_ELEMENTS);				    
	      WeakPointerVector<Element >& rE2 = rConditionGeometry[1].GetValue(NEIGHBOUR_ELEMENTS);

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
			      if( (   rConditionGeometry[0].Id() == rElementGeom[lpofa(1,i)].Id() 
				      && rConditionGeometry[1].Id() == rElementGeom[lpofa(2,i)].Id() ) || 
				  (   rConditionGeometry[0].Id() == rElementGeom[lpofa(2,i)].Id() 
				      && rConditionGeometry[1].Id() == rElementGeom[lpofa(1,i)].Id() ) )
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
      if(mrModelPart.Conditions(MeshId).size()>0)
	totalcond = mrModelPart.Conditions(MeshId).size();
			  

      if(counter == totalcond){
	if( mEchoLevel > 1 )
	  std::cout<<"   Condition Masters (mesh "<<MeshId<<"): LOCATED ["<<counter<<"]"<<std::endl;
	found=true;
      }
      else{
	if( mEchoLevel > 1 )
	  std::cout<<"   Condition Masters (mesh "<<MeshId<<"): not LOCATED ["<<counter-totalcond<<"]"<<std::endl;
	found=false;
      }
			
      return found;

    }


    void PrintSkin (ModelPart::IndexType MeshId=0)
    {
      //PRINT SKIN:		
      std::cout<<" CONDITIONS: geometry nodes ("<<mrModelPart.Conditions(MeshId).size()<<")"<<std::endl;

      ConditionsContainerType& rCond = mrModelPart.Conditions(MeshId);
      for(ConditionsContainerType::iterator ic = rCond.begin(); ic!= rCond.end(); ic++)
	{
			
	  Geometry< Node<3> >& rConditionGeometry = ic->GetGeometry();
	  std::cout<<"["<<ic->Id()<<"]:"<<std::endl;
	  //ic->PrintInfo(std::cout);
	  std::cout<<"( ";
	  for(unsigned int i = 0; i < rConditionGeometry.size(); i++)
	    {
	      std::cout<< rConditionGeometry[i].Id()<<", ";
	    }
	  std::cout<<" ): ";

	  ic->GetValue(MASTER_ELEMENTS)[0].PrintInfo(std::cout);
							
	  std::cout<<std::endl;

	}
      std::cout<<std::endl;

    }


    bool FindNodeInCondition(Geometry< Node<3> >& rConditionGeometry,Geometry< Node<3> >& rGeometry , boost::numeric::ublas::matrix<unsigned int>& lpofa, boost::numeric::ublas::vector<unsigned int>& lnofa, unsigned int& iface)
    {
      
      // not equivalent geometry sizes for boundary conditions:
      if( rConditionGeometry.size() != lnofa[iface] )
	return false;
      
      // line boundary condition:
      if( lnofa[iface] == 2 )
	{
	  if( rConditionGeometry[0].Id() == rGeometry[lpofa(1,iface)].Id()  ||
	      rConditionGeometry[1].Id() == rGeometry[lpofa(2,iface)].Id()  || 
	      rConditionGeometry[0].Id() == rGeometry[lpofa(2,iface)].Id()  ||
	      rConditionGeometry[1].Id() == rGeometry[lpofa(1,iface)].Id()  )
	  {	 
	    return true;
	  }
	  else
	  {
	    return false;
	  }
	    
	}
      
      //3D faces:
      if(  lnofa[iface] == 3 )
	{
	  if( rConditionGeometry[0].Id() == rGeometry[lpofa(1,iface)].Id() || 
	      rConditionGeometry[1].Id() == rGeometry[lpofa(2,iface)].Id() ||
	      rConditionGeometry[2].Id() == rGeometry[lpofa(3,iface)].Id() || 
	      rConditionGeometry[0].Id() == rGeometry[lpofa(3,iface)].Id() || 
	      rConditionGeometry[1].Id() == rGeometry[lpofa(1,iface)].Id() ||
	      rConditionGeometry[2].Id() == rGeometry[lpofa(2,iface)].Id() ||
	      rConditionGeometry[0].Id() == rGeometry[lpofa(2,iface)].Id() ||
	      rConditionGeometry[1].Id() == rGeometry[lpofa(3,iface)].Id() ||
	      rConditionGeometry[2].Id() == rGeometry[lpofa(1,iface)].Id()  )
	  {
	    return true;
	  }
	  else
	  {
	    return false;
	  }
	  
	}

      if(  lnofa[iface] > 3 )
	{
	  KRATOS_THROW_ERROR( std::logic_error, "Wrong Condition Number of Face Nodes",*this );
	}

      return false;
  
    }


    bool FindCondition(Geometry< Node<3> >& rConditionGeometry,Geometry< Node<3> >& rGeometry , boost::numeric::ublas::matrix<unsigned int>& lpofa, boost::numeric::ublas::vector<unsigned int>& lnofa, unsigned int& iface)
    {
      
      // not equivalent geometry sizes for boundary conditions:
      if( rConditionGeometry.size() != lnofa[iface] )
	return false;
      
      // line boundary condition:
      if( lnofa[iface] == 2 )
	{
	  if( (   rConditionGeometry[0].Id() == rGeometry[lpofa(1,iface)].Id() 
		  && rConditionGeometry[1].Id() == rGeometry[lpofa(2,iface)].Id() ) || 
	      (   rConditionGeometry[0].Id() == rGeometry[lpofa(2,iface)].Id() 
		  && rConditionGeometry[1].Id() == rGeometry[lpofa(1,iface)].Id() ) )
	  {	 
	    return true;
	  }
	  else
	  {
	    return false;
	  }
	    
	}
      
      //3D faces:
      if(  lnofa[iface] == 3 )
	{
	  if( (   rConditionGeometry[0].Id() == rGeometry[lpofa(1,iface)].Id() 
		  && rConditionGeometry[1].Id() == rGeometry[lpofa(2,iface)].Id()
		  && rConditionGeometry[2].Id() == rGeometry[lpofa(3,iface)].Id() ) || 
	      (   rConditionGeometry[0].Id() == rGeometry[lpofa(3,iface)].Id() 
		  && rConditionGeometry[1].Id() == rGeometry[lpofa(1,iface)].Id()
		  && rConditionGeometry[2].Id() == rGeometry[lpofa(2,iface)].Id() ) ||
	      (   rConditionGeometry[0].Id() == rGeometry[lpofa(2,iface)].Id() 
		  && rConditionGeometry[1].Id() == rGeometry[lpofa(3,iface)].Id()
		  && rConditionGeometry[2].Id() == rGeometry[lpofa(1,iface)].Id() ) )
	  {
	    return true;
	  }
	  else
	  {
	    return false;
	  }
	  
	}

      if(  lnofa[iface] > 3 )
	{
	  KRATOS_THROW_ERROR( std::logic_error, "Wrong Condition Number of Face Nodes",*this );
	}

      return false;
  
    }


    bool FindConditionID(Geometry< Node<3> >& rConditionGeometry, int& MeshId)
    {

      //check if the conditions belongs to the MeshId checking the nodes Id
      for(unsigned int i=0; i<rConditionGeometry.size(); i++)
	{
	  for(ModelPart::NodesContainerType::const_iterator in = mrModelPart.NodesBegin(MeshId); in!=mrModelPart.NodesEnd(MeshId); in++)
	    {			
	      if( rConditionGeometry[i].Id() == in->Id() )
		return true;
	    }
	}

      return false;
    }


    bool UniqueSkinSearch( int MeshId = 0 )
    {

      if( mEchoLevel > 0 ){
	std::cout<<" [ SET BOUNDARY CONDITIONS : "<<std::endl;
	std::cout<<"   Initial Conditions : "<<mrModelPart.Conditions(MeshId).size()<<" [MESH:"<<MeshId<<"]"<<std::endl;
      }

      //properties to be used in the generation
      int number_properties = mrModelPart.NumberOfProperties();
      Properties::Pointer properties = mrModelPart.GetMesh().pGetProperties(number_properties-1);
			
      //reset the boundary flag in all nodes
      for(ModelPart::NodesContainerType::const_iterator in = mrModelPart.NodesBegin(MeshId); in!=mrModelPart.NodesEnd(MeshId); in++)
	{
	  in->Reset(BOUNDARY);
	}


      //swap conditions for a temporary use
      ModelPart::ConditionsContainerType TemporaryConditions;
      TemporaryConditions.reserve(mrModelPart.Conditions(MeshId).size());
      TemporaryConditions.swap(mrModelPart.Conditions(MeshId));

      //set consecutive ids in the mesh conditions
      unsigned int ConditionId=1;
      for(ModelPart::ConditionsContainerType::iterator ic = TemporaryConditions.begin(); ic!= TemporaryConditions.end(); ic++)
	{
	  // remeshing rebuild
	  Geometry< Node<3> >& rConditionGeometry = ic->GetGeometry();
	  for( unsigned int i=0; i<rConditionGeometry.size(); i++ )
	    {
	      if( rConditionGeometry[i].Is(TO_ERASE)){
		ic->Set(TO_ERASE);
		break;
	      }
	    }
	  // remeshing rebuild

	  ic->SetId(ConditionId);
	  ConditionId++;
	}

     
      //control the previous mesh conditions
      std::vector<int> PreservedConditions( TemporaryConditions.size() );
      std::fill( PreservedConditions.begin(), PreservedConditions.end(), 0 );
		

      ModelPart::ElementsContainerType::iterator elements_begin  = mrModelPart.ElementsBegin(MeshId);
      ModelPart::ElementsContainerType::iterator elements_end    = mrModelPart.ElementsEnd(MeshId);
    
      ConditionId=0;
      for(ModelPart::ElementsContainerType::iterator ie = elements_begin; ie != elements_end ; ie++)
	{
	  
	  Geometry< Node<3> >& rGeometry = ie->GetGeometry();
	  
	  if( rGeometry.FacesNumber() >= 3 ){ //3 or 4

	    /*each face is opposite to the corresponding node number so in 2D triangle
	      0 ----- 1 2
	      1 ----- 2 0
	      2 ----- 0 1
	    */

	    /*each face is opposite to the corresponding node number so in 3D tetrahedron
	      0 ----- 1 2 3
	      1 ----- 2 0 3
	      2 ----- 0 1 3
	      3 ----- 0 2 1
	    */

	    //finding boundaries and creating the "skin"
	    //
	    //********************************************************************

	    boost::numeric::ublas::matrix<unsigned int> lpofa; //connectivities of points defining faces
	    boost::numeric::ublas::vector<unsigned int> lnofa; //number of points defining faces
	 
	    // WeakPointerVector<Element >& rE = ie->GetValue(NEIGHBOUR_ELEMENTS);
	    
	    //get matrix nodes in faces
	    rGeometry.NodesInFaces(lpofa);
	    rGeometry.NumberNodesInFaces(lnofa);
	    

	    int Id=ie->Id() -1 ;

	    ModelPart::ElementsContainerType::iterator element_neighbour;

	    //Get the standard ReferenceCondition
	    const Condition & rReferenceCondition = mrRemesh.GetReferenceCondition();
	    
	    //loop on element faces
	    for(unsigned int iface = 0; iface<rGeometry.FacesNumber(); iface++)
	      {

		unsigned int NumberNodesInFace = lnofa[iface];

		int index = mrRemesh.NeighbourList[Id][iface];
		
		if(index > 0)
		  {
		    index = mrRemesh.PreservedElements[index-1];
		    
		    //check if this element exists:: if not found-> returns the last element
		    element_neighbour = (mrModelPart.Elements(MeshId)).find( elements_begin->Id() + index-1 ); 
		  }
		else
		  {
		    element_neighbour = elements_end;
		  }

	      
		if( element_neighbour == elements_end )
		  {
		    
		    //if no neighbour is present => the face is free surface
		    for(unsigned int j=1; j<=NumberNodesInFace; j++)
		      {
			rGeometry[lpofa(j,iface)].Set(BOUNDARY);
		      }
	
	
		    //Get the correct ReferenceCondition
		    Condition::Pointer pBoundaryCondition;
		    bool condition_found = false;
		    bool point_condition = false; 

		    bool inserted = false;
		    for(ModelPart::ConditionsContainerType::iterator ic = TemporaryConditions.begin(); ic!= TemporaryConditions.end(); ic++)
		      {
			Geometry< Node<3> >& rConditionGeometry = ic->GetGeometry();

			if( ic->IsNot(TO_ERASE) ){

			  if( ic->IsNot(CONTACT) ){
			  
			    if(ic->Is(NEW_ENTITY)){
			      inserted = false;
			    }
			    else{
			      // remeshing rebuild
			      for( unsigned int i=0; i<rConditionGeometry.size(); i++ )
				{
				  if( rConditionGeometry[i].Is(TO_ERASE)){
				    inserted = true;
				    break;
				  }
				}
			      // remeshing rebuild
			    }

			    if( !inserted ){
			      
			      if( PreservedConditions[ic->Id()-1] == 0 ){
			      
				condition_found = this->FindCondition(rConditionGeometry,rGeometry,lpofa,lnofa,iface);
				
				if( condition_found ){
				
				  pBoundaryCondition = (*(ic.base())); //accessing boost::shared_ptr  get() to obtain the raw pointer
				  PreservedConditions[ic->Id()-1] += 1; //add each time is used

				  if( rConditionGeometry.PointsNumber() == 1 )
				    point_condition = true;

				  //break;
				}
			      }
			    
			    }
			    else{
			    
			      if( PreservedConditions[ic->Id()-1] < 2 ){
				
				condition_found = this->FindNodeInCondition(rConditionGeometry,rGeometry,lpofa,lnofa,iface);
				
				if( condition_found ){
				
				  pBoundaryCondition = (*(ic.base())); //accessing boost::shared_ptr  get() to obtain the raw pointer
				  PreservedConditions[ic->Id()-1] += 1; //add each time is used

				  if( rConditionGeometry.PointsNumber() == 1 )
				    point_condition = true;
				  
				  //break;
				}	

				// std::cout<<" INSERTED COND "<<ic->Id()<<std::endl;
			      }			
			  
			    }
			    			    			    			    
			  }
			  else{
			    
			    PreservedConditions[ic->Id()-1] += 1;  //will not be restored
			    //std::cout<<" Condition Contact "<<ic->Id()<<std::endl;

			  }
			
			}
			  

			if(condition_found==true){
			  // std::cout<<" Condition Found:  "<<ic->Id()<<" ("<<ic->GetGeometry()[0].Id()<<", "<<ic->GetGeometry()[1].Id()<<") == ("<<rGeom[lpofa(1,i)].Id()<<" "<<rGeom[lpofa(2,i)].Id()<<") ->  Used : "<<PreservedConditions[ic->Id()-1]<<" times "<<std::endl;
			  break;
			}
		      }


		    // Set new conditions:  start
		    if( !point_condition ){

		      //1.- create geometry: points array and geometry type
		      
		      Condition::NodesArrayType        FaceNodes;
		      Condition::GeometryType::Pointer ConditionVertices;
		      
		      FaceNodes.reserve(NumberNodesInFace);

		      for(unsigned int j=1; j<=NumberNodesInFace; j++)
			{
			  FaceNodes.push_back(rGeometry(lpofa(j,iface)));
			}

		      ConditionId +=1;

		      //Create a condition
		      Condition::Pointer p_cond;
		      if(condition_found)
			{
			  p_cond = pBoundaryCondition->Clone(ConditionId, FaceNodes);
		      
			  p_cond->Data() =pBoundaryCondition->Data();
		      
			  WeakPointerVector< Element >& MasterElements = p_cond->GetValue(MASTER_ELEMENTS);
			  MasterElements.push_back( Element::WeakPointer( *(ie.base()) ) );
			  p_cond->SetValue(MASTER_ELEMENTS,MasterElements);

			  //p_cond->GetValue(MASTER_NODES).push_back( Node<3>::WeakPointer( rGeometry(lpofa(0,i)) ) );			
			  WeakPointerVector< Node<3> >& MasterNodes = p_cond->GetValue(MASTER_NODES);
			  MasterNodes.push_back( Node<3>::WeakPointer( rGeometry(lpofa(NumberNodesInFace,iface)) ) );
			  p_cond->SetValue(MASTER_NODES,MasterNodes);

			}
		      else
			{
		  
			  if( mEchoLevel > 0 )
			    std::cout<<"   NOT FOUND CONDITION :: CREATED-> ["<<ConditionId<<"] ("<<FaceNodes[0].Id()<<","<<FaceNodes[1].Id()<<")"<<std::endl;
			  
			  //std::cout<<" ReferenceCondition "<<rReferenceCondition<<std::endl;

			  p_cond = rReferenceCondition.Create(ConditionId, FaceNodes, properties);
		      
			  //if a condition is created new nodes must be labeled TO_REFINE
			  for(unsigned int j=0; j<NumberNodesInFace; j++)
			    {
			      FaceNodes[j].Set(TO_REFINE);
			    }

			  MeshDataTransferUtilities TransferUtilities;

			  TransferUtilities.InitializeBoundaryData(p_cond, mrRemesh.TransferVariables); 

			  WeakPointerVector< Element >& MasterElements = p_cond->GetValue(MASTER_ELEMENTS);
			  MasterElements.push_back( Element::WeakPointer( *(ie.base()) ) );
			  p_cond->SetValue(MASTER_ELEMENTS,MasterElements);

			  //p_cond->GetValue(MASTER_NODES).push_back( Node<3>::WeakPointer( rGeometry(lpofa(0,i)) ) );			
			  WeakPointerVector< Node<3> >& MasterNodes = p_cond->GetValue(MASTER_NODES);
			  MasterNodes.push_back( Node<3>::WeakPointer( rGeometry(lpofa(NumberNodesInFace,iface)) ) );
			  p_cond->SetValue(MASTER_NODES,MasterNodes);

			}


		      mrModelPart.AddCondition(p_cond, MeshId);
		    }
		    // Set new conditions: end
		  }
	      }

	  }
	}

      this->AddOtherConditions(TemporaryConditions, PreservedConditions, ConditionId, MeshId);
	
      return true;
    }

    bool AddOtherConditions(ModelPart::ConditionsContainerType& rTemporaryConditions, std::vector<int>& PreservedConditions, unsigned int& rConditionId, int MeshId = 0 )
    {
      //add all previous conditions not found in the skin search are added:
      for(ModelPart::ConditionsContainerType::iterator ic = rTemporaryConditions.begin(); ic!= rTemporaryConditions.end(); ic++)
	{		    

	  bool node_not_preserved = false;
	  bool condition_not_preserved = false;

	  if( PreservedConditions[ic->Id()-1] == 0 ){

	    Geometry< Node<3> >& rGeometry = ic->GetGeometry();
	    
	    Condition::NodesArrayType FaceNodes;

	    FaceNodes.reserve(rGeometry.size() );

	    for(unsigned int j=0; j<rGeometry.size(); j++)
	      {
		FaceNodes.push_back(rGeometry(j));
		if( FaceNodes[j].Is(TO_ERASE) || FaceNodes[j].Is(TO_REFINE) )
		  node_not_preserved = true;
	      }

	    if( ic->Is(TO_ERASE) )
	      condition_not_preserved = true;

	    if(node_not_preserved == true || condition_not_preserved == true)
	      continue;

	    PreservedConditions[ic->Id()-1] += 1;

	    rConditionId +=1;

	    Condition::Pointer p_cond = ic->Clone(rConditionId, FaceNodes);
	    p_cond->Data() = ic->Data();

	    mrModelPart.AddCondition(p_cond,MeshId);
	    //mrModelPart.Conditions(MeshId).push_back(ic->Clone(rConditionId,FaceNodes));

	    if( mEchoLevel > 0 ){
	      std::cout<<" Temporal Condition Not Set "<<ic->Id()<<"("<<ic->GetGeometry()[0].Id()<<","<<ic->GetGeometry()[1].Id()<<")"<<std::endl;
	      std::cout<<" Push Back Not Set Conditions "<<rConditionId<<"("<<FaceNodes[0].Id()<<","<<FaceNodes[1].Id()<<")"<<std::endl;
	    }

	  }
	}


      //control if all previous conditions have been added:
      bool all_assigned = true;
      for(unsigned int i=0; i<PreservedConditions.size(); i++)
	{
	  if( PreservedConditions[i] == 0 )
	    all_assigned = false;
	}


      if( mEchoLevel >= 1 ){

	std::cout<<"   New Conditions : "<<mrModelPart.NumberOfConditions(MeshId)<<"] [MESH:"<<MeshId<<"]"<<std::endl;

	if(all_assigned == true)
	  std::cout<<"   Boundary Conditions RELOCATED "<<std::endl;
	else
	  std::cout<<"   Boundary Conditions NOT relocated "<<std::endl;
      }

      return all_assigned;

    }


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
    ReconstructMeshBoundaryProcess& operator=(ReconstructMeshBoundaryProcess const& rOther);

    /// Copy constructor.
    //ReconstructMeshBoundaryProcess(ReconstructMeshBoundaryProcess const& rOther);


    ///@}

  }; // Class ReconstructMeshBoundaryProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    ReconstructMeshBoundaryProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const ReconstructMeshBoundaryProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // KRATOS_RECONSTRUCT_MESH_BOUNDARY_PROCESS_H_INCLUDED  defined 
