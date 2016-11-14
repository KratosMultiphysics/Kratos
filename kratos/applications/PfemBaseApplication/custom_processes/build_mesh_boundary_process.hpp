//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_BUILD_MESH_BOUNDARY_PROCESS_H_INCLUDED )
#define  KRATOS_BUILD_MESH_BOUNDARY_PROCESS_H_INCLUDED


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
#include "pfem_base_application_variables.h"

///VARIABLES used:
//Data:     MASTER_ELEMENTS(set), MASTER_NODES(set), NEIGHBOUR_ELEMENTS
//StepData: RIGID_WALL
//Flags:    (checked) CONTACT
//          (set)     BOUNDARY(set)
//          (modified)  
//          (reset)   
//(set):=(set in this process)

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
  class BuildMeshBoundaryProcess
    : public Process
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BuildMeshBoundaryProcess
    KRATOS_CLASS_POINTER_DEFINITION( BuildMeshBoundaryProcess );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BuildMeshBoundaryProcess(ModelPart& rModelPart,
			     int MeshId = 0,
			     int EchoLevel = 0)
      : mrModelPart(rModelPart)
    { 
      mMeshId = MeshId;
      mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~BuildMeshBoundaryProcess()
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
      KRATOS_TRY
	
      bool success=false;

      boost::timer auxiliary;

      unsigned int NumberOfMeshes=mrModelPart.NumberOfMeshes();
	
		
      if( mMeshId == 0 ){

	unsigned int start=0;
	if(NumberOfMeshes>1) 
	  start=1;
			  
	for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
	  {
			    
	    if( mEchoLevel >= 1 )
	      std::cout<<" [ Skin Search on Mesh["<<MeshId<<"] ]"<<std::endl;

	    //success=SkinSearch(MeshId);
	    success=UniqueSkinSearch(MeshId);
			    
	    if(!success)
	      {
		std::cout<<"  ERROR:  Skin Search FAILED on mesh : ["<<MeshId<<"] "<<std::endl;
	      }
	    else
	      {
		if( mEchoLevel >= 1 )
		  std::cout<<" [ Search performed in Time = "<<auxiliary.elapsed()<<" ]"<<std::endl;
		//PrintSkin(MeshId);
	      }
	  }
      }
      else{
	
	if( mEchoLevel >= 1 )
	  std::cout<<" [ Skin Search on Mesh["<<mMeshId<<"] ]"<<std::endl;

	//success=SkinSearch(MeshId);
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
      }
	
      if( NumberOfMeshes > 1 )
	SetGlobalConditions();

      
      if( mEchoLevel >= 1 )
	std::cout<<"  ::[SET BOUNDARY NORMALS]:: "<<std::endl;

      //ComputeBoundaryNormals BoundUtils;
      BoundaryNormalsCalculationUtilities BoundaryComputation;
      BoundaryComputation.CalculateWeightedBoundaryNormals(mrModelPart, mEchoLevel);

	
      if( mEchoLevel >= 1 )
	std::cout<<"   -> Boundary Normals Computed <- "<<std::endl;

     KRATOS_CATCH(" ")
      
    }



   bool SearchConditionMasters(int MeshId = 0)
    {
      
      KRATOS_TRY
	
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

      KRATOS_CATCH(" ")    

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
      return "BuildMeshBoundaryProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
      rOStream << "BuildMeshBoundaryProcess";
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

    ModelPart& mrModelPart;

    int mMeshId;
    
    int mEchoLevel;


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    //**************************************************************************
    //**************************************************************************


    virtual bool UniqueSkinSearch( int MeshId = 0 )
    {
      KRATOS_TRY
      
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
	 
	    WeakPointerVector<Element >& rE = ie->GetValue(NEIGHBOUR_ELEMENTS);
	    
	    //get matrix nodes in faces
	    rGeometry.NodesInFaces(lpofa);
	    rGeometry.NumberNodesInFaces(lnofa);
	    
	    //loop on neighbour elements of an element
	    unsigned int iface=0;
	    for(WeakPointerVector< Element >::iterator ne = rE.begin(); ne!=rE.end(); ne++)
	      {
		unsigned int NumberNodesInFace = lnofa[iface];

		if (ne->Id() == ie->Id())
		  {
		    //if no neighbour is present => the face is free surface
		    for(unsigned int j=1; j<=NumberNodesInFace; j++)
		      {
			rGeometry[lpofa(j,iface)].Set(BOUNDARY);
		      }

		    //1.- create geometry: points array and geometry type
		    Condition::NodesArrayType        FaceNodes;
		    Condition::GeometryType::Pointer ConditionVertices;
		      
		    FaceNodes.reserve(NumberNodesInFace);

		    for(unsigned int j=1; j<=NumberNodesInFace; j++)
		      {
			FaceNodes.push_back(rGeometry(lpofa(j,iface)));
		      }
				    
							
		    if( NumberNodesInFace == 2 ){					  
		      ConditionVertices = Condition::GeometryType::Pointer(new Line2D2< Node<3> >(FaceNodes) );
		    }
		    else if ( NumberNodesInFace == 3 ){
		      ConditionVertices = Condition::GeometryType::Pointer(new Triangle3D3< Node<3> >(FaceNodes) );
		    }

		    ConditionId +=1;
		    
		    //Create a composite condition
		    CompositeCondition::Pointer p_cond = CompositeCondition::Pointer(new CompositeCondition(ConditionId,ConditionVertices,properties) ); 

		    bool condition_found = false;
		    bool point_condition = false;
					       
		    // Search for existing conditions: start
		    for(ModelPart::ConditionsContainerType::iterator ic = TemporaryConditions.begin(); ic!= TemporaryConditions.end(); ic++)
		      {
			Geometry< Node<3> >& rConditionGeometry = ic->GetGeometry();
						
			condition_found = this->FindCondition(rConditionGeometry,rGeometry,lpofa,lnofa,iface);

			if( condition_found ){

			  p_cond->AddChild(*(ic.base()));

			  PreservedConditions[ic->Id()-1] += 1;		  

			  if( rConditionGeometry.PointsNumber() == 1 )
			    point_condition = true;
			}
			
		      }
		    // Search for existing conditions: end

		    if( !point_condition ){
		      // usually one MasterElement and one MasterNode in 2D; in 3D can be more than one -> it has to be extended to other 3D geometries
		      //p_cond->GetValue(MASTER_ELEMENTS).push_back( Element::WeakPointer( *(ie.base()) ) );
		      WeakPointerVector< Element >& MasterElements = p_cond->GetValue(MASTER_ELEMENTS);
		      MasterElements.push_back( Element::WeakPointer( *(ie.base()) ) );
		      p_cond->SetValue(MASTER_ELEMENTS,MasterElements);

		      //p_cond->GetValue(MASTER_NODES).push_back( Node<3>::WeakPointer( rGeometry(lpofa(0,i)) ) );	
		      WeakPointerVector< Node<3> >& MasterNodes = p_cond->GetValue(MASTER_NODES);
		      MasterNodes.push_back( Node<3>::WeakPointer( rGeometry(lpofa(0,iface)) ) );
		      p_cond->SetValue(MASTER_NODES,MasterNodes);
		    }

		    mrModelPart.AddCondition(Condition::Pointer(p_cond), MeshId);
		    //mrModelPart.Conditions(MeshId).push_back(Condition::Pointer(p_cond));

		    // Set new conditions: end
		    
		  }
					
		iface+=1;
	      }

	  }
	}


      this->AddOtherConditions(TemporaryConditions, PreservedConditions, ConditionId, MeshId);
	
      return true;

      KRATOS_CATCH(" ")
     
    }

    //**************************************************************************
    //**************************************************************************

    bool FindCondition(Geometry< Node<3> >& rConditionGeometry,Geometry< Node<3> >& rGeometry , boost::numeric::ublas::matrix<unsigned int>& lpofa, boost::numeric::ublas::vector<unsigned int>& lnofa, unsigned int& iface)
    {
      KRATOS_TRY
      
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

     KRATOS_CATCH(" ")
      
    }


    //**************************************************************************
    //**************************************************************************


    bool FindConditionID(Geometry< Node<3> >& rConditionGeometry, int& MeshId)
    {
      KRATOS_TRY
      
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

      KRATOS_CATCH(" ")	   
    }


    //**************************************************************************
    //**************************************************************************


    void PrintSkin (ModelPart::IndexType MeshId=0)
    {

      KRATOS_TRY
	
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

      KRATOS_CATCH(" ")
    }

    //**************************************************************************
    //**************************************************************************


    virtual bool AddOtherConditions(ModelPart::ConditionsContainerType& rTemporaryConditions, std::vector<int>& PreservedConditions, unsigned int& rConditionId, int MeshId = 0 )
    {

      KRATOS_TRY
	
      //add all previous conditions not found in the skin search are added:
      for(ModelPart::ConditionsContainerType::iterator ic = rTemporaryConditions.begin(); ic!= rTemporaryConditions.end(); ic++)
	{	

	  if( PreservedConditions[ic->Id()-1] == 0 ){

	    Geometry< Node<3> >& rGeometry = ic->GetGeometry();
				
	    if( FindConditionID(rGeometry, MeshId) ){

	      Condition::NodesArrayType FaceNodes;

	      FaceNodes.reserve(rGeometry.size() );

	      for(unsigned int j=0; j<rGeometry.size(); j++)
		{
		  FaceNodes.push_back(rGeometry(j));
		}

	      PreservedConditions[ic->Id()-1] += 1;

	      rConditionId +=1;


	      Condition::Pointer p_cond = ic->Clone(rConditionId, FaceNodes);
	      p_cond->Data() = ic->Data();
	      
	      mrModelPart.AddCondition(p_cond,MeshId);
	      //mrModelPart.Conditions(MeshId).push_back(ic->Clone(rConditionId,FaceNodes));

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

	std::cout<<"    New Conditions: "<<mrModelPart.NumberOfConditions(MeshId)<<"] [MESH:"<<MeshId<<"]"<<std::endl;

	if(all_assigned == true)
	  std::cout<<"   Boundary Conditions RELOCATED "<<std::endl;
	else
	  std::cout<<"   Boundary Conditions NOT relocated "<<std::endl;
      }

      return all_assigned;
      
      KRATOS_CATCH(" ")

    }

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

    //**************************************************************************
    //**************************************************************************


    bool SkinSearch( int MeshId = 0 )
    {

      KRATOS_TRY
      
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
	 
	    WeakPointerVector<Element >& rE = ie->GetValue(NEIGHBOUR_ELEMENTS);
	    
	    //get matrix nodes in faces
	    rGeometry.NodesInFaces(lpofa);
	    rGeometry.NumberNodesInFaces(lnofa);
	    
	    //loop on neighbour elements of an element
	    unsigned int iface=0;
	    for(WeakPointerVector< Element >::iterator ne = rE.begin(); ne!=rE.end(); ne++)
	      {
		unsigned int NumberNodesInFace = lnofa[iface];

		if (ne->Id() == ie->Id())
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
						
		    // Search for existing conditions: start
		    for(ModelPart::ConditionsContainerType::iterator ic = TemporaryConditions.begin(); ic!= TemporaryConditions.end(); ic++)
		      {
			Geometry< Node<3> >& rConditionGeometry = ic->GetGeometry();
						
			condition_found = this->FindCondition(rConditionGeometry,rGeometry,lpofa,lnofa,iface);

			if( condition_found ){

			  pBoundaryCondition = (*(ic.base())); //accessing boost::shared_ptr  get() to obtain the raw pointer

			  PreservedConditions[ic->Id()-1] += 1;		  

			  if( rConditionGeometry.PointsNumber() == 1 )
			    point_condition = true;
						
			  break;
			}
			
		      }
		    // Search for existing conditions: end


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
				   						
		      if( NumberNodesInFace == 2 ){					  
			ConditionVertices = Condition::GeometryType::Pointer(new Line2D2< Node<3> >(FaceNodes) );
		      }
		      else if ( NumberNodesInFace == 3 ){
			ConditionVertices = Condition::GeometryType::Pointer(new Triangle3D3< Node<3> >(FaceNodes) );
		      }

		      ConditionId +=1;

		      //Create a condition
		      Condition::Pointer p_cond;
		      if(condition_found)
			{
			  p_cond = pBoundaryCondition->Clone(ConditionId, FaceNodes);
			}
		      else
			{
			  p_cond = Condition::Pointer(new Condition(ConditionId,ConditionVertices,properties) ); 
			  //p_cond = mr_reference_condition.Create(ConditionId, ConditionVertices, properties);
			}
 
		      //usually one MasterElement and one MasterNode in 2D in 3D can be more than one

		      //p_cond->GetValue(MASTER_ELEMENTS).push_back( Element::WeakPointer( *(ie.base()) ) );
		      WeakPointerVector< Element >& MasterElements = p_cond->GetValue(MASTER_ELEMENTS);
		      MasterElements.push_back( Element::WeakPointer( *(ie.base()) ) );
		      p_cond->SetValue(MASTER_ELEMENTS,MasterElements);

		      //p_cond->GetValue(MASTER_NODES).push_back( Node<3>::WeakPointer( rGeometry(lpofa(0,i)) ) );			
		      WeakPointerVector< Node<3> >& MasterNodes = p_cond->GetValue(MASTER_NODES);
		      MasterNodes.push_back( Node<3>::WeakPointer( rGeometry(lpofa(0,iface)) ) );
		      p_cond->SetValue(MASTER_NODES,MasterNodes);

		      mrModelPart.AddCondition(p_cond, MeshId);
		      //mrModelPart.Conditions(MeshId).push_back(p_cond);

		    }
		    // Set new conditions: end

		  }
					
		iface+=1;
	      }

	  }
	}


      this->AddOtherConditions(TemporaryConditions, PreservedConditions, ConditionId, MeshId);
	
      return true;

      KRATOS_CATCH(" ")
      
    }


    //**************************************************************************
    //**************************************************************************

    void SetGlobalConditions()
    {

      KRATOS_TRY

      if( mEchoLevel >= 1 ){
	std::cout<<" [MESH:0]: "<<std::endl;
	std::cout<<" [OLD TOTAL CONDITIONS: "<<mrModelPart.NumberOfConditions()<<"] "<<std::endl;
      }

      //contact conditions are located on Mesh_0
      ModelPart::ConditionsContainerType KeepConditions;

      unsigned int condId=1;
      unsigned int start=0;  
      unsigned int NumberOfMeshes=mrModelPart.NumberOfMeshes();
      if(NumberOfMeshes>1) 
	start=1;

      for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
	{
	  for(ModelPart::ConditionsContainerType::iterator i_cond = mrModelPart.ConditionsBegin(MeshId) ; i_cond != mrModelPart.ConditionsEnd(MeshId) ; i_cond++)
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


      for(ModelPart::ConditionsContainerType::iterator i_cond = mrModelPart.ConditionsBegin(); i_cond!= mrModelPart.ConditionsEnd(); i_cond++)
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
      
      mrModelPart.Conditions().swap(KeepConditions);

      if( mEchoLevel >= 1 )
	std::cout<<" [NEW TOTAL CONDITIONS: "<<mrModelPart.NumberOfConditions()<<"] "<<std::endl;

      KRATOS_CATCH(" ")     
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
    BuildMeshBoundaryProcess& operator=(BuildMeshBoundaryProcess const& rOther);

    /// Copy constructor.
    //BuildMeshBoundaryProcess(BuildMeshBoundaryProcess const& rOther);


    ///@}

  }; // Class BuildMeshBoundaryProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    BuildMeshBoundaryProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const BuildMeshBoundaryProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // KRATOS_BUILD_MESH_BOUNDARY_PROCESS_H_INCLUDED  defined 
