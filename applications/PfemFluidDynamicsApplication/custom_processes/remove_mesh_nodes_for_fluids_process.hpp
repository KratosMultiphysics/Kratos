//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:                   AFranci $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                October 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_REMOVE_MESH_NODES_FOR_FLUIDS_PROCESS_H_INCLUDED )
#define  KRATOS_REMOVE_MESH_NODES_FOR_FLUIDS_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "includes/model_part.h"
#include "custom_utilities/mesh_error_calculation_utilities.hpp"
#include "custom_utilities/modeler_utilities.hpp"

#include "pfem_fluid_dynamics_application_variables.h"

///VARIABLES used:
//Data:     NORMAL, MASTER_NODES, NEIGHBOUR_NODES, NEIGBOUR_ELEMENTS
//StepData: MEAN_ERROR, CONTACT_FORCE
//Flags:    (checked) TO_ERASE, BOUNDARY, STRUCTURE, TO_SPLIT, CONTACT, NEW_ENTITY, BLOCKED
//          (set)     TO_ERASE(conditions,nodes)(set), NEW_ENTITY(conditions,nodes)(set), BLOCKED(nodes)->locally, VISITED(nodes)(set)
//          (modified)  
//          (reset)   BLOCKED->locally
//(set):=(set in this process)

namespace Kratos
{

///@name Kratos Classes
///@{

/// Remove Mesh Nodes Process for 2D and 3D cases
/** The process labels the nodes to be erased (TO_ERASE)
    if they are too close (mRemoveOnDistance == true)
    if the error of the patch they belong is very small (REMOVE_NODES_ON_ERROR)
    In the interior of the domain or in the boundary (REMOVE_BOUNDARY_NODES) ...

    Additional treatment of the nonconvex boundaries is also going to erase nodes.

    At the end of the execution nodes are cleaned (only in the current mesh)
    If boundary nodes are removed, conditions must be build again (new conditions are build in the current mesh)   
*/

class RemoveMeshNodesForFluidsProcess
  : public Process
{
public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of Process
  KRATOS_CLASS_POINTER_DEFINITION( RemoveMeshNodesForFluidsProcess );

  typedef ModelPart::ConditionType         ConditionType;
  typedef ModelPart::PropertiesType       PropertiesType;
  typedef ConditionType::GeometryType       GeometryType;
  typedef Bucket<3, Node<3>, std::vector<Node<3>::Pointer>, Node<3>::Pointer, std::vector<Node<3>::Pointer>::iterator, std::vector<double>::iterator > BucketType;
  typedef Tree< KDTreePartition<BucketType> >                          KdtreeType; //Kdtree
  typedef ModelPart::MeshType::GeometryType::PointsArrayType      PointsArrayType;
  ///@}
  ///@name Life Cycle
  ///@{

    /// Default constructor.
    RemoveMeshNodesForFluidsProcess(ModelPart& rModelPart,
			   ModelerUtilities::MeshingParameters& rRemeshingParameters,
			   int EchoLevel) 
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    {
      std::cout<<" remove_mesh_nodes_process_for_fluids "<<std::endl;

      mMeshId = mrRemesh.MeshId;
      mEchoLevel = EchoLevel;
    }


  /// Destructor.
  virtual ~RemoveMeshNodesForFluidsProcess() {}


  ///@}
  ///@name Operators
  ///@{

  /// This operator is provided to call the process as a function and simply calls the Execute method.
  void operator()()
  {
    Execute();
  }


  ///@}
  ///@name Operations
  ///@{


  /// Execute method is used to execute the Process algorithms.
  virtual void Execute()
  {

    KRATOS_TRY

      if( mEchoLevel > 0 ){
	std::cout<<" [ REMOVE CLOSE NODES: "<<std::endl;
	//std::cout<<"   Nodes before erasing : "<<mrModelPart.Nodes(mMeshId).size()<<std::endl;
      }

    double RemovedConditions = mrModelPart.NumberOfConditions(mMeshId);
    double NumberOfNodes = mrModelPart.NumberOfNodes(mMeshId);

    bool any_node_removed      = false;
    bool any_condition_removed = false;
    
    int error_nodes_removed    = 0;
    int inside_nodes_removed   = 0;
    int boundary_nodes_removed = 0;


    //if the remove_node switch is activated, we check if the nodes got too close
    if (mrRemesh.Refine->RemovingOptions.Is(ModelerUtilities::REMOVE_NODES))
      {
	std::cout<<" REMOVE_NODES is TRUE "<<std::endl;
	bool any_node_removed_on_error = false;
	////////////////////////////////////////////////////////////
	if (mrRemesh.Refine->RemovingOptions.Is(ModelerUtilities::REMOVE_NODES_ON_ERROR))	      
	  {
	    std::cout<<" REMOVE_NODES_ON_ERROR is TRUE "<<std::endl;
	  
	    any_node_removed_on_error = RemoveNodesOnError(error_nodes_removed); //2D and 3D
	  }
	//////////////////////////////////////////////////////////// 
	std::cout<<"error_nodes_removed :"<<error_nodes_removed<<std::endl;
	bool any_convex_condition_removed = false;

	bool any_node_removed_on_distance = false;
	////////////////////////////////////////////////////////////
	if (mrRemesh.Refine->RemovingOptions.Is(ModelerUtilities::REMOVE_NODES_ON_DISTANCE))	      
	  {
	    std::cout<<" REMOVE_NODES_ON_DISTANCE is TRUE "<<std::endl;
	    // double  MeanRadius=0;
	    any_node_removed_on_distance = RemoveNodesOnDistance(inside_nodes_removed, boundary_nodes_removed, any_condition_removed); 
	  }
	// REMOVE ON DISTANCE
	////////////////////////////////////////////////////////////
	  
	if(any_node_removed_on_error || any_node_removed_on_distance)
	  any_node_removed = true;

	if(any_convex_condition_removed || any_condition_removed)
	  any_condition_removed = true;
	  

	if(any_node_removed)
	  this->CleanRemovedNodes(mrModelPart,mMeshId);

	if(any_condition_removed){
	  //Clean Conditions
	  ModelPart::ConditionsContainerType RemoveConditions;

	  //id = 0;
	  for(ModelPart::ConditionsContainerType::iterator ic = mrModelPart.ConditionsBegin(mMeshId); ic!= mrModelPart.ConditionsEnd(mMeshId); ic++)
	    {

	      if(ic->IsNot(TO_ERASE)){
		//id+=1;
		RemoveConditions.push_back(*(ic.base()));
		//RemoveConditions.back().SetId(id);
	      }
	      else{
		std::cout<<"   Condition RELEASED:"<<ic->Id()<<std::endl;
	      }
	    }

          mrModelPart.Conditions(mMeshId).swap(RemoveConditions);

	}


      }

      
    // number of removed nodes:
    mrRemesh.Info->RemovedNodes = NumberOfNodes - mrModelPart.NumberOfNodes(mMeshId);
    int distance_remove =  inside_nodes_removed + boundary_nodes_removed;
      
    RemovedConditions -= mrModelPart.NumberOfConditions(mMeshId);
      
    if( mEchoLevel > 0 ){
      std::cout<<"  remove_mesh_nodes_process_for_fluids  [ CONDITIONS ( removed : "<<RemovedConditions<<" ) ]"<<std::endl;
      std::cout<<"   [ NODES      ( removed : "<<mrRemesh.Info->RemovedNodes<<" ) ]"<<std::endl;
      std::cout<<"   [ Error(removed: "<<error_nodes_removed<<"); Distance(removed: "<<distance_remove<<"; inside: "<<inside_nodes_removed<<"; boundary: "<<boundary_nodes_removed<<") ]"<<std::endl;


      //std::cout<<"   Nodes after  erasing : "<<mrModelPart.Nodes(mMeshId).size()<<std::endl;
      std::cout<<"   REMOVE CLOSE NODES ]; "<<std::endl;
    }

    KRATOS_CATCH(" ")
      }


  /// this function is designed for being called at the beginning of the computations
  /// right after reading the model and the groups
  virtual void ExecuteInitialize()
  {
  }

  /// this function is designed for being execute once before the solution loop but after all of the
  /// solvers where built
  virtual void ExecuteBeforeSolutionLoop()
  {
  }

  /// this function will be executed at every time step BEFORE performing the solve phase
  virtual void ExecuteInitializeSolutionStep()
  {	
  }

  /// this function will be executed at every time step AFTER performing the solve phase
  virtual void ExecuteFinalizeSolutionStep()
  {
  }

  /// this function will be executed at every time step BEFORE  writing the output
  virtual void ExecuteBeforeOutputStep()
  {
  }

  /// this function will be executed at every time step AFTER writing the output
  virtual void ExecuteAfterOutputStep()
  {
  }

  /// this function is designed for being called at the end of the computations
  /// right after reading the model and the groups
  virtual void ExecuteFinalize()
  {
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
    return "RemoveMeshNodesForFluidsProcess";
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream& rOStream) const
  {
    rOStream << "RemoveMeshNodesForFluidsProcess";
  }

  /// Print object's data.
  virtual void PrintData(std::ostream& rOStream) const
  {
  }


  ///@}
  ///@name Friends
  ///@{

  ///@}


private:
  ///@name Static Member Variables
  ///@{

  ///@}
  ///@name Static Member Variables
  ///@{
  ModelPart& mrModelPart;
 
  ModelerUtilities::MeshingParameters& mrRemesh;

  ModelerUtilities mModelerUtilities;  

  ModelPart::IndexType mMeshId; 

  int mEchoLevel;

  ///@}
  ///@name Un accessible methods
  ///@{

  //**************************************************************************
  //**************************************************************************

  void CleanRemovedNodes(ModelPart& rModelPart,ModelPart::IndexType MeshId)
  {
    KRATOS_TRY

      //MESH 0 total domain mesh
      ModelPart::NodesContainerType temporal_nodes;
    temporal_nodes.reserve(rModelPart.Nodes(MeshId).size());
	
    temporal_nodes.swap(rModelPart.Nodes(MeshId));
	
    for(ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; i_node++)
      {
	if( i_node->IsNot(TO_ERASE) ){
	  (rModelPart.Nodes(MeshId)).push_back(*(i_node.base()));	
//	  if(i_node->Id()>24880){
//	  std::cout<<" NOT ERASED NODE: "<<i_node->Id()<<"  COORDINATES "<<i_node->Coordinates()<<std::endl;
//		  }	 
	
	}
	// else{
	//     // std::cout<<"   ERASED NODE: "<<i_node->Id()<<"  COORDINATES "<<i_node->Coordinates()<<std::endl;
	//   if( i_node->Is(BOUNDARY) )
	//     std::cout<<"   BOUNDARY NODE RELEASED "<<i_node->Id()<<std::endl;
	// }
      }
	
    rModelPart.Nodes(MeshId).Sort();
	
	
    KRATOS_CATCH( "" )
      }

  //**************************************************************************
  //**************************************************************************

  bool RemoveNodesOnError(int& error_removed_nodes)
  {
    KRATOS_TRY
	 
    //***SIZES :::: parameters do define the tolerance in mesh size: 
    double size_for_criterion_error   = 2.0 * mrRemesh.Refine->CriticalRadius; //compared with mean node radius
       
    bool any_node_removed = false;

    MeshErrorCalculationUtilities MeshErrorDistribution;
    MeshErrorDistribution.SetEchoLevel(mEchoLevel);

    std::vector<double> NodalError;
    std::vector<int>    nodes_ids;

	      
    MeshErrorDistribution.NodalErrorCalculation(mrModelPart,NodalError,nodes_ids,mMeshId,mrRemesh.Refine->GetErrorVariable());

    for(ModelPart::NodesContainerType::const_iterator in = mrModelPart.NodesBegin(mMeshId); in != mrModelPart.NodesEnd(mMeshId); in++)
      {

	WeakPointerVector<Node<3> >& rN = in->GetValue(NEIGHBOUR_NODES);
	int erased_nodes =0;
	for(unsigned int i = 0; i < rN.size(); i++)
	  {
	    if(rN[i].Is(TO_ERASE))
	      erased_nodes += 1;
	  }
	   
	   
	if( in->IsNot(BOUNDARY) &&  in->IsNot(STRUCTURE) && erased_nodes < 1 )
	  {
	    double& MeanError = in->FastGetSolutionStepValue(MEAN_ERROR);
	    MeanError = NodalError[nodes_ids[in->Id()]];
	       
	    WeakPointerVector<Element >& neighb_elems = in->GetValue(NEIGHBOUR_ELEMENTS);
	    double mean_node_radius = 0;
	    for(WeakPointerVector< Element >::iterator ne = neighb_elems.begin(); ne!=neighb_elems.end(); ne++)
	      {
		mean_node_radius+= mModelerUtilities.CalculateElementRadius(ne->GetGeometry()); //Triangle 2D, Tetrahedron 3D
		//mean_node_radius+= mModelerUtilities.CalculateTriangleRadius(ne->GetGeometry());
		//mean_node_radius+= mModelerUtilities.CalculateTetrahedronRadius(ne->GetGeometry());
	      }
	       
	    mean_node_radius /= double(neighb_elems.size());
	       
	    if(NodalError[nodes_ids[in->Id()]] < mrRemesh.Refine->ReferenceError && mean_node_radius < size_for_criterion_error)
	      {
		//std::cout<<"   Energy : node remove ["<<in->Id()<<"] : "<<NodalError[nodes_ids[in->Id()]]<<std::endl;
		//std::cout<<"   mean_node_radius "<<mean_node_radius<<" < "<<size_for_criterion_error<<" size_for_criterion_error"<<std::endl;
		in->Set(TO_ERASE);
		any_node_removed = true;
		error_removed_nodes++;
	      }
	  }
      }
       
    return any_node_removed;

    KRATOS_CATCH(" ")

      }



  bool RemoveNodesOnDistance(int& inside_nodes_removed, int& boundary_nodes_removed, bool& any_condition_removed)
  {
    KRATOS_TRY
	 

    //***SIZES :::: parameters do define the tolerance in mesh size: 
    double initialMeanRadius=0;
    initialMeanRadius=mrRemesh.Refine->InitialRadius;
    std::cout<<"initialMeanRadius "<<initialMeanRadius<<std::endl;
    double size_for_distance_inside       = 1.25  * mrRemesh.Refine->CriticalRadius;//compared to element radius
    double size_for_distance_boundary     = 1.25  * size_for_distance_inside; //compared to element radius
    double size_for_wall_tip_contact_side = 0.15 * mrRemesh.Refine->CriticalSide;
    size_for_distance_inside       = 0.6 * initialMeanRadius;//compared to element radius
    size_for_distance_boundary     = 0.6 * initialMeanRadius; //compared to element radius
    size_for_wall_tip_contact_side = 0.15 * mrRemesh.Refine->CriticalSide;

 
    bool derefine_wall_tip_contact = false;

    bool any_node_removed = false;

    //bucket size definition:
    unsigned int bucket_size = 20;
       
    //create the list of the nodes to be check during the search
    std::vector<Node<3>::Pointer> list_of_nodes;
    list_of_nodes.reserve(mrModelPart.NumberOfNodes(mMeshId));
    for(ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin(mMeshId) ; i_node != mrModelPart.NodesEnd(mMeshId) ; i_node++)
      {
	(list_of_nodes).push_back(*(i_node.base()));
      }
       
    KdtreeType nodes_tree(list_of_nodes.begin(),list_of_nodes.end(), bucket_size);

    ////////////////////////////////////////////////////////////
       
    //all of the nodes in this list will be preserved
    unsigned int num_neighbours = 100;
       
    std::vector<Node<3>::Pointer> neighbours         (num_neighbours);
    std::vector<double>           neighbour_distances(num_neighbours);
       
       
    //radius means the distance, if the distance between two nodes is closer to radius -> mark for removing
    double radius=0;
    Node<3> work_point(0,0.0,0.0,0.0);
    unsigned int n_points_in_radius;
       
       
    for(ModelPart::NodesContainerType::const_iterator in = mrModelPart.NodesBegin(mMeshId); in != mrModelPart.NodesEnd(mMeshId); in++)
      {
	bool on_contact_tip = false;
	bool contact_active = false;
	   
	if( in->SolutionStepsDataHas(CONTACT_FORCE) ){
	  array_1d<double, 3 > & ContactForceNormal  = in->FastGetSolutionStepValue(CONTACT_FORCE);
	  if(norm_2(ContactForceNormal)>0)
	    contact_active = true;
	}

	if(contact_active || in->Is(TO_SPLIT) || in->Is(CONTACT) )
	  on_contact_tip = true;				  

	if( in->IsNot(NEW_ENTITY) )
	  {
	    radius = size_for_distance_inside;

	    work_point[0]=in->X();
	    work_point[1]=in->Y();
	    work_point[2]=in->Z();
	    if(in->Is(FREE_SURFACE)){// it must be more difficult to erase a free_surface node, otherwise, lot of volume is lost
	      // radius = 0.5  * initialMeanRadius;//compared with element radius
	      radius = 0.3  * initialMeanRadius;//compared with element radius
       
	      // radius*=0.85;
	    }
	    n_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, neighbours.begin(),neighbour_distances.begin(), num_neighbours);

	    if (n_points_in_radius>1)
	      {
		// std::cout<<"     Points in Radius "<< n_points_in_radius<<" radius "<<radius<<std::endl;

		//if( in->IsNot(STRUCTURE) ) {//MEANS DOFS FIXED

		// if ( in->IsNot(BOUNDARY) && in->IsNot(FREE_SURFACE) && in->IsNot(RIGID) && in->IsNot(ISOLATED))
		if ( in->IsNot(RIGID))
		  {
		    if( mrRemesh.Refine->RemovingOptions.Is(ModelerUtilities::REMOVE_NODES_ON_DISTANCE) ){
		      //look if we are already erasing any of the other nodes
		      unsigned int contact_nodes = 0;
		      unsigned int erased_nodes = 0;
		      for(std::vector<Node<3>::Pointer>::iterator nn=neighbours.begin(); nn!=neighbours.begin() + n_points_in_radius ; nn++)
			{
			  if( (*nn)->Is(BOUNDARY) && (*nn)->Is(CONTACT) )
			    contact_nodes += 1;

			  if( (*nn)->Is(TO_ERASE) )
			    erased_nodes += 1;
			}

		      if( erased_nodes < 1 && contact_nodes < 1){ //we release the node if no other nodes neighbours are being erased
			in->Set(TO_ERASE);
			// std::cout<<"     Distance Criterion Node ["<<in->Id()<<"] TO_ERASE "<<std::endl;
			any_node_removed = true;
			inside_nodes_removed++;
			//distance_remove++;
		      }

		    }

		  }
		// else if ( (mrRemesh.Refine->RemovingOptions.Is(ModelerUtilities::REMOVE_BOUNDARY_NODES) && mrRemesh.Refine->RemovingOptions.Is(ModelerUtilities::REMOVE_BOUNDARY_NODES_ON_DISTANCE)) && (in)->IsNot(TO_ERASE)) //boundary nodes will be removed if they get REALLY close to another boundary node (0.2(=extra_factor) * h_factor)
		else     {

		  // std::cout<<"  Remove close boundary nodes: Candidate ["<<in->Id()<<"]"<<std::endl;

		  //here we loop over the neighbouring nodes and if there are nodes
		  //with BOUNDARY flag and closer than 0.2*nodal_h from our node, we remove the node we are considering
		  unsigned int k = 0;
		  unsigned int counter = 0;
		  for(std::vector<Node<3>::Pointer>::iterator nn=neighbours.begin(); nn!=neighbours.begin() + n_points_in_radius ; nn++)
		    {
		      bool nn_on_contact_tip = false;
		      bool contact_active = false;
	   
		      if( (*nn)->SolutionStepsDataHas(CONTACT_FORCE) ){
			array_1d<double, 3 > & ContactForceNormal  = (*nn)->FastGetSolutionStepValue(CONTACT_FORCE);
			if(norm_2(ContactForceNormal)>0)
			  contact_active = true;
		      }
			   
		      if(contact_active || (*nn)->Is(TO_SPLIT) || (*nn)->Is(CONTACT) )
			nn_on_contact_tip = true;				  

		      //std::cout<<" radius * extra_factor "<<(extra_factor*radius)<<" >? "<<neighbour_distances[k]<<std::endl;
		      if ( (*nn)->Is(BOUNDARY) && !nn_on_contact_tip && neighbour_distances[k] < size_for_distance_boundary && neighbour_distances[k] > 0.0 )
			{
			  //KRATOS_WATCH( neighbours_distances[k] )
			  if((*nn)->IsNot(TO_ERASE)){
			    counter += 1;
			  }
			}

		      if ( (*nn)->Is(BOUNDARY) && nn_on_contact_tip && neighbour_distances[k] < size_for_wall_tip_contact_side ) {
			if ( (*nn)->IsNot(TO_ERASE)) { 
			  counter += 1;
			}
		      }


		      k++;
		    }

		  if(counter > 1 && in->IsNot(NEW_ENTITY) && !on_contact_tip ){ //Can be inserted in the boundary refine
		    in->Set(TO_ERASE);
		    std::cout<<"     Removed Boundary Node ["<<in->Id()<<"] on Distance "<<std::endl;
		    any_node_removed = true;
		    boundary_nodes_removed++;
		    //distance_remove ++;
		  }
		  else if ( counter > 2 && in->IsNot(NEW_ENTITY) && on_contact_tip && derefine_wall_tip_contact) {
		    in->Set(TO_ERASE);
		    std::cout << "     Removing a TIP POINT due to that criterion [" << in->Id() << "]" << std::endl;
		    any_node_removed = true;
		    boundary_nodes_removed++;
		  }

		}

		//}

	      }

	  }	
      }

    uint erased_nodes=0;
    for(ModelPart::ElementsContainerType::const_iterator ie = mrModelPart.ElementsBegin(mMeshId);
	ie != mrModelPart.ElementsEnd(mMeshId); ie++)
      {
	const unsigned int dimension = ie->GetGeometry().WorkingSpaceDimension();
	uint rigidNodes=0;
	//coordinates
	for(unsigned int i=0; i<ie->GetGeometry().size(); i++)
	  {	
	    if(ie->GetGeometry()[i].Is(RIGID)){
	      rigidNodes++;
	    }
	  }
	if(dimension==2){

	  if(rigidNodes>1){
	    // std::cout<<"erased_nodes "<<erased_nodes<<std::endl;
	    double distWall=0;
	    double safetyCoefficient2D=0.5;
	    
	    double elementVolume=ie->GetGeometry().Volume();
	    double CriticalVolume=0.1*mrRemesh.Refine->MeanVolume;
	    if(elementVolume<CriticalVolume){
	      for(unsigned int i=0; i<ie->GetGeometry().size(); i++)
		{	
		  if(ie->GetGeometry()[i].IsNot(RIGID) && ie->GetGeometry()[i].IsNot(TO_ERASE)){
		    ie->GetGeometry()[i].Set(TO_ERASE);
		    std::cout<<"erase this layer node because it may be potentially dangerous and pass through the solid contour"<<std::endl;
		    erased_nodes += 1;
		    break;
		  }
		}
	   
	    }
	    array_1d<double,3> Edges(3,0.0);
	    array_1d<uint,3> FirstEdgeNode(3,0);
	    array_1d<uint,3> SecondEdgeNode(3,0);
	    double WallCharacteristicDistance=0;
	    array_1d<double,2> CoorDifference(3,0.0);
	    CoorDifference = ie->GetGeometry()[1].Coordinates() - ie->GetGeometry()[0].Coordinates();
	    double SquaredLength = CoorDifference[0]*CoorDifference[0] + CoorDifference[1]*CoorDifference[1];
	    Edges[0]=sqrt(SquaredLength);
	    FirstEdgeNode[0]=0;
	    SecondEdgeNode[0]=1;
	    if(ie->GetGeometry()[0].Is(RIGID) && ie->GetGeometry()[1].Is(RIGID)){
	      WallCharacteristicDistance=Edges[0];
	    }
	    uint Counter=0;
	    for (uint i = 2; i < ie->GetGeometry().size(); i++){
	      for(uint j = 0; j < i; j++)
		{
		  CoorDifference = ie->GetGeometry()[i].Coordinates() - ie->GetGeometry()[j].Coordinates();
		  SquaredLength = CoorDifference[0]*CoorDifference[0] + CoorDifference[1]*CoorDifference[1];
		  Counter+=1;
		  Edges[Counter]=sqrt(SquaredLength);
		  FirstEdgeNode[Counter]=j;
		  SecondEdgeNode[Counter]=i;
		  if(ie->GetGeometry()[i].Is(RIGID) && ie->GetGeometry()[j].Is(RIGID) && Edges[Counter]>WallCharacteristicDistance ){
		    WallCharacteristicDistance=Edges[Counter];
		  }
		}

	    }
	    distWall=WallCharacteristicDistance;

	    for (uint i = 0; i < 3; i++){
	      if(((ie->GetGeometry()[FirstEdgeNode[i]].Is(RIGID) && !ie->GetGeometry()[SecondEdgeNode[i]].Is(RIGID)) ||
		  (ie->GetGeometry()[SecondEdgeNode[i]].Is(RIGID) && !ie->GetGeometry()[FirstEdgeNode[i]].Is(RIGID))) &&
		 ie->GetGeometry()[FirstEdgeNode[i]].IsNot(TO_ERASE) &&
		 ie->GetGeometry()[SecondEdgeNode[i]].IsNot(TO_ERASE)&&
		 Edges[i]<safetyCoefficient2D*distWall){
		if(ie->GetGeometry()[FirstEdgeNode[i]].IsNot(RIGID) && ie->GetGeometry()[FirstEdgeNode[i]].IsNot(TO_ERASE)){
		  ie->GetGeometry()[FirstEdgeNode[i]].Set(TO_ERASE);
		  // std::cout<<"erase node 0"<<std::endl;
		  erased_nodes += 1;
		}else if(ie->GetGeometry()[SecondEdgeNode[i]].IsNot(RIGID) && ie->GetGeometry()[SecondEdgeNode[i]].IsNot(TO_ERASE)){
		  ie->GetGeometry()[SecondEdgeNode[i]].Set(TO_ERASE);
		  // std::cout<<"erase node 1"<<std::endl;
		  erased_nodes += 1;
		}

	      }

	    }


	  }



	}else if(dimension==3){
	  // std::cout<<"boundary_nodes_removed "<<std::endl;

	  if(rigidNodes>1){
	    // std::cout<<"erased_nodes "<<erased_nodes<<std::endl;
	    double distWall=0;
	    double safetyCoefficient3D=0.6;
	    
	    double elementVolume=ie->GetGeometry().Volume();
	    double CriticalVolume=0.1*mrRemesh.Refine->MeanVolume;
	    if(elementVolume<CriticalVolume){
	      for(unsigned int i=0; i<ie->GetGeometry().size(); i++)
		{
		  if(ie->GetGeometry()[i].IsNot(RIGID) && ie->GetGeometry()[i].IsNot(TO_ERASE)){
		    ie->GetGeometry()[i].Set(TO_ERASE);
		    std::cout<<"erase this layer node because it may be potentially dangerous and pass through the solid contour"<<std::endl;
		    erased_nodes += 1;
		    break;
		  }
		}
	   
	    }
	    array_1d<double,6> Edges(6,0.0);
	    array_1d<uint,6> FirstEdgeNode(6,0);
	    array_1d<uint,6> SecondEdgeNode(6,0);
	    double WallCharacteristicDistance=0;
	    array_1d<double,3> CoorDifference(3,0.0);
	    CoorDifference = ie->GetGeometry()[1].Coordinates() - ie->GetGeometry()[0].Coordinates();
	    double SquaredLength = CoorDifference[0]*CoorDifference[0] +
	      CoorDifference[1]*CoorDifference[1] +
	      CoorDifference[2]*CoorDifference[2];
	    Edges[0]=sqrt(SquaredLength);
	    FirstEdgeNode[0]=0;
	    SecondEdgeNode[0]=1;
	    if(ie->GetGeometry()[0].Is(RIGID) && ie->GetGeometry()[1].Is(RIGID)){
	      WallCharacteristicDistance=Edges[0];
	    }
	    uint Counter=0;
	    for (uint i = 2; i < ie->GetGeometry().size(); i++){
	      for(uint j = 0; j < i; j++)
	    	{
	    	  CoorDifference = ie->GetGeometry()[i].Coordinates() - ie->GetGeometry()[j].Coordinates();
	    	  SquaredLength = CoorDifference[0]*CoorDifference[0] +
	    	    CoorDifference[1]*CoorDifference[1] +
	    	    CoorDifference[2]*CoorDifference[2];
	    	  Counter+=1;
	    	  Edges[Counter]=sqrt(SquaredLength);
	    	  FirstEdgeNode[Counter]=j;
	    	  SecondEdgeNode[Counter]=i;
	    	  // if(ie->GetGeometry()[i].Is(RIGID) && ie->GetGeometry()[j].Is(RIGID) && Edges[Counter]>WallCharacteristicDistance ){
	    	  if(ie->GetGeometry()[i].Is(RIGID) && ie->GetGeometry()[j].Is(RIGID) && WallCharacteristicDistance==0 ){
	    	    WallCharacteristicDistance=Edges[Counter];
	    	  }
	    	}

	    }
	    distWall=WallCharacteristicDistance;


	    for (uint i = 0; i < 6; i++){
	      if(((ie->GetGeometry()[FirstEdgeNode[i]].Is(RIGID) && !ie->GetGeometry()[SecondEdgeNode[i]].Is(RIGID)) ||
	    	  (ie->GetGeometry()[SecondEdgeNode[i]].Is(RIGID) && !ie->GetGeometry()[FirstEdgeNode[i]].Is(RIGID))) &&
	    	 ie->GetGeometry()[FirstEdgeNode[i]].IsNot(TO_ERASE) &&
	    	 ie->GetGeometry()[SecondEdgeNode[i]].IsNot(TO_ERASE)&&
	    	 Edges[i]<safetyCoefficient3D*distWall){
	    	// if(ie->GetGeometry()[FirstEdgeNode[i]].IsNot(RIGID) && ie->GetGeometry()[FirstEdgeNode[i]].IsNot(TO_ERASE) && ie->GetGeometry()[FirstEdgeNode[i]].IsNot(FREE_SURFACE)){
	    	if(ie->GetGeometry()[FirstEdgeNode[i]].IsNot(RIGID) && ie->GetGeometry()[FirstEdgeNode[i]].IsNot(TO_ERASE)){

	    	  ie->GetGeometry()[FirstEdgeNode[i]].Set(TO_ERASE);
	    	  // std::cout<<"erase node 0"<<std::endl;
	    	  erased_nodes += 1;
	    	// }else if(ie->GetGeometry()[SecondEdgeNode[i]].IsNot(RIGID) && ie->GetGeometry()[SecondEdgeNode[i]].IsNot(TO_ERASE) && ie->GetGeometry()[SecondEdgeNode[i]].IsNot(FREE_SURFACE)){
	    	}else if(ie->GetGeometry()[SecondEdgeNode[i]].IsNot(RIGID) && ie->GetGeometry()[SecondEdgeNode[i]].IsNot(TO_ERASE)){
	    	  ie->GetGeometry()[SecondEdgeNode[i]].Set(TO_ERASE);
	    	  // std::cout<<"erase node 1"<<std::endl;
	    	  erased_nodes += 1;
	    	}

	      }

	    }



	  }

	}

      }
 

    if(erased_nodes>0){
      std::cout<<"layer_nodes_removed "<<erased_nodes<<std::endl;
      any_node_removed = true;
    }
    //Build boundary after removing boundary nodes due distance criterion
    std::cout<<"boundary_nodes_removed "<<boundary_nodes_removed<<std::endl;
    std::cout<<"inside_nodes_removed "<<inside_nodes_removed<<std::endl;

    return any_node_removed;
       
    KRATOS_CATCH(" ")

      } 

    /// Assignment operator.
    RemoveMeshNodesForFluidsProcess& operator=(RemoveMeshNodesForFluidsProcess const& rOther);


    /// this function is a private function


    /// Copy constructor.
    //Process(Process const& rOther);


    ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  RemoveMeshNodesForFluidsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const RemoveMeshNodesForFluidsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_REMOVE_MESH_NODES_FOR_FLUIDS_PROCESS_H_INCLUDED  defined 


