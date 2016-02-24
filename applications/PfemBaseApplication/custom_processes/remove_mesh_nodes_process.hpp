//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_REMOVE_MESH_NODES_PROCESS_H_INCLUDED )
#define  KRATOS_REMOVE_MESH_NODES_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "includes/model_part.h"
#include "custom_utilities/mesh_error_calculation_utilities.hpp"
#include "custom_utilities/modeler_utilities.hpp"

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

/// Remove Mesh Nodes Process
/** The process labels the nodes to be erased (TO_ERASE)
    if they are too close (mRemoveOnDistance == true)
    if the error of the patch they belong is very small (REMOVE_NODES_ON_ERROR)
    In the interior of the domain or in the boundary (REMOVE_BOUNDARY_NODES) ...

    Additional treatment of the nonconvex boundaries is also going to erase nodes.

    At the end of the execution nodes are cleaned (only in the current mesh)
    If boundary nodes are removed, conditions must be build again (new conditions are build in the current mesh)   
*/

class RemoveMeshNodesProcess
  : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION( RemoveMeshNodesProcess );

    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;

    typedef Bucket<3, Node<3>, std::vector<Node<3>::Pointer>, Node<3>::Pointer, std::vector<Node<3>::Pointer>::iterator, std::vector<double>::iterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> >                          KdtreeType; //Kdtree
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RemoveMeshNodesProcess(ModelPart& rModelPart,
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
    virtual ~RemoveMeshNodesProcess() {}


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
	  bool any_node_removed_on_error = false;
	  ////////////////////////////////////////////////////////////
	  if (mrRemesh.Refine->RemovingOptions.Is(ModelerUtilities::REMOVE_NODES_ON_ERROR))	      
	    {
	      any_node_removed_on_error = RemoveNodesOnError(error_nodes_removed);
	    }
	  //////////////////////////////////////////////////////////// 

	  bool any_convex_condition_removed = false;
	  ////////////////////////////////////////////////////////////
	  if (mrRemesh.Refine->RemovingOptions.Is(ModelerUtilities::REMOVE_BOUNDARY_NODES))
	    {
	      any_convex_condition_removed = RemoveNonConvexBoundary();
	    }
	  //////////////////////////////////////////////////////////// 


	  bool any_node_removed_on_distance = false;
	  ////////////////////////////////////////////////////////////
	  if (mrRemesh.Refine->RemovingOptions.Is(ModelerUtilities::REMOVE_NODES_ON_DISTANCE) || mrRemesh.Refine->RemovingOptions.Is(ModelerUtilities::REMOVE_BOUNDARY_NODES_ON_DISTANCE))	      
	    {
	      any_node_removed_on_distance = RemoveNodesOnDistance(inside_nodes_removed, boundary_nodes_removed, any_condition_removed);
	    }
	  // REMOVE ON DISTANCE
	  ////////////////////////////////////////////////////////////
	  
	  if(any_node_removed_on_error || any_node_removed_on_distance)
	    any_node_removed = true;

	  if(any_convex_condition_removed || any_condition_removed)
	    any_condition_removed = true;
	  

	  if(any_node_removed)
	    mModelerUtilities.CleanRemovedNodes(mrModelPart,mMeshId);

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
	std::cout<<"   [ CONDITIONS ( removed : "<<RemovedConditions<<" ) ]"<<std::endl;
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
        return "RemoveMeshNodesProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "RemoveMeshNodesProcess";
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
		   mean_node_radius+= mModelerUtilities.CalculateTriangleRadius(ne->GetGeometry());
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

    //**************************************************************************
    //**************************************************************************

    bool RemoveNodesOnDistance(int& inside_nodes_removed, int& boundary_nodes_removed, bool& any_condition_removed)
    {
       KRATOS_TRY
	 

       //***SIZES :::: parameters do define the tolerance in mesh size: 
       double size_for_distance_inside       = 1.0  * mrRemesh.Refine->CriticalRadius; //compared with element radius
       double size_for_distance_boundary     = 1.5  * size_for_distance_inside; //compared with element radius
       double size_for_wall_tip_contact_side = 0.15 * mrRemesh.Refine->CriticalSide;
 
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
	   array_1d<double, 3 > & ContactForceNormal  = in->FastGetSolutionStepValue(CONTACT_FORCE);

	   if(norm_2(ContactForceNormal)>0 || in->Is(TO_SPLIT) || in->Is(CONTACT) )
	     on_contact_tip = true;				  

	   if( in->IsNot(NEW_ENTITY) )
	     {
	       radius = size_for_distance_inside;

	       work_point[0]=in->X();
	       work_point[1]=in->Y();
	       work_point[2]=in->Z();

	       n_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, neighbours.begin(),neighbour_distances.begin(), num_neighbours);

	       if (n_points_in_radius>1)
		 {
		   //std::cout<<"     Points in Radius "<< n_points_in_radius<<" radius "<<radius<<std::endl;

		   //if( in->IsNot(STRUCTURE) ) {//MEANS DOFS FIXED

		   if ( in->IsNot(BOUNDARY) )
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
			   //std::cout<<"     Distance Criterion Node ["<<in->Id()<<"] TO_ERASE "<<std::endl;
			   any_node_removed = true;
			   inside_nodes_removed++;
			   //distance_remove++;
			 }

		       }

		     }
		   else if ( (mrRemesh.Refine->RemovingOptions.Is(ModelerUtilities::REMOVE_BOUNDARY_NODES) && mrRemesh.Refine->RemovingOptions.Is(ModelerUtilities::REMOVE_BOUNDARY_NODES_ON_DISTANCE)) && (in)->IsNot(TO_ERASE)) //boundary nodes will be removed if they get REALLY close to another boundary node (0.2(=extra_factor) * h_factor)
		     {

		       //std::cout<<"  Remove close boundary nodes: Candidate ["<<in->Id()<<"]"<<std::endl;

		       //here we loop over the neighbouring nodes and if there are nodes
		       //with BOUNDARY flag and closer than 0.2*nodal_h from our node, we remove the node we are considering
		       unsigned int k = 0;
		       unsigned int counter = 0;
		       for(std::vector<Node<3>::Pointer>::iterator nn=neighbours.begin(); nn!=neighbours.begin() + n_points_in_radius ; nn++)
			 {
			   bool nn_on_contact_tip = false;
			   array_1d<double, 3 > & ContactForceNormal  = (*nn)->FastGetSolutionStepValue(CONTACT_FORCE);

			   if(norm_2(ContactForceNormal)>0 || (*nn)->Is(TO_SPLIT) || (*nn)->Is(CONTACT) )
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
			 //std::cout<<"     Removed Boundary Node ["<<in->Id()<<"] on Distance "<<std::endl;
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

       //Build boundary after removing boundary nodes due distance criterion
       if(boundary_nodes_removed){
	 any_condition_removed = RebuildBoundary();
       }
       //Build boundary after removing boundary nodes due distance criterion
       

       return any_node_removed;
       
       KRATOS_CATCH(" ")

    }

    //**************************************************************************
    //**************************************************************************

    bool RebuildBoundary()
    {
      KRATOS_TRY
	   
     bool any_condition_removed = false;

      std::vector<std::vector<Condition::Pointer> > node_shared_conditions(mrModelPart.NumberOfNodes()+1); //all domain nodes
      
      for(ModelPart::ConditionsContainerType::iterator ic = mrModelPart.ConditionsBegin(mMeshId); ic!= mrModelPart.ConditionsEnd(mMeshId); ic++)
	{	 
	  if(ic->IsNot(NEW_ENTITY) && ic->IsNot(TO_ERASE)){
	    Geometry< Node<3> >& rConditionGeom = ic->GetGeometry();
	    for(unsigned int i=0; i<rConditionGeom.size(); i++){
	      //std::cout<<"["<<ic->Id()<<"] i "<<i<<" condition "<<rConditionGeom[i].Id()<<std::endl;
	      if(rConditionGeom[i].Is(TO_ERASE)){
		if( mEchoLevel > 0 )
		  std::cout<<"     Released node condition ["<<rConditionGeom[i].Id()<<"]: WARNING "<<std::endl;
	      }
	      
	      node_shared_conditions[rConditionGeom[i].Id()].push_back(*(ic.base()));	  
	    }
	  }
	  
	}
      
      
      //nodes
      int i=0,j=0;
      unsigned int initial_cond_size = mrModelPart.Conditions().size()+1; //total model part node size
      unsigned int id = 1;
      unsigned int new_id = 0;
      
      for(ModelPart::NodesContainerType::const_iterator in = mrModelPart.NodesBegin(mMeshId); in != mrModelPart.NodesEnd(mMeshId); in++) 
	{
	  
	  if( in->Is(BOUNDARY) && in->IsNot(BLOCKED) && in->IsNot(NEW_ENTITY) && in->Is(TO_ERASE) ){

	    unsigned int nodeId = in->Id();

	    if(node_shared_conditions[nodeId].size()>=2){

	      // std::cout<<"     nodeId "<<nodeId<<std::endl;
	      if(node_shared_conditions[nodeId][0]->IsNot(TO_ERASE) && node_shared_conditions[nodeId][1]->IsNot(TO_ERASE)){

		if(node_shared_conditions[nodeId][0]->GetGeometry()[0].Id() == in->Id()){
		  i = 1;
		  j = 0;
		}
		else{
		  i = 0;
		  j = 1;
		}


		Geometry< Node<3> >& rConditionGeom1 = node_shared_conditions[nodeId][i]->GetGeometry();
		Geometry< Node<3> >& rConditionGeom2 = node_shared_conditions[nodeId][j]->GetGeometry();

		//node in id Node1;

		Node<3> & Node0 = rConditionGeom1[0]; // other node in condition [1]
		Node<3> & Node2 = rConditionGeom2[1]; // other node in condition [2]

		node_shared_conditions[nodeId][i]->Set(TO_ERASE); //release condition [1]
		node_shared_conditions[nodeId][j]->Set(TO_ERASE); //release condition [2]

		//condition will be removed
		any_condition_removed = true;

		Condition::Pointer NewCond = node_shared_conditions[nodeId][i];

		Node0.Set(BLOCKED);
		Node0.Set(ModelerUtilities::ENGAGED_NODES);

		Node2.Set(BLOCKED);
		Node2.Set(ModelerUtilities::ENGAGED_NODES);

		//create new condition Node0-NodeB
		Condition::NodesArrayType face;
		face.reserve(2);

		face.push_back(rConditionGeom1(0));
		face.push_back(rConditionGeom2(1));

		new_id = initial_cond_size + id;
		//properties to be used in the generation
		Condition::Pointer pcond       = NewCond->Clone(new_id, face);
		// std::cout<<"     ID"<<id<<" 1s "<<pcond1->GetGeometry()[0].Id()<<" "<<pcond1->GetGeometry()[1].Id()<<std::endl;

		pcond->Set(NEW_ENTITY);

		//std::cout<<"     Condition INSERTED (Id: "<<new_id<<") ["<<rConditionGeom1[0].Id()<<", "<<rConditionGeom2[1].Id()<<"] "<<std::endl;

		pcond->SetValue(NORMAL, NewCond->GetValue(NORMAL) );
		pcond->SetValue(MASTER_NODES, NewCond->GetValue(MASTER_NODES) );

		MeshDataTransferUtilities TransferUtilities;
		TransferUtilities.TransferBoundaryData(pcond, NewCond, mrRemesh.TransferVariables); 

		(mrModelPart.Conditions(mMeshId)).push_back(pcond);

		id +=1;
	      }

	    }
	  }

	}

      for(ModelPart::NodesContainerType::const_iterator in = mrModelPart.NodesBegin(mMeshId); in != mrModelPart.NodesEnd(mMeshId); in++)
	{
	  in->Reset(BLOCKED);
	}

      return any_condition_removed;

      KRATOS_CATCH(" ")
    }

    //**************************************************************************
    //**************************************************************************

    bool RemoveNonConvexBoundary()
    {
      KRATOS_TRY

      if( mEchoLevel > 0 ){
	std::cout<<"   [ REMOVE NON CONVEX BOUNDARY : "<<std::endl;
	//std::cout<<"     Starting Conditions : "<<mrModelPart.Conditions(mMeshId).size()<<std::endl;
      }

      double RemovedConditions = mrModelPart.NumberOfConditions(mMeshId);

      //***SIZES :::: parameters do define the tolerance in mesh size: 
      double critical_angle        = -120;
      double size_for_side_normal  =  mrRemesh.Refine->CriticalRadius;
      
      
      std::vector<std::vector<Condition::Pointer> > node_shared_conditions(mrModelPart.NumberOfNodes()+1); //all domain nodes
      
      //std::cout<<"     Shared Conditions Size "<<node_shared_conditions.size()<<std::endl;
      
      for(ModelPart::ConditionsContainerType::iterator ic = mrModelPart.ConditionsBegin(mMeshId); ic!= mrModelPart.ConditionsEnd(mMeshId); ic++)
	{	 
	  if(ic->IsNot(NEW_ENTITY) && ic->IsNot(TO_ERASE)){
	    Geometry< Node<3> >& rConditionGeom = ic->GetGeometry();
	    for(unsigned int i=0; i<rConditionGeom.size(); i++){
	      //std::cout<<"["<<ic->Id()<<"] i "<<i<<" condition "<<rConditionGeom[i].Id()<<std::endl;
	      if(rConditionGeom[i].Is(TO_ERASE))
		std::cout<<"     WARNING: Released node condition "<<std::endl;

	      node_shared_conditions[rConditionGeom[i].Id()].push_back(*(ic.base()));	  
	    }
	  }
	}

      //std::cout<<"     Node Shared Conditions (Pair of Condition Nodes) is now set "<<std::endl;
      
      //angles 
      double condition_angle = 0;
      
      //vector of the neighbour conditions
      array_1d<double,3> S1;
      array_1d<double,3> S2;
      S1.clear();
      S2.clear();
      
      //normals of the neighbour conditions
      array_1d<double,3> N1;
      array_1d<double,3> N2;
      N1.clear();
      N2.clear();
      
      //nodes
      int i=0,j=0;
      
      //condition id and size
      unsigned int initial_cond_size = mrModelPart.Conditions().size()+1; //total model part node size
      unsigned int id = 1;
      unsigned int new_id = 0;
      int RemovedNodes =0;
      
      for(ModelPart::NodesContainerType::const_iterator in = mrModelPart.NodesBegin(mMeshId); in != mrModelPart.NodesEnd(mMeshId); in++)
	{

	  if( in->Is(BOUNDARY) && in->IsNot(BLOCKED) && in->IsNot(NEW_ENTITY) )
	    {
	      unsigned int nodeId = in->Id();

	      if(node_shared_conditions[nodeId].size()>=2){

		// std::cout<<"     nodeId "<<nodeId<<std::endl;
		if(node_shared_conditions[nodeId][0]->IsNot(TO_ERASE) && node_shared_conditions[nodeId][1]->IsNot(TO_ERASE)){
		
		  if(node_shared_conditions[nodeId][0]->GetGeometry()[0].Id() == in->Id()){
		    i = 1;
		    j = 0;
		  }
		  else{
		    i = 0;
		    j = 1;
		  }
	      
		  
		  //Node1*  neighbour conditions in 2D:   (Node0) ---[1]--- (Node1*)

		  //normal condition [1]
		  N1 = node_shared_conditions[nodeId][i]->GetValue(NORMAL);
		  //normal condition [2]
		  N2 = node_shared_conditions[nodeId][j]->GetValue(NORMAL);
	      
		  // std::cout<<"     N1 "<<N1<<std::endl;
		  // std::cout<<"     N2 "<<N2<<std::endl;
		

		  Geometry< Node<3> >& rConditionGeom1 = node_shared_conditions[nodeId][i]->GetGeometry();
		  Geometry< Node<3> >& rConditionGeom2 = node_shared_conditions[nodeId][j]->GetGeometry();
	      
		  //node in id Node1;

		  Node<3> & Node0 = rConditionGeom1[0]; // other node in condition [1]
		  Node<3> & Node2 = rConditionGeom2[1]; // other node in condition [2]


		  // std::cout<<"     Node0: "<<rConditionGeom1[0].Id()<<" Node 1: "<<rConditionGeom1[1].Id()<<std::endl;
		  // std::cout<<"     Node1: "<<rConditionGeom2[0].Id()<<" Node 2: "<<rConditionGeom2[1].Id()<<std::endl;

		  //segment condition [1]
		  S1[0] = rConditionGeom1[1].X() - rConditionGeom1[0].X();
		  S1[1] = rConditionGeom1[1].Y() - rConditionGeom1[0].Y();
	      
		  if(norm_2(S1)!=0)
		    S1/=norm_2(S1);

		  //segment condition [2]
		  S2[0] = rConditionGeom2[1].X() - rConditionGeom2[0].X();
		  S2[1] = rConditionGeom2[1].Y() - rConditionGeom2[0].Y();

		  if(norm_2(S2)!=0)
		    S2/=norm_2(S2);
		  
		  // std::cout<<"     S1 "<<S1<<std::endl;
		  // std::cout<<"     S2 "<<S2<<std::endl;


		  bool remove_S1 = false;
		  if(norm_2(S1)<size_for_side_normal)
		    remove_S1 = true;

		  bool remove_S2 = false;
		  if(norm_2(S2)<size_for_side_normal)
		    remove_S2 = true;

		  if(remove_S1 || remove_S2){
		    
		    node_shared_conditions[nodeId][i]->Set(TO_ERASE); //release condition [1]
		    node_shared_conditions[nodeId][j]->Set(TO_ERASE); //release condition [2]
		    in->Set(TO_ERASE);    //release Node1*

		    Condition::Pointer NewCond = node_shared_conditions[nodeId][i];
		    
		    Node0.Set(BLOCKED);
		    Node2.Set(BLOCKED);

		    //create new condition Node0-NodeB
		    Condition::NodesArrayType face;
		    face.reserve(2);

		    face.push_back(rConditionGeom1(0));
		    face.push_back(rConditionGeom2(1));
		
		    new_id = initial_cond_size + id;
		    //properties to be used in the generation
		    Condition::Pointer pcond       = NewCond->Clone(new_id, face);
		    // std::cout<<"     ID"<<id<<" 1s "<<pcond1->GetGeometry()[0].Id()<<" "<<pcond1->GetGeometry()[1].Id()<<std::endl;

		    pcond->Set(NEW_ENTITY);

		    //std::cout<<"     Condition INSERTED (Id: "<<new_id<<") ["<<rConditionGeom1[0].Id()<<", "<<rConditionGeom2[1].Id()<<"] "<<std::endl;

		    pcond->SetValue(NORMAL, NewCond->GetValue(NORMAL) );
		    pcond->SetValue(MASTER_NODES, NewCond->GetValue(MASTER_NODES) );

		    MeshDataTransferUtilities TransferUtilities;
		    TransferUtilities.TransferBoundaryData(pcond, NewCond, mrRemesh.TransferVariables); 

		    (mrModelPart.Conditions(mMeshId)).push_back(pcond);

		    RemovedNodes += 1;
		    id +=1;

		   
		  }
		  else{

		    double projection_sides   = inner_prod(S1,S2);
		    double projection_normals = inner_prod(N1,N2);
		    double relative_angle = 0;

		    if(projection_normals!=0)
		      relative_angle = projection_sides/projection_normals;
		  
		    if(relative_angle<=1 && relative_angle>=-1 )
		      condition_angle = (180.0/3.14159) * std::acos(relative_angle);
	    
		    if(inner_prod(S1,N2)<0) 
		      condition_angle *=(-1);

		    // std::cout<<"     projection_sides "<<projection_sides<<std::endl;
		    // std::cout<<"     projection_normals "<<projection_normals<<std::endl;
		    // std::cout<<"     relative_angle "<<relative_angle<<std::endl;
		    // std::cout<<"     condition_angle "<<condition_angle<<" critical_angle "<<critical_angle<<std::endl;
		  

		    if( condition_angle < -40 ){		    
		      // std::cout<<"     B NODE "<<in->Id()<<std::endl;
		      // std::cout<<"     projection_sides "<<projection_sides<<std::endl;
		      // std::cout<<"     projection_normals "<<projection_normals<<std::endl;
		      // std::cout<<"     relative_angle "<<relative_angle<<std::endl;
		      // std::cout<<"     condition_angle "<<condition_angle<<" critical_angle "<<critical_angle<<std::endl;
		      in->Set(VISITED);

		      Node0.Set(VISITED);
		      Node2.Set(VISITED);
		  
		    }

		    if(condition_angle<critical_angle){
		

		      //Path of neighbour conditions in 2D:   (NodeA) ---[0]--- (Node0) ---[1]--- (Node1*) ---[2]--- (Node2) ---[2]--- (Node2) ---[3]--- (NodeB)
		    
		      //realease positions:
		      node_shared_conditions[nodeId][i]->Set(TO_ERASE); //release condition [1]
		      node_shared_conditions[nodeId][j]->Set(TO_ERASE); //release condition [2]

		      in->Set(TO_ERASE);    //release Node1*
		      Node2.Set(TO_ERASE);  //release Node2

		      if( mEchoLevel > 0 ){
			std::cout<<"     Node Release/Modify  i "<<in->Id()<<std::endl;
			std::cout<<"     Node Release/Modify  j "<<Node2.Id()<<std::endl;
		      }

		      //set Node0 to a new position (between 0 and 2)
		      Node0.X() = 0.5 * ( Node0.X() + Node2.X() );
		      Node0.Y() = 0.5 * ( Node0.Y() + Node2.Y() );
		      Node0.Z() = 0.5 * ( Node0.Z() + Node2.Z() );

		      //assign data to dofs
		      unsigned int buffer_size = Node0.GetBufferSize();
		      unsigned int step_data_size = mrModelPart.GetNodalSolutionStepDataSize();

		      for(unsigned int step = 0; step<buffer_size; step++)
			{
			  //getting the data of the solution step
			  double* step_data = Node0.SolutionStepData().Data(step);

			  double* node0_data = Node0.SolutionStepData().Data(step);
			  double* node1_data = Node0.SolutionStepData().Data(step);

			  //copying this data in the position of the vector we are interested in
			  for(unsigned int j= 0; j<step_data_size; j++)
			    {
			      step_data[j] = 0.5*node0_data[j] + 0.5*node1_data[j];
			    }
			}
			
		      //recover the original position of the node
		      const array_1d<double,3>& disp = Node0.FastGetSolutionStepValue(DISPLACEMENT);
		      Node0.X0() = Node0.X() - disp[0];
		      Node0.Y0() = Node0.Y() - disp[1];
		      Node0.Z0() = Node0.Z() - disp[2];
	
		      //search shared condition of Node0 and Node A
		      if(node_shared_conditions[Node0.Id()][0]->Id() == Node0.Id()){
			i = 1;
		      }
		      else{
			i = 0;
		      }
		
		      Geometry< Node<3> >& rConditionGeom0 = node_shared_conditions[Node0.Id()][i]->GetGeometry();
		      Node<3> & NodeA = rConditionGeom0[0];

		      //search shared condition of Node2 and Node B
		      if(node_shared_conditions[Node2.Id()][0]->Id() == Node2.Id()){
			i = 0;
		      }
		      else{
			i = 1;
		      }
		
		      //New conditions profile in 2D:  (NodeA) ---[0]--- (Node0**) ---[3]--- (NodeB)   where (Node0**) is (Node0) in another position

		      Condition::Pointer NewCond = node_shared_conditions[Node2.Id()][i];
		      NewCond->Set(TO_ERASE);
		      Geometry< Node<3> >& rConditionGeom3 = NewCond->GetGeometry();
		      Node<3> & NodeB = rConditionGeom3[1];

		      NodeA.Set(ModelerUtilities::ENGAGED_NODES);
		      NodeB.Set(ModelerUtilities::ENGAGED_NODES);

		      Node0.Set(ModelerUtilities::ENGAGED_NODES);
		
		
		      //create new condition Node0-NodeB
		      Condition::NodesArrayType face;
		      face.reserve(2);

		      face.push_back(rConditionGeom1(0));
		      face.push_back(rConditionGeom3(1));
		
		      new_id = initial_cond_size + id;
		      //properties to be used in the generation
		      Condition::Pointer pcond       = NewCond->Clone(new_id, face);
		      // std::cout<<" ID"<<id<<" 1s "<<pcond1->GetGeometry()[0].Id()<<" "<<pcond1->GetGeometry()[1].Id()<<std::endl;

		      pcond->Set(NEW_ENTITY);

		      if( mEchoLevel > 0 ){
			std::cout<<"     Condition INSERTED (Id: "<<new_id<<") ["<<rConditionGeom1[0].Id()<<", "<<rConditionGeom3[1].Id()<<"] "<<std::endl;
		      }

		      rConditionGeom1[0].Set(TO_ERASE,false);  // do not release Node1
		      rConditionGeom3[1].Set(TO_ERASE,false);  // do not release Node2


		      pcond->SetValue(NORMAL, NewCond->GetValue(NORMAL) );
		      pcond->SetValue(MASTER_NODES, NewCond->GetValue(MASTER_NODES) );

		      MeshDataTransferUtilities TransferUtilities;
		      TransferUtilities.TransferBoundaryData(pcond, NewCond, mrRemesh.TransferVariables); 

		      (mrModelPart.Conditions(mMeshId)).push_back(pcond);

		      RemovedNodes += 1;
		      id +=1;
		
		    }
		  }
		}
	      }
	    }
	}

      for(ModelPart::NodesContainerType::const_iterator in = mrModelPart.NodesBegin(mMeshId); in != mrModelPart.NodesEnd(mMeshId); in++)
	{
	  in->Reset(BLOCKED);
	}
	      
    
      RemovedConditions = mrModelPart.Conditions(mMeshId).size() - RemovedConditions;

      if( mEchoLevel > 0 ){
	std::cout<<"     [ CONDITIONS ( removed : "<<RemovedConditions<<" ) ]"<<std::endl;
	std::cout<<"     [ NODES      ( removed : "<<RemovedNodes<<" ) ]"<<std::endl;
    
	std::cout<<"     Ending   Conditions : "<<mrModelPart.Conditions(mMeshId).size()<<"  (Removed nodes: "<< RemovedNodes<<" ) "<<std::endl;
	std::cout<<"     REMOVE NON CONVEX BOUNDARY ]; "<<std::endl;
      }

      if(RemovedNodes)
	return true;
      else
	return false;
      
      KRATOS_CATCH(" ")

    }


    /// Assignment operator.
    RemoveMeshNodesProcess& operator=(RemoveMeshNodesProcess const& rOther);


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
                                  RemoveMeshNodesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const RemoveMeshNodesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_REMOVE_MESH_NODES_PROCESS_H_INCLUDED  defined 


