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
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_processes/mesher_process.hpp"

#include "pfem_application_variables.h"

///VARIABLES used:
//Data:     NORMAL, MASTER_NODES, NEIGHBOUR_NODES, NEIGBOUR_ELEMENTS
//StepData: MEAN_ERROR
//Flags:    (checked) TO_ERASE, BOUNDARY, STRUCTURE, TO_SPLIT, NEW_ENTITY, BLOCKED
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

class RemoveFluidNodesMesherProcess
  : public MesherProcess
{
public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of Process
  KRATOS_CLASS_POINTER_DEFINITION( RemoveFluidNodesMesherProcess );

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
    RemoveFluidNodesMesherProcess(ModelPart& rModelPart,
			   MesherUtilities::MeshingParameters& rRemeshingParameters,
			   int EchoLevel)
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    {
      mEchoLevel = EchoLevel;
    }


  /// Destructor.
  virtual ~RemoveFluidNodesMesherProcess() {}


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
  void Execute() override
  {

    KRATOS_TRY

      if( mEchoLevel > 1 ){
	std::cout<<" [ REMOVE CLOSE NODES: "<<std::endl;
      }

    double RemovedConditions = mrModelPart.NumberOfConditions();
    // double NumberOfNodes = mrModelPart.NumberOfNodes();

    bool any_node_removed      = false;
    bool any_condition_removed = false;

    int error_nodes_removed    = 0;
    int inside_nodes_removed   = 0;
    int boundary_nodes_removed = 0;

    //if the remove_node switch is activated, we check if the nodes got too close
    if (mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_NODES))
      {
	if( mEchoLevel > 1 )
	  std::cout<<" REMOVE_NODES is TRUE "<<std::endl;
	bool any_node_removed_on_error = false;
	////////////////////////////////////////////////////////////
	if (mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_NODES_ON_ERROR))
	  {
	    if( mEchoLevel > 1 )
	      std::cout<<" REMOVE_NODES_ON_ERROR is TRUE "<<std::endl;

	    any_node_removed_on_error = RemoveNodesOnError(error_nodes_removed); //2D and 3D
	  }
	////////////////////////////////////////////////////////////
	if( mEchoLevel > 1 )
	  std::cout<<"error_nodes_removed :"<<error_nodes_removed<<std::endl;
	bool any_convex_condition_removed = false;

	bool any_node_removed_on_distance = false;
	////////////////////////////////////////////////////////////
	if (mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_NODES_ON_DISTANCE))
	  {
	    if( mEchoLevel > 1 )
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
	  this->CleanRemovedNodes(mrModelPart);

	if(any_condition_removed){
	  //Clean Conditions
	  ModelPart::ConditionsContainerType RemoveConditions;

	  //id = 0;
	  for(ModelPart::ConditionsContainerType::iterator ic = mrModelPart.ConditionsBegin(); ic!= mrModelPart.ConditionsEnd(); ++ic)
	    {

	      if(ic->IsNot(TO_ERASE)){
		//id+=1;
		RemoveConditions.push_back(*(ic.base()));
		//RemoveConditions.back().SetId(id);
	      }
	      else{
		if( mEchoLevel > 1 )
		  std::cout<<"   Condition RELEASED:"<<ic->Id()<<std::endl;
	      }
	    }

          mrModelPart.Conditions().swap(RemoveConditions);

	}


      }


    // number of removed nodes:
    // mrRemesh.Info->RemovedNodes = NumberOfNodes - mrModelPart.NumberOfNodes();
    mrRemesh.Info->RemovedNodes +=  error_nodes_removed + inside_nodes_removed + boundary_nodes_removed;
    int distance_remove =  inside_nodes_removed + boundary_nodes_removed;

    RemovedConditions -= mrModelPart.NumberOfConditions();

    if( mEchoLevel > 1 ){
      std::cout<<"  remove_mesh_nodes_process_for_fluids  [ CONDITIONS ( removed : "<<RemovedConditions<<" ) ]"<<std::endl;
      std::cout<<"   [ NODES      ( removed : "<<mrRemesh.Info->RemovedNodes<<" ) ]"<<std::endl;
      std::cout<<"   [ Error(removed: "<<error_nodes_removed<<"); Distance(removed: "<<distance_remove<<"; inside: "<<inside_nodes_removed<<"; boundary: "<<boundary_nodes_removed<<") ]"<<std::endl;


      //std::cout<<"   Nodes after  erasing : "<<mrModelPart.Nodes().size()<<std::endl;
      std::cout<<"   REMOVE CLOSE NODES ]; "<<std::endl;
    }

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
  std::string Info() const override
  {
    return "RemoveFluidNodesMesherProcess";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "RemoveFluidNodesMesherProcess";
  }

  /// Print object's data.
  void PrintData(std::ostream& rOStream) const override
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

  MesherUtilities::MeshingParameters& mrRemesh;

  MesherUtilities mMesherUtilities;

  int mEchoLevel;

  ///@}
  ///@name Un accessible methods
  ///@{

  //**************************************************************************
  //**************************************************************************

  void CleanRemovedNodes(ModelPart& rModelPart)
  {
    KRATOS_TRY

      //MESH 0 total domain mesh
      ModelPart::NodesContainerType temporal_nodes;
    temporal_nodes.reserve(rModelPart.Nodes().size());

    temporal_nodes.swap(rModelPart.Nodes());

    for(ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; ++i_node)
      {
	if( i_node->IsNot(TO_ERASE) ){
	  (rModelPart.Nodes()).push_back(*(i_node.base()));
	}
      }

    rModelPart.Nodes().Sort();


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


    MeshErrorDistribution.NodalErrorCalculation(mrModelPart,NodalError,nodes_ids,mrRemesh.Refine->GetErrorVariable());

    for(ModelPart::NodesContainerType::const_iterator in = mrModelPart.NodesBegin(); in != mrModelPart.NodesEnd(); ++in)
      {

	WeakPointerVector<Node<3> >& rN = in->GetValue(NEIGHBOUR_NODES);
	int erased_nodes =0;
	for(unsigned int i = 0; i < rN.size(); ++i)
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
	    for(WeakPointerVector< Element >::iterator ne = neighb_elems.begin(); ne!=neighb_elems.end(); ++ne)
	      {
		mean_node_radius+= mMesherUtilities.CalculateElementRadius(ne->GetGeometry()); //Triangle 2D, Tetrahedron 3D
		//mean_node_radius+= mMesherUtilities.CalculateTriangleRadius(ne->GetGeometry());
		//mean_node_radius+= mMesherUtilities.CalculateTetrahedronRadius(ne->GetGeometry());
	      }

	    mean_node_radius /= double(neighb_elems.size());

	    if(NodalError[nodes_ids[in->Id()]] < mrRemesh.Refine->ReferenceError && mean_node_radius < size_for_criterion_error)
	      {
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

      const unsigned int dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

    //***SIZES :::: parameters do define the tolerance in mesh size:
    double initialMeanRadius=0;
    initialMeanRadius=mrRemesh.Refine->InitialRadius;

    if( mEchoLevel > 1 )
      std::cout<<"initialMeanRadius "<<initialMeanRadius<<std::endl;
    double size_for_distance_inside       = 0.6 * initialMeanRadius;//compared to element radius
    double size_for_distance_boundary     = 0.6 * initialMeanRadius; //compared to element radius
    double size_for_wall_tip_contact_side = 0.15 * mrRemesh.Refine->CriticalSide;

    // if(dimension==3){
    //   size_for_distance_inside       = 0.5 * initialMeanRadius;//compared to element radius
    //   size_for_distance_boundary     = 0.5 * initialMeanRadius; //compared to element radius
    //   size_for_wall_tip_contact_side = 0.15 * mrRemesh.Refine->CriticalSide;
    // }

    bool any_node_removed = false;

    //bucket size definition:
    unsigned int bucket_size = 20;

    //create the list of the nodes to be check during the search
    std::vector<Node<3>::Pointer> list_of_nodes;
    list_of_nodes.reserve(mrModelPart.NumberOfNodes());
    for(ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin() ; i_node != mrModelPart.NodesEnd() ; ++i_node)
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


    for(ModelPart::NodesContainerType::const_iterator in = mrModelPart.NodesBegin(); in != mrModelPart.NodesEnd(); ++in)
      {
	if(in->Is(TO_ERASE)){
	  any_node_removed = true;
	}

	if( in->IsNot(NEW_ENTITY) &&  in->IsNot(INLET) )
	// if( in->IsNot(NEW_ENTITY) )
	  {
	    radius = size_for_distance_inside;

	    work_point[0]=in->X();
	    work_point[1]=in->Y();
	    work_point[2]=in->Z();
	    unsigned int freeSurfaceNeighNodes=0;
	    // unsigned int rigidNeighNodes=0;

	    if(in->Is(FREE_SURFACE)){// it must be more difficult to erase a free_surface node, otherwise, lot of volume is lost
	      radius = 0.5  * initialMeanRadius;//compared with element radius
	      // radius = 0.4  * initialMeanRadius;//compared with element radius
	      WeakPointerVector< Node < 3 > >& neighb_nodes = in->GetValue(NEIGHBOUR_NODES);
	      unsigned int countRigid=0;
	      for (WeakPointerVector< Node <3> >::iterator nn = neighb_nodes.begin();nn != neighb_nodes.end(); ++nn)
		{
		  if(nn->Is(RIGID) || nn->Is(SOLID)){
		    countRigid++;
		  }
		}
	      if(countRigid==neighb_nodes.size()){
		radius = 0.15  * initialMeanRadius;
	      }
	    }else{
	      WeakPointerVector< Node < 3 > >& neighb_nodes = in->GetValue(NEIGHBOUR_NODES);
	      for (WeakPointerVector< Node <3> >::iterator nn = neighb_nodes.begin();nn != neighb_nodes.end(); ++nn)
		{
		  if(nn->Is(FREE_SURFACE)){
		    freeSurfaceNeighNodes++;
		  }
		  // if(nn->Is(RIGID)){
		  //   rigidNeighNodes++;
		  // }
		}
	    }

	    if(in->Is(INLET)){
	      radius = 0.3  * initialMeanRadius;//compared with element radius
	    }
	    n_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, neighbours.begin(),neighbour_distances.begin(), num_neighbours);

	    if (n_points_in_radius>1)
	      {

		if (  in->IsNot(INLET) && in->IsNot(RIGID) && in->IsNot(SOLID) && in->IsNot(ISOLATED) )
		  {
		    if( mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_NODES_ON_DISTANCE) ){

                      // if the node is close to the surface, do not erase it, move to a mean (laplacian) position
		       if (in->IsNot(FREE_SURFACE) && in->IsNot(RIGID) && freeSurfaceNeighNodes==dimension){

		      	WeakPointerVector< Node < 3 > >& neighb_nodes = in->GetValue(NEIGHBOUR_NODES);
		      	array_1d<double,3> sumOfCoordinates=in->Coordinates();
                        array_1d<double,3> sumOfCurrentDisplacements=in->FastGetSolutionStepValue(DISPLACEMENT,0);
		      	array_1d<double,3> sumOfCurrentVelocities=in->FastGetSolutionStepValue(VELOCITY,0);
		      	array_1d<double,3> sumOfPreviousVelocities=in->FastGetSolutionStepValue(VELOCITY,1);
		      	double sumOfPressures=in->FastGetSolutionStepValue(PRESSURE,0);
		      	double counter=1.0;

		      	for (WeakPointerVector< Node <3> >::iterator nn = neighb_nodes.begin();nn != neighb_nodes.end(); ++nn)
		      	  {
		      	    counter+=1.0;
		      	    noalias(sumOfCoordinates)+=nn->Coordinates();
                            noalias(sumOfCurrentDisplacements)+=nn->FastGetSolutionStepValue(DISPLACEMENT,0);
		      	    noalias(sumOfCurrentVelocities)+=nn->FastGetSolutionStepValue(VELOCITY,0);
		      	    noalias(sumOfPreviousVelocities)+=nn->FastGetSolutionStepValue(VELOCITY,1);
		      	    sumOfPressures+=nn->FastGetSolutionStepValue(PRESSURE,0);
		      	  }
		      	in->X()  = sumOfCoordinates[0]/counter;
		      	in->Y()  = sumOfCoordinates[1]/counter;
		      	in->X0() = in->X()-sumOfCurrentDisplacements[0]/counter;
		      	in->Y0() = in->Y()-sumOfCurrentDisplacements[1]/counter;
		      	in->FastGetSolutionStepValue(DISPLACEMENT_X,0)=sumOfCurrentDisplacements[0]/counter;
		      	in->FastGetSolutionStepValue(DISPLACEMENT_Y,0)=sumOfCurrentDisplacements[1]/counter;
		      	in->FastGetSolutionStepValue(DISPLACEMENT_X,1)=sumOfCurrentDisplacements[0]/counter;
		      	in->FastGetSolutionStepValue(DISPLACEMENT_Y,1)=sumOfCurrentDisplacements[1]/counter;
		      	in->FastGetSolutionStepValue(VELOCITY_X,0)=sumOfCurrentVelocities[0]/counter;
		      	in->FastGetSolutionStepValue(VELOCITY_Y,0)=sumOfCurrentVelocities[1]/counter;
		      	in->FastGetSolutionStepValue(VELOCITY_X,1)=sumOfPreviousVelocities[0]/counter;
		      	in->FastGetSolutionStepValue(VELOCITY_Y,1)=sumOfPreviousVelocities[1]/counter;
		      	in->FastGetSolutionStepValue(PRESSURE,0)=sumOfPressures/counter;
		      	if(dimension==3){
		      	  in->Z() =sumOfCoordinates[2]/counter;
		      	  in->Z0() =in->Z()-sumOfCurrentDisplacements[2]/counter;
		      	  in->FastGetSolutionStepValue(DISPLACEMENT_Z,0)=sumOfCurrentDisplacements[2]/counter;
		      	  in->FastGetSolutionStepValue(DISPLACEMENT_Z,1)=sumOfCurrentDisplacements[2]/counter;
		      	  in->FastGetSolutionStepValue(VELOCITY_Z,0)=sumOfCurrentVelocities[2]/counter;
		      	  in->FastGetSolutionStepValue(VELOCITY_Z,1)=sumOfPreviousVelocities[2]/counter;
		      	}

		      } else{
			//look if we are already erasing any of the other nodes
			unsigned int erased_nodes = 0;
			for(std::vector<Node<3>::Pointer>::iterator nn=neighbours.begin(); nn!=neighbours.begin() + n_points_in_radius ; ++nn)
			  {
			    if( (*nn)->Is(TO_ERASE) )
			      erased_nodes += 1;
			  }

			if( erased_nodes < 1 ){ //we release the node if no other nodes neighbours are being erased
			  in->Set(TO_ERASE);
			  any_node_removed = true;
			  inside_nodes_removed++;
			  //distance_remove++;
			}
		      }
		    }

		  }

		else if ( in->IsNot(INLET) )  {

		  // std::cout<<"  Remove close boundary nodes: Candidate ["<<in->Id()<<"]"<<std::endl;
		  //here we loop over the neighbouring nodes and if there are nodes
		  //with BOUNDARY flag and closer than 0.2*nodal_h from our node, we remove the node we are considering
		  unsigned int k = 0;
		  unsigned int counter = 0;
		  for(std::vector<Node<3>::Pointer>::iterator nn=neighbours.begin(); nn!=neighbours.begin() + n_points_in_radius ; ++nn)
                  {

                    if ( (*nn)->Is(BOUNDARY) &&  neighbour_distances[k] < size_for_distance_boundary && neighbour_distances[k] > 0.0 )
                    {
                      if((*nn)->IsNot(TO_ERASE)){
                        counter += 1;
                      }
                    }

                    if ( (*nn)->Is(BOUNDARY) && neighbour_distances[k] < size_for_wall_tip_contact_side ) {
                      if ( (*nn)->IsNot(TO_ERASE)) {
                        counter += 1;
                      }
                    }

                    k++;
                  }

		  if(counter > 1 && in->IsNot(RIGID) && in->IsNot(SOLID) && in->IsNot(NEW_ENTITY)){ //Can be inserted in the boundary refine
		    in->Set(TO_ERASE);
		    if( mEchoLevel > 1 )
		      std::cout<<"     Removed Boundary Node ["<<in->Id()<<"] on Distance "<<std::endl;
		    any_node_removed = true;
		    boundary_nodes_removed++;
		    //distance_remove ++;
		  }
		  else if ( counter > 2 && in->IsNot(RIGID) && in->IsNot(SOLID) && in->IsNot(NEW_ENTITY) ) {
		    in->Set(TO_ERASE);
		    if( mEchoLevel > 1 )
		      std::cout << "     Removing a TIP POINT due to that criterion [" << in->Id() << "]" << std::endl;
		    any_node_removed = true;
		    boundary_nodes_removed++;
		  }

		}


	      }
	  }
      }

    unsigned int erased_nodes=0;
    for(ModelPart::ElementsContainerType::const_iterator ie = mrModelPart.ElementsBegin();
	ie != mrModelPart.ElementsEnd(); ++ie)
    {
      unsigned int rigidNodes=0;
      //coordinates
      for(unsigned int i=0; i<ie->GetGeometry().size(); ++i)
      {
        if((ie->GetGeometry()[i].Is(RIGID) && ie->GetGeometry()[i].IsNot(INLET)) || ie->GetGeometry()[i].Is(SOLID)){
          rigidNodes++;
        }
      }


      if(dimension==2){
        if(rigidNodes==2)
          EraseCriticalNodes2D(ie->GetGeometry(),erased_nodes,inside_nodes_removed);
      }else if(dimension==3){
        if(rigidNodes>1)
          EraseCriticalNodes3D(ie->GetGeometry(),erased_nodes,inside_nodes_removed);
      }


    }


    if(erased_nodes>0){
      if( mEchoLevel > 1 )
	std::cout<<"layer_nodes_removed "<<erased_nodes<<std::endl;
      any_node_removed = true;
    }
    //Build boundary after removing boundary nodes due distance criterion
    if( mEchoLevel > 1 ){
      std::cout<<"boundary_nodes_removed "<<boundary_nodes_removed<<std::endl;
      std::cout<<"inside_nodes_removed "<<inside_nodes_removed<<std::endl;
    }
    return any_node_removed;

    KRATOS_CATCH(" ")

      }


  void EraseCriticalNodes2D( Element::GeometryType& Element, unsigned int &erased_nodes,int& inside_nodes_removed)
  {

    KRATOS_TRY

      // std::cout<<"erased_nodes "<<erased_nodes<<std::endl;
    double safetyCoefficient2D=0.5;
    double elementVolume=Element.Area();

    // ////////  it erases nodes in very small elements /////////
    // double criticalVolume=0.1*mrRemesh.Refine->MeanVolume;
    // criticalVolume=0;
    // if(elementVolume<criticalVolume){
    //   for(unsigned int i=0; i<Element.size(); ++i)
    // 	{
    // 	  if(Element[i].IsNot(RIGID) && Element[i].IsNot(SOLID) && Element[i].IsNot(TO_ERASE)){
    // 	    Element[i].Set(TO_ERASE);
    // 	    if( mEchoLevel > 1 )
    // 	      std::cout<<"erase this layer node because it may be potentially dangerous and pass through the solid contour"<<std::endl;
    // 	    erased_nodes += 1;
    // 	    inside_nodes_removed++;
    // 	    break;
    // 	  }
    // 	}

    // }


    array_1d<double,3> Edges(3,0.0);
    array_1d<unsigned int,3> FirstEdgeNode(3,0);
    array_1d<unsigned int,3> SecondEdgeNode(3,0);
    double wallLength=0;
    // array_1d<double,3> CoorDifference(3,0.0);

    // ////////  to compute the length of the wall edge /////////
    // noalias(CoorDifference) = Element[1].Coordinates() - Element[0].Coordinates();
    array_1d<double,3> CoorDifference= Element[1].Coordinates() - Element[0].Coordinates();
    double SquaredLength = CoorDifference[0]*CoorDifference[0] + CoorDifference[1]*CoorDifference[1];
    Edges[0]=sqrt(SquaredLength);
    FirstEdgeNode[0]=0;
    SecondEdgeNode[0]=1;
    if(Element[0].Is(RIGID)  && Element[1].Is(RIGID) ){
      wallLength=Edges[0];
    }
    unsigned int counter=0;
    for (unsigned int i = 2; i < Element.size(); ++i){
      for(unsigned int j = 0; j < i; ++j)
	{
	  noalias(CoorDifference) = Element[i].Coordinates() - Element[j].Coordinates();
	  SquaredLength = CoorDifference[0]*CoorDifference[0] + CoorDifference[1]*CoorDifference[1];
	  counter+=1;
	  Edges[counter]=sqrt(SquaredLength);
	  FirstEdgeNode[counter]=j;
	  SecondEdgeNode[counter]=i;
	  if(Element[i].Is(RIGID) && Element[j].Is(RIGID) && Edges[counter]>wallLength ){
	    wallLength=Edges[counter];
	  }
	}

    }

    ////////  to compare the triangle height to wall edge length /////////
    for (unsigned int i = 0; i < Element.size(); ++i){
      if(Element[i].IsNot(RIGID) && Element[i].IsNot(TO_ERASE) && Element[i].IsNot(SOLID) && Element[i].IsNot(ISOLATED)){
    	double height=elementVolume*2.0/wallLength;

	//////it is evident when a freesurface particle in touch with wall is erased --> reduce the safety coeff
	if(Element[i].Is(FREE_SURFACE)){
	  WeakPointerVector< Node < 3 > >& neighb_nodes = Element[i].GetValue(NEIGHBOUR_NODES);
	  unsigned int countRigid=0;
	  unsigned int countFreeSurface=0;
	  for (WeakPointerVector< Node <3> >::iterator nn = neighb_nodes.begin();nn != neighb_nodes.end(); ++nn)
	      {
		if(nn->Is(RIGID) || nn->Is(SOLID)){
		  countRigid++;
		}
		if(nn->Is(FREE_SURFACE) && nn->IsNot(RIGID) && nn->IsNot(SOLID)){
		  countFreeSurface++;
		}
	      }
	  if((countRigid+countFreeSurface)==neighb_nodes.size() && countRigid>0){
	    safetyCoefficient2D=0.25;
	  }
	}

	////// if the node is very close to the wall is erased in any case
    	if(height<(0.5*safetyCoefficient2D*wallLength)){
	  Element[i].Set(TO_ERASE);
	  erased_nodes += 1;
	  inside_nodes_removed++;
	}

	// // if the node is near to the wall but not too close, if possible, it is not erased but just moved in the middle of its largest edge (not shared with a wall node)
	// else if(height<safetyCoefficient2D*wallLength){
	//   bool eraseNode=true;
	//   eraseNode=CheckForMovingLayerNodes(Element[i],wallLength);

	//     if(eraseNode==true){
	//       // std::cout<<"I will erase this node because too close to neighbour nodes "<<std::endl;
	//       std::cout<<"(distances:  "<<height<<" vs "<<wallLength<<")"<<std::endl;
	//     Element[i].Set(TO_ERASE);
	//     erased_nodes += 1;
	//     inside_nodes_removed++;
	//   }
	// }


      }

    }

    // ////////  to compare the non-wall length to wall edge length /////////
    // for (unsigned int i = 0; i < 3; ++i){
    //   if(((Element[FirstEdgeNode[i]].Is(RIGID) && Element[SecondEdgeNode[i]].IsNot(RIGID)) ||
    // 	  (Element[SecondEdgeNode[i]].Is(RIGID) && Element[FirstEdgeNode[i]].IsNot(RIGID))) &&
    // 	 Element[FirstEdgeNode[i]].IsNot(TO_ERASE) &&
    // 	 Element[SecondEdgeNode[i]].IsNot(TO_ERASE)&&
    // 	 Edges[i]<safetyCoefficient2D*wallLength){
    // 	if(Element[FirstEdgeNode[i]].IsNot(RIGID) && Element[FirstEdgeNode[i]].IsNot(SOLID) && Element[FirstEdgeNode[i]].IsNot(TO_ERASE)){
    // 	  Element[FirstEdgeNode[i]].Set(TO_ERASE);
    // 	  erased_nodes += 1;
    // 	  inside_nodes_removed++;
    // 	}else if(Element[SecondEdgeNode[i]].IsNot(RIGID) && Element[SecondEdgeNode[i]].IsNot(SOLID) && Element[SecondEdgeNode[i]].IsNot(TO_ERASE)){
    // 	  Element[SecondEdgeNode[i]].Set(TO_ERASE);
    // 	  erased_nodes += 1;
    // 	  inside_nodes_removed++;
    // 	}

    //   }

    // }
    KRATOS_CATCH( "" )


      }


  void EraseCriticalNodes3D( Element::GeometryType& Element, unsigned int &erased_nodes,int& inside_nodes_removed)
  {

    KRATOS_TRY

    double safetyCoefficient3D=0.6;
    // double safetyCoefficient3D=0.7;


    double elementVolume=Element.Volume();
    double criticalVolume=0.1*mrRemesh.Refine->MeanVolume;
    if(elementVolume<criticalVolume){
      for(unsigned int i=0; i<Element.size(); ++i)
	{
	  if(Element[i].IsNot(RIGID) && Element[i].IsNot(SOLID) && Element[i].IsNot(TO_ERASE)){
	    Element[i].Set(TO_ERASE);
	    if( mEchoLevel > 1)
	      std::cout<<"erase this layer node because it may be potentially dangerous and pass through the solid contour"<<std::endl;
	    erased_nodes += 1;
	    inside_nodes_removed++;
	    break;
	  }
	}

    }
    array_1d<double,6> Edges(6,0.0);
    array_1d<unsigned int,6> FirstEdgeNode(6,0);
    array_1d<unsigned int,6> SecondEdgeNode(6,0);
    double wallLength=0;
    double minimumLength=0;
    // array_1d<double,3> CoorDifference(3,0.0);

    // ////////  to compute the length of the wall edge /////////
    // CoorDifference = Element[1].Coordinates() - Element[0].Coordinates();
    array_1d<double,3> CoorDifference= Element[1].Coordinates() - Element[0].Coordinates();

    double SquaredLength = CoorDifference[0]*CoorDifference[0] +
      CoorDifference[1]*CoorDifference[1] +
      CoorDifference[2]*CoorDifference[2];
    Edges[0]=sqrt(SquaredLength);
    FirstEdgeNode[0]=0;
    SecondEdgeNode[0]=1;
    if(Element[0].Is(RIGID) && Element[1].Is(RIGID)){
      wallLength=Edges[0];
    }
    if((Element[0].Is(RIGID) && Element[1].IsNot(RIGID)) ||
       (Element[1].Is(RIGID) && Element[0].IsNot(RIGID)) ){
      minimumLength=Edges[0];
    }
    unsigned int counter=0;
    for (unsigned int i = 2; i < Element.size(); ++i){
      for(unsigned int j = 0; j < i; ++j)
	{
	  noalias(CoorDifference) = Element[i].Coordinates() - Element[j].Coordinates();
	  // CoorDifference = Element[i].Coordinates() - Element[j].Coordinates();
	  SquaredLength = CoorDifference[0]*CoorDifference[0] +
	    CoorDifference[1]*CoorDifference[1] +
	    CoorDifference[2]*CoorDifference[2];
	  counter+=1;
	  Edges[counter]=sqrt(SquaredLength);
	  FirstEdgeNode[counter]=j;
	  SecondEdgeNode[counter]=i;
	  if(Element[i].Is(RIGID) && Element[j].Is(RIGID) && wallLength==0 ){
	    wallLength=Edges[counter];
	  }
	  if(((Element[i].Is(RIGID) && Element[j].IsNot(RIGID)) ||
	     (Element[j].Is(RIGID) && Element[i].IsNot(RIGID)) ) &&
	     (Edges[counter]<minimumLength || minimumLength==0)){
	    minimumLength=Edges[counter];
	  }
	}

    }

    ////////  to avoid the elimanation of isolated free-surface-rigid elements /////////
    for (unsigned int i = 0; i < Element.size(); ++i){
      if(Element[i].IsNot(RIGID) && Element[i].IsNot(TO_ERASE) && Element[i].IsNot(SOLID) && Element[i].IsNot(ISOLATED)){
	//////it is evident when a freesurface particle in touch with wall is erased --> reduce the safety coeff
	if(Element[i].Is(FREE_SURFACE)){
	  WeakPointerVector< Node < 3 > >& neighb_nodes = Element[i].GetValue(NEIGHBOUR_NODES);
	  unsigned int countRigid=0;
	  unsigned int countFreeSurface=0;
	  for (WeakPointerVector< Node <3> >::iterator nn = neighb_nodes.begin();nn != neighb_nodes.end(); ++nn)
	    {
	      if(nn->Is(RIGID) || nn->Is(SOLID)){
		countRigid++;
	      }
	      if(nn->Is(FREE_SURFACE) && nn->IsNot(RIGID) && nn->IsNot(SOLID)){
		countFreeSurface++;
	      }
	    }
	  if((countRigid+countFreeSurface)==neighb_nodes.size() && countRigid>0){
	    safetyCoefficient3D=0.25;
	  }
	}
      }

    }

    // ////////// ////////// ////////// ////////// ////////// ////////
    // if(minimumLength<(0.5*safetyCoefficient3D*wallLength)){
    //   std::cout<<"(1. minimumLength:  "<<minimumLength<<" vs "<<wallLength<<")"<<std::endl;

    //   for (unsigned int i = 0; i < Element.size(); ++i){
    // 	if(Element[i].IsNot(RIGID) && Element[i].IsNot(TO_ERASE) && Element[i].IsNot(SOLID) && Element[i].IsNot(ISOLATED)){

    // 	  Element[i].Set(TO_ERASE);
    // 	  inside_nodes_removed++;
    // 	  erased_nodes += 1;
    // 	}
    //   }
    // }
    // else if(minimumLength<safetyCoefficient3D*wallLength){

    //   std::cout<<"(2. minimumLength:  "<<minimumLength<<" vs "<<wallLength<<")"<<std::endl;

    //   for (unsigned int i = 0; i < Element.size(); ++i){
    // 	if(Element[i].IsNot(RIGID) && Element[i].IsNot(TO_ERASE) && Element[i].IsNot(SOLID) && Element[i].IsNot(ISOLATED)){
    // 	  bool eraseNode=true;
    // 	  eraseNode=CheckForMovingLayerNodes(Element[i],wallLength);

    // 	  if(eraseNode==true){
    // 	    std::cout<<"I will erase this node because too close to neighbour nodes "<<std::endl;
    // 	    Element[i].Set(TO_ERASE);
    // 	    erased_nodes += 1;
    // 	    inside_nodes_removed++;
    // 	  }
    // 	}
    //   }
    // }
    // ////////// ////////// ////////// ////////// ////////// ////////





    // ////////  to compare the non-wall length to wall edge length /////////
    for (unsigned int i = 0; i < Edges.size(); ++i){
      if(((Element[FirstEdgeNode[i]].Is(RIGID) && Element[SecondEdgeNode[i]].IsNot(RIGID)) ||
    	  (Element[SecondEdgeNode[i]].Is(RIGID) && Element[FirstEdgeNode[i]].IsNot(RIGID))) &&
    	 Element[FirstEdgeNode[i]].IsNot(TO_ERASE) &&
    	 Element[SecondEdgeNode[i]].IsNot(TO_ERASE)&&
    	 Edges[i]<safetyCoefficient3D*wallLength){
    	if(Element[FirstEdgeNode[i]].IsNot(RIGID) && Element[FirstEdgeNode[i]].IsNot(SOLID) && Element[FirstEdgeNode[i]].IsNot(TO_ERASE) && Element[FirstEdgeNode[i]].IsNot(ISOLATED)){

    	  Element[FirstEdgeNode[i]].Set(TO_ERASE);
    	  inside_nodes_removed++;
    	  erased_nodes += 1;
    	}else if(Element[SecondEdgeNode[i]].IsNot(RIGID) && Element[SecondEdgeNode[i]].IsNot(SOLID) && Element[SecondEdgeNode[i]].IsNot(TO_ERASE) && Element[SecondEdgeNode[i]].IsNot(ISOLATED)){
    	  Element[SecondEdgeNode[i]].Set(TO_ERASE);
    	  inside_nodes_removed++;
    	  erased_nodes += 1;
    	}

      }

    }



    KRATOS_CATCH( "" )

      }


  bool CheckForMovingLayerNodes( Node < 3 >&  CheckedNode,const double wallLength)
  {
    KRATOS_TRY
      const unsigned int dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

    WeakPointerVector< Node < 3 > >& neighb_nodes = CheckedNode.GetValue(NEIGHBOUR_NODES);
    bool eraseNode=true;
    double maxSquaredDistance=0;
    WeakPointerVector< Node < 3 > >::iterator j = neighb_nodes.begin();
    for (WeakPointerVector< Node <3> >::iterator nn = neighb_nodes.begin();nn != neighb_nodes.end(); ++nn)
      {
	if(nn->IsNot(RIGID) && nn->IsNot(SOLID)){
	  // std::cout<<"neigh coordinates: "<<nn->X()<<" "<<nn->Y()<<std::endl;
	  array_1d<double,3>  CoorNeighDifference=CheckedNode.Coordinates()-nn->Coordinates();
	  double squaredDistance=CoorNeighDifference[0]*CoorNeighDifference[0]+CoorNeighDifference[1]*CoorNeighDifference[1];
	  if(dimension==3){
	    squaredDistance+=CoorNeighDifference[2]*CoorNeighDifference[2];
	  }
	  if(squaredDistance>maxSquaredDistance){
	    // std::cout<<"(distances:  "<<squaredDistance<<" vs "<<maxSquaredDistance<<")"<<std::endl;
	    maxSquaredDistance=squaredDistance;
	    j=nn;
	  }
	}
      }
    //I have looked for the biggest edge for moving there the layer node
    double maxNeighDistance=sqrt(maxSquaredDistance);
    if(maxNeighDistance>wallLength && wallLength>0){
      for (WeakPointerVector< Node<3> >::iterator nn = neighb_nodes.begin();nn != neighb_nodes.end(); ++nn)
	{
	  if(nn==j){

	    unsigned int idMaster = CheckedNode.GetId();
	    unsigned int idSlave = j->GetId();
	    InterpolateFromTwoNodes(idMaster,idMaster,idSlave);
	    std::vector<double> NewCoordinates(3);
	    NewCoordinates[0] = (CheckedNode.X()+j->X())*0.5;
	    NewCoordinates[1] = (CheckedNode.Y()+j->Y())*0.5;
	    CheckedNode.X() =NewCoordinates[0];
	    CheckedNode.Y() =NewCoordinates[1];
	    CheckedNode.X0() =NewCoordinates[0];
	    CheckedNode.Y0() =NewCoordinates[1];
	    CheckedNode.FastGetSolutionStepValue(DISPLACEMENT_X,0)=0;
	    CheckedNode.FastGetSolutionStepValue(DISPLACEMENT_Y,0)=0;
	    CheckedNode.FastGetSolutionStepValue(DISPLACEMENT_X,1)=0;
	    CheckedNode.FastGetSolutionStepValue(DISPLACEMENT_Y,1)=0;
	    if(dimension==3){
	      NewCoordinates[2] = (CheckedNode.Z()+j->Z())*0.5;
	      CheckedNode.Z() =NewCoordinates[2];
	      CheckedNode.Z0() =NewCoordinates[2];
	      CheckedNode.FastGetSolutionStepValue(DISPLACEMENT_Z,0)=0;
	      CheckedNode.FastGetSolutionStepValue(DISPLACEMENT_Z,1)=0;
	    }
	    // std::cout<<"new coordinates: "<<CheckedNode.X()<<" "<<CheckedNode.Y()<<std::endl;
	    // std::cout<<"(distances:  "<<maxNeighDistance<<" vs "<<wallLength<<")"<<std::endl;
	  }
	  eraseNode=false;
	}
    }
    return eraseNode;

    KRATOS_CATCH( "" )

      }



  void InterpolateFromTwoNodes(unsigned int idMaster,unsigned int idSlave1,unsigned int idSlave2)
  {

    KRATOS_TRY

      Node<3>::Pointer MasterNode = mrModelPart.pGetNode(idMaster);
      Node<3>::Pointer SlaveNode1 = mrModelPart.pGetNode(idSlave1);
      Node<3>::Pointer SlaveNode2 = mrModelPart.pGetNode(idSlave2);

      VariablesList& rVariablesList = mrModelPart.GetNodalSolutionStepVariablesList();

      unsigned int buffer_size = MasterNode->GetBufferSize();

    for(VariablesList::const_iterator i_variable =  rVariablesList.begin();  i_variable != rVariablesList.end() ; ++i_variable)
      {
	//std::cout<<" name "<<i_variable->Name()<<std::endl;
	//std::cout<<" type "<<typeid(*i_variable).name()<<std::endl;
	std::string variable_name = i_variable->Name();
	if(KratosComponents<Variable<double> >::Has(variable_name))
	  {
	    //std::cout<<"double"<<std::endl;
	    Variable<double> variable = KratosComponents<Variable<double> >::Get(variable_name);
	    for(unsigned int step = 0; step<buffer_size; ++step)
	      {
		//getting the data of the solution step
		double& node_data = MasterNode->FastGetSolutionStepValue(variable, step);

		double node0_data = SlaveNode1->FastGetSolutionStepValue(variable, step);
		double node1_data = SlaveNode2->FastGetSolutionStepValue(variable, step);

		node_data = (0.5*node0_data + 0.5*node1_data);

	      }
	  }
	else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name))
	  {
	    //std::cout<<"array1d"<<std::endl;
	    Variable<array_1d<double, 3> > variable = KratosComponents<Variable<array_1d<double, 3> > >::Get(variable_name);
	    for(unsigned int step = 0; step<buffer_size; ++step)
	      {
		//getting the data of the solution step
		array_1d<double, 3>& node_data = MasterNode->FastGetSolutionStepValue(variable, step);

		const array_1d<double, 3>& node0_data = SlaveNode1->FastGetSolutionStepValue(variable, step);
		const array_1d<double, 3>& node1_data = SlaveNode2->FastGetSolutionStepValue(variable, step);
		noalias(node_data) = (0.5*node0_data + 0.5*node1_data);
		// node_data = (0.5*node0_data + 0.5*node1_data);
	      }

	  }
	else if(KratosComponents<Variable<int > >::Has(variable_name))
	  {
	    //std::cout<<"int"<<std::endl;
	    //NO INTERPOLATION
	  }
	else if(KratosComponents<Variable<bool > >::Has(variable_name))
	  {
	    //std::cout<<"bool"<<std::endl;
	    //NO INTERPOLATION
	  }
	else if(KratosComponents<Variable<Matrix > >::Has(variable_name))
	  {
	    //std::cout<<"Matrix"<<std::endl;
	    Variable<Matrix > variable = KratosComponents<Variable<Matrix > >::Get(variable_name);
	    for(unsigned int step = 0; step<buffer_size; ++step)
	      {
		//getting the data of the solution step
		Matrix& node_data = MasterNode->FastGetSolutionStepValue(variable, step);

		Matrix& node0_data = SlaveNode1->FastGetSolutionStepValue(variable, step);
		Matrix& node1_data = SlaveNode2->FastGetSolutionStepValue(variable, step);

		if( node_data.size1() > 0 && node_data.size2() ){
		  if( node_data.size1() == node0_data.size1() && node_data.size2() == node0_data.size2() &&
		      node_data.size1() == node1_data.size1() && node_data.size2() == node1_data.size2() ) {
		    noalias(node_data) = (0.5*node0_data + 0.5*node1_data);
		    // node_data = (0.5*node0_data + 0.5*node1_data);
		  }
		}
	      }

	  }
	else if(KratosComponents<Variable<Vector > >::Has(variable_name))
	  {
	    //std::cout<<"Vector"<<std::endl;
	    Variable<Vector >variable = KratosComponents<Variable<Vector > >::Get(variable_name);
	    for(unsigned int step = 0; step<buffer_size; ++step)
	      {
		//getting the data of the solution step
		Vector& node_data = MasterNode->FastGetSolutionStepValue(variable, step);

		Vector& node0_data = SlaveNode1->FastGetSolutionStepValue(variable, step);
		Vector& node1_data = SlaveNode2->FastGetSolutionStepValue(variable, step);

		if( node_data.size() > 0 ){
		  if( node_data.size() == node0_data.size() &&
		      node_data.size() == node1_data.size()) {
		    noalias(node_data) = (0.5*node0_data + 0.5*node1_data);
		    // node_data = (0.5*node0_data + 0.5*node1_data);
		  }
		}
	      }
	  }

      }

    KRATOS_CATCH( "" )

      }




    /// Assignment operator.
    RemoveFluidNodesMesherProcess& operator=(RemoveFluidNodesMesherProcess const& rOther);


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
                                  RemoveFluidNodesMesherProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const RemoveFluidNodesMesherProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_REMOVE_MESH_NODES_FOR_FLUIDS_PROCESS_H_INCLUDED  defined
