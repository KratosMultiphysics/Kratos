//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_REMOVE_NODES_MESHER_PROCESS_H_INCLUDED )
#define  KRATOS_REMOVE_NODES_MESHER_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "includes/model_part.h"
#include "custom_utilities/mesh_error_calculation_utilities.hpp"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_processes/mesher_process.hpp"

///VARIABLES used:
//Data:     NORMAL, MASTER_NODES, NEIGHBOUR_NODES, NEIGBOUR_ELEMENTS
//StepData: MEAN_ERROR, CONTACT_FORCE
//Flags:    (checked) TO_ERASE, BOUNDARY,  CONTACT, NEW_ENTITY, BLOCKED
//          (set)     TO_ERASE(conditions,nodes)(set), NEW_ENTITY(conditions,nodes)(set), BLOCKED(nodes)(set), INSIDE(nodes)(set)
//          (modified)
//          (reset)
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


/**NOTE: as kratos uses shared pointers, the main list of nodes can be changed.
   the element contains the geometry which has pointers to nodes, so that nodes are kept in memory and their information can be recovered.
   that means also that if the element information is needed it is kept until the end when new elements are created and replace the old ones.
*/


class RemoveNodesMesherProcess
    : public MesherProcess
{
 public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of Process
  KRATOS_CLASS_POINTER_DEFINITION( RemoveNodesMesherProcess );

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
  RemoveNodesMesherProcess(ModelPart& rModelPart,
			   MesherUtilities::MeshingParameters& rRemeshingParameters,
			   int EchoLevel)
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
  {
    mEchoLevel = EchoLevel;
  }


  /// Destructor.
  virtual ~RemoveNodesMesherProcess() {}


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

        if( mEchoLevel > 0 ){
          std::cout<<" [ REMOVE CLOSE NODES: "<<std::endl;
          //std::cout<<"   Nodes before erasing : "<<mrModelPart.Nodes().size()<<std::endl;
        }

    if( mrModelPart.Name() != mrRemesh.SubModelPartName )
      std::cout<<" ModelPart Supplied do not corresponds to the Meshing Domain: ("<<mrModelPart.Name()<<" != "<<mrRemesh.SubModelPartName<<")"<<std::endl;


    mrRemesh.Refine->Info.BodyNodesRemoved.Initialize();
    mrRemesh.Refine->Info.BoundaryNodesRemoved.Initialize();

    double RemovedConditions = mrModelPart.NumberOfConditions();
    double NumberOfNodes     = mrModelPart.NumberOfNodes();

    bool any_node_removed      = false;
    bool any_condition_removed = false;

    //if the remove_node switch is activated, we check if the nodes got too close
    if (mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_NODES) || (mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_BOUNDARY_NODES)) )
    {
      bool any_node_removed_on_error = false;
      ////////////////////////////////////////////////////////////
      if (mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_NODES_ON_ERROR) && mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_NODES) )
      {
        any_node_removed_on_error = this->RemoveNodesOnError(mrModelPart, mrRemesh.Refine->Info.BodyNodesRemoved.on_error); //2D and 3D
      }
      ////////////////////////////////////////////////////////////

      bool any_convex_condition_removed = false;
      ////////////////////////////////////////////////////////////
      if (mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_BOUNDARY_NODES))
      {
        any_convex_condition_removed = this->RemoveNonConvexBoundary(mrModelPart); //2D only
      }
      ////////////////////////////////////////////////////////////


      bool any_node_removed_on_distance = false;
      ////////////////////////////////////////////////////////////
      if (mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_NODES_ON_DISTANCE) || mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_BOUNDARY_NODES_ON_DISTANCE))
      {
        any_node_removed_on_distance = this->RemoveNodesOnDistance(mrModelPart, mrRemesh.Refine->Info.BodyNodesRemoved.on_distance, mrRemesh.Refine->Info.BoundaryNodesRemoved.on_distance, any_condition_removed); //2D only (RebuildBoundary is only 2D)
      }
      // REMOVE ON DISTANCE
      ////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////
      // REMOVE CONTACT NODES (and boundary near the contact)
      if ( mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_BOUNDARY_NODES_ON_DISTANCE) )
      {
        bool any_node_removed_on_contact = this->RemoveNodesOnContact( mrModelPart, mrRemesh.Refine->Info.BodyNodesRemoved.on_distance, mrRemesh.Refine->Info.BoundaryNodesRemoved.on_distance, any_condition_removed);
        if ( any_node_removed_on_contact || any_node_removed_on_distance)
          any_node_removed_on_distance = true;
      }
      ////////////////////////////////////////////////////////////

      if(any_node_removed_on_error || any_node_removed_on_distance)
        any_node_removed = true;

      if(any_convex_condition_removed || any_condition_removed)
        any_condition_removed = true;


      if(any_node_removed)
        this->CleanRemovedNodes(mrModelPart);

      if(any_condition_removed){
        //Clean Conditions
        ModelPart::ConditionsContainerType PreservedConditions;

        //id = 0;
        for(ModelPart::ConditionsContainerType::iterator ic = mrModelPart.ConditionsBegin(); ic!= mrModelPart.ConditionsEnd(); ++ic)
        {

          if(ic->IsNot(TO_ERASE)){
            //id+=1;
            PreservedConditions.push_back(*(ic.base()));
            //PreservedConditions.back().SetId(id);
          }
          else{
            std::cout<<"   Condition RELEASED:"<<ic->Id()<<std::endl;
          }
        }

        mrModelPart.Conditions().swap(PreservedConditions);

      }


    }


    // number of removed nodes:
    mrRemesh.Info->RemovedNodes = NumberOfNodes - mrModelPart.NumberOfNodes();
    RemovedConditions -= mrModelPart.NumberOfConditions();

    if( mEchoLevel > 0 ){
      std::cout<<"   [ CONDITIONS ( removed : "<<RemovedConditions<<" ) ]"<<std::endl;
      std::cout<<"   [ NODES      ( removed : "<<mrRemesh.Info->RemovedNodes<<" ) ]"<<std::endl;

      if( this->mEchoLevel >=1 ){
        mrRemesh.Refine->Info.BodyNodesRemoved.EchoStats();
        mrRemesh.Refine->Info.BoundaryNodesRemoved.EchoStats();
      }

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
    return "RemoveNodesMesherProcess";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "RemoveNodesMesherProcess";
  }

  /// Print object's data.
  void PrintData(std::ostream& rOStream) const override
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

  MesherUtilities::MeshingParameters& mrRemesh;

  MesherUtilities mMesherUtilities;

  int mEchoLevel;

  ///@}
  ///@name Protected Operators
  ///@{
  ///@}
  ///@name Protected Operations
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
      else{
        if( i_node->Is(BOUNDARY) )
          if( mEchoLevel > 0 )
            std::cout<<"   BOUNDARY NODE RELEASED "<<i_node->Id()<<std::endl;
      }
    }

    rModelPart.Nodes().Sort();


    KRATOS_CATCH( "" )
  }


  //**************************************************************************
  //**************************************************************************

  virtual bool RemoveNodesOnDistance(ModelPart& rModelPart, unsigned int& inside_nodes_removed, unsigned int& boundary_nodes_removed, bool& any_condition_removed)
  {
    KRATOS_TRY


    //***SIZES :::: parameters do define the tolerance in mesh size:
    double size_for_distance_inside       = 1.0 * mrRemesh.Refine->CriticalRadius; //compared with element radius
    double size_for_distance_boundary     = 1.5 * size_for_distance_inside; //compared with element radius


    bool any_node_removed = false;

    //bucket size definition:
    unsigned int bucket_size = 20;

    //create the list of the nodes to be check during the search
    std::vector<Node<3>::Pointer> list_of_nodes;
    list_of_nodes.reserve(rModelPart.NumberOfNodes());
    for(ModelPart::NodesContainerType::iterator i_node = rModelPart.NodesBegin() ; i_node != rModelPart.NodesEnd() ; ++i_node)
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


    for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(); in != rModelPart.NodesEnd(); ++in)
    {
      bool on_contact_tip = false;
      bool contact_active = false;

      if( in->SolutionStepsDataHas(CONTACT_FORCE) ){
        array_1d<double, 3 > & ContactForceNormal  = in->FastGetSolutionStepValue(CONTACT_FORCE);
        if(norm_2(ContactForceNormal)>0)
          contact_active = true;
      }

      if(contact_active || in->Is(CONTACT) )
        on_contact_tip = true;

      if( in->IsNot(NEW_ENTITY) && in->IsNot(INLET) && in->IsNot(RIGID) && in->IsNot(TO_ERASE) )
      {
        radius = size_for_distance_inside;

        work_point[0]=in->X();
        work_point[1]=in->Y();
        work_point[2]=in->Z();

        n_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, neighbours.begin(),neighbour_distances.begin(), num_neighbours);

        if (n_points_in_radius>1)
        {
          //std::cout<<"     Points in Radius "<< n_points_in_radius<<" radius "<<radius<<std::endl;
          if( in->IsNot(BOUNDARY) )
          {
            if( mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_NODES_ON_DISTANCE) ){

              if( !this->CheckEngagedNode((*in),neighbours,neighbour_distances,n_points_in_radius) ){ //we release the node if no other nodes neighbours are being erased
                in->Set(TO_ERASE);
                any_node_removed = true;
                inside_nodes_removed++;
                //distance_remove++;
                //std::cout<<"     Distance Criterion Node ["<<in->Id()<<"] TO_ERASE "<<std::endl;
              }

            }

          }
          else if ( (mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_BOUNDARY_NODES) && mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_BOUNDARY_NODES_ON_DISTANCE)) ) //boundary nodes will be removed if they get REALLY close to another boundary node (0.2(=extra_factor) * h_factor)
          {

            //std::cout<<"  Remove close boundary nodes: Candidate ["<<in->Id()<<"]"<<std::endl;

            //here we loop over the neighbouring nodes and if there are nodes
            //with BOUNDARY flag and closer than 0.2*nodal_h from our node, we remove the node we are considering
            unsigned int k = 0;
            unsigned int counter = 0;
            for(std::vector<Node<3>::Pointer>::iterator nn=neighbours.begin(); nn!=neighbours.begin() + n_points_in_radius ; ++nn)
            {
              bool nn_on_contact_tip = false;
              bool contact_active = false;

              if( (*nn)->SolutionStepsDataHas(CONTACT_FORCE) ){
                array_1d<double, 3 > & ContactForceNormal  = (*nn)->FastGetSolutionStepValue(CONTACT_FORCE);
                if(norm_2(ContactForceNormal)>0)
                  contact_active = true;
              }

              if(contact_active || (*nn)->Is(CONTACT) )
                nn_on_contact_tip = true;

              //std::cout<<" radius * extra_factor "<<(extra_factor*radius)<<" >? "<<neighbour_distances[k]<<std::endl;
              if ( (*nn)->Is(BOUNDARY) && !nn_on_contact_tip && neighbour_distances[k] < size_for_distance_boundary && neighbour_distances[k] > 0.0 )
              {
                //KRATOS_WATCH( neighbours_distances[k] )
                if((*nn)->IsNot(TO_ERASE)){
                  counter += 1;
                }
              }


              k++;
            }

            if(counter > 1 && !on_contact_tip ){ //Can be inserted in the boundary refine
              in->Set(TO_ERASE);
              //std::cout<<"     Removed Boundary Node ["<<in->Id()<<"] on Distance "<<std::endl;
              any_node_removed = true;
              boundary_nodes_removed++;
              //distance_remove ++;
            }

          }

          //}

        }

      }
    }

    //Build boundary after removing boundary nodes due distance criterion
    if(boundary_nodes_removed){
      any_condition_removed = RebuildBoundary(rModelPart);
    }
    //Build boundary after removing boundary nodes due distance criterion


    return any_node_removed;

    KRATOS_CATCH(" ")

  }

  //**************************************************************************
  //**************************************************************************

  virtual bool CheckEngagedNode(Node<3>& rNode, std::vector<Node<3>::Pointer>& rNeighbours, std::vector<double>& rNeighbourDistances, unsigned int& rn_points_in_radius)
  {
    KRATOS_TRY

    //look if we are already erasing any of the other nodes
    bool engaged_node = false;

    //trying to remove the node in the shorter distance anyway
    double min_distance = std::numeric_limits<double>::max();
    unsigned int counter = 0;
    unsigned int closest_node = 0;

    for(std::vector<double>::iterator nd=rNeighbourDistances.begin(); nd!=rNeighbourDistances.begin()+rn_points_in_radius; ++nd)
    {
      if( (*nd) < min_distance && (*nd) > 0.0 ){
        min_distance = (*nd);
        closest_node = counter;
      }
      //std::cout<<" Distances "<<counter<<" "<< (*nd) <<std::endl;
      ++counter;
    }

    if( rNeighbours[closest_node]->Is(TO_ERASE) && rNeighbours[closest_node]->IsNot(BOUNDARY) ){
      engaged_node = true;
    }
    else{
      counter = 0;
      unsigned int erased_node = 0;
      for(std::vector<Node<3>::Pointer>::iterator nn=rNeighbours.begin(); nn!=rNeighbours.begin()+rn_points_in_radius; ++nn)
      {
        if( (*nn)->IsNot(BOUNDARY) && (*nn)->Is(TO_ERASE) )
        {
          erased_node = counter;
          engaged_node = true;
          break;
        }
        ++counter;
      }
      if( engaged_node ){
        // if the distance is 5 times smaller remove anyway
        if( rNeighbourDistances[closest_node] * 5 <  rNeighbourDistances[erased_node] ){
          engaged_node = false;
          //std::cout<<"     Distance Criterion Node ["<<in->Id()<<"] TO_ERASE "<<closest_node<<" "<<erased_node<<std::endl;
        }
      }
    }

    return engaged_node;

    KRATOS_CATCH(" ")
  }

  //**************************************************************************
  //**************************************************************************

  virtual bool RemoveNodesOnError(ModelPart& rModelPart, unsigned int& error_removed_nodes)
  {
    KRATOS_TRY

    //***SIZES :::: parameters do define the tolerance in mesh size:
    double size_for_criterion_error   = 2.0 * mrRemesh.Refine->CriticalRadius; //compared with mean node radius

    bool any_node_removed = false;

    MeshErrorCalculationUtilities MeshErrorDistribution;
    MeshErrorDistribution.SetEchoLevel(mEchoLevel);

    std::vector<double> NodalError(rModelPart.NumberOfNodes()+1);

    unsigned int number_of_nodes = 0;

    if(mrRemesh.InputInitializedFlag)
      number_of_nodes = mrRemesh.NodeMaxId+1;
    else
      number_of_nodes = MesherUtilities::GetMaxNodeId(rModelPart)+1;

    std::vector<int> nodes_ids(number_of_nodes);


    MeshErrorDistribution.NodalErrorCalculation(rModelPart,NodalError,nodes_ids,mrRemesh.Refine->GetErrorVariable());

    for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(); in != rModelPart.NodesEnd(); ++in)
    {

      WeakPointerVector<Node<3> >& rN = in->GetValue(NEIGHBOUR_NODES);
      int erased_nodes =0;
      for(unsigned int i = 0; i < rN.size(); ++i)
      {
        if(rN[i].Is(TO_ERASE))
          erased_nodes += 1;
      }


      if( in->IsNot(BOUNDARY) && erased_nodes < 1 )
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

  virtual bool RemoveNodesOnContact(ModelPart& rModelPart, unsigned int& inside_nodes_removed, unsigned int& boundary_nodes_removed, bool& any_condition_removed)
  {
    KRATOS_TRY

    // sizes
    double size_for_distance_inside       = 1.0  * mrRemesh.Refine->CriticalRadius;
    double size_for_distance_boundary     = 1.5  * size_for_distance_inside;

    double size_for_wall_tip_contact_side = (0.5 * mrRemesh.Refine->CriticalSide); // the distance to refine
    size_for_wall_tip_contact_side *= 0.5; // half the distance
    size_for_wall_tip_contact_side =  size_for_wall_tip_contact_side * size_for_wall_tip_contact_side;  // the return of  search in radius is distance^2

    double SF = 0.80;
    size_for_wall_tip_contact_side *= SF*SF;


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
      bool on_contact_tip = false;
      array_1d<double, 3 > & ContactForceNormal  = in->FastGetSolutionStepValue(CONTACT_FORCE);

      if(norm_2(ContactForceNormal)>0 || in->Is(CONTACT) )
        on_contact_tip = true;

      bool on_contact_tip_strict = false;
      if (norm_2(ContactForceNormal) > 0)
        on_contact_tip_strict = true;

      if( in->IsNot(NEW_ENTITY) )
      {
        radius = size_for_distance_inside;

        work_point[0]=in->X();
        work_point[1]=in->Y();
        work_point[2]=in->Z();

        n_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, neighbours.begin(),neighbour_distances.begin(), num_neighbours);

        if (n_points_in_radius>1)
        {

          if ( in->IsNot(BOUNDARY) )
          {

            if( mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_NODES_ON_DISTANCE) ){
              //look if we are already erasing any of the other nodes
              unsigned int contact_nodes = 0;
              unsigned int erased_nodes = 0;
              unsigned int near_to_contact_nodes = 0;
              unsigned int kk = 0;
              for(std::vector<Node<3>::Pointer>::iterator nn=neighbours.begin(); nn!=neighbours.begin() + n_points_in_radius ; ++nn)
              {
                if( (*nn)->Is(BOUNDARY) && (*nn)->Is(CONTACT) )
                  contact_nodes += 1;

                if( (*nn)->Is(TO_ERASE) )
                  erased_nodes += 1;

                // to remove a node that is very close to a contact node (two times the safety factor)
                if ( (*nn)->Is(BOUNDARY) && (*nn)->Is(CONTACT) && (neighbour_distances[kk] < SF*size_for_wall_tip_contact_side) )
                  near_to_contact_nodes += 1;

                kk++;

              } // end for neighbours

              if ( erased_nodes < 1 && contact_nodes < 2 && near_to_contact_nodes == 1) // we release node if it is very very near a contact
              {
                // to remove an interior node to is to close to a contacting node
                in->Set(TO_ERASE);
                any_node_removed = true;
                inside_nodes_removed++;
                std::cout <<"   RemovingC0, an interior node very very near to a (possibly) contacting node " << in->Id() << std::endl;
                std::cout <<"      X: " << in->X() << " Y: " << in->Y() << std::endl;
              } // end else if

            }

          }
          else if ( (mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_BOUNDARY_NODES) && mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_BOUNDARY_NODES_ON_DISTANCE)) && (in)->IsNot(TO_ERASE))
            //boundary nodes will be removed if they get REALLY close to another boundary node (0.2(=extra_factor) * h_factor)
          {

            //std::cout<<"  Remove close boundary nodes: Candidate ["<<in->Id()<<"]"<<std::endl;

            //here we loop over the neighbouring nodes and if there are nodes
            //with BOUNDARY flag and closer than 0.2*nodal_h from our node, we remove the node we are considering
            unsigned int k = 0;
            unsigned int counterC2 = 0, counterC3 = 0;
            for(std::vector<Node<3>::Pointer>::iterator nn=neighbours.begin(); nn!=neighbours.begin() + n_points_in_radius ; ++nn)
            {


              bool nn_on_contact_tip = false;
              array_1d<double, 3 > & ContactForceNormal  = (*nn)->FastGetSolutionStepValue(CONTACT_FORCE);

              if(norm_2(ContactForceNormal)>0 || (*nn)->Is(CONTACT) )
                nn_on_contact_tip = true;

              bool nn_on_contact_tip_strict = false;
              if (norm_2(ContactForceNormal)>0)
                nn_on_contact_tip_strict = true;

              if ( (*nn)->Is(BOUNDARY)  && neighbour_distances[k] > 0.0 && (*nn)->IsNot(TO_ERASE) )
              {
                //KRATOS_WATCH( neighbours_distances[k] )
                if ( neighbour_distances[k] < size_for_wall_tip_contact_side ) {
                  if ( nn_on_contact_tip_strict && (*nn)->IsNot(NEW_ENTITY) ) {
                    counterC2 += 1;
                  }
                  if ( nn_on_contact_tip && nn_on_contact_tip_strict && neighbour_distances[k] < SF*size_for_distance_boundary ) {
                    counterC3 += 1;
                  }
                }
              }

              k++;
            }  // end for each neighbour

            if ( counterC2 > 1 && in->IsNot(NEW_ENTITY) && on_contact_tip_strict) {
              in->Set(TO_ERASE);
              std::cout << "     RemovingC2: three contacting nodes where close, removing the middle one [" << in->Id() << "]" << std::endl;
              any_node_removed = true;
              boundary_nodes_removed++;
              std::cout << "      X: " << in->X() << " Y: " << in->Y() << std::endl;
            }
            else if ( counterC3 > 0 && in->IsNot(NEW_ENTITY) && on_contact_tip && !on_contact_tip_strict ) {

              in->Set(TO_ERASE);
              any_node_removed = true;
              boundary_nodes_removed++;
              std::cout << "    RemovingC3: a non_contacting_node was to close to a contacting. removing the non_contacting " << in->Id() << std::endl;
              std::cout << "      X: " << in->X() << " Y: " << in->Y() << std::endl;

            }
          }
        }

      }
    }

    // New loop to see if just two contacting nodes are very near. (it has to be done after the others to not to remove a pair)

    if ( mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_BOUNDARY_NODES) &&  mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_BOUNDARY_NODES_ON_DISTANCE) )
    {
      for (ModelPart::NodesContainerType::const_iterator in = mrModelPart.NodesBegin(); in != mrModelPart.NodesEnd(); ++in)
      {

        bool on_contact_tip = false;
        array_1d<double, 3 > & ContactForceNormal  = in->FastGetSolutionStepValue(CONTACT_FORCE);

        if(norm_2(ContactForceNormal)>0 || in->Is(CONTACT) )
          on_contact_tip = true;

        bool on_contact_tip_strict = false;
        if (norm_2(ContactForceNormal) > 0)
          on_contact_tip_strict = true;

        if ( in->IsNot(NEW_ENTITY) &&  in->IsNot(TO_ERASE) && in->Is(BOUNDARY) && on_contact_tip )
        {
          radius = size_for_distance_inside;

          work_point[0]=in->X();
          work_point[1]=in->Y();
          work_point[2]=in->Z();

          n_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, neighbours.begin(),neighbour_distances.begin(), num_neighbours);

          if (n_points_in_radius>1)
          {
            unsigned int k = 0;
            unsigned int counterC4 = 0;
            for (std::vector<Node<3>::Pointer>::iterator  nn = neighbours.begin(); nn!=neighbours.begin() + n_points_in_radius; ++nn)
            {
              bool nn_on_contact_tip = false;
              array_1d<double, 3 > & ContactForceNormal  = (*nn)->FastGetSolutionStepValue(CONTACT_FORCE);

              // alternative definition
              if(norm_2(ContactForceNormal)>0 || (*nn)->Is(CONTACT) )
                nn_on_contact_tip = true;
              bool nn_on_contact_tip_strict = false;
              if (norm_2(ContactForceNormal) > 0)
                nn_on_contact_tip_strict = true;

              if ( (*nn)->IsNot(NEW_ENTITY) && (*nn)->IsNot(TO_ERASE) && (*nn)->Is(BOUNDARY) && neighbour_distances[k] > 0.0)
              {
                if ( nn_on_contact_tip && nn_on_contact_tip_strict && neighbour_distances[k] < SF*SF*size_for_wall_tip_contact_side)
                {
                  counterC4 += 1;
                  std::cout << " THIS IS THE CONTRARY NODE: " << (*nn)->X() << " " << (*nn)->Y() << std::endl;
                  std::cout << " THIS IS THE CONTRARY FORCE: " << ContactForceNormal << std::endl;
                  std::cout << " module " << norm_2(ContactForceNormal) << std::endl;
                }
              } // first if for C4 Condi
              k++;
            } // end for all Neighbours

            if (counterC4 > 0 && in->IsNot(NEW_ENTITY) && on_contact_tip && on_contact_tip_strict )
            {
              in->Set(TO_ERASE);
              any_node_removed = true;
              boundary_nodes_removed++;
              std::cout << "    RemovingC4: two contacting nodes are very very near, removing one " << in->Id() << std::endl;
              std::cout << "      X: " << in->X() << " Y: " << in->Y() << std::endl;
            }
          }
        }

      }
    }


    //Build boundary after removing boundary nodes due distance criterion
    if(boundary_nodes_removed){
      any_condition_removed = RebuildBoundary(rModelPart);
    }
    //Build boundary after removing boundary nodes due distance criterion


    return any_node_removed;
    KRATOS_CATCH("")
  }

  //**************************************************************************
  //**************************************************************************

  bool RebuildBoundary(ModelPart& rModelPart)
  {
    KRATOS_TRY

        bool any_condition_removed = false;

    unsigned int number_of_nodes = 0;
    if( !rModelPart.IsSubModelPart() )
      number_of_nodes = rModelPart.NumberOfNodes()+1;
    else
      number_of_nodes = rModelPart.GetParentModelPart()->NumberOfNodes()+1;

    std::vector<std::vector<Condition::Pointer> > node_shared_conditions(number_of_nodes); //all domain nodes

    for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(); ic!= rModelPart.ConditionsEnd(); ++ic)
    {
      if(ic->IsNot(NEW_ENTITY) && ic->IsNot(TO_ERASE)){
        Geometry< Node<3> >& rConditionGeom = ic->GetGeometry();
        for(unsigned int i=0; i<rConditionGeom.size(); ++i){
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

    unsigned int initial_cond_size = 0;
    if( !rModelPart.IsSubModelPart() )
      initial_cond_size = rModelPart.NumberOfConditions()+1; //total model part conditions size
    else
      initial_cond_size = rModelPart.GetParentModelPart()->NumberOfConditions()+1;

    unsigned int id = 1;
    unsigned int new_id = 0;

    for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(); in != rModelPart.NodesEnd(); ++in)
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
            if ( pcond->Is(TO_ERASE) ) {
              pcond->Reset(TO_ERASE); // due to the new cloning
            }


            //std::cout<<"     Condition INSERTED (Id: "<<new_id<<") ["<<rConditionGeom1[0].Id()<<", "<<rConditionGeom2[1].Id()<<"] "<<std::endl;

            pcond->SetValue(NORMAL, NewCond->GetValue(NORMAL) );
            pcond->SetValue(MASTER_NODES, NewCond->GetValue(MASTER_NODES) );

            MeshDataTransferUtilities TransferUtilities;
            TransferUtilities.TransferBoundaryData(pcond, NewCond, *(mrRemesh.Transfer));

            (rModelPart.Conditions()).push_back(pcond);

            id +=1;
          }

        } else {
          in->Set(TO_ERASE, false);
          std::cout << "FINALLY NOT Removing " << in->Id() << std::endl;
        }
      }

    }

    // for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(); in != rModelPart.NodesEnd(); ++in)
    // 	{
    // 	  in->Reset(BLOCKED);
    // 	}

    return any_condition_removed;

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

  ///@name Private Static Member Variables
  ///@{
  ///@}
  ///@name Private Static Member Variables
  ///@{
  ///@}
  ///@name Private Operators
  ///@{
  ///@}
  ///@name Private Operations
  ///@{



  //**************************************************************************
  //**************************************************************************

  bool RemoveNonConvexBoundary(ModelPart& rModelPart)
  {
    KRATOS_TRY

    if( mEchoLevel > 0 ){
      std::cout<<"   [ REMOVE NON CONVEX BOUNDARY : "<<std::endl;
      //std::cout<<"     Starting Conditions : "<<rModelPart.Conditions().size()<<std::endl;
    }

    double RemovedConditions = rModelPart.NumberOfConditions();

    //***SIZES :::: parameters do define the tolerance in mesh size:
    double critical_angle        = -120;
    double size_for_side_normal  =  mrRemesh.Refine->CriticalRadius;


    unsigned int number_of_nodes = 0;
    if( !rModelPart.IsSubModelPart() )
      number_of_nodes = rModelPart.NumberOfNodes()+1;
    else
      number_of_nodes = rModelPart.GetParentModelPart()->NumberOfNodes()+1;


    std::vector<std::vector<Condition::Pointer> > node_shared_conditions(number_of_nodes); //all domain nodes

    //std::cout<<"     Shared Conditions Size "<<node_shared_conditions.size()<<std::endl;

    for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(); ic!= rModelPart.ConditionsEnd(); ++ic)
    {
      if(ic->IsNot(NEW_ENTITY) && ic->IsNot(TO_ERASE)){
        Geometry< Node<3> >& rConditionGeom = ic->GetGeometry();
        for(unsigned int i=0; i<rConditionGeom.size(); ++i){
          //std::cout<<"["<<ic->Id()<<"] i "<<i<<" condition "<<rConditionGeom[i].Id()<<std::endl;
          if(rConditionGeom[i].Is(TO_ERASE))
            std::cout<<"     WARNING: Released node condition "<<std::endl;

          node_shared_conditions[rConditionGeom[i].Id()].push_back(*(ic.base()));
        }
      }
    }

    //std::cout<<"     Node Shared Conditions (Pair of Condition Nodes) is now set "<<std::endl;


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

    unsigned int initial_cond_size = 0;
    if( !rModelPart.IsSubModelPart() )
      initial_cond_size = rModelPart.NumberOfConditions()+1; //total model part conditions size
    else
      initial_cond_size = rModelPart.GetParentModelPart()->NumberOfConditions()+1;


    unsigned int id = 1;
    unsigned int new_id = 0;
    int RemovedNodes =0;

    for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(); in != rModelPart.NodesEnd(); ++in)
    {

      //angles
      double condition_angle = 0;

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
              TransferUtilities.TransferBoundaryData(pcond, NewCond,*(mrRemesh.Transfer));

              (rModelPart.Conditions()).push_back(pcond);

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

                //review this flag it is reused in laplacian smoothing....
                in->Set(INSIDE);

                Node0.Set(INSIDE);
                Node2.Set(INSIDE);

              }

              if(condition_angle<critical_angle){


                //Path of neighbour conditions in 2D:   (NodeA) ---[0]--- (Node0) ---[1]--- (Node1*) ---[2]--- (Node2) --- [3]--- (NodeB)

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
                VariablesList& variables_list = rModelPart.GetNodalSolutionStepVariablesList();

                PointsArrayType  PointsArray;
                PointsArray.push_back( *(in.base()) );
                //PointsArray.push_back( &Node2 );
                PointsArray.push_back( rConditionGeom2(1) );

                Geometry<Node<3> > geom( PointsArray );

                std::vector<double> ShapeFunctionsN(2);
                std::fill( ShapeFunctionsN.begin(), ShapeFunctionsN.end(), 0.0 );
                ShapeFunctionsN[0] = 0.5;
                ShapeFunctionsN[1] = 0.5;

                MeshDataTransferUtilities DataTransferUtilities;
                DataTransferUtilities.Interpolate2Nodes( geom, ShapeFunctionsN, variables_list, Node0);

                // unsigned int buffer_size = Node0.GetBufferSize();
                // unsigned int step_data_size = rModelPart.GetNodalSolutionStepDataSize();

                // for(unsigned int step = 0; step<buffer_size; ++step)
                // 	{
                // 	  //getting the data of the solution step
                // 	  double* step_data = Node0.SolutionStepData().Data(step);

                // 	  double* node0_data = Node0.SolutionStepData().Data(step);
                // 	  double* node1_data = Node2.SolutionStepData().Data(step);

                // 	  //copying this data in the position of the vector we are interested in
                // 	  for(unsigned int j= 0; j<step_data_size; ++j)
                // 	    {
                // 	      step_data[j] = 0.5*node0_data[j] + 0.5*node1_data[j];
                // 	    }
                // 	}

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

                Geometry<Node<3> >& rConditionGeom0 = node_shared_conditions[Node0.Id()][i]->GetGeometry();
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
                Geometry<Node<3> >& rConditionGeom3 = NewCond->GetGeometry();
                Node<3> & NodeB = rConditionGeom3[1];

                NodeA.Set(BLOCKED);
                NodeB.Set(BLOCKED);
                Node0.Set(BLOCKED);

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
                TransferUtilities.TransferBoundaryData(pcond, NewCond, *(mrRemesh.Transfer) );

                (rModelPart.Conditions()).push_back(pcond);

                RemovedNodes += 1;
                id +=1;

              }
            }
          }
        }
      }
    }

    for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(); in != rModelPart.NodesEnd(); ++in)
    {
      in->Reset(BLOCKED);
    }


    RemovedConditions = rModelPart.Conditions().size() - RemovedConditions;

    if( mEchoLevel > 0 ){
      std::cout<<"     [ CONDITIONS ( removed : "<<RemovedConditions<<" ) ]"<<std::endl;
      std::cout<<"     [ NODES      ( removed : "<<RemovedNodes<<" ) ]"<<std::endl;

      std::cout<<"     Ending   Conditions : "<<rModelPart.Conditions().size()<<"  (Removed nodes: "<< RemovedNodes<<" ) "<<std::endl;
      std::cout<<"     REMOVE NON CONVEX BOUNDARY ]; "<<std::endl;
    }

    if(RemovedNodes)
      return true;
    else
      return false;

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
  RemoveNodesMesherProcess& operator=(RemoveNodesMesherProcess const& rOther);

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
                                  RemoveNodesMesherProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const RemoveNodesMesherProcess& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_REMOVE_NODES_MESHER_PROCESS_H_INCLUDED  defined
