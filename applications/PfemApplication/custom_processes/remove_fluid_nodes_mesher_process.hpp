//
//   Project Name:        KratosPfemFluidApplication $
//   Created by:          $Author:       JMCarbonell $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:           July 2018 $
//   Revision:            $Revision:             0.0 $
//
//

#if !defined(KRATOS_REMOVE_FLUID_NODES_MESHER_PROCESS_H_INCLUDED )
#define  KRATOS_REMOVE_FLUID_NODES_MESHER_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "custom_processes/remove_nodes_mesher_process.hpp"

///VARIABLES used:
//Data:     NORMAL, MASTER_NODES, NEIGHBOUR_NODES, NEIGBOUR_ELEMENTS
//StepData: MEAN_ERROR
//Flags:    (checked) TO_ERASE, BOUNDARY, STRUCTURE, NEW_ENTITY, BLOCKED
//          (set)     TO_ERASE(conditions,nodes)(set), NEW_ENTITY(conditions,nodes)(set), BLOCKED(nodes)->locally, VISITED(nodes)(set)
//          (modified)
//          (reset)   BLOCKED->locally
//(set):=(set in this process)

namespace Kratos
{

///@name Kratos Classes
///@{

/// Remove Fluid Nodes Process for 2D and 3D cases

class RemoveFluidNodesMesherProcess
    : public RemoveNodesMesherProcess
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
      : RemoveNodesMesherProcess(rModelPart,rRemeshingParameters,EchoLevel)
  {
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

 protected:

  ///@name Protected static Member Variables
  ///@{
  ///@}
  ///@name Protected member Variables
  ///@{
  ///@name Protected Operators
  ///@{
  ///@}
  ///@name Protected Operations
  ///@{
  ///@}

  bool RemoveNodesOnDistance(ModelPart& rModelPart, unsigned int& inside_nodes_removed, unsigned int& boundary_nodes_removed, bool& any_condition_removed) override
  {
    KRATOS_TRY

    const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

    //***SIZES :::: parameters do define the tolerance in mesh size:
    double size_for_distance_inside       = 2.0 * mrRemesh.Refine->CriticalRadius; //compared with element radius
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
    unsigned int num_neighbours = 20;

    std::vector<Node<3>::Pointer> neighbours         (num_neighbours);
    std::vector<double>           neighbour_distances(num_neighbours);


    //radius means the distance, if the distance between two nodes is closer to radius -> mark for removing
    double radius=0;
    Node<3> work_point(0,0.0,0.0,0.0);
    unsigned int n_points_in_radius;


    for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(); in != rModelPart.NodesEnd(); ++in)
    {
      if(in->Is(TO_ERASE)){
        any_node_removed = true;
        std::cout<<" TO_ERASE "<<in->Id()<<" "<<in->Coordinates()<<std::endl;
      }

      if( in->IsNot(NEW_ENTITY) && in->IsNot(BLOCKED) && in->IsNot(SOLID) && in->IsNot(TO_ERASE) )
      {
        radius = size_for_distance_inside;

        work_point[0]=in->X();
        work_point[1]=in->Y();
        work_point[2]=in->Z();

        unsigned int FreeSurfaceNeighbours=0;
        unsigned int RigidNeighbours=0;
        WeakPointerVector< Node < 3 > >& NeighbourNodes = in->GetValue(NEIGHBOUR_NODES);
        unsigned int NumberOfNeighbourNodes = NeighbourNodes.size();
        for(WeakPointerVector< Node <3> >::iterator nn = NeighbourNodes.begin();nn != NeighbourNodes.end(); ++nn)
        {
          if(nn->Is(BLOCKED) || nn->Is(SOLID)){
            ++RigidNeighbours;
          }
          if(nn->Is(FREE_SURFACE)){
            ++FreeSurfaceNeighbours;
          }
        }

        if(in->Is(FREE_SURFACE)){ // it must be more difficult to erase a free_surface node, otherwise, lot of volume is lost

          if( RigidNeighbours == NumberOfNeighbourNodes ){
            radius*=0.25;
          }
          else{
            radius*=0.85;
          }

        }

        n_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, neighbours.begin(),neighbour_distances.begin(), num_neighbours);

        if (n_points_in_radius>1)
        {
          //std::cout<<"     Points in Radius "<< n_points_in_radius<<" radius "<<radius<<std::endl;
          if ( in->IsNot(BOUNDARY) )
          {
            if( this->mrRemesh.Refine->RemovingOptions.Is(MesherUtilities::REMOVE_NODES_ON_DISTANCE) ){

              // if the node is close to the surface, do not erase it, move to a mean (laplacian) position
              if(in->IsNot(FREE_SURFACE) && FreeSurfaceNeighbours>=dimension){
                this->MoveInsideNode((*in));
              }
              else{
                if( !this->CheckEngagedNode((*in),neighbours,neighbour_distances,n_points_in_radius) ){ //we release the node if no other nodes neighbours are being erased
                  in->Set(TO_ERASE);
                  any_node_removed = true;
                  ++inside_nodes_removed;
                  //std::cout<<"     Distance Criterion Node ["<<in->Id()<<"] TO_ERASE "<<in->Coordinates()<<std::endl;
                }
              }
            }

          }
          else{

            // std::cout<<"  Remove close boundary nodes: Candidate ["<<in->Id()<<"]"<<std::endl;
            bool engaged_node = false;
            unsigned int counter = 0;
            for(std::vector<Node<3>::Pointer>::iterator nn=neighbours.begin(); nn!=neighbours.begin() + n_points_in_radius ; ++nn)
            {

              if ( (*nn)->Is(BOUNDARY) && (neighbour_distances[counter] < 2.0 * size_for_distance_boundary) && (neighbour_distances[counter] > 0.0) )
              {
                if((*nn)->Is(TO_ERASE) || (*nn)->Is(BLOCKED)){
                  engaged_node = true;
                  break;
                }
              }

              ++counter;
            }

            if(!engaged_node && in->IsNot(BLOCKED)){ //Can be inserted in the boundary refine
              in->Set(TO_ERASE);
              ++boundary_nodes_removed;
              //std::cout<<"     Removed Boundary Node ["<<in->Id()<<"] on Distance "<<std::endl;
            }

          }


        }
      }
    }

    if( boundary_nodes_removed > 0 )
      this->MoveBoundaries(rModelPart,boundary_nodes_removed);

    if( boundary_nodes_removed > 0 )
      any_node_removed = true;


    bool critical_nodes_removed = false;
    critical_nodes_removed = this->EraseCriticalNodes(rModelPart,inside_nodes_removed);

    if( any_node_removed || critical_nodes_removed )
      any_node_removed = true;


    //Build boundary after removing boundary nodes due distance criterion
    if( this->mEchoLevel > 0 ){
      std::cout<<"boundary_nodes_removed "<<boundary_nodes_removed<<std::endl;
      std::cout<<"inside_nodes_removed "<<inside_nodes_removed<<std::endl;
      std::cout<<"critical_nodes_removed "<<critical_nodes_removed<<std::endl;
    }

    //Build boundary after removing boundary nodes due distance criterion
    // if(boundary_nodes_removed){
    //   std::cout<<"     Rebuild boundary needed after release on Distance "<<std::endl;
    //   //any_condition_removed = this->RebuildBoundary(rModelPart);
    // }

    return any_node_removed;

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

  void MoveBoundaries(ModelPart& rModelPart, unsigned int& boundary_nodes_removed)
  {

    KRATOS_TRY

    for(ModelPart::ConditionsContainerType::const_iterator i_cond = rModelPart.ConditionsBegin(); i_cond != rModelPart.ConditionsEnd(); ++i_cond)
    {
      Condition::GeometryType& rGeometry = i_cond->GetGeometry();
      unsigned int NumberOfVertices = rGeometry.size();

      unsigned int counter = 0;
      int id = -1;
      for(unsigned int i=0; i<NumberOfVertices; ++i)
      {
        if(rGeometry[i].Is(TO_ERASE)){
          id = i;
          ++counter;
        }

      }

      if(counter==1 && id>=0){

        if( this->MoveBoundaryNode(rGeometry[id]) ){
          rGeometry[id].Set(TO_ERASE,false);
          --boundary_nodes_removed;
        }

      }
    }

    KRATOS_CATCH( "" )
  }

  //**************************************************************************
  //**************************************************************************

  bool CheckApproachingPoint(Node<3>& rNode, Node<3>& rEdgeNode)
  {
    KRATOS_TRY

    bool approaching_point = false;

    array_1d<double,3> EdgeVelocity = rEdgeNode.FastGetSolutionStepValue(VELOCITY);

    array_1d<double,3> Direction = rEdgeNode.FastGetSolutionStepValue(NORMAL);
    if( norm_2(Direction) )
      Direction /= norm_2(Direction);

    array_1d<double,3> Distance = (rEdgeNode.Coordinates()-rNode.Coordinates());
    double distance = inner_prod(Distance, Direction);
    Distance = distance * Direction;

    array_1d<double,3> VelocityDirection = rNode.FastGetSolutionStepValue(VELOCITY);

    double velocity = norm_2(VelocityDirection);
    if(velocity!=0)
       VelocityDirection/=velocity;

    //different velocity directions
    if( inner_prod( EdgeVelocity, VelocityDirection ) < 0 ){
      //node velocity direction towards edge
      if( inner_prod( Distance, VelocityDirection ) > 0 )
        approaching_point = true;
      else
        approaching_point = false;
    } //same velocity directions
    else if( velocity > 0.1 * norm_2(EdgeVelocity) ){
      //node velocity direction towards edge
      if( inner_prod( Distance, VelocityDirection ) > 0 )
        approaching_point = true;
      else
        approaching_point = false;
    }

    return approaching_point;

    KRATOS_CATCH( "" )
  }

  //**************************************************************************
  //**************************************************************************

  bool CheckApproachingEdge(Node<3>& rNode, Node<3>& rEdgeNodeA, Node<3>& rEdgeNodeB)
  {
    KRATOS_TRY

    bool approaching_edge = false;

    array_1d<double,3> EdgeVelocity = 0.5 * (rEdgeNodeA.FastGetSolutionStepValue(VELOCITY)+rEdgeNodeB.FastGetSolutionStepValue(VELOCITY));
    array_1d<double,3> MidPoint = 0.5 * (rEdgeNodeA.Coordinates() + rEdgeNodeB.Coordinates());

    array_1d<double,3> Direction;
    this->GetDirectionToEdge(Direction,rEdgeNodeA,rEdgeNodeB);

    array_1d<double,3> Distance = (MidPoint-rNode.Coordinates());
    double distance = inner_prod(Distance, Direction);
    Distance = distance * Direction;

    array_1d<double,3> VelocityDirection = rNode.FastGetSolutionStepValue(VELOCITY);

    double velocity = norm_2(VelocityDirection);
    if(velocity!=0)
       VelocityDirection/=velocity;

    //different velocity directions
    if( inner_prod( EdgeVelocity, VelocityDirection ) < 0 ){
      //node velocity direction towards edge
      if( inner_prod( Distance, VelocityDirection ) > 0 )
        approaching_edge = true;
      else
        approaching_edge = false;
    } //same velocity directions
    else if( velocity > 0.1 * norm_2(EdgeVelocity) ){
      //node velocity direction towards edge
      if( inner_prod( Distance, VelocityDirection ) > 0 )
        approaching_edge = true;
      else
        approaching_edge = false;
    }

    return approaching_edge;

    KRATOS_CATCH( "" )
  }


  //**************************************************************************
  //**************************************************************************

  bool CheckApproachingFace(Node<3>& rNode, Node<3>& rEdgeNodeA, Node<3>& rEdgeNodeB, Node<3>& rEdgeNodeC)
  {
    KRATOS_TRY

    bool approaching_face = false;

    array_1d<double,3> EdgeVelocity = (1.0/3.0) * (rEdgeNodeA.FastGetSolutionStepValue(VELOCITY)+rEdgeNodeB.FastGetSolutionStepValue(VELOCITY)+rEdgeNodeC.FastGetSolutionStepValue(VELOCITY));
    array_1d<double,3> MidPoint = (1.0/3.0) * (rEdgeNodeA.Coordinates() + rEdgeNodeB.Coordinates() + rEdgeNodeC.Coordinates());

    array_1d<double,3> Direction;
    this->GetDirectionToFace(Direction,rEdgeNodeA,rEdgeNodeB,rEdgeNodeC);

    array_1d<double,3> Distance = (MidPoint-rNode.Coordinates());
    double distance = inner_prod(Distance, Direction);
    Distance = distance * Direction;

    array_1d<double,3> VelocityDirection = rNode.FastGetSolutionStepValue(VELOCITY);

    double velocity = norm_2(VelocityDirection);
    if(velocity!=0)
       VelocityDirection/=velocity;

    //different velocity directions
    if( inner_prod( EdgeVelocity, VelocityDirection ) < 0 ){
      //node velocity direction towards edge
      if( inner_prod( Distance, VelocityDirection ) > 0 )
        approaching_face = true;
      else
        approaching_face = false;
    } //same velocity directions
    else if( velocity > 0.1 * norm_2(EdgeVelocity) ){
      //node velocity direction towards edge
      if( inner_prod( Distance, VelocityDirection ) > 0 )
        approaching_face = true;
      else
        approaching_face = false;
    }

    return approaching_face;

    KRATOS_CATCH( "" )
  }

  //**************************************************************************
  //**************************************************************************

  double GetDistanceToNode(Node<3>& rNode, Node<3>& rEdgeNode, array_1d<double,3>& rDirection)
  {
    KRATOS_TRY

    array_1d<double,3> Direction = rEdgeNode.FastGetSolutionStepValue(NORMAL);
    if( norm_2(Direction) )
      Direction /= norm_2(Direction);

    array_1d<double,3> Distance = (rEdgeNode.Coordinates()-rNode.Coordinates());
    double distance = fabs(inner_prod(Distance, Direction));

    return distance;

    KRATOS_CATCH( "" )
  }

  //**************************************************************************
  //**************************************************************************

  double GetDistanceToEdge(Node<3>& rNode, Node<3>& rEdgeNodeA, Node<3>& rEdgeNodeB, array_1d<double,3>& rDirection)
  {
    KRATOS_TRY

    array_1d<double,3> MidPoint  = 0.5 * (rEdgeNodeA.Coordinates() + rEdgeNodeB.Coordinates()) ;

    this->GetDirectionToEdge(rDirection,rEdgeNodeA,rEdgeNodeB);

    double distance = fabs(inner_prod((MidPoint-rNode.Coordinates()), rDirection));

    return distance;

    KRATOS_CATCH( "" )
  }

  //**************************************************************************
  //**************************************************************************

  double GetDistanceToFace(Node<3>& rNode, Node<3>& rEdgeNodeA, Node<3>& rEdgeNodeB, Node<3>& rEdgeNodeC, array_1d<double,3>& rDirection)
  {
    KRATOS_TRY

    array_1d<double,3> MidPoint  = (1.0/3.0) * (rEdgeNodeA.Coordinates() + rEdgeNodeB.Coordinates() + rEdgeNodeC.Coordinates() ) ;

    this->GetDirectionToFace(rDirection,rEdgeNodeA,rEdgeNodeB,rEdgeNodeC);

    double distance = fabs(inner_prod((MidPoint-rNode.Coordinates()), rDirection));

    return distance;

    KRATOS_CATCH( "" )
  }


  //**************************************************************************
  //**************************************************************************

  void GetDirectionToEdge(array_1d<double,3>& rDirection, Node<3>& rEdgeNodeA, Node<3>& rEdgeNodeB)
  {
    //get wall direction from normals:
    if( rEdgeNodeA.FastGetSolutionStepValue(SHRINK_FACTOR) == 1 && rEdgeNodeB.FastGetSolutionStepValue(SHRINK_FACTOR) == 1 ){
      rDirection = 0.5 * (rEdgeNodeA.FastGetSolutionStepValue(NORMAL) + rEdgeNodeB.FastGetSolutionStepValue(NORMAL) ) ;
    }
    else{
      if( rEdgeNodeA.FastGetSolutionStepValue(SHRINK_FACTOR) == 1 ){
        rDirection = rEdgeNodeA.FastGetSolutionStepValue(NORMAL);
      }
      else if( rEdgeNodeB.FastGetSolutionStepValue(SHRINK_FACTOR) == 1 ){
        rDirection = rEdgeNodeB.FastGetSolutionStepValue(NORMAL);
      }
      else{
        rDirection = 0.5 * (rEdgeNodeA.FastGetSolutionStepValue(NORMAL) + rEdgeNodeB.FastGetSolutionStepValue(NORMAL) ) ;
      }
    }

    if( norm_2(rDirection) )
      rDirection /= norm_2(rDirection);

  }

  //**************************************************************************
  //**************************************************************************

  void GetDirectionToFace(array_1d<double,3>& rDirection, Node<3>& rEdgeNodeA, Node<3>& rEdgeNodeB, Node<3>& rEdgeNodeC)
  {

    array_1d<double,3> VectorB = rEdgeNodeB.Coordinates() - rEdgeNodeA.Coordinates();
    array_1d<double,3> VectorC = rEdgeNodeC.Coordinates() - rEdgeNodeA.Coordinates();

    MathUtils<double>::CrossProduct(rDirection,VectorB,VectorC);

    if( norm_2(rDirection) )
      rDirection /= norm_2(rDirection);

  }

  //**************************************************************************
  //**************************************************************************

  bool EraseCriticalNodes(ModelPart& rModelPart, unsigned int& inside_nodes_removed)
  {

    KRATOS_TRY

    bool any_node_removed = false;
    unsigned int erased_nodes=0;

    ModelPart::NodesContainerType LayerNodes;

    MesherUtilities MesherUtils;
    double MaxRelativeVelocity = 1.5; //arbitrary value, will depend on time step (AF)

    for(ModelPart::ElementsContainerType::const_iterator ie = rModelPart.ElementsBegin(); ie != rModelPart.ElementsEnd(); ++ie)
    {
      Element::GeometryType& rGeometry = ie->GetGeometry();
      const unsigned int NumberOfNodes = rGeometry.size();

      bool wall_boundary = false;
      for(unsigned int i=0; i<NumberOfNodes; ++i)
      {
        if(rGeometry[i].Is(BLOCKED) || rGeometry[i].Is(SOLID)){
          wall_boundary = true;
          break;
        }
      }

      if( wall_boundary ){

        bool speedy_approach = MesherUtils.CheckRelativeVelocities(rGeometry, MaxRelativeVelocity);

        if( speedy_approach ){

          for(unsigned int i=0; i<NumberOfNodes; ++i)
          {
            if((rGeometry[i].IsNot(BLOCKED) && rGeometry[i].IsNot(SOLID))){
              LayerNodes.push_back(rGeometry(i));
            }
          }

        }
      }
    }

    const int nnodes = LayerNodes.size();
    // const unsigned int& dimension = rModelPart.GetProcessInfo()[SPACE_DIMENSION];

    if (nnodes != 0)
    {
      ModelPart::NodesContainerType::iterator it_begin = LayerNodes.begin();

      unsigned int inside_nodes_removed_accum = 0;
      #pragma omp parallel for reduction(+:inside_nodes_removed_accum,erased_nodes)
      for (int i = 0; i < nnodes; ++i)
      {
        ModelPart::NodesContainerType::iterator it = it_begin + i;

        double MinimumDistance = std::numeric_limits<double>::max();
        array_1d<double,3> Direction;
        // PointsArrayType Vertices(dimension);
        double distance = 0;
        unsigned int face = 0;

        WeakPointerVector<Element>& rNE = it->GetValue(NEIGHBOUR_ELEMENTS);

 	for(WeakPointerVector<Element>::iterator ie = rNE.begin(); ie!=rNE.end(); ++ie)
        {
          WeakPointerVector<Element>& rE = ie->GetValue(NEIGHBOUR_ELEMENTS);

          DenseMatrix<unsigned int> lpofa; //connectivities of points defining faces

          distance = 0;
          face = 0;
          for(WeakPointerVector<Element>::iterator je = rE.begin(); je!=rE.end(); ++je)
          {
            if (ie->Id() == je->Id()){

              GeometryType& rGeometry = ie->GetGeometry();

              rGeometry.NodesInFaces(lpofa);

              const unsigned int NumberOfNodes = rGeometry.size();

              unsigned int wall_boundary = 0;
              std::vector<unsigned int> wall_nodes;
              for(unsigned int j=1; j<NumberOfNodes; ++j)
              {
                if( rGeometry[lpofa(j,face)].Is(BLOCKED) || rGeometry[lpofa(j,face)].Is(SOLID) ){
                  ++wall_boundary;
                  wall_nodes.push_back(j);
                }
              }

              if( wall_boundary == NumberOfNodes-1 ){ //node to face criterion

                if( NumberOfNodes == 3 ){

                  if( this->CheckApproachingEdge(*it,rGeometry[lpofa(1,face)], rGeometry[lpofa(2,face)]) ){
                    distance = this->GetDistanceToEdge(*it, rGeometry[lpofa(1,face)], rGeometry[lpofa(2,face)], Direction);
                    if( distance < MinimumDistance ){
                      MinimumDistance = distance;
                      // Vertices[0] = rGeometry[lpofa(1,face)];
                      // Vertices[1] = rGeometry[lpofa(2,face)];
                    }
                  }

                }
                else if( NumberOfNodes == 4 ){

                  if( this->CheckApproachingFace(*it,rGeometry[lpofa(1,face)], rGeometry[lpofa(2,face)], rGeometry[lpofa(3,face)]) ){
                    distance = this->GetDistanceToFace(*it, rGeometry[lpofa(1,face)], rGeometry[lpofa(2,face)], rGeometry[lpofa(3,face)], Direction);
                    if( distance < MinimumDistance ){
                      MinimumDistance = distance;
                      // Vertices[0] = rGeometry[lpofa(1,face)];
                      // Vertices[1] = rGeometry[lpofa(2,face)];
                      // Vertices[2] = rGeometry[lpofa(3,face)];
                    }
                  }
                }
              }
              else if( wall_boundary == NumberOfNodes-2 ){ //node to edge criterion

                if( NumberOfNodes == 3 ){

                  if( this->CheckApproachingPoint(*it,rGeometry[lpofa(wall_nodes.front(),face)]) ){
                    distance = this->GetDistanceToNode(*it, rGeometry[lpofa(wall_nodes.front(),face)], Direction);
                    if( distance < MinimumDistance ){
                      MinimumDistance = distance;
                    }
                  }

                }
                else if( NumberOfNodes == 4 ){

                  if( this->CheckApproachingEdge(*it,rGeometry[lpofa(wall_nodes.front(),face)],rGeometry[lpofa(wall_nodes.back(),face)]) ){
                    distance = this->GetDistanceToEdge(*it, rGeometry[lpofa(wall_nodes.front(),face)], rGeometry[lpofa(wall_nodes.back(),face)], Direction);
                    if( distance < MinimumDistance ){
                      MinimumDistance = distance;
                    }
                  }
                }

              }
            }
            ++face;
          }
        }

        if( MinimumDistance < 0.25 * mrRemesh.Refine->CriticalRadius ){
          it->Set(TO_ERASE);
          ++erased_nodes;
          ++inside_nodes_removed_accum;
        }
        else if( MinimumDistance < 1.5 * mrRemesh.Refine->CriticalRadius ){
          distance = (1.5 * mrRemesh.Refine->CriticalRadius - MinimumDistance);
          this->MoveLayerNode(*it, Direction, distance);
        }

      }

	  inside_nodes_removed += inside_nodes_removed_accum;

    }
    if(erased_nodes>0){
      //std::cout<<" Layer Nodes removed "<<erased_nodes<<std::endl;
      any_node_removed = true;
    }

    return any_node_removed;

    KRATOS_CATCH( "" )
  }


  //**************************************************************************
  //**************************************************************************

  bool MoveBoundaryNode(Node<3>& rNode)
  {

    KRATOS_TRY

    bool moved_node = false;
    //std::cout<<" Boundary to Move Pre ["<<rNode.Id()<<"] "<<rNode.Coordinates()<<std::endl;
    unsigned int FreeSurfaceNodes = 0;
    WeakPointerVector< Node < 3 > >& NeighbourNodes = rNode.GetValue(NEIGHBOUR_NODES);
    WeakPointerVector< Node < 3 > > FreeNeighbours;
    for(WeakPointerVector< Node <3> >::iterator nn = NeighbourNodes.begin(); nn != NeighbourNodes.end(); ++nn)
    {
      if( nn->Is(FREE_SURFACE) ){
        FreeNeighbours.push_back(*(nn.base()));
            ++FreeSurfaceNodes;
      }
    }

    if( FreeSurfaceNodes == 2 )
    {
      array_1d<double,3> MidPoint  = 0.5 * (FreeNeighbours.front().Coordinates()+FreeNeighbours.back().Coordinates());
      array_1d<double,3> Direction = (FreeNeighbours.front().Coordinates()-FreeNeighbours.back().Coordinates());

      if(norm_2(Direction))
        Direction/=norm_2(Direction);

      array_1d<double,3> Displacement = inner_prod( (MidPoint-rNode.Coordinates()), Direction ) * Direction;
      noalias(rNode.Coordinates()) += Displacement;
      noalias(rNode.FastGetSolutionStepValue(DISPLACEMENT))   += Displacement;
      noalias(rNode.FastGetSolutionStepValue(DISPLACEMENT,1)) += Displacement;

      for(WeakPointerVector< Node < 3 > >::iterator fsn = FreeNeighbours.begin(); fsn != FreeNeighbours.end(); ++fsn)
      {
        noalias(rNode.FastGetSolutionStepValue(VELOCITY))       += fsn->FastGetSolutionStepValue(VELOCITY);
        noalias(rNode.FastGetSolutionStepValue(VELOCITY,1))     += fsn->FastGetSolutionStepValue(VELOCITY,1);
        noalias(rNode.FastGetSolutionStepValue(ACCELERATION))   += fsn->FastGetSolutionStepValue(ACCELERATION);
        noalias(rNode.FastGetSolutionStepValue(ACCELERATION,1)) += fsn->FastGetSolutionStepValue(ACCELERATION,1);
        rNode.FastGetSolutionStepValue(PRESSURE)                += fsn->FastGetSolutionStepValue(PRESSURE);
        rNode.FastGetSolutionStepValue(PRESSURE_VELOCITY)       += fsn->FastGetSolutionStepValue(PRESSURE_VELOCITY);
        rNode.FastGetSolutionStepValue(PRESSURE_VELOCITY,1)     += fsn->FastGetSolutionStepValue(PRESSURE_VELOCITY,1);
      }


      double quotient = 1.0/double(FreeSurfaceNodes+1);
      rNode.FastGetSolutionStepValue(VELOCITY)       *= quotient;
      rNode.FastGetSolutionStepValue(VELOCITY,1)     *= quotient;
      rNode.FastGetSolutionStepValue(ACCELERATION)   *= quotient;
      rNode.FastGetSolutionStepValue(ACCELERATION,1) *= quotient;
      rNode.FastGetSolutionStepValue(PRESSURE)            *= quotient;
      rNode.FastGetSolutionStepValue(PRESSURE_VELOCITY)   *= quotient;
      rNode.FastGetSolutionStepValue(PRESSURE_VELOCITY,1) *= quotient;

      moved_node = true;

    }
    else if(FreeSurfaceNodes > 2) {

      array_1d<double,3> MidPoint;
      noalias(MidPoint) = ZeroVector(3);
      double quotient = 1.0/double(FreeSurfaceNodes);
      for(WeakPointerVector< Node < 3 > >::iterator fsn = FreeNeighbours.begin(); fsn != FreeNeighbours.end(); ++fsn)
      {
        MidPoint += fsn->Coordinates();

        noalias(rNode.FastGetSolutionStepValue(VELOCITY))       += fsn->FastGetSolutionStepValue(VELOCITY);
        noalias(rNode.FastGetSolutionStepValue(VELOCITY,1))     += fsn->FastGetSolutionStepValue(VELOCITY,1);
        noalias(rNode.FastGetSolutionStepValue(ACCELERATION))   += fsn->FastGetSolutionStepValue(ACCELERATION);
        noalias(rNode.FastGetSolutionStepValue(ACCELERATION,1)) += fsn->FastGetSolutionStepValue(ACCELERATION,1);
        rNode.FastGetSolutionStepValue(PRESSURE)                += fsn->FastGetSolutionStepValue(PRESSURE);
        rNode.FastGetSolutionStepValue(PRESSURE_VELOCITY)       += fsn->FastGetSolutionStepValue(PRESSURE_VELOCITY);
        rNode.FastGetSolutionStepValue(PRESSURE_VELOCITY,1)     += fsn->FastGetSolutionStepValue(PRESSURE_VELOCITY,1);
      }
      MidPoint *= quotient;
      array_1d<double,3> Normal = rNode.FastGetSolutionStepValue(NORMAL);

      if(norm_2(Normal))
        Normal/=norm_2(Normal);

      array_1d<double,3> Displacement = (MidPoint-rNode.Coordinates()) - inner_prod( (MidPoint-rNode.Coordinates()), Normal ) * Normal;
      noalias(rNode.Coordinates()) += Displacement;
      noalias(rNode.FastGetSolutionStepValue(DISPLACEMENT)) += Displacement;
      noalias(rNode.FastGetSolutionStepValue(DISPLACEMENT,1)) += Displacement;

      quotient = 1.0/double(FreeSurfaceNodes+1);
      rNode.FastGetSolutionStepValue(VELOCITY)       *= quotient;
      rNode.FastGetSolutionStepValue(VELOCITY,1)     *= quotient;
      rNode.FastGetSolutionStepValue(ACCELERATION)   *= quotient;
      rNode.FastGetSolutionStepValue(ACCELERATION,1) *= quotient;
      rNode.FastGetSolutionStepValue(PRESSURE)            *= quotient;
      rNode.FastGetSolutionStepValue(PRESSURE_VELOCITY)   *= quotient;
      rNode.FastGetSolutionStepValue(PRESSURE_VELOCITY,1) *= quotient;

      moved_node = true;
    }
    else{
      std::cout<<" Boundary node with only one FREE_SURFACE neighbour "<<std::endl;
    }

    //std::cout<<" Boundary to Move Post ["<<rNode.Id()<<"] "<<rNode.Coordinates()<<std::endl;

    return moved_node;

    KRATOS_CATCH( "" )
  }

  //**************************************************************************
  //**************************************************************************

  void MoveInsideNode(Node<3>& rNode)
  {

    KRATOS_TRY

    WeakPointerVector< Node < 3 > >& NeighbourNodes = rNode.GetValue(NEIGHBOUR_NODES);
    unsigned int NumberOfNeighbourNodes = NeighbourNodes.size();

    //std::cout<<" Moved Node Pre ["<<rNode.Id()<<"] Displacement"<<rNode.FastGetSolutionStepValue(DISPLACEMENT)<<" Position "<<rNode.Coordinates()<<" Initial Position "<<rNode.GetInitialPosition()<<std::endl;

    //array_1d<double,3> CurrentPosition = rNode.Coordinates();

    for(WeakPointerVector< Node <3> >::iterator nn = NeighbourNodes.begin();nn != NeighbourNodes.end(); ++nn)
    {
      noalias(rNode.Coordinates())                            += nn->Coordinates();
      noalias(rNode.FastGetSolutionStepValue(DISPLACEMENT))   += nn->FastGetSolutionStepValue(DISPLACEMENT);
      noalias(rNode.FastGetSolutionStepValue(DISPLACEMENT,1)) += nn->FastGetSolutionStepValue(DISPLACEMENT,1);
      noalias(rNode.FastGetSolutionStepValue(VELOCITY))       += nn->FastGetSolutionStepValue(VELOCITY);
      noalias(rNode.FastGetSolutionStepValue(VELOCITY,1))     += nn->FastGetSolutionStepValue(VELOCITY,1);
      noalias(rNode.FastGetSolutionStepValue(ACCELERATION))   += nn->FastGetSolutionStepValue(ACCELERATION);
      noalias(rNode.FastGetSolutionStepValue(ACCELERATION,1)) += nn->FastGetSolutionStepValue(ACCELERATION,1);
      rNode.FastGetSolutionStepValue(PRESSURE)                += nn->FastGetSolutionStepValue(PRESSURE);
      rNode.FastGetSolutionStepValue(PRESSURE_VELOCITY)       += nn->FastGetSolutionStepValue(PRESSURE_VELOCITY);
      rNode.FastGetSolutionStepValue(PRESSURE_VELOCITY,1)     += nn->FastGetSolutionStepValue(PRESSURE_VELOCITY,1);
    }

    double quotient = 1.0/double(NumberOfNeighbourNodes+1);

    rNode.Coordinates() *= quotient;
    rNode.FastGetSolutionStepValue(DISPLACEMENT)   *= quotient;
    rNode.FastGetSolutionStepValue(DISPLACEMENT,1) *= quotient;
    //rNode.FastGetSolutionStepValue(DISPLACEMENT)   += rNode.Coordinates()-CurrentPosition;
    //rNode.FastGetSolutionStepValue(DISPLACEMENT,1) += rNode.Coordinates()-CurrentPosition;
    rNode.GetInitialPosition() = (rNode.Coordinates() - rNode.FastGetSolutionStepValue(DISPLACEMENT));
    rNode.FastGetSolutionStepValue(VELOCITY)       *= quotient;
    rNode.FastGetSolutionStepValue(VELOCITY,1)     *= quotient;
    rNode.FastGetSolutionStepValue(ACCELERATION)   *= quotient;
    rNode.FastGetSolutionStepValue(ACCELERATION,1) *= quotient;
    rNode.FastGetSolutionStepValue(PRESSURE)            *= quotient;
    rNode.FastGetSolutionStepValue(PRESSURE_VELOCITY)   *= quotient;
    rNode.FastGetSolutionStepValue(PRESSURE_VELOCITY,1) *= quotient;

    //std::cout<<" Moved Node Post ["<<rNode.Id()<<"] Displacement"<<rNode.FastGetSolutionStepValue(DISPLACEMENT)<<" Position "<<rNode.Coordinates()<<" Initial Position "<<rNode.GetInitialPosition()<<std::endl;

    KRATOS_CATCH( "" )
  }

  //**************************************************************************
  //**************************************************************************

  void MoveLayerNode(Node<3>& rNode, const array_1d<double,3>& rDirection, const double& rDistance)
  {

    KRATOS_TRY

    //std::cout<<" Moved Node Pre ["<<rNode.Id()<<"] Displacement"<<rNode.FastGetSolutionStepValue(DISPLACEMENT)<<" Position "<<rNode.Coordinates()<<std::endl;

    const array_1d<double,3>& VelocityDirection = rNode.FastGetSolutionStepValue(VELOCITY);

    double sign = 1;
    if( inner_prod(VelocityDirection,rDirection) > 0 )
      sign*=(-1);

    noalias(rNode.Coordinates()) += sign * rDistance * rDirection;
    rNode.FastGetSolutionStepValue(DISPLACEMENT)   += sign * rDistance * rDirection;
    rNode.FastGetSolutionStepValue(DISPLACEMENT,1) += sign * rDistance * rDirection;

    //std::cout<<" Moved Node Post ["<<rNode.Id()<<"] Displacement"<<rNode.FastGetSolutionStepValue(DISPLACEMENT)<<" Position "<<rNode.Coordinates()<<std::endl;

    KRATOS_CATCH( "" )
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

#endif // KRATOS_REMOVE_FLUID_NODES_MESHER_PROCESS_H_INCLUDED  defined
