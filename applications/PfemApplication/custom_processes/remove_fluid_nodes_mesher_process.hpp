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
      }

      if( in->IsNot(NEW_ENTITY) &&  in->IsNot(INLET) && in->IsNot(RIGID) && in->IsNot(SOLID) && in->IsNot(TO_ERASE) )
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
          if(nn->Is(RIGID) || nn->Is(SOLID)){
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
                this->MoveNodeToMeanPosition((*in));
              }
              else{
                if( !this->CheckEngagedNode((*in),neighbours,neighbour_distances,n_points_in_radius) ){ //we release the node if no other nodes neighbours are being erased
                  in->Set(TO_ERASE);
                  any_node_removed = true;
                  ++inside_nodes_removed;
                  //std::cout<<"     Distance Criterion Node ["<<in->Id()<<"] TO_ERASE "<<std::endl;
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
                if((*nn)->Is(TO_ERASE)){
                  engaged_node = true;
                  break;
                }
              }

              ++counter;
            }

            if(!engaged_node){ //Can be inserted in the boundary refine
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
    critical_nodes_removed = EraseCriticalNodes(rModelPart,inside_nodes_removed);

    if( any_node_removed || critical_nodes_removed )
      any_node_removed = true;


    //Build boundary after removing boundary nodes due distance criterion
    if( this->mEchoLevel > 0 ){
      std::cout<<"boundary_nodes_removed "<<boundary_nodes_removed<<std::endl;
      std::cout<<"inside_nodes_removed "<<inside_nodes_removed<<std::endl;
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

        if( this->MoveBoundaryNodeToMeanPosition(rGeometry[id]) ){
          rGeometry[id].Set(TO_ERASE,false);
          --boundary_nodes_removed;
        }

      }
    }

    KRATOS_CATCH( "" )
  }

  //**************************************************************************
  //**************************************************************************

  bool MoveBoundaryNodeToMeanPosition(Node<3>& rNode)
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

  void MoveNodeToMeanPosition(Node<3>& rNode)
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

  void MoveNodeToNonRigidMeanPosition(Node<3>& rNode, Node<3>& rEdgeNodeA, Node<3>& rEdgeNodeB, const double& rMeanEdgeLength)
  {

    KRATOS_TRY

    array_1d<double,3> VelocityDirection = rNode.FastGetSolutionStepValue(VELOCITY);

    double velocity = norm_2(VelocityDirection);
    if(velocity!=0)
       VelocityDirection/=velocity;

    double advance  = norm_2(rNode.FastGetSolutionStepValue(DISPLACEMENT)-rNode.FastGetSolutionStepValue(DISPLACEMENT,1));
    double distance = this->GetDistanceToEdge(rNode,rEdgeNodeA,rEdgeNodeB);

    VelocityDirection *= (-1)*((rMeanEdgeLength-distance)+advance);

    //std::cout<<" Moved Node Pre ["<<rNode.Id()<<"] Displacement"<<rNode.FastGetSolutionStepValue(DISPLACEMENT)<<" Position "<<rNode.Coordinates()<<" Initial Position "<<rNode.GetInitialPosition()<<std::endl;

    noalias(rNode.Coordinates()) += VelocityDirection;
    rNode.FastGetSolutionStepValue(DISPLACEMENT)   += VelocityDirection;
    rNode.FastGetSolutionStepValue(DISPLACEMENT,1) += VelocityDirection;
    //rNode.GetInitialPosition() = (rNode.Coordinates() - rNode.FastGetSolutionStepValue(DISPLACEMENT));


    //std::cout<<" Moved Node Post ["<<rNode.Id()<<"] Displacement"<<rNode.FastGetSolutionStepValue(DISPLACEMENT)<<" Position "<<rNode.Coordinates()<<" Initial Position "<<rNode.GetInitialPosition()<<std::endl;

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

    array_1d<double,3> Distance = (MidPoint-rNode.Coordinates());

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

    //check approach magnitude
    if( approaching_edge == true) {

      array_1d<double,3> Direction;

      this->GetDirectionToEdge(Direction,rEdgeNodeA,rEdgeNodeB);

      double distance = fabs(inner_prod((MidPoint-rNode.Coordinates()), Direction));

      double advance  = fabs(inner_prod((rNode.FastGetSolutionStepValue(DISPLACEMENT)-rNode.FastGetSolutionStepValue(DISPLACEMENT,1)), Direction));

      if( 3*advance < distance )
        approaching_edge = false;

      //check also velocity direction: (when normal well calculated maybe it can be deactivated)
      distance = fabs(inner_prod((MidPoint-rNode.Coordinates()), VelocityDirection));

      advance  = fabs(inner_prod((rNode.FastGetSolutionStepValue(DISPLACEMENT)-rNode.FastGetSolutionStepValue(DISPLACEMENT,1)), VelocityDirection));

      if( 3*advance > distance )
        approaching_edge = true;
      //check also velocity direction: (optional)

      //std::cout<<" aproaching edge: "<<approaching_edge<<" (distance: "<< distance <<" advance: "<<advance<<")"<<std::endl;
    }

    return approaching_edge;

    KRATOS_CATCH( "" )
  }


  //**************************************************************************
  //**************************************************************************

  double GetDistanceToEdge(Node<3>& rNode, Node<3>& rEdgeNodeA, Node<3>& rEdgeNodeB)
  {
    KRATOS_TRY

    array_1d<double,3> MidPoint  = 0.5 * (rEdgeNodeA.Coordinates() + rEdgeNodeB.Coordinates()) ;

    array_1d<double,3> Direction;
    this->GetDirectionToEdge(Direction,rEdgeNodeA,rEdgeNodeB);

    double distance = fabs(inner_prod((MidPoint-rNode.Coordinates()), Direction));

    //check also velocity direction: (when normal well calculated maybe it can be deactivated)
    // Direction = rNode.FastGetSolutionStepValue(VELOCITY);

    // double modulus = norm_2(Direction);
    // if(modulus!=0)
    //    Direction/=modulus;

    // modulus = fabs(inner_prod((MidPoint-rNode.Coordinates()), Direction));

    // if( modulus < distance )
    //   distance = modulus;
    //check also velocity direction: (optional)

    return distance;

    KRATOS_CATCH( "" )
  }

  //**************************************************************************
  //**************************************************************************

  void GetDirectionToEdge(array_1d<double,3>& Direction, Node<3>& rEdgeNodeA, Node<3>& rEdgeNodeB)
  {
    //get wall direction from normals:
     if( rEdgeNodeA.FastGetSolutionStepValue(SHRINK_FACTOR) == 1 && rEdgeNodeB.FastGetSolutionStepValue(SHRINK_FACTOR) == 1 ){
      Direction = 0.5 * (rEdgeNodeA.FastGetSolutionStepValue(NORMAL) + rEdgeNodeB.FastGetSolutionStepValue(NORMAL) ) ;
    }
    else{
      if( rEdgeNodeA.FastGetSolutionStepValue(SHRINK_FACTOR) == 1 ){
        Direction = rEdgeNodeA.FastGetSolutionStepValue(NORMAL);
      }
      else if( rEdgeNodeB.FastGetSolutionStepValue(SHRINK_FACTOR) == 1 ){
        Direction = rEdgeNodeB.FastGetSolutionStepValue(NORMAL);
      }
      else{
        Direction = 0.5 * (rEdgeNodeA.FastGetSolutionStepValue(NORMAL) + rEdgeNodeB.FastGetSolutionStepValue(NORMAL) ) ;
      }
    }

    if( norm_2(Direction) )
      Direction /= norm_2(Direction);

  }


  //**************************************************************************
  //**************************************************************************

  bool EraseCriticalNodes(ModelPart& rModelPart, unsigned int& inside_nodes_removed)
  {

    KRATOS_TRY

    bool any_node_removed = false;
    unsigned int erased_nodes=0;

    for(ModelPart::ElementsContainerType::const_iterator ie = rModelPart.ElementsBegin(); ie != rModelPart.ElementsEnd(); ++ie)
    {
      std::vector<Node<3>::Pointer> EdgeNodes;

      Element::GeometryType& rGeometry = ie->GetGeometry();
      //coordinates
      for(unsigned int i=0; i<rGeometry.size(); ++i)
      {
        if((rGeometry[i].Is(RIGID) || rGeometry[i].Is(SOLID)) && rGeometry[i].IsNot(INLET)){
          EdgeNodes.push_back(rGeometry(i));
        }
      }

      if( EdgeNodes.size() >= 2 ){

        const double dimension = ie->GetGeometry().WorkingSpaceDimension();

        double CriticalVolume = 0;
        double ElementSize    = 0;
        double ElementHeight  = 0;
        double MeanEdgeLength = 0;

        double RigidEdgeLength = std::numeric_limits<double>::min();
        if(EdgeNodes.size() == 2){
          RigidEdgeLength = norm_2(EdgeNodes.front()->Coordinates() - EdgeNodes.back()->Coordinates());
        }

        double VolumeTolerance = 1e-2*pow(4.0*mrRemesh.Refine->CriticalRadius,dimension);

        //set distances
        if( dimension == 2 ){
          ElementSize    = ie->GetGeometry().Area();
          CriticalVolume = 0.5 * mrRemesh.Refine->CriticalRadius * mrRemesh.Refine->CriticalRadius;
          VolumeTolerance *= 1.15e-2;
        }
        else if( dimension == 3 ){
          ElementSize    = ie->GetGeometry().Volume();
          CriticalVolume = 0.5 * mrRemesh.Refine->CriticalRadius * mrRemesh.Refine->CriticalRadius * mrRemesh.Refine->CriticalRadius;
          VolumeTolerance *= 0.25;
        }

        MesherUtilities MesherUtils;
        double MaxRelativeVelocity = 1.5; //arbitrary value, will depend on time step (AF)

        double VolumeChange = 0;

        bool speedy_approach = MesherUtils.CheckRelativeVelocities(rGeometry, MaxRelativeVelocity);
        bool volume_decrease = MesherUtils.CheckVolumeDecrease(rGeometry, dimension, VolumeTolerance, VolumeChange);

        //compare the element height to wall edge length
        for (unsigned int i = 0; i < rGeometry.size(); ++i)
        {

          if( rGeometry[i].IsNot(RIGID) && rGeometry[i].IsNot(SOLID) && rGeometry[i].IsNot(TO_ERASE) && rGeometry[i].IsNot(ISOLATED) && rGeometry[i].IsNot(INLET) ){

            //set distances
            bool approaching_edge = false;
            std::vector<Node<3>::Pointer> SelectedEdgeNodes(2);

            if( EdgeNodes.size() == 2 ){
              SelectedEdgeNodes = EdgeNodes;
              approaching_edge = this->CheckApproachingEdge(rGeometry[i],*EdgeNodes.front(),*EdgeNodes.back());
            }
            else{

              //select larger edge where the node is approaching
              for(unsigned int j=0; j<EdgeNodes.size()-1; ++j)
              {
                for(unsigned int k=j+1; k<EdgeNodes.size(); ++k)
                {
                  double Length = norm_2(EdgeNodes[k]->Coordinates() - EdgeNodes[j]->Coordinates());
                  if( Length > RigidEdgeLength && this->CheckApproachingEdge(rGeometry[i],*EdgeNodes[k],*EdgeNodes[j]) ){
                    RigidEdgeLength = Length;
                    SelectedEdgeNodes[0] = EdgeNodes[k];
                    SelectedEdgeNodes[1] = EdgeNodes[j];
                    approaching_edge = true;
                  }
                }
              }
            }

            // if the node is moving towards to the wall
            if( approaching_edge ){

              ElementHeight  = 0;
              MeanEdgeLength = 0;

              if( dimension == 2 ){
                //ElementHeight  = 2.0 * ElementSize / RigidEdgeLength; (AF)
                //MeanEdgeLength = 0.5 * RigidEdgeLength;
                ElementHeight  = 2.0 * this->GetDistanceToEdge(rGeometry[i],*SelectedEdgeNodes[0],*SelectedEdgeNodes[1]);
                MeanEdgeLength = 0.5 * RigidEdgeLength;
              }
              else if( dimension == 3 ){
                ElementHeight  = 2.0 * this->GetDistanceToEdge(rGeometry[i],*SelectedEdgeNodes[0],*SelectedEdgeNodes[1]);
                MeanEdgeLength = 0.5 * RigidEdgeLength;
              }

              //avoid erasing freesurface particle in touch with wall (penalize length)
              if(rGeometry[i].Is(FREE_SURFACE)){

                WeakPointerVector<Node<3> >& NeighbourNodes = rGeometry[i].GetValue(NEIGHBOUR_NODES);
                unsigned int NumberOfNeighbours = NeighbourNodes.size();
                unsigned int BoundaryWallNodes=0;
                unsigned int FreeSurfaceNodes=0;
                for (WeakPointerVector< Node <3> >::iterator nn = NeighbourNodes.begin();nn != NeighbourNodes.end(); ++nn)
                {
                  if(nn->Is(RIGID) || nn->Is(SOLID)){
                    ++BoundaryWallNodes;
                  }
                  else if(nn->Is(FREE_SURFACE)){
                    ++FreeSurfaceNodes;
                  }
                }
                if( (BoundaryWallNodes+FreeSurfaceNodes) == NumberOfNeighbours && BoundaryWallNodes>0){
                  MeanEdgeLength*=0.5;
                }
              }
              // if the node is very close to the wall or the volume is critically small, it is erased
              if( ElementHeight < MeanEdgeLength && (speedy_approach && volume_decrease) ){

                if( 4 * ElementHeight < MeanEdgeLength ){
                  //std::cout<<" LAYER NODE ERASE DISTANCE ["<<rGeometry[i].Id()<<"] volumes:"<<ElementSize<<" < "<<CriticalVolume<<" heights: "<<ElementHeight<<" < "<<MeanEdgeLength<<" speedy "<<speedy_approach<<" volume_decrease "<<volume_decrease<<std::endl;
                  rGeometry[i].Set(TO_ERASE);
                  ++erased_nodes;
                  ++inside_nodes_removed;
                }
                else{
                  if( rGeometry[i].IsNot(FREE_SURFACE) ){
                    this->MoveNodeToNonRigidMeanPosition(rGeometry[i],*SelectedEdgeNodes[0],*SelectedEdgeNodes[1],MeanEdgeLength);
                  }
                  else{
                    if( !this->MoveBoundaryNodeToMeanPosition(rGeometry[i]) ){
                      rGeometry[i].Set(TO_ERASE);
                      ++erased_nodes;
                      ++inside_nodes_removed;
                    }
                  }
                }

              }
              else if( ElementSize<CriticalVolume ){
                //std::cout<<" LAYER NODE ERASE VOLUME ["<<rGeometry[i].Id()<<"] volumes:"<<ElementSize<<" < "<<CriticalVolume<<" heights: "<<ElementHeight<<" < "<<MeanEdgeLength<<" speedy "<<speedy_approach<<" volume_decrease "<<volume_decrease<<std::endl;
                // if( rGeometry[i].Is(FREE_SURFACE) )
                //   std::cout<<" LAYER NODE FREE_SURFACE ERASED "<<std::endl;
                rGeometry[i].Set(TO_ERASE);
                ++erased_nodes;
                ++inside_nodes_removed;
              }

            }

          }

        }

      }
    }

    if(erased_nodes>0){
      std::cout<<" Layer Nodes removed "<<erased_nodes<<std::endl;
      any_node_removed = true;
    }

    return any_node_removed;

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
