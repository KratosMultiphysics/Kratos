//
//   Project Name:        KratosPfemFluidApplication $
//   Created by:          $Author:       JMCarbonell $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:           July 2018 $
//   Revision:            $Revision:             0.0 $
//
//

#if !defined(KRATOS_REMOVE_MESH_NODES_FOR_FLUIDS_PROCESS_H_INCLUDED )
#define  KRATOS_REMOVE_MESH_NODES_FOR_FLUIDS_PROCESS_H_INCLUDED


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
              if(in->IsNot(FREE_SURFACE) && FreeSurfaceNeighbours==dimension){
                this->MoveNodeToMeanPosition((*in));
              }
              else{
                if( !this->CheckEngagedNode((*in),neighbours,neighbour_distances,n_points_in_radius) ){ //we release the node if no other nodes neighbours are being erased
                  in->Set(TO_ERASE);
                  any_node_removed = true;
                  inside_nodes_removed++;
                  //distance_remove++;
                  //std::cout<<"     Distance Criterion Node ["<<in->Id()<<"] TO_ERASE "<<std::endl;
                }
              }
            }

          }
          else if( dimension== 2 ){

            // std::cout<<"  Remove close boundary nodes: Candidate ["<<in->Id()<<"]"<<std::endl;
            bool engaged_node = false;
            unsigned int counter = 0;
            for(std::vector<Node<3>::Pointer>::iterator nn=neighbours.begin(); nn!=neighbours.begin() + n_points_in_radius ; ++nn)
            {

              if ( (*nn)->Is(BOUNDARY) && neighbour_distances[counter] < 1.5 * size_for_distance_boundary && neighbour_distances[counter] > 0.0 )
              {
                if((*nn)->Is(TO_ERASE)){
                  engaged_node = true;
                }
              }

              ++counter;
            }

            if(!engaged_node){ //Can be inserted in the boundary refine
              in->Set(TO_ERASE);
              //any_node_removed = true;
              boundary_nodes_removed++;
              //distance_remove ++;
              std::cout<<"     Removed Boundary Node ["<<in->Id()<<"] on Distance "<<std::endl;
            }

          }


        }
      }
    }

    if( dimension== 2 ){
      if( boundary_nodes_removed > 0 )
        this->MoveBoundaryNodesToMeanPosition(rModelPart,boundary_nodes_removed);

      if( boundary_nodes_removed > 0 )
        any_node_removed = true;
    }

    bool critical_nodes_removed = false;
    if( dimension== 2 )
      critical_nodes_removed = EraseCriticalNodes2D(rModelPart,inside_nodes_removed);

    if( any_node_removed || critical_nodes_removed )
      any_node_removed = true;

    // if( dimension ==3 ){
    //   unsigned int erased_nodes=0;
    //   for(ModelPart::ElementsContainerType::const_iterator ie = rModelPart.ElementsBegin(); ie != rModelPart.ElementsEnd(); ++ie)
    //   {
    //     unsigned int WallNodes=0;
    //     //coordinates
    //     for(unsigned int i=0; i<ie->GetGeometry().size(); ++i)
    //     {
    //       if((ie->GetGeometry()[i].Is(RIGID) && ie->GetGeometry()[i].IsNot(INLET)) || ie->GetGeometry()[i].Is(SOLID)){
    //         WallNodes++;
    //       }
    //     }

    //     if(WallNodes>1)
    //       EraseCriticalNodes3D(rModelPart,ie->GetGeometry(),erased_nodes,inside_nodes_removed);
    //   }
    // }

     //Build boundary after removing boundary nodes due distance criterion
    if( this->mEchoLevel > 1 ){
      std::cout<<"boundary_nodes_removed "<<boundary_nodes_removed<<std::endl;
      std::cout<<"inside_nodes_removed "<<inside_nodes_removed<<std::endl;
    }

    //Build boundary after removing boundary nodes due distance criterion
    if(boundary_nodes_removed){
      std::cout<<"     Rebuild boundary needed after release on Distance "<<std::endl;
      //any_condition_removed = this->RebuildBoundary(rModelPart);
    }

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

  void MoveBoundaryNodesToMeanPosition(ModelPart& rModelPart, unsigned int& boundary_nodes_removed)
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

      if(counter == 1 && id>=0){

        //std::cout<<" Boundary to Move Pre ["<<rGeometry[id].Id()<<"] "<<rGeometry[id].Coordinates()<<std::endl;
        unsigned int FreeSurfaceNodes = 0;
        WeakPointerVector< Node < 3 > >& NeighbourNodes = rGeometry[id].GetValue(NEIGHBOUR_NODES);
        std::vector<array_1d<double,3> > FreeNeighbours;
        for(WeakPointerVector< Node <3> >::iterator nn = NeighbourNodes.begin(); nn != NeighbourNodes.end(); ++nn)
        {
          if(nn->Is(FREE_SURFACE) && nn->IsNot(TO_ERASE)){
            ++FreeSurfaceNodes;
            FreeNeighbours.push_back(nn->Coordinates());
          }
        }

        if( FreeSurfaceNodes == 2 )
        {
          array_1d<double,3> MidPoint  = 0.5 * (FreeNeighbours.front()+FreeNeighbours.back());
          array_1d<double,3> Direction = (FreeNeighbours.front()-FreeNeighbours.back());

          if(norm_2(Direction))
            Direction/=norm_2(Direction);

          array_1d<double,3> Displacement = inner_prod( (MidPoint-rGeometry[id].Coordinates()), Direction ) * Direction;
          noalias(rGeometry[id].Coordinates()) += Displacement;
          noalias(rGeometry[id].FastGetSolutionStepValue(DISPLACEMENT)) += Displacement;
          noalias(rGeometry[id].FastGetSolutionStepValue(DISPLACEMENT,1)) += Displacement;

          rGeometry[id].Set(TO_ERASE,false);
          --boundary_nodes_removed;

        }
        //std::cout<<" Boundary to Move Post ["<<rGeometry[id].Id()<<"] "<<rGeometry[id].Coordinates()<<std::endl;
      }

    }

    KRATOS_CATCH( "" )
  }

  //**************************************************************************
  //**************************************************************************

  void MoveNodeToMeanPosition(Node<3>& rNode)
  {

    KRATOS_TRY

    WeakPointerVector< Node < 3 > >& NeighbourNodes = rNode.GetValue(NEIGHBOUR_NODES);
    unsigned int NumberOfNeighbourNodes = NeighbourNodes.size();

    // std::cout<<" Moved Node Pre ["<<rNode.Id()<<"] Displacement"<<rNode.FastGetSolutionStepValue(DISPLACEMENT)<<" Position "<<rNode.Coordinates()<<" Initial Position "<<rNode.GetInitialPosition()<<std::endl;

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
    rNode.GetInitialPosition() = (rNode.Coordinates() - rNode.FastGetSolutionStepValue(DISPLACEMENT));
    rNode.FastGetSolutionStepValue(VELOCITY)       *= quotient;
    rNode.FastGetSolutionStepValue(VELOCITY,1)     *= quotient;
    rNode.FastGetSolutionStepValue(ACCELERATION)   *= quotient;
    rNode.FastGetSolutionStepValue(ACCELERATION,1) *= quotient;
    rNode.FastGetSolutionStepValue(PRESSURE)            *= quotient;
    rNode.FastGetSolutionStepValue(PRESSURE_VELOCITY)   *= quotient;
    rNode.FastGetSolutionStepValue(PRESSURE_VELOCITY,1) *= quotient;

    // std::cout<<" Moved Node Post ["<<rNode.Id()<<"] Displacement"<<rNode.FastGetSolutionStepValue(DISPLACEMENT)<<" Position "<<rNode.Coordinates()<<" Initial Position "<<rNode.GetInitialPosition()<<std::endl;

    KRATOS_CATCH( "" )
  }

  //**************************************************************************
  //**************************************************************************

  bool EraseCriticalNodes2D(ModelPart& rModelPart, unsigned int& inside_nodes_removed)
  {

    KRATOS_TRY

    bool any_node_removed = false;
    unsigned int erased_nodes=0;

    for(ModelPart::ElementsContainerType::const_iterator ie = rModelPart.ElementsBegin(); ie != rModelPart.ElementsEnd(); ++ie)
    {
      unsigned int WallNodes=0;
      std::vector<Node<3>::Pointer> EdgeNodes;

      Element::GeometryType& rGeometry = ie->GetGeometry();
      //coordinates
      for(unsigned int i=0; i<rGeometry.size(); ++i)
      {
        if((rGeometry[i].Is(RIGID) && rGeometry[i].IsNot(INLET)) || rGeometry[i].Is(SOLID)){
          EdgeNodes.push_back(rGeometry(i));
          ++WallNodes;
        }
      }

      if( WallNodes == 2 ){

        double RigidEdgeLength = norm_2(EdgeNodes.front()->Coordinates() - EdgeNodes.back()->Coordinates());
        double ElementArea     = ie->GetGeometry().Area();

        MesherUtilities MesherUtils;
        const double MaxRelativeVelocity = 1.5; //arbitrary value, will depend on time step (AF)
        bool speedy_approach = MesherUtils.CheckRelativeVelocities(rGeometry, MaxRelativeVelocity);

        const double dimension = ie->GetGeometry().WorkingSpaceDimension();
        bool volume_decrease = MesherUtils.CheckVolumeDecrease(rGeometry, dimension);


        //compare the triangle height to wall edge length
        for (unsigned int i = 0; i < rGeometry.size(); ++i){

          if(rGeometry[i].IsNot(RIGID) && rGeometry[i].IsNot(TO_ERASE) && rGeometry[i].IsNot(SOLID) && rGeometry[i].IsNot(ISOLATED)){

            double ElementHeight  = 2.0 * ElementArea / RigidEdgeLength;
            double MeanEdgeLength = 0.5 * RigidEdgeLength;

            //it is evident when a freesurface particle in touch with wall is erased --> reduce the safety coeff
            if(rGeometry[i].Is(FREE_SURFACE)){

              WeakPointerVector< Node < 3 > >& NeighbourNodes = rGeometry[i].GetValue(NEIGHBOUR_NODES);
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

            // if the node is very close to the wall and moving towards it is erased
            if( (ElementHeight < MeanEdgeLength) && (speedy_approach && volume_decrease) ){
              rGeometry[i].Set(TO_ERASE);
              ++erased_nodes;
              ++inside_nodes_removed;
            }

          }

        }

      }
    }

    if(erased_nodes>0){
      if( this->mEchoLevel > 1 )
	std::cout<<"layer_nodes_removed "<<erased_nodes<<std::endl;
      any_node_removed = true;
    }

    return any_node_removed;

    KRATOS_CATCH( "" )
  }


  //**************************************************************************
  //**************************************************************************

  void EraseCriticalNodes3D(ModelPart& rModelPart, Element::GeometryType& Element, unsigned int &erased_nodes, unsigned int& inside_nodes_removed)
  {

    KRATOS_TRY

    double safetyCoefficient3D=0.6;
    // double safetyCoefficient3D=0.7;


    double elementVolume=Element.Volume();
    double criticalVolume=0.1*this->mrRemesh.Refine->MeanVolume;
    if(elementVolume<criticalVolume){
      for(unsigned int i=0; i<Element.size(); ++i)
      {
        if(Element[i].IsNot(RIGID) && Element[i].IsNot(SOLID) && Element[i].IsNot(TO_ERASE)){
          if( norm_2(Element[i].FastGetSolutionStepValue(VELOCITY)) > 1e-3)
            Element[i].Set(TO_ERASE);
          if( this->mEchoLevel > 1)
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


    // ////////  to compare the non-wall length to wall edge length /////////
    for (unsigned int i = 0; i < Edges.size(); ++i){
      if(((Element[FirstEdgeNode[i]].Is(RIGID) && Element[SecondEdgeNode[i]].IsNot(RIGID)) ||
    	  (Element[SecondEdgeNode[i]].Is(RIGID) && Element[FirstEdgeNode[i]].IsNot(RIGID))) &&
    	 Element[FirstEdgeNode[i]].IsNot(TO_ERASE) &&
    	 Element[SecondEdgeNode[i]].IsNot(TO_ERASE)&&
    	 Edges[i]<safetyCoefficient3D*wallLength){
    	if(Element[FirstEdgeNode[i]].IsNot(RIGID) && Element[FirstEdgeNode[i]].IsNot(SOLID) && Element[FirstEdgeNode[i]].IsNot(TO_ERASE) && Element[FirstEdgeNode[i]].IsNot(ISOLATED)){
          if( norm_2(Element[FirstEdgeNode[i]].FastGetSolutionStepValue(VELOCITY)) > 1e-3)
            Element[FirstEdgeNode[i]].Set(TO_ERASE);
    	  inside_nodes_removed++;
    	  erased_nodes += 1;
    	}else if(Element[SecondEdgeNode[i]].IsNot(RIGID) && Element[SecondEdgeNode[i]].IsNot(SOLID) && Element[SecondEdgeNode[i]].IsNot(TO_ERASE) && Element[SecondEdgeNode[i]].IsNot(ISOLATED)){
          if( norm_2(Element[SecondEdgeNode[i]].FastGetSolutionStepValue(VELOCITY)) > 1e-3)
            Element[SecondEdgeNode[i]].Set(TO_ERASE);
    	  inside_nodes_removed++;
    	  erased_nodes += 1;
    	}

      }

    }

    KRATOS_CATCH( "" )

  }

  //**************************************************************************
  //**************************************************************************

  bool CheckForMovingLayerNodes(ModelPart& rModelPart, Node<3>& CheckedNode, const double wallLength)
  {
    KRATOS_TRY

    const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

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
          InterpolateFromTwoNodes(rModelPart,idMaster,idMaster,idSlave);
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

  //**************************************************************************
  //**************************************************************************


  void InterpolateFromTwoNodes(ModelPart& rModelPart, unsigned int idMaster,unsigned int idSlave1,unsigned int idSlave2)
  {

    KRATOS_TRY

    Node<3>::Pointer MasterNode = rModelPart.pGetNode(idMaster);
    Node<3>::Pointer SlaveNode1 = rModelPart.pGetNode(idSlave1);
    Node<3>::Pointer SlaveNode2 = rModelPart.pGetNode(idSlave2);

    VariablesList& rVariablesList = rModelPart.GetNodalSolutionStepVariablesList();

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

#endif // KRATOS_REMOVE_MESH_NODES_FOR_FLUIDS_PROCESS_H_INCLUDED  defined
