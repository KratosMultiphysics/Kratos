//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:          AFranci $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:          July 2018 $
//   Revision:            $Revision:            0.0 $
//
//


#if !defined(KRATOS_INSERT_FLUID_NODES_MESHER_PROCESS_H_INCLUDED )
#define  KRATOS_INSERT_FLUID_NODES_MESHER_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "includes/model_part.h"
#include "custom_utilities/mesh_error_calculation_utilities.hpp"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_processes/mesher_process.hpp"


namespace Kratos
{

///@name Kratos Classes
///@{

/// Refine Mesh Elements Process 2D and 3D
/** Inserts as many nodes as were released in the removal process
    the nodes are inserted in the edges of the largest elements
*/

class InsertFluidNodesMesherProcess
    : public MesherProcess
{
 public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of Process
  KRATOS_CLASS_POINTER_DEFINITION( InsertFluidNodesMesherProcess );

  typedef ModelPart::NodeType                   NodeType;
  typedef ModelPart::ConditionType         ConditionType;
  typedef ModelPart::PropertiesType       PropertiesType;
  typedef ConditionType::GeometryType       GeometryType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  InsertFluidNodesMesherProcess(ModelPart& rModelPart,
                              MesherUtilities::MeshingParameters& rRemeshingParameters,
                              int EchoLevel)
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
  {
    mEchoLevel = EchoLevel;
  }


  /// Destructor.
  virtual ~InsertFluidNodesMesherProcess() {}


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

    if( mEchoLevel > 1 )
      std::cout<<" [ INSERT NEW NODES for homomgeneous mesh: "<<std::endl;

    if( mrModelPart.Name() != mrRemesh.SubModelPartName )
      std::cout<<" ModelPart Supplied do not corresponds to the Meshing Domain: ("<<mrModelPart.Name()<<" != "<<mrRemesh.SubModelPartName<<")"<<std::endl;

    unsigned int ElementsToRefine = mrRemesh.Info->RemovedNodes;

    if( ElementsToRefine > 0 )
    {
      std::vector<array_1d<double,3> >         NewPositions(ElementsToRefine);
      std::vector<array_1d< unsigned int,4 > > NodeIdsToInterpolate(ElementsToRefine);
      std::vector<Node<3>::DofsContainerType > NewDofs(ElementsToRefine);
      std::vector<double >                     LargestVolumes(ElementsToRefine);
      std::fill(LargestVolumes.begin(), LargestVolumes.end(), -1 );

      unsigned int NodesToRefine=0;

      for(ModelPart::ElementsContainerType::const_iterator i_elem = mrModelPart.ElementsBegin(); i_elem != mrModelPart.ElementsEnd(); ++i_elem)
      {

        const unsigned int dimension = i_elem->GetGeometry().WorkingSpaceDimension();

        if(dimension==2){
          SelectEdgeToRefine2D(i_elem->GetGeometry(),NewPositions,LargestVolumes,NodeIdsToInterpolate,NewDofs,NodesToRefine,ElementsToRefine);
        } else if(dimension==3){
          SelectEdgeToRefine3D(i_elem->GetGeometry(),NewPositions,LargestVolumes,NodeIdsToInterpolate,NewDofs,NodesToRefine,ElementsToRefine);
        }

      }

      mrRemesh.Info->InsertedNodes = NodesToRefine;

      // Here only splits as vertices as removed nodes ( without cheking the longest vertices )
      this->CreateAndAddNewNodes(NewPositions,NodeIdsToInterpolate,NewDofs,ElementsToRefine);

    }

    if( mEchoLevel > 1 )
      std::cout<<"   INSERT NEW NODES ]; ("<<mrRemesh.Info->InsertedNodes<<")"<<std::endl;


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
    return "InsertFluidNodesMesherProcess";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "InsertFluidNodesMesherProcess";
  }

  /// Print object's data.s
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
  ///@name Private Operators
  ///@{


  ///@}
  ///@name Private Operations
  ///@{

  //**************************************************************************
  //**************************************************************************

  void FlagsCounter(GeometryType& rGeometry, unsigned int& rigid_nodes,unsigned int& freesurface_nodes,unsigned int& inlet_nodes,unsigned int& nodes_to_split, unsigned int& nodes_to_erase)
  {
    const unsigned int NumberOfNodes = rGeometry.size();
    for(unsigned int i=0; i<NumberOfNodes; ++i)
    {
      if(rGeometry[i].Is(RIGID)){
        ++rigid_nodes;
      }
      if(rGeometry[i].Is(FREE_SURFACE)){
        ++freesurface_nodes;
      }
      if(rGeometry[i].Is(INLET)){
        ++inlet_nodes;
      }
      if(rGeometry[i].Is(TO_SPLIT)){
        ++nodes_to_split;
      }
      if(rGeometry[i].Is(TO_ERASE)){
        ++nodes_to_erase;
      }

    }

  }

  //**************************************************************************
  //**************************************************************************

  void SelectEdgeToRefine2D( Element::GeometryType& rGeometry,
			     std::vector<array_1d<double,3> >& rNewPositions,
			     std::vector<double >& rLargestVolumes,
			     std::vector<array_1d< unsigned int,4 > >& rNodeIdsToInterpolate,
			     std::vector<Node<3>::DofsContainerType >& rNewDofs,
			     unsigned int& rNodesToRefine,
			     const unsigned int& rElementsToRefine)
  {
    KRATOS_TRY

    const unsigned int NumberOfNodes = rGeometry.size();

    unsigned int rigid_nodes=0;
    unsigned int freesurface_nodes=0;
    unsigned int inlet_nodes=0;
    unsigned int nodes_to_split=0;
    unsigned int nodes_to_erase=0;

    this->FlagsCounter(rGeometry, rigid_nodes, freesurface_nodes, inlet_nodes, nodes_to_split, nodes_to_erase);

    bool any_node_to_erase=false;
    if( nodes_to_erase > 0 )
      any_node_to_erase = true;

    double critical_edge_length=5.0*mrRemesh.Refine->CriticalRadius;

    double length_tolerance=1.5;
    double penalization=1.0;

    if(rigid_nodes>1){
      // penalization=0.7;
      penalization=0.8;
      if(inlet_nodes>0){
      	penalization=0.9;
      }
    }

    double ElementalVolume =  rGeometry.Area();

    array_1d<double,3> Edges(3,0.0);
    array_1d<unsigned int,3> FirstEdgeNode(3,0);
    array_1d<unsigned int,3> SecondEdgeNode(3,0);

    double WallCharacteristicDistance=0;
    array_1d<double,3> CoorDifference;
    noalias(CoorDifference) = rGeometry[1].Coordinates() - rGeometry[0].Coordinates();

    double SquaredLength = CoorDifference[0]*CoorDifference[0] + CoorDifference[1]*CoorDifference[1];
    Edges[0]=sqrt(SquaredLength);
    FirstEdgeNode[0]=0;
    SecondEdgeNode[0]=1;
    if(rGeometry[0].Is(RIGID) && rGeometry[1].Is(RIGID)){
      WallCharacteristicDistance=Edges[0];
    }
    unsigned int Counter=0;
    for (unsigned int i = 2; i < NumberOfNodes; ++i)
    {
      for(unsigned int j = 0; j < i; ++j)
      {
        noalias(CoorDifference) = rGeometry[i].Coordinates() - rGeometry[j].Coordinates();
        SquaredLength = CoorDifference[0]*CoorDifference[0] + CoorDifference[1]*CoorDifference[1];
        Counter+=1;
        Edges[Counter]=sqrt(SquaredLength);
        FirstEdgeNode[Counter]=j;
        SecondEdgeNode[Counter]=i;
        if(rGeometry[i].Is(RIGID) && rGeometry[j].Is(RIGID) && Edges[Counter]>WallCharacteristicDistance ){
          WallCharacteristicDistance=Edges[Counter];
        }
      }

    }

    bool dangerousElement=false;

    if(rigid_nodes>1){

      for (unsigned int i = 0; i < 3; ++i)
      {
	if((Edges[i]<WallCharacteristicDistance*length_tolerance && (rGeometry[FirstEdgeNode[i]].Is(RIGID) || rGeometry[SecondEdgeNode[i]].Is(RIGID))) ||
	   (rGeometry[FirstEdgeNode[i]].Is(RIGID) && rGeometry[SecondEdgeNode[i]].Is(RIGID) )){
	  Edges[i]=0;

	}
	if((rGeometry[FirstEdgeNode[i]].Is(FREE_SURFACE) || rGeometry[FirstEdgeNode[i]].Is(RIGID))  &&
	   (rGeometry[SecondEdgeNode[i]].Is(FREE_SURFACE)|| rGeometry[SecondEdgeNode[i]].Is(RIGID))){
	  Edges[i]=0;
	}
      }

    }

    if((Edges[0]==0 && Edges[1]==0 && Edges[2]==0) || rigid_nodes==3){
      dangerousElement=true;
    }

    if(dangerousElement==false && any_node_to_erase==false && nodes_to_split<2){

      unsigned int maxCount=3;
      double LargestEdge=0;

      for(unsigned int i=0; i<3; ++i)
      {
        if(Edges[i]>LargestEdge){
          maxCount=i;
          LargestEdge=Edges[i];
        }
      }

      if(rNodesToRefine<rElementsToRefine && LargestEdge>critical_edge_length){

	array_1d<double,3> NewPosition;
	noalias(NewPosition)= 0.5*(rGeometry[FirstEdgeNode[maxCount]].Coordinates()+rGeometry[SecondEdgeNode[maxCount]].Coordinates());

	rNodeIdsToInterpolate[rNodesToRefine][0]=rGeometry[FirstEdgeNode[maxCount]].GetId();
	rNodeIdsToInterpolate[rNodesToRefine][1]=rGeometry[SecondEdgeNode[maxCount]].GetId();

        rGeometry[FirstEdgeNode[maxCount]].Set(TO_SPLIT);
        rGeometry[SecondEdgeNode[maxCount]].Set(TO_SPLIT);

	if(rGeometry[SecondEdgeNode[maxCount]].IsNot(RIGID)){
	  rNewDofs[rNodesToRefine]=rGeometry[SecondEdgeNode[maxCount]].GetDofs();
	}else if(rGeometry[FirstEdgeNode[maxCount]].IsNot(RIGID)){
	  rNewDofs[rNodesToRefine]=rGeometry[FirstEdgeNode[maxCount]].GetDofs();
	}else{
	  std::cout<<"CAUTION! THIS IS A WALL EDGE"<<std::endl;
	}

	rLargestVolumes[rNodesToRefine]=ElementalVolume;
	rNewPositions[rNodesToRefine]=NewPosition;
	rNodesToRefine++;

      }
      else if( freesurface_nodes<3 && rigid_nodes<3 ){

	ElementalVolume*=penalization;
	for(unsigned int nn= 0; nn< rElementsToRefine; ++nn)
        {
          if(ElementalVolume>rLargestVolumes[nn]){

            bool suitableElement=true;
            if(maxCount<3 && LargestEdge>critical_edge_length){
              array_1d<double,3> NewPosition;
              noalias(NewPosition) = 0.5*(rGeometry[FirstEdgeNode[maxCount]].Coordinates()+rGeometry[SecondEdgeNode[maxCount]].Coordinates());

              for(unsigned int j= 0; j< rElementsToRefine; ++j)
              {
                if(rNewPositions[j][0]==NewPosition[0] && rNewPositions[j][1]==NewPosition[1]){
                  suitableElement=false;
                }
              }

              if(suitableElement==true){
                rNodeIdsToInterpolate[nn][0]=rGeometry[FirstEdgeNode[maxCount]].GetId();
                rNodeIdsToInterpolate[nn][1]=rGeometry[SecondEdgeNode[maxCount]].GetId();

                rGeometry[FirstEdgeNode[maxCount]].Set(TO_SPLIT);
                rGeometry[SecondEdgeNode[maxCount]].Set(TO_SPLIT);

                if(rGeometry[SecondEdgeNode[maxCount]].IsNot(RIGID)){
                  rNewDofs[nn]=rGeometry[SecondEdgeNode[maxCount]].GetDofs();
                }else if(rGeometry[FirstEdgeNode[maxCount]].IsNot(RIGID)){
                  rNewDofs[nn]=rGeometry[FirstEdgeNode[maxCount]].GetDofs();
                }else{
                  std::cout<<"CAUTION! THIS IS A WALL EDGE"<<std::endl;
                }
                rLargestVolumes[nn]=ElementalVolume;
                rNewPositions[nn]=NewPosition;
              }

            }

            break;
          }
        }

      }
    }

    // reset TO_SPLIT
    for(ModelPart::NodesContainerType::const_iterator in = mrModelPart.NodesBegin(); in != mrModelPart.NodesEnd(); ++in)
    {
      if(in->Is(TO_SPLIT))
        in->Set(TO_SPLIT,false);
    }

    KRATOS_CATCH( "" )
  }

  //**************************************************************************
  //**************************************************************************


  void SelectEdgeToRefine3D( GeometryType& rGeometry,
			     std::vector<array_1d<double,3> >& rNewPositions,
			     std::vector<double >& rLargestVolumes,
			     std::vector<array_1d< unsigned int,4 > >& rNodeIdsToInterpolate,
			     std::vector<Node<3>::DofsContainerType >& rNewDofs,
			     unsigned int& rNodesToRefine,
			     const unsigned int& rElementsToRefine)
  {
    KRATOS_TRY

    const unsigned int NumberOfNodes = rGeometry.size();

    unsigned int rigid_nodes=0;
    unsigned int freesurface_nodes=0;
    unsigned int inlet_nodes=0;
    unsigned int nodes_to_split=0;
    unsigned int nodes_to_erase=0;

    this->FlagsCounter(rGeometry, rigid_nodes, freesurface_nodes, inlet_nodes, nodes_to_split, nodes_to_erase);

    bool any_node_to_erase=false;
    if( nodes_to_erase > 0 )
      any_node_to_erase = true;

    double critical_edge_length=5.0*mrRemesh.Refine->CriticalRadius;
    double length_tolerance=1.6;
    double penalization=1.0;

    if(rigid_nodes>2){
      penalization=0.7;
      if(inlet_nodes>0){
	penalization=0.9;
      }
    }


    double ElementalVolume =  rGeometry.Volume();

    array_1d<double,6> Edges(6,0.0);
    array_1d<unsigned int,6> FirstEdgeNode(6,0);
    array_1d<unsigned int,6> SecondEdgeNode(6,0);
    double WallCharacteristicDistance=0;
    array_1d<double,3> CoorDifference;
    noalias(CoorDifference) = rGeometry[1].Coordinates() - rGeometry[0].Coordinates();

    double SquaredLength = CoorDifference[0]*CoorDifference[0] + CoorDifference[1]*CoorDifference[1]  + CoorDifference[2]*CoorDifference[2];
    Edges[0]=sqrt(SquaredLength);
    FirstEdgeNode[0]=0;
    SecondEdgeNode[0]=1;

    if(rGeometry[0].Is(RIGID) && rGeometry[1].Is(RIGID)){
      WallCharacteristicDistance=Edges[0];
    }

    unsigned int Counter=0;
    for (unsigned int i = 2; i < NumberOfNodes; ++i)
    {
      for(unsigned int j = 0; j < i; ++j)
      {
        noalias(CoorDifference) = rGeometry[i].Coordinates() - rGeometry[j].Coordinates();
        SquaredLength = CoorDifference[0]*CoorDifference[0] + CoorDifference[1]*CoorDifference[1] + CoorDifference[2]*CoorDifference[2];
        Counter+=1;
        Edges[Counter]=sqrt(SquaredLength);
        FirstEdgeNode[Counter]=j;
        SecondEdgeNode[Counter]=i;
        if(rGeometry[i].Is(RIGID) && rGeometry[j].Is(RIGID) && Edges[Counter]>WallCharacteristicDistance ){
          WallCharacteristicDistance=Edges[Counter];
        }
      }

    }
    //Edges connectivity: Edges[0]=d01, Edges[1]=d20, Edges[2]=d21, Edges[3]=d30, Edges[4]=d31, Edges[5]=d32
    bool dangerousElement=false;

    if(rigid_nodes>1){

      for (unsigned int i = 0; i < 6; ++i)
      {
	if((Edges[i]<WallCharacteristicDistance*length_tolerance && (rGeometry[FirstEdgeNode[i]].Is(RIGID) || rGeometry[SecondEdgeNode[i]].Is(RIGID))) ||
	   (rGeometry[FirstEdgeNode[i]].Is(RIGID) && rGeometry[SecondEdgeNode[i]].Is(RIGID) )){
	  Edges[i]=0;
	}

	if((rGeometry[FirstEdgeNode[i]].Is(FREE_SURFACE) || rGeometry[FirstEdgeNode[i]].Is(RIGID))  &&
	   (rGeometry[SecondEdgeNode[i]].Is(FREE_SURFACE)|| rGeometry[SecondEdgeNode[i]].Is(RIGID))){
	  Edges[i]=0;
	}
      }

    }
    else if(rigid_nodes==1){

      if(rGeometry[0].Is(RIGID)){
	Edges[0]=0;
	Edges[1]=0;
	Edges[3]=0;
      }
      if(rGeometry[1].Is(RIGID)){
	Edges[0]=0;
	Edges[2]=0;
	Edges[4]=0;
      }
      if(rGeometry[2].Is(RIGID)){
	Edges[1]=0;
	Edges[2]=0;
	Edges[5]=0;
      }
      if(rGeometry[3].Is(RIGID)){
	Edges[3]=0;
	Edges[4]=0;
	Edges[5]=0;
      }
    }

    if((Edges[0]==0 && Edges[1]==0 && Edges[2]==0 && Edges[3]==0 && Edges[4]==0 && Edges[5]==0) || rigid_nodes>2){
      dangerousElement=true;
    }

    //just to fill the vector
    if(dangerousElement==false && any_node_to_erase==false && nodes_to_split<2){

      unsigned int maxCount=6;
      double LargestEdge=0;

      for(unsigned int i=0; i<6; ++i)
      {
        if(Edges[i]>LargestEdge){
          maxCount=i;
          LargestEdge=Edges[i];
        }
      }

      if(rNodesToRefine<rElementsToRefine && LargestEdge>critical_edge_length){
	array_1d<double,3> NewPosition;
        noalias(NewPosition) = 0.5*(rGeometry[FirstEdgeNode[maxCount]].Coordinates()+rGeometry[SecondEdgeNode[maxCount]].Coordinates());

	rNodeIdsToInterpolate[rNodesToRefine][0]=rGeometry[FirstEdgeNode[maxCount]].GetId();
	rNodeIdsToInterpolate[rNodesToRefine][1]=rGeometry[SecondEdgeNode[maxCount]].GetId();

        rGeometry[FirstEdgeNode[maxCount]].Set(TO_SPLIT);
        rGeometry[SecondEdgeNode[maxCount]].Set(TO_SPLIT);

	if(rGeometry[SecondEdgeNode[maxCount]].IsNot(RIGID)){
	  rNewDofs[rNodesToRefine]=rGeometry[SecondEdgeNode[maxCount]].GetDofs();
	}else if(rGeometry[FirstEdgeNode[maxCount]].IsNot(RIGID)){
	  rNewDofs[rNodesToRefine]=rGeometry[FirstEdgeNode[maxCount]].GetDofs();
	}else{
	  std::cout<<"CAUTION! THIS IS A WALL EDGE"<<std::endl;
	}
	rLargestVolumes[rNodesToRefine]=ElementalVolume;
	rNewPositions[rNodesToRefine]=NewPosition;
	rNodesToRefine++;

      }
      else if(freesurface_nodes<4 && rigid_nodes<4){

	ElementalVolume*=penalization;
	for(unsigned int nn= 0; nn< rElementsToRefine; ++nn)
        {
          if(ElementalVolume>rLargestVolumes[nn]){


            bool suitableElement=true;

            if(maxCount<6 && LargestEdge>critical_edge_length){
              array_1d<double,3> NewPosition;
              noalias(NewPosition) = 0.5 * (rGeometry[FirstEdgeNode[maxCount]].Coordinates()+rGeometry[SecondEdgeNode[maxCount]].Coordinates());

              for(unsigned int j= 0; j< rElementsToRefine; ++j)
              {
                if(rNewPositions[j][0]==NewPosition[0] && rNewPositions[j][1]==NewPosition[1] && rNewPositions[j][2]==NewPosition[2]){
                  suitableElement=false; //this is a repeated node, I have already choose this from another element
                }
              }

              if(suitableElement==true){

                rNodeIdsToInterpolate[nn][0]=rGeometry[FirstEdgeNode[maxCount]].GetId();
                rNodeIdsToInterpolate[nn][1]=rGeometry[SecondEdgeNode[maxCount]].GetId();

                rGeometry[FirstEdgeNode[maxCount]].Set(TO_SPLIT);
                rGeometry[SecondEdgeNode[maxCount]].Set(TO_SPLIT);

                if(rGeometry[SecondEdgeNode[maxCount]].IsNot(RIGID)){
                  rNewDofs[nn]=rGeometry[SecondEdgeNode[maxCount]].GetDofs();
                }else if(rGeometry[FirstEdgeNode[maxCount]].IsNot(RIGID)){
                  rNewDofs[nn]=rGeometry[FirstEdgeNode[maxCount]].GetDofs();
                }else{
                  std::cout<<"CAUTION! THIS IS A WALL EDGE"<<std::endl;
                }
                rLargestVolumes[nn]=ElementalVolume;
                rNewPositions[nn]=NewPosition;
              }

            }

            break;
          }
        }
      }
    }

    // reset TO_SPLIT
    for(ModelPart::NodesContainerType::const_iterator in = mrModelPart.NodesBegin(); in != mrModelPart.NodesEnd(); ++in)
    {
      if(in->Is(TO_SPLIT))
        in->Set(TO_SPLIT,false);
    }

    KRATOS_CATCH( "" )
  }


  void CreateAndAddNewNodes(std::vector<array_1d<double,3> >& rNewPositions,
			    std::vector<array_1d< unsigned int,4 > >& rNodeIdsToInterpolate,
			    std::vector<Node<3>::DofsContainerType >& rNewDofs,
			    const unsigned int& rElementsToRefine)
  {
    KRATOS_TRY

    const unsigned int dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

    std::vector<Node<3>::Pointer > list_of_new_nodes;
    const unsigned int NodeIdParent = MesherUtilities::GetMaxNodeId( *(mrModelPart.GetParentModelPart()) );
    const unsigned int NodeId = MesherUtilities::GetMaxNodeId(mrModelPart);

    unsigned int Id = NodeIdParent + 1; //total model part node size

    if(NodeId>NodeIdParent){
      Id = NodeId + 1;
      std::cout<<"initial_node_size  "<<Id<<std::endl;
    }

    //assign data to dofs
    VariablesList& rVariablesList = mrModelPart.GetNodalSolutionStepVariablesList();

    for(unsigned int nn= 0; nn< rNewPositions.size(); ++nn)
    {

      double  x = rNewPositions[nn][0];
      double  y = rNewPositions[nn][1];
      double  z = 0;
      if(dimension==3)
        z=rNewPositions[nn][2];

      //create a new node
      Node<3>::Pointer pnode = mrModelPart.CreateNewNode(Id,x,y,z);
      ++Id;

      //to control the inserted nodes
      // pnode->Set(MODIFIED);
      //std::cout<<" Insert new node "<<pnode->Id()<<std::endl;

      pnode->Set(NEW_ENTITY,true); //not boundary
      list_of_new_nodes.push_back( pnode );

      // commented JMC august 15-18
      // if(mrRemesh.InputInitializedFlag){
      //   mrRemesh.NodalPreIds.push_back( pnode->Id() );
      //   pnode->SetId(Id);
      // }

      //giving model part variables list to the node
      pnode->SetSolutionStepVariablesList(&rVariablesList);

      //set buffer size
      pnode->SetBufferSize(mrModelPart.GetBufferSize());

      Node<3>::DofsContainerType& Reference_dofs = rNewDofs[nn];

      //generating the dofs
      for(Node<3>::DofsContainerType::iterator iii = Reference_dofs.begin(); iii != Reference_dofs.end(); ++iii)
      {
        Node<3>::DofType& rDof = *iii;
        Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );

        //(p_new_dof)->FreeDof();
      }

      Geometry<Node<3> >::PointsArrayType  PointsArray;
      PointsArray.push_back( mrModelPart.pGetNode(rNodeIdsToInterpolate[nn][0]) );
      PointsArray.push_back( mrModelPart.pGetNode(rNodeIdsToInterpolate[nn][1]) );

      Geometry<Node<3> > LineGeometry( PointsArray );

      std::vector<double> ShapeFunctionsN(2);
      std::fill( ShapeFunctionsN.begin(), ShapeFunctionsN.end(), 0.0 );
      ShapeFunctionsN[0] = 0.5;
      ShapeFunctionsN[1] = 0.5;

      MeshDataTransferUtilities DataTransferUtilities;
      DataTransferUtilities.Interpolate2Nodes( LineGeometry, ShapeFunctionsN, rVariablesList, *pnode);

      if( PointsArray[0].Is(FREE_SURFACE) && PointsArray[1].Is(FREE_SURFACE) )
        pnode->Set(FREE_SURFACE);

      if( PointsArray[0].Is(BOUNDARY) && PointsArray[1].Is(BOUNDARY) )
        pnode->Set(BOUNDARY);

    }


    //set the coordinates to the original value
    const array_1d<double,3> ZeroNormal(3,0.0);
    for(std::vector<Node<3>::Pointer>::iterator it =  list_of_new_nodes.begin(); it!=list_of_new_nodes.end(); ++it)
    {
      const array_1d<double,3>& displacement = (*it)->FastGetSolutionStepValue(DISPLACEMENT);
      (*it)->X0() = (*it)->X() - displacement[0];
      (*it)->Y0() = (*it)->Y() - displacement[1];
      (*it)->Z0() = (*it)->Z() - displacement[2];

      (*it)->Set(FLUID);
      (*it)->Set(ACTIVE);

      // std::cout<<" New Node [ "<<(*it)->Id()<<"]: Displacement "<<(*it)->FastGetSolutionStepValue(DISPLACEMENT)<<" Position "<<(*it)->Coordinates()<<std::endl;

      //correct contact_normal interpolation
      if( (*it)->SolutionStepsDataHas(CONTACT_FORCE) )
        noalias((*it)->GetSolutionStepValue(CONTACT_FORCE)) = ZeroNormal;

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
  InsertFluidNodesMesherProcess& operator=(InsertFluidNodesMesherProcess const& rOther);


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
                                  InsertFluidNodesMesherProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const InsertFluidNodesMesherProcess& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_INSERT_FLUID_NODES_MESHER_PROCESS_H_INCLUDED  defined
