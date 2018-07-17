//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:                   AFranci $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                October 2016 $
//   Revision:            $Revision:                     0.0 $
//
//


#if !defined(KRATOS_INSERT_NEW_NODES_MESHER_PROCESS_H_INCLUDED )
#define  KRATOS_INSERT_NEW_NODES_MESHER_PROCESS_H_INCLUDED


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
//Data:
//StepData: CONTACT_FORCE, DISPLACEMENT
//Flags:    (checked)
//          (set)
//          (modified)
//          (reset)
//(set):=(set in this process)

namespace Kratos
{

///@name Kratos Classes
///@{

/// Refine Mesh Elements Process 2D and 3D
/** The process labels the nodes to be refined (TO_REFINE)
    if the ThresholdVariable  is larger than a ReferenceThreshold
*/

class InsertNewNodesMesherProcess
    : public MesherProcess
{
 public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of Process
  KRATOS_CLASS_POINTER_DEFINITION( InsertNewNodesMesherProcess );

  typedef ModelPart::NodeType                   NodeType;
  typedef ModelPart::ConditionType         ConditionType;
  typedef ModelPart::PropertiesType       PropertiesType;
  typedef ConditionType::GeometryType       GeometryType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  InsertNewNodesMesherProcess(ModelPart& rModelPart,
                              MesherUtilities::MeshingParameters& rRemeshingParameters,
                              int EchoLevel)
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
  {
    mEchoLevel = EchoLevel;
  }


  /// Destructor.
  virtual ~InsertNewNodesMesherProcess() {}


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
          std::cout<<" [ GENERATE NEW NODES for homomgeneous mesh: "<<std::endl;

    if( mrModelPart.Name() != mrRemesh.SubModelPartName )
      std::cout<<" ModelPart Supplied do not corresponds to the Meshing Domain: ("<<mrModelPart.Name()<<" != "<<mrRemesh.SubModelPartName<<")"<<std::endl;

    const ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();
    double currentTime = rCurrentProcessInfo[TIME];
    double timeInterval = rCurrentProcessInfo[DELTA_TIME];
    if(currentTime<2*timeInterval){
      mrRemesh.Info->RemovedNodes=0;
      if( mEchoLevel > 1 )
	std::cout<<" First meshes: I repare the mesh without adding new nodes"<<std::endl;
    }

    int ElementsToRefine=0;
    ElementsToRefine=mrRemesh.Info->RemovedNodes;

    if(ElementsToRefine>0 && mEchoLevel > 1)
      std::cout<<" I will look for "<<ElementsToRefine <<" new nodes"<<std::endl;


    if(ElementsToRefine>0 )
    {
      std::vector<array_1d<double,3> > NewPositions;
      std::vector<double > BiggestVolumes;
      std::vector<array_1d< unsigned int,4 > > NodesIDToInterpolate;
      std::vector<Node<3>::DofsContainerType > NewDofs;

      int CountNodes=0;

      NewPositions.resize(ElementsToRefine);
      BiggestVolumes.resize(ElementsToRefine);
      NodesIDToInterpolate.resize(ElementsToRefine);
      NewDofs.resize(ElementsToRefine);

      for(int nn= 0; nn< ElementsToRefine; ++nn)
      {
        BiggestVolumes[nn]=-1.0;
      }

      ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin();
      // const unsigned int nds = element_begin->GetGeometry().size();
      for(ModelPart::ElementsContainerType::const_iterator ie = element_begin; ie != mrModelPart.ElementsEnd(); ++ie)
      {

        const unsigned int dimension = ie->GetGeometry().WorkingSpaceDimension();

        //////// choose the right (big and safe) elements to refine and compute the new node position and variables ////////
        if(dimension==2){
          SelectEdgeToRefine2D(ie->GetGeometry(),NewPositions,BiggestVolumes,NodesIDToInterpolate,NewDofs,CountNodes,ElementsToRefine);
        } else if(dimension==3){
          SelectEdgeToRefine3D(ie->GetGeometry(),NewPositions,BiggestVolumes,NodesIDToInterpolate,NewDofs,CountNodes,ElementsToRefine);
        }

      }// elements loop


      mrRemesh.Info->RemovedNodes -=ElementsToRefine;
      if(CountNodes<ElementsToRefine){
        mrRemesh.Info->RemovedNodes +=ElementsToRefine-CountNodes;
        NewPositions.resize(CountNodes);
        BiggestVolumes.resize(CountNodes);
        NodesIDToInterpolate.resize(CountNodes);
        NewDofs.resize(CountNodes);
      }


      CreateAndAddNewNodes(NewPositions,NodesIDToInterpolate,NewDofs,ElementsToRefine);

    }//if ElementsToRefine>0

    mrRemesh.InputInitializedFlag=false;


    if( mEchoLevel > 1 )
      std::cout<<"   GENERATE NEW NODES ]; "<<std::endl;

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
    return "InsertNewNodesMesherProcess";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "InsertNewNodesMesherProcess";
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


  void SelectEdgeToRefine2D( Element::GeometryType& Element,
			     std::vector<array_1d<double,3> >& NewPositions,
			     std::vector<double >& BiggestVolumes,
			     std::vector<array_1d< unsigned int,4 > >& NodesIDToInterpolate,
			     std::vector<Node<3>::DofsContainerType >& NewDofs,
			     int &CountNodes,
			     int ElementsToRefine)
  {
    KRATOS_TRY

        const unsigned int nds = Element.size();

    unsigned int rigidNodes=0;
    unsigned int freesurfaceNodes=0;
    unsigned int inletNodes=0;
    bool toEraseNodeFound=false;

    for(unsigned int pn=0; pn<nds; ++pn)
    {
      if(Element[pn].Is(RIGID)){
        rigidNodes++;
      }
      if(Element[pn].Is(TO_ERASE)){
        toEraseNodeFound=true;
      }
      if(Element[pn].Is(FREE_SURFACE)){
        freesurfaceNodes++;
      }
      if(Element[pn].Is(INLET)){
        inletNodes++;
      }
    }


    double  limitEdgeLength=1.4*mrRemesh.Refine->CriticalRadius;
    double safetyCoefficient2D=1.5;
    double penalization=1.0;
    if(rigidNodes>1){
      // penalization=0.7;
      penalization=0.8;
      if(inletNodes>0){
      	penalization=0.9;
      }
    }

    double ElementalVolume =  Element.Area();

    array_1d<double,3> Edges(3,0.0);
    array_1d<unsigned int,3> FirstEdgeNode(3,0);
    array_1d<unsigned int,3> SecondEdgeNode(3,0);
    double WallCharacteristicDistance=0;
    array_1d<double,3> CoorDifference = Element[1].Coordinates() - Element[0].Coordinates();
    // array_1d<double,3> CoorDifference(3,0.0);
    // noalias(CoorDifference) = Element[1].Coordinates() - Element[0].Coordinates();
    // CoorDifference = Element[1].Coordinates() - Element[0].Coordinates();
    double SquaredLength = CoorDifference[0]*CoorDifference[0] + CoorDifference[1]*CoorDifference[1];
    Edges[0]=sqrt(SquaredLength);
    FirstEdgeNode[0]=0;
    SecondEdgeNode[0]=1;
    if(Element[0].Is(RIGID) && Element[1].Is(RIGID)){
      WallCharacteristicDistance=Edges[0];
    }
    unsigned int Counter=0;
    for (unsigned int i = 2; i < nds; ++i){
      for(unsigned int j = 0; j < i; ++j)
      {
        noalias(CoorDifference) = Element[i].Coordinates() - Element[j].Coordinates();
        // CoorDifference = Element[i].Coordinates() - Element[j].Coordinates();
        SquaredLength = CoorDifference[0]*CoorDifference[0] + CoorDifference[1]*CoorDifference[1];
        Counter+=1;
        Edges[Counter]=sqrt(SquaredLength);
        FirstEdgeNode[Counter]=j;
        SecondEdgeNode[Counter]=i;
        if(Element[i].Is(RIGID) && Element[j].Is(RIGID) && Edges[Counter]>WallCharacteristicDistance ){
          WallCharacteristicDistance=Edges[Counter];
        }
      }

    }

    bool dangerousElement=false;
    if(rigidNodes>1){
      for (unsigned int i = 0; i < 3; ++i){
	if((Edges[i]<WallCharacteristicDistance*safetyCoefficient2D && (Element[FirstEdgeNode[i]].Is(RIGID) || Element[SecondEdgeNode[i]].Is(RIGID))) ||
	   (Element[FirstEdgeNode[i]].Is(RIGID) && Element[SecondEdgeNode[i]].Is(RIGID) )){
	  Edges[i]=0;

	}
	// if(Element[FirstEdgeNode[i]].Is(FREE_SURFACE) && Element[SecondEdgeNode[i]].Is(FREE_SURFACE)){
	//   Edges[i]=0;
	//   // Edges[i]*=penalizationFreeSurface;
	// }
	if((Element[FirstEdgeNode[i]].Is(FREE_SURFACE) || Element[FirstEdgeNode[i]].Is(RIGID))  &&
	   (Element[SecondEdgeNode[i]].Is(FREE_SURFACE)|| Element[SecondEdgeNode[i]].Is(RIGID))){
	  Edges[i]=0;
	}
      }

    }
    if((Edges[0]==0 && Edges[1]==0 && Edges[2]==0) || rigidNodes==3){
      dangerousElement=true;
    }

    if(dangerousElement==false && toEraseNodeFound==false){

      // array_1d<double,3> NewPosition(3,0.0);
      unsigned int maxCount=3;
      double LargestEdge=0;

      for(unsigned int i=0; i<3; ++i)
      {
        if(Edges[i]>LargestEdge){
          maxCount=i;
          LargestEdge=Edges[i];
        }
      }

      if(CountNodes<ElementsToRefine && LargestEdge>limitEdgeLength){
	array_1d<double,3> NewPosition=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
	// noalias(NewPosition)=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
	// NewPosition=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
	NodesIDToInterpolate[CountNodes][0]=Element[FirstEdgeNode[maxCount]].GetId();
	NodesIDToInterpolate[CountNodes][1]=Element[SecondEdgeNode[maxCount]].GetId();
	if(Element[SecondEdgeNode[maxCount]].IsNot(RIGID)){
	  NewDofs[CountNodes]=Element[SecondEdgeNode[maxCount]].GetDofs();
	}else if(Element[FirstEdgeNode[maxCount]].IsNot(RIGID)){
	  NewDofs[CountNodes]=Element[FirstEdgeNode[maxCount]].GetDofs();
	}else{
	  std::cout<<"CAUTION! THIS IS A WALL EDGE"<<std::endl;
	}
	BiggestVolumes[CountNodes]=ElementalVolume;
	NewPositions[CountNodes]=NewPosition;
	CountNodes++;
      }else  if (freesurfaceNodes<3 && rigidNodes<3 ){


	ElementalVolume*=penalization;
	for(int nn= 0; nn< ElementsToRefine; ++nn)
        {
          if(ElementalVolume>BiggestVolumes[nn]){


            bool suitableElement=true;
            if(maxCount<3 && LargestEdge>limitEdgeLength){
              array_1d<double,3> NewPosition=(Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
              // noalias(NewPosition)=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
              // NewPosition=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
              for(int j= 0; j< ElementsToRefine; ++j)
              {
                if(NewPositions[j][0]==NewPosition[0] && NewPositions[j][1]==NewPosition[1]){
                  suitableElement=false;
                }
              }
              if(suitableElement==true){
                NodesIDToInterpolate[nn][0]=Element[FirstEdgeNode[maxCount]].GetId();
                NodesIDToInterpolate[nn][1]=Element[SecondEdgeNode[maxCount]].GetId();
                if(Element[SecondEdgeNode[maxCount]].IsNot(RIGID)){
                  NewDofs[nn]=Element[SecondEdgeNode[maxCount]].GetDofs();
                }else if(Element[FirstEdgeNode[maxCount]].IsNot(RIGID)){
                  NewDofs[nn]=Element[FirstEdgeNode[maxCount]].GetDofs();
                }else{
                  std::cout<<"CAUTION! THIS IS A WALL EDGE"<<std::endl;
                }
                BiggestVolumes[nn]=ElementalVolume;
                NewPositions[nn]=NewPosition;
              }

            }

            break;
          }
        }

      }
    }


    KRATOS_CATCH( "" )

        }




  void SelectEdgeToRefine3D( Element::GeometryType& Element,
			     std::vector<array_1d<double,3> >& NewPositions,
			     std::vector<double >& BiggestVolumes,
			     std::vector<array_1d< unsigned int,4 > >& NodesIDToInterpolate,
			     std::vector<Node<3>::DofsContainerType >& NewDofs,
			     int &CountNodes,
			     int ElementsToRefine)
  {
    KRATOS_TRY



    const unsigned int nds = Element.size();

    unsigned int rigidNodes=0;
    unsigned int freesurfaceNodes=0;
    unsigned int inletNodes=0;
    bool toEraseNodeFound=false;

    for(unsigned int pn=0; pn<nds; ++pn)
    {
      if(Element[pn].Is(RIGID)){
        rigidNodes++;
      }
      if(Element[pn].Is(TO_ERASE)){
        toEraseNodeFound=true;
      }
      if(Element[pn].Is(FREE_SURFACE)){
        freesurfaceNodes++;
      }
      if(Element[pn].Is(INLET)){
        inletNodes++;
      }
    }

    double   limitEdgeLength=1.25*mrRemesh.Refine->CriticalRadius;
    double safetyCoefficient3D=1.6;
    double penalization=1.0;
    if(rigidNodes>2){
      penalization=0.7;
      if(inletNodes>0){
	penalization=0.9;
      }
    }
    // if(freesurfaceNodes>2){
    //   penalization=0.6;
    // }

    double ElementalVolume =  Element.Volume();

    array_1d<double,6> Edges(6,0.0);
    array_1d<unsigned int,6> FirstEdgeNode(6,0);
    array_1d<unsigned int,6> SecondEdgeNode(6,0);
    double WallCharacteristicDistance=0;
    array_1d<double,3> CoorDifference  = Element[1].Coordinates() - Element[0].Coordinates();
    // array_1d<double,3> CoorDifference(3,0.0);
    // noalias(CoorDifference) = Element[1].Coordinates() - Element[0].Coordinates();
    // CoorDifference = Element[1].Coordinates() - Element[0].Coordinates();
    double SquaredLength = CoorDifference[0]*CoorDifference[0] + CoorDifference[1]*CoorDifference[1]  + CoorDifference[2]*CoorDifference[2];
    Edges[0]=sqrt(SquaredLength);
    FirstEdgeNode[0]=0;
    SecondEdgeNode[0]=1;
    if(Element[0].Is(RIGID) && Element[1].Is(RIGID)){
      WallCharacteristicDistance=Edges[0];
    }
    unsigned int Counter=0;
    for (unsigned int i = 2; i < nds; ++i){
      for(unsigned int j = 0; j < i; ++j)
      {
        noalias(CoorDifference) = Element[i].Coordinates() - Element[j].Coordinates();
        // CoorDifference = Element[i].Coordinates() - Element[j].Coordinates();
        SquaredLength = CoorDifference[0]*CoorDifference[0] + CoorDifference[1]*CoorDifference[1] + CoorDifference[2]*CoorDifference[2];
        Counter+=1;
        Edges[Counter]=sqrt(SquaredLength);
        FirstEdgeNode[Counter]=j;
        SecondEdgeNode[Counter]=i;
        if(Element[i].Is(RIGID) && Element[j].Is(RIGID) && Edges[Counter]>WallCharacteristicDistance ){
          WallCharacteristicDistance=Edges[Counter];
        }
      }

    }
    //Edges connectivity: Edges[0]=d01, Edges[1]=d20, Edges[2]=d21, Edges[3]=d30, Edges[4]=d31, Edges[5]=d32
    bool dangerousElement=false;
    if(rigidNodes>1){
      for (unsigned int i = 0; i < 6; ++i){
	if((Edges[i]<WallCharacteristicDistance*safetyCoefficient3D && (Element[FirstEdgeNode[i]].Is(RIGID) || Element[SecondEdgeNode[i]].Is(RIGID))) ||
	   (Element[FirstEdgeNode[i]].Is(RIGID) && Element[SecondEdgeNode[i]].Is(RIGID) )){
	  Edges[i]=0;
	}
	// if(Element[FirstEdgeNode[i]].Is(FREE_SURFACE) && Element[SecondEdgeNode[i]].Is(FREE_SURFACE)){
	//   Edges[i]=0;
	// }
	if((Element[FirstEdgeNode[i]].Is(FREE_SURFACE) || Element[FirstEdgeNode[i]].Is(RIGID))  &&
	   (Element[SecondEdgeNode[i]].Is(FREE_SURFACE)|| Element[SecondEdgeNode[i]].Is(RIGID))){
	  Edges[i]=0;
	}
      }

    }else if(rigidNodes==1){
      if(Element[0].Is(RIGID)){
	Edges[0]=0;
	Edges[1]=0;
	Edges[3]=0;
      }
      if(Element[1].Is(RIGID)){
	Edges[0]=0;
	Edges[2]=0;
	Edges[4]=0;
      }
      if(Element[2].Is(RIGID)){
	Edges[1]=0;
	Edges[2]=0;
	Edges[5]=0;
      }
      if(Element[3].Is(RIGID)){
	Edges[3]=0;
	Edges[4]=0;
	Edges[5]=0;
      }
    }

    if((Edges[0]==0 && Edges[1]==0 && Edges[2]==0 && Edges[3]==0 && Edges[4]==0 && Edges[5]==0) || rigidNodes>2){
      dangerousElement=true;
    }

    //just to fill the vector
    if(dangerousElement==false && toEraseNodeFound==false){

      // array_1d<double,3> NewPosition(3,0.0);
      unsigned int maxCount=6;
      double LargestEdge=0;

      for(unsigned int i=0; i<6; ++i)
      {
        if(Edges[i]>LargestEdge){
          maxCount=i;
          LargestEdge=Edges[i];
        }
      }

      if(CountNodes<ElementsToRefine && LargestEdge>limitEdgeLength){
	array_1d<double,3> NewPosition= (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
	// noalias(NewPosition)=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
	// NewPosition=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
	NodesIDToInterpolate[CountNodes][0]=Element[FirstEdgeNode[maxCount]].GetId();
	NodesIDToInterpolate[CountNodes][1]=Element[SecondEdgeNode[maxCount]].GetId();
	if(Element[SecondEdgeNode[maxCount]].IsNot(RIGID)){
	  NewDofs[CountNodes]=Element[SecondEdgeNode[maxCount]].GetDofs();
	}else if(Element[FirstEdgeNode[maxCount]].IsNot(RIGID)){
	  NewDofs[CountNodes]=Element[FirstEdgeNode[maxCount]].GetDofs();
	}else{
	  std::cout<<"CAUTION! THIS IS A WALL EDGE"<<std::endl;
	}
	BiggestVolumes[CountNodes]=ElementalVolume;
	NewPositions[CountNodes]=NewPosition;
	CountNodes++;
      }else  if (freesurfaceNodes<4 && rigidNodes<4){

	ElementalVolume*=penalization;
	for(int nn= 0; nn< ElementsToRefine; ++nn)
        {
          if(ElementalVolume>BiggestVolumes[nn]){


            bool suitableElement=true;

            if(maxCount<6 && LargestEdge>limitEdgeLength){
              array_1d<double,3> NewPosition= (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
              // noalias(NewPosition)=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
              // NewPosition=    (Element[FirstEdgeNode[maxCount]].Coordinates()+Element[SecondEdgeNode[maxCount]].Coordinates())*0.5;
              for(int j= 0; j< ElementsToRefine; ++j)
              {
                if(NewPositions[j][0]==NewPosition[0] && NewPositions[j][1]==NewPosition[1] && NewPositions[j][2]==NewPosition[2]){
                  suitableElement=false; //this is a repeated node, I have already choose this from another element
                }
              }
              if(suitableElement==true){
                NodesIDToInterpolate[nn][0]=Element[FirstEdgeNode[maxCount]].GetId();
                NodesIDToInterpolate[nn][1]=Element[SecondEdgeNode[maxCount]].GetId();
                if(Element[SecondEdgeNode[maxCount]].IsNot(RIGID)){
                  NewDofs[nn]=Element[SecondEdgeNode[maxCount]].GetDofs();
                }else if(Element[FirstEdgeNode[maxCount]].IsNot(RIGID)){
                  NewDofs[nn]=Element[FirstEdgeNode[maxCount]].GetDofs();
                }else{
                  std::cout<<"CAUTION! THIS IS A WALL EDGE"<<std::endl;
                }
                BiggestVolumes[nn]=ElementalVolume;
                NewPositions[nn]=NewPosition;
              }

            }

            break;
          }
        }
      }
    }

    KRATOS_CATCH( "" )

  }


  void CreateAndAddNewNodes(std::vector<array_1d<double,3> >& NewPositions,
			    std::vector<array_1d< unsigned int,4 > >& NodesIDToInterpolate,
			    std::vector<Node<3>::DofsContainerType >& NewDofs,
			    int ElementsToRefine)
  {
    KRATOS_TRY

    const unsigned int dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

    std::vector<Node<3>::Pointer > list_of_new_nodes;
    double NodeIdParent = MesherUtilities::GetMaxNodeId( *(mrModelPart.GetParentModelPart()) );
    double NodeId = MesherUtilities::GetMaxNodeId(mrModelPart);

    unsigned int initial_node_size =NodeIdParent + 1 + ElementsToRefine; //total model part node size

    if(NodeId>NodeIdParent){
      initial_node_size =NodeId + 1 + ElementsToRefine;
      std::cout<<"initial_node_size  "<<initial_node_size<<std::endl;
    }

    //assign data to dofs
    VariablesList& rVariablesList = mrModelPart.GetNodalSolutionStepVariablesList();

    for(unsigned int nn= 0; nn< NewPositions.size(); ++nn)
    {

      unsigned int id = initial_node_size + nn;

      double  x = NewPositions[nn][0];
      double  y = NewPositions[nn][1];
      double  z = 0;
      if(dimension==3)
        z=NewPositions[nn][2];

      //create a new node
      Node<3>::Pointer pnode = mrModelPart.CreateNewNode(id,x,y,z);

      pnode->Set(NEW_ENTITY); //not boundary
      list_of_new_nodes.push_back( pnode );

      if(mrRemesh.InputInitializedFlag){
        mrRemesh.NodalPreIds.push_back( pnode->Id() );
        pnode->SetId(id);
      }

      //giving model part variables list to the node
      pnode->SetSolutionStepVariablesList(&rVariablesList);

      //set buffer size
      pnode->SetBufferSize(mrModelPart.GetBufferSize());

      Node<3>::DofsContainerType& Reference_dofs = NewDofs[nn];

       //generating the dofs
      for(Node<3>::DofsContainerType::iterator iii = Reference_dofs.begin(); iii != Reference_dofs.end(); ++iii)
      {
        Node<3>::DofType& rDof = *iii;
        Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );

        //(p_new_dof)->FreeDof();
      }

      PointsArrayType  PointsArray;
      PointsArray.push_back( mrModelPart.pGetNode(NodesIDToInterpolate[nn][0]) );
      PointsArray.push_back( mrModelPart.pGetNode(NodesIDToInterpolate[nn][1]) );

      Geometry<Node<3> > LineGeometry( PointsArray );

      std::vector<double> ShapeFunctionsN(2);
      std::fill( ShapeFunctionsN.begin(), ShapeFunctionsN.end(), 0.0 );
      ShapeFunctionsN[0] = 0.5;
      ShapeFunctionsN[1] = 0.5;

      MeshDataTransferUtilities DataTransferUtilities;
      DataTransferUtilities.Interpolate2Nodes( LineGeometry, ShapeFunctionsN, rVariablesList, *pnode);

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
  InsertNewNodesMesherProcess& operator=(InsertNewNodesMesherProcess const& rOther);


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
                                  InsertNewNodesMesherProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const InsertNewNodesMesherProcess& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_INSERT_NEW_NODES_MESHER_PROCESS_H_INCLUDED  defined
