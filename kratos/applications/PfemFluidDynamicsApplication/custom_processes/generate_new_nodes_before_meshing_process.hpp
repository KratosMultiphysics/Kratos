//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:                   AFranci $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                October 2016 $
//   Revision:            $Revision:                     0.0 $
//
//


#if !defined(KRATOS_GENERATE_NEW_NODES_BEFORE_MESHING_PROCESS_H_INCLUDED )
#define  KRATOS_GENERATE_NEW_NODES_BEFORE_MESHING_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes 
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "includes/model_part.h"
#include "custom_utilities/mesh_error_calculation_utilities.hpp"
#include "custom_utilities/modeler_utilities.hpp"

///VARIABLES used:
//Data:      
//StepData: DOMAIN_LABEL, CONTACT_FORCE, DISPLACEMENT
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

class GenerateNewNodesBeforeMeshingProcess
  : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION( GenerateNewNodesBeforeMeshingProcess );

  typedef ModelPart::NodeType                   NodeType;
  typedef ModelPart::ConditionType         ConditionType;
  typedef ModelPart::PropertiesType       PropertiesType;
  typedef ConditionType::GeometryType       GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GenerateNewNodesBeforeMeshingProcess(ModelPart& rModelPart,
					 ModelerUtilities::MeshingParameters& rRemeshingParameters,
					 int EchoLevel) 
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    {
      std::cout<<" generate_new_nodes_process_for_homogeneous_mesh CONSTRUCTOR "<<std::endl;

      mMeshId = mrRemesh.MeshId;
      mEchoLevel = EchoLevel;
    }


    /// Destructor.
    virtual ~GenerateNewNodesBeforeMeshingProcess() {}


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

      if( mEchoLevel > 0 )
	std::cout<<" [ GENERATE NEW NODES for homomgeneous mesh: "<<std::endl;

      if( mrModelPart.Name() != mrRemesh.SubModelPartName )
	std::cout<<" ModelPart Supplied do not corresponds to the Meshing Domain: ("<<mrModelPart.Name()<<" != "<<mrRemesh.SubModelPartName<<")"<<std::endl;

      SelectNewNodesForHomogeneousMeshes();
      std::cout<<"               I will generate "<<mrRemesh.Info->RemovedNodes <<" new nodes"<<std::endl;

      if( mEchoLevel > 0 )
	std::cout<<"   GENERATE NEW NODES ]; "<<std::endl;

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
        return "GenerateNewNodesBeforeMeshingProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GenerateNewNodesBeforeMeshingProcess";
    }

    /// Print object's data.s
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
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{
  
    void ProjectVariablesToNewNodes(std::vector<Node<3>::Pointer >& rlist_of_new_nodes)			    
    {

      KRATOS_TRY

      //defintions for spatial search
      typedef Node<3>                                  PointType;
      typedef Node<3>::Pointer                  PointPointerType;
      typedef std::vector<PointPointerType>          PointVector;
      typedef PointVector::iterator                PointIterator;
      typedef std::vector<double>                 DistanceVector;
      typedef std::vector<double>::iterator     DistanceIterator;
      
      typedef Bucket<3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > BucketType;
      typedef Tree< KDTreePartition<BucketType> >     KdtreeType; //Kdtree
      //defintions for spatial search
      
      unsigned int  bucket_size = 20;
      KdtreeType    NodesTree(rlist_of_new_nodes.begin(),rlist_of_new_nodes.end(),bucket_size);

      //Find out where the new nodes belong to:
      std::vector<double> ShapeFunctionsN;    
      std::vector<VariablesListDataValueContainer> VariablesListVector(rlist_of_new_nodes.size());
      
      VariablesList&  variables_list = mrModelPart.GetNodalSolutionStepVariablesList();
      
      //find the center and "radius" of the element
      double  radius = 0;
      Node<3> center(0,0.0,0.0,0.0);
      
      unsigned int MaximumNumberOfPointsInRadius = rlist_of_new_nodes.size();
      std::vector<Node<3>::Pointer> PointsInRadius (MaximumNumberOfPointsInRadius);
      std::vector<double>  PointsInRadiusDistances (MaximumNumberOfPointsInRadius);
      
      //geometry
      ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin(mMeshId);	  
      const unsigned int nds = element_begin->GetGeometry().size();

      std::vector<std::vector<double> > ElementPointCoordinates(nds);
      std::vector<double> PointCoordinates(3);
      std::fill( PointCoordinates.begin(), PointCoordinates.end(), 0.0 );
      std::fill( ElementPointCoordinates.begin(), ElementPointCoordinates.end(), PointCoordinates );
      
      for(ModelPart::ElementsContainerType::const_iterator ie = mrModelPart.ElementsBegin(mMeshId);
	ie != mrModelPart.ElementsEnd(mMeshId); ie++)
	{

	  //coordinates
	  for(unsigned int i=0; i<ie->GetGeometry().size(); i++)
	    {	
	      PointCoordinates[0] = ie->GetGeometry()[i].X();
	      PointCoordinates[1] = ie->GetGeometry()[i].Y();
	      PointCoordinates[2] = ie->GetGeometry()[i].Z();

	      ElementPointCoordinates[i] = PointCoordinates;
	    }
	    
	  std::fill( PointCoordinates.begin(), PointCoordinates.end(), 0.0 );
	  MeshDataTransferUtilities DataTransferUtilities;
	  DataTransferUtilities.CalculateCenterAndSearchRadius( ElementPointCoordinates, PointCoordinates, radius );

	  //find all of the new nodes within the radius
	  center.X() = PointCoordinates[0];
	  center.Y() = PointCoordinates[1];
	  center.Z() = PointCoordinates[2];

	  double Radius = radius * 1.01;
	  int NumberOfPointsInRadius = NodesTree.SearchInRadius (center, Radius, PointsInRadius.begin(), PointsInRadiusDistances.begin(),  MaximumNumberOfPointsInRadius);


	  //check if inside and eventually interpolate
	  for(std::vector<Node<3>::Pointer>::iterator it_found = PointsInRadius.begin(); it_found != (PointsInRadius.begin() + NumberOfPointsInRadius) ; it_found++)
	    {
		
	      PointCoordinates[0] = (*it_found)->X();
	      PointCoordinates[1] = (*it_found)->Y();
	      PointCoordinates[2] = (*it_found)->Z();

	      bool is_inside = false;
	      is_inside = ModelerUtilities::CalculatePosition( ElementPointCoordinates, PointCoordinates, ShapeFunctionsN );

	      if(is_inside == true)
		{
		  double alpha = 1; //1 to interpolate, 0 to leave the original data
		  DataTransferUtilities.Interpolate( ie->GetGeometry(), ShapeFunctionsN, variables_list, *(it_found), alpha );		    
		}
	    }
	}
      
      KRATOS_CATCH( "" )

    }




  void SelectNewNodesForHomogeneousMeshes()
  {
    KRATOS_TRY

 
      const unsigned int dimension = mrModelPart.ElementsBegin(mMeshId)->GetGeometry().WorkingSpaceDimension();


    int ElementsToRefine=0;
    ElementsToRefine=mrRemesh.Info->RemovedNodes;
    // mrRemesh.Info->RemovedNodes=ElementsToRefine;
    std::cout<<" I will find "<<ElementsToRefine <<" new nodes"<<std::endl;

    std::vector<array_1d<double,3> > NewPositions;
    std::vector<array_1d<double,3> > NewVelocities;
    std::vector<array_1d<double,3> > NewAccelerations;
    std::vector<double > NewPressures;
    std::vector<double > BiggestVolumes;
    std::vector<bool > midPointInterpolation;
    int CountNodes=0;
    double athird=0.3333333333333333;
    NewPositions.resize(ElementsToRefine);
    NewVelocities.resize(ElementsToRefine);
    NewAccelerations.resize(ElementsToRefine);
    NewPressures.resize(ElementsToRefine);
    NewPressures.resize(ElementsToRefine);
    midPointInterpolation.resize(ElementsToRefine);
    BiggestVolumes.resize(ElementsToRefine);
  
    // Node<3>::DofsContainerType& reference_dofs= mrModelPart.Nodes(mMeshId).front().GetDofs();
    Node<3>::DofsContainerType& reference_dofs = (mrModelPart.NodesBegin())->GetDofs();

    bool DofsFound=false;

    if(ElementsToRefine>0 )
      {
	ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin(mMeshId);	  
	ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin(mMeshId);
	const unsigned int nds = element_begin->GetGeometry().size();
	int count=0;
	for(ModelPart::ElementsContainerType::const_iterator ie = element_begin; ie != mrModelPart.ElementsEnd(mMeshId); ie++)
	  {
		
	    const unsigned int dimension = ie->GetGeometry().WorkingSpaceDimension();

	    ModelerUtilities ModelerUtils;
	    double ElementalVolume = 0;	
	    uint rigidNodes=0;
	    uint freesurfaceNodes=0;
	    if(dimension==2)
	      {
		array_1d<double,2> NewPosition;
		array_1d<double,2> NewVelocity;
		array_1d<double,2> NewAcceleration;
		double NewPressure;
		for(uint pn=0; pn<nds; pn++)
		  {
		    if(ie->GetGeometry()[pn].Is(RIGID)){
		      rigidNodes++;
		    }
		  //   else if(DofsFound==false){
		  //   // reference_dofs = ie->GetGeometry()[pn].GetDofs();
		  //   DofsFound=true;
		  // }
		    if(ie->GetGeometry()[pn].Is(FREE_SURFACE))
		      freesurfaceNodes++;
		  }
	  		  
		// std::cout<<"  (2D) the rigid nodes are "<< rigidNodes<<std::endl;

		ElementalVolume =  ie->GetGeometry().Area();

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
		for (uint i = 2; i < nds; i++){
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

		double safetyCoefficient2D=1.5;
		bool dangerousElement=false;
		if(rigidNodes>1){
		  for (uint i = 0; i < 3; i++){
		    if(Edges[i]<WallCharacteristicDistance*safetyCoefficient2D && (ie->GetGeometry()[FirstEdgeNode[i]].Is(RIGID) || ie->GetGeometry()[SecondEdgeNode[i]].Is(RIGID))){
		      Edges[i]=0;
		    }

		  }

		}
		if((Edges[0]==0 && Edges[1]==0 && Edges[2]==0) || rigidNodes==3){
		  dangerousElement=true;
		}

		if(dangerousElement==false){
		  count++;
		  if(count<ElementsToRefine){
		    uint maxCount=3;
		    double LargestEdge=0;
		    for(uint i=0; i<3; i++)
		      {
			if(Edges[i]>LargestEdge){
			  maxCount=i;
			  LargestEdge=Edges[i];
			}
		      }
		    if(maxCount<3){
		    NewPosition=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].Coordinates()                         +ie->GetGeometry()[SecondEdgeNode[maxCount]].Coordinates())*0.5;
		    NewVelocity=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(VELOCITY))*0.5;
		    NewAcceleration=(ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(ACCELERATION))*0.5;
		    NewPressure=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(PRESSURE))*0.5;
		    }else{
		      NewPosition=    (ie->GetGeometry()[0].Coordinates()                         +ie->GetGeometry()[1].Coordinates()                         +ie->GetGeometry()[2].Coordinates())*athird;
		      NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY))*athird;
		      NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION))*athird;
		      NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE))*athird;
		    }

		    BiggestVolumes[CountNodes]=ElementalVolume;
		    NewPositions[CountNodes]=NewPosition;
		    NewVelocities[CountNodes]=NewVelocity;
		    NewAccelerations[CountNodes]=NewAcceleration;
		    NewPressures[CountNodes]=NewPressure;
		    if(freesurfaceNodes<2){
		      midPointInterpolation[CountNodes]=false;
		    }else{
		      midPointInterpolation[CountNodes]=true;
		    }
		    CountNodes++;
		  }else if (freesurfaceNodes<2){
		    for(int nn= 0; nn< ElementsToRefine; nn++)
		      {
			if(ElementalVolume>BiggestVolumes[nn]){
			  uint maxCount=3;
			  double LargestEdge=0;
			  for(uint i=0; i<3; i++)
			    {
			      if(Edges[i]>LargestEdge){
				maxCount=i;
				LargestEdge=Edges[i];
			      }
			    }
			  if(maxCount<3){
			  NewPosition=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].Coordinates()                         +ie->GetGeometry()[SecondEdgeNode[maxCount]].Coordinates())*0.5;
			  NewVelocity=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(VELOCITY))*0.5;
			  NewAcceleration=(ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(ACCELERATION))*0.5;
			  NewPressure=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(PRESSURE))*0.5;
			  }else{
			    NewPosition=    (ie->GetGeometry()[0].Coordinates()                         +ie->GetGeometry()[1].Coordinates()                         +ie->GetGeometry()[2].Coordinates())*athird;
			    NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY))*athird;
			    NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION))*athird;
			    NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE))*athird;
			  }

			  BiggestVolumes[nn]=ElementalVolume;
			  NewPositions[nn]=NewPosition;
			  NewVelocities[nn]=NewVelocity;
			  NewAccelerations[nn]=NewAcceleration;
			  NewPressures[nn]=NewPressure;
			  if(freesurfaceNodes<2){
			    midPointInterpolation[CountNodes]=false;
			  }else{
			    midPointInterpolation[CountNodes]=true;
			  }
			  break;
			}
		      }
		  }
		}
	      }else if(dimension==3){
	      array_1d<double,3> NewPosition;
	      array_1d<double,3> NewVelocity;
	      array_1d<double,3> NewAcceleration;
	      double NewPressure;
	      for(uint pn=0; pn<nds; pn++)
		{
		  if(ie->GetGeometry()[pn].Is(RIGID)){
		    rigidNodes++;
		  }else if(DofsFound==false){
		    reference_dofs = ie->GetGeometry()[pn].GetDofs();
		    DofsFound=true;
		  }
		  if(ie->GetGeometry()[pn].Is(FREE_SURFACE))
		    freesurfaceNodes++;
		}

	      ElementalVolume =  ie->GetGeometry().Volume();

	      // std::cout<<"  (3D) the rigid nodes are "<< rigidNodes<<std::endl;

	      array_1d<double,6> Edges(6,0.0);
	      array_1d<uint,6> FirstEdgeNode(6,0);
	      array_1d<uint,6> SecondEdgeNode(6,0);
	      double WallCharacteristicDistance=0;
	      array_1d<double,3> CoorDifference(3,0.0);
	      CoorDifference = ie->GetGeometry()[1].Coordinates() - ie->GetGeometry()[0].Coordinates();
	      double SquaredLength = CoorDifference[0]*CoorDifference[0] + CoorDifference[1]*CoorDifference[1]  + CoorDifference[2]*CoorDifference[2];
	      Edges[0]=sqrt(SquaredLength);
	      FirstEdgeNode[0]=0;
	      SecondEdgeNode[0]=1;
	      if(ie->GetGeometry()[0].Is(RIGID) && ie->GetGeometry()[1].Is(RIGID)){
		WallCharacteristicDistance=Edges[0];
	      }
	      uint Counter=0;
	      for (uint i = 2; i < nds; i++){
		for(uint j = 0; j < i; j++)
		  {
		    CoorDifference = ie->GetGeometry()[i].Coordinates() - ie->GetGeometry()[j].Coordinates();
		    SquaredLength = CoorDifference[0]*CoorDifference[0] + CoorDifference[1]*CoorDifference[1] + CoorDifference[2]*CoorDifference[2];
		    Counter+=1;
		    Edges[Counter]=sqrt(SquaredLength);
		    FirstEdgeNode[Counter]=j;
		    SecondEdgeNode[Counter]=i;
		    if(ie->GetGeometry()[i].Is(RIGID) && ie->GetGeometry()[j].Is(RIGID) && Edges[Counter]>WallCharacteristicDistance ){
		      WallCharacteristicDistance=Edges[Counter];
		    }
		  }

	      }
	      //Edges connectivity: Edges[0]=d01, Edges[1]=d20, Edges[2]=d21, Edges[3]=d30, Edges[4]=d31, Edges[5]=d32

	      // double distWall=0;
	      double safetyCoefficient3D=1.8;
	      bool dangerousElement=false;
	      if(rigidNodes>1){
		for (uint i = 0; i < 6; i++){
		  if(Edges[i]<WallCharacteristicDistance*safetyCoefficient3D && (ie->GetGeometry()[FirstEdgeNode[i]].Is(RIGID) || ie->GetGeometry()[SecondEdgeNode[i]].Is(RIGID))){
		    Edges[i]=0;
		  }

		}

	      }

	      if(rigidNodes==1){
		if(ie->GetGeometry()[0].Is(RIGID)){
		  Edges[0]=0;
		  Edges[1]=0;
		  Edges[3]=0;
		}
		if(ie->GetGeometry()[1].Is(RIGID)){
		  Edges[0]=0;
		  Edges[2]=0;
		  Edges[4]=0;
		}
		if(ie->GetGeometry()[2].Is(RIGID)){
		  Edges[1]=0;
		  Edges[2]=0;
		  Edges[5]=0;
		}
		if(ie->GetGeometry()[3].Is(RIGID)){
		  Edges[3]=0;
		  Edges[4]=0;
		  Edges[5]=0;
		}
	      }

	      if((Edges[0]==0 && Edges[1]==0 && Edges[2]==0 && Edges[3]==0 && Edges[4]==0 && Edges[5]==0) || rigidNodes>2){
		dangerousElement=true;
	      }

	      //just to fill the vector
	      if(dangerousElement==false){
		count++;
		if(count<ElementsToRefine){
		  uint maxCount=6;
		  double LargestEdge=0;
		  for(uint i=0; i<6; i++)
		    {
		      if(Edges[i]>LargestEdge){
			maxCount=i;
			LargestEdge=Edges[i];
		      }
		    }
		  if(maxCount<6){
		  NewPosition=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].Coordinates()                         +ie->GetGeometry()[SecondEdgeNode[maxCount]].Coordinates())*0.5;
		  NewVelocity=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(VELOCITY))*0.5;
		  NewAcceleration=(ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(ACCELERATION))*0.5;
		  NewPressure=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(PRESSURE))*0.5;
		  }else{
		    NewPosition=    (ie->GetGeometry()[0].Coordinates()+
				     ie->GetGeometry()[1].Coordinates()+
				     ie->GetGeometry()[2].Coordinates()+
				     ie->GetGeometry()[3].Coordinates())*0.25;
		    NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)+
				     ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY)+
				     ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY)+
				     ie->GetGeometry()[3].FastGetSolutionStepValue(VELOCITY))*0.25;
		    NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+
				     ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION)+
				     ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION)+
				     ie->GetGeometry()[3].FastGetSolutionStepValue(ACCELERATION))*0.25;
		    NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)+
				     ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)+
				     ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)+
				     ie->GetGeometry()[3].FastGetSolutionStepValue(PRESSURE))*0.25;
		  }


		  BiggestVolumes[CountNodes]=ElementalVolume;
		  NewPositions[CountNodes]=NewPosition;
		  NewVelocities[CountNodes]=NewVelocity;
		  NewAccelerations[CountNodes]=NewAcceleration;
		  NewPressures[CountNodes]=NewPressure;
		  if(freesurfaceNodes<2){
		    midPointInterpolation[CountNodes]=false;
		  }else{
		    midPointInterpolation[CountNodes]=true;
		  }
		  CountNodes++;
		}else if (freesurfaceNodes<3){
		  for(int nn= 0; nn< ElementsToRefine; nn++)
		    {
		      if(ElementalVolume>BiggestVolumes[nn]){
			uint maxCount=6;
			double LargestEdge=0;
			for(uint i=0; i<6; i++)
			  {
			    if(Edges[i]>LargestEdge){
			      maxCount=i;
			      LargestEdge=Edges[i];
			    }
			  }
			if(maxCount<6){
			  NewPosition=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].Coordinates()                         +ie->GetGeometry()[SecondEdgeNode[maxCount]].Coordinates())*0.5;
			  NewVelocity=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(VELOCITY))*0.5;
			  NewAcceleration=(ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(ACCELERATION))*0.5;
			  NewPressure=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(PRESSURE))*0.5;
			}else{
			  NewPosition=    (ie->GetGeometry()[0].Coordinates()+
					   ie->GetGeometry()[1].Coordinates()+
					   ie->GetGeometry()[2].Coordinates()+
					   ie->GetGeometry()[3].Coordinates())*0.25;
			  NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)+
					   ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY)+
					   ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY)+
					   ie->GetGeometry()[3].FastGetSolutionStepValue(VELOCITY))*0.25;
			  NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+
					   ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION)+
					   ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION)+
					   ie->GetGeometry()[3].FastGetSolutionStepValue(ACCELERATION))*0.25;
			  NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)+
					   ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)+
					   ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)+
					   ie->GetGeometry()[3].FastGetSolutionStepValue(PRESSURE))*0.25;
			}

		
			BiggestVolumes[nn]=ElementalVolume;
			NewPositions[nn]=NewPosition;
			NewVelocities[nn]=NewVelocity;
			NewAccelerations[nn]=NewAcceleration;
			NewPressures[nn]=NewPressure;
			if(freesurfaceNodes<2){
			  midPointInterpolation[CountNodes]=false;
			}else{
			  midPointInterpolation[CountNodes]=true;
			}
			break;
		      }
		    }
		}
	      }
	    }


    	  }	 
      }

    std::vector<Node<3>::Pointer > list_of_new_nodes;
    std::vector<NodeType::Pointer> list_of_NEW_nodes;
    // int idNode=mrModelPart.NumberOfNodes(mMeshId);

    unsigned int initial_node_size = mrModelPart.Nodes().size()+1+ElementsToRefine; //total model part node size

    // int max_id =  mrModelPart.Nodes().back().Id();

    // if(initial_node_size>max_id){
    //   std::cout<<"ATTENTION!!! INITIAL_NODE_SIZE > MAX_ID !!!"<<std::endl;
    //   initial_node_size=max_id; 
    // }
    unsigned int ID = ModelerUtilities::GetMaxNodeId(mrModelPart) + 1;
    std::cout<<"ID ="<<initial_node_size;
    // initial_node_size =ID;
    std::cout<<"              other ID ="<<ID<<std::endl;
    
    //assign data to dofs
    VariablesList& VariablesList = mrModelPart.GetNodalSolutionStepVariablesList();

    for(uint nn= 0; nn< NewPositions.size(); nn++)
      {

    	unsigned int id = initial_node_size + nn ;

        std::cout<<"NEW POINT ("<<id<<") COORDINATES x,y: "<<NewPositions[nn][0]<<" "<<NewPositions[nn][1];
        // std::cout<<"NEW INTERPOLATED VELOCITY: "<<NewVelocities[nn][0]<<" "<<NewVelocities[nn][1]<<" "<<NewVelocities[nn][2]<<std::endl;
        // std::cout<<"NEW INTERPOLATED ACCEL: "<<NewAccelerations[nn][0]<<" "<<NewAccelerations[nn][1]<<" "<<NewAccelerations[nn][2];
        // std::cout<<"NEW INTERPOLATED PRES: "<<NewPressures[nn];
        std::cout<<" ELEMENT VOLUME: "<<BiggestVolumes[nn]<<std::endl;


    	double  x = NewPositions[nn][0];
    	double  y = NewPositions[nn][1];
    	double  z = 0; 
    	if(dimension==3)
    	  z=NewPositions[nn][2];

	
    	Node<3>::Pointer pnode = mrModelPart.CreateNewNode(id,x,y,z);
    	pnode->Set(NEW_ENTITY); //not boundary
    	// pnode->SetBufferSize(mrModelPart.NodesBegin(mMeshId)->GetBufferSize() );
	list_of_new_nodes.push_back( pnode );
	if(mrRemesh.InputInitializedFlag){
	  mrRemesh.NodalPreIds.push_back( pnode->Id() );
	  //pnode->SetId(mrRemesh.NodalPreIds.size()-1);
	  pnode->SetId(nn+1);
	}

	// //giving model part variables list to the node
	pnode->SetSolutionStepVariablesList(&VariablesList);
	      
	// //set buffer size
	pnode->SetBufferSize(mrModelPart.GetBufferSize());

    	if(mMeshId!=0)
    	  mrModelPart.AddNode(pnode,mMeshId);


	//generating the dofs
	// Node<3>::DofsContainerType& reference_dofs = mrModelPart.Nodes(mMeshId).front().GetDofs();
	if(DofsFound==false){
	  reference_dofs = (mrModelPart.NodesBegin())->GetDofs();
	}
	// std::cout<<" node.front "<< mrModelPart.Nodes(mMeshId).front() <<std::endl;
	for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
	  {
	    Node<3>::DofType& rDof = *iii;
	    Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
	    (p_new_dof)->FreeDof();
	  }
	


	// //////////////////////NEW PART ///////////////////////////////////////////////////
	// NodeType::Pointer pNode = boost::make_shared< NodeType >( ID, x, y, z );
	// //set new id
	// if(mrRemesh.InputInitializedFlag){
	//   mrRemesh.NodalPreIds.push_back( ID );
	//   pNode->SetId(nn+1);
	//   if( ID > mrRemesh.MaxNodeIdNumber )
	//     mrRemesh.MaxNodeIdNumber = ID;		
	// }
      
	// //giving model part variables list to the node
	// pNode->SetSolutionStepVariablesList(&VariablesList);
	      
	// //set buffer size
	// pNode->SetBufferSize(mrModelPart.GetBufferSize());

	// //generating the dofs
	// for(Node<3>::DofsContainerType::iterator i_dof = ReferenceDofs.begin(); i_dof != ReferenceDofs.end(); i_dof++)
	//   {
	//     NodeType::DofType& rDof = *i_dof;
	//     NodeType::DofType::Pointer pNewDof = pNode->pAddDof( rDof );

	//     (pNewDof)->FreeDof();
	//   }
	            	          
	// list_of_NEW_nodes.push_back( pNode );

	// ID++;
	// //////////////////////NEW PART ///////////////////////////////////////////////////

      
      }


     if( list_of_new_nodes.size() > 0)
	ProjectVariablesToNewNodes( list_of_new_nodes );

     //set the coordinates to the original value
      const array_1d<double,3> ZeroNormal(3,0.0);
      uint count=0;
      for(std::vector<Node<3>::Pointer>::iterator it =  list_of_new_nodes.begin(); it!=list_of_new_nodes.end(); it++)
	{
	  const array_1d<double,3>& displacement = (*it)->FastGetSolutionStepValue(DISPLACEMENT);
	  (*it)->X0() = (*it)->X() - displacement[0];
	  (*it)->Y0() = (*it)->Y() - displacement[1];
	  (*it)->Z0() = (*it)->Z() - displacement[2];

	  // midPointInterpolation[count]=true;
	  if(midPointInterpolation[count]==true){
	    (*it)->FastGetSolutionStepValue(VELOCITY)=NewVelocities[count];
	    (*it)->FastGetSolutionStepValue(ACCELERATION)=NewAccelerations[count];
	    (*it)->FastGetSolutionStepValue(PRESSURE)=NewPressures[count];
	    array_1d<double,3> NewInterpolatedVelocity=(*it)->FastGetSolutionStepValue(VELOCITY);
	    array_1d<double,3> NewInterpolatedAcc=(*it)->FastGetSolutionStepValue(ACCELERATION);
	    double NewInterpolatedPre=(*it)->FastGetSolutionStepValue(PRESSURE);
	    std::cout<<"..... CHECK NEW INTERPOLATED VELOCITY: "<<NewInterpolatedVelocity[0]<<" "<<NewInterpolatedVelocity[1]<<" "<<NewInterpolatedVelocity[2];
	    std::cout<<" NEW INTERPOLATED ACCEL: "<<NewInterpolatedAcc[0]<<" "<<NewInterpolatedAcc[1]<<" "<<NewInterpolatedAcc[2];
	    std::cout<<" NEW pressure: "<<NewInterpolatedPre<<std::endl;
	  }
	  count++;



	  double density=(*it)->FastGetSolutionStepValue(DENSITY);
	  if(density==0){
	    std::cout<<"...........................  the density of this node was 0 now is 1000 "<<std::endl;
	    (*it)->FastGetSolutionStepValue(DENSITY)=1000.0;
	  }

	//correct contact_normal interpolation
	  if( (*it)->SolutionStepsDataHas(CONTACT_FORCE) )
	    noalias((*it)->GetSolutionStepValue(CONTACT_FORCE)) = ZeroNormal;
		    
	  (*it)->SetValue(DOMAIN_LABEL,mMeshId);
	  
	}

      mrRemesh.InputInitializedFlag=false;


    //mModelerUtilities.SetNodes(mrModelPart,mrRemesh,mMeshId);

    // if(!mrRemesh.InputInitializedFlag){
    //   mrRemesh.InputInitializedFlag=true;
    // }

    KRATOS_CATCH( "" )

      }



  ///////////////////////////// INITIAL BASIC VERSION ////////////////////
  // void SelectNewNodesForHomogeneousMeshes()

  // {
  //   KRATOS_TRY

  //     ////////////////////////////////////////////////  new part  /////////////////////////////////////////////////
 
  //     const unsigned int dimension = mrModelPart.ElementsBegin(mMeshId)->GetGeometry().WorkingSpaceDimension();


  //   int ElementsToRefine=0;
  //   ElementsToRefine=mrRemesh.Info->RemovedNodes;
  //   // mrRemesh.Info->RemovedNodes=ElementsToRefine;
  //   std::cout<<" I will find "<<ElementsToRefine <<" new nodes"<<std::endl;

  //   std::vector<array_1d<double,3> > NewPositions;
  //   std::vector<array_1d<double,3> > NewVelocities;
  //   std::vector<array_1d<double,3> > NewAccelerations;
  //   std::vector<double > NewPressures;
  //   std::vector<double > BiggestVolumes;
  //   std::vector<bool > midPointInterpolation;
  //   int CountNodes=0;
  //   double athird=0.3333333333333333;
  //   NewPositions.resize(ElementsToRefine);
  //   NewVelocities.resize(ElementsToRefine);
  //   NewAccelerations.resize(ElementsToRefine);
  //   NewPressures.resize(ElementsToRefine);
  //   NewPressures.resize(ElementsToRefine);
  //   midPointInterpolation.resize(ElementsToRefine);
  //   BiggestVolumes.resize(ElementsToRefine);
  
  //   if(ElementsToRefine>0 )
  //     {
  // 	// int& OutNumberOfElements  = mrRemesh.OutMesh.GetNumberOfElements();
  // 	// int* OutElementList       = mrRemesh.OutMesh.GetElementList();

  // 	ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin(mMeshId);	  
  // 	ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin(mMeshId);
  // 	const unsigned int nds = element_begin->GetGeometry().size();
  // 	int count=0;
  // 	for(ModelPart::ElementsContainerType::const_iterator ie = element_begin; ie != mrModelPart.ElementsEnd(mMeshId); ie++)
  // 	  {
		
  // 	    const unsigned int dimension = ie->GetGeometry().WorkingSpaceDimension();

  // 	    ModelerUtilities ModelerUtils;
  // 	    double ElementalVolume = 0;	
  // 	    uint rigidNodes=0;
  // 	    uint freesurfaceNodes=0;
  // 	    if(dimension==2)
  // 	      {
  // 		array_1d<double,2> NewPosition;
  // 		array_1d<double,2> NewVelocity;
  // 		array_1d<double,2> NewAcceleration;
  // 		double NewPressure;
  // 		for(uint pn=0; pn<nds; pn++)
  // 		  {
  // 		    if(ie->GetGeometry()[pn].Is(RIGID))
  // 		      rigidNodes++;
  // 		    if(ie->GetGeometry()[pn].Is(FREE_SURFACE))
  // 		      freesurfaceNodes++;
  // 		  }
	  		  
		
  // 		ElementalVolume =  ie->GetGeometry().Area();
	
  // 		array_1d<double,3> Edges(3,0.0);
  // 		array_1d<uint,3> FirstEdgeNode(3,0);
  // 		array_1d<uint,3> SecondEdgeNode(3,0);
  // 		double WallCharacteristicDistance=0;
  // 		array_1d<double,3> CoorDifference(3,0.0);
  // 		CoorDifference = ie->GetGeometry()[1].Coordinates() - ie->GetGeometry()[0].Coordinates();
  // 		double SquaredLength = CoorDifference[0]*CoorDifference[0] + CoorDifference[1]*CoorDifference[1];
  // 		Edges[0]=sqrt(SquaredLength);
  // 		FirstEdgeNode[0]=0;
  // 		SecondEdgeNode[0]=1;
  // 		if(ie->GetGeometry()[0].Is(RIGID) && ie->GetGeometry()[1].Is(RIGID)){
  // 		  WallCharacteristicDistance=Edges[0];
  // 		}
  // 		uint Counter=0;
  // 		for (uint i = 2; i < nds; i++){
  // 		  for(uint j = 0; j < i; j++)
  // 		    {
  // 		      CoorDifference = ie->GetGeometry()[i].Coordinates() - ie->GetGeometry()[j].Coordinates();
  // 		      SquaredLength = CoorDifference[0]*CoorDifference[0] + CoorDifference[1]*CoorDifference[1];
  // 		      Counter+=1;
  // 		      Edges[Counter]=sqrt(SquaredLength);
  // 		      FirstEdgeNode[Counter]=j;
  // 		      SecondEdgeNode[Counter]=i;
  // 		      if(ie->GetGeometry()[i].Is(RIGID) && ie->GetGeometry()[j].Is(RIGID) && Edges[Counter]>WallCharacteristicDistance ){
  // 			WallCharacteristicDistance=Edges[Counter];
  // 		      }
  // 		    }

  // 		}


  // 		double safetyCoefficient2D=1.5;
  // 		bool dangerousElement=false;
  // 		if(rigidNodes>1){

  // 		    if(Edges[0]<WallCharacteristicDistance*safetyCoefficient2D && (ie->GetGeometry()[0].Is(RIGID) || ie->GetGeometry()[1].Is(RIGID))){
  // 		      Edges[0]=0;
  // 		    }
  // 		    if(Edges[1]<WallCharacteristicDistance*safetyCoefficient2D && (ie->GetGeometry()[0].Is(RIGID) || ie->GetGeometry()[2].Is(RIGID))){
  // 		      Edges[1]=0;
  // 		    }
  // 		    if(Edges[2]<WallCharacteristicDistance*safetyCoefficient2D && (ie->GetGeometry()[1].Is(RIGID) || ie->GetGeometry()[2].Is(RIGID))){
  // 		      Edges[2]=0;
  // 		    }
  // 		  }

  // 		}
  // 		if((Edges[0]==0 && Edges[1]==0 && Edges[2]==0) || rigidNodes==3){
  // 		  dangerousElement=true;
  // 		}

  // 		if(dangerousElement==false){
  // 		  count++;
  // 		  if(count<ElementsToRefine){

  // 		    if(Edges[0]>Edges[1] && Edges[0]>Edges[2]){
  // 		      NewPosition=    (ie->GetGeometry()[0].Coordinates()                         +ie->GetGeometry()[1].Coordinates())*0.5;
  // 		      NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY))*0.5;
  // 		      NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION))*0.5;
  // 		      NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE))*0.5;
  // 		    }else if(Edges[1]>Edges[0] && Edges[1]>Edges[2]){
  // 		      NewPosition=    (ie->GetGeometry()[0].Coordinates()                         +ie->GetGeometry()[2].Coordinates())*0.5;
  // 		      NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY))*0.5;
  // 		      NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION))*0.5;
  // 		      NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE))*0.5;
  // 		    }else if(Edges[2]>Edges[0] && Edges[2]>Edges[1]){
  // 		      NewPosition=    (ie->GetGeometry()[1].Coordinates()                         +ie->GetGeometry()[2].Coordinates())*0.5;
  // 		      NewVelocity=    (ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY))*0.5;
  // 		      NewAcceleration=(ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION))*0.5;
  // 		      NewPressure=    (ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE))*0.5;
  // 		    }else{
  // 		      NewPosition=    (ie->GetGeometry()[0].Coordinates()                         +ie->GetGeometry()[1].Coordinates()                         +ie->GetGeometry()[2].Coordinates())*athird;
  // 		      NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY))*athird;
  // 		      NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION))*athird;
  // 		      NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE))*athird;
  // 		    }

  // 		    BiggestVolumes[CountNodes]=ElementalVolume;
  // 		    NewPositions[CountNodes]=NewPosition;
  // 		    NewVelocities[CountNodes]=NewVelocity;
  // 		    NewAccelerations[CountNodes]=NewAcceleration;
  // 		    NewPressures[CountNodes]=NewPressure;
  // 		    if(freesurfaceNodes<2){
  // 		      midPointInterpolation[CountNodes]=false;
  // 		    }else{
  // 		      midPointInterpolation[CountNodes]=true;
  // 		    }
  // 		    CountNodes++;
  // 		  }else if (freesurfaceNodes<2){
  // 		    for(int nn= 0; nn< ElementsToRefine; nn++)
  // 		      {
  // 			if(ElementalVolume>BiggestVolumes[nn]){
  // 			  if(Edges[0]>Edges[1] && Edges[0]>Edges[2]){
  // 			    NewPosition=    (ie->GetGeometry()[0].Coordinates()                         +ie->GetGeometry()[1].Coordinates())*0.5;
  // 			    NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY))*0.5;
  // 			    NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION))*0.5;
  // 			    NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE))*0.5;
  // 			  }else if(Edges[1]>Edges[0] && Edges[1]>Edges[2]){
  // 			    NewPosition=    (ie->GetGeometry()[0].Coordinates()                         +ie->GetGeometry()[2].Coordinates())*0.5;
  // 			    NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY))*0.5;
  // 			    NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION))*0.5;
  // 			    NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE))*0.5;
  //      			  }else if(Edges[2]>Edges[0] && Edges[2]>Edges[1]){
  // 			    NewPosition=    (ie->GetGeometry()[1].Coordinates()                         +ie->GetGeometry()[2].Coordinates())*0.5;
  // 			    NewVelocity=    (ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY))*0.5;
  // 			    NewAcceleration=(ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION))*0.5;
  // 			    NewPressure=    (ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE))*0.5;
  // 			  }else{
  // 			    NewPosition=    (ie->GetGeometry()[0].Coordinates()                         +ie->GetGeometry()[1].Coordinates()                         +ie->GetGeometry()[2].Coordinates())*athird;
  // 			    NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY))*athird;
  // 			    NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION))*athird;
  // 			    NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE))*athird;
  // 			  }

  // 			  BiggestVolumes[nn]=ElementalVolume;
  // 			  NewPositions[nn]=NewPosition;
  // 			  NewVelocities[nn]=NewVelocity;
  // 			  NewAccelerations[nn]=NewAcceleration;
  // 			  NewPressures[nn]=NewPressure;
  // 			  if(freesurfaceNodes<2){
  // 			    midPointInterpolation[CountNodes]=false;
  // 			  }else{
  // 			    midPointInterpolation[CountNodes]=true;
  // 			  }
  // 			  break;
  // 			}
  // 		      }
  // 		  }
  // 		}
  // 	      }else if(dimension==3){
  // 	      array_1d<double,3> NewPosition;
  // 	      array_1d<double,3> NewVelocity;
  // 	      array_1d<double,3> NewAcceleration;
  // 	      double NewPressure;
  // 	      for(uint pn=0; pn<nds; pn++)
  // 		{
  // 		  if(ie->GetGeometry()[pn].Is(RIGID))
  // 		    rigidNodes++;
  // 		  if(ie->GetGeometry()[pn].Is(FREE_SURFACE))
  // 		    freesurfaceNodes++;
  // 		}

  // 	      ElementalVolume =  ie->GetGeometry().Volume();
  // 	      // std::cout<<"ElementalVolume is"<<ElementalVolume<<std::endl;
  // 	      double dist01=0;
  // 	      double dist02=0;
  // 	      double dist03=0;
  // 	      double dist12=0;
  // 	      double dist13=0;
  // 	      double dist23=0;
  // 	      dist01=sqrt((ie->GetGeometry()[0].X()-ie->GetGeometry()[1].X())*(ie->GetGeometry()[0].X()-ie->GetGeometry()[1].X())+
  // 			  (ie->GetGeometry()[0].Y()-ie->GetGeometry()[1].Y())*(ie->GetGeometry()[0].Y()-ie->GetGeometry()[1].Y())+
  // 			  (ie->GetGeometry()[0].Z()-ie->GetGeometry()[1].Z())*(ie->GetGeometry()[0].Z()-ie->GetGeometry()[1].Z()));

  // 	      dist02=sqrt((ie->GetGeometry()[0].X()-ie->GetGeometry()[2].X())*(ie->GetGeometry()[0].X()-ie->GetGeometry()[2].X())+
  // 			  (ie->GetGeometry()[0].Y()-ie->GetGeometry()[2].Y())*(ie->GetGeometry()[0].Y()-ie->GetGeometry()[2].Y())+
  // 			  (ie->GetGeometry()[0].Z()-ie->GetGeometry()[2].Z())*(ie->GetGeometry()[0].Z()-ie->GetGeometry()[2].Z()));

  // 	      dist03=sqrt((ie->GetGeometry()[0].X()-ie->GetGeometry()[3].X())*(ie->GetGeometry()[0].X()-ie->GetGeometry()[3].X())+
  // 			  (ie->GetGeometry()[0].Y()-ie->GetGeometry()[3].Y())*(ie->GetGeometry()[0].Y()-ie->GetGeometry()[3].Y())+
  // 			  (ie->GetGeometry()[0].Z()-ie->GetGeometry()[3].Z())*(ie->GetGeometry()[0].Z()-ie->GetGeometry()[3].Z()));

  // 	      dist12=sqrt((ie->GetGeometry()[1].X()-ie->GetGeometry()[2].X())*(ie->GetGeometry()[1].X()-ie->GetGeometry()[2].X())+
  // 			  (ie->GetGeometry()[1].Y()-ie->GetGeometry()[2].Y())*(ie->GetGeometry()[1].Y()-ie->GetGeometry()[2].Y())+
  // 			  (ie->GetGeometry()[1].Z()-ie->GetGeometry()[2].Z())*(ie->GetGeometry()[1].Z()-ie->GetGeometry()[2].Z()));

  // 	      dist13=sqrt((ie->GetGeometry()[1].X()-ie->GetGeometry()[3].X())*(ie->GetGeometry()[1].X()-ie->GetGeometry()[3].X())+
  // 			  (ie->GetGeometry()[1].Y()-ie->GetGeometry()[3].Y())*(ie->GetGeometry()[1].Y()-ie->GetGeometry()[3].Y())+
  // 			  (ie->GetGeometry()[1].Z()-ie->GetGeometry()[3].Z())*(ie->GetGeometry()[1].Z()-ie->GetGeometry()[3].Z()));
	   
  // 	      dist23=sqrt((ie->GetGeometry()[2].X()-ie->GetGeometry()[3].X())*(ie->GetGeometry()[2].X()-ie->GetGeometry()[3].X())+
  // 			  (ie->GetGeometry()[2].Y()-ie->GetGeometry()[3].Y())*(ie->GetGeometry()[2].Y()-ie->GetGeometry()[3].Y())+
  // 			  (ie->GetGeometry()[2].Z()-ie->GetGeometry()[3].Z())*(ie->GetGeometry()[2].Z()-ie->GetGeometry()[3].Z()));

  // 	      double distWall=0;
  // 	      double safetyCoefficient3D=1.8;
  // 	      bool dangerousElement=false;
  // 	      if(rigidNodes>1){
  // 		if(ie->GetGeometry()[0].Is(RIGID) && ie->GetGeometry()[1].Is(RIGID)){
  // 		  distWall=sqrt((ie->GetGeometry()[0].X()-ie->GetGeometry()[1].X())*(ie->GetGeometry()[0].X()-ie->GetGeometry()[1].X())+
  // 				(ie->GetGeometry()[0].Y()-ie->GetGeometry()[1].Y())*(ie->GetGeometry()[0].Y()-ie->GetGeometry()[1].Y())+
  // 				(ie->GetGeometry()[0].Z()-ie->GetGeometry()[1].Z())*(ie->GetGeometry()[0].Z()-ie->GetGeometry()[1].Z()));

  // 		}else if(ie->GetGeometry()[0].Is(RIGID) && ie->GetGeometry()[2].Is(RIGID)){
  // 		  distWall=sqrt((ie->GetGeometry()[0].X()-ie->GetGeometry()[2].X())*(ie->GetGeometry()[0].X()-ie->GetGeometry()[2].X())+
  // 				(ie->GetGeometry()[0].Y()-ie->GetGeometry()[2].Y())*(ie->GetGeometry()[0].Y()-ie->GetGeometry()[2].Y())+
  // 				(ie->GetGeometry()[0].Z()-ie->GetGeometry()[2].Z())*(ie->GetGeometry()[0].Z()-ie->GetGeometry()[2].Z()));

  // 		}else if(ie->GetGeometry()[0].Is(RIGID) && ie->GetGeometry()[3].Is(RIGID)){
  // 		  distWall=sqrt((ie->GetGeometry()[0].X()-ie->GetGeometry()[3].X())*(ie->GetGeometry()[0].X()-ie->GetGeometry()[3].X())+
  // 				(ie->GetGeometry()[0].Y()-ie->GetGeometry()[3].Y())*(ie->GetGeometry()[0].Y()-ie->GetGeometry()[3].Y())+
  // 				(ie->GetGeometry()[0].Z()-ie->GetGeometry()[3].Z())*(ie->GetGeometry()[0].Z()-ie->GetGeometry()[3].Z()));

  // 		}else if(ie->GetGeometry()[1].Is(RIGID) && ie->GetGeometry()[2].Is(RIGID)){
  // 		  distWall=sqrt((ie->GetGeometry()[1].X()-ie->GetGeometry()[2].X())*(ie->GetGeometry()[1].X()-ie->GetGeometry()[2].X())+
  // 				(ie->GetGeometry()[1].Y()-ie->GetGeometry()[2].Y())*(ie->GetGeometry()[1].Y()-ie->GetGeometry()[2].Y())+
  // 				(ie->GetGeometry()[1].Z()-ie->GetGeometry()[2].Z())*(ie->GetGeometry()[1].Z()-ie->GetGeometry()[2].Z()));

  // 		}else if(ie->GetGeometry()[1].Is(RIGID) && ie->GetGeometry()[3].Is(RIGID)){
  // 		  distWall=sqrt((ie->GetGeometry()[1].X()-ie->GetGeometry()[3].X())*(ie->GetGeometry()[1].X()-ie->GetGeometry()[3].X())+
  // 				(ie->GetGeometry()[1].Y()-ie->GetGeometry()[3].Y())*(ie->GetGeometry()[1].Y()-ie->GetGeometry()[3].Y())+
  // 				(ie->GetGeometry()[1].Z()-ie->GetGeometry()[3].Z())*(ie->GetGeometry()[1].Z()-ie->GetGeometry()[3].Z()));

  // 		}else if(ie->GetGeometry()[2].Is(RIGID) && ie->GetGeometry()[3].Is(RIGID)){
  // 		  distWall=sqrt((ie->GetGeometry()[2].X()-ie->GetGeometry()[3].X())*(ie->GetGeometry()[2].X()-ie->GetGeometry()[3].X())+
  // 				(ie->GetGeometry()[2].Y()-ie->GetGeometry()[3].Y())*(ie->GetGeometry()[2].Y()-ie->GetGeometry()[3].Y())+
  // 				(ie->GetGeometry()[2].Z()-ie->GetGeometry()[3].Z())*(ie->GetGeometry()[2].Z()-ie->GetGeometry()[3].Z()));

  // 		}

  // 		if(dist01<distWall*safetyCoefficient3D && (ie->GetGeometry()[0].Is(RIGID) || ie->GetGeometry()[1].Is(RIGID))){
  // 		  dist01=0;
  // 		}
  // 		if(dist02<distWall*safetyCoefficient3D && (ie->GetGeometry()[0].Is(RIGID) || ie->GetGeometry()[2].Is(RIGID))){
  // 		  dist02=0;
  // 		}
  // 		if(dist03<distWall*safetyCoefficient3D && (ie->GetGeometry()[0].Is(RIGID) || ie->GetGeometry()[3].Is(RIGID))){
  // 		  dist03=0;
  // 		}
  // 		if(dist12<distWall*safetyCoefficient3D && (ie->GetGeometry()[1].Is(RIGID) || ie->GetGeometry()[2].Is(RIGID))){
  // 		  dist12=0;
  // 		}
  // 		if(dist13<distWall*safetyCoefficient3D && (ie->GetGeometry()[1].Is(RIGID) || ie->GetGeometry()[3].Is(RIGID))){
  // 		  dist13=0;
  // 		}
  // 		if(dist23<distWall*safetyCoefficient3D && (ie->GetGeometry()[2].Is(RIGID) || ie->GetGeometry()[3].Is(RIGID))){
  // 		  dist23=0;
  // 		}

  // 	      }

  // 	      if(rigidNodes==1){
  // 		if(ie->GetGeometry()[0].Is(RIGID)){
  // 		  dist01=0;
  // 		  dist02=0;
  // 		  dist03=0;
  // 		}
  // 		if(ie->GetGeometry()[1].Is(RIGID)){
  // 		  dist01=0;
  // 		  dist12=0;
  // 		  dist13=0;
  // 		}
  // 		if(ie->GetGeometry()[2].Is(RIGID)){
  // 		  dist02=0;
  // 		  dist12=0;
  // 		  dist23=0;
  // 		}
  // 		if(ie->GetGeometry()[3].Is(RIGID)){
  // 		  dist03=0;
  // 		  dist13=0;
  // 		  dist23=0;
  // 		}
  // 	      }

  // 	      if((dist01==0 && dist02==0 && dist03==0 && dist12==0 && dist13==0 && dist23==0) || rigidNodes>2){
  // 		dangerousElement=true;
  // 	      }

  // 	      //just to fill the vector
  // 	      if(dangerousElement==false){
  // 		count++;
  // 		if(count<ElementsToRefine){

  // 		  if(dist01>dist02 && dist01>dist03 && dist01>dist12 && dist01>dist13 && dist01>dist23){
  // 		    NewPosition=    (ie->GetGeometry()[0].Coordinates()                         +ie->GetGeometry()[1].Coordinates())*0.5;
  // 		    NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY))*0.5;
  // 		    NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION))*0.5;
  // 		    NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE))*0.5;
  // 		  }else if(dist02>dist01 && dist02>dist03 && dist02>dist12 && dist02>dist13 && dist02>dist23){
  // 		    NewPosition=    (ie->GetGeometry()[0].Coordinates()                         +ie->GetGeometry()[2].Coordinates())*0.5;
  // 		    NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY))*0.5;
  // 		    NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION))*0.5;
  // 		    NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE))*0.5;
  // 		  }else if(dist03>dist01 && dist03>dist02 && dist03>dist12 && dist03>dist13 && dist03>dist23){
  // 		    NewPosition=    (ie->GetGeometry()[0].Coordinates()                         +ie->GetGeometry()[3].Coordinates())*0.5;
  // 		    NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[3].FastGetSolutionStepValue(VELOCITY))*0.5;
  // 		    NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[3].FastGetSolutionStepValue(ACCELERATION))*0.5;
  // 		    NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[3].FastGetSolutionStepValue(PRESSURE))*0.5;
  // 		  }else if(dist12>dist01 && dist12>dist02 && dist12>dist03 && dist12>dist13 && dist12>dist23){
  // 		    NewPosition=    (ie->GetGeometry()[1].Coordinates()                         +ie->GetGeometry()[2].Coordinates())*0.5;
  // 		    NewVelocity=    (ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY))*0.5;
  // 		    NewAcceleration=(ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION))*0.5;
  // 		    NewPressure=    (ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE))*0.5;
  // 		  }else if(dist13>dist01 && dist13>dist02 && dist13>dist03 && dist13>dist12 && dist13>dist23){
  // 		    NewPosition=    (ie->GetGeometry()[1].Coordinates()                         +ie->GetGeometry()[3].Coordinates())*0.5;
  // 		    NewVelocity=    (ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[3].FastGetSolutionStepValue(VELOCITY))*0.5;
  // 		    NewAcceleration=(ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[3].FastGetSolutionStepValue(ACCELERATION))*0.5;
  // 		    NewPressure=    (ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[3].FastGetSolutionStepValue(PRESSURE))*0.5;
  // 		  }else if(dist23>dist01 && dist23>dist02 && dist23>dist03 && dist23>dist12 && dist23>dist13){
  // 		    NewPosition=    (ie->GetGeometry()[2].Coordinates()                         +ie->GetGeometry()[3].Coordinates())*0.5;
  // 		    NewVelocity=    (ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[3].FastGetSolutionStepValue(VELOCITY))*0.5;
  // 		    NewAcceleration=(ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[3].FastGetSolutionStepValue(ACCELERATION))*0.5;
  // 		    NewPressure=    (ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[3].FastGetSolutionStepValue(PRESSURE))*0.5;
  // 		  }else{
  // 		    std::cout<<"dist01 "<<dist01<<" dist02 "<<dist02<<" dist12 "<<dist12<<std::endl;
  // 		    NewPosition=    (ie->GetGeometry()[0].Coordinates()+
  // 				     ie->GetGeometry()[1].Coordinates()+
  // 				     ie->GetGeometry()[2].Coordinates()+
  // 				     ie->GetGeometry()[3].Coordinates())*0.25;
  // 		    NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)+
  // 				     ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY)+
  // 				     ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY)+
  // 				     ie->GetGeometry()[3].FastGetSolutionStepValue(VELOCITY))*0.25;
  // 		    NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+
  // 				     ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION)+
  // 				     ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION)+
  // 				     ie->GetGeometry()[3].FastGetSolutionStepValue(ACCELERATION))*0.25;
  // 		    NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)+
  // 				     ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)+
  // 				     ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)+
  // 				     ie->GetGeometry()[3].FastGetSolutionStepValue(PRESSURE))*0.25;
  // 		  }
  // 		  BiggestVolumes[CountNodes]=ElementalVolume;
  // 		  NewPositions[CountNodes]=NewPosition;
  // 		  NewVelocities[CountNodes]=NewVelocity;
  // 		  NewAccelerations[CountNodes]=NewAcceleration;
  // 		  NewPressures[CountNodes]=NewPressure;
  // 		  if(freesurfaceNodes<2){
  // 		    midPointInterpolation[CountNodes]=false;
  // 		  }else{
  // 		    midPointInterpolation[CountNodes]=true;
  // 		  }
  // 		  CountNodes++;
  // 		}else if (freesurfaceNodes<3){
  // 		  for(int nn= 0; nn< ElementsToRefine; nn++)
  // 		    {
  // 		      if(ElementalVolume>BiggestVolumes[nn]){
  // 			if(dist01>dist02 && dist01>dist03 && dist01>dist12 && dist01>dist13 && dist01>dist23){
  // 			  NewPosition=    (ie->GetGeometry()[0].Coordinates()                         +ie->GetGeometry()[1].Coordinates())*0.5;
  // 			  NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY))*0.5;
  // 			  NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION))*0.5;
  // 			  NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE))*0.5;
  // 			}else if(dist02>dist01 && dist02>dist03 && dist02>dist12 && dist02>dist13 && dist02>dist23){
  // 			  NewPosition=    (ie->GetGeometry()[0].Coordinates()                         +ie->GetGeometry()[2].Coordinates())*0.5;
  // 			  NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY))*0.5;
  // 			  NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION))*0.5;
  // 			  NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE))*0.5;
  // 			}else if(dist03>dist01 && dist03>dist02 && dist03>dist12 && dist03>dist13 && dist03>dist23){
  // 			  NewPosition=    (ie->GetGeometry()[0].Coordinates()                         +ie->GetGeometry()[3].Coordinates())*0.5;
  // 			  NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[3].FastGetSolutionStepValue(VELOCITY))*0.5;
  // 			  NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[3].FastGetSolutionStepValue(ACCELERATION))*0.5;
  // 			  NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[3].FastGetSolutionStepValue(PRESSURE))*0.5;
  // 			}else if(dist12>dist01 && dist12>dist02 && dist12>dist03 && dist12>dist13 && dist12>dist23){
  // 			  NewPosition=    (ie->GetGeometry()[1].Coordinates()                         +ie->GetGeometry()[2].Coordinates())*0.5;
  // 			  NewVelocity=    (ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY))*0.5;
  // 			  NewAcceleration=(ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION))*0.5;
  // 			  NewPressure=    (ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE))*0.5;
  // 			}else if(dist13>dist01 && dist13>dist02 && dist13>dist03 && dist13>dist12 && dist13>dist23){
  // 			  NewPosition=    (ie->GetGeometry()[1].Coordinates()                         +ie->GetGeometry()[3].Coordinates())*0.5;
  // 			  NewVelocity=    (ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[3].FastGetSolutionStepValue(VELOCITY))*0.5;
  // 			  NewAcceleration=(ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[3].FastGetSolutionStepValue(ACCELERATION))*0.5;
  // 			  NewPressure=    (ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[3].FastGetSolutionStepValue(PRESSURE))*0.5;
  // 			}else if(dist23>dist01 && dist23>dist02 && dist23>dist03 && dist23>dist12 && dist23>dist13){
  // 			  NewPosition=    (ie->GetGeometry()[2].Coordinates()                         +ie->GetGeometry()[3].Coordinates())*0.5;
  // 			  NewVelocity=    (ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[3].FastGetSolutionStepValue(VELOCITY))*0.5;
  // 			  NewAcceleration=(ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[3].FastGetSolutionStepValue(ACCELERATION))*0.5;
  // 			  NewPressure=    (ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[3].FastGetSolutionStepValue(PRESSURE))*0.5;
  // 			}else{
  // 			  std::cout<<"dist01 "<<dist01<<" dist02 "<<dist02<<" dist12 "<<dist12<<std::endl;
  // 			  NewPosition=    (ie->GetGeometry()[0].Coordinates()+
  // 					   ie->GetGeometry()[1].Coordinates()+
  // 					   ie->GetGeometry()[2].Coordinates()+
  // 					   ie->GetGeometry()[3].Coordinates())*0.25;
  // 			  NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)+
  // 					   ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY)+
  // 					   ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY)+
  // 					   ie->GetGeometry()[3].FastGetSolutionStepValue(VELOCITY))*0.25;
  // 			  NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+
  // 					   ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION)+
  // 					   ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION)+
  // 					   ie->GetGeometry()[3].FastGetSolutionStepValue(ACCELERATION))*0.25;
  // 			  NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)+
  // 					   ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)+
  // 					   ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)+
  // 					   ie->GetGeometry()[3].FastGetSolutionStepValue(PRESSURE))*0.25;
  // 			}
  // 			// std::cout<<"NewPosition[2]"<<NewPosition[2]<<" "<<ie->GetGeometry()[0].Z()<<", "<<ie->GetGeometry()[1].Z()<<", "<<ie->GetGeometry()[2].Z()<<", "<<ie->GetGeometry()[3].Z()<<" volume:"<<ElementalVolume<<std::endl;
  // 			BiggestVolumes[nn]=ElementalVolume;
  // 			NewPositions[nn]=NewPosition;
  // 			NewVelocities[nn]=NewVelocity;
  // 			NewAccelerations[nn]=NewAcceleration;
  // 			NewPressures[nn]=NewPressure;
  // 			if(freesurfaceNodes<2){
  // 			  midPointInterpolation[CountNodes]=false;
  // 			}else{
  // 			  midPointInterpolation[CountNodes]=true;
  // 			}
  // 			break;
  // 		      }
  // 		    }
  // 		}
  // 	      }
  // 	    }


  //   	  }	 
  //     }

  //   std::vector<Node<3>::Pointer > list_of_new_nodes;
  //   // int idNode=mrModelPart.NumberOfNodes(mMeshId);

  //   unsigned int initial_node_size = mrModelPart.Nodes().size()+1; //total model part node size

  //   for(uint nn= 0; nn< NewPositions.size(); nn++)
  //     {
  //   	unsigned int id = initial_node_size + nn ;

  //   	double  x = NewPositions[nn][0];
  //   	double  y = NewPositions[nn][1];
  //   	double  z = 0; 
  //   	if(dimension==3)
  //   	  z=NewPositions[nn][2];

  //   	Node<3>::Pointer pnode = mrModelPart.CreateNewNode(id,x,y,z);

  //   	pnode->Set(NEW_ENTITY); //not boundary

  //   	pnode->SetBufferSize(mrModelPart.NodesBegin(mMeshId)->GetBufferSize() );
  // 	list_of_new_nodes.push_back( pnode );

  // 	if(mrRemesh.InputInitializedFlag){
  // 	  mrRemesh.NodalPreIds.push_back( pnode->Id() );
  // 	  pnode->SetId(mrRemesh.NodalPreIds.size()-1);
  // 	}

  //   	if(mMeshId!=0)
  //   	  mrModelPart.AddNode(pnode,mMeshId);

  // 	//generating the dofs
  // 	Node<3>::DofsContainerType& reference_dofs = (mrModelPart.NodesBegin(mMeshId))->GetDofs();
  // 	for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
  // 	  {
  // 	    Node<3>::DofType& rDof = *iii;
  // 	    Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );

  // 	    (p_new_dof)->FreeDof();
  // 	  }
      
  //     }


  //   if( list_of_new_nodes.size() > 0)
  //     ProjectVariablesToNewNodes( list_of_new_nodes );

  //   //set the coordinates to the original value
  //   const array_1d<double,3> ZeroNormal(3,0.0);
  //   uint count=0;
  //   for(std::vector<Node<3>::Pointer>::iterator it =  list_of_new_nodes.begin(); it!=list_of_new_nodes.end(); it++)
  //     {
  // 	const array_1d<double,3>& displacement = (*it)->FastGetSolutionStepValue(DISPLACEMENT);
  // 	(*it)->X0() = (*it)->X() - displacement[0];
  // 	(*it)->Y0() = (*it)->Y() - displacement[1];
  // 	(*it)->Z0() = (*it)->Z() - displacement[2];

  // 	// midPointInterpolation[count]=true;
  // 	if(midPointInterpolation[count]==true){
  // 	  (*it)->FastGetSolutionStepValue(VELOCITY)=NewVelocities[count];
  // 	  (*it)->FastGetSolutionStepValue(ACCELERATION)=NewAccelerations[count];
  // 	  (*it)->FastGetSolutionStepValue(PRESSURE)=NewPressures[count];
  // 	  array_1d<double,3> NewInterpolatedVelocity=(*it)->FastGetSolutionStepValue(VELOCITY);
  // 	  array_1d<double,3> NewInterpolatedAcc=(*it)->FastGetSolutionStepValue(ACCELERATION);
  // 	  double NewInterpolatedPre=(*it)->FastGetSolutionStepValue(PRESSURE);
  // 	}
  // 	count++;

  // 	double density=(*it)->FastGetSolutionStepValue(DENSITY);
  // 	if(density==0){
  // 	  std::cout<<"..................................  the density of this node was 0 now is 1000 "<<std::endl;
  // 	  (*it)->FastGetSolutionStepValue(DENSITY)=1000.0;
  // 	}

  // 	//correct contact_normal interpolation
  // 	if( (*it)->SolutionStepsDataHas(CONTACT_FORCE) )
  // 	  noalias((*it)->GetSolutionStepValue(CONTACT_FORCE)) = ZeroNormal;
		    
  // 	(*it)->SetValue(DOMAIN_LABEL,mMeshId);
	  
  //     }








  //   mrRemesh.InputInitializedFlag=false;




  //   //mModelerUtilities.SetNodes(mrModelPart,mrRemesh,mMeshId);

  //   // if(!mrRemesh.InputInitializedFlag){
  //   //   mrRemesh.InputInitializedFlag=true;
  //   // }

  //   KRATOS_CATCH( "" )

  //     }














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
    GenerateNewNodesBeforeMeshingProcess& operator=(GenerateNewNodesBeforeMeshingProcess const& rOther);


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
                                  GenerateNewNodesBeforeMeshingProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GenerateNewNodesBeforeMeshingProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_GENERATE_NEW_NODES_BEFORE_MESHING_PROCESS_H_INCLUDED  defined 


