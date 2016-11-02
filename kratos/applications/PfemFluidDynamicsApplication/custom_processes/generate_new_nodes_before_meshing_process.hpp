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
      // std::cout<<"     I will generate "<<mrRemesh.Info->RemovedNodes <<" new nodes"<<std::endl;

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
    double viscosity=0;
    double density=0;
    double bulkModulus=0;
    array_1d<double, 3 > VolumeAcceleration(3,0.0);
    bool fluidNodeFound=false;
    int CountNodes=0;
    double athird=0.3333333333333333;
    NewPositions.resize(ElementsToRefine);
    NewVelocities.resize(ElementsToRefine);
    NewAccelerations.resize(ElementsToRefine);
    NewPressures.resize(ElementsToRefine);
    NewPressures.resize(ElementsToRefine);
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


	    ///////////////// FOR 2D CASE: choose the right (big and safe) elements to refine and compute the new node  position and variables //////////////////////
	    if(dimension==2)
	      {

		for(uint pn=0; pn<nds; pn++)
		  {
		    if(ie->GetGeometry()[pn].Is(RIGID)){
		      rigidNodes++;
		    }
		    if(fluidNodeFound==false && ie->GetGeometry()[pn].Is(FLUID) && !ie->GetGeometry()[pn].Is(RIGID)){
		      viscosity=ie->GetGeometry()[pn].FastGetSolutionStepValue(VISCOSITY);
		      density=ie->GetGeometry()[pn].FastGetSolutionStepValue(DENSITY);
		      bulkModulus=ie->GetGeometry()[pn].FastGetSolutionStepValue(BULK_MODULUS);
		      VolumeAcceleration = ie->GetGeometry()[pn].FastGetSolutionStepValue(VOLUME_ACCELERATION);
		      fluidNodeFound=true;
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
		array_1d<double,2> CoorDifference(2,0.0);
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
		  // count++;
		  if(count<ElementsToRefine){

		    array_1d<double,2> NewPosition(2,0.0);
		    array_1d<double,2> NewVelocity(2,0.0);
		    array_1d<double,2> NewAcceleration(2,0.0);
		    double NewPressure=0;
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
		    }
		    else{
		      // for problematic elements
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
		    // std::cout<<NewPositions[CountNodes][0]<<" "<<NewPositions[CountNodes][1]<<"_volume_ "<<ElementalVolume<<std::endl;

		    CountNodes++;
		  }else if (freesurfaceNodes<2 && rigidNodes<2){

		    for(int nn= 0; nn< ElementsToRefine; nn++)
		      {

			if(ElementalVolume>BiggestVolumes[nn]){

			  array_1d<double,2> NewPosition(2,0.0);
			  array_1d<double,2> NewVelocity(2,0.0);
			  array_1d<double,2> NewAcceleration(2,0.0);
			  double NewPressure=0;
			  uint maxCount=3;
			  double LargestEdge=0;

			  for(uint i=0; i<3; i++)
			    {
			      if(Edges[i]>LargestEdge){
				maxCount=i;
				LargestEdge=Edges[i];
			      }
			    }
			  bool suitableElement=true;
			  if(maxCount<3){
			    NewPosition=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].Coordinates()                         +ie->GetGeometry()[SecondEdgeNode[maxCount]].Coordinates())*0.5;			    
			    NewVelocity=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(VELOCITY))*0.5;
			    NewAcceleration=(ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(ACCELERATION))*0.5;
			    NewPressure=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(PRESSURE))*0.5;
			  }
			  else  if(rigidNodes==1 && freesurfaceNodes>0){
			    // suitableElement=false;
			    // for problematic elements
			    NewPosition=    (ie->GetGeometry()[0].Coordinates()                         +ie->GetGeometry()[1].Coordinates()                         +ie->GetGeometry()[2].Coordinates())*athird;
			    NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY))*athird;
			    NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION))*athird;
			    NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE))*athird;
			  }
			  for(int j= 0; j< ElementsToRefine; j++)
			    {
			      if(NewPositions[j][0]==NewPosition[0] && NewPositions[j][1]==NewPosition[1]){
				suitableElement=false;
				// std::cout<<"this node has been already used, skip this element"<<std::endl;
			      }
			    }
			  if(suitableElement==true){
			    BiggestVolumes[nn]=ElementalVolume;
			    NewPositions[nn]=NewPosition;
			    NewVelocities[nn]=NewVelocity;
			    NewAccelerations[nn]=NewAcceleration;
			    NewPressures[nn]=NewPressure;
			  }
			  break;
			}
		      }
		  }
		  count++;
		}
	       
	      }
	    ///////////////// FOR 3D CASE: choose the right (big and safe) elements to refine and compute the new node  position and variables //////////////////////

	    else if(dimension==3){

	      for(uint pn=0; pn<nds; pn++)
		{
		  if(ie->GetGeometry()[pn].Is(RIGID)){
		    rigidNodes++;
		  }else if(DofsFound==false){
		    reference_dofs = ie->GetGeometry()[pn].GetDofs();
		    DofsFound=true;
		  }
		  if(fluidNodeFound==false && ie->GetGeometry()[pn].Is(FLUID) && !ie->GetGeometry()[pn].Is(RIGID)){
		    viscosity=ie->GetGeometry()[pn].FastGetSolutionStepValue(VISCOSITY);
		    density=ie->GetGeometry()[pn].FastGetSolutionStepValue(DENSITY);
		    bulkModulus=ie->GetGeometry()[pn].FastGetSolutionStepValue(BULK_MODULUS);
		    VolumeAcceleration = ie->GetGeometry()[pn].FastGetSolutionStepValue(VOLUME_ACCELERATION);
		    fluidNodeFound=true;
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

	      }else if(rigidNodes==1){
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
		// count++;
		if(count<ElementsToRefine){ 

		  array_1d<double,3> NewPosition(3,0.0);
		  array_1d<double,3> NewVelocity(3,0.0);
		  array_1d<double,3> NewAcceleration(3,0.0);
		  double NewPressure=0;
		  uint maxCount=6;
		  double LargestEdge=0;

		  for(uint i=0; i<6; i++)
		    {
		      if(Edges[i]>LargestEdge){
			maxCount=i;
			LargestEdge=Edges[i];
		      }
		    }
		  // if(maxCount<6){
		  NewPosition=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].Coordinates()                         +ie->GetGeometry()[SecondEdgeNode[maxCount]].Coordinates())*0.5;
		  NewVelocity=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(VELOCITY))*0.5;
		  NewAcceleration=(ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(ACCELERATION))*0.5;
		  NewPressure=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(PRESSURE))*0.5;
		  // }
		  // else{
		  //   // //for problematic elements
		  //   // NewPosition=    (ie->GetGeometry()[0].Coordinates()+
		  //   // 		     ie->GetGeometry()[1].Coordinates()+
		  //   // 		     ie->GetGeometry()[2].Coordinates()+
		  //   // 		     ie->GetGeometry()[3].Coordinates())*0.25;
		  //   // NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)+
		  //   // 		     ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY)+
		  //   // 		     ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY)+
		  //   // 		     ie->GetGeometry()[3].FastGetSolutionStepValue(VELOCITY))*0.25;
		  //   // NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+
		  //   // 		     ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION)+
		  //   // 		     ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION)+
		  //   // 		     ie->GetGeometry()[3].FastGetSolutionStepValue(ACCELERATION))*0.25;
		  //   // NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)+
		  //   // 		     ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)+
		  //   // 		     ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)+
		  //   // 		     ie->GetGeometry()[3].FastGetSolutionStepValue(PRESSURE))*0.25;
		  // }
		  BiggestVolumes[CountNodes]=ElementalVolume;
		  NewPositions[CountNodes]=NewPosition;
		  NewVelocities[CountNodes]=NewVelocity;
		  NewAccelerations[CountNodes]=NewAcceleration;
		  NewPressures[CountNodes]=NewPressure;
		  CountNodes++;
		}else if (freesurfaceNodes<3){
		  for(int nn= 0; nn< ElementsToRefine; nn++)
		    {
		      if(ElementalVolume>BiggestVolumes[nn]){

			array_1d<double,3> NewPosition(3,0.0);
			array_1d<double,3> NewVelocity(3,0.0);
			array_1d<double,3> NewAcceleration(3,0.0);
			double NewPressure=0;
			uint maxCount=6;
			double LargestEdge=0;
			
			for(uint i=0; i<6; i++)
			  {
			    if(Edges[i]>LargestEdge){
			      maxCount=i;
			      LargestEdge=Edges[i];
			    }
			  }
			bool suitableElement=true;

			if(maxCount<6){

			  NewPosition=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].Coordinates()                         +ie->GetGeometry()[SecondEdgeNode[maxCount]].Coordinates())*0.5;
			  NewVelocity=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(VELOCITY)    +ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(VELOCITY))*0.5;
			  NewAcceleration=(ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(ACCELERATION)+ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(ACCELERATION))*0.5;
			  NewPressure=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].FastGetSolutionStepValue(PRESSURE)    +ie->GetGeometry()[SecondEdgeNode[maxCount]].FastGetSolutionStepValue(PRESSURE))*0.5;
			}
			else{
			  // std::cout<<"problematic element"<<std::endl;
			  suitableElement=false;		
			  // NewPosition=    (ie->GetGeometry()[0].Coordinates()+
			  // 		   ie->GetGeometry()[1].Coordinates()+
			  // 		   ie->GetGeometry()[2].Coordinates()+
			  // 		   ie->GetGeometry()[3].Coordinates())*0.25;
			  // NewVelocity=    (ie->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)+
			  // 		   ie->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY)+
			  // 		   ie->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY)+
			  // 		   ie->GetGeometry()[3].FastGetSolutionStepValue(VELOCITY))*0.25;
			  // NewAcceleration=(ie->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION)+
			  // 		   ie->GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION)+
			  // 		   ie->GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION)+
			  // 		   ie->GetGeometry()[3].FastGetSolutionStepValue(ACCELERATION))*0.25;
			  // NewPressure=    (ie->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)+
			  // 		   ie->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE)+
			  // 		   ie->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE)+
			  // 		   ie->GetGeometry()[3].FastGetSolutionStepValue(PRESSURE))*0.25;
			}

			for(int j= 0; j< ElementsToRefine; j++)
			  {
			    if(NewPositions[j][0]==NewPosition[0] && NewPositions[j][1]==NewPosition[1] && NewPositions[j][2]==NewPosition[2]){
			      suitableElement=false; //this is a repeated node, I have already choose this from another element
			    }
			  }
			if(suitableElement==true){
			  BiggestVolumes[nn]=ElementalVolume;
			  NewPositions[nn]=NewPosition;
			  NewVelocities[nn]=NewVelocity;
			  NewAccelerations[nn]=NewAcceleration;
			  NewPressures[nn]=NewPressure;
			}
			break;
		      }
		    }
		}
		count++;
	      }
	    }


    	  }	 
      }



    std::vector<Node<3>::Pointer > list_of_new_nodes;
    // int idNode=mrModelPart.NumberOfNodes(mMeshId);

    unsigned int initial_node_size = mrModelPart.Nodes().size()+1+ElementsToRefine; //total model part node size

    // int max_id =  mrModelPart.Nodes().back().Id();

    // if(initial_node_size>max_id){
    //   std::cout<<"ATTENTION!!! INITIAL_NODE_SIZE > MAX_ID !!!"<<std::endl;
    //   initial_node_size=max_id; 
    // }
    // unsigned int ID = ModelerUtilities::GetMaxNodeId(mrModelPart) + 1;

    //assign data to dofs
    VariablesList& VariablesList = mrModelPart.GetNodalSolutionStepVariablesList();

    for(uint nn= 0; nn< NewPositions.size(); nn++)
      {

    	unsigned int id = initial_node_size + nn ;

        // std::cout<<"NEW POINT ("<<id<<") COORDINATES x,y: "<<NewPositions[nn][0]<<" "<<NewPositions[nn][1];
        // std::cout<<"NEW INTERPOLATED VELOCITY: "<<NewVelocities[nn][0]<<" "<<NewVelocities[nn][1]<<" "<<NewVelocities[nn][2]<<std::endl;
        // std::cout<<"NEW INTERPOLATED ACCEL: "<<NewAccelerations[nn][0]<<" "<<NewAccelerations[nn][1]<<" "<<NewAccelerations[nn][2];
        // std::cout<<"NEW INTERPOLATED PRES: "<<NewPressures[nn];
        // std::cout<<" ELEMENT VOLUME: "<<BiggestVolumes[nn]<<std::endl;


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

	(*it)->FastGetSolutionStepValue(VELOCITY)=NewVelocities[count];
	(*it)->FastGetSolutionStepValue(ACCELERATION)=NewAccelerations[count];
	(*it)->FastGetSolutionStepValue(PRESSURE)=NewPressures[count];
	count++;

	(*it)->FastGetSolutionStepValue(DENSITY)=density;
	(*it)->FastGetSolutionStepValue(VISCOSITY)=viscosity;
	(*it)->FastGetSolutionStepValue(BULK_MODULUS)=bulkModulus;
	if((*it)->SolutionStepsDataHas(VOLUME_ACCELERATION) ){
	  (*it)->FastGetSolutionStepValue(VOLUME_ACCELERATION)=VolumeAcceleration;
	}else{
	  std::cout<<" this node has not volume acceleration!!!! "<<std::endl;
	}
	if(density==0){
	  std::cout<<" the density of this new node is 0 !!!! "<<std::endl;
	}
	if(viscosity==0){
	  std::cout<<" the viscosity of this new node is 0 !!!!  "<<std::endl;
	}
	if(bulkModulus==0){
	  std::cout<<" the bulkModulus of this new node is 0 !!!!  "<<std::endl;
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


