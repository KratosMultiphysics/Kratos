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



 
 void SelectNewNodesForHomogeneousMeshes()

  {
    KRATOS_TRY

 
      const unsigned int dimension = mrModelPart.ElementsBegin(mMeshId)->GetGeometry().WorkingSpaceDimension();


    int ElementsToRefine=0;
    ElementsToRefine=mrRemesh.Info->RemovedNodes;
    std::cout<<" I will find "<<ElementsToRefine <<" new nodes"<<std::endl;
    std::vector<array_1d<double,3> > NewPositions;
    std::vector<double > BiggestVolumes;
    std::vector<array_1d< unsigned int,4 > > NodesIDToInterpolate;
    std::vector<bool > ElementalInterpolation;

    int CountNodes=0;
    // double athird=0.3333333333333333;
    NewPositions.resize(ElementsToRefine);
    BiggestVolumes.resize(ElementsToRefine);
    NodesIDToInterpolate.resize(ElementsToRefine);
    ElementalInterpolation.resize(ElementsToRefine);

    Node<3>::DofsContainerType& reference_dofs = (mrModelPart.NodesBegin())->GetDofs();

    bool DofsFound=false;

    if(ElementsToRefine>0 )
      {
 	ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin(mMeshId);	  
 	//ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin(mMeshId);
 	const unsigned int nds = element_begin->GetGeometry().size();
 	int count=0;
 	for(ModelPart::ElementsContainerType::const_iterator ie = element_begin; ie != mrModelPart.ElementsEnd(mMeshId); ie++)
 	  {
		
 	    const unsigned int dimension = ie->GetGeometry().WorkingSpaceDimension();

 	    ModelerUtilities ModelerUtils;
 	    double ElementalVolume = 0;	
 	    uint rigidNodes=0;
 	    uint freesurfaceNodes=0;
 	    bool toEraseNodeFound=false;

 	    ///////////////// FOR 2D CASE: choose the right (big and safe) elements to refine and compute the new node  position and variables //////////////////////
 	    if(dimension==2)
 	      {

 		for(uint pn=0; pn<nds; pn++)
 		  {
 		    if(ie->GetGeometry()[pn].Is(RIGID)){
 		      rigidNodes++;
 		    }
 		    if(ie->GetGeometry()[pn].Is(TO_ERASE)){
 		      toEraseNodeFound=true;
 		    }

 		    if(ie->GetGeometry()[pn].Is(FREE_SURFACE))
 		      freesurfaceNodes++;
 		  }
	  		  

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

 		if(dangerousElement==false && toEraseNodeFound==false){
 		  // count++;
 		  if(count<ElementsToRefine){

 		    array_1d<double,2> NewPosition(2,0.0);
 		    uint maxCount=3;
 		    double LargestEdge=0;

 		    for(uint i=0; i<3; i++)
 		      {
 			if(Edges[i]>LargestEdge){
 			  maxCount=i;
 			  LargestEdge=Edges[i];
 			}
 		      }
 		    NewPosition=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].Coordinates()+ie->GetGeometry()[SecondEdgeNode[maxCount]].Coordinates())*0.5;
 		    NodesIDToInterpolate[CountNodes][0]=ie->GetGeometry()[FirstEdgeNode[maxCount]].GetId();
 		    NodesIDToInterpolate[CountNodes][1]=ie->GetGeometry()[SecondEdgeNode[maxCount]].GetId();
 		    ElementalInterpolation[CountNodes]=false;
 		    // std::cout<<"...1st NodeToInterpol "<< NodesIDToInterpolate[CountNodes][0]<<std::endl;
 		    // std::cout<<"...2nd NodeToInterpol "<< NodesIDToInterpolate[CountNodes][1]<<std::endl;
 		    BiggestVolumes[CountNodes]=ElementalVolume;
 		    NewPositions[CountNodes]=NewPosition;
 		    CountNodes++;
 		  }else if (freesurfaceNodes<2 && rigidNodes<2){

 		    for(int nn= 0; nn< ElementsToRefine; nn++)
 		      {

 			if(ElementalVolume>BiggestVolumes[nn]){

 			  array_1d<double,2> NewPosition(2,0.0);
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
 			    NewPosition=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].Coordinates()+ie->GetGeometry()[SecondEdgeNode[maxCount]].Coordinates())*0.5;	
 			    for(int j= 0; j< ElementsToRefine; j++)
 			      {
 				if(NewPositions[j][0]==NewPosition[0] && NewPositions[j][1]==NewPosition[1]){
 				  suitableElement=false;
 				}
 			      }	
 			    if(suitableElement==true){	    
 			      NodesIDToInterpolate[nn][0]=ie->GetGeometry()[FirstEdgeNode[maxCount]].GetId();
 			      NodesIDToInterpolate[nn][1]=ie->GetGeometry()[SecondEdgeNode[maxCount]].GetId();
 			      ElementalInterpolation[nn]=false;
 			    }

 			  }
 			  else{
 			    suitableElement=false;
 			    std::cout<<"ATTENTION!!! suitableElement=false"<<std::endl;
 			  }
 			  if(suitableElement==true){
 			    BiggestVolumes[nn]=ElementalVolume;
 			    NewPositions[nn]=NewPosition;
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
 		  if(ie->GetGeometry()[pn].Is(TO_ERASE)){
 		    toEraseNodeFound=true;
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
 	      if(dangerousElement==false && toEraseNodeFound==false){
 		// count++;
 		if(count<ElementsToRefine){ 

 		  array_1d<double,3> NewPosition(3,0.0);
 		  uint maxCount=6;
 		  double LargestEdge=0;

 		  for(uint i=0; i<6; i++)
 		    {
 		      if(Edges[i]>LargestEdge){
 			maxCount=i;
 			LargestEdge=Edges[i];
 		      }
 		    }
 		  NewPosition=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].Coordinates()                         +ie->GetGeometry()[SecondEdgeNode[maxCount]].Coordinates())*0.5;
 		  NodesIDToInterpolate[CountNodes][0]=ie->GetGeometry()[FirstEdgeNode[maxCount]].GetId();
 		  NodesIDToInterpolate[CountNodes][1]=ie->GetGeometry()[SecondEdgeNode[maxCount]].GetId();
 		  ElementalInterpolation[CountNodes]=false;
 		  BiggestVolumes[CountNodes]=ElementalVolume;
 		  NewPositions[CountNodes]=NewPosition;
 		  CountNodes++;
 		}else if (freesurfaceNodes<3 && rigidNodes<3){
 		  for(int nn= 0; nn< ElementsToRefine; nn++)
 		    {
 		      if(ElementalVolume>BiggestVolumes[nn]){

 			array_1d<double,3> NewPosition(3,0.0);
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

 			  NewPosition=    (ie->GetGeometry()[FirstEdgeNode[maxCount]].Coordinates()+ie->GetGeometry()[SecondEdgeNode[maxCount]].Coordinates())*0.5;
 			  for(int j= 0; j< ElementsToRefine; j++)
 			    {
 			      if(NewPositions[j][0]==NewPosition[0] && NewPositions[j][1]==NewPosition[1] && NewPositions[j][2]==NewPosition[2]){
 				suitableElement=false; //this is a repeated node, I have already choose this from another element
 			      }
 			    }
 			  if(suitableElement==true){	    
 			    NodesIDToInterpolate[nn][0]=ie->GetGeometry()[FirstEdgeNode[maxCount]].GetId();
 			    NodesIDToInterpolate[nn][1]=ie->GetGeometry()[SecondEdgeNode[maxCount]].GetId();
 			    ElementalInterpolation[nn]=false;
 			  }
 			}
 			else{
 			  suitableElement=false;		
 			  std::cout<<"ATTENTION!!! suitableElement=false"<<std::endl;
 			}

 			if(suitableElement==true){
 			  BiggestVolumes[nn]=ElementalVolume;
 			  NewPositions[nn]=NewPosition;
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

    unsigned int initial_node_size = mrModelPart.Nodes().size()+1+ElementsToRefine; //total model part node size

    //assign data to dofs
    VariablesList& VariablesList = mrModelPart.GetNodalSolutionStepVariablesList();

    for(uint nn= 0; nn< NewPositions.size(); nn++)
      {

    	unsigned int id = initial_node_size + nn ;

        // std::cout<<"NEW POINT ("<<id<<") COORDINATES x,y: "<<NewPositions[nn][0]<<" "<<NewPositions[nn][1]<<std::endl;

 	// std::cout<<"interpolNodes: "<<NodesIDToInterpolate[nn][0]<<" "<<NodesIDToInterpolate[nn][1]<<std::endl;
 	// std::cout<<"true or false? "<<ElementalInterpolation[nn]<<std::endl;

    	double  x = NewPositions[nn][0];
    	double  y = NewPositions[nn][1];
    	double  z = 0; 
    	if(dimension==3)
    	  z=NewPositions[nn][2];

	
    	Node<3>::Pointer pnode = mrModelPart.CreateNewNode(id,x,y,z);
    	pnode->Set(NEW_ENTITY); //not boundary
 	list_of_new_nodes.push_back( pnode );
 	if(mrRemesh.InputInitializedFlag){
 	  mrRemesh.NodalPreIds.push_back( pnode->Id() );
 	  pnode->SetId(nn+1);
 	}

 	// //giving model part variables list to the node
 	pnode->SetSolutionStepVariablesList(&VariablesList);
	      
 	// //set buffer size
 	pnode->SetBufferSize(mrModelPart.GetBufferSize());

 	if(mMeshId!=0)
 	  mrModelPart.AddNode(pnode,mMeshId);

 	//generating the dofs
 	if(DofsFound==false){
 	  reference_dofs = (mrModelPart.NodesBegin())->GetDofs();
 	}
 	for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
 	  {
 	    Node<3>::DofType& rDof = *iii;
 	    Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
 	    (p_new_dof)->FreeDof();
 	  }

 	Node<3>::Pointer SlaveNode1 = mrModelPart.pGetNode(NodesIDToInterpolate[nn][0]);
 	Node<3>::Pointer SlaveNode2 = mrModelPart.pGetNode(NodesIDToInterpolate[nn][1]);	
 	InterpolateFromTwoNodes(pnode,SlaveNode1,SlaveNode2,VariablesList);
      
      }


    //set the coordinates to the original value
    const array_1d<double,3> ZeroNormal(3,0.0);
    for(std::vector<Node<3>::Pointer>::iterator it =  list_of_new_nodes.begin(); it!=list_of_new_nodes.end(); it++)
      {
 	const array_1d<double,3>& displacement = (*it)->FastGetSolutionStepValue(DISPLACEMENT);
 	(*it)->X0() = (*it)->X() - displacement[0];
 	(*it)->Y0() = (*it)->Y() - displacement[1];
 	(*it)->Z0() = (*it)->Z() - displacement[2];

 	// count++;
 	(*it)->Set(FLUID);
 	(*it)->Set(ACTIVE);
 	// std::cout<<"velocity_x "<<(*it)->FastGetSolutionStepValue(VELOCITY_X,0)<<std::endl;
 	// std::cout<<"velocity_x "<<(*it)->FastGetSolutionStepValue(VELOCITY_X,1)<<std::endl;
 	// std::cout<<"velocity_x "<<(*it)->FastGetSolutionStepValue(VELOCITY_X,2)<<std::endl;
 	// std::cout<<"pressure "<<(*it)->FastGetSolutionStepValue(PRESSURE,0)<<std::endl;
 	// std::cout<<"pressure "<<(*it)->FastGetSolutionStepValue(PRESSURE,1)<<std::endl;
 	// std::cout<<"pressure "<<(*it)->FastGetSolutionStepValue(PRESSURE,2)<<std::endl;
 	// std::cout<<"acc "<<(*it)->FastGetSolutionStepValue(ACCELERATION_X,0)<<std::endl;
 	// std::cout<<"acc "<<(*it)->FastGetSolutionStepValue(ACCELERATION_X,1)<<std::endl;
 	// std::cout<<"acc "<<(*it)->FastGetSolutionStepValue(ACCELERATION_X,2)<<std::endl;
 	// std::cout<<"bulkModulus "<<(*it)->FastGetSolutionStepValue(BULK_MODULUS)<<std::endl;
 	// std::cout<<"density "<<(*it)->FastGetSolutionStepValue(DENSITY)<<std::endl;
 	// std::cout<<"viscosity "<<(*it)->FastGetSolutionStepValue(VISCOSITY)<<std::endl;
 	//correct contact_normal interpolation
 	if( (*it)->SolutionStepsDataHas(CONTACT_FORCE) )
 	  noalias((*it)->GetSolutionStepValue(CONTACT_FORCE)) = ZeroNormal;
		    
 	(*it)->SetValue(DOMAIN_LABEL,mMeshId);
	  
      }

    mrRemesh.InputInitializedFlag=false;


    KRATOS_CATCH( "" )

      }


  void InterpolateFromTwoNodes( Node<3>::Pointer MasterNode,Node<3>::Pointer SlaveNode1,Node<3>::Pointer SlaveNode2,VariablesList& rVariablesList)
  { 
	  
    KRATOS_TRY

      unsigned int buffer_size = MasterNode->GetBufferSize();

     
    for(VariablesList::const_iterator i_variable =  rVariablesList.begin();  i_variable != rVariablesList.end() ; i_variable++)
      {
	//std::cout<<" name "<<i_variable->Name()<<std::endl;
	//std::cout<<" type "<<typeid(*i_variable).name()<<std::endl;
	std::string variable_name = i_variable->Name();
	if(KratosComponents<Variable<double> >::Has(variable_name))
	  {
	    //std::cout<<"double"<<std::endl;
	    Variable<double> variable = KratosComponents<Variable<double> >::Get(variable_name);
	    for(unsigned int step = 0; step<buffer_size; step++)
	      {
		//getting the data of the solution step
		double& node_data = MasterNode->FastGetSolutionStepValue(variable, step);
		  
		double& node0_data = SlaveNode1->FastGetSolutionStepValue(variable, step);
		double& node1_data = SlaveNode2->FastGetSolutionStepValue(variable, step);
		  
		node_data = (0.5*node0_data + 0.5*node1_data);
		  
	      }
	  }
	else if(KratosComponents<Variable<array_1d<double, 3> > >::Has(variable_name))
	  {
	    //std::cout<<"array1d"<<std::endl;
	    Variable<array_1d<double, 3> > variable = KratosComponents<Variable<array_1d<double, 3> > >::Get(variable_name);
	    for(unsigned int step = 0; step<buffer_size; step++)
	      {
		//getting the data of the solution step
		array_1d<double, 3>& node_data = MasterNode->FastGetSolutionStepValue(variable, step);
		  
		array_1d<double, 3>& node0_data = SlaveNode1->FastGetSolutionStepValue(variable, step);
		array_1d<double, 3>& node1_data = SlaveNode2->FastGetSolutionStepValue(variable, step);
		  
		node_data = (0.5*node0_data + 0.5*node1_data);		  
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
	    for(unsigned int step = 0; step<buffer_size; step++)
	      {
		//getting the data of the solution step
		Matrix& node_data = MasterNode->FastGetSolutionStepValue(variable, step);
		  
		Matrix& node0_data = SlaveNode1->FastGetSolutionStepValue(variable, step);
		Matrix& node1_data = SlaveNode2->FastGetSolutionStepValue(variable, step);
		  
		if( node_data.size1() > 0 && node_data.size2() ){
		  if( node_data.size1() == node0_data.size1() && node_data.size2() == node0_data.size2() &&
		      node_data.size1() == node1_data.size1() && node_data.size2() == node1_data.size2() ) {
		      
		    node_data = (0.5*node0_data + 0.5*node1_data);	       
		  }
		}
	      }

	  }
	else if(KratosComponents<Variable<Vector > >::Has(variable_name))
	  {
	    //std::cout<<"Vector"<<std::endl;
	    Variable<Vector >variable = KratosComponents<Variable<Vector > >::Get(variable_name);
	    for(unsigned int step = 0; step<buffer_size; step++)
	      {
		//getting the data of the solution step
		Vector& node_data = MasterNode->FastGetSolutionStepValue(variable, step);
		  
		Vector& node0_data = SlaveNode1->FastGetSolutionStepValue(variable, step);
		Vector& node1_data = SlaveNode2->FastGetSolutionStepValue(variable, step);
		  
		if( node_data.size() > 0 ){
		  if( node_data.size() == node0_data.size() &&
		      node_data.size() == node1_data.size()) {
		      
		    node_data = (0.5*node0_data + 0.5*node1_data);	       
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


