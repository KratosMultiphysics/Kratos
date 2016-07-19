//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_GENERATE_NEW_NODES_FOR_HOMOGENEOUS_MESH_PROCESS_H_INCLUDED )
#define  KRATOS_GENERATE_NEW_NODES_FOR_HOMOGENEOUS_MESH_PROCESS_H_INCLUDED


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

class GenerateNewNodesForHomogeneousMeshProcess
  : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION( GenerateNewNodesForHomogeneousMeshProcess );

    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GenerateNewNodesForHomogeneousMeshProcess(ModelPart& rModelPart,
					 ModelerUtilities::MeshingParameters& rRemeshingParameters,
					 ModelPart::IndexType MeshId,
					 int EchoLevel) 
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    {
      std::cout<<" generate_new_nodes_process_for_homogeneous_mesh CONSTRUCTOR "<<std::endl;

      mMeshId = MeshId;
      mEchoLevel = EchoLevel;
    }


    /// Destructor.
    virtual ~GenerateNewNodesForHomogeneousMeshProcess() {}


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
	std::cout<<" [ GENERATE NEW NODES: "<<std::endl;

      SelectNewNodesForHomogeneousMeshes();
      std::cout<<"                    I will generate "<<mrRemesh.Info->RemovedNodes <<" new nodes"<<std::endl;

      //Find out where the new nodes belong to:

      //creating an auxiliary list for the new nodes
      std::vector<Node<3>::Pointer > list_of_new_nodes;
      //std::vector<int> local_ids;
      
      //node to get the DOFs from
      Node<3>::DofsContainerType& reference_dofs = (mrModelPart.NodesBegin(mMeshId))->GetDofs();
      
      unsigned int initial_node_size = mrModelPart.Nodes().size()+1; //total model part node size
      
      ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin(mMeshId);	  
      const unsigned int dimension = element_begin->GetGeometry().WorkingSpaceDimension();

      double* OutPointList = mrRemesh.OutMesh.GetPointList();
      double* InPointList = mrRemesh.InMesh.GetPointList();

      int& InNumberOfPoints = mrRemesh.InMesh.GetNumberOfPoints();
      int& OutNumberOfPoints = mrRemesh.OutMesh.GetNumberOfPoints();
     
      //if points were added, new nodes must be added to ModelPart
      int j = 0;
      if (InNumberOfPoints > OutNumberOfPoints)
	{
	  for(int i = OutNumberOfPoints; i<InNumberOfPoints; i++)
	    {
	      unsigned int id = initial_node_size + j ;
	      int base = i*dimension;
	    
	      double& x = InPointList[base];
	      double& y = InPointList[base+1];
	      double  z = 0; 
	      if(dimension==3)
		z=OutPointList[base+2];
	    
	      //std::cout<<" domain node id "<<id<<" local id "<<i+1<<std::endl;
	      //std::cout<<" node creation position ("<<x<<", "<<y<<")"<<std::endl;
	      Node<3>::Pointer pnode = mrModelPart.CreateNewNode(id,x,y,z);
		
	      pnode->Set(NEW_ENTITY); //not boundary

	      //set to the main mesh (Mesh 0) to avoid problems in the NodalPreIds (number of nodes: change) in other methods
	      pnode->SetBufferSize(mrModelPart.NodesBegin(mMeshId)->GetBufferSize() );

	      list_of_new_nodes.push_back( pnode );

	      if(mrRemesh.InputInitializedFlag){
		mrRemesh.NodalPreIds.push_back( pnode->Id() );
		pnode->SetId(i+1);
	      }
	      //local_ids.push_back(i+1);

	      //set to the main mesh (Mesh 0) to avoid problems in the NodalPreIds (number of nodes: change) in other methods
	      if(mMeshId!=0)
		mrModelPart.AddNode(pnode,mMeshId);

	      //generating the dofs
	      for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
		{
		  Node<3>::DofType& rDof = *iii;
		  Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );

		  (p_new_dof)->FreeDof();
		}
		
	      j++;

	    }

	}
      
      //Inserted nodes
      mrRemesh.Info->InsertedNodes = OutNumberOfPoints-InNumberOfPoints;

      if( mEchoLevel > 0 )
	std::cout <<"   [ GENERATED NODES: ( added: " << mrRemesh.Info->InsertedNodes <<" ) ]"<<std::endl;
	

      //if points were added project variables
      if( list_of_new_nodes.size() > 0)
	ProjectVariablesToNewNodes( list_of_new_nodes );

      //set the coordinates to the original value
      const array_1d<double,3> ZeroNormal(3,0.0);
      
      for(std::vector<Node<3>::Pointer>::iterator it =  list_of_new_nodes.begin(); it!=list_of_new_nodes.end(); it++)
	{
	  const array_1d<double,3>& displacement = (*it)->FastGetSolutionStepValue(DISPLACEMENT);
	  (*it)->X0() = (*it)->X() - displacement[0];
	  (*it)->Y0() = (*it)->Y() - displacement[1];
	  (*it)->Z0() = (*it)->Z() - displacement[2];
	//correct contact_normal interpolation
	  if( (*it)->SolutionStepsDataHas(CONTACT_FORCE) )
	    noalias((*it)->GetSolutionStepValue(CONTACT_FORCE)) = ZeroNormal;
		    
	  (*it)->SetValue(DOMAIN_LABEL,mMeshId);
	  
	}

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
        return "GenerateNewNodesForHomogeneousMeshProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GenerateNewNodesForHomogeneousMeshProcess";
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

      ////////////////////////////////////////////////  new part  /////////////////////////////////////////////////
 
      const unsigned int dimension = mrModelPart.ElementsBegin(mMeshId)->GetGeometry().WorkingSpaceDimension();

    int ElementsToRefine=0;
    ElementsToRefine=mrRemesh.Info->RemovedNodes;
    std::vector<array_1d<double,3> > NewPoints;
    std::vector<double > BiggestVolumes;
    int CountNodes=0;
    NewPoints.resize(ElementsToRefine);
    BiggestVolumes.resize(ElementsToRefine);
  
    if(ElementsToRefine>0 )
      {
	int& OutNumberOfElements  = mrRemesh.OutMesh.GetNumberOfElements();
	int* OutElementList       = mrRemesh.OutMesh.GetElementList();

	ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin(mMeshId);	  
	ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin(mMeshId);
	const unsigned int nds = element_begin->GetGeometry().size();
	ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes(mMeshId);

    	for(int el = 0; el< OutNumberOfElements; el++)
    	  {

    	    if(mrRemesh.PreservedElements[el])
    	      {
		ModelerUtilities ModelerUtils;
		double ElementalVolume = 0;	
		uint rigidNodes=0;
		uint freesurfaceNodes=0;
		if(dimension==2)
		  {
		    Geometry<Node<3> > vertices;
		    array_1d<double,2> NewPoint;

		    for(uint pn=0; pn<nds; pn++)
		      {
			vertices.push_back(rNodes(OutElementList[el*nds+pn]));
			if(vertices[pn].Is(RIGID))
			  rigidNodes++;
			if(vertices[pn].Is(FREE_SURFACE))
			  freesurfaceNodes++;
		      }
	  		  

		    ElementalVolume = ModelerUtilities::CalculateTriangleArea(vertices[0].X(),vertices[0].Y(),
									      vertices[1].X(),vertices[1].Y(),
									      vertices[2].X(),vertices[2].Y());

		    //just to fill the vector
		    if(el<ElementsToRefine){
		      NewPoint[0]=(vertices[0].X()+vertices[1].X()+vertices[2].X())*0.33333333;
		      NewPoint[1]=(vertices[0].Y()+vertices[1].Y()+vertices[2].Y())*0.33333333;
		      BiggestVolumes[CountNodes]=ElementalVolume;
		      NewPoints[CountNodes]=NewPoint;
		      CountNodes++;
		    }else if (freesurfaceNodes<2 && rigidNodes<2){
		      for(int nn= 0; nn< ElementsToRefine; nn++)
			{
			  if(ElementalVolume>BiggestVolumes[nn]){
			    NewPoint[0]=(vertices[0].X()+vertices[1].X()+vertices[2].X())*0.33333333;
			    NewPoint[1]=(vertices[0].Y()+vertices[1].Y()+vertices[2].Y())*0.33333333;
			    BiggestVolumes[nn]=ElementalVolume;
			    NewPoints[nn]=NewPoint;
			    break;
			  }
			}
		    }
		  }else if(dimension==3){
		  std::cout<<"!!!!!!!!!!!!!!!!!!! FUNCTION TO TEST !!!!!!!!!!!!!!!!"<<std::endl;
		  std::cout<<"!!!!!!!!!!!!!!!!! ONLY 2D VERSION TESTED !!!!!!!!!!!!!!!"<<std::endl;
		  Geometry<Node<3> > vertices;
		  array_1d<double,2> NewPoint;

		  for(uint pn=0; pn<nds; pn++)
		    {
		      vertices.push_back(rNodes(OutElementList[el*nds+pn]));
		      if(vertices[pn].Is(RIGID))
			rigidNodes++;
		      if(vertices[pn].Is(FREE_SURFACE))
			freesurfaceNodes++;
		    }
	  		  

		  ElementalVolume = ModelerUtilities::CalculateTetrahedronVolume(vertices[0].X(),vertices[0].Y(),vertices[0].Z(),
										 vertices[1].X(),vertices[1].Y(),vertices[1].Z(),
										 vertices[2].X(),vertices[2].Y(),vertices[2].Z(),
										 vertices[3].X(),vertices[3].Y(),vertices[3].Z());

		  //just to fill the vector
		  if(el<ElementsToRefine){
		    NewPoint[0]=(vertices[0].X()+vertices[1].X()+vertices[2].X()+vertices[3].X())*0.25;
		    NewPoint[1]=(vertices[0].Y()+vertices[1].Y()+vertices[2].Y()+vertices[3].Y())*0.25;
		    NewPoint[2]=(vertices[0].Z()+vertices[1].Z()+vertices[2].Z()+vertices[3].Z())*0.25;
		    BiggestVolumes[CountNodes]=ElementalVolume;
		    NewPoints[CountNodes]=NewPoint;
		    CountNodes++;
		  }else if (freesurfaceNodes<2 && rigidNodes<2){
		    for(int nn= 0; nn< ElementsToRefine; nn++)
		      {
			if(ElementalVolume>BiggestVolumes[nn]){
			  NewPoint[0]=(vertices[0].X()+vertices[1].X()+vertices[2].X()+vertices[3].X())*0.25;
			  NewPoint[1]=(vertices[0].Y()+vertices[1].Y()+vertices[2].Y()+vertices[3].Y())*0.25;
			  NewPoint[2]=(vertices[0].Z()+vertices[1].Z()+vertices[2].Z()+vertices[3].Z())*0.25;
			  BiggestVolumes[nn]=ElementalVolume;
			  NewPoints[nn]=NewPoint;
			  break;
			}
		      }
		  }
		}
    	      }	  

    	  }	 
      }

    for(int nn= 0; nn< ElementsToRefine; nn++)
      {
        std::cout<<"NEW POINT COORDINATES: "<<NewPoints[nn][0]<<" "<<NewPoints[nn][1];
        std::cout<<"ELEMENT VOLUME: "<<BiggestVolumes[nn]<<std::endl;
      }


    //input mesh: NODES    
    mrRemesh.InMesh.CreatePointList((mrModelPart.Nodes(mMeshId).size()+ElementsToRefine), dimension);
    
    double* OutPointList    = mrRemesh.OutMesh.GetPointList();
    int& OutNumberOfPoints  = mrRemesh.OutMesh.GetNumberOfPoints();

    double* InPointList    = mrRemesh.InMesh.GetPointList();
    int& InNumberOfPoints  = mrRemesh.InMesh.GetNumberOfPoints();


    int base = 0;
    for(int i = 0; i<(OutNumberOfPoints); i++)
      {
	for(unsigned int j=0; j<dimension; j++)
	  InPointList[base+j] = OutPointList[base+j];
	base+=dimension;
      }

    // int* OutElementList    = mrRemesh.OutMesh.GetElementList();
    // int& OutNumberOfElements  = mrRemesh.OutMesh.GetNumberOfElements();

    // int* InElementList         = mrRemesh.InMesh.GetElementList();
    // int& InNumberOfElements   = mrRemesh.InMesh.GetNumberOfElements();

    // for(int i = 0; i<(OutNumberOfElements); i++)
    //   {
    // 	InElementList[i] = OutElementList[i];
    //   }

    // if(!mrRemesh.InputInitializedFlag){

    //   if((int)mrRemesh.NodalPreIds.size() != OutNumberOfPoints)
    // 	mrRemesh.NodalPreIds.resize(OutNumberOfPoints+1);

    //   std::fill( mrRemesh.NodalPreIds.begin(), mrRemesh.NodalPreIds.end(), 0 );
    // }
	 
    //writing the points coordinates in a vector and reordening the Id's
    ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin(mMeshId);

    std::cout<<"mrModelPart.Nodes(mMeshId).size() "<<mrModelPart.Nodes(mMeshId).size()<<std::endl;

    int direct = 1;

    int addedNodes=0;
    for(int i = OutNumberOfPoints; i<InNumberOfPoints; i++)
      {
	  
  	if(i>=int(mrModelPart.Nodes(mMeshId).size())){
	  for(unsigned int j=0; j<dimension; j++){
	    InPointList[base+j]   = NewPoints[addedNodes][j];
	  }
  	  std::cout<<"--------------------- NewPoints[addedNodes][0] "<<NewPoints[addedNodes][0]<<std::endl;
  	  addedNodes++;
  	} 
	    
  	base+=dimension;
  	direct+=1;
 
      }


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
    GenerateNewNodesForHomogeneousMeshProcess& operator=(GenerateNewNodesForHomogeneousMeshProcess const& rOther);


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
                                  GenerateNewNodesForHomogeneousMeshProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GenerateNewNodesForHomogeneousMeshProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_GENERATE_NEW_NODES_FOR_HOMOGENEOUS_MESH_PROCESS_H_INCLUDED  defined 


