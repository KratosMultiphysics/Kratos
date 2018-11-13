//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_GENERATE_NEW_NODES_PROCESS_H_INCLUDED )
#define  KRATOS_GENERATE_NEW_NODES_PROCESS_H_INCLUDED


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

class GenerateNewNodesProcess
  : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION( GenerateNewNodesProcess );

    typedef ModelPart::NodeType                   NodeType;
    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GenerateNewNodesProcess(ModelPart& rModelPart,
			    ModelerUtilities::MeshingParameters& rRemeshingParameters,
			    int EchoLevel) 
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    {
      mEchoLevel = EchoLevel;
    }


    /// Destructor.
    virtual ~GenerateNewNodesProcess() {}


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

		if( mrModelPart.Name() != mrRemesh.SubModelPartName )
			std::cout<<" ModelPart Supplied do not corresponds to the Meshing Domain: ("<<mrModelPart.Name()<<" != "<<mrRemesh.SubModelPartName<<")"<<std::endl;

		//Find out where the new nodes belong to:

		//creating an auxiliary list for the new nodes
		std::vector<NodeType::Pointer > list_of_new_nodes;
      
		this->GenerateNewNodes(mrModelPart, list_of_new_nodes);
      
		if( mEchoLevel > 0 )
			std::cout <<"   [ GENERATED NODES: ( added: " << mrRemesh.Info->InsertedNodes <<" ) ]"<<std::endl;

		//project variables to new nodes from mesh elements
		if( list_of_new_nodes.size() > 0)
			this->ProjectVariablesToNewNodes( mrModelPart, list_of_new_nodes );

		//set node variables
		for(std::vector<NodeType::Pointer>::iterator it =  list_of_new_nodes.begin(); it!=list_of_new_nodes.end(); it++)
		{
			this->SetNewNodeVariables(mrModelPart,(*it));
		}

		this->SetNodesToModelPart(mrModelPart, list_of_new_nodes);
      
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
        return "GenerateNewNodesProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GenerateNewNodesProcess";
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

    int mEchoLevel;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{
     

    //*******************************************************************************************
    //*******************************************************************************************

    void GenerateNewNodes(ModelPart& rModelPart, std::vector<NodeType::Pointer>& list_of_nodes)
    {
		KRATOS_TRY

		NodeType::Pointer pNode;

		//center
		double xc = 0;
		double yc = 0;
		double zc = 0;
      
		//assign data to dofs
		NodeType::DofsContainerType& ReferenceDofs = rModelPart.Nodes().front().GetDofs();

		VariablesList& VariablesList = rModelPart.GetNodalSolutionStepVariablesList();

		unsigned int id = ModelerUtilities::GetMaxNodeId(rModelPart) + 1;

		// get dimension of the problem
		ModelPart::ElementsContainerType::iterator element_begin = rModelPart.ElementsBegin();	  
		const unsigned int dimension = element_begin->GetGeometry().WorkingSpaceDimension();

		// get mesh data from mesher
		double* OutPointList = mrRemesh.OutMesh.GetPointList();

		int& InNumberOfPoints  = mrRemesh.InMesh.GetNumberOfPoints();
		int& OutNumberOfPoints = mrRemesh.OutMesh.GetNumberOfPoints();

		// check if mesher created new nodes (inside the domain, no boundary nodes)
		if (OutNumberOfPoints > InNumberOfPoints)
		{
			// loop over new nodes that need to be added to the problem
			for(int i = InNumberOfPoints; i<OutNumberOfPoints; i++)
			{
				int base = i*dimension;
				// assign position from OutPointList
				xc = OutPointList[base];
				yc = OutPointList[base+1];
				zc = 0; 
				if(dimension==3)
					zc=OutPointList[base+2];

				//create a new node
				pNode = boost::make_shared< NodeType >( id, xc, yc, zc );

				//set new id
				if(mrRemesh.InputInitializedFlag)
				{
					mrRemesh.NodalPreIds.push_back( id );
					pNode->SetId(i+1);
					if( id > mrRemesh.NodeMaxId )
						mrRemesh.NodeMaxId = id;		
				}
		  
				//giving model part variables list to the node
				pNode->SetSolutionStepVariablesList(&VariablesList);
	      
				//set buffer size
				pNode->SetBufferSize(rModelPart.GetBufferSize());

				//generating the dofs
				for(Node<3>::DofsContainerType::iterator i_dof = ReferenceDofs.begin(); i_dof != ReferenceDofs.end(); i_dof++)
				{
					NodeType::DofType& rDof = *i_dof;
					NodeType::DofType::Pointer pNewDof = pNode->pAddDof( rDof );

					(pNewDof)->FreeDof();
				}
				
				list_of_nodes.push_back(pNode);
	      
				id++;
			}	
		}

		//Inserted nodes
		mrRemesh.Info->InsertedNodes = OutNumberOfPoints-InNumberOfPoints;
      
		KRATOS_CATCH( "" )
	}
  
    //**************************************************************************
    //**************************************************************************

    void ProjectVariablesToNewNodes(ModelPart& rModelPart, std::vector<NodeType::Pointer >& list_of_new_nodes)			    
    {

		KRATOS_TRY

		//defintions for spatial search
		typedef Node<3>                                  PointType;
		typedef Node<3>::Pointer                  PointPointerType;
		typedef std::vector<PointPointerType>          PointVector;
		typedef PointVector::iterator                PointIterator;
		//typedef std::vector<double>                 DistanceVector;
		typedef std::vector<double>::iterator     DistanceIterator;

		typedef Bucket<3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > BucketType;
		typedef Tree< KDTreePartition<BucketType> >     KdtreeType; //Kdtree
		//defintions for spatial search

		unsigned int  bucket_size = 20;
		KdtreeType    NodesTree(list_of_new_nodes.begin(),list_of_new_nodes.end(),bucket_size);

		//Find out where the new nodes belong to:
		std::vector<double> ShapeFunctionsN;    
		std::vector<VariablesListDataValueContainer> VariablesListVector(list_of_new_nodes.size());

		VariablesList&  variables_list = rModelPart.GetNodalSolutionStepVariablesList();

		//find the center and "radius" of the element
		double  radius = 0;
		Node<3> center(0,0.0,0.0,0.0);

		unsigned int MaximumNumberOfPointsInRadius = list_of_new_nodes.size();
		std::vector<Node<3>::Pointer> PointsInRadius (MaximumNumberOfPointsInRadius);
		std::vector<double>  PointsInRadiusDistances (MaximumNumberOfPointsInRadius);

		//geometry
		ModelPart::ElementsContainerType::iterator element_begin = rModelPart.ElementsBegin();	  
		const unsigned int nds = element_begin->GetGeometry().size();

		std::vector<std::vector<double> > ElementPointCoordinates(nds);
		std::vector<double> PointCoordinates(3);
		std::fill( PointCoordinates.begin(), PointCoordinates.end(), 0.0 );
		std::fill( ElementPointCoordinates.begin(), ElementPointCoordinates.end(), PointCoordinates );

		//loop over elements
		for(ModelPart::ElementsContainerType::const_iterator ie = rModelPart.ElementsBegin(); ie != rModelPart.ElementsEnd(); ie++)
		{
			//write element points coordinates into ElementPointCoordinates
			for(unsigned int i=0; i<ie->GetGeometry().size(); i++)
			{	
				PointCoordinates[0] = ie->GetGeometry()[i].X();
				PointCoordinates[1] = ie->GetGeometry()[i].Y();
				PointCoordinates[2] = ie->GetGeometry()[i].Z();

				ElementPointCoordinates[i] = PointCoordinates;
			}

			std::fill( PointCoordinates.begin(), PointCoordinates.end(), 0.0 );
			MeshDataTransferUtilities DataTransferUtilities;
			//writes the center point of ElementPointCoordintes into PointCoordintes; updated the largest radius
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
				//check if new node (PointCoordinates) is inside ElementPointCoordinates and calculates shape functions vector at PointCoordinates
				is_inside = ModelerUtilities::CalculatePosition( ElementPointCoordinates, PointCoordinates, ShapeFunctionsN );

				if(is_inside == true)
				{
					//do the interpolation based on shape functions vector for new node
					double alpha = 1; //1 to interpolate, 0 to leave the original data
					DataTransferUtilities.Interpolate( ie->GetGeometry(), ShapeFunctionsN, variables_list, *(it_found), alpha );		    
				}
			}
		}

		KRATOS_CATCH( "" )
    }


      //*******************************************************************************************
    //*******************************************************************************************

    virtual void SetNewNodeVariables(ModelPart& rModelPart, NodeType::Pointer& pNode)
    {
		KRATOS_TRY

		//set model part
		pNode->SetValue(MODEL_PART_NAME,rModelPart.Name());
     
		//set original position
		const array_1d<double,3>& Displacement = pNode->FastGetSolutionStepValue(DISPLACEMENT);
		pNode->X0() = pNode->X() - Displacement[0];
		pNode->Y0() = pNode->Y() - Displacement[1];
		pNode->Z0() = pNode->Z() - Displacement[2];
      
		//set contact force
		pNode->FastGetSolutionStepValue(CONTACT_FORCE).clear();
      
		KRATOS_CATCH( "" )
    }


    //*******************************************************************************************
    //*******************************************************************************************
 
    void SetNodesToModelPart(ModelPart& rModelPart, std::vector<NodeType::Pointer>& list_of_nodes)
    {
		KRATOS_TRY

		if(list_of_nodes.size())
		{
			//add new conditions: ( SOLID body model part )
			for(std::vector<NodeType::Pointer>::iterator i_node = list_of_nodes.begin(); i_node!= list_of_nodes.end(); i_node++)
			{
				rModelPart.Nodes().push_back(*(i_node));
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
    GenerateNewNodesProcess& operator=(GenerateNewNodesProcess const& rOther);


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
                                  GenerateNewNodesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GenerateNewNodesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_GENERATE_NEW_NODES_PROCESS_H_INCLUDED  defined 


