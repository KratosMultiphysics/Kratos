//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_REFINE_MESH_ELEMENTS_IN_EDGES_PROCESS_H_INCLUDED )
#define  KRATOS_REFINE_MESH_ELEMENTS_IN_EDGES_PROCESS_H_INCLUDED


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
//StepData: 
//Flags:    (checked) 
//          (set)     TO_SPLIT / VISITED
//          (modified)  
//          (reset)   TO_SPLIT / VISITED
//(set):=(set in this process)

namespace Kratos
{

///@name Kratos Classes
///@{

/// Refine Mesh Elements Process 2D and 3D
/** The process labels the nodes to be refined (TO_REFINE)
    if the ThresholdVariable  is larger than a ReferenceThreshold
*/

class RefineMeshElementsInEdgesProcess
  : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION( RefineMeshElementsInEdgesProcess );

    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RefineMeshElementsInEdgesProcess(ModelPart& rModelPart,
				     ModelerUtilities::MeshingParameters& rRemeshingParameters,
				     int EchoLevel) 
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    {
      mEchoLevel = EchoLevel;
    }


    /// Destructor.
    virtual ~RefineMeshElementsInEdgesProcess() {}


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

		if( ( mrRemesh.Refine->RefiningOptions.Is(ModelerUtilities::REFINE_ADD_NODES) || mrRemesh.Refine->RefiningOptions.Is(ModelerUtilities::REFINE_INSERT_NODES) ) )
		{
			//0.- Clean locally used flags
			this->CleanUsedFlags(mrModelPart);
			
			//1.- Select Elements TO_SPLIT (edge elements with all nodes as boundary-free surface)
			ModelPart::ElementsContainerType   BoundaryEdgedElements;
			ModelPart::ConditionsContainerType BoundaryEdgedConditions;
			
			this->SelectFullBoundaryEdgedElements(mrModelPart, BoundaryEdgedElements,BoundaryEdgedConditions); 
			
			//2.- Select Inside Faces to refine
			std::vector<Geometry< Node<3> > > ListOfFacesToSplit;
			this->SelectFacesToSplit(BoundaryEdgedElements,BoundaryEdgedConditions,ListOfFacesToSplit);
		
			//3.- Create and insert new nodes
			std::vector<Node<3>::Pointer>  ListOfNewNodes;
			this->GenerateNewNodes(mrModelPart,ListOfNewNodes,ListOfFacesToSplit);

			//4.- Insert new nodes to model part
			this->SetNodesToModelPart(mrModelPart, ListOfNewNodes);
			
			//0.- Clean locally used flags
			this->CleanUsedFlags(mrModelPart);
			
		}
      
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
        return "RefineMeshElementsInEdgesProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "RefineMeshElementsInEdgesProcess";
    }

    /// Print object's data.
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
    ///@name Un accessible methods
    ///@{

    //**************************************************************************
    //**************************************************************************

    void CleanUsedFlags(ModelPart& rModelPart)
    {
		KRATOS_TRY

		for(ModelPart::ElementsContainerType::const_iterator i_elem = rModelPart.ElementsBegin(); i_elem != rModelPart.ElementsEnd(); i_elem++)
		{
			i_elem->Set(TO_SPLIT,false);
		}

		for(ModelPart::NodesContainerType::const_iterator i_node = rModelPart.NodesBegin(); i_node != rModelPart.NodesEnd(); i_node++)
		{
			i_node->Set(TO_SPLIT,false);
		}
      
		KRATOS_CATCH( "" )
    }
  
    //**************************************************************************
    //**************************************************************************

    void SelectFullBoundaryEdgedElements(ModelPart& rModelPart,
					 ModelPart::ElementsContainerType& rBoundaryEdgedElements,
					 ModelPart::ConditionsContainerType& rBoundaryEdgedConditions)
    {
		KRATOS_TRY
	
		bool is_full_boundary = false;
		// loop over elements of model part
		// mark full boundary elements and its node TO_SPLIT
		for(ModelPart::ElementsContainerType::iterator i_elem = rModelPart.ElementsBegin(); i_elem != rModelPart.ElementsEnd(); i_elem++)
		{
			Geometry< Node<3> >& rGeometry = i_elem->GetGeometry();

			is_full_boundary = true;	    
			for(unsigned int i=0; i<rGeometry.size(); i++)
			{
				if( rGeometry[i].IsNot(BOUNDARY) ){
					is_full_boundary = false;
					break;
				}
			}

			if( is_full_boundary )
			{
				rBoundaryEdgedElements.push_back(*(i_elem.base()));
				i_elem->Set(TO_SPLIT,true);

				for(unsigned int i=0; i<rGeometry.size(); i++)
				{
					rGeometry[i].Set(TO_SPLIT,true);
				}
			}
		}

		// loop over conditions of the model part
		// add full boundary conditions (all nodes are already set TO_SPLIT) to rBoundaryEdgedConditions
		for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.ConditionsBegin(); i_cond != rModelPart.ConditionsEnd(); i_cond++)
		{
			if( i_cond->Is(BOUNDARY) ){
		
				Geometry< Node<3> >& rGeometry = i_cond->GetGeometry();
				is_full_boundary = true;	    
				for(unsigned int i=0; i<rGeometry.size(); i++)
				{
					if( rGeometry[i].IsNot(TO_SPLIT) ){
						is_full_boundary = false;
						break;
					}
				}
				
				if( is_full_boundary )
					rBoundaryEdgedConditions.push_back(*(i_cond.base()));
					
			}
		}
  
		KRATOS_CATCH( "" )
    }


    //**************************************************************************
    //**************************************************************************

    void SelectFacesToSplit(ModelPart::ElementsContainerType& rBoundaryEdgedElements,
			    ModelPart::ConditionsContainerType& rBoundaryEdgedConditions,
			    std::vector<Geometry< Node<3> > >& rListOfFacesToSplit)
    {
		KRATOS_TRY

		boost::numeric::ublas::matrix<unsigned int> lpofa; //connectivities of points defining faces
        boost::numeric::ublas::vector<unsigned int> lnofa; //number of points defining faces

		// set all condition flags to VISITED == false
		for(ModelPart::ConditionsContainerType::iterator i_cond = rBoundaryEdgedConditions.begin(); i_cond != rBoundaryEdgedConditions.end(); i_cond++)
		{    
			i_cond->Set(VISITED,false);
		}
	
		// loop over full boundary elements
		for(ModelPart::ElementsContainerType::iterator i_elem = rBoundaryEdgedElements.begin(); i_elem != rBoundaryEdgedElements.end(); i_elem++)
		{
			Geometry< Node<3> >& rGeometry = i_elem->GetGeometry();
	    
			rGeometry.NodesInFaces(lpofa);
			rGeometry.NumberNodesInFaces(lnofa);

			bool condition_found = false;
			
			// loop over faces of boundary element
			for(unsigned int i_face=0; i_face < rGeometry.FacesNumber(); i_face++)
			{
				condition_found = false;
				// loop over full boundary conditions
				for(ModelPart::ConditionsContainerType::iterator i_cond = rBoundaryEdgedConditions.begin(); i_cond != rBoundaryEdgedConditions.end(); i_cond++)
				{
					if( i_cond->IsNot(VISITED) )
					{
						i_cond->Set(VISITED,true);
						Geometry< Node<3> >& rConditionGeometry = i_cond->GetGeometry();

						ModelerUtilities ModelerUtils;
						condition_found = ModelerUtils.FindCondition(rConditionGeometry,rGeometry,lpofa,lnofa,i_face);
		    
						if( condition_found ){
							i_cond->Set(VISITED,true);
							break;
						}
					}
				}

				if( condition_found == false )
				{

					unsigned int NumberNodesInFace = lnofa[i_face];
					Condition::NodesArrayType    FaceNodes;
					FaceNodes.reserve(NumberNodesInFace);

					for(unsigned int j=1; j<=NumberNodesInFace; j++)
					{
						FaceNodes.push_back(rGeometry(lpofa(j,i_face)));
					}
		  
					Geometry<Node<3> > InsideFace(FaceNodes);
					rListOfFacesToSplit.push_back(InsideFace);
					break;
				}
			}
		}

		std::cout<<" rEdgeElements "<< rBoundaryEdgedElements.size()<<std::endl;
		std::cout<<" rEdgeConditions "<< rBoundaryEdgedConditions.size()<<std::endl;
		std::cout<<" FacesToSplit "<<rListOfFacesToSplit.size()<<std::endl;
	
		for(ModelPart::ConditionsContainerType::iterator i_cond = rBoundaryEdgedConditions.begin(); i_cond != rBoundaryEdgedConditions.end(); i_cond++)
		{
			if( i_cond->Is(VISITED) )
				i_cond->Set(VISITED,false);
			else
				std::cout<<" Some Edge Conditions Not VISITED :: something is wrong "<<std::endl;
		}
    
		KRATOS_CATCH( "" )
    }

    //*******************************************************************************************
    //*******************************************************************************************

	void GenerateNewNodes(ModelPart& rModelPart,
			std::vector<Node<3>::Pointer>& rListOfNewNodes,
			std::vector<Geometry<Node<3> > >& rListOfFacesToSplit)
    {
		KRATOS_TRY

		ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

		MeshDataTransferUtilities DataTransferUtilities;
      
		Node<3>::Pointer pNode;

		//center
		double xc = 0;
		double yc = 0;
		double zc = 0;
      
		//radius
		double radius = 0;

		//assign data to dofs
		Node<3>::DofsContainerType& ReferenceDofs = rModelPart.Nodes().front().GetDofs();

		VariablesList& VariablesList = rModelPart.GetNodalSolutionStepVariablesList();

		std::vector<double> ShapeFunctionsN;
      
		unsigned int id = MathUtils<int>::Max(ModelerUtilities::GetMaxNodeId(*(rModelPart.GetParentModelPart())), ModelerUtilities::GetMaxNodeId(rModelPart))+1;

		unsigned int size  = 0;
		unsigned int count = 0;

		for (std::vector<Geometry<Node<3> > >::iterator i_face = rListOfFacesToSplit.begin() ; i_face != rListOfFacesToSplit.end(); ++i_face)
		{
	  
			size = i_face->size();	           

			ShapeFunctionsN.resize(size);
	      
	  
			if( size == 2 )	      	
				DataTransferUtilities.CalculateCenterAndSearchRadius( (*i_face)[0].X(), (*i_face)[0].Y(), (*i_face)[1].X(), (*i_face)[1].Y(), xc,yc,radius);
	  
			if( size == 3 )
				DataTransferUtilities.CalculateCenterAndSearchRadius( (*i_face)[0].X(), (*i_face)[0].Y(), (*i_face)[0].Z(), (*i_face)[1].X(), (*i_face)[1].Y(), (*i_face)[1].Z(), (*i_face)[2].X(), (*i_face)[2].Y(), (*i_face)[2].Z(), xc,yc,zc,radius);

			//create a new node
			pNode = boost::make_shared< Node<3> >( id, xc, yc, zc );

			//giving model part variables list to the node
			pNode->SetSolutionStepVariablesList(&VariablesList);
			  
			//set buffer size
			pNode->SetBufferSize(rModelPart.GetBufferSize());
		  
			//generating the dofs
			for(Node<3>::DofsContainerType::iterator i_dof = ReferenceDofs.begin(); i_dof != ReferenceDofs.end(); i_dof++)
			{
				Node<3>::DofType& rDof = *i_dof;
				Node<3>::DofType::Pointer pNewDof = pNode->pAddDof( rDof );

				count = 0;
				for( unsigned int i = 0; i<size; i++ )
				{
					if((*i_face)[i].IsFixed(rDof.GetVariable()))
						count++;
				}
		    
				if( count == size )
					(pNewDof)->FixDof();
				else
					(pNewDof)->FreeDof();
			}
	      
			std::fill(ShapeFunctionsN.begin(), ShapeFunctionsN.end(), 1.0/double(size));            

			double alpha = 1;
			DataTransferUtilities.Interpolate( (*i_face), ShapeFunctionsN, VariablesList, pNode, alpha );

			//set flags
			pNode->Set(NEW_ENTITY);

			//set variables
			this->SetNewNodeVariables(rModelPart, pNode);
	      	          
			rListOfNewNodes.push_back(pNode);

			std::cout<<" New node added: ID "<<pNode->Id()<<" ["<<pNode->X()<<", "<<pNode->Y()<<", "<<pNode->Z()<<"]"<<std::endl;

			id++;     
	    }

		KRATOS_CATCH( "" )
    }


    //*******************************************************************************************
    //*******************************************************************************************

    void SetNewNodeVariables(ModelPart& rModelPart, Node<3>::Pointer& pNode)
    {
      KRATOS_TRY
      
		//set model part
		pNode->SetValue(MODEL_PART_NAME,rModelPart.Name());

		//set nodal_h
		pNode->FastGetSolutionStepValue(NODAL_H) = mrRemesh.Refine->CriticalSide*2.0;

		//set original position
		const array_1d<double,3>& Displacement = pNode->FastGetSolutionStepValue(DISPLACEMENT);
		pNode->X0() = pNode->X() - Displacement[0];
		pNode->Y0() = pNode->Y() - Displacement[1];
		pNode->Z0() = pNode->Z() - Displacement[2];
      
		//reset contact force
		pNode->FastGetSolutionStepValue(CONTACT_FORCE).clear();

		KRATOS_CATCH( "" )
    }

  
    //*******************************************************************************************
    //*******************************************************************************************
 
    void SetNodesToModelPart(ModelPart& rModelPart, std::vector<Node<3>::Pointer>& rListOfNewNodes)
    {
		KRATOS_TRY

		if(rListOfNewNodes.size())
		{
			//add new conditions: ( SOLID body model part )
			for(std::vector<Node<3>::Pointer>::iterator i_node = rListOfNewNodes.begin(); i_node!= rListOfNewNodes.end(); i_node++)
			{
				rModelPart.Nodes().push_back(*(i_node));
			}
			std::cout<<"+++ModelerUtilities::GetMaxNodeId(*(rModelPart.GetParentModelPart())) "<<ModelerUtilities::GetMaxNodeId(*(rModelPart.GetParentModelPart()))<<std::endl;
			std::cout<<"+++ModelerUtilities::GetMaxNodeId(*(rModelPart)) "<<ModelerUtilities::GetMaxNodeId(rModelPart)<<std::endl;
		}

		KRATOS_CATCH( "" )
    }  


    /// Assignment operator.
    RefineMeshElementsInEdgesProcess& operator=(RefineMeshElementsInEdgesProcess const& rOther);


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
                                  RefineMeshElementsInEdgesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const RefineMeshElementsInEdgesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_REFINE_MESH_ELEMENTS_IN_EDGES_PROCESS_H_INCLUDED  defined 


