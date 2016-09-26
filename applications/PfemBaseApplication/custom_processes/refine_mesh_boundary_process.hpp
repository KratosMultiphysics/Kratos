//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_REFINE_MESH_BOUNDARY_PROCESS_H_INCLUDED )
#define  KRATOS_REFINE_MESH_BOUNDARY_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "includes/model_part.h"
#include "custom_utilities/modeler_utilities.hpp"

///VARIABLES used:
//Data:     DOMAIN_LABEL(nodes)(set)
//StepData: NODAL_H, NORMAL, CONTACT_FORCE, DISPLACEMENT
//Flags:    (checked) BOUNDARY, TO_SPLIT
//          (set)     BOUNDARY(nodes), TO_ERASE(conditions), NEW_ENTITY(conditions,nodes)(set), TO_SPLIT(conditions)->locally
//          (modified)  
//          (reset)   TO_SPLIT
//(set):=(set in this process)

namespace Kratos
{

///@name Kratos Classes
///@{

/// Refine Mesh Boundary Process
/** The process labels the boundary conditions (TO_SPLIT)
    Dependencies: RemoveMeshNodesProcess.Execute()  is needed as a previous step
    
    Determines if new conditions must be inserted in boundary.
    If boundary must to be kept (CONSTRAINED), 
    New conditions will be rebuild (splitting the old ones and inserting new nodes)
    Old conditions will be erased at the end.
    
*/

class RefineMeshBoundaryProcess
  : public Process
 {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION( RefineMeshBoundaryProcess );

    typedef ModelPart::NodeType                   NodeType;
    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;
    typedef PointerVector<NodeType>        PointsArrayType;
  
    typedef PointerVectorSet<ConditionType, IndexedObject> ConditionsContainerType;
    typedef ConditionsContainerType::iterator                    ConditionIterator;
    typedef ConditionsContainerType::const_iterator      ConditionConstantIterator;
  
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RefineMeshBoundaryProcess(ModelPart& rModelPart,
			      ModelerUtilities::MeshingParameters& rRemeshingParameters,
			      int EchoLevel) 
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    {     
      mMeshId = mrRemesh.MeshId;
      mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~RefineMeshBoundaryProcess() {}


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

      if( this->mEchoLevel > 0 ){
        std::cout<<" [ REFINE BOUNDARY : "<<std::endl;
	//std::cout<<"   Nodes and Conditions : "<<mrModelPart.Nodes(mMeshId).size()<<", "<<mrModelPart.Conditions(mMeshId).size()<<std::endl;
      }


      if( mrModelPart.Name() != mrRemesh.SubModelPartName )
	std::cout<<" ModelPart Supplied do not corresponds to the Meshing Domain: ("<<mrModelPart.Name()<<" != "<<mrRemesh.SubModelPartName<<")"<<std::endl;

      
      mrRemesh.Info->InsertedConditions    = mrModelPart.NumberOfConditions(mMeshId);
      mrRemesh.Info->InsertedBoundaryNodes = mrModelPart.NumberOfNodes(mMeshId);


     //if the insert switches are activated, we check if the boundaries got too coarse
     if( (mrRemesh.Refine->RefiningOptions.Is(ModelerUtilities::REFINE_INSERT_NODES) || mrRemesh.Refine->RefiningOptions.Is(ModelerUtilities::REFINE_ADD_NODES)) && mrRemesh.Refine->RefiningOptions.Is(ModelerUtilities::REFINE_BOUNDARY) )
     {

        std::vector<NodeType::Pointer>  list_of_nodes;
        std::vector<ConditionType::Pointer> list_of_conditions;

	unsigned int conditions_size = mrModelPart.NumberOfConditions();
	
	list_of_nodes.reserve(conditions_size);
	list_of_conditions.reserve(conditions_size);
	    
	this->SelectBoundaryToRefine(mrModelPart); //conditions (TO_REFINE)  contact_nodes (TO_SPLIT)

	this->GenerateNewNodes(mrModelPart, list_of_nodes, list_of_conditions); //points (NEW_ENTITY) 
	  
	this->GenerateNewConditions(mrModelPart, list_of_nodes, list_of_conditions);// new conditions(NEW_ENTITY)  //old conditions (TO_ERASE)

	this->SetNodesToModelPart(mrModelPart, list_of_nodes);
	
	this->SetConditionsToModelPart(mrModelPart, list_of_conditions);

	//new conditions are added to model part and later added to condition model parts via (NEW_ENTITY) and (MODEL_PART_NAME)

	
     } // REFINE END;


     mrRemesh.Info->InsertedConditions    = mrModelPart.NumberOfConditions(mMeshId)-mrRemesh.Info->InsertedConditions;
     mrRemesh.Info->InsertedBoundaryNodes = mrModelPart.NumberOfNodes(mMeshId)-mrRemesh.Info->InsertedBoundaryNodes;

     if( this->mEchoLevel > 0 ){
        std::cout<<"   [ CONDITIONS ( inserted : "<<mrRemesh.Info->InsertedConditions<<" ) ]"<<std::endl;
        std::cout<<"   [ NODES      ( inserted : "<<mrRemesh.Info->InsertedBoundaryNodes<<" ) ]"<<std::endl;
        std::cout<<"   REFINE BOUNDARY ]; "<<std::endl;
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
        return "RefineMeshBoundaryProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "RefineMeshBoundaryProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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

    ModelPart& mrModelPart;
 
    ModelerUtilities::MeshingParameters& mrRemesh;

    ModelerUtilities mModelerUtilities;
  
    ModelPart::IndexType mMeshId; 

    int mEchoLevel;
  
    ///@}
    ///@name Protected Operators
    ///@{

    bool RefineOnThreshold(ConditionType::Pointer& pCondition, ProcessInfo& rCurrentProcessInfo, double& critical_size)
    {
      KRATOS_TRY
	           	      
      if( pCondition->GetValue(MASTER_ELEMENTS).size() > 0 ){
	
	Element::ElementType& MasterElement = pCondition->GetValue(MASTER_ELEMENTS).back();
	
	std::vector<double> Value;
	
	MasterElement.GetValueOnIntegrationPoints(mrRemesh.Refine->GetThresholdVariable(),Value,rCurrentProcessInfo);
		
	//calculate threshold value (plastic power)
	double threshold_value = 0;
	
	for(std::vector<double>::iterator v = Value.begin(); v!=Value.end(); v++)
	  threshold_value += *v;

	threshold_value /= double(Value.size());
	threshold_value *= MasterElement.GetGeometry().Area();
	
	//calculate condition length
	double face_size = mModelerUtilities.CalculateBoundarySize(pCondition->GetGeometry());
	
	if( threshold_value > mrRemesh.Refine->ReferenceThreshold * MasterElement.GetGeometry().Area() && face_size > critical_size )
	  return true;
      }

      return false;
      
      KRATOS_CATCH( "" )
    }

    //*******************************************************************************************
    //*******************************************************************************************


    bool RefineOnDistance(ConditionType::Pointer& pCondition, double& critical_size)
    {
      KRATOS_TRY
      	      
      //calculate condition length
      double face_size = mModelerUtilities.CalculateBoundarySize(pCondition->GetGeometry());
	
      if( face_size > critical_size )
	return true;
      
      return false;
      
      KRATOS_CATCH( "" )
    }

    //*******************************************************************************************
    //*******************************************************************************************
   
    void SelectBoundaryToRefine(ModelPart& rModelPart)
    {

      KRATOS_TRY
	     
      ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

      //***SIZES :::: parameters do define the tolerance in mesh size: 

      //THRESHOLD VALUE INSERT
      double size_for_threshold_face  = 1.50 * mrRemesh.Refine->CriticalSide; 

      //DISTANCE VALUE INSERT
      double size_for_boundary_face   = 3.50 * mrRemesh.Refine->CriticalSide;


      bool threshold_insert = false;
      bool face_size_insert = false;
      

      //LOOP TO CONSIDER ALL SUBDOMAIN CONDITIONS

      bool refine_candidate = false;
      for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.ConditionsBegin(mMeshId); i_cond!= rModelPart.ConditionsEnd(mMeshId); i_cond++)
	{

	  refine_candidate = false;
	  
	  if( mrRemesh.Options.Is(ModelerUtilities::CONSTRAINED) ){
	    if( i_cond->Is(BOUNDARY) ) //ONLY SET TO THE BOUNDARY SKIN CONDITIONS (CompositeCondition)
	      refine_candidate = true;
	    else
	      refine_candidate = false;
	  }
	  else{
	    refine_candidate = true; 
	  }


	  if( refine_candidate ){
	    if (mrRemesh.Refine->RefiningBoxSetFlag == true ){
	      refine_candidate = mModelerUtilities.CheckConditionInBox(*(i_cond.base()), *(mrRemesh.Refine->RefiningBox), rCurrentProcessInfo);
	    }
	  }


	  if( refine_candidate ){

	    threshold_insert = false;
	    face_size_insert = false;


	    //double condition_radius = 0;
	    if( i_cond->IsNot(TO_ERASE) ){

	      if ( mrRemesh.Refine->RefiningOptions.Is(ModelerUtilities::REFINE_BOUNDARY_ON_THRESHOLD) )
		threshold_insert = this->RefineOnThreshold(*(i_cond.base()), rCurrentProcessInfo, size_for_threshold_face);
	      
					
	      if ( mrRemesh.Refine->RefiningOptions.Is(ModelerUtilities::REFINE_BOUNDARY_ON_DISTANCE) )
		face_size_insert = this->RefineOnDistance(*(i_cond.base()), size_for_boundary_face);

	      if( threshold_insert || face_size_insert )
		i_cond->Set(TO_REFINE);

	    }
	    else{
	      if( this->mEchoLevel > 0 )
		std::cout<<" Condition "<<i_cond->Id()<<" TO_ERASE "<<std::endl;
	    }

	  }
        }


      KRATOS_CATCH( "" )
    }


    //*******************************************************************************************
    //*******************************************************************************************

    void GenerateNewNodes(ModelPart& rModelPart, std::vector<NodeType::Pointer>& list_of_nodes, std::vector<ConditionType::Pointer>& list_of_conditions)
    {
      KRATOS_TRY

      MeshDataTransferUtilities DataTransferUtilities;
      
      NodeType::Pointer pNode;

      //center
      double xc = 0;
      double yc = 0;
      double zc = 0;
      
      //radius
      double radius = 0;

      //assign data to dofs
      NodeType::DofsContainerType& ReferenceDofs = rModelPart.Nodes(mMeshId).front().GetDofs();

      VariablesList& VariablesList = rModelPart.GetNodalSolutionStepVariablesList();

      
      std::vector<double> ShapeFunctionsN;
      
      unsigned int id = ModelerUtilities::GetMaxNodeId(rModelPart) + 1;

      unsigned int size  = 0;
      unsigned int count = 0;
      
      for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.ConditionsBegin(mMeshId); i_cond!= rModelPart.ConditionsEnd(mMeshId); i_cond++)
	{
	  if( i_cond->Is(TO_REFINE) )
	    {	            

	      Geometry< Node<3> >& rGeometry = i_cond->GetGeometry();

	      size = rGeometry.size();	           

	      ShapeFunctionsN.resize(size);
	      
		
	      if( size == 2 )	      	
		DataTransferUtilities.CalculateCenterAndSearchRadius( rGeometry[0].X(), rGeometry[0].Y(), 
								      rGeometry[1].X(), rGeometry[1].Y(),
								      xc,yc,radius);
		
	     	      
	      if( size == 3 )
		DataTransferUtilities.CalculateCenterAndSearchRadius( rGeometry[0].X(), rGeometry[0].Y(), rGeometry[0].Z(),
								      rGeometry[1].X(), rGeometry[1].Y(), rGeometry[1].Z(),
								      rGeometry[2].X(), rGeometry[2].Y(), rGeometry[2].Z(),
								      xc,yc,zc,radius);


	      //create a new node
	      pNode = boost::make_shared< NodeType >( id, xc, yc, zc );

	      //giving model part variables list to the node
	      pNode->SetSolutionStepVariablesList(&VariablesList);
	      
	      //set buffer size
	      pNode->SetBufferSize(rModelPart.GetBufferSize());

	      //generating the dofs
	      for(Node<3>::DofsContainerType::iterator i_dof = ReferenceDofs.begin(); i_dof != ReferenceDofs.end(); i_dof++)
		{
		  NodeType::DofType& rDof = *i_dof;
		  NodeType::DofType::Pointer pNewDof = pNode->pAddDof( rDof );

		  count = 0;
		  for( unsigned int i = 0; i<size; i++ )
		    {
		      if(rGeometry[i].IsFixed(rDof.GetVariable()))
			count++;
		    }

		  if( count == size )
		    (pNewDof)->FixDof();
		  else
		    (pNewDof)->FreeDof();
		}
	      
	      std::fill(ShapeFunctionsN.begin(), ShapeFunctionsN.end(), 1.0/double(size));            

	      double alpha = 1;
	      DataTransferUtilities.Interpolate( rGeometry, ShapeFunctionsN, VariablesList, pNode, alpha );

	      //set flags
	      pNode->Set(NEW_ENTITY);
	      pNode->Set(BOUNDARY);

	      //set variables
	      this->SetNewNodeVariables(rModelPart, *(i_cond.base()), pNode);
	      	          
	      list_of_nodes.push_back(pNode);
	      list_of_conditions.push_back(*(i_cond.base()));
	      
	      id++;     
	      
	    }
			
	}
      
      KRATOS_CATCH( "" )
    }


    //*******************************************************************************************
    //*******************************************************************************************

    virtual void SetNewNodeVariables(ModelPart& rModelPart, ConditionType::Pointer& pCondition, NodeType::Pointer& pNode)
    {
      KRATOS_TRY
	           	      
      //set variables:
      Geometry< Node<3> >& rGeometry = pCondition->GetGeometry();
      
	
      //set model part
      pNode->SetValue(MODEL_PART_NAME,rModelPart.Name());

      //set nodal_h
      pNode->FastGetSolutionStepValue(NODAL_H) = mrRemesh.Refine->CriticalSide;
      
      //set normal
      noalias(pNode->FastGetSolutionStepValue(NORMAL)) = pCondition->GetValue(NORMAL);
      
      //set original position
      const array_1d<double,3>& Displacement = pNode->FastGetSolutionStepValue(DISPLACEMENT);
      pNode->X0() = pNode->X() - Displacement[0];
      pNode->Y0() = pNode->Y() - Displacement[1];
      pNode->Z0() = pNode->Z() - Displacement[2];
      
      //set contact force
      unsigned int count = 0;
      for( unsigned int i = 0; i<rGeometry.size(); i++ )
	{
	  if( norm_2(rGeometry[i].FastGetSolutionStepValue(CONTACT_FORCE)) == 0 )
	    count++;
	}
      
      if( count )
	pNode->FastGetSolutionStepValue(CONTACT_FORCE).clear();

      
      KRATOS_CATCH( "" )
    }

  
    //*******************************************************************************************
    //*******************************************************************************************
  
    void GenerateNewConditions(ModelPart& rModelPart, std::vector<NodeType::Pointer >& list_of_nodes, std::vector<ConditionType::Pointer>& list_of_conditions)
    {
      KRATOS_TRY

      std::vector<ConditionType::Pointer> list_of_new_conditions;
	
      unsigned int id = ModelerUtilities::GetMaxConditionId(rModelPart) + 1;

      ConditionType::Pointer pCondition;
      
      int size = 0;
      
      unsigned int counter = 0;
      
      for(std::vector<ConditionType::Pointer>::iterator i_cond = list_of_conditions.begin(); i_cond!= list_of_conditions.end(); i_cond++)
	{
	  Geometry< Node<3> >& rGeometry = (*i_cond)->GetGeometry();

	  size = rGeometry.size();
	  
	  PointsArrayType Nodes(size);
		
	  if( size == 2 ){	      	

	    //new condition 1
	    Nodes(0) = rGeometry(0);
	    Nodes(1) = list_of_nodes[counter];
	    
	    pCondition = (*i_cond)->Clone(id, Nodes);

	    //set flags
	    pCondition->Set(NEW_ENTITY);

	    SetNewConditionVariables((*i_cond), pCondition);

	    id++;
	    
	    list_of_new_conditions.push_back(pCondition);

	    //new condition 2
	    Nodes(0) = list_of_nodes[counter];
	    Nodes(1) = rGeometry(1);
	    
	    pCondition = (*i_cond)->Clone(id, Nodes);

	    //set flags
	    pCondition->Set(NEW_ENTITY);

	    //set variables
	    this->SetNewConditionVariables((*i_cond), pCondition);

	    id++;
	    
	    list_of_new_conditions.push_back(pCondition);
	  }
		    
	  if( size == 3 ){

	    //new condition 1
	    Nodes(0) = rGeometry(0);
	    Nodes(1) = rGeometry(1);
	    Nodes(2) = list_of_nodes[counter];
	    
	    pCondition = (*i_cond)->Clone(id, Nodes);

	    //set flags
	    pCondition->Set(NEW_ENTITY);

	    SetNewConditionVariables((*i_cond), pCondition);

	    id++;
	    
	    list_of_new_conditions.push_back(pCondition);

	    //new condition 2
	    Nodes(0) = rGeometry(1);
	    Nodes(1) = rGeometry(2);
	    Nodes(2) = list_of_nodes[counter];
	    
	    
	    pCondition = (*i_cond)->Clone(id, Nodes);

	    //set flags
	    pCondition->Set(NEW_ENTITY);

	    //set variables
	    SetNewConditionVariables((*i_cond), pCondition);

	    id++;
	    
	    list_of_new_conditions.push_back(pCondition);


	    //new condition 3
	    Nodes(0) = rGeometry(2);
	    Nodes(1) = rGeometry(0);
	    Nodes(2) = list_of_nodes[counter];
	    
	    
	    pCondition = (*i_cond)->Clone(id, Nodes);

	    //set flags
	    pCondition->Set(NEW_ENTITY);

	    //set variables
	    SetNewConditionVariables((*i_cond), pCondition);

	    id++;
	    
	    list_of_new_conditions.push_back(pCondition);

	    
	  }

	  // once the condition is refined set to erase
	  (*i_cond)->Set(TO_ERASE);
	  (*i_cond)->Set(TO_REFINE, false);

	  ConditionsContainerType& ChildrenConditions = (*i_cond)->GetValue(CHILDREN_CONDITIONS);

	  for (ConditionConstantIterator cn = ChildrenConditions.begin() ; cn != ChildrenConditions.end(); ++cn)
	    {
	      cn->Set(TO_ERASE);
	    }

	  
	  counter++;
	    
	}

      //update the list of old conditions with the list of new conditions
      list_of_conditions = list_of_new_conditions;
	      
      KRATOS_CATCH( "" )
    }


      
    //*******************************************************************************************
    //*******************************************************************************************

    virtual void SetNewConditionVariables(ConditionType::Pointer& pOldCondition, ConditionType::Pointer& pNewCondition)
    {
      KRATOS_TRY
	           	      
      //set variables
      pNewCondition->SetValue( MASTER_NODES          , pOldCondition->GetValue(MASTER_NODES)         );
      pNewCondition->SetValue( NORMAL                , pOldCondition->GetValue(NORMAL)               );
      pNewCondition->SetValue( CAUCHY_STRESS_VECTOR  , pOldCondition->GetValue(CAUCHY_STRESS_VECTOR) );
      pNewCondition->SetValue( DEFORMATION_GRADIENT  , pOldCondition->GetValue(DEFORMATION_GRADIENT) );

      KRATOS_CATCH( "" )
    }
      

    //*******************************************************************************************
    //*******************************************************************************************
 
    void SetNodesToModelPart(ModelPart& rModelPart, std::vector<NodeType::Pointer>& list_of_nodes)
    {
      KRATOS_TRY

      if(list_of_nodes.size()){

	//add new conditions: ( SOLID body model part )
	for(std::vector<NodeType::Pointer>::iterator i_node = list_of_nodes.begin(); i_node!= list_of_nodes.end(); i_node++)
	  {
	    rModelPart.Nodes().push_back(*(i_node));
	  }
	
      }

      KRATOS_CATCH( "" )
    }
    //*******************************************************************************************
    //*******************************************************************************************
 
    void SetConditionsToModelPart(ModelPart& rModelPart, std::vector<ConditionType::Pointer>& list_of_conditions)
    {
      KRATOS_TRY

      if(list_of_conditions.size()){

	//clear erased conditions: ( SOLID body model part )
	this->CleanModelPartConditions(rModelPart);
  
	//add new conditions: ( SOLID body model part )
	for(std::vector<ConditionType::Pointer>::iterator i_cond = list_of_conditions.begin(); i_cond!= list_of_conditions.end(); i_cond++)
	  {
	    rModelPart.Conditions().push_back(*(i_cond));
	  }
	
      }

      //renumerate conditions
      // unsigned int id=1;
      // for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.ConditionsBegin(); i_cond!= rModelPart.ConditionsEnd(); i_cond++)
      // 	{
      // 	  i_cond->SetId(id);
      // 	  id++;
      // 	}
      
      
      KRATOS_CATCH( "" )
    }
       

    //*******************************************************************************************
    //*******************************************************************************************

    void CleanModelPartConditions(ModelPart& rModelPart)
    {
      KRATOS_TRY

	
      //clean old conditions (TO_ERASE) and add new conditions (NEW_ENTITY)
      ModelPart::ConditionsContainerType PreservedConditions;
      PreservedConditions.reserve(rModelPart.Conditions().size());
      PreservedConditions.swap(rModelPart.Conditions());
      
      for(ModelPart::ConditionsContainerType::iterator i_cond = PreservedConditions.begin(); i_cond!= PreservedConditions.end(); i_cond++)
    	{
    	  if(i_cond->IsNot(TO_ERASE))
    	    rModelPart.Conditions().push_back(*(i_cond.base()));
    	}
	  	
      KRATOS_CATCH( "" )
    }
  
    ///@}
    ///@name Protected Operations
    ///@{


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
    ///@name Private Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{
  
    ///@}
    ///@name Private Operations
    ///@{

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
    RefineMeshBoundaryProcess& operator=(RefineMeshBoundaryProcess const& rOther);


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
                                  RefineMeshBoundaryProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const RefineMeshBoundaryProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_REFINE_MESH_BOUNDARY_PROCESS_H_INCLUDED  defined 


