//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_MODEL_START_END_MESHING_PROCESS_H_INCLUDED )
#define  KRATOS_MODEL_START_END_MESHING_PROCESS_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_utilities/modeler_utilities.hpp"
#include "custom_processes/nodal_neighbours_search_process.hpp"
#include "custom_processes/build_model_part_boundary_process.hpp"

#include "pfem_base_application_variables.h"

///VARIABLES used:
//Data:     
//StepData: 
//Flags:    (checked) 
//          (set)     
//          (modified)  
//          (reset)   


namespace Kratos
{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{
  typedef  ModelPart::ConditionType                                ConditionType;
  typedef  ModelPart::NodesContainerType                      NodesContainerType;
  typedef  ModelPart::ElementsContainerType                ElementsContainerType;
  typedef  ModelPart::ConditionsContainerType            ConditionsContainerType;
  typedef  ModelPart::MeshType::GeometryType::PointsArrayType    PointsArrayType;

  typedef PointerVectorSet<ConditionType, IndexedObject> ConditionsContainerType;
  typedef ConditionsContainerType::iterator                    ConditionIterator;
  typedef ConditionsContainerType::const_iterator      ConditionConstantIterator;
  
  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{

  ///@}
  ///@name Kratos Classes
  ///@{

  /// Short class definition.
  /** Detail class definition.
   */
  class ModelStartEndMeshingProcess
    : public Process
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ModelStartEndMeshingProcess
    KRATOS_CLASS_POINTER_DEFINITION( ModelStartEndMeshingProcess );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ModelStartEndMeshingProcess(ModelPart& rMainModelPart,
				Flags Options,
				int EchoLevel = 0)
      : mrMainModelPart(rMainModelPart)
    { 
      mOptions   = Options;
      mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~ModelStartEndMeshingProcess()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
      Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    virtual void Execute()
    {

    };

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    virtual void ExecuteInitialize()
    {

      KRATOS_TRY

      //Sort Conditions
      this->SortModelPartConditions();

      KRATOS_CATCH(" ")
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    virtual void ExecuteFinalize()
    {
      KRATOS_TRY
      
      // Restore meshes coherency
      //this->BuildModelPartMeshes();
      this->BuildModelPartStructure();

      // Perform searches for next processes (i.e. contact search)
      this->PerformModelSearches();

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
    virtual std::string Info() const
    {
      return "ModelStartEndMeshingProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
      rOStream << "ModelStartEndMeshingProcess";
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

    ModelPart& mrMainModelPart;

    Flags mOptions;

    int mEchoLevel;


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    void SortModelPartConditions()
    {

      KRATOS_TRY

      //Sort Conditions
      unsigned int consecutive_index = 1;
      for(ModelPart::ConditionsContainerType::iterator ic = mrMainModelPart.ConditionsBegin(); ic!=mrMainModelPart.ConditionsEnd(); ic++)
	(ic)->SetId(consecutive_index++);

      mrMainModelPart.Conditions().Sort();
      // mrMainModelPart.Conditions().Unique();


      KRATOS_CATCH(" ")

    }


    //*******************************************************************************************
    //*******************************************************************************************

    void SortModelPartElements()
    {

      KRATOS_TRY


      // Sort Elements
      unsigned int consecutive_index = 1;
      for(ModelPart::ElementsContainerType::iterator ie = mrMainModelPart.ElementsBegin(); ie!=mrMainModelPart.ElementsEnd(); ie++)
	ie->SetId(consecutive_index++);

      mrMainModelPart.Elements().Sort();
      // mrMainModelPart.Elements().Unique();

      KRATOS_CATCH(" ")

     }


    //*******************************************************************************************
    //*******************************************************************************************
    
    void SortModelPartNodes()
    {

      KRATOS_TRY


      //Sort Nodes, set STRUCTURE nodes at end
      unsigned int consecutive_index = 1;
      unsigned int reverse_index = mrMainModelPart.Nodes().size();
      for(ModelPart::NodesContainerType::iterator in = mrMainModelPart.NodesBegin(); in!=mrMainModelPart.NodesEnd(); in++)
	{
	  if(in->IsNot(STRUCTURE) ){
	    in->SetId(consecutive_index++);
	  }
	  else{
	    in->SetId(reverse_index--);
	  }
	}

      mrMainModelPart.Nodes().Sort();
      // mrMainModelPart.Nodes().Unique();

      KRATOS_CATCH(" ")

     }

    
    //*******************************************************************************************
    //*******************************************************************************************
    
    void BuildModelPartMeshes()
    {

      KRATOS_TRY

      //Once all meshes are build, the main mesh Id=0 must be reassigned
      unsigned int NumberOfMeshes = mrMainModelPart.NumberOfMeshes();

      if(NumberOfMeshes>1){
	this->BuildTotalMesh(mrMainModelPart, mEchoLevel);
      }
      else{
	this->CleanMeshFlags(mrMainModelPart,0);
      }


      KRATOS_CATCH(" ")

    }


    //*******************************************************************************************
    //*******************************************************************************************
    
    void BuildModelPartStructure()
    {

      KRATOS_TRY

      //Once all model parts are build, the main model part must be reconstructed coherently
      unsigned int NumberOfSubModelParts=mrMainModelPart.NumberOfSubModelParts();

      if(NumberOfSubModelParts>0){
	this->BuildTotalModelPart(mrMainModelPart, mEchoLevel);
      }
      else{
	this->CleanMeshFlags(mrMainModelPart,0);
      }


      KRATOS_CATCH(" ")

    }

    //*******************************************************************************************
    //*******************************************************************************************
    
    virtual void PerformModelSearches()
    {

      KRATOS_TRY
       
	//NODAL NEIGHBOURS SEARCH
	NodalNeighboursSearchProcess FindNeighbours(mrMainModelPart);
	FindNeighbours.Execute();

	//NODAL_H SEARCH	    
	//FindNodalHProcess FindNodalH(mrMainModelPart);
	//FindNodalH.Execute();

	//CONDITIONS MASTER_ELEMENTS and MASTER_NODES SEARCH
	BuildModelPartBoundaryProcess BuildBoundaryProcess(mrMainModelPart, mrMainModelPart.Name(), mEchoLevel);
	BuildBoundaryProcess.SearchConditionMasters(mrMainModelPart);

	//BOUNDARY NORMALS SEARCH and SHRINKAGE FACTOR
	BoundaryNormalsCalculationUtilities BoundaryComputation;
	BoundaryComputation.CalculateWeightedBoundaryNormals(mrMainModelPart, mEchoLevel);   


      KRATOS_CATCH(" ")

     }


    //*******************************************************************************************
    //*******************************************************************************************
    
    void ChangeTipElements()
    {

      KRATOS_TRY

	
	//CHANGE ELEMENTS WITH ALL NODES IN BOUNDARY  (TIP ELEMENTS)
	// for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
	//   {
		
	//     if(Meshes[MeshId].Is( ModelerUtilities::REMESH )){
	//       ChangeTipElementsUtilities TipElements;
	//       if(mpMeshingVariables->AvoidTipElementsFlag){
	// 	TipElements.SwapDiagonals(mrMainModelPart,MeshId);
	//       }
	//     }
	//   }	
	
      KRATOS_CATCH(" ")

     }

    //*******************************************************************************************
    //*******************************************************************************************

    virtual void BuildTotalModelPart(ModelPart& rModelPart, int EchoLevel)
    {

      KRATOS_TRY

      //Mesh Id=0
      
      if( EchoLevel > 0 )
	std::cout<<"   [ START MODEL PART ["<<rModelPart.Name()<<"] [Elems=:"<<rModelPart.NumberOfElements()<<"|Nodes="<<rModelPart.NumberOfNodes()<<"|Conds="<<rModelPart.NumberOfConditions()<<"] ] "<<std::endl;      
      
      rModelPart.Nodes().clear();
      rModelPart.Elements().clear();

      //contact conditions are located on Mesh_0
      ModelPart::ConditionsContainerType PreservedConditions;
      
      unsigned int nodeId=1;
      unsigned int elemId=1;
      unsigned int condId=1;

      this->BuildBodyModelParts(rModelPart, PreservedConditions, nodeId, elemId, condId);

      this->BuildBoundaryModelParts(rModelPart,PreservedConditions, nodeId, elemId, condId);
           
      this->BuildContactModelParts(rModelPart, PreservedConditions, nodeId, elemId, condId);
          
      //now set new conditions
      rModelPart.Conditions().swap(PreservedConditions);
      
      //Sort
      rModelPart.Nodes().Sort();
      rModelPart.Elements().Sort();
      rModelPart.Conditions().Sort();     
      
      //Unique
      rModelPart.Nodes().Unique();
      rModelPart.Elements().Unique();
      rModelPart.Conditions().Unique();
      
      //Sort Again to have coherent numeration for nodes (mesh with shared nodes)
      unsigned int consecutive_index = 1;
      for(ModelPart::NodesContainerType::iterator in = rModelPart.NodesBegin() ; in != rModelPart.NodesEnd() ; in++)
	in->SetId(consecutive_index++);
           
      this->BuildComputingDomain(rModelPart, EchoLevel);

      if( EchoLevel > 0 )
	std::cout<<"   [ END MODEL PART ["<<rModelPart.Name()<<"] [Elems=:"<<rModelPart.NumberOfElements()<<"|Nodes="<<rModelPart.NumberOfNodes()<<"|Conds="<<rModelPart.NumberOfConditions()<<"] ] "<<std::endl;      
      

      KRATOS_CATCH(" ")

    }


    //*******************************************************************************************
    //*******************************************************************************************
 
    void BuildBodyModelParts(ModelPart& rModelPart, ModelPart::ConditionsContainerType& rPreservedConditions, unsigned int& rNodeId, unsigned int& rElemId, unsigned int& rCondId)
    {
      KRATOS_TRY
	     
      //Add Fluid Bodies modelparts to main modelpart  flags: ( FLUID / NOT_ACTIVE / NOT_BOUNDARY )
      //Add Solid Bodies modelparts to main modelpart  flags: ( SOLID / NOT_ACTIVE / NOT_BOUNDARY )
      //Add Rigid Bodies modelparts to main modelpart  flags: ( RIGID / NOT_ACTIVE / NOT_BOUNDARY )

      const array_1d<double,3> ZeroNormal(3,0.0);
	
      for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin() ; i_mp!=rModelPart.SubModelPartsEnd(); i_mp++)
	{

	  bool add_to_main_model_part = false;
	  
	  if( i_mp->IsNot(ACTIVE) && i_mp->IsNot(BOUNDARY) ){ //only the domains (no computing, no boundary)

	    if( i_mp->Is(SOLID) || i_mp->Is(FLUID) || i_mp->Is(RIGID) || i_mp->IsNot(CONTACT) )
	      add_to_main_model_part = true;
	  }
	  
	  
	  if( add_to_main_model_part  ){ 

    
	    if( mEchoLevel > 0 )
	      std::cout<<"    [ SUBMODEL PART ["<<i_mp->Name()<<"] [Elems="<<i_mp->NumberOfElements()<<"|Nodes="<<i_mp->NumberOfNodes()<<"|Conds="<<i_mp->NumberOfConditions()<<"] ";//Clean Nodes when redefining the main model part:
	    ModelPart::NodesContainerType temporal_nodes;
	    temporal_nodes.reserve(i_mp->Nodes().size());
	    temporal_nodes.swap(i_mp->Nodes());

	    
	    if( !i_mp->NumberOfElements() ){ // usually rigid domains
	      
	      for(ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; i_node++)
		{
		  if( i_node->IsNot(TO_ERASE) ){
		    (i_mp->Nodes()).push_back(*(i_node.base()));
		    (rModelPart.Nodes()).push_back(*(i_node.base()));	
		    rModelPart.Nodes().back().SetId(rNodeId);
		    rNodeId+=1;
		  }
		}
	    }
	    else{
	    
	      for(ModelPart::ElementsContainerType::iterator i_elem = i_mp->ElementsBegin() ; i_elem != i_mp->ElementsEnd() ; i_elem++)
		{
		  if( i_elem->IsNot(TO_ERASE) ){ //at this point any element must be TO_ERASE		
		  
		    PointsArrayType& vertices=i_elem->GetGeometry().Points();
		    for(unsigned int i=0; i<vertices.size(); i++)
		      {
			vertices[i].Set(BLOCKED);
		      }
		  
		    (rModelPart.Elements()).push_back(*(i_elem.base()));	
		    rModelPart.Elements().back().SetId(rElemId);
		    rElemId+=1;
		  
		  }
		  else{
		    std::cout<<" ELEMENT TO_ERASE must be RELEASED in a previous stage "<<std::endl;
		  }
		}

	      for(ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; i_node++)
		{
		  if(i_node->Is(BLOCKED) || i_node->Is(RIGID)){ //all nodes belonging to an element had been blocked previously
		      
		    i_node->Reset(ISOLATED);   //reset isolated
		    i_node->Reset(TO_REFINE);  //reset if was labeled to refine (to not duplicate boundary conditions)
		    i_node->Reset(BLOCKED); 
		      
		    if( i_node->IsNot(TO_ERASE) ){

		      (i_mp->Nodes()).push_back(*(i_node.base()));
		      (rModelPart.Nodes()).push_back(*(i_node.base()));	
		      rModelPart.Nodes().back().SetId(rNodeId);
		      rNodeId+=1;
			
		    }

		  }
		  else{

		    i_node->Set(ISOLATED);
		    i_node->Reset(TO_REFINE);  //reset if was labeled to refine (to not duplicate boundary conditions)
		    i_node->Reset(BLOCKED); 

		    if( mOptions.Is(ModelerUtilities::KEEP_ISOLATED_NODES) && i_node->IsNot(TO_ERASE) ){

		      (i_mp->Nodes()).push_back(*(i_node.base()));
		      (rModelPart.Nodes()).push_back(*(i_node.base()));	
		      rModelPart.Nodes().back().SetId(rNodeId);
		      rNodeId+=1;

		    }
		    
		  }
		
		  if(i_node->Is(BOUNDARY)){

		    if(i_node->IsNot(RIGID))
		      i_node->Set(FREE_SURFACE);

		  }
		  else{

		    if( i_node->SolutionStepsDataHas(CONTACT_FORCE) )
		      noalias(i_node->GetSolutionStepValue(CONTACT_FORCE)) = ZeroNormal;
		  }
				
		}

	    }

	    for(ModelPart::ConditionsContainerType::iterator i_cond = i_mp->ConditionsBegin() ; i_cond != i_mp->ConditionsEnd() ; i_cond++)
	      {
		
		if( i_cond->IsNot(TO_ERASE) ){
		  i_cond->Reset(TO_REFINE);  //reset if was labeled to refine (to not duplicate boundary conditions)
		  rPreservedConditions.push_back(*(i_cond.base()));
		  rPreservedConditions.back().SetId(rCondId);
		  rCondId+=1;	
		}
		
	      }
	    
	    if( mEchoLevel > 0 )
	      std::cout<<" / [Elems="<<i_mp->NumberOfElements()<<"|Nodes="<<i_mp->NumberOfNodes()<<"|Conds="<<i_mp->NumberOfConditions()<<"] ] "<<std::endl;

	  }
	  
	}
	
      KRATOS_CATCH(" ")
    }

    
    //*******************************************************************************************
    //*******************************************************************************************
 
    void BuildBoundaryModelParts(ModelPart& rModelPart, ModelPart::ConditionsContainerType& rPreservedConditions, unsigned int& rNodeId, unsigned int& rElemId, unsigned int& rCondId)
    {
      
      KRATOS_TRY

      unsigned int body_model_part_conditions = rPreservedConditions.size();
      
      if(body_model_part_conditions > 0){

	//add new conditions: ( BOUNDARY model parts )
	for(ModelPart::SubModelPartIterator i_mp= mrMainModelPart.SubModelPartsBegin() ; i_mp!=mrMainModelPart.SubModelPartsEnd(); i_mp++)
	  {
	    if( i_mp->Is(BOUNDARY) ){ //boundary model part

	      if( mEchoLevel > 0 )
		std::cout<<"    [ SUBMODEL PART ["<<i_mp->Name()<<"]  initial [Elems="<<i_mp->NumberOfElements()<<"|Nodes="<<i_mp->NumberOfNodes()<<"|Conds="<<i_mp->NumberOfConditions()<<"] ]"<<std::endl;

	      this->CleanModelPartConditions(*i_mp);
		
	      for(ModelPart::ConditionsContainerType::iterator i_cond = rPreservedConditions.begin(); i_cond!= rPreservedConditions.end(); i_cond++)
		{   
		  ConditionsContainerType& ChildrenConditions = i_cond->GetValue(CHILDREN_CONDITIONS);

		  //this conditions are cloned, then the id has no coherence, must be renumbered at the end of the assignation
		  for (ConditionConstantIterator cn = ChildrenConditions.begin() ; cn != ChildrenConditions.end(); ++cn)
		    {
		      if( cn->GetValue(MODEL_PART_NAME) == i_mp->Name() ){
			i_mp->Conditions().push_back(*(cn.base()));

			if( i_cond->Is(NEW_ENTITY) ){
			  for(unsigned int i=0; i<i_cond->GetGeometry().size(); i++)
			    {
			      if( i_cond->GetGeometry()[i].Is(NEW_ENTITY) ){
				(i_mp->Nodes()).push_back(i_cond->GetGeometry()(i));
				//i_cond->GetGeometry()[i].Reset(NEW_ENTITY); //reset if was new 
			      }
			    }
			  //i_cond->Reset(NEW_ENTITY); //reset here if the condition is inserted
			}
			  
		      }

		    }
		  
		}
	    }	      

	  } 
	  
      }

      // Set new nodes to the dirichlet sub model parts (works in 2D. not shure in 3D).  Must be reviewed
      for(ModelPart::SubModelPartIterator i_model_part = mrMainModelPart.SubModelPartsBegin() ; i_model_part != mrMainModelPart.SubModelPartsEnd();  i_model_part++)
	{
	  if( i_model_part->IsNot(BOUNDARY) &&  i_model_part->IsNot(ACTIVE) && i_model_part->IsNot(RIGID) ){
	    
	    for(ModelPart::NodesContainerType::iterator i_node = i_model_part->NodesBegin() ; i_node != i_model_part->NodesEnd(); i_node++ )
	      {

		if( i_node->Is(BOUNDARY) && i_node->Is(NEW_ENTITY) ){

		  // Generate a list of neighbour nodes
		  unsigned int NodeId = i_node->Id();
		  
		  std::vector<int> list_of_neighbour_nodes; 

		  for( ModelPart::ConditionsContainerType::iterator j_cond = rPreservedConditions.begin(); j_cond != rPreservedConditions.end(); j_cond++)
		    {

		      bool node_belongs_to_condition = false; 
		      Geometry< Node<3 > > & rjGeom = j_cond->GetGeometry();

		      if ( j_cond->Is(NEW_ENTITY) ){
			for ( unsigned int j = 0; j < rjGeom.size() ; j++) {
			  if ( rjGeom[j].Id() == NodeId) {
			    node_belongs_to_condition = true; 
			    break;
			  }
			}

			if (node_belongs_to_condition){
			  for (unsigned int j = 0; j < rjGeom.size() ; j++)  {
			    list_of_neighbour_nodes.push_back( rjGeom[j].Id() );
			  }
			}

		      }
		    }
		  
		  if(list_of_neighbour_nodes.size() == 0){
		    std::cout << " something wierd, this new node does not have any new neighbour: " << NodeId << std::endl;
		    // aqui falta un continue o algu ( no un break)
		    continue;
		  }

		  
		  // unique and sort
		  std::sort(list_of_neighbour_nodes.begin(), list_of_neighbour_nodes.end() );
		  std::vector<int>::iterator new_end = std::unique( list_of_neighbour_nodes.begin(), list_of_neighbour_nodes.end() );
		  list_of_neighbour_nodes.resize( std::distance( list_of_neighbour_nodes.begin(), new_end) ); 

		  for (ModelPart::SubModelPartIterator i_mp = mrMainModelPart.SubModelPartsBegin(); i_mp != mrMainModelPart.SubModelPartsEnd(); i_mp++)
		    {
		      if ( i_mp->Is(BOUNDARY) && i_mp->IsNot(CONTACT) && (i_mp->NumberOfConditions() == 0) )
			{
			  unsigned int counter = 0;

			  for (unsigned int ii = 0; ii < list_of_neighbour_nodes.size(); ii++)
			    {
			      unsigned int target = list_of_neighbour_nodes[ii];

			      for (ModelPart::NodesContainerType::iterator i_node = i_mp->NodesBegin(); i_node != i_mp->NodesEnd(); i_node++)
				{
				  if ( i_node->Id() == target) {
				    counter++;
				  }
				}
			    }

			  if ( counter == list_of_neighbour_nodes.size()-1)
			    i_mp->Nodes().push_back( *(i_node.base() ) );
			}
		    }

		}
	      }
	    

	  }

	}

      //reset NEW_ENTITIES in conditions
      for( ModelPart::ConditionsContainerType::iterator j_cond = rPreservedConditions.begin(); j_cond != rPreservedConditions.end(); j_cond++)
	{
	  j_cond->Reset(NEW_ENTITY);
	  for(unsigned int j=0; j<j_cond->GetGeometry().size(); j++)
	    {
	      if( j_cond->GetGeometry()[j].Is(NEW_ENTITY) ){
		j_cond->GetGeometry()[j].Reset(NEW_ENTITY); //reset if was new
	      }
	    }
	}

      
      //add new nodes: ( BOUNDARY model parts ) and remove erased nodes
      for(ModelPart::SubModelPartIterator i_mp= mrMainModelPart.SubModelPartsBegin() ; i_mp!=mrMainModelPart.SubModelPartsEnd(); i_mp++)
	{
	  if( i_mp->Is(BOUNDARY) ){ //boundary model part

	    ModelPart::NodesContainerType temporal_nodes;
	    temporal_nodes.reserve(i_mp->Nodes().size());
	    temporal_nodes.swap(i_mp->Nodes());

	    for(ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; i_node++)
	      {
		if( i_node->IsNot(TO_ERASE) )
		  (i_mp->Nodes()).push_back(*(i_node.base()));
	      }

	    if( mEchoLevel > 0 )
	      std::cout<<"    [ SUBMODEL PART ["<<i_mp->Name()<<"]  final [Elems="<<i_mp->NumberOfElements()<<"|Nodes="<<i_mp->NumberOfNodes()<<"|Conds="<<i_mp->NumberOfConditions()<<"] ] "<<std::endl;
	  }
	}


      //add boundary domain conditions to preserved conditions
      for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin() ; i_mp!=rModelPart.SubModelPartsEnd(); i_mp++)
	{
	  if( i_mp->Is(BOUNDARY) ){ //boundary model part

	    for(ModelPart::ConditionsContainerType::iterator i_cond = i_mp->ConditionsBegin() ; i_cond != i_mp->ConditionsEnd() ; i_cond++)
	      {
		if( i_cond->IsNot(TO_ERASE) ){
		  i_cond->Reset(NEW_ENTITY); //reset here if the condition is inserted
		  rPreservedConditions.push_back(*(i_cond.base()));
		  rPreservedConditions.back().SetId(rCondId);
		  rCondId+=1;	
		}
	      }
	  }
	  
	}      
	
      KRATOS_CATCH( "" )
    }



    //*******************************************************************************************
    //*******************************************************************************************
 
    void BuildContactModelParts(ModelPart& rModelPart, ModelPart::ConditionsContainerType& rPreservedConditions, unsigned int& rNodeId, unsigned int& rElemId, unsigned int& rCondId)
    {
      KRATOS_TRY

      //Add Contact modelparts to main modelpart  flags: ( CONTACT ) in contact model parts keep only nodes and contact conditions // after that a contact search will be needed

      //if contact condition has the same geometry size as an elements printing ids will coincide,
      //renumber conditions with rElemId instead of rCondId :: in order to ensure it check maximun and apply it
      unsigned int rContactId = rCondId;
      if( rElemId > rCondId )
	rContactId = rElemId;
       

      for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin() ; i_mp!=rModelPart.SubModelPartsEnd(); i_mp++)
	{
  
	  if( i_mp->Is(CONTACT) ){ //keep only contact conditions
	    
	    if( mEchoLevel > 0 )
	      std::cout<<"    [ SUBMODEL PART ["<<i_mp->Name()<<"] [Elems="<<i_mp->NumberOfElements()<<"|Nodes="<<i_mp->NumberOfNodes()<<"|Conds="<<i_mp->NumberOfConditions()<<"] ";

	    i_mp->Elements().clear();
	    	    
	    //Clean Nodes when redefining the main model part:
	    ModelPart::NodesContainerType temporal_nodes;
	    temporal_nodes.reserve(i_mp->Nodes().size());
	    temporal_nodes.swap(i_mp->Nodes());
	  
	    
	    for(ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; i_node++)
	      {
		if( i_node->IsNot(TO_ERASE) ){
		  (i_mp->Nodes()).push_back(*(i_node.base()));
		}
	      }
	    
	    ModelPart::ConditionsContainerType temporal_conditions;
	    temporal_conditions.reserve(i_mp->Conditions().size());
	    temporal_conditions.swap(i_mp->Conditions());
	    	    
	    for(ModelPart::ConditionsContainerType::iterator i_cond = temporal_conditions.begin() ; i_cond != temporal_conditions.end() ; i_cond++)
	      {
		
		if( i_cond->Is(CONTACT) ){ //keep only contact conditions

		  if( i_cond->IsNot(TO_ERASE) ){ //it can not be to erase
		    
		    (i_mp->Conditions()).push_back(*(i_cond.base()));
		    rPreservedConditions.push_back(*(i_cond.base()));
		    rPreservedConditions.back().SetId(rContactId);
		    rContactId+=1;
		    
		  }
		  
		}
		
	      }	  

	    if( mEchoLevel > 0 )
	      std::cout<<" / [Elems="<<i_mp->NumberOfElements()<<"|Nodes="<<i_mp->NumberOfNodes()<<"|Conds="<<i_mp->NumberOfConditions()<<"] ] "<<std::endl;


	  }
	  
	}

      
      if( rElemId > rCondId )
	rElemId = rContactId;
      else
	rCondId = rContactId;
      

      KRATOS_CATCH(" ")
    }

    
    //*******************************************************************************************
    //*******************************************************************************************

    void CleanModelPartConditions(ModelPart& rModelPart)
    {
      
      KRATOS_TRY

      if( rModelPart.Is(BOUNDARY) )
	rModelPart.Conditions().clear();
      
      //clean old conditions (TO_ERASE) and add new conditions (NEW_ENTITY)
      // ModelPart::ConditionsContainerType PreservedConditions;
      // PreservedConditions.reserve(rModelPart.Conditions().size());
      // PreservedConditions.swap(rModelPart.Conditions());
      
      // for(ModelPart::ConditionsContainerType::iterator i_cond = PreservedConditions.begin(); i_cond!= PreservedConditions.end(); i_cond++)
      // 	{
      // 	  if(i_cond->IsNot(TO_ERASE))
      // 	    rModelPart.Conditions().push_back(*(i_cond.base()));
      // 	}


	
      KRATOS_CATCH( "" )
    }
    

    //*******************************************************************************************
    //*******************************************************************************************

    virtual void BuildComputingDomain (ModelPart& rModelPart, int EchoLevel)
    {
      KRATOS_TRY
	
      std::string ComputingModelPartName;
      for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin() ; i_mp!=rModelPart.SubModelPartsEnd(); i_mp++)
	{
	  if( i_mp->Is(ACTIVE) ){ //computing_domain
	    ComputingModelPartName = i_mp->Name();
	  }
	}


      ModelPart& rComputingModelPart = rModelPart.GetSubModelPart(ComputingModelPartName);

      rComputingModelPart.Nodes().clear();
      rComputingModelPart.Elements().clear();
      rComputingModelPart.Conditions().clear();

      //add all needed computing entities (elements, nodes, conditions)

      for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin() ; i_mp!=rModelPart.SubModelPartsEnd(); i_mp++)
	{
	  if( (i_mp->IsNot(BOUNDARY) && i_mp->IsNot(ACTIVE)) || (i_mp->Is(RIGID)) ){ 

	    if( i_mp->IsNot(CONTACT)  ){
	      
	      for(ModelPart::NodesContainerType::iterator i_node = i_mp->NodesBegin() ; i_node != i_mp->NodesEnd() ; i_node++)
		{
		  (rComputingModelPart.Nodes()).push_back(*(i_node.base()));
		}


	      for(ModelPart::ConditionsContainerType::iterator i_cond = i_mp->ConditionsBegin() ; i_cond != i_mp->ConditionsEnd() ; i_cond++)
		{
		  (rComputingModelPart.Conditions()).push_back(*(i_cond.base()));
		}

	      for(ModelPart::ElementsContainerType::iterator i_elem = i_mp->ElementsBegin() ; i_elem != i_mp->ElementsEnd() ; i_elem++)
		{
		  (rComputingModelPart.Elements()).push_back(*(i_elem.base()));	
		}
	      
	    }
	   
	  }
	}

      //add all contact conditions 
      
      for(ModelPart::SubModelPartIterator i_mp= rModelPart.SubModelPartsBegin() ; i_mp!=rModelPart.SubModelPartsEnd(); i_mp++)
	{
	  if( i_mp->Is(CONTACT)  ){
	      
	    for(ModelPart::ConditionsContainerType::iterator i_cond = i_mp->ConditionsBegin() ; i_cond != i_mp->ConditionsEnd() ; i_cond++)
	      {
		if( i_cond->Is(CONTACT) )
		  (rComputingModelPart.Conditions()).push_back(*(i_cond.base()));
	      }
	      
	  }
	   
	}
      
      
      //Sort
      // rComputingModelPart.Nodes().Sort();
      // rComputingModelPart.Elements().Sort();
      // rComputingModelPart.Conditions().Sort();     
      
      //Unique
      // rComputingModelPart.Nodes().Unique();
      // rComputingModelPart.Elements().Unique();
      // rComputingModelPart.Conditions().Unique();
      
      if( EchoLevel > 0 )
	std::cout<<"    [ SUBMODEL PART ["<<rComputingModelPart.Name()<<"] [Elems="<<rComputingModelPart.NumberOfElements()<<"|Nodes="<<rComputingModelPart.NumberOfNodes()<<"|Conds="<<rComputingModelPart.NumberOfConditions()<<"] ] "<<std::endl;

 
      KRATOS_CATCH(" ")
    }

    //*******************************************************************************************
    //*******************************************************************************************

    void BuildTotalMesh (ModelPart& rModelPart, int EchoLevel)
    {

      KRATOS_TRY

      //Mesh Id=0
      
      if( EchoLevel > 0 )
	std::cout<<"   [ START MESH [Id=0] [Elems=:"<<rModelPart.NumberOfElements()<<"|Nodes="<<rModelPart.NumberOfNodes()<<"|Conds="<<rModelPart.NumberOfConditions()<<"] ] "<<std::endl;      
      
      rModelPart.Nodes().clear();
      rModelPart.Elements().clear();

      //contact conditions are located on Mesh_0
      ModelPart::ConditionsContainerType KeepConditions;
      

      //std::cout<<" [ Number of Meshes "<<rModelPart.GetMeshes().size()-1<<" ]"<<std::endl;
      
      unsigned int nodeId=1;
      unsigned int elemId=1;
      unsigned int condId=1;
      
      unsigned int start=0;
      unsigned int NumberOfMeshes=rModelPart.NumberOfMeshes();
      if(NumberOfMeshes>1) 
	start=1;
      

      for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
	{
	  if( EchoLevel > 0 )
	    std::cout<<"    [ CHILD MESH [Id:"<<MeshId<<"] [Elems="<<rModelPart.NumberOfElements(MeshId)<<"|Nodes="<<rModelPart.NumberOfNodes(MeshId)<<"|Conds="<<rModelPart.NumberOfConditions(MeshId)<<"] ] "<<std::endl;


	  //Clean Nodes when redefining the total mesh:
	  const array_1d<double,3> ZeroNormal(3,0.0);
	  ModelPart::NodesContainerType temporal_nodes;
	  temporal_nodes.reserve(rModelPart.Nodes(MeshId).size());
	  temporal_nodes.swap(rModelPart.Nodes(MeshId));

	  if (!rModelPart.NumberOfElements(MeshId)){	  
	    for(ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; i_node++)
	      {
     		(rModelPart.Nodes(MeshId)).push_back(*(i_node.base()));
		(rModelPart.Nodes()).push_back(*(i_node.base()));	
		rModelPart.Nodes().back().SetId(nodeId);
		nodeId+=1;
	      }
	  }
	  else{

	    for(ModelPart::ElementsContainerType::iterator i_elem = rModelPart.ElementsBegin(MeshId) ; i_elem != rModelPart.ElementsEnd(MeshId) ; i_elem++)
	      {
		PointsArrayType& vertices=i_elem->GetGeometry().Points();
		for(unsigned int i=0; i<vertices.size(); i++)
		  {
		    vertices[i].Set(BLOCKED);
		  }
		
		(rModelPart.Elements()).push_back(*(i_elem.base()));	
		rModelPart.Elements().back().SetId(elemId);
		elemId+=1;
	      }
	    
	    
	    for(ModelPart::NodesContainerType::iterator i_node = temporal_nodes.begin() ; i_node != temporal_nodes.end() ; i_node++)
	      {
		//i_node->PrintInfo(std::cout);
		//std::cout<<std::endl;

		if(i_node->Is(BLOCKED) || i_node->Is(RIGID))
		  {
		    i_node->Reset(ISOLATED);   //reset isolated
		    i_node->Reset(NEW_ENTITY); //reset if was new 
		    i_node->Reset(TO_REFINE);  //reset if was labeled to refine (to not duplicate boundary conditions)
		    i_node->Reset(BLOCKED); 
		    
		    (rModelPart.Nodes(MeshId)).push_back(*(i_node.base()));
		    (rModelPart.Nodes()).push_back(*(i_node.base()));	
		    rModelPart.Nodes().back().SetId(nodeId);
		    nodeId+=1;

		  }
		else{

		  i_node->Set(ISOLATED);
		  i_node->Reset(NEW_ENTITY); //reset if was new 
		  i_node->Reset(TO_REFINE);  //reset if was labeled to refine (to not duplicate boundary conditions)
		  i_node->Reset(BLOCKED); 

		  if( mOptions.Is(ModelerUtilities::KEEP_ISOLATED_NODES) ){

		    (rModelPart.Nodes(MeshId)).push_back(*(i_node.base()));
		    (rModelPart.Nodes()).push_back(*(i_node.base()));	
		    rModelPart.Nodes().back().SetId(nodeId);
		    nodeId+=1;

		  }
		    
		}
		
		if(i_node->Is(BOUNDARY)){

		  if(i_node->IsNot(RIGID))
		    i_node->Set(FREE_SURFACE);

		}
		else{

		  if( i_node->SolutionStepsDataHas(CONTACT_FORCE) )
		    noalias(i_node->GetSolutionStepValue(CONTACT_FORCE)) = ZeroNormal;
		}
				
	      }
	  }
	  //rModelPart.Nodes(MeshId).Sort();  

	  for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.ConditionsBegin(MeshId) ; i_cond != rModelPart.ConditionsEnd(MeshId) ; i_cond++)
	    {
	      i_cond->Reset(NEW_ENTITY); //reset here if the node is inserted
	      KeepConditions.push_back(*(i_cond.base()));
	      KeepConditions.back().SetId(condId);
	      condId+=1;	
	    }
	}


      for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.ConditionsBegin(); i_cond!= rModelPart.ConditionsEnd(); i_cond++)
	{
	  if(i_cond->Is(CONTACT)){
	    KeepConditions.push_back(*(i_cond.base()));
	    KeepConditions.back().SetId(condId);
	    condId+=1;
	  }
	}
      
      rModelPart.Conditions().swap(KeepConditions);
      
      
      //Sort
      // rModelPart.Nodes().Sort();
      // rModelPart.Elements().Sort();
      // rModelPart.Conditions().Sort();     
      
      //Unique
      // rModelPart.Nodes().Unique();
      // rModelPart.Elements().Unique();
      // rModelPart.Conditions().Unique();
      
      //Sort Again to have coherent numeration for nodes (mesh with shared nodes)
      unsigned int consecutive_index = 1;
      for(ModelPart::NodesContainerType::iterator in = rModelPart.NodesBegin(0) ; in != rModelPart.NodesEnd(0) ; in++)
	in->SetId(consecutive_index++);
      
      if( EchoLevel > 0 )
	std::cout<<"   [ END MESH [Id=0] [Elems=:"<<rModelPart.NumberOfElements()<<"|Nodes="<<rModelPart.NumberOfNodes()<<"|Conds="<<rModelPart.NumberOfConditions()<<"] ] "<<std::endl;      
      

      KRATOS_CATCH(" ")

   }
    

  //*******************************************************************************************
  //*******************************************************************************************

  void CleanMeshFlags(ModelPart& rModelPart, ModelPart::IndexType MeshId)
  {
    
    KRATOS_TRY

    for(ModelPart::NodesContainerType::const_iterator i_node = rModelPart.NodesBegin(MeshId); i_node != rModelPart.NodesEnd(MeshId); i_node++)
      {

	i_node->Reset(NEW_ENTITY); //reset here if the node is labeled as insert 
	i_node->Reset(TO_REFINE);  //reset here if the node is labeled as refine (to not duplicate bo

      }

    for(ModelPart::ConditionsContainerType::iterator i_cond = rModelPart.ConditionsBegin(MeshId) ; i_cond != rModelPart.ConditionsEnd(MeshId) ; i_cond++)
      {
	i_cond->Reset(NEW_ENTITY); //reset here if the node is inserted
      }

    KRATOS_CATCH(" ")
  }


    //*******************************************************************************************
    //*******************************************************************************************


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
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
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
    ModelStartEndMeshingProcess& operator=(ModelStartEndMeshingProcess const& rOther);

    /// Copy constructor.
    //ModelStartEndMeshingProcess(ModelStartEndMeshingProcess const& rOther);


    ///@}

  }; // Class ModelStartEndMeshingProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    ModelStartEndMeshingProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const ModelStartEndMeshingProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // KRATOS_MODEL_START_END_MESHING_PROCESS_H_INCLUDED  defined 
