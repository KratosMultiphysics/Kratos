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
#include "custom_processes/build_mesh_boundary_process.hpp"

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
  typedef  ModelPart::NodesContainerType NodesContainerType;
  typedef  ModelPart::ElementsContainerType ElementsContainerType;
  typedef  ModelPart::ConditionsContainerType ConditionsContainerType;
  typedef  ModelPart::MeshType::GeometryType::PointsArrayType PointsArrayType;
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
    ModelStartEndMeshingProcess(ModelPart& rModelPart,
				Flags Options,
				int EchoLevel = 0)
      : mrModelPart(rModelPart)
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
      
      // Restore mehses coherency
      this->BuildModelPartMeshes();

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

    ModelPart& mrModelPart;

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
      for(ModelPart::ConditionsContainerType::iterator ic = mrModelPart.ConditionsBegin(); ic!=mrModelPart.ConditionsEnd(); ic++)
	(ic)->SetId(consecutive_index++);

      mrModelPart.Conditions().Sort();
      mrModelPart.Conditions().Unique();


      KRATOS_CATCH(" ")

    }


    //*******************************************************************************************
    //*******************************************************************************************

    void SortModelPartElements()
    {

      KRATOS_TRY


      // Sort Elements
      unsigned int consecutive_index = 1;
      for(ModelPart::ElementsContainerType::iterator ie = mrModelPart.ElementsBegin(); ie!=mrModelPart.ElementsEnd(); ie++)
	ie->SetId(consecutive_index++);

      mrModelPart.Elements().Sort();
      mrModelPart.Elements().Unique();

      KRATOS_CATCH(" ")

     }


    //*******************************************************************************************
    //*******************************************************************************************
    
    void SortModelPartNodes()
    {

      KRATOS_TRY


      //Sort Nodes, set STRUCTURE nodes at end
      unsigned int consecutive_index = 1;
      unsigned int reverse_index = mrModelPart.Nodes().size();
      for(ModelPart::NodesContainerType::iterator in = mrModelPart.NodesBegin(); in!=mrModelPart.NodesEnd(); in++)
	{
	  if(in->IsNot(STRUCTURE) ){
	    in->SetId(consecutive_index++);
	  }
	  else{
	    in->SetId(reverse_index--);
	  }
	}

      mrModelPart.Nodes().Sort();
      mrModelPart.Nodes().Unique();

      KRATOS_CATCH(" ")

     }

    
    //*******************************************************************************************
    //*******************************************************************************************
    
    void BuildModelPartMeshes()
    {

      KRATOS_TRY

      //Once all meshes are build, the main mesh Id=0 must be reassigned
      ModelerUtilities ModelerUtils;

      unsigned int NumberOfMeshes = mrModelPart.NumberOfMeshes();

      if(NumberOfMeshes>1){
	this->BuildTotalMesh(mrModelPart, mEchoLevel);
      }
      else{
	this->CleanMeshFlags(mrModelPart,0);
      }


      KRATOS_CATCH(" ")

     }


    
    //*******************************************************************************************
    //*******************************************************************************************
    
    virtual void PerformModelSearches()
    {

      KRATOS_TRY
       
	//NODAL NEIGHBOURS SEARCH
	NodalNeighboursSearchProcess FindNeighbours(mrModelPart);
	FindNeighbours.Execute();

	//NODAL_H SEARCH	    
	//FindNodalHProcess FindNodalH(mrModelPart);
	//FindNodalH.Execute();

	//CONDITIONS MASTER_ELEMENTS and MASTER_NODES SEARCH
	BuildMeshBoundaryProcess BuildBoundaryProcess(mrModelPart);
	BuildBoundaryProcess.SearchConditionMasters();

	//BOUNDARY NORMALS SEARCH and  SHRINKAGE FACTOR
	BoundaryNormalsCalculationUtilities BoundaryComputation;
	BoundaryComputation.CalculateMeshBoundaryNormals(mrModelPart, mEchoLevel);   


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
	// 	TipElements.SwapDiagonals(mrModelPart,MeshId);
	//       }
	//     }
	//   }	
	
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
		  //std::cout<<" NOT ENGAGED NODE "<<i_node->Id()<<std::endl;
		  i_node->Set(ISOLATED);		 	  
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
      rModelPart.Nodes().Sort();
      rModelPart.Elements().Sort();
      rModelPart.Conditions().Sort();     
      
      //Unique
      rModelPart.Nodes().Unique();
      rModelPart.Elements().Unique();
      rModelPart.Conditions().Unique();
      
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

  void CleanMeshFlags(ModelPart& rModelPart,ModelPart::IndexType MeshId)
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
