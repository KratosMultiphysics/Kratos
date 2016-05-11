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
#include <string>
#include <iostream>


// External includes
#include <boost/timer.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "custom_utilities/modeler_utilities.hpp"

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
      unsigned int consecutive_index = 1;
      for(ModelPart::ConditionsContainerType::iterator ic = mrModelPart.ConditionsBegin(); ic!=mrModelPart.ConditionsEnd(); ic++)
	(ic)->SetId(consecutive_index++);

      mrModelPart.Conditions().Sort();
      mrModelPart.Conditions().Unique();

      // //Sort Elements
      // consecutive_index = 1;
      // for(ModelPart::ElementsContainerType::iterator ie = mrModelPart.ElementsBegin(); ie!=mrModelPart.ElementsEnd(); ie++)
      //   ie->SetId(consecutive_index++);

      // mrModelPart.Elements().Sort();
      // mrModelPart.Elements().Unique();

      // //Sort Nodes, set STRUCTURE NODES (RIGID TOOL NODES) AT END
      // consecutive_index = 1;
      // unsigned int reverse_index = mrModelPart.Nodes().size();
      // for(ModelPart::NodesContainerType::iterator in = mrModelPart.NodesBegin(); in!=mrModelPart.NodesEnd(); in++){
      //   if(in->IsNot(STRUCTURE) ){
      //     in->SetId(consecutive_index++);
      //   }
      //   else{
      //     in->SetId(reverse_index--);
      //   }
      // }

      // mrModelPart.Nodes().Sort();
      // mrModelPart.Nodes().Unique();

      //check that contact forces are not cleared
      // for(ModelPart::NodesContainerType::iterator in = mrModelPart.NodesBegin(); in!=mrModelPart.NodesEnd(); in++){

      // 	if(in->Is(BOUNDARY) ){
      // 	  array_1d<double, 3 > &ContactForce = in->FastGetSolutionStepValue(CONTACT_FORCE);
      // 	  std::cout<<" ContactForce["<<in->Id()<<"] "<<ContactForce<<std::endl;
      // 	}
	
      // }


      KRATOS_CATCH(" ")
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    virtual void ExecuteFinalize()
    {
      KRATOS_TRY

      //Once all meshes are build, the main mesh Id=0 must be reassigned
      ModelerUtilities ModelerUtils;

      MeshDataTransferUtilities TransferUtilities;

      unsigned int NumberOfMeshes = mrModelPart.NumberOfMeshes();

      if(NumberOfMeshes>1){
	ModelerUtils.BuildTotalMesh(mrModelPart, mEchoLevel);
      }
      else{
	ModelerUtils.CleanMeshFlags(mrModelPart,0);
      }
	
      unsigned int start = 0;
      if(NumberOfMeshes>1) 
	start=1;
      
      bool remesh_performed =true; //false; //deactivate searches
	   
      if(remesh_performed){


	//CHANGE ELEMENTS WITH ALL NODES IN BOUNDARY  (TIP ELEMENTS)
	for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
	  {
		
	    //if(Meshes[MeshId].Is( ModelerUtilities::REMESH )){
	      // ChangeTipElementsUtilities TipElements;
	      // if(mpMeshingVariables->AvoidTipElementsFlag){
	      //   //TipElements.SwapDiagonals(mrModelPart,MeshId);
	      // }
	    //}
	  }		
	    
	    
	//NEIGHBOURS SEARCH
	NodalNeighboursSearchProcess FindNeighbours(mrModelPart);
	FindNeighbours.Execute();

	//NODAL_H SEARCH	    
	//FindNodalHProcess FindNodalH(mrModelPart);
	//FindNodalH.Execute();

	//CONDITIONS MASTER_ELEMENTS and MASTER_NODES SEARCH
	BuildMeshBoundaryProcess BuildBoundaryProcess(mrModelPart,2);
	for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
	  {
	    BuildBoundaryProcess.SearchConditionMasters(MeshId);
	  }


	//BOUNDARY NORMALS SEARCH // SHRINKAGE FACTOR
	//ComputeBoundaryNormals BoundUtils;
	BoundaryNormalsCalculationUtilities BoundaryComputation;
	for(unsigned int MeshId=start; MeshId<NumberOfMeshes; MeshId++)
	  {
	    BoundaryComputation.CalculateMeshBoundaryNormals(mrModelPart, MeshId, mEchoLevel);
	  }



      }

      // if(!out_buffer_active)
      //   std::cout.rdbuf(buffer);

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


    ///@}
    ///@name Protected Operators
    ///@{


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
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    ModelPart& mrModelPart;

    Flags mOptions;

    int mEchoLevel;

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
