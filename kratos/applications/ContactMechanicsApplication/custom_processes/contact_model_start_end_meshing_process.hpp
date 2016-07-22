//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_CONTACT_MODEL_START_END_MESHING_PROCESS_H_INCLUDED )
#define  KRATOS_CONTACT_MODEL_START_END_MESHING_PROCESS_H_INCLUDED


// System includes

// External includes

// Project includes

#include "includes/model_part.h"
#include "custom_utilities/modeler_utilities.hpp"
#include "custom_processes/model_start_end_meshing_process.hpp"

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
  class KRATOS_API(CONTACT_MECHANICS_APPLICATION) ContactModelStartEndMeshingProcess
    : public ModelStartEndMeshingProcess
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ContactModelStartEndMeshingProcess
    KRATOS_CLASS_POINTER_DEFINITION( ContactModelStartEndMeshingProcess );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ContactModelStartEndMeshingProcess(ModelPart& rModelPart,
				       Flags Options,
				       int EchoLevel = 0)
      : ModelStartEndMeshingProcess(rModelPart, Options, EchoLevel)
    { 
    }

    /// Destructor.
    virtual ~ContactModelStartEndMeshingProcess()
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

      // to keep history, maybe must be cleaned after remeshing
      this->ClearContactConditions();

      //Sort Conditions
      // this->SortModelPartConditions(); //sorted in ClearContactConditions()
     
      //Update Boundary Normals before Contact Search   
      //(needed when meshing of the domains is not performed:
      // normal directions change with mesh movement)
      BoundaryNormalsCalculationUtilities BoundaryComputation;
      BoundaryComputation.CalculateBoundaryNormals(mrModelPart, mEchoLevel);

      
      KRATOS_CATCH(" ")
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    virtual void ExecuteFinalize()
    {
      KRATOS_TRY
      
      // Renumerate conditions
      this->RenumerateConditions();

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
      return "ContactModelStartEndMeshingProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
      rOStream << "ContactModelStartEndMeshingProcess";
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

    void RenumerateConditions()
    {

      KRATOS_TRY

	unsigned int consecutive_index = 0;
	//Renumerate conditions to add in the end of the Elements array (for writing purposes in the ID)
	//mrModelPart.Elements().Sort();
	int LastElementId   = (mrModelPart.Elements().end()-1)->Id();
        int LastConditionId = (mrModelPart.Conditions().end()-1)->Id();
	if(LastElementId>LastConditionId){
	  consecutive_index = LastElementId+1;
	}
	else{
	  consecutive_index = LastConditionId+1;
	}
	
	for(ModelPart::ConditionsContainerType::iterator it = mrModelPart.ConditionsBegin(); it!=mrModelPart.ConditionsEnd(); it++){
	  if(it->Is(CONTACT)){
	    it->SetId(consecutive_index); 
	    consecutive_index++;
	  }
	}

      KRATOS_CATCH(" ")

    }

   
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

    void ClearContactConditions()
    {

      KRATOS_TRY

      //*******************************************************************
      //clearing contact conditions
      //
	
      std::cout<<" Total Conditions BEFORE: ["<<mrModelPart.Conditions().size()<<"]"<<std::endl;

      ModelPart::ConditionsContainerType PreservedConditions;

      int contact_Id=0;
      for(ModelPart::ConditionsContainerType::iterator ic = mrModelPart.ConditionsBegin(); ic!= mrModelPart.ConditionsEnd(); ic++)
	{

	  if(ic->IsNot(CONTACT)){
	    contact_Id+=1;
	    PreservedConditions.push_back(*(ic.base()));
	    PreservedConditions.back().SetId(contact_Id);
	  }
	}
      
      mrModelPart.Conditions().swap(PreservedConditions);
	      
      mrModelPart.Conditions().Sort();
      mrModelPart.Conditions().Unique();

      std::cout<<" Total Conditions CLEAN: ["<<mrModelPart.Conditions().size()<<"]"<<std::endl;
      
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
    ContactModelStartEndMeshingProcess& operator=(ContactModelStartEndMeshingProcess const& rOther);

    /// Copy constructor.
    //ContactModelStartEndMeshingProcess(ContactModelStartEndMeshingProcess const& rOther);


    ///@}

  }; // Class ContactModelStartEndMeshingProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    ContactModelStartEndMeshingProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const ContactModelStartEndMeshingProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // CONTACT_KRATOS_MODEL_START_END_MESHING_PROCESS_H_INCLUDED  defined 
