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
    ContactModelStartEndMeshingProcess(ModelPart& rMainModelPart,
				       Flags Options,
				       int EchoLevel = 0)
      : ModelStartEndMeshingProcess(rMainModelPart, Options, EchoLevel)
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

      this->ClearContactConditions();

      this->ClearContactFlags();

      //CONDITIONS MASTER_ELEMENTS and MASTER_NODES SEARCH
      if( mrMainModelPart.GetProcessInfo()[IS_RESTARTED] == true ){
      	BuildModelPartBoundaryProcess BuildBoundaryProcess(mrMainModelPart, mrMainModelPart.Name(), mEchoLevel);
      	BuildBoundaryProcess.SearchConditionMasters();
      }
      
      //Update Boundary Normals before Contact Search   
      //(needed when meshing of the domains is not performed:
      // normal directions change with mesh movement)
      BoundaryNormalsCalculationUtilities BoundaryComputation;
      BoundaryComputation.CalculateWeightedBoundaryNormals(mrMainModelPart, mEchoLevel);

      //mrMainModelPart.Conditions().Sort();
      //mrMainModelPart.Conditions().Unique();

      
      KRATOS_CATCH(" ")
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    virtual void ExecuteFinalize()
    {
      KRATOS_TRY
      
      // Add contact conditions
      this->AddContactConditions();

      // Restore contact nodal flags
      this->RestoreContactFlags();
      
      // Renumerate conditions
      this->RenumerateConditions();

      // Clear Contact Forces
      this->ClearContactForces();

      // Clear Contact Normals
      this->ClearContactNormals();

      
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
      //mrMainModelPart.Elements().Sort();
      int LastElementId   = mrMainModelPart.Elements().back().Id();
      int LastConditionId = mrMainModelPart.Conditions().back().Id();
      if(LastElementId>LastConditionId){
	consecutive_index = LastElementId+1;
      }
      else{
	consecutive_index = LastConditionId+1;
      }
	
      for(ModelPart::ConditionsContainerType::iterator it = mrMainModelPart.ConditionsBegin(); it!=mrMainModelPart.ConditionsEnd(); it++){
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



    //**************************************************************************
    //**************************************************************************

    void ClearContactForces()
    {
      KRATOS_TRY

     
      ModelPart::NodesContainerType& rNodes = mrMainModelPart.Nodes();

      // create contact condition for rigid and deformable bodies
      for(ModelPart::NodesContainerType::ptr_iterator nd = rNodes.ptr_begin(); nd != rNodes.ptr_end(); ++nd)
	{
	  if( (*nd)->Is(BOUNDARY) && (*nd)->IsNot(CONTACT) ){
	    array_1d<double, 3> & ContactForce  = (*nd)->FastGetSolutionStepValue(CONTACT_FORCE);
	    ContactForce.clear();
	  }

	}
                
      KRATOS_CATCH( "" )
	
    }
    
    //**************************************************************************
    //**************************************************************************

    void ClearContactNormals()
    {
      KRATOS_TRY

     
      ModelPart::NodesContainerType& rNodes = mrMainModelPart.Nodes();

      // create contact condition for rigid and deformable bodies
      for(ModelPart::NodesContainerType::ptr_iterator nd = rNodes.ptr_begin(); nd != rNodes.ptr_end(); ++nd)
	{
	  if( (*nd)->Is(BOUNDARY) && (*nd)->IsNot(CONTACT) ){
	    array_1d<double, 3> & ContactNormal  = (*nd)->FastGetSolutionStepValue(CONTACT_NORMAL);
	    ContactNormal.clear();
	  }

	}
                
      KRATOS_CATCH( "" )
	
    }
    
    //**************************************************************************
    //**************************************************************************

    void ClearContactFlags ( )
    {
     KRATOS_TRY

      for(ModelPart::NodesContainerType::iterator i_node = mrMainModelPart.NodesBegin(); i_node!= mrMainModelPart.NodesEnd(); i_node++)
	{
	  if( i_node->Is(CONTACT) ){
	    i_node->Set(CONTACT,false);
	  }
	}

      KRATOS_CATCH( "" )
    }


    //**************************************************************************
    //**************************************************************************

    void RestoreContactFlags ( )
    {

      KRATOS_TRY

      for(ModelPart::ConditionsContainerType::iterator i_cond = mrMainModelPart.ConditionsBegin(); i_cond!= mrMainModelPart.ConditionsEnd(); i_cond++)
	{
	  if( i_cond->Is(CONTACT) ){
	    for(unsigned int i=0; i<i_cond->GetGeometry().size(); i++)
	      {
		i_cond->GetGeometry()[i].Set(CONTACT,true);
	      }
	  }
	}

      KRATOS_CATCH( "" )
    }
    
    //**************************************************************************
    //**************************************************************************

    void AddContactConditions()
    {

      KRATOS_TRY

      //Add contact conditions from contact domain
      std::string ModelPartName;
      for(ModelPart::SubModelPartIterator i_mp= mrMainModelPart.SubModelPartsBegin(); i_mp!=mrMainModelPart.SubModelPartsEnd(); i_mp++)
	{
	  if(i_mp->Is(CONTACT))
	    ModelPartName = i_mp->Name();
	}
      
      ModelPart& ContactModelPart = mrMainModelPart.GetSubModelPart(ModelPartName);

      //Add contact conditions to computing domain
      for(ModelPart::SubModelPartIterator i_mp= mrMainModelPart.SubModelPartsBegin(); i_mp!=mrMainModelPart.SubModelPartsEnd(); i_mp++)
	{
	  if(i_mp->Is(SOLID) && i_mp->Is(ACTIVE))
	    ModelPartName = i_mp->Name();
	}

      
      AddContactConditions(ContactModelPart, mrMainModelPart.GetSubModelPart(ModelPartName));

      //Add contact conditions to  main domain (if added in computing domaing with AddCondition, are automatically added to the main domain, else use push_back)
      AddContactConditions(ContactModelPart, mrMainModelPart);

      KRATOS_CATCH( "" )
	
    }

    //**************************************************************************
    //**************************************************************************

    void AddContactConditions(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart)
    {

      KRATOS_TRY

      //*******************************************************************
      //adding contact conditions
      //
	
      if( mEchoLevel >= 1 ){
	std::cout<<" ["<<rDestinationModelPart.Name()<<" :: CONDITIONS [OLD:"<<rDestinationModelPart.NumberOfConditions();
      }

      for(ModelPart::ConditionsContainerType::iterator ic = rOriginModelPart.ConditionsBegin(); ic!= rOriginModelPart.ConditionsEnd(); ic++)
	{

	  if(ic->Is(CONTACT))
	    rDestinationModelPart.AddCondition(*(ic.base()));
	  
	}
      
      if( mEchoLevel >= 1 ){
	std::cout<<" / NEW:"<<rDestinationModelPart.NumberOfConditions()<<"] "<<std::endl;
      }
            
      KRATOS_CATCH( "" )
	
    }




    //**************************************************************************
    //**************************************************************************

    void ClearContactConditions()
    {

      KRATOS_TRY

      //Clear contact conditions from computing domain
      std::string ModelPartName;
      for(ModelPart::SubModelPartIterator i_mp= mrMainModelPart.SubModelPartsBegin(); i_mp!=mrMainModelPart.SubModelPartsEnd(); i_mp++)
	{
	  if(i_mp->Is(SOLID) && i_mp->Is(ACTIVE))
	    ModelPartName = i_mp->Name();
	}
      
      ClearContactConditions(mrMainModelPart.GetSubModelPart(ModelPartName));

      //Clear contact conditions from the main domain
      ClearContactConditions(mrMainModelPart);

      KRATOS_CATCH( "" )
	
    }
    

    //**************************************************************************
    //**************************************************************************

    void ClearContactConditions(ModelPart& rModelPart)
    {

      KRATOS_TRY

      //*******************************************************************
      //clearing contact conditions
      //
	
      if( mEchoLevel >= 1 ){
	std::cout<<" ["<<rModelPart.Name()<<" :: CONDITIONS [OLD:"<<rModelPart.NumberOfConditions();
      }

      ModelPart::ConditionsContainerType PreservedConditions;

      for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(); ic!= rModelPart.ConditionsEnd(); ic++)
	{

	  if(ic->IsNot(CONTACT) && ic->GetGeometry().size() > 1){
	    PreservedConditions.push_back(*(ic.base()));
	  }
	}
      
      rModelPart.Conditions().swap(PreservedConditions);
	      
      //rModelPart.Conditions().Sort();
      //rModelPart.Conditions().Unique();

      if( mEchoLevel >= 1 ){
	std::cout<<" / NEW:"<<rModelPart.NumberOfConditions()<<"] "<<std::endl;
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

#endif // KRATOS_CONTACT_MODEL_START_END_MESHING_PROCESS_H_INCLUDED  defined 
