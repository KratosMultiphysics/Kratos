//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_COMPUTING_MODEL_PART_TRANSFER_PROCESS_H_INCLUDED)
#define  KRATOS_COMPUTING_MODEL_PART_TRANSFER_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Process to transfer model part entities to a computing submodelpart

class ComputingModelPartTransferProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputingModelPartTransferProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputingModelPartTransferProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    ComputingModelPartTransferProcess(ModelPart& model_part,
                                      Parameters rParameters
                                      ) : Process() , mrModelPart(model_part)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
            {
               "reference_element_type": "Element2D3N",
               "reference_contact_condition_type": "CompositeCondition2D2N"
            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mpReferenceElement = &(KratosComponents<Element>::Get(rParameters["reference_element_type"].GetString()));
        mpReferenceContactCondition = &(KratosComponents<Condition>::Get(rParameters["reference_contact_condition_type"].GetString()));

        KRATOS_CATCH("")
    }


    /// Destructor.
    virtual ~ComputingModelPartTransferProcess() {}


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


    /// Execute method is used to execute the ComputingModelPartTransferProcess algorithms.
    void Execute() override
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
      KRATOS_TRY
          
      this->CreateAndTransferEntities(mrModelPart->GetParentModelPart(), mrModelPart, *mpReferenceElement, *mpReferenceContactConditions);
 
      KRATOS_CATCH("")
    }

    /// this function is designed for being execute once before the solution loop but after all of the
    /// solvers where built
    void ExecuteBeforeSolutionLoop() override
    {
    }


    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
      KRATOS_TRY

      this->CreateAndTransferEntities(mrModelPart->GetParentModelPart(), mrModelPart, *mpReferenceElement, *mpReferenceContactConditions);

      KRATOS_CATCH("")
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {
    }


    /// this function will be executed at every time step BEFORE  writing the output
    void ExecuteBeforeOutputStep() override
    {
    }


    /// this function will be executed at every time step AFTER writing the output
    void ExecuteAfterOutputStep() override
    {
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    void ExecuteFinalize() override
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
    std::string Info() const override
    {
        return "ComputingModelPartTransferProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputingModelPartTransferProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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

    /// Copy constructor.
    ComputingModelPartTransferProcess(ComputingModelPartTransferProcess const& rOther);

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
  
    const Element* mpReferenceElement;

    const Condition* mpReferenceContactCondition;

    ///@}
    ///@name Private Operators
    ///@{
    ///@}
    ///@name Private Operations
    ///@{

    //**********************************************************************************************
    //**********************************************************************************************
    ///This function fills the @param DestinationModelPart using the data obtained from @param  OriginModelPart
    ///the elements and conditions of the DestinationModelPart part use the same connectivity (and id) as the
    ///OriginModelPart but their type is determined by @param rReferenceElement and @param rReferenceBoundaryCondition

    void CreateAndTransferEntities(ModelPart& OriginModelPart,
                                   ModelPart& DestinationModelPart,
                                   Element const& rReferenceElement,
                                   Condition const& rReferenceContactCondition)
    {
        KRATOS_TRY

        //nodes and conditions previously assigned:
        if(DestinationModelPart.Nodes().size() == 0)
          KRATOS_ERROR << " Destination model part with no nodes " << std::endl;
        
        //clear presious conditions
        DestinationModelPart.Conditions().clear();
        
        //clear previous elements
        DestinationModelPart.Elements().clear();

        //assigning Properties
        DestinationModelPart.SetProperties( OriginModelPart.pProperties() );

        //generating the elements
        for (ModelPart::ElementsContainerType::iterator ie = OriginModelPart.ElementsBegin(); ie != OriginModelPart.ElementsEnd(); ie++)
        {
	  if( ie->GetGeometry().size()>2 && ie->IsNot(RIGID) ){
	  
            Properties::Pointer properties = ie->pGetProperties();
            Element::Pointer p_element = rReferenceElement.Create(ie->Id(), ie->GetGeometry(), properties);
	        
	    p_element->Set(THERMAL);
	    
	    //set mechanical element as pointer
	    WeakPointerVector< Element > MasterElements;
	    MasterElements.push_back(Element::WeakPointer( *(ie.base()) ) );

	    p_element->SetValue(MASTER_ELEMENTS, MasterElements);

            DestinationModelPart.Elements().push_back(p_element);
	    
	  }
        }

        //generating the conditions
        for (ModelPart::ConditionsContainerType::iterator ic = OriginModelPart.ConditionsBegin(); ic != OriginModelPart.ConditionsEnd(); ic++)
        {
          if(ic->Is(CONTACT))
          {

            if( ic->GetGeometry().size() >= 2 ){
		
              Properties::Pointer properties = ic->pGetProperties();
		  
              Condition::Pointer p_condition = rReferenceBoundaryCondition.Create(ic->Id(), ic->GetGeometry(), properties);
		  
              //set mechanical variables to contact conditions:
              p_condition->AssignFlags( *ic ); //it overwrites THERMAL and CONTACT
              p_condition->Data() = ic->Data();
		  
              p_condition->Set(ACTIVE);
              p_condition->Set(THERMAL);
		  	  
              DestinationModelPart.Conditions().push_back(p_condition);
		  
              //set only thermal conditions from skin boundary conditions:
            }
		
          }

        }
        
        //generating tables
	DestinationModelPart.Tables() = OriginModelPart.Tables();

        Communicator::Pointer pComm = OriginModelPart.GetCommunicator().Create();
        DestinationModelPart.SetCommunicator(pComm);

	if( GetEchoLevel() >= 1 )
	  std::cout<<" DestinationModelPart "<<DestinationModelPart<<std::endl;

        KRATOS_CATCH("")
    }
  
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    ComputingModelPartTransferProcess& operator=(ComputingModelPartTransferProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class ComputingModelPartTransferProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ComputingModelPartTransferProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ComputingModelPartTransferProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_COMPUTING_MODEL_PART_TRANSFER_PROCESS_H_INCLUDED  defined
