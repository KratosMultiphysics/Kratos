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
#include "processes/transfer_entitites_between_model_parts_process.hpp"

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
               "model_part_name": "computing_sub_model_part",
               "element_transfer_flags": [],
               "composite_boundary_conditions": True,
               "create_new_entities": False,
               "creation_options" :{
                     "entity_creation_flags" [],
                     "reference_element_type": "ElementType",
                     "reference_contact_condition_type": "ContactConditionType"
               }
            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        for(unsigned int i=0; i<rParameters["element_transfer_flags"].size(); ++i)
        {
          mTransferFlags.push_back(rParameters["element_transfer_flags"][i].GetString());
        }

        mCompositeConditions = rParameters["composite_boundary_conditions"].GetBool();
        mCreateNewEntities = rParameters["create_new_entitites"].GetBool();
        if( mCreateNewEntities ){
          for(unsigned int i=0; i<rParameters["creation_options"]["entity_creation_flags"].size(); ++i)
          {
            mCreationFlags.push_back(rParameters["creation_options"]["entity_creation_flags"][i].GetString());
          }
          mpReferenceElement = &(KratosComponents<Element>::Get(rParameters["creation_options"]["reference_element_type"].GetString()));
          mpReferenceContactCondition = &(KratosComponents<Condition>::Get(rParameters["creation_options"]["reference_contact_condition_type"].GetString()));
        }

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
      KRATOS_TRY

      if( mCreateNewEntities ){
        this->CreateAndTransferEntities(mrModelPart->GetParentModelPart(), mrModelPart, mTransferFlags, mCreationFlags, *mpReferenceElement, *mpReferenceContactConditions);
      }
      else{
        this->TransferEntities(mrModelPart->GetParentModelPart(), mrModelPart, mTransferFlags);
      }

      KRATOS_CATCH("")
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

    bool mCreateNewEntities;

    bool mCompositeConditions;
  
    const std::vector<Flags> mTransferFlags;
  
    const std::vector<Flags> mCreationFlags;
  
    ///@}
    ///@name Private Operators
    ///@{
    ///@}
    ///@name Private Operations
    ///@{

    //**********************************************************************************************
    //**********************************************************************************************
    ///This function fills the @param DestinationModelPart using the data obtained from @param rOriginModelPart
    ///the elements and conditions of the rDestinationModelPart part use the same connectivity (and id) as the
    ///OriginModelPart but their type is determined by @param rReferenceElement and @param rReferenceBoundaryCondition

    void CreateAndTransferEntities(ModelPart& rOriginModelPart,
                                   ModelPart& rDestinationModelPart,
                                   const std::vector<Flags>& rTransferFlags,
                                   const std::vector<Flags>& rCreationFlags,
                                   Element const& rReferenceElement,
                                   Condition const& rReferenceContactCondition,
                                   bool& rCompositeConditions)
    {
        KRATOS_TRY

        //clear previous nodes
        std::string Nodes = "Nodes";
        rDestinationModelPart.Nodes().clear();
        TransferEntitiesBetweenModelPartsProcess NodesTransferProcess = TransferEntitiesBetweenModelPartsProcess(rDestinationModelPart, rOriginModelPart, Nodes);
        NodesTransferProcess.Execute();

        //clear previous elements
        rDestinationModelPart.Elements().clear();

        //generating the elements
        for (ModelPart::ElementsContainerType::iterator ie = rOriginModelPart.ElementsBegin(); ie != rOriginModelPart.ElementsEnd(); ++ie)
        {
          if (this->MatchFlags(*(ie.base()),rCreationFlags))
          {
            Properties::Pointer pProperties = ie->pGetProperties();
            Element::Pointer pElement = rReferenceElement.Create(ie->Id(), ie->GetGeometry(), pProperties);

	    //set origin element as pointer
	    WeakPointerVector< Element > MasterElements;
	    MasterElements.push_back(Element::WeakPointer( *(ie.base()) ) );

	    pElement->SetValue(MASTER_ELEMENTS, MasterElements);

            rDestinationModelPart.Elements().push_back(pElement);
	  }
        }


        //clear previous conditions
        std::string Conditions = "Conditions";
        rDestinationModelPart.Conditions().clear();

        if( rCompositeConditions ){
          this->TransferBoundaryConditions(rOriginModelPart, rDestinationModelPart);
        }
        else{
          TransferEntitiesBetweenModelPartsProcess ConditionsTransferProcess = TransferEntitiesBetweenModelPartsProcess(rDestinationModelPart, rOriginModelPart, Conditions, rTransferFlags);
          ConditionsTransferProcess.Execute();
        }
        
        //generating contact conditions
        for (ModelPart::ConditionsContainerType::iterator ic = rOriginModelPart.ConditionsBegin(); ic != rOriginModelPart.ConditionsEnd(); ++ic)
        {
          if(ic->Is(CONTACT))
          {
            Properties::Pointer pProperties = ic->pGetProperties();
            Condition::Pointer pCondition = rReferenceContactCondition.Create(ic->Id(), ic->GetGeometry(), pProperties);

            //set mechanical variables to contact conditions:
            pCondition->Data() = ic->Data();

            rDestinationModelPart.Conditions().push_back(pCondition);
          }
        }
               
        //assigning Properties
        rDestinationModelPart.SetProperties( rOriginModelPart.pProperties() );

        //assigning Tables
	rDestinationModelPart.Tables() = rOriginModelPart.Tables();
 
        Communicator::Pointer pComm = OriginModelPart.GetCommunicator().Create();
        rDestinationModelPart.SetCommunicator(pComm);

	if( GetEchoLevel() >= 1 )
	  std::cout<<" DestinationModelPart "<<rDestinationModelPart<<std::endl;

        KRATOS_CATCH("")
    }


    void TransferEntities(ModelPart& rOriginModelPart,
                          ModelPart& rDestinationModelPart,
                          const std::vector<Flags>& rTransferFlags,
                          bool& rCompositeConditions)
    {
        KRATOS_TRY

        //clear previous nodes
        std::string Nodes = "Nodes";
        rDestinationModelPart.Nodes().clear();
        TransferEntitiesBetweenModelPartsProcess NodesTransferProcess = TransferEntitiesBetweenModelPartsProcess(rDestinationModelPart, rOriginModelPart, Nodes);
        NodesTransferProcess.Execute();

        //clear previous elements
        std::string Elements = "Elements";
        rDestinationModelPart.Elements().clear();
        TransferEntitiesBetweenModelPartsProcess ElementsTransferProcess = TransferEntitiesBetweenModelPartsProcess(rDestinationModelPart, rOriginModelPart, Elements, rTransferFlags);
        ElementsTransferProcess.Execute();

        //clear previous conditions
        std::string Conditions = "Conditions";
        rDestinationModelPart.Conditions().clear();
        if( rCompositeConditions ){
          this->TransferBoundaryConditions(rOriginModelPart, rDestinationModelPart);
        }
        else{
          TransferEntitiesBetweenModelPartsProcess ConditionsTransferProcess = TransferEntitiesBetweenModelPartsProcess(rDestinationModelPart, rOriginModelPart, Conditions, rTransferFlags);
          ConditionsTransferProcess.Execute();
        }

        //assigning Properties
        rDestinationModelPart.SetProperties( rOriginModelPart.pProperties() );

        //assigning Tables
	rDestinationModelPart.Tables() = rOriginModelPart.Tables();

        Communicator::Pointer pComm = rOriginModelPart.GetCommunicator().Create();
        rDestinationModelPart.SetCommunicator(pComm);

	if( GetEchoLevel() >= 1 )
	  std::cout<<" DestinationModelPart "<<rDestinationModelPart<<std::endl;
        
        KRATOS_CATCH("")
    }


    void TransferBoundaryConditions(ModelPart& rOriginModelPart,
                                    ModelPart& rDestinationModelPart)
    {
        KRATOS_TRY

        std::vector<Flags> BoundaryFlags;
        BoundaryFlags.push_back(BOUNDARY);
        std::string Conditions = "Conditions";
        TransferEntitiesBetweenModelPartsProcess ConditionsTransferProcess = TransferEntitiesBetweenModelPartsProcess(rDestinationModelPart, rOriginModelPart, Conditions, BoundaryFlags);
        ConditionsTransferProcess.Execute();

        KRATOS_CATCH("")
    }

  
    bool MatchFlags(const Element::Pointer& pElement, const std::vector<Flags>& rFlags)
    {

      for(unsigned int i = 0; i<rFlags.size(); i++)
	{
	  if( pElement->IsNot(rFlags[i]) )
	    return false;
	}

      return true;	  
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
