//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_TRANSFER_COMPUTING_MODEL_PART_ENTITIES_PROCESS_H_INCLUDED)
#define  KRATOS_TRANSFER_COMPUTING_MODEL_PART_ENTITIES_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_processes/transfer_entities_between_model_parts_process.hpp"
#include "solid_mechanics_application_variables.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Process to transfer model part entities to a computing submodelpart

class TransferComputingModelPartEntitiesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TransferComputingModelPartEntitiesProcess
    KRATOS_CLASS_POINTER_DEFINITION(TransferComputingModelPartEntitiesProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    TransferComputingModelPartEntitiesProcess(ModelPart& rModelPart,
                                      Parameters rParameters
                                      ) : Process(), mrModelPart(rModelPart)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
            {
               "model_part_name": "computing_sub_model_part",
               "transfer_flags": [],
               "composite_conditions": False,
               "generate_entities": []
            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        for(unsigned int i=0; i<rParameters["entities_transfer_flags"].size(); ++i)
        {
          mTransferFlags.push_back(KratosComponents<Flags>::Get(rParameters["entities_transfer_flags"][i].GetString()));
        }

        mCompositeConditions = rParameters["composite_conditions"].GetBool();


        for(unsigned int i=0; i<rParameters["generate_entities"].size(); ++i)
        {
          Parameters EntitiesGeneration = rParameters["generate_entities"][i];
          this->AddToEntitiesGenerationList(EntitiesGeneration);
        }

        KRATOS_CATCH("")
    }


    /// Destructor.
    ~TransferComputingModelPartEntitiesProcess() override {}


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


    /// Execute method is used to execute the TransferComputingModelPartEntitiesProcess algorithms.
    void Execute() override
    {
      KRATOS_TRY

      if( mElementGenerationList.size() != 0 || mConditionGenerationList.size() != 0 ){
        this->GenerateAndTransferEntities(*(mrModelPart.GetParentModelPart()), mrModelPart, mElementGenerationList, mConditionGenerationList, mCompositeConditions);
      }
      else{
        this->TransferEntities(*(mrModelPart.GetParentModelPart()), mrModelPart, mTransferFlags, mCompositeConditions);
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
        return "TransferComputingModelPartEntitiesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "TransferComputingModelPartEntitiesProcess";
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
    TransferComputingModelPartEntitiesProcess(TransferComputingModelPartEntitiesProcess const& rOther);

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


    template<class EntityType>
    class EntityGeneration
    {
     private:
      // member variables
      const EntityType* mpEntityType;

     public:

      //variables
      std::vector<Flags> OriginEntityFlags;

      // set methods
      void SetEntityType(const std::string& rName) {mpEntityType = &(KratosComponents<EntityType>::Get(rName));};

      void SetEntityType(const EntityType& rEntityType) { mpEntityType = &rEntityType; };

      // get methods
      const EntityType& GetEntityType() const { return *mpEntityType; };
    };


    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;

    bool mCompositeConditions;

    std::vector<Flags>  mTransferFlags;

    std::vector<EntityGeneration<Element> >  mElementGenerationList;

    std::vector<EntityGeneration<Condition> >  mConditionGenerationList;


    ///@}
    ///@name Private Operators
    ///@{
    ///@}
    ///@name Private Operations
    ///@{


    //**********************************************************************************************
    //**********************************************************************************************

    void AddToEntitiesGenerationList(Parameters& rParameters)
    {
        Parameters default_parameters( R"(
            {
                "entity_type": "Element",
                "creation_options" :{
                     "original_flags" : [],
                     "entity_kratos_type": "ElementType"
                }

            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        if( rParameters["entity_type"].GetString() == "Element" ){
          EntityGeneration<Element> Entity;
          for(unsigned int i=0; i<rParameters["creation_options"]["original_flags"].size(); ++i)
          {
            Entity.OriginEntityFlags.push_back(KratosComponents<Flags>::Get(rParameters["creation_options"]["original_flags"][i].GetString()));
          }
          Entity.SetEntityType(rParameters["creation_options"]["entity_kratos_type"].GetString());
          mElementGenerationList.push_back( Entity );
        }
        else if( rParameters["entity_type"].GetString() == "Condition" ){
          EntityGeneration<Condition> Entity;
          for(unsigned int i=0; i<rParameters["creation_options"]["original_flags"].size(); ++i)
          {
            Entity.OriginEntityFlags.push_back(KratosComponents<Flags>::Get(rParameters["creation_options"]["original_flags"][i].GetString()));
          }
          Entity.SetEntityType(rParameters["creation_options"]["entity_kratos_type"].GetString());
          mConditionGenerationList.push_back( Entity );
        }
    }

    //**********************************************************************************************
    //**********************************************************************************************

    ///This function fills the @param DestinationModelPart using the data obtained from @param rOriginModelPart
    ///the elements and conditions of the rDestinationModelPart part use the same connectivity (and id) as the
    ///OriginModelPart but their type is determined by @param rReferenceElement and @param rReferenceBoundaryCondition

    void GenerateAndTransferEntities(ModelPart& rOriginModelPart,
                                     ModelPart& rDestinationModelPart,
                                     const std::vector<EntityGeneration<Element> >&  rElementGenerationList,
                                     const std::vector<EntityGeneration<Condition> >&  rConditionGenerationList,
                                     bool& rCompositeConditions)
    {
        KRATOS_TRY

        //clear previous nodes
        std::string Nodes = "Nodes";
        rDestinationModelPart.Nodes().clear();
        TransferEntitiesBetweenModelPartsProcess NodesTransferProcess(rDestinationModelPart, rOriginModelPart, Nodes);
        NodesTransferProcess.Execute();

        //clear previous elements
        rDestinationModelPart.Elements().clear();

        //generating the elements
        for(std::vector<EntityGeneration<Element> >::const_iterator eg = rElementGenerationList.begin(); eg != rElementGenerationList.end(); ++eg)
        {
          for (ModelPart::ElementsContainerType::iterator ie = rOriginModelPart.ElementsBegin(); ie != rOriginModelPart.ElementsEnd(); ++ie)
          {
            if (this->MatchFlags(*(ie.base()),eg->OriginEntityFlags))
            {
            Properties::Pointer pProperties = ie->pGetProperties();
            Element::Pointer pElement = eg->GetEntityType().Create(ie->Id(), ie->GetGeometry(), pProperties);

	    //set origin element as pointer
	    WeakPointerVector< Element > MasterElements;
	    MasterElements.push_back(Element::WeakPointer( *(ie.base()) ) );

	    pElement->SetValue(MASTER_ELEMENTS, MasterElements);

            rDestinationModelPart.Elements().push_back(pElement);
            }
          }
        }

        //clear previous conditions
        std::string Conditions = "Conditions";
        rDestinationModelPart.Conditions().clear();

        if( rCompositeConditions ){
          this->TransferBoundaryConditions(rOriginModelPart, rDestinationModelPart);
        }

        for(std::vector<EntityGeneration<Condition> >::const_iterator cg = rConditionGenerationList.begin(); cg != rConditionGenerationList.end(); ++cg)
        {
          //generating contact conditions
          for (ModelPart::ConditionsContainerType::iterator ic = rOriginModelPart.ConditionsBegin(); ic != rOriginModelPart.ConditionsEnd(); ++ic)
          {
            if (this->MatchFlags(*(ic.base()),cg->OriginEntityFlags))
            {
              Properties::Pointer pProperties = ic->pGetProperties();
              Condition::Pointer pCondition = cg->GetEntityType().Create(ic->Id(), ic->GetGeometry(), pProperties);

              //set mechanical variables to contact conditions:
              pCondition->Data() = ic->Data();

              rDestinationModelPart.Conditions().push_back(pCondition);
            }
          }
        }

        //assigning Properties
        rDestinationModelPart.SetProperties( rOriginModelPart.pProperties() );

        //assigning Tables
	rDestinationModelPart.Tables() = rOriginModelPart.Tables();

        Communicator::Pointer pComm = rOriginModelPart.GetCommunicator().Create();
        rDestinationModelPart.SetCommunicator(pComm);

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
        TransferEntitiesBetweenModelPartsProcess NodesTransferProcess(rDestinationModelPart, rOriginModelPart, Nodes);
        NodesTransferProcess.Execute();

        //clear previous elements
        std::string Elements = "Elements";
        rDestinationModelPart.Elements().clear();
        TransferEntitiesBetweenModelPartsProcess ElementsTransferProcess(rDestinationModelPart, rOriginModelPart, Elements, rTransferFlags);
        ElementsTransferProcess.Execute();

        //clear previous conditions
        std::string Conditions = "Conditions";
        rDestinationModelPart.Conditions().clear();
        if( rCompositeConditions ){
          this->TransferBoundaryConditions(rOriginModelPart, rDestinationModelPart);
        }
        else{
          TransferEntitiesBetweenModelPartsProcess ConditionsTransferProcess(rDestinationModelPart, rOriginModelPart, Conditions, rTransferFlags);
          ConditionsTransferProcess.Execute();
        }

        //assigning Properties
        rDestinationModelPart.SetProperties( rOriginModelPart.pProperties() );

        //assigning Tables
	rDestinationModelPart.Tables() = rOriginModelPart.Tables();

        Communicator::Pointer pComm = rOriginModelPart.GetCommunicator().Create();
        rDestinationModelPart.SetCommunicator(pComm);

        KRATOS_CATCH("")
    }


    void TransferBoundaryConditions(ModelPart& rOriginModelPart,
                                    ModelPart& rDestinationModelPart)
    {
        KRATOS_TRY

        std::vector<Flags> BoundaryFlags;
        BoundaryFlags.push_back(BOUNDARY);
        std::string Conditions = "Conditions";
        TransferEntitiesBetweenModelPartsProcess ConditionsTransferProcess(rDestinationModelPart, rOriginModelPart, Conditions, BoundaryFlags);
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


    bool MatchFlags(const Condition::Pointer& pCondition, const std::vector<Flags>& rFlags)
    {

      for(unsigned int i = 0; i<rFlags.size(); i++)
	{
	  if( pCondition->IsNot(rFlags[i]) )
	    return false;
	}

      return true;
    }

    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    TransferComputingModelPartEntitiesProcess& operator=(TransferComputingModelPartEntitiesProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class TransferComputingModelPartEntitiesProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  TransferComputingModelPartEntitiesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const TransferComputingModelPartEntitiesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_TRANSFER_COMPUTING_MODEL_PART_ENTITIES_PROCESS_H_INCLUDED  defined
