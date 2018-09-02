//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_TRANSFER_SOLVING_MODEL_PART_ENTITIES_PROCESS_H_INCLUDED)
#define  KRATOS_TRANSFER_SOLVING_MODEL_PART_ENTITIES_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_processes/transfer_entities_between_model_parts_process.hpp"
#include "solid_mechanics_application_variables.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Process to transfer model part entities to a solving submodelpart

class TransferSolvingModelPartEntitiesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3>   NodeType;
    typedef Condition ConditionType;
    typedef Element   ElementType;

    /// Pointer definition of TransferSolvingModelPartEntitiesProcess
    KRATOS_CLASS_POINTER_DEFINITION(TransferSolvingModelPartEntitiesProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    TransferSolvingModelPartEntitiesProcess(ModelPart& rModelPart,
                                            Parameters rParameters
                                            ) : Process(), mrModelPart(rModelPart)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
            {
               "model_part_name": "new_computing_domain",
               "assign_flags": [],
               "composite_conditions": false,
               "transfer_entities": [],
               "generate_entities": []
            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mSolvingModelPartName = rParameters["model_part_name"].GetString();
        if(!mrModelPart.HasSubModelPart(mSolvingModelPartName))
          mrModelPart.CreateSubModelPart(mSolvingModelPartName);


        for(unsigned int i = 0; i<rParameters["assign_flags"].size(); i++)
          mrModelPart.GetSubModelPart(mSolvingModelPartName).Set(KratosComponents<Flags>::Get(rParameters["assign_flags"][i].GetString()));

        mCompositeConditions = rParameters["composite_conditions"].GetBool();

        for(unsigned int i=0; i<rParameters["transfer_entities"].size(); ++i)
        {
          Parameters EntitiesTransfer = rParameters["transfer_entities"][i];
          this->AddToEntitiesTransferList(EntitiesTransfer);
        }

        for(unsigned int i=0; i<rParameters["generate_entities"].size(); ++i)
        {
          Parameters EntitiesGeneration = rParameters["generate_entities"][i];
          this->AddToEntitiesGenerationList(EntitiesGeneration);
        }

        KRATOS_CATCH("")
    }


    /// Destructor.
    ~TransferSolvingModelPartEntitiesProcess() override {}


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


    /// Execute method is used to execute the TransferSolvingModelPartEntitiesProcess algorithms.
    void Execute() override
    {
      KRATOS_TRY

      ModelPart& SolvingModelPart = mrModelPart.GetSubModelPart(mSolvingModelPartName);

      //clear previous nodes
      SolvingModelPart.Nodes().clear();
      //clear previous elements
      SolvingModelPart.Elements().clear();
      //clear previous conditions
      SolvingModelPart.Conditions().clear();

      if( mElementGenerationList.size() != 0 || mConditionGenerationList.size() != 0 ){
        this->GenerateAndTransferEntities(SolvingModelPart, mElementGenerationList, mConditionGenerationList);
      }

      if( mNodeTransferList.size()!=0 || mElementTransferList.size() != 0 || mConditionTransferList.size() != 0 ){
        this->TransferEntities(SolvingModelPart, mNodeTransferList, mElementTransferList, mConditionTransferList);
      }

      if( mCompositeConditions ){
        this->TransferBoundaryConditions(SolvingModelPart,mrModelPart);
      }

      //assigning Properties
      SolvingModelPart.SetProperties( mrModelPart.pProperties() );

      //assigning Tables
      SolvingModelPart.Tables() = mrModelPart.Tables();

      Communicator::Pointer pComm = mrModelPart.GetCommunicator().Create();
      SolvingModelPart.SetCommunicator(pComm);

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
        return "TransferSolvingModelPartEntitiesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "TransferSolvingModelPartEntitiesProcess";
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
    TransferSolvingModelPartEntitiesProcess(TransferSolvingModelPartEntitiesProcess const& rOther);

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
    class EntityTransfer
    {
     private:
      // member variables
      const EntityType* mpEntityType;

     public:

      //variables
      std::vector<ModelPart*>  OriginModelParts;
      std::vector<Flags>  TransferFlags;
      std::vector<Flags>  AssignFlags;

      // set methods
      void SetEntityType(const std::string& rName) {mpEntityType = &(KratosComponents<EntityType>::Get(rName));};

      void SetEntityType(const EntityType& rEntityType) {mpEntityType = &rEntityType;};

      // get methods
      const EntityType& GetEntityType() const {return *mpEntityType;};
    };


    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;

    std::string mSolvingModelPartName;

    bool mCompositeConditions;

    std::vector<EntityTransfer<NodeType> >  mNodeTransferList;

    std::vector<EntityTransfer<ElementType> >  mElementTransferList;

    std::vector<EntityTransfer<ConditionType> >  mConditionTransferList;

    std::vector<EntityTransfer<ElementType> >  mElementGenerationList;

    std::vector<EntityTransfer<ConditionType> >  mConditionGenerationList;

    ///@}
    ///@name Private Operators
    ///@{
    ///@}
    ///@name Private Operations
    ///@{


    //**********************************************************************************************
    //**********************************************************************************************

    void AddToEntitiesTransferList(Parameters& rParameters)
    {
        Parameters default_parameters( R"(
            {
                "origin_model_parts_list": [],
                "entity_type": "Element",
                "transfer_flags" : [],
                "assign_flags" : []
            }  )" );

        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        if( rParameters["entity_type"].GetString() == "Node" ){
          EntityTransfer<NodeType> Entity;
          for( unsigned int i=0; i<rParameters["origin_model_parts_list"].size(); ++i)
          {
            Entity.OriginModelParts.push_back(&mrModelPart.GetSubModelPart(rParameters["origin_model_parts_list"][i].GetString()));
          }
          for(unsigned int i=0; i<rParameters["transfer_flags"].size(); ++i)
          {
            Entity.TransferFlags.push_back(KratosComponents<Flags>::Get(rParameters["transfer_flags"][i].GetString()));
          }
          for(unsigned int i=0; i<rParameters["assign_flags"].size(); ++i)
          {
            Entity.AssignFlags.push_back(KratosComponents<Flags>::Get(rParameters["assign_flags"][i].GetString()));
          }
          mNodeTransferList.push_back( Entity );
        }
        else if( rParameters["entity_type"].GetString() == "Element" ){
          EntityTransfer<ElementType> Entity;
          for( unsigned int i=0; i<rParameters["origin_model_parts_list"].size(); ++i)
          {
            Entity.OriginModelParts.push_back(&mrModelPart.GetSubModelPart(rParameters["origin_model_parts_list"][i].GetString()));
          }
          for(unsigned int i=0; i<rParameters["transfer_flags"].size(); ++i)
          {
            Entity.TransferFlags.push_back(KratosComponents<Flags>::Get(rParameters["transfer_flags"][i].GetString()));
          }
          for(unsigned int i=0; i<rParameters["assign_flags"].size(); ++i)
          {
            Entity.AssignFlags.push_back(KratosComponents<Flags>::Get(rParameters["assign_flags"][i].GetString()));
          }
          mElementTransferList.push_back( Entity );
        }
        else if( rParameters["entity_type"].GetString() == "Condition" )
        {
          EntityTransfer<ConditionType> Entity;
          for( unsigned int i=0; i<rParameters["origin_model_parts_list"].size(); ++i)
          {
            Entity.OriginModelParts.push_back(&mrModelPart.GetSubModelPart(rParameters["origin_model_parts_list"][i].GetString()));
          }
          for(unsigned int i=0; i<rParameters["transfer_flags"].size(); ++i)
          {
            Entity.TransferFlags.push_back(KratosComponents<Flags>::Get(rParameters["transfer_flags"][i].GetString()));
          }
          for(unsigned int i=0; i<rParameters["assign_flags"].size(); ++i)
          {
            Entity.AssignFlags.push_back(KratosComponents<Flags>::Get(rParameters["assign_flags"][i].GetString()));
          }
          mConditionTransferList.push_back( Entity );
        }

    }


    //**********************************************************************************************
    //**********************************************************************************************

    void AddToEntitiesGenerationList(Parameters& rParameters)
    {
        Parameters default_parameters( R"(
            {
                "origin_model_parts_list": [],
                "entity_type": "Element",
                "transfer_flags": [],
                "assign_flags": [],
                "entity_kratos_type": "ElementType"
            }  )" );

        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        if( rParameters["entity_type"].GetString() == "Element" ){
          EntityTransfer<ElementType> Entity;
          for( unsigned int i=0; i<rParameters["origin_model_parts_list"].size(); ++i)
          {
            Entity.OriginModelParts.push_back(&mrModelPart.GetSubModelPart(rParameters["origin_model_parts_list"][i].GetString()));
          }
          for(unsigned int i=0; i<rParameters["transfer_flags"].size(); ++i)
          {
            Entity.TransferFlags.push_back(KratosComponents<Flags>::Get(rParameters["transfer_flags"][i].GetString()));
          }
          for(unsigned int i=0; i<rParameters["assign_flags"].size(); ++i)
          {
            Entity.AssignFlags.push_back(KratosComponents<Flags>::Get(rParameters["assign_flags"][i].GetString()));
          }
          Entity.SetEntityType(rParameters["entity_kratos_type"].GetString());
          mElementGenerationList.push_back( Entity );
        }
        else if( rParameters["entity_type"].GetString() == "Condition" )
        {
          EntityTransfer<ConditionType> Entity;
          for( unsigned int i=0; i<rParameters["origin_model_parts_list"].size(); ++i)
          {
            Entity.OriginModelParts.push_back(&mrModelPart.GetSubModelPart(rParameters["origin_model_parts_list"][i].GetString()));
          }
          for(unsigned int i=0; i<rParameters["transfer_flags"].size(); ++i)
          {
            Entity.TransferFlags.push_back(KratosComponents<Flags>::Get(rParameters["transfer_flags"][i].GetString()));
          }
          for(unsigned int i=0; i<rParameters["assign_flags"].size(); ++i)
          {
            Entity.AssignFlags.push_back(KratosComponents<Flags>::Get(rParameters["assign_flags"][i].GetString()));
          }
          Entity.SetEntityType(rParameters["entity_kratos_type"].GetString());
          mConditionGenerationList.push_back( Entity );
        }
    }

    //**********************************************************************************************
    //**********************************************************************************************

    ///This function fills the @param DestinationModelPart using the data obtained from @param rOriginModelPart
    ///the elements and conditions of the rDestinationModelPart part use the same connectivity (and id) as the
    ///OriginModelPart but their type is determined by @param rReferenceElement and @param rReferenceBoundaryCondition

    void GenerateAndTransferEntities(ModelPart& rDestinationModelPart,
                                     std::vector<EntityTransfer<ElementType> >& rElementGenerationList,
                                     std::vector<EntityTransfer<ConditionType> >& rConditionGenerationList)
    {
        KRATOS_TRY

        //generate elements
        for(std::vector<EntityTransfer<Element> >::iterator i_entity = rElementGenerationList.begin(); i_entity != rElementGenerationList.end(); ++i_entity)
        {
          for(std::vector<ModelPart*>::iterator i_part = i_entity->OriginModelParts.begin(); i_part != i_entity->OriginModelParts.end(); ++i_part)
          {
            for(ModelPart::ElementsContainerType::iterator i_elem = (*i_part)->ElementsBegin(); i_elem != (*i_part)->ElementsEnd(); ++i_elem)
            {
              if (this->MatchFlags(*(i_elem.base()),i_entity->TransferFlags))
              {
                Properties::Pointer pProperties = i_elem->pGetProperties();
                Element::Pointer pElement = i_entity->GetEntityType().Create(i_elem->Id(), i_elem->GetGeometry(), pProperties);

                //set origin element as pointer
                WeakPointerVector< Element > MasterElements;
                MasterElements.push_back(Element::WeakPointer( *(i_elem.base()) ) );

                pElement->SetValue(MASTER_ELEMENTS, MasterElements);

                rDestinationModelPart.Elements().push_back(pElement);
              }
            }
          }
        }

        //generate conditions
        for(std::vector<EntityTransfer<Condition> >::iterator i_entity = rConditionGenerationList.begin(); i_entity != rConditionGenerationList.end(); ++i_entity)
        {
          for(std::vector<ModelPart*>::iterator i_part = i_entity->OriginModelParts.begin(); i_part != i_entity->OriginModelParts.end(); ++i_part)
          {
            //generate contact conditions
            for (ModelPart::ConditionsContainerType::iterator i_cond = (*i_part)->ConditionsBegin(); i_cond != (*i_part)->ConditionsEnd(); ++i_cond)
            {
              if (this->MatchFlags(*(i_cond.base()),i_entity->TransferFlags))
              {
                Properties::Pointer pProperties = i_cond->pGetProperties();
                Condition::Pointer pCondition = i_entity->GetEntityType().Create(i_cond->Id(), i_cond->GetGeometry(), pProperties);

                //set mechanical variables to contact conditions:
                pCondition->Data() = i_cond->Data();
               
                rDestinationModelPart.Conditions().push_back(pCondition);
              }
            }
          }
        }

        KRATOS_CATCH("")
    }


    void TransferEntities(ModelPart& rDestinationModelPart,
                          std::vector<EntityTransfer<NodeType> >& rNodeTransferList,
                          std::vector<EntityTransfer<ElementType> >&  rElementTransferList,
                          std::vector<EntityTransfer<ConditionType> >&  rConditionTransferList)
    {
        KRATOS_TRY


        std::string Nodes = "Nodes";

        //transfer nodes
        for(std::vector<EntityTransfer<NodeType> >::iterator i_entity = rNodeTransferList.begin(); i_entity != rNodeTransferList.end(); ++i_entity)
        {
          for(std::vector<ModelPart*>::iterator i_part = i_entity->OriginModelParts.begin(); i_part != i_entity->OriginModelParts.end(); ++i_part)
          {
            TransferEntitiesBetweenModelPartsProcess NodesTransferProcess(rDestinationModelPart, *(*i_part), Nodes, i_entity->TransferFlags, i_entity->AssignFlags);
            NodesTransferProcess.Execute();
          }
        }

        std::string Elements = "Elements";

        //transfer elements
        for(std::vector<EntityTransfer<ElementType> >::iterator i_entity = rElementTransferList.begin(); i_entity != rElementTransferList.end(); ++i_entity)
        {
          for(std::vector<ModelPart*>::iterator i_part = i_entity->OriginModelParts.begin(); i_part != i_entity->OriginModelParts.end(); ++i_part)
          {
            TransferEntitiesBetweenModelPartsProcess ElementsTransferProcess(rDestinationModelPart, *(*i_part), Elements, i_entity->TransferFlags, i_entity->AssignFlags);
            ElementsTransferProcess.Execute();
          }
        }

        std::string Conditions = "Conditions";

        //transfer conditions
        for(std::vector<EntityTransfer<ConditionType> >::iterator i_entity = rConditionTransferList.begin(); i_entity != rConditionTransferList.end(); ++i_entity)
        {
          for(std::vector<ModelPart*>::iterator i_part = i_entity->OriginModelParts.begin(); i_part != i_entity->OriginModelParts.end(); ++i_part)
          {
            TransferEntitiesBetweenModelPartsProcess ConditionsTransferProcess(rDestinationModelPart, *(*i_part), Conditions, i_entity->TransferFlags, i_entity->AssignFlags);
            ConditionsTransferProcess.Execute();
          }
        }

        KRATOS_CATCH("")
    }


    void TransferBoundaryConditions(ModelPart& rDestinationModelPart,
                                    ModelPart& rOriginModelPart)
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
    TransferSolvingModelPartEntitiesProcess& operator=(TransferSolvingModelPartEntitiesProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class TransferSolvingModelPartEntitiesProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  TransferSolvingModelPartEntitiesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const TransferSolvingModelPartEntitiesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_TRANSFER_SOLVING_MODEL_PART_ENTITIES_PROCESS_H_INCLUDED  defined
