// System includes
#include <string>
#include <vector>
#include <sstream>

// External includes

// Project includes
#include "includes/define.h"
#include "modeler/combine_model_part_modeler.h"
#include "utilities/variable_utils.h"

namespace Kratos
{

CombineModelPartModeler::CombineModelPartModeler(
    Model& rModel,
    Parameters ModelerParameters
    ) : Modeler(rModel, ModelerParameters),
        mpModel(&rModel),
        mParameters(ModelerParameters)
{
    mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

/***********************************************************************************/
/***********************************************************************************/

Modeler::Pointer CombineModelPartModeler::Create(
    Model& rModel,
    const Parameters ModelParameters
    ) const
{
    return Kratos::make_shared<CombineModelPartModeler>(rModel, ModelParameters);
}

/***********************************************************************************/
/***********************************************************************************/

void CombineModelPartModeler::SetupModelPart()
{
    KRATOS_TRY;
    const auto& r_new_model_part_name = mParameters["combined_model_part_name"].GetString();

    auto& r_combined_model_part = mpModel->HasModelPart(r_new_model_part_name) ?
        mpModel->GetModelPart(r_new_model_part_name) : mpModel->CreateModelPart(r_new_model_part_name);

    this->ResetModelPart(r_combined_model_part);

    this->CheckOriginModelPartsAndAssignRoot();

    this->CopyCommonData(r_combined_model_part);

    this->DuplicateMesh();

    this->CreateCommunicators();

    this->PopulateCommunicators();

    this->CreateSubModelParts();

    KRATOS_CATCH("Failure in SetupModelPart");
}

/***********************************************************************************/
/***********************************************************************************/

void CombineModelPartModeler::ResetModelPart(ModelPart& rCombinedModelPart) const
{
    VariableUtils().SetFlag(TO_ERASE, true, rCombinedModelPart.Nodes());
    VariableUtils().SetFlag(TO_ERASE, true, rCombinedModelPart.Elements());
    VariableUtils().SetFlag(TO_ERASE, true, rCombinedModelPart.Conditions());

    rCombinedModelPart.RemoveNodesFromAllLevels(TO_ERASE);
    rCombinedModelPart.RemoveElementsFromAllLevels(TO_ERASE);
    rCombinedModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
}


/***********************************************************************************/
/***********************************************************************************/

void CombineModelPartModeler::CheckOriginModelPartsAndAssignRoot()
{
    KRATOS_TRY;

    Parameters model_part_list = mParameters["model_part_list"];

    KRATOS_ERROR_IF(model_part_list.size() == 0) <<
        "\"model_part_list\" is empty!\n";

    const std::string& first_origin_model_pat = model_part_list[0]["origin_model_part"].GetString();
    mpOriginRootModelPart =  &(mpModel->GetModelPart(first_origin_model_pat).GetRootModelPart());

    for (unsigned int i = 1; i < model_part_list.size(); i++) {
        ModelPart& r_origin_model_part = mpModel->GetModelPart(model_part_list[i]["origin_model_part"].GetString());
        if (r_origin_model_part.GetRootModelPart().FullName() != mpOriginRootModelPart->FullName()) {
            KRATOS_ERROR << "The origin model part \"" << r_origin_model_part.FullName() <<
                "\" does not have the same root as the rest of origin model parts: \"" <<
                mpOriginRootModelPart->FullName() << "\".\n";
        }
    }

    KRATOS_CATCH("Failure in CheckOriginModelPartsAndAssignRoot");

}

/***********************************************************************************/
/***********************************************************************************/

void CombineModelPartModeler::CopyCommonData(
    ModelPart& rCombinedModelPart) const
{
    KRATOS_TRY;

    if (!rCombinedModelPart.IsSubModelPart()) {
        rCombinedModelPart.SetNodalSolutionStepVariablesList(mpOriginRootModelPart->pGetNodalSolutionStepVariablesList());
        rCombinedModelPart.SetBufferSize( mpOriginRootModelPart->GetBufferSize() );
    }

    rCombinedModelPart.SetProcessInfo( mpOriginRootModelPart->pGetProcessInfo() );
    rCombinedModelPart.PropertiesArray() = mpOriginRootModelPart->PropertiesArray();
    rCombinedModelPart.Tables() = mpOriginRootModelPart->Tables();
    KRATOS_CATCH("Failure in CopyCommonData");

}

/***********************************************************************************/
/***********************************************************************************/

void CombineModelPartModeler::DuplicateMesh() const
{
    KRATOS_TRY;
    const std::string element_name = mParameters["element_name"].GetString();
    const std::string condition_name = mParameters["condition_name"].GetString();

    // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
    // so that an error is thrown if they don't exist
    KRATOS_ERROR_IF(element_name != "" && !KratosComponents<Element>::Has(element_name)) << "Element name not found in KratosComponents< Element > -- name is " << element_name << std::endl;
    KRATOS_ERROR_IF(condition_name != "" && !KratosComponents<Condition>::Has(condition_name)) << "Condition name not found in KratosComponents< Condition > -- name is " << condition_name << std::endl;

    const auto& reference_condition = KratosComponents<Condition>::Get(condition_name);
    const auto& reference_element = KratosComponents<Element>::Get(element_name);

    Parameters model_part_list = mParameters["model_part_list"];

    for (unsigned int i = 0; i < model_part_list.size(); i++) {
        ModelPart& r_origin_model_part = mpModel->GetModelPart(model_part_list[i]["origin_model_part"].GetString());
        const std::string& destination_model_part_name = model_part_list[i]["destination_model_part"].GetString();
        ModelPart* p_destination_model_part;
        if (mpModel->HasModelPart(destination_model_part_name)) {
            p_destination_model_part = &mpModel->GetModelPart(destination_model_part_name);
        } else {
            p_destination_model_part = &mpModel->CreateModelPart(destination_model_part_name);
        }
        p_destination_model_part->AddNodes(r_origin_model_part.NodesBegin(), r_origin_model_part.NodesEnd());
        this->DuplicateElements(r_origin_model_part, *p_destination_model_part, reference_element);
        this->DuplicateConditions(r_origin_model_part, *p_destination_model_part, reference_condition);
    }
    KRATOS_CATCH("Failure in DuplicateElementsAndConditions");
}

/***********************************************************************************/
/***********************************************************************************/

void CombineModelPartModeler::CreateCommunicators()
{
    KRATOS_TRY;
    Parameters model_part_list = mParameters["model_part_list"];

    for (unsigned int i = 0; i < mParameters["model_part_list"].size(); i++) {
        ModelPart& r_origin_model_part = mpModel->GetModelPart(model_part_list[i]["origin_model_part"].GetString());
        ModelPart& r_destination_model_part = mpModel->GetModelPart(model_part_list[i]["destination_model_part"].GetString());
        Communicator& r_reference_comm = r_origin_model_part.GetCommunicator();
        if (r_origin_model_part.GetCommunicator().IsDistributed()) {
            Communicator::Pointer p_destination_comm = r_reference_comm.Create();
            p_destination_comm->SetNumberOfColors(r_reference_comm.GetNumberOfColors());
            p_destination_comm->NeighbourIndices() = r_reference_comm.NeighbourIndices();
            r_destination_model_part.SetCommunicator( p_destination_comm );

            ModelPart* p_current_model_part = &r_destination_model_part;
            while (p_current_model_part->IsSubModelPart()) {
                p_current_model_part = &(p_current_model_part->GetParentModelPart());
                Communicator& r_parent_reference_comm = mpOriginRootModelPart->GetCommunicator();
                Communicator::Pointer p_parent_destination_comm = r_parent_reference_comm.Create();
                p_parent_destination_comm->SetNumberOfColors(r_parent_reference_comm.GetNumberOfColors());
                p_parent_destination_comm->NeighbourIndices() = r_parent_reference_comm.NeighbourIndices();
                p_current_model_part->SetCommunicator( p_parent_destination_comm );

            }
        }
    }
    KRATOS_CATCH("Failure in CreateCommunicators");
}

/***********************************************************************************/
/***********************************************************************************/

void CombineModelPartModeler::PopulateCommunicators()
{
    Parameters model_part_list = mParameters["model_part_list"];

    for (unsigned int i = 0; i < mParameters["model_part_list"].size(); i++) {
        ModelPart& r_origin_model_part = mpModel->GetModelPart(model_part_list[i]["origin_model_part"].GetString());
        ModelPart& r_destination_model_part = mpModel->GetModelPart(model_part_list[i]["destination_model_part"].GetString());
        ModelPart* p_current_model_part = &r_destination_model_part;
        Communicator& r_origin_communicator = r_origin_model_part.GetCommunicator();
        if (r_origin_model_part.GetCommunicator().IsDistributed()) {
            bool keep_going_up = true;
            while (keep_going_up) {
                Communicator& r_destination_communicator = p_current_model_part->GetCommunicator();
                PopulateLocalMesh(
                    r_origin_communicator,
                    r_destination_communicator,
                    *p_current_model_part
                );
                keep_going_up = p_current_model_part->IsSubModelPart();
                p_current_model_part = &p_current_model_part->GetParentModelPart();
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CombineModelPartModeler::DuplicateElements(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart,
    const Element& rReferenceElement
    ) const
{
    // Generate the elements
    ModelPart::ElementsContainerType temp_elements;
    temp_elements.reserve(rOriginModelPart.NumberOfElements());
    for (auto i_elem = rOriginModelPart.ElementsBegin(); i_elem != rOriginModelPart.ElementsEnd(); ++i_elem) {
        if (!rDestinationModelPart.GetRootModelPart().HasElement(i_elem->Id())) {
            Properties::Pointer properties = i_elem->pGetProperties();

            // Reuse the geometry of the old element (to save memory)
            Element::Pointer p_element = rReferenceElement.Create(i_elem->Id(), i_elem->pGetGeometry(), properties);
            temp_elements.push_back(p_element);
        } else {
            temp_elements.push_back(rDestinationModelPart.GetRootModelPart().pGetElement(i_elem->Id()));
        }
    }

    rDestinationModelPart.AddElements(temp_elements.begin(), temp_elements.end());
}

/***********************************************************************************/
/***********************************************************************************/

void CombineModelPartModeler::DuplicateConditions(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart,
    const Condition& rReferenceBoundaryCondition
    ) const
{
    // Generate the conditions
    ModelPart::ConditionsContainerType temp_conditions;
    temp_conditions.reserve(rOriginModelPart.NumberOfConditions());
    for (auto i_cond = rOriginModelPart.ConditionsBegin(); i_cond != rOriginModelPart.ConditionsEnd(); ++i_cond) {
        if (!rDestinationModelPart.GetRootModelPart().HasCondition(i_cond->Id())) {
            Properties::Pointer properties = i_cond->pGetProperties();
            // Reuse the geometry of the old element (to save memory)
            Condition::Pointer p_condition = rReferenceBoundaryCondition.Create(i_cond->Id(), i_cond->pGetGeometry(), properties);

            temp_conditions.push_back(p_condition);
        } else {
            temp_conditions.push_back(rDestinationModelPart.GetRootModelPart().pGetCondition(i_cond->Id()));
        }
    }

    rDestinationModelPart.AddConditions(temp_conditions.begin(), temp_conditions.end());
}

/***********************************************************************************/
/***********************************************************************************/

void CombineModelPartModeler::PopulateLocalMesh(
    Communicator& rReferenceComm,
    Communicator& rDestinationComm,
    ModelPart& rDestinationModelPart
) const
{
    if (rReferenceComm.IsDistributed()) {

        ModelPart::NodesContainerType& rDestinationLocalNodes = rDestinationComm.LocalMesh().Nodes();
        rDestinationLocalNodes.reserve(rReferenceComm.LocalMesh().NumberOfNodes());
        auto local_nodes_mutable_pass = rDestinationLocalNodes.GetMutablePass();
        for (auto p_node = rReferenceComm.LocalMesh().Nodes().ptr_begin(); p_node != rReferenceComm.LocalMesh().Nodes().ptr_end(); ++p_node) {
            local_nodes_mutable_pass.push_back(*p_node);
        }

        ModelPart::NodesContainerType& rDestinationInterfaceNodes = rDestinationComm.InterfaceMesh().Nodes();
        rDestinationInterfaceNodes.reserve(rReferenceComm.InterfaceMesh().NumberOfNodes());
        auto interface_nodes_mutable_pass = rDestinationInterfaceNodes.GetMutablePass();
        for (auto p_node = rReferenceComm.InterfaceMesh().Nodes().ptr_begin(); p_node != rReferenceComm.InterfaceMesh().Nodes().ptr_end(); ++p_node) {
            interface_nodes_mutable_pass.push_back(*p_node);
        }

        ModelPart::NodesContainerType& rDestinationGhostNodes = rDestinationComm.GhostMesh().Nodes();
        rDestinationGhostNodes.reserve(rReferenceComm.GhostMesh().NumberOfNodes());
        auto ghost_nodes_mutable_pass = rDestinationGhostNodes.GetMutablePass();
        for (auto p_node = rReferenceComm.GhostMesh().Nodes().ptr_begin(); p_node != rReferenceComm.GhostMesh().Nodes().ptr_end(); ++p_node) {
            ghost_nodes_mutable_pass.push_back(*p_node);
        }

        ModelPart::ConditionsContainerType& rDestinationLocalConditions = rDestinationComm.LocalMesh().Conditions();
        rDestinationLocalConditions.reserve(rDestinationModelPart.NumberOfConditions());
        auto conditions_mutable_pass = rDestinationLocalConditions.GetMutablePass();
        for (auto p_cond = rDestinationModelPart.Conditions().ptr_begin(); p_cond != rDestinationModelPart.Conditions().ptr_end(); ++p_cond) {
            conditions_mutable_pass.push_back(*p_cond);
        }

        ModelPart::ElementsContainerType& rDestinationLocalElements = rDestinationComm.LocalMesh().Elements();
        rDestinationLocalElements.reserve(rDestinationModelPart.NumberOfElements());
        auto elements_mutable_pass = rDestinationLocalElements.GetMutablePass();
        for (auto p_elem = rDestinationModelPart.Elements().ptr_begin(); p_elem != rDestinationModelPart.Elements().ptr_end(); ++p_elem) {
            elements_mutable_pass.push_back(*p_elem);
        }
    } else {
        ModelPart::NodesContainerType& rDestinationLocalNodes = rDestinationComm.LocalMesh().Nodes();
        rDestinationLocalNodes.reserve(rDestinationModelPart.NumberOfNodes());
        auto local_nodes_mutable_pass = rDestinationLocalNodes.GetMutablePass();
        for (auto p_node = rDestinationModelPart.Nodes().ptr_begin(); p_node != rDestinationModelPart.Nodes().ptr_end(); ++p_node) {
            local_nodes_mutable_pass.push_back(*p_node);
        }

        ModelPart::ConditionsContainerType& rDestinationLocalConditions = rDestinationComm.LocalMesh().Conditions();
        rDestinationLocalConditions.reserve(rDestinationModelPart.NumberOfConditions());
        auto conditions_mutable_pass = rDestinationLocalConditions.GetMutablePass();
        for (auto p_cond = rDestinationModelPart.Conditions().ptr_begin(); p_cond != rDestinationModelPart.Conditions().ptr_end(); ++p_cond) {
            conditions_mutable_pass.push_back(*p_cond);
        }

        ModelPart::ElementsContainerType& rDestinationLocalElements = rDestinationComm.LocalMesh().Elements();
        rDestinationLocalElements.reserve(rDestinationModelPart.NumberOfElements());
        auto elements_mutable_pass = rDestinationLocalElements.GetMutablePass();
        for (auto p_elem = rDestinationModelPart.Elements().ptr_begin(); p_elem != rDestinationModelPart.Elements().ptr_end(); ++p_elem) {
            elements_mutable_pass.push_back(*p_elem);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CombineModelPartModeler::DuplicateCommunicatorData(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart) const
{
    Communicator& rReferenceComm = rOriginModelPart.GetCommunicator();
    if (rReferenceComm.IsDistributed()) {
        Communicator::Pointer pDestinationComm = rReferenceComm.Create();
        pDestinationComm->SetNumberOfColors( rReferenceComm.GetNumberOfColors() );
        pDestinationComm->NeighbourIndices() = rReferenceComm.NeighbourIndices();
        rDestinationModelPart.SetCommunicator( pDestinationComm );

        auto& rDestinationComm = rDestinationModelPart.GetCommunicator();
        PopulateLocalMesh(
            rReferenceComm,
            rDestinationComm,
            rDestinationModelPart
        );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CombineModelPartModeler::CreateSubModelParts()
{
     Parameters model_part_list = mParameters["model_part_list"];

    for (unsigned int i = 0; i < mParameters["model_part_list"].size(); i++) {
        ModelPart& r_origin_model_part = mpModel->GetModelPart(model_part_list[i]["origin_model_part"].GetString());
        ModelPart& r_destination_model_part = mpModel->GetModelPart(model_part_list[i]["destination_model_part"].GetString());
        this->DuplicateSubModelParts(r_origin_model_part, r_destination_model_part);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CombineModelPartModeler::DuplicateSubModelParts(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart) const
{
    for (auto i_part = rOriginModelPart.SubModelPartsBegin(); i_part != rOriginModelPart.SubModelPartsEnd(); ++i_part) {
        if(!rDestinationModelPart.HasSubModelPart(i_part->Name())) {
            rDestinationModelPart.CreateSubModelPart(i_part->Name());
        }

        ModelPart& destination_part = rDestinationModelPart.GetSubModelPart(i_part->Name());

        destination_part.AddNodes(i_part->NodesBegin(), i_part->NodesEnd());

        std::vector<ModelPart::IndexType> ids;
        ids.reserve(i_part->Elements().size());

        // Execute only if we created elements in the destination
        if (rDestinationModelPart.NumberOfElements() > 0)
        {
            //adding by index
            for(auto it=i_part->ElementsBegin(); it!=i_part->ElementsEnd(); ++it) {
                if (!destination_part.HasElement(it->Id())) {
                    ids.push_back(it->Id());
                }
            }
            destination_part.AddElements(ids, 0); //adding by index
        }

        // Execute only if we created conditions in the destination
        if (rDestinationModelPart.NumberOfConditions() > 0)
        {
            ids.clear();
            for(auto it=i_part->ConditionsBegin(); it!=i_part->ConditionsEnd(); ++it) {
                if (!destination_part.HasCondition(it->Id())) {
                    ids.push_back(it->Id());
                }
            }
            destination_part.AddConditions(ids, 0);
        }

        // Duplicate the Communicator for this SubModelPart
        this->DuplicateCommunicatorData(*i_part, destination_part);

        // Recursively call this function to duplicate any child SubModelParts
        this->DuplicateSubModelParts(*i_part, destination_part);
    }
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters CombineModelPartModeler::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"({
        "combined_model_part_name" : "",
        "model_part_list"          : [],
        "element_name"             : "Element3D4N",
        "condition_name"           : "SurfaceCondition3D3N"
    })");
    return default_parameters;
}

} // namespace Kratos
