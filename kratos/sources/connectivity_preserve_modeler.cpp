//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#include "modeler/connectivity_preserve_modeler.h"
#include "utilities/variable_utils.h"

namespace Kratos
{

// Public methods //////////////////////////////////////////////////////////////

void ConnectivityPreserveModeler::GenerateModelPart(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart,
    Element const& rReferenceElement,
    Condition const& rReferenceBoundaryCondition)
{
    KRATOS_TRY;

    this->CheckVariableLists(rOriginModelPart, rDestinationModelPart);

    this->ResetModelPart(rDestinationModelPart);

    this->CopyCommonData(rOriginModelPart, rDestinationModelPart);

    this->DuplicateElements(rOriginModelPart, rDestinationModelPart, rReferenceElement);

    this->DuplicateConditions(rOriginModelPart, rDestinationModelPart, rReferenceBoundaryCondition);

    this->DuplicateCommunicatorData(rOriginModelPart,rDestinationModelPart);

    this->DuplicateSubModelParts(rOriginModelPart, rDestinationModelPart);

    KRATOS_CATCH("");
}

void ConnectivityPreserveModeler::GenerateModelPart(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart,
    Element const& rReferenceElement)
{
    KRATOS_TRY;

    this->CheckVariableLists(rOriginModelPart, rDestinationModelPart);

    this->ResetModelPart(rDestinationModelPart);

    this->CopyCommonData(rOriginModelPart, rDestinationModelPart);

    this->DuplicateElements(rOriginModelPart, rDestinationModelPart, rReferenceElement);

    this->DuplicateCommunicatorData(rOriginModelPart,rDestinationModelPart);

    this->DuplicateSubModelParts(rOriginModelPart, rDestinationModelPart);

    KRATOS_CATCH("");
}

void ConnectivityPreserveModeler::GenerateModelPart(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart,
    Condition const& rReferenceCondition)
{
    KRATOS_TRY;

    this->CheckVariableLists(rOriginModelPart, rDestinationModelPart);

    this->ResetModelPart(rDestinationModelPart);

    this->CopyCommonData(rOriginModelPart, rDestinationModelPart);

    this->DuplicateConditions(rOriginModelPart, rDestinationModelPart, rReferenceCondition);

    this->DuplicateCommunicatorData(rOriginModelPart,rDestinationModelPart);

    this->DuplicateSubModelParts(rOriginModelPart, rDestinationModelPart);

    KRATOS_CATCH("");
}

// Private methods /////////////////////////////////////////////////////////////
void ConnectivityPreserveModeler::CheckVariableLists(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart) const
{
    //check that the variable lists are matching
    const auto& r_destination_variable_list = rDestinationModelPart.GetNodalSolutionStepVariablesList();
    const auto& r_origin_variable_list = rOriginModelPart.GetNodalSolutionStepVariablesList();

    for (const auto& var : r_destination_variable_list) {
        KRATOS_WARNING_IF("VARIABLE LIST MISMATCH - ", !r_origin_variable_list.Has(var))
            << "Variable: " << var << " is in rDestinationModelPart variables "
            << "but not in the rOriginModelPart variables" << std::endl;
    }

    for (const auto& var : r_origin_variable_list) {
        KRATOS_WARNING_IF("VARIABLE LIST MISMATCH - ", !r_destination_variable_list.Has(var))
            << "Variable: " << var << " is in rOriginModelPart variables "
            << "but not in the rDestinationModelPart variables" << std::endl;
    }
}

void ConnectivityPreserveModeler::ResetModelPart(ModelPart& rDestinationModelPart) const
{
    VariableUtils().SetFlag(TO_ERASE, true, rDestinationModelPart.Nodes());
    VariableUtils().SetFlag(TO_ERASE, true, rDestinationModelPart.Elements());
    VariableUtils().SetFlag(TO_ERASE, true, rDestinationModelPart.Conditions());

    rDestinationModelPart.RemoveNodesFromAllLevels(TO_ERASE);
    rDestinationModelPart.RemoveElementsFromAllLevels(TO_ERASE);
    rDestinationModelPart.RemoveConditionsFromAllLevels(TO_ERASE);
}

void ConnectivityPreserveModeler::CopyCommonData(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart) const
{
    // Do not try to change some of the things we need to change if the destination is a SubModelPart
    if( rDestinationModelPart.IsSubModelPart() ) {
        KRATOS_ERROR_IF_NOT(rOriginModelPart.GetNodalSolutionStepVariablesList() == rDestinationModelPart.GetNodalSolutionStepVariablesList())
            << "Attempting to change the SolutionStepVariablesList of the Destination Model Part, which is a SubModelPart." << std::endl
            << "Aborting, since this would break its parent ModelPart." << std::endl;

        KRATOS_ERROR_IF_NOT(rDestinationModelPart.GetBufferSize() == rOriginModelPart.GetBufferSize())
            << "Attempting to change the BufferSize of the Destination Model Part, which is a SubModelPart." << std::endl
            << "Aborting, since this would break its parent ModelPart." << std::endl;
    } else {
        rDestinationModelPart.SetNodalSolutionStepVariablesList(rOriginModelPart.pGetNodalSolutionStepVariablesList());
        rDestinationModelPart.SetBufferSize( rOriginModelPart.GetBufferSize() );
    }

    // These should be safe for SubModelParts
    rDestinationModelPart.SetProcessInfo( rOriginModelPart.pGetProcessInfo() );
    rDestinationModelPart.SetProperties( rOriginModelPart.pProperties() );
    rDestinationModelPart.Tables() = rOriginModelPart.Tables();

    // Assign the nodes to the new model part
    rDestinationModelPart.AddNodes(rOriginModelPart.NodesBegin(), rOriginModelPart.NodesEnd());
}

void ConnectivityPreserveModeler::DuplicateElements(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart,
    const Element& rReferenceElement) const
{
    // Generate the elements
    ModelPart::ElementsContainerType& r_elements = rDestinationModelPart.Elements();
    r_elements.reserve(rOriginModelPart.NumberOfElements());

    // Calculates the maximum element id available in destination model part
    // In here, passing original model part communicator makes sense, since
    // it is only used to MaxAll of a double value (destination model part communicators are still not constructed).
    const IndexType max_element_id =
        GetMaxId(rDestinationModelPart.IsSubModelPart()
                     ? rDestinationModelPart.GetRootModelPart().Elements()
                     : r_elements,
                 rOriginModelPart.GetCommunicator());

    for (auto i_elem = rOriginModelPart.ElementsBegin(); i_elem != rOriginModelPart.ElementsEnd(); ++i_elem) {
        Properties::Pointer properties = i_elem->pGetProperties();

        // Reuse the geometry of the old element (to save memory)
        Element::Pointer p_element = rReferenceElement.Create(i_elem->Id() + max_element_id, i_elem->pGetGeometry(), properties);

        r_elements.push_back(p_element);
    }
}

void ConnectivityPreserveModeler::DuplicateConditions(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart,
    const Condition& rReferenceBoundaryCondition) const
{
    // Generate the conditions
    ModelPart::ConditionsContainerType& r_conditions = rDestinationModelPart.Conditions();
    r_conditions.reserve(rOriginModelPart.NumberOfConditions());

    // Calculates the maximum condition id available in destination model part
    // In here, passing original model part communicator makes sense, since
    // it is only used to MaxAll of a double value (destination model part communicators are still not constructed).
    const IndexType max_condition_id =
        GetMaxId(rDestinationModelPart.IsSubModelPart()
                     ? rDestinationModelPart.GetRootModelPart().Conditions()
                     : r_conditions,
                 rOriginModelPart.GetCommunicator());

    for (auto i_cond = rOriginModelPart.ConditionsBegin(); i_cond != rOriginModelPart.ConditionsEnd(); ++i_cond) {
        Properties::Pointer properties = i_cond->pGetProperties();

        // Reuse the geometry of the old element (to save memory)
        Condition::Pointer p_condition = rReferenceBoundaryCondition.Create(i_cond->Id() + max_condition_id, i_cond->pGetGeometry(), properties);

        r_conditions.push_back(p_condition);
    }
}

void ConnectivityPreserveModeler::DuplicateCommunicatorData(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart) const
{
    /* Create a new communicator for rDestinationModelPart and fill it with the information of the original one
     * Only "general" information and node lists are copied, element and condition lists will be created later
     * using the new elements.
     */
    Communicator& rReferenceComm = rOriginModelPart.GetCommunicator();
    Communicator::Pointer pDestinationComm = rReferenceComm.Create();
    pDestinationComm->SetNumberOfColors( rReferenceComm.GetNumberOfColors() );
    pDestinationComm->NeighbourIndices() = rReferenceComm.NeighbourIndices();

    if (rReferenceComm.IsDistributed()) {
        pDestinationComm->LocalMesh().SetNodes( rReferenceComm.LocalMesh().pNodes() );
        pDestinationComm->InterfaceMesh().SetNodes( rReferenceComm.InterfaceMesh().pNodes() );
        pDestinationComm->GhostMesh().SetNodes( rReferenceComm.GhostMesh().pNodes() );
        for (unsigned int i = 0; i < rReferenceComm.GetNumberOfColors(); i++) {
            pDestinationComm->pLocalMesh(i)->SetNodes( rReferenceComm.pLocalMesh(i)->pNodes() );
            pDestinationComm->pInterfaceMesh(i)->SetNodes( rReferenceComm.pInterfaceMesh(i)->pNodes() );
            pDestinationComm->pGhostMesh(i)->SetNodes( rReferenceComm.pGhostMesh(i)->pNodes() );
        }

        // All elements are passed as local elements to the new communicator
        ModelPart::ElementsContainerType& rDestinationLocalElements = pDestinationComm->LocalMesh().Elements();
        rDestinationLocalElements.clear();
        rDestinationLocalElements.reserve(rDestinationModelPart.NumberOfElements());
        for (auto i_elem = rDestinationModelPart.Elements().ptr_begin(); i_elem != rDestinationModelPart.Elements().ptr_end(); ++i_elem) {
            rDestinationLocalElements.push_back(*i_elem);
        }

        // Do the same for Conditions
        ModelPart::ConditionsContainerType& rDestinationLocalConditions = pDestinationComm->LocalMesh().Conditions();
        rDestinationLocalConditions.clear();
        rDestinationLocalConditions.reserve(rDestinationModelPart.NumberOfConditions());
        for (auto i_cond = rDestinationModelPart.Conditions().ptr_begin(); i_cond != rDestinationModelPart.Conditions().ptr_end(); ++i_cond) {
            rDestinationLocalConditions.push_back(*i_cond);
        }
    } else {
        pDestinationComm->SetLocalMesh(rDestinationModelPart.pGetMesh());
    }

    rDestinationModelPart.SetCommunicator( pDestinationComm );
}

void ConnectivityPreserveModeler::DuplicateSubModelParts(
ModelPart& rOriginModelPart,
ModelPart& rDestinationModelPart) const
{
    // If we copy root model part into a submodel part, then it will create element/condition whith ids which may
    // conflict with destination model part root model part ids
    // or we need to add all the submodel part elements and conditions in the origin model part as well
    if (!rDestinationModelPart.IsSubModelPart())
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
                for(auto it=i_part->ElementsBegin(); it!=i_part->ElementsEnd(); ++it)
                    ids.push_back(it->Id());
                destination_part.AddElements(ids, 0); //adding by index
            }

            // Execute only if we created conditions in the destination
            if (rDestinationModelPart.NumberOfConditions() > 0)
            {
                ids.clear();
                for(auto it=i_part->ConditionsBegin(); it!=i_part->ConditionsEnd(); ++it)
                    ids.push_back(it->Id());
                destination_part.AddConditions(ids, 0);
            }

            // Duplicate the Communicator for this SubModelPart
            this->DuplicateCommunicatorData(*i_part, destination_part);

            // Recursively call this function to duplicate any child SubModelParts
            this->DuplicateSubModelParts(*i_part, destination_part);
        }
    }
}

}
