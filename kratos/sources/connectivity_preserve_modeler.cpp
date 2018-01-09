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

namespace Kratos
{

// Public methods //////////////////////////////////////////////////////////////

ConnectivityPreserveModeler::ConnectivityPreserveModeler():
    Modeler()
{
}


ConnectivityPreserveModeler::~ConnectivityPreserveModeler()
{
}


void ConnectivityPreserveModeler::GenerateModelPart(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart,
    Element const& rReferenceElement,
    Condition const& rReferenceBoundaryCondition)
{
    KRATOS_TRY;

    this->ResetModelPart(rDestinationModelPart);

    this->CopyCommonData(rOriginModelPart, rDestinationModelPart);

    this->DuplicateElements(rOriginModelPart, rDestinationModelPart, rReferenceElement);

    this->DuplicateConditions(rOriginModelPart, rDestinationModelPart, rReferenceBoundaryCondition);

    this->DuplicateCommunicatorData(rOriginModelPart,rDestinationModelPart);

    this->DuplicateSubModelParts(rOriginModelPart, rDestinationModelPart);

    KRATOS_CATCH("");
}

// Private methods /////////////////////////////////////////////////////////////

void ConnectivityPreserveModeler::ResetModelPart(ModelPart &rDestinationModelPart)
{
    for(auto it = rDestinationModelPart.NodesBegin(); it != rDestinationModelPart.NodesEnd(); it++)
        it->Set(TO_ERASE);
    rDestinationModelPart.RemoveNodes(TO_ERASE);

    for(auto it = rDestinationModelPart.ElementsBegin(); it != rDestinationModelPart.ElementsEnd(); it++)
        it->Set(TO_ERASE);
    rDestinationModelPart.RemoveElements(TO_ERASE);

    for(auto it = rDestinationModelPart.ElementsBegin(); it != rDestinationModelPart.ElementsEnd(); it++)
        it->Set(TO_ERASE);
    rDestinationModelPart.RemoveConditions(TO_ERASE);
}


void ConnectivityPreserveModeler::CopyCommonData(
    ModelPart &rOriginModelPart,
    ModelPart &rDestinationModelPart)
{
    // Do not try to change some of the things we need to change if the destination is a SubModelPart
    if( rDestinationModelPart.IsSubModelPart() )
    {
        if( !(rOriginModelPart.GetNodalSolutionStepVariablesList() == rDestinationModelPart.GetNodalSolutionStepVariablesList()) )
        {
            KRATOS_ERROR << "Attempting to change the SolutionStepVariablesList of the Destination Model Part, which is a SubModelPart." << std::endl
                         << "Aborting, since this would break its parent ModelPart." << std::endl;
        }

        if( rDestinationModelPart.GetBufferSize() != rOriginModelPart.GetBufferSize() )
        {
            KRATOS_ERROR << "Attempting to change the BufferSize of the Destination Model Part, which is a SubModelPart." << std::endl
                         << "Aborting, since this would break its parent ModelPart." << std::endl;
        }
    }
    else
    {
        rDestinationModelPart.GetNodalSolutionStepVariablesList() = rOriginModelPart.GetNodalSolutionStepVariablesList();
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
    ModelPart &rOriginModelPart,
    ModelPart &rDestinationModelPart,
    Element const &rReferenceElement)
{
    // Generate the elements
    ModelPart::ElementsContainerType temp_elements;
    temp_elements.reserve(rOriginModelPart.NumberOfElements());
    for (auto i_elem = rOriginModelPart.ElementsBegin(); i_elem != rOriginModelPart.ElementsEnd(); ++i_elem)
    {
        Properties::Pointer properties = i_elem->pGetProperties();
        Element::Pointer p_element = rReferenceElement.Create(i_elem->Id(), i_elem->GetGeometry(), properties);

        // Reuse the geometry of the old element (to save memory)
        p_element->pGetGeometry() = i_elem->pGetGeometry();

        temp_elements.push_back(p_element);
    }

    rDestinationModelPart.AddElements(temp_elements.begin(), temp_elements.end());
}


void ConnectivityPreserveModeler::DuplicateConditions(
        ModelPart &rOriginModelPart,
        ModelPart &rDestinationModelPart,
        Condition const &rReferenceBoundaryCondition)
{
    // Generate the conditions
    ModelPart::ConditionsContainerType temp_conditions;
    temp_conditions.reserve(rOriginModelPart.NumberOfConditions());
    for (auto i_cond = rOriginModelPart.ConditionsBegin(); i_cond != rOriginModelPart.ConditionsEnd(); ++i_cond)
    {
        Properties::Pointer properties = i_cond->pGetProperties();
        Condition::Pointer p_condition = rReferenceBoundaryCondition.Create(i_cond->Id(), i_cond->GetGeometry(), properties);

        // Reuse the geometry of the old element (to save memory)
        p_condition->pGetGeometry() = i_cond->pGetGeometry();

        temp_conditions.push_back(p_condition);
    }

    rDestinationModelPart.AddConditions(temp_conditions.begin(), temp_conditions.end());
}

void ConnectivityPreserveModeler::DuplicateCommunicatorData(
    ModelPart &rOriginModelPart,
    ModelPart &rDestinationModelPart)
{
    /* Create a new communicator for rDestinationModelPart and fill it with the information of the original one
     * Only "general" information and node lists are copied, element and condition lists will be created later
     * using the new elements.
     */
    Communicator& rReferenceComm = rOriginModelPart.GetCommunicator();
    Communicator::Pointer pDestinationComm = rReferenceComm.Create();
    pDestinationComm->SetNumberOfColors( rReferenceComm.GetNumberOfColors() );
    pDestinationComm->NeighbourIndices() = rReferenceComm.NeighbourIndices();
    pDestinationComm->LocalMesh().SetNodes( rReferenceComm.LocalMesh().pNodes() );
    pDestinationComm->InterfaceMesh().SetNodes( rReferenceComm.InterfaceMesh().pNodes() );
    pDestinationComm->GhostMesh().SetNodes( rReferenceComm.GhostMesh().pNodes() );
    for (unsigned int i = 0; i < rReferenceComm.GetNumberOfColors(); i++)
    {
        pDestinationComm->pLocalMesh(i)->SetNodes( rReferenceComm.pLocalMesh(i)->pNodes() );
        pDestinationComm->pInterfaceMesh(i)->SetNodes( rReferenceComm.pInterfaceMesh(i)->pNodes() );
        pDestinationComm->pGhostMesh(i)->SetNodes( rReferenceComm.pGhostMesh(i)->pNodes() );
    }

    // This is a dirty hack to detect if the communicator is a Communicator or an MPICommunicator
    // Note that downcasting would not work here because MPICommunicator is not compiled in non-MPI builds
    bool is_mpi = ( rOriginModelPart.pElements() == rReferenceComm.LocalMesh().pElements() );

    if (is_mpi)
    {
        // All elements are passed as local elements to the new communicator
        ModelPart::ElementsContainerType& rDestinationLocalElements = pDestinationComm->LocalMesh().Elements();
        rDestinationLocalElements.clear();
        rDestinationLocalElements.reserve(rDestinationModelPart.NumberOfElements());
        for (auto i_elem = rDestinationModelPart.Elements().ptr_begin(); i_elem != rDestinationModelPart.Elements().ptr_end(); ++i_elem)
        {
            rDestinationLocalElements.push_back(*i_elem);
        }

        // Do the same for Conditions
        ModelPart::ConditionsContainerType& rDestinationLocalConditions = pDestinationComm->LocalMesh().Conditions();
        rDestinationLocalConditions.clear();
        rDestinationLocalConditions.reserve(rDestinationModelPart.NumberOfConditions());
        for (auto i_cond = rDestinationModelPart.Conditions().ptr_begin(); i_cond != rDestinationModelPart.Conditions().ptr_end(); ++i_cond)
        {
            rDestinationLocalConditions.push_back(*i_cond);
        }
    }
    else
    {
        pDestinationComm->LocalMesh().pElements() = rDestinationModelPart.pElements();
        pDestinationComm->LocalMesh().pConditions() = rDestinationModelPart.pConditions();
    }

    rDestinationModelPart.SetCommunicator( pDestinationComm );
}


void ConnectivityPreserveModeler::DuplicateSubModelParts(
    ModelPart &rOriginModelPart,
    ModelPart &rDestinationModelPart)
{
    for(auto i_part = rOriginModelPart.SubModelPartsBegin(); i_part != rOriginModelPart.SubModelPartsEnd(); ++i_part)
    {
        ModelPart& destination_part = *(rDestinationModelPart.CreateSubModelPart(i_part->Name()));

        destination_part.AddNodes(i_part->NodesBegin(), i_part->NodesEnd());

        destination_part.AddElements(i_part->ElementsBegin(), i_part->ElementsEnd());

        destination_part.AddConditions(i_part->ConditionsBegin(), i_part->ConditionsEnd());

        // Duplicate the Communicator for this SubModelPart
        this->DuplicateCommunicatorData(*i_part, destination_part);

        // Recursively call this function to duplicate any child SubModelParts
        this->DuplicateSubModelParts(*i_part, destination_part);
    }
}

}
