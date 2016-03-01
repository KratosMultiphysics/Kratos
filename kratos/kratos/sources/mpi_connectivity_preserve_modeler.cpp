//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ \.  
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#include "includes/node.h"
#include "includes/properties.h"
#include "includes/communicator.h"


#include "modeler/mpi_connectivity_preserve_modeler.h"

namespace Kratos
{

MPIConnectivityPreserveModeler::MPIConnectivityPreserveModeler():
    Modeler()
{
}

MPIConnectivityPreserveModeler::~MPIConnectivityPreserveModeler()
{

}

void MPIConnectivityPreserveModeler::GenerateModelPart(ModelPart &rOriginModelPart,
                                                       ModelPart &rDestinationModelPart,
                                                       const Element &rReferenceElement,
                                                       const Condition &rReferenceCondition)
{
    // Copy general ModelPart properites
    rDestinationModelPart.GetNodalSolutionStepVariablesList() = rOriginModelPart.GetNodalSolutionStepVariablesList();
    rDestinationModelPart.SetBufferSize( rOriginModelPart.GetBufferSize() );
    rDestinationModelPart.SetProcessInfo( rOriginModelPart.pGetProcessInfo() );
    rDestinationModelPart.SetProperties( rOriginModelPart.pProperties() );

    // Copy tables
    rDestinationModelPart.Tables() = rOriginModelPart.Tables();

    // Copy the node list so that both model parts share the same nodes
    rDestinationModelPart.SetNodes( rOriginModelPart.pNodes() );

    /* Create a new communicator for rDestinationModelPart and fill it with the information of the original one
     * Only "general" information and node lists are passed, element and condition lists will be created later
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

    rDestinationModelPart.SetCommunicator( pDestinationComm );

    // Reset element container and create new elements
    rDestinationModelPart.Elements().clear();
    rDestinationModelPart.Elements().reserve(rOriginModelPart.NumberOfElements());

    for (ModelPart::ElementsContainerType::iterator iEl = rOriginModelPart.ElementsBegin();
         iEl != rOriginModelPart.ElementsEnd(); iEl++)
    {
        Properties::Pointer pProp = iEl->pGetProperties();
        Element::Pointer pElem = rReferenceElement.Create(iEl->Id(),iEl->GetGeometry(),pProp);
        rDestinationModelPart.Elements().push_back(pElem);
    }

    // All elements are passed as local elements to the new communicator
    ModelPart::ElementsContainerType& rDestinationLocalElements = pDestinationComm->LocalMesh().Elements();
    rDestinationLocalElements.clear();
    rDestinationLocalElements.reserve(rDestinationModelPart.NumberOfElements());
    for (ModelPart::ElementsContainerType::ptr_iterator iEl = rDestinationModelPart.Elements().ptr_begin();
         iEl != rDestinationModelPart.Elements().ptr_end(); iEl++)
    {
        rDestinationLocalElements.push_back(*iEl);
    }

    // Reset condition container and create new conditions
    rDestinationModelPart.Conditions().clear();
    rDestinationModelPart.Conditions().reserve(rOriginModelPart.NumberOfConditions());
    for (ModelPart::ConditionsContainerType::iterator iCo = rOriginModelPart.ConditionsBegin();
         iCo != rOriginModelPart.ConditionsEnd(); iCo++)
    {
        Properties::Pointer pProp = iCo->pGetProperties();
        Condition::Pointer pCond = rReferenceCondition.Create(iCo->Id(),iCo->GetGeometry(),pProp);
        rDestinationModelPart.Conditions().push_back(pCond);
    }

    // Create new communicator local condition list
    ModelPart::ConditionsContainerType& rDestinationLocalConditions = pDestinationComm->LocalMesh().Conditions();
    rDestinationLocalConditions.clear();
    rDestinationLocalConditions.reserve(rDestinationModelPart.NumberOfConditions());
    for (ModelPart::ConditionsContainerType::ptr_iterator iCo = rDestinationModelPart.Conditions().ptr_begin();
         iCo != rDestinationModelPart.Conditions().ptr_end(); iCo++)
    {
        rDestinationLocalConditions.push_back(*iCo);
    }
}

}
