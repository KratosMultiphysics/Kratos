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

// System includes

// External includes

// Project includes
#include "modeler/connectivity_preserve_modeler.h"
#include "utilities/variable_utils.h"

namespace Kratos
{

// Public methods //////////////////////////////////////////////////////////////

ConnectivityPreserveModeler::ConnectivityPreserveModeler(Model& rModel, Parameters Settings):
    Modeler(Settings),
    mpModel(&rModel)
{
    mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}


Modeler::Pointer ConnectivityPreserveModeler::Create(Model& rModel, const Parameters Settings) const
{
    return Kratos::make_shared<ConnectivityPreserveModeler>(rModel, Settings);
}


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

void ConnectivityPreserveModeler::GenerateModelPart(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart)
{
    KRATOS_TRY;

    this->CheckVariableLists(rOriginModelPart, rDestinationModelPart);

    this->ResetModelPart(rDestinationModelPart);

    this->CopyCommonData(rOriginModelPart, rDestinationModelPart);

    this->DuplicateEntities(rOriginModelPart, rDestinationModelPart);

    this->DuplicateCommunicatorData(rOriginModelPart,rDestinationModelPart);

    this->DuplicateSubModelParts(rOriginModelPart, rDestinationModelPart);

    KRATOS_CATCH("");
}


void ConnectivityPreserveModeler::SetupModelPart()
{
    ModelPart& r_origin_model_part = mpModel->GetModelPart(mParameters["origin_model_part_name"].GetString());

    const std::string destination_model_part_name = mParameters["destination_model_part_name"].GetString();
    if (!mpModel->HasModelPart(destination_model_part_name)) {
        mpModel->CreateModelPart(destination_model_part_name, r_origin_model_part.GetBufferSize());
    }
    ModelPart& r_destination_model_part = mpModel->GetModelPart(destination_model_part_name);

    const std::string element_name = mParameters["reference_element"].GetString();
    const std::string condition_name = mParameters["reference_condition"].GetString();

    if (element_name != "") {
        if (condition_name != "") {
            GenerateModelPart(
                r_origin_model_part,
                r_destination_model_part,
                KratosComponents<Element>::Get(element_name),
                KratosComponents<Condition>::Get(condition_name));
        } else {
            GenerateModelPart(
                r_origin_model_part,
                r_destination_model_part,
                KratosComponents<Element>::Get(element_name));
        }
    } else if (condition_name != "") {
        GenerateModelPart(
            r_origin_model_part,
            r_destination_model_part,
            KratosComponents<Condition>::Get(condition_name));
    } else {
        KRATOS_INFO("ConnectivityPreserveModeler") << "Generating destination model part with default Kratos elements and conditions." << std::endl;
        GenerateModelPart(r_origin_model_part, r_destination_model_part);
    }
}


const Parameters ConnectivityPreserveModeler::GetDefaultParameters() const
{
    return Parameters(R"({
        "origin_model_part_name"        : "undefined_origin_model_part_name",
        "destination_model_part_name"   : "undefined_destination_model_part_name",
        "reference_element"             : "",
        "reference_condition"           : ""
    })");
}


std::string ConnectivityPreserveModeler::Info() const
{
    return "ConnectivityPreserveModeler";
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
    rDestinationModelPart.PropertiesArray() = rOriginModelPart.PropertiesArray();
    rDestinationModelPart.Tables() = rOriginModelPart.Tables();

    // Assign the nodes to the new model part
    rDestinationModelPart.AddNodes(rOriginModelPart.NodesBegin(), rOriginModelPart.NodesEnd());

    // Assign the geometries to the new model part
    rDestinationModelPart.Geometries() = rOriginModelPart.Geometries();
}

void ConnectivityPreserveModeler::DuplicateElements(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart,
    const Element& rReferenceElement) const
{
    // Generate the elements
    ModelPart::ElementsContainerType temp_elements;
    temp_elements.reserve(rOriginModelPart.NumberOfElements());
    for (auto i_elem = rOriginModelPart.ElementsBegin(); i_elem != rOriginModelPart.ElementsEnd(); ++i_elem) {
        Properties::Pointer properties = i_elem->pGetProperties();

        // Reuse the geometry of the old element (to save memory)
        Element::Pointer p_element = rReferenceElement.Create(i_elem->Id(), i_elem->pGetGeometry(), properties);

        temp_elements.push_back(p_element);
    }

    rDestinationModelPart.AddElements(temp_elements.begin(), temp_elements.end());
}

void ConnectivityPreserveModeler::DuplicateConditions(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart,
    const Condition& rReferenceBoundaryCondition) const
{
    // Generate the conditions
    ModelPart::ConditionsContainerType temp_conditions;
    temp_conditions.reserve(rOriginModelPart.NumberOfConditions());
    for (auto i_cond = rOriginModelPart.ConditionsBegin(); i_cond != rOriginModelPart.ConditionsEnd(); ++i_cond) {
        Properties::Pointer properties = i_cond->pGetProperties();

        // Reuse the geometry of the old element (to save memory)
        Condition::Pointer p_condition = rReferenceBoundaryCondition.Create(i_cond->Id(), i_cond->pGetGeometry(), properties);

        temp_conditions.push_back(p_condition);
    }

    rDestinationModelPart.AddConditions(temp_conditions.begin(), temp_conditions.end());
}

void ConnectivityPreserveModeler::DuplicateEntities(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart) const
{
    // Set the geometries to entities maps
    const auto& r_elems_map = SetGeometryElementMap();
    const auto& r_conds_map = SetGeometryConditionMap();

    // Generate the elements
    ModelPart::ElementsContainerType temp_elements;
    temp_elements.reserve(rOriginModelPart.NumberOfElements());
    temp_elements.reserve(rOriginModelPart.NumberOfElements());
    for (auto it_elem = rOriginModelPart.ElementsBegin(); it_elem != rOriginModelPart.ElementsEnd(); ++it_elem) {
        // Create the new element reusing the geometry of the old one (to save memory)
        const auto& r_ref_elem = r_elems_map.at(it_elem->GetGeometry().GetGeometryType());
        auto p_element = r_ref_elem.Create(it_elem->Id(), it_elem->pGetGeometry(), it_elem->pGetProperties());
        // Add the new element to the temporary container
        temp_elements.push_back(p_element);
    }
    rDestinationModelPart.AddElements(temp_elements.begin(), temp_elements.end());

    // Generate the conditions
    ModelPart::ConditionsContainerType temp_conditions;
    temp_conditions.reserve(rOriginModelPart.NumberOfConditions());
    for (auto it_cond = rOriginModelPart.ConditionsBegin(); it_cond != rOriginModelPart.ConditionsEnd(); ++it_cond) {
        // Create the new condition reusing the geometry of the old one (to save memory)
        const auto& r_ref_cond = r_conds_map.at(it_cond->GetGeometry().GetGeometryType());
        Condition::Pointer p_condition = r_ref_cond.Create(it_cond->Id(), it_cond->pGetGeometry(), it_cond->pGetProperties());
        // Add the new condition to the temporary container
        temp_conditions.push_back(p_condition);
    }
    rDestinationModelPart.AddConditions(temp_conditions.begin(), temp_conditions.end());
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

ConnectivityPreserveModeler::GeometryElementMapType ConnectivityPreserveModeler::SetGeometryElementMap() const
{
    GeometryElementMapType geom_elem_map = {
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4, KratosComponents<Element>::Get("Element2D4N")},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8, KratosComponents<Element>::Get("Element2D8N")},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9, KratosComponents<Element>::Get("Element2D9N")},
        {GeometryData::KratosGeometryType::Kratos_Triangle2D3, KratosComponents<Element>::Get("Element2D3N")},
        {GeometryData::KratosGeometryType::Kratos_Triangle2D6, KratosComponents<Element>::Get("Element2D6N")},
        {GeometryData::KratosGeometryType::Kratos_Hexahedra3D8, KratosComponents<Element>::Get("Element3D8N")},
        {GeometryData::KratosGeometryType::Kratos_Hexahedra3D20, KratosComponents<Element>::Get("Element3D20N")},
        {GeometryData::KratosGeometryType::Kratos_Hexahedra3D20, KratosComponents<Element>::Get("Element3D27N")},
        {GeometryData::KratosGeometryType::Kratos_Prism3D6, KratosComponents<Element>::Get("Element3D6N")},
        {GeometryData::KratosGeometryType::Kratos_Prism3D15, KratosComponents<Element>::Get("Element3D15N")},
        {GeometryData::KratosGeometryType::Kratos_Pyramid3D5, KratosComponents<Element>::Get("Element3D5N")},
        {GeometryData::KratosGeometryType::Kratos_Pyramid3D13, KratosComponents<Element>::Get("Element3D13N")},
        {GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10, KratosComponents<Element>::Get("Element3D10N")},
        {GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4, KratosComponents<Element>::Get("Element3D4N")},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4, KratosComponents<Element>::Get("SurfaceElement3D4N")},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8, KratosComponents<Element>::Get("SurfaceElement3D8N")},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9, KratosComponents<Element>::Get("SurfaceElement3D9N")},
        {GeometryData::KratosGeometryType::Kratos_Triangle3D3, KratosComponents<Element>::Get("SurfaceElement3D3N")},
        {GeometryData::KratosGeometryType::Kratos_Triangle3D6, KratosComponents<Element>::Get("SurfaceElement3D6N")},
        {GeometryData::KratosGeometryType::Kratos_Line2D2, KratosComponents<Element>::Get("LineElement2D2N")},
        {GeometryData::KratosGeometryType::Kratos_Line2D3, KratosComponents<Element>::Get("LineElement2D3N")},
        {GeometryData::KratosGeometryType::Kratos_Line3D2, KratosComponents<Element>::Get("LineElement3D2N")},
        {GeometryData::KratosGeometryType::Kratos_Line3D3, KratosComponents<Element>::Get("LineElement3D3N")},
        {GeometryData::KratosGeometryType::Kratos_Point2D, KratosComponents<Element>::Get("PointElement2D1N")},
        {GeometryData::KratosGeometryType::Kratos_Point3D, KratosComponents<Element>::Get("PointElement3D1N")},
    };

        return geom_elem_map;
}

ConnectivityPreserveModeler::GeometryConditionMapType ConnectivityPreserveModeler::SetGeometryConditionMap() const
{
    GeometryConditionMapType geom_cond_map = {
        {GeometryData::KratosGeometryType::Kratos_Line2D2, KratosComponents<Condition>::Get("LineCondition2D2N")},
        {GeometryData::KratosGeometryType::Kratos_Line2D3, KratosComponents<Condition>::Get("LineCondition2D3N")},
        {GeometryData::KratosGeometryType::Kratos_Line3D2, KratosComponents<Condition>::Get("LineCondition3D2N")},
        {GeometryData::KratosGeometryType::Kratos_Line3D3, KratosComponents<Condition>::Get("LineCondition3D3N")},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4, KratosComponents<Condition>::Get("SurfaceCondition3D4N")},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8, KratosComponents<Condition>::Get("SurfaceCondition3D8N")},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9, KratosComponents<Condition>::Get("SurfaceCondition3D9N")},
        {GeometryData::KratosGeometryType::Kratos_Triangle3D3, KratosComponents<Condition>::Get("SurfaceCondition3D3N")},
        {GeometryData::KratosGeometryType::Kratos_Triangle3D6, KratosComponents<Condition>::Get("SurfaceCondition3D6N")},
        {GeometryData::KratosGeometryType::Kratos_Point2D, KratosComponents<Condition>::Get("PointCondition2D1N")},
        {GeometryData::KratosGeometryType::Kratos_Point3D, KratosComponents<Condition>::Get("PointCondition3D1N")},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4, KratosComponents<Condition>::Get("PrismCondition2D4N")},
        {GeometryData::KratosGeometryType::Kratos_Prism3D6, KratosComponents<Condition>::Get("PrismCondition3D6N")},
    };

    return geom_cond_map;
}

}
