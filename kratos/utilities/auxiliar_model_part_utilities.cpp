//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "includes/key_hash.h"
#include "utilities/auxiliar_model_part_utilities.h"
#include "utilities/parallel_utilities.h"
#include "variable_utils.h"
#include "containers/model.h"

namespace Kratos
{
void AuxiliarModelPartUtilities::CopySubModelPartStructure(const ModelPart& rModelPartToCopyFromIt, ModelPart& rModelPartToCopyIntoIt)
{
    for (auto& r_old_sub_model_part : rModelPartToCopyFromIt.SubModelParts()) {
        auto& r_new_sub_model_part = rModelPartToCopyIntoIt.CreateSubModelPart(r_old_sub_model_part.Name());
        if (r_old_sub_model_part.NumberOfSubModelParts() > 0) {
            CopySubModelPartStructure(r_old_sub_model_part, r_new_sub_model_part);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RecursiveEnsureModelPartOwnsProperties(const bool RemovePreviousProperties)
{
    // First we do in this model part
    EnsureModelPartOwnsProperties(RemovePreviousProperties);

    // Now we do in submodelparts
    for (auto& r_old_sub_model_part : mrModelPart.SubModelParts()) {
        AuxiliarModelPartUtilities(r_old_sub_model_part).RecursiveEnsureModelPartOwnsProperties(RemovePreviousProperties);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::EnsureModelPartOwnsProperties(const bool RemovePreviousProperties)
{
    // First we clear the properties if we want so
    if (RemovePreviousProperties) {
        mrModelPart.GetMesh(0).pProperties()->clear();
    }

    // The list of properties
    std::unordered_set<Properties::Pointer, IndexedObjectPointerHasher<Properties::Pointer>, IndexedObjectPointerComparator<Properties::Pointer>> list_of_properties;

    // Iterating over the elements
    auto& r_elements_array = mrModelPart.Elements();
    const auto it_elem_begin= r_elements_array.begin();
    const int number_of_elements = static_cast<int>(r_elements_array.size());

    // Iterating over the conditions
    auto& r_conditions_array = mrModelPart.Conditions();
    const auto it_cond_begin= r_conditions_array.begin();
    const int number_of_conditions = static_cast<int>(r_conditions_array.size());

    #pragma omp parallel
    {
        // The list of properties
        std::unordered_set<Properties::Pointer, IndexedObjectPointerHasher<Properties::Pointer>, IndexedObjectPointerComparator<Properties::Pointer>> buffer_list_of_properties;

        #pragma omp for schedule(dynamic, 512) nowait
        for (int i = 0; i < number_of_elements; ++i) {
            auto it_elem = it_elem_begin + i;

            Properties::Pointer p_prop = it_elem->pGetProperties();

            if (buffer_list_of_properties.find(p_prop) == buffer_list_of_properties.end()) {
                buffer_list_of_properties.insert(p_prop);
            }
        }

        #pragma omp for schedule(dynamic, 512) nowait
        for (int i = 0; i < number_of_conditions; ++i) {
            auto it_cond = it_cond_begin + i;

            Properties::Pointer p_prop = it_cond->pGetProperties();
            if (buffer_list_of_properties.find(p_prop) == buffer_list_of_properties.end()) {
                buffer_list_of_properties.insert(p_prop);
            }
        }

        // Combine buffers together
        #pragma omp critical
        {
            list_of_properties.insert(buffer_list_of_properties.begin(),buffer_list_of_properties.end());
        }
    }

    // Add properties to respective model parts
    for (const auto& p_prop : list_of_properties) {
        if (!mrModelPart.HasProperties(p_prop->Id())) {
            mrModelPart.AddProperties(p_prop);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementAndBelongings(
    IndexType ElementId, Flags IdentifierFlag, IndexType ThisIndex)
{
    auto& r_array_nodes = mrModelPart.Nodes(ThisIndex);

    VariableUtils().SetFlag(IdentifierFlag, true, r_array_nodes);

    block_for_each(
        mrModelPart.Elements(ThisIndex), [&IdentifierFlag,ElementId]( ModelPart::ElementType& rElement )
        {
            if (rElement.Id() != ElementId)
                for (auto& r_node : rElement.GetGeometry())
                    r_node.Set(IdentifierFlag, false);
        }
    );

    bool condition_to_remove;
    for (auto& cond : mrModelPart.Conditions(ThisIndex)) {
        condition_to_remove = true;
        for (auto& node : cond.GetGeometry()) {
            if (node.IsNot(IdentifierFlag)) {
                condition_to_remove = false;
                break;
            }
            if (condition_to_remove) cond.Set(IdentifierFlag);
        }
    }

    mrModelPart.RemoveElement(ElementId, ThisIndex);
    mrModelPart.RemoveConditions(IdentifierFlag);
    mrModelPart.RemoveNodes(IdentifierFlag);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementAndBelongings(Element& ThisElement, const Flags IdentifierFlag , IndexType ThisIndex)
{
    RemoveElementAndBelongings(ThisElement.Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementAndBelongings(Element::Pointer pThisElement, const Flags IdentifierFlag, IndexType ThisIndex)
{
    RemoveElementAndBelongings(pThisElement->Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementAndBelongingsFromAllLevels(IndexType ElementId, const Flags IdentifierFlag, IndexType ThisIndex)
{
    if (mrModelPart.IsSubModelPart()) {
        AuxiliarModelPartUtilities aux_utility = AuxiliarModelPartUtilities(mrModelPart.GetParentModelPart());
        aux_utility.RemoveElementAndBelongings(ElementId, IdentifierFlag, ThisIndex);
    } else {
        RemoveElementAndBelongings(ElementId, IdentifierFlag, ThisIndex);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementAndBelongingsFromAllLevels(Element& ThisElement, const Flags IdentifierFlag, IndexType ThisIndex)
{
    RemoveElementAndBelongingsFromAllLevels(ThisElement.Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementAndBelongingsFromAllLevels(Element::Pointer pThisElement, const Flags IdentifierFlag, IndexType ThisIndex)
{
    RemoveElementAndBelongingsFromAllLevels(pThisElement->Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementsAndBelongings(Flags IdentifierFlag)
{
    //loop over all the meshes
    VariableUtils variable_utils;
    auto& meshes = mrModelPart.GetMeshes();
    for(auto i_mesh = meshes.begin() ; i_mesh != meshes.end() ; i_mesh++) {
        auto& r_array_nodes = i_mesh->Nodes();
        variable_utils.SetFlag(IdentifierFlag, true, r_array_nodes);

        block_for_each(
            i_mesh->Elements(),
            [&IdentifierFlag](Element& rElement)
            {
                if (rElement.IsNot(IdentifierFlag))
                    for (auto& r_node : rElement.GetGeometry())
                        r_node.Set(IdentifierFlag, false);
            }
        );

        bool condition_to_remove;
        for (auto& cond : i_mesh->Conditions()) {
            condition_to_remove = true;
            for (auto& node : cond.GetGeometry()) {
                if (node.IsNot(IdentifierFlag)) {
                    condition_to_remove = false;
                    break;
                }
                if (condition_to_remove) cond.Set(IdentifierFlag);
            }
        }
    }

    mrModelPart.RemoveElements(IdentifierFlag);
    mrModelPart.RemoveConditions(IdentifierFlag);
    mrModelPart.RemoveNodes(IdentifierFlag);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveElementsAndBelongingsFromAllLevels(const Flags IdentifierFlag)
{
    ModelPart& root_model_part = mrModelPart.GetRootModelPart();
    AuxiliarModelPartUtilities aux_utility = AuxiliarModelPartUtilities(root_model_part);
    aux_utility.RemoveElementsAndBelongings(IdentifierFlag);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionAndBelongings(IndexType ConditionId, Flags IdentifierFlag, IndexType ThisIndex)
{
    auto& r_array_nodes = mrModelPart.Nodes(ThisIndex);
    VariableUtils().SetFlag(IdentifierFlag, true, r_array_nodes);

    block_for_each(
        mrModelPart.Conditions(ThisIndex),
        [&IdentifierFlag,ConditionId](ModelPart::ConditionType& rCondition)
        {
            if (rCondition.Id() != ConditionId)
                for (auto& r_node : rCondition.GetGeometry())
                    r_node.Set(IdentifierFlag, false);
        }
    );
    bool element_to_remove;
    for (auto& elem : mrModelPart.Elements(ThisIndex)) {
        element_to_remove = true;
        for (auto& node : elem.GetGeometry()) {
            if (node.IsNot(IdentifierFlag)) {
                element_to_remove = false;
                break;
            }
            if (element_to_remove) elem.Set(IdentifierFlag);
        }
    }

    mrModelPart.RemoveCondition(ConditionId, ThisIndex);
    mrModelPart.RemoveElements(IdentifierFlag);
    mrModelPart.RemoveNodes(IdentifierFlag);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionAndBelongings(Condition& ThisCondition, const Flags IdentifierFlag , IndexType ThisIndex)
{
    RemoveConditionAndBelongings(ThisCondition.Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionAndBelongings(Condition::Pointer pThisCondition, const Flags IdentifierFlag, IndexType ThisIndex)
{
    RemoveConditionAndBelongings(pThisCondition->Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionAndBelongingsFromAllLevels(IndexType ConditionId, const Flags IdentifierFlag, IndexType ThisIndex)
{
    if (mrModelPart.IsSubModelPart()) {
        AuxiliarModelPartUtilities aux_utility = AuxiliarModelPartUtilities(mrModelPart.GetParentModelPart());
        aux_utility.RemoveConditionAndBelongings(ConditionId, IdentifierFlag, ThisIndex);
    } else {
        RemoveConditionAndBelongings(ConditionId, IdentifierFlag, ThisIndex);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionAndBelongingsFromAllLevels(Condition& ThisCondition, const Flags IdentifierFlag, IndexType ThisIndex)
{
    RemoveConditionAndBelongingsFromAllLevels(ThisCondition.Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionAndBelongingsFromAllLevels(Condition::Pointer pThisCondition, const Flags IdentifierFlag, IndexType ThisIndex)
{
    RemoveConditionAndBelongingsFromAllLevels(pThisCondition->Id(), IdentifierFlag, ThisIndex);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionsAndBelongings(Flags IdentifierFlag)
{
    //loop over all the meshes
    VariableUtils variable_utils;
    auto& meshes = mrModelPart.GetMeshes();
    for(auto i_mesh = meshes.begin() ; i_mesh != meshes.end() ; i_mesh++) {
        auto& r_array_nodes = i_mesh->Nodes();
        variable_utils.SetFlag(IdentifierFlag, true, r_array_nodes);

        block_for_each(
            i_mesh->Conditions(),
            [&IdentifierFlag](ModelPart::ConditionType& rCondition)
            {
                if (rCondition.IsNot(IdentifierFlag))
                    for (auto& r_node : rCondition.GetGeometry())
                        r_node.Set(IdentifierFlag, false);
            }
        );

        bool element_to_remove;
        for (auto& elem : i_mesh->Elements()) {
            element_to_remove = true;
            for (auto& node : elem.GetGeometry()) {
                if (node.IsNot(IdentifierFlag)) {
                    element_to_remove = false;
                    break;
                }
                if (element_to_remove) elem.Set(IdentifierFlag);
            }
        }
    }

    mrModelPart.RemoveConditions(IdentifierFlag);
    mrModelPart.RemoveElements(IdentifierFlag);
    mrModelPart.RemoveNodes(IdentifierFlag);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveConditionsAndBelongingsFromAllLevels(const Flags IdentifierFlag)
{
    ModelPart& root_model_part = mrModelPart.GetRootModelPart();
    AuxiliarModelPartUtilities aux_utility = AuxiliarModelPartUtilities(root_model_part);
    aux_utility.RemoveConditionsAndBelongings(IdentifierFlag);
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::RemoveOrphanNodesFromSubModelParts()
{
    VariableUtils variable_utils;
    for (auto& r_sub_model_part : mrModelPart.SubModelParts()) {
        if (r_sub_model_part.NumberOfNodes() > 0) {
            auto& r_array_nodes = r_sub_model_part.Nodes();
            variable_utils.SetFlag(TO_ERASE, true, r_array_nodes);
            // Check orphans nodes
            if (r_sub_model_part.NumberOfElements() > 0 || r_sub_model_part.NumberOfConditions() > 0 || r_sub_model_part.NumberOfGeometries() > 0) {
                for (auto& r_elem : r_sub_model_part.Elements()) {
                    auto& r_geometry = r_elem.GetGeometry();
                    for (auto& r_node : r_geometry) {
                        r_node.Set(TO_ERASE, false);
                    }
                }
                for (auto& r_cond : r_sub_model_part.Conditions()) {
                    auto& r_geometry = r_cond.GetGeometry();
                    for (auto& r_node : r_geometry) {
                        r_node.Set(TO_ERASE, false);
                    }
                }
                const auto& r_geometries = r_sub_model_part.Geometries();
                for (auto it_geom = r_geometries.begin(); it_geom != r_geometries.end(); ++it_geom) {
                    auto& r_geometry = *((it_geom.base())->second);
                    for (auto& r_node : r_geometry) {
                        r_node.Set(TO_ERASE, false);
                    }
                }
            }
            r_sub_model_part.RemoveNodes(TO_ERASE);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

std::unordered_map<IndexType, std::vector<IndexType>> AuxiliarModelPartUtilities::RetrieveElementsNeighbourElementsIds()
{
    // Initialize the container to store the IDs of neighboring elements for each element
    std::unordered_map<IndexType, std::vector<IndexType>> elements_neighbours_elements_ids;

    // Determine dimension
    unsigned int dim = 0;
    if (mrModelPart.NumberOfElements() > 0) {
        dim = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
    }

    // Loop over each element in the model part
    for (auto& r_elem : mrModelPart.Elements()) {
        // Retrieve the neighboring elements for the current element
        auto& r_neighbours = r_elem.GetValue(NEIGHBOUR_ELEMENTS);

        // Check if the number of neighboring elements matches the expected dimension
        KRATOS_ERROR_IF_NOT(r_neighbours.size() == dim) << "Condition " << r_elem.Id() << " has not correct size solution for NEIGHBOUR_ELEMENTS " << r_neighbours.size() << " vs " << dim << std::endl;

        // Create a vector to store the IDs of neighboring elements
        std::vector<IndexType> solution;
        solution.reserve(dim);

        // Loop over each dimension and store the ID of the neighboring element
        for (unsigned int i_dim = 0; i_dim < dim; ++i_dim) {
            auto p_neighbour = r_neighbours(i_dim);
            if (p_neighbour.get()) {
                solution.push_back(p_neighbour->Id());
            }
        }

        // Shrink to fit
        solution.shrink_to_fit();

        // Insert the IDs of neighboring elements into the map with the ID of the current element as the key
        elements_neighbours_elements_ids.insert({r_elem.Id(), solution});
    }

    // Return the map containing the IDs of neighboring elements for each element
    return elements_neighbours_elements_ids;
}

/***********************************************************************************/
/***********************************************************************************/

std::unordered_map<IndexType, std::vector<IndexType>> AuxiliarModelPartUtilities::RetrieveConditionsNeighbourConditionsIds()
{
    // Initialize the container to store the IDs of neighboring conditions for each condition
    std::unordered_map<IndexType, std::vector<IndexType>> conditions_neighbours_conditions_ids;

    // Determine dimension
    unsigned int dim = 0;
    if (mrModelPart.NumberOfConditions() > 0) {
        dim = mrModelPart.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();
    }

    // Loop over each condition in the model part
    for (auto& r_cond : mrModelPart.Conditions()) {
        // Retrieve the neighboring conditions for the current condition
        auto& r_neighbours = r_cond.GetValue(NEIGHBOUR_CONDITIONS);

        // Check if the number of neighboring conditions matches the expected dimension
        KRATOS_ERROR_IF_NOT(r_neighbours.size() == dim) << "Condition " << r_cond.Id() << " has not correct size solution for NEIGHBOUR_CONDITIONS " << r_neighbours.size() << " vs " << dim << std::endl;

        // Create a vector to store the IDs of neighboring conditions
        std::vector<IndexType> solution;
        solution.reserve(dim);

        // Loop over each dimension and store the ID of the neighboring condition
        for (unsigned int i_dim = 0; i_dim < dim; ++i_dim) {
            auto p_neighbour = r_neighbours(i_dim);
            if (p_neighbour.get()) {
                solution.push_back(p_neighbour->Id());
            }
        }

        // Shrink to fit
        solution.shrink_to_fit();

        // Insert the IDs of neighboring conditions into the map with the ID of the current condition as the key
        conditions_neighbours_conditions_ids.insert({r_cond.Id(), solution});
    }

    // Return the map containing the IDs of neighboring conditions for each condition
    return conditions_neighbours_conditions_ids;
}

/***********************************************************************************/
/***********************************************************************************/

ModelPart& AuxiliarModelPartUtilities::DeepCopyModelPart(
    const std::string& rNewModelPartName,
    Model* pModel
    )
{
    // First we check if the model is not null
    if (pModel == nullptr)
        pModel = &mrModelPart.GetModel();

    // We create the new model part
    ModelPart& r_model_part = pModel->CreateModelPart(rNewModelPartName);

    /// We will copy the member variables of the model part one by one

    // We copy the buffer size (direct copy)
    r_model_part.SetBufferSize(mrModelPart.GetBufferSize());

    // We copy the process info
    r_model_part.GetProcessInfo() = mrModelPart.GetProcessInfo();

    // We copy the list of variables (using copy constructor)
    r_model_part.SetNodalSolutionStepVariablesList(Kratos::make_intrusive<VariablesList>(mrModelPart.GetNodalSolutionStepVariablesList()));
    r_model_part.SetNodalSolutionStepVariablesList();

    // We copy the tables, first using the copy constructor, and then reassigning each table so it doesn't point to the original one
    const auto& r_reference_tables = mrModelPart.Tables();
    auto& r_tables = r_model_part.Tables();
    r_tables.SetMaxBufferSize(r_reference_tables.GetMaxBufferSize());
    r_tables.SetSortedPartSize(r_reference_tables.GetSortedPartSize());
    const auto& r_reference_tables_container = r_reference_tables.GetContainer();
    auto& r_tables_container = r_tables.GetContainer();
    const IndexType number_tables = r_reference_tables_container.size();
    r_tables_container.resize(number_tables);
    const auto it_tab_begin = r_reference_tables_container.begin();
    IndexPartition<std::size_t>(number_tables).for_each([&it_tab_begin,&r_tables_container](std::size_t i) {
        auto it_tab = it_tab_begin + i;
        const auto index = it_tab->first;
        const auto& r_pointer_table = it_tab->second;
        r_tables_container[i] = std::pair<std::size_t,Table<double,double>::Pointer>(index, Kratos::make_shared<Table<double,double>>(*r_pointer_table));
    });

    // We copy the meshes (here is the heavy work)
    // NOTE: From the mesh I am not going to copy neither the Flags, neither the DataValueContainer, as those are unused and I think it is needed to open a discussion about clean up of the code and remove those derivations (multiple derivations have problems of overhead https://isocpp.org/wiki/faq/multiple-inheritance)
    // RecursiveEnsureModelPartOwnsProperties(); //NOTE: To be activated in case people doesn't create the model parts properly and the properties are not created in the model part before assigning tho the elements and conditions. For the moment I would not activate it because I don't like to patronize the code with this kind of stuff. 

    // Copy properties, first using the copy constructor, and then reassigning each table so it doesn't point to the original one
    const auto& r_reference_properties = mrModelPart.rProperties();
    auto& r_properties= r_model_part.rProperties();
    r_properties.SetMaxBufferSize(r_reference_properties.GetMaxBufferSize());
    r_properties.SetSortedPartSize(r_reference_properties.GetSortedPartSize());
    const auto& r_reference_properties_container = r_reference_properties.GetContainer();
    auto& r_properties_container = r_properties.GetContainer();
    const IndexType number_properties = r_reference_properties_container.size();
    r_properties_container.resize(number_properties);
    const auto it_prop_begin = r_reference_properties_container.begin();
    IndexPartition<std::size_t>(number_properties).for_each([&it_prop_begin,&r_properties_container](std::size_t i) {
        auto it_prop = it_prop_begin + i;
        r_properties_container[i] = Kratos::make_shared<Properties>(*(*it_prop));
    });

    // Copy nodes
    const auto& r_reference_nodes = mrModelPart.Nodes();
    auto& r_nodes = r_model_part.Nodes();
    r_nodes.SetMaxBufferSize(r_reference_nodes.GetMaxBufferSize());
    r_nodes.SetSortedPartSize(r_reference_nodes.GetSortedPartSize());
    const auto& r_reference_nodes_container = r_reference_nodes.GetContainer();
    auto& r_nodes_container = r_nodes.GetContainer();
    const IndexType number_nodes = r_reference_nodes_container.size();
    r_nodes_container.resize(number_nodes);
    const auto it_node_begin = r_reference_nodes_container.begin();
    IndexPartition<std::size_t>(number_nodes).for_each([&it_node_begin,&r_nodes_container](std::size_t i) {
        auto it_node = it_node_begin + i;
        r_nodes_container[i] = (*it_node)->Clone();
    });

    // First, before copy, we create a database of pointers of the geometries
    PointerVector<Node> points_geometry;
    std::unordered_map<Geometry<Node>::Pointer,Geometry<Node>::Pointer> geometry_pointers_database;

    // The database of elements
    const auto& r_reference_elements = mrModelPart.Elements();
    for (auto& r_elem : r_reference_elements) {
        const auto& p_old_geometry = r_elem.pGetGeometry();
        if (geometry_pointers_database.find(p_old_geometry) == geometry_pointers_database.end()) {
            const auto& p_old_points = p_old_geometry->Points();
            if (points_geometry.size() != p_old_points.size()) {
                points_geometry.resize(p_old_points.size());
            }
            for (IndexType i = 0; i < p_old_points.size(); ++i) {
                points_geometry(i) = r_nodes(p_old_points[i].Id());
            }
            const IndexType id = p_old_geometry->IsIdSelfAssigned() ? 0 : p_old_geometry->Id();
            auto p_new_geometry = p_old_geometry->Create(id,points_geometry);
            p_new_geometry->GetData() = p_old_geometry->GetData();
            geometry_pointers_database[p_old_geometry] = p_new_geometry;
        }
    }

    // The database of conditions
    const auto& r_reference_conditions = mrModelPart.Conditions();
    for (auto& r_cond : r_reference_conditions) {
        const auto& p_old_geometry = r_cond.pGetGeometry();
        if (geometry_pointers_database.find(p_old_geometry) == geometry_pointers_database.end()) {
            const auto& p_old_points = p_old_geometry->Points();
            if (points_geometry.size() != p_old_points.size()) {
                points_geometry.resize(p_old_points.size());
            }
            for (IndexType i = 0; i < p_old_points.size(); ++i) {
                points_geometry(i) = r_nodes(p_old_points[i].Id());
            }
            const IndexType id = p_old_geometry->IsIdSelfAssigned() ? 0 : p_old_geometry->Id();
            auto p_new_geometry = p_old_geometry->Create(id,points_geometry);
            p_new_geometry->GetData() = p_old_geometry->GetData();
            geometry_pointers_database[p_old_geometry] = p_new_geometry;
        }
    }

    // The database of geometries
    const auto& r_reference_geometries = mrModelPart.Geometries();
    for (auto it_geom = r_reference_geometries.begin(); it_geom != r_reference_geometries.end(); ++it_geom) {
        auto p_old_geometry = (it_geom.base())->second;
        if (geometry_pointers_database.find(p_old_geometry) == geometry_pointers_database.end()) {
            const auto& p_old_points = p_old_geometry->Points();
            if (points_geometry.size() != p_old_points.size()) {
                points_geometry.resize(p_old_points.size());
            }
            for (IndexType i = 0; i < p_old_points.size(); ++i) {
                points_geometry(i) = r_nodes(p_old_points[i].Id());
            }
            const IndexType id = p_old_geometry->IsIdSelfAssigned() ? 0 : p_old_geometry->Id();
            auto p_new_geometry = p_old_geometry->Create(id,points_geometry);
            p_new_geometry->GetData() = p_old_geometry->GetData();
            geometry_pointers_database[p_old_geometry] = p_new_geometry;
        }
    }

    // Copy elements
    auto& r_elements = r_model_part.Elements();
    DeepCopyEntities(r_model_part, r_elements, r_reference_elements, geometry_pointers_database);

    // Copy conditions
    auto& r_conditions = r_model_part.Conditions();
    DeepCopyEntities(r_model_part, r_conditions, r_reference_conditions, geometry_pointers_database);

    // Copy constraints
    // NOTE: Constraints cannot be deep copied (most of the information in implemented in derived classes as the LinearMasterSlaveConstraint), therefore we will use the Clone method of the constraint to create a new one.
    const auto& r_reference_constraints = mrModelPart.MasterSlaveConstraints();
    auto& r_constraints = r_model_part.MasterSlaveConstraints();
    r_constraints.SetMaxBufferSize(r_reference_constraints.GetMaxBufferSize());
    r_constraints.SetSortedPartSize(r_reference_constraints.GetSortedPartSize());
    const auto& r_reference_const_container = r_reference_constraints.GetContainer();
    auto& r_const_container = r_constraints.GetContainer();
    const IndexType number_constraints = r_reference_const_container.size();
    r_const_container.resize(number_constraints);
    const auto it_const_begin = r_reference_const_container.begin();
    IndexPartition<std::size_t>(number_constraints).for_each([&it_const_begin,&r_const_container](std::size_t i) {
        auto it_const = it_const_begin + i;
        r_const_container[i] = (*it_const)->Clone((*it_const)->Id());
    });

    // We copy the geometries
    for (auto it_geom = r_reference_geometries.begin(); it_geom != r_reference_geometries.end(); ++it_geom) {
        auto p_old_geometry = (it_geom.base())->second;
        r_model_part.AddGeometry(geometry_pointers_database[p_old_geometry]);
    }

    // We copy the communicator (using copy constructor)
    r_model_part.SetCommunicator(Kratos::make_shared<Communicator>(mrModelPart.GetCommunicator()));

    // We cannot copy the parent model part as it will break the concept of deep copy, which a priori assumes this is the parent model part, so nothing to do here

    // We copy the sub model parts 
    // NOTE: It is assumed that the submodelparts that working only with Id of the different entities will be enough, as we have ensured to copy everything, including the ids
    DeepCopySubModelPart(mrModelPart, r_model_part);

    // Finally the Model is set in the initial creation

    return r_model_part;
}

/***********************************************************************************/
/***********************************************************************************/

void AuxiliarModelPartUtilities::DeepCopySubModelPart(const ModelPart& rOldModelPart, ModelPart& rNewModelPart)
{
    if (rOldModelPart.NumberOfSubModelParts() == 0) {
        return;
    } else {
        for (auto& r_old_sub_model_part : rOldModelPart.SubModelParts()) {
            const std::string& r_old_sub_model_part_name = r_old_sub_model_part.Name();
            ModelPart& r_new_sub_model_part = rNewModelPart.CreateSubModelPart(r_old_sub_model_part_name);

            // Now reassigning the entities (first nodes, then elements, then conditions, then constraints)
            const IndexType number_of_nodes = r_old_sub_model_part.NumberOfNodes();
            std::vector<IndexType> index_list(number_of_nodes);

            // Iterate in the nodes
            auto& r_nodes_array = r_old_sub_model_part.Nodes();
            const auto it_node_begin = r_nodes_array.begin();
            IndexPartition<std::size_t>(number_of_nodes).for_each([&it_node_begin,&index_list](std::size_t i) {
                auto it_node = it_node_begin + i;
                index_list[i] = it_node->Id();
            });
            r_new_sub_model_part.AddNodes(index_list);

            // Iterate in the elements
            const IndexType number_of_elements = r_old_sub_model_part.NumberOfElements();
            index_list.resize(number_of_elements);
            auto& r_elements_array = r_old_sub_model_part.Elements();
            const auto it_elem_begin = r_elements_array.begin();
            IndexPartition<std::size_t>(number_of_elements).for_each([&it_elem_begin,&index_list](std::size_t i) {
                auto it_elem = it_elem_begin + i;
                index_list[i] = it_elem->Id();
            });
            r_new_sub_model_part.AddElements(index_list);

            // Iterate in the conditions
            const IndexType number_of_conditions = r_old_sub_model_part.NumberOfConditions();
            index_list.resize(number_of_conditions);
            auto& r_conditions_array = r_old_sub_model_part.Conditions();
            const auto it_cond_begin = r_conditions_array.begin();
            IndexPartition<std::size_t>(number_of_conditions).for_each([&it_cond_begin,&index_list](std::size_t i) {
                auto it_cond = it_cond_begin + i;
                index_list[i] = it_cond->Id();
            });
            r_new_sub_model_part.AddConditions(index_list);


            // Iterate in the constraints
            const IndexType number_of_master_slave_constraints = r_old_sub_model_part.NumberOfConditions();
            index_list.resize(number_of_master_slave_constraints);
            auto& r_master_slave_constraints_array = r_old_sub_model_part.MasterSlaveConstraints();
            const auto it_const_begin = r_master_slave_constraints_array.begin();
            IndexPartition<std::size_t>(number_of_master_slave_constraints).for_each([&it_const_begin,&index_list](std::size_t i) {
                auto it_const = it_const_begin + i;
                index_list[i] = it_const->Id();
            });
            r_new_sub_model_part.AddMasterSlaveConstraints(index_list);

            // TODO: Properties. The problem with the properties is that to the contrary to the entities, it is not guaranteed that a property from a sub model part will be present in the parent model part.

            // Iterate in the geometries
            const IndexType number_of_geometries = r_old_sub_model_part.NumberOfGeometries();
            index_list.resize(number_of_geometries);
            auto& r_geometries_array = r_old_sub_model_part.Geometries();
            const auto it_geo_begin = r_geometries_array.begin();
            IndexPartition<std::size_t>(number_of_geometries).for_each([&it_geo_begin,&index_list](std::size_t i) {
                auto it_geo = it_geo_begin;
                for (std::size_t j = 0; j < i; ++j) ++it_geo; // TODO: Redefine the iterators adaptors to accept arbitrary integers
                index_list[i] = it_geo->Id();
            });
            r_new_sub_model_part.AddGeometries(index_list);

            // Finally we do a loop over the submodelparts of the submodelpart to copy them (this is done recursively, so the copy will be done until there are no more submodelparts)
            DeepCopySubModelPart(r_old_sub_model_part, r_new_sub_model_part);
        }
    }
}

}  // namespace Kratos.
