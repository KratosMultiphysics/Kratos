//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes


// External includes


// Project includes
#include "containers/pointer_vector_set.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "rom_application_variables.h"
#include "rom_auxiliary_utilities.h"

namespace Kratos
{

namespace
{
    std::map<GeometryData::KratosGeometryType, std::string> AuxiliaryGeometryToConditionMap {
        {GeometryData::KratosGeometryType::Kratos_Line2D2, "LineCondition2D2N"},
        {GeometryData::KratosGeometryType::Kratos_Triangle3D3, "SurfaceCondition3D3N"},
        {GeometryData::KratosGeometryType::Kratos_Triangle3D6, "SurfaceCondition3D6N"},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4, "SurfaceCondition3D4N"},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8, "SurfaceCondition3D8N"},
        {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9, "SurfaceCondition3D9N"},
    };
}

void RomAuxiliaryUtilities::SetHRomComputingModelPart(
    const Parameters HRomWeights,
    const ModelPart& rOriginModelPart,
    ModelPart& rHRomComputingModelPart)
{
    // Ensure that the provided destination model part is empty
    rHRomComputingModelPart.Clear();

    // Auxiliary containers to save the entities involved in the HROM mesh
    // Note that we use a set for the nodes to make sure that the same node is not added by more than one element/condition
    NodesPointerSetType hrom_nodes_set;
    std::vector<Element::Pointer> hrom_elems_vect;
    std::vector<Condition::Pointer> hrom_conds_vect;

    const auto& r_elem_weights = HRomWeights["Elements"];
    hrom_elems_vect.reserve(rOriginModelPart.NumberOfElements());
    for (auto it = r_elem_weights.begin(); it != r_elem_weights.end(); ++it) {
        // Get element from origin mesh
        const IndexType elem_id = stoi(it.name());
        const auto p_elem = rOriginModelPart.pGetElement(elem_id + 1); //FIXME: WHY THIS +1?

        // Add the element to the auxiliary container and to the main HROM model part
        hrom_elems_vect.push_back(p_elem);
        rHRomComputingModelPart.AddElement(p_elem);

        // Add the element nodes to the auxiliary set and to the main HROM model part
        const auto& r_geom = p_elem->GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            NodeType::Pointer p_node = r_geom(i_node);
            hrom_nodes_set.insert(hrom_nodes_set.end(), p_node);
            rHRomComputingModelPart.AddNode(p_node);
        }
    }
    hrom_elems_vect.shrink_to_fit();

    const auto& r_cond_weights = HRomWeights["Conditions"];
    hrom_conds_vect.reserve(rOriginModelPart.NumberOfConditions());
    for (auto it = r_cond_weights.begin(); it != r_cond_weights.end(); ++it) {
        // Get the condition from origin mesh
        const IndexType cond_id = stoi(it.name());
        auto p_cond = rOriginModelPart.pGetCondition(cond_id + 1); //FIXME: WHY THIS +1?

        // Add the condition to the auxiliary container and to the main HROM model part
        hrom_conds_vect.push_back(p_cond);
        rHRomComputingModelPart.AddCondition(p_cond);

        // Add the condition nodes to the auxiliary set and to the main HROM model part
        const auto& r_geom = p_cond->GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            auto p_node = r_geom(i_node);
            hrom_nodes_set.insert(hrom_nodes_set.end(), p_node);
            rHRomComputingModelPart.AddNode(p_node);
        }
    }
    hrom_conds_vect.shrink_to_fit();

    //TODO: ADD MPC'S

    // Add properties to the HROM mesh
    // Note that we add all the properties although some of them might note be used in the HROM mesh
    auto& r_root_model_part = const_cast<ModelPart&>(rOriginModelPart).GetRootModelPart();
    auto& r_properties = r_root_model_part.rProperties();
    for (auto it_p_prop = r_properties.ptr_begin(); it_p_prop < r_properties.ptr_end(); ++it_p_prop) {
        rHRomComputingModelPart.AddProperties(*it_p_prop);
    }

    // Create and fill the HROM calculation sub model parts
    for (auto& r_orig_sub_mp : rOriginModelPart.SubModelParts()) {
        RecursiveHRomModelPartCreation(hrom_nodes_set, hrom_elems_vect, hrom_conds_vect, r_orig_sub_mp, rHRomComputingModelPart);
    }
}

void RomAuxiliaryUtilities::RecursiveHRomModelPartCreation(
    const NodesPointerSetType& rNodesSet,
    const std::vector<Element::Pointer>& rElementsVector,
    const std::vector<Condition::Pointer>& rConditionsVector,
    const ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart)
{
    // Emulate the origin submodelpart hierarchy
    auto& r_hrom_sub_mp = rDestinationModelPart.CreateSubModelPart(rOriginModelPart.Name());

    // Add nodes
    std::vector<IndexType> aux_node_ids;
    aux_node_ids.reserve(rOriginModelPart.NumberOfNodes());
    for (const auto& r_node : rOriginModelPart.Nodes()) {
        if (rNodesSet.find(r_node.Id()) != rNodesSet.end()) {
            aux_node_ids.push_back(r_node.Id());
        }
    }
    r_hrom_sub_mp.AddNodes(aux_node_ids);

    // Add elements
    std::vector<IndexType> aux_elem_ids;
    aux_elem_ids.reserve(rOriginModelPart.NumberOfElements());
    for (const auto& r_elem : rOriginModelPart.Elements()) {
        auto is_found = [&r_elem](Element::Pointer p_elem){return r_elem.Id() == p_elem->Id();};
        if (std::find_if(rElementsVector.begin(), rElementsVector.end(), is_found) != rElementsVector.end()) {
            aux_elem_ids.push_back(r_elem.Id());
        }
    }
    r_hrom_sub_mp.AddElements(aux_elem_ids);

    // Add conditions
    std::vector<IndexType> aux_cond_ids;
    aux_cond_ids.reserve(rOriginModelPart.NumberOfConditions());
    for (const auto& r_cond : rOriginModelPart.Conditions()) {
        auto is_found = [&r_cond](Condition::Pointer p_cond){return r_cond.Id() == p_cond->Id();};
        if (std::find_if(rConditionsVector.begin(), rConditionsVector.end(), is_found) != rConditionsVector.end()) {
            aux_cond_ids.push_back(r_cond.Id());
        }
    }
    r_hrom_sub_mp.AddConditions(aux_cond_ids);

    // Add properties
    auto& r_properties = const_cast<ModelPart&>(rOriginModelPart).rProperties();
    for (auto it_p_prop = r_properties.ptr_begin(); it_p_prop < r_properties.ptr_end(); ++it_p_prop) {
        r_hrom_sub_mp.AddProperties(*it_p_prop);
    }

    //TODO: ADD MPCs

    // Recursive addition
    for (auto& r_orig_sub_mp : rOriginModelPart.SubModelParts()) {
        RecursiveHRomModelPartCreation(rNodesSet, rElementsVector, rConditionsVector, r_orig_sub_mp, r_hrom_sub_mp);
    }
}

//TODO: Make it thin walled and beam compatible
void RomAuxiliaryUtilities::SetHRomVolumetricVisualizationModelPart(
    const ModelPart& rOriginModelPart,
    ModelPart& rHRomVisualizationModelPart)
{
    // Create a map for the potential skin entities
    // Key is a sorted vector with the face ids
    // Value is a tuple with a bool indicating if the entity is repeated (first) with a pointer to the origin face entity to be cloned (second)
    ElementFacesMapType element_faces_map;

    // Find the volumetric body skin
    std::vector<IndexType> bd_ids;
    GeometryType::GeometriesArrayType boundary_entities;
    for (const auto& r_elem : rOriginModelPart.Elements()) {
        // Get the geometry face ids
        const auto& r_geom = r_elem.GetGeometry();
        boundary_entities = r_geom.GenerateBoundariesEntities();

        // Loop the boundary entities
        const SizeType n_bd_entities = boundary_entities.size();
        for (IndexType i_bd_entity = 0; i_bd_entity < n_bd_entities; ++i_bd_entity) {
            // Get the boundary entity geometry
            const auto& r_bd_entity_geom = boundary_entities[i_bd_entity];

            // Set an auxiliary array with the sorted ids to be used as key
            SizeType n_nodes_bd = r_bd_entity_geom.PointsNumber();
            if (bd_ids.size() != n_nodes_bd) {
                bd_ids.resize(n_nodes_bd);
            }
            IndexType i = 0;
            for (const auto& r_node : r_bd_entity_geom) {
                bd_ids[i++] = r_node.Id();
            }
            std::sort(bd_ids.begin(), bd_ids.end());

            // Search for the current boundary entity ids
            // If not added, do the first insert in the map with a false value of the repeated flag
            // If already added, modify the existent value to flag the entity as repeated
            auto p_bd_geom = boundary_entities(i_bd_entity);
            auto it_search = element_faces_map.find(bd_ids);
            if (it_search == element_faces_map.end()) {
                auto value = std::make_pair(false, p_bd_geom);
                element_faces_map.insert(std::make_pair(bd_ids, value));
            } else {
                auto value = std::make_pair(true, p_bd_geom);
                it_search->second = value;
            }
        }
    }

    // Filter the skin entities from the face entities in the map
    NodesPointerSetType skin_nodes_set;
    std::vector<GeometryPointerType> skin_geom_prototypes;
    for (auto& r_map_entry : element_faces_map) {
        const auto value = r_map_entry.second;
        // Note that the first pair value indicates if the face entity is repeated (interior)
        if (!std::get<0>(value)) {
            // Add current boundary face to the prototypes list
            auto p_bd_geom_prot = std::get<1>(value);
            skin_geom_prototypes.push_back(p_bd_geom_prot);

            // Add current boundary face nodes to the auxiliary set
            const auto& r_geom = *p_bd_geom_prot;
            const SizeType n_face_nodes = r_geom.PointsNumber();
            for (IndexType i_node = 0; i_node < n_face_nodes; ++i_node) {
                auto p_node = r_geom(i_node);
                skin_nodes_set.insert(skin_nodes_set.end(), p_node);
            }
        }
    }

    // Add missing nodes to the HROM main model part
    std::vector<IndexType> skin_nodes_ids;
    skin_nodes_ids.reserve(skin_nodes_set.size());
    for (auto it_p_node = skin_nodes_set.ptr_begin(); it_p_node != skin_nodes_set.ptr_end(); ++it_p_node) {
        skin_nodes_ids.push_back((*it_p_node)->Id());
        rHRomVisualizationModelPart.AddNode(*it_p_node);
    }

    // Add entities to the HROM visualization model part
    std::sort(skin_nodes_ids.begin(), skin_nodes_ids.end());
    rHRomVisualizationModelPart.AddNodes(skin_nodes_ids);

    // Create fake conditions for the HROM visualization
    IndexType max_cond_id = (rHRomVisualizationModelPart.GetRootModelPart().ConditionsEnd()-1)->Id();
    const IndexType max_prop_id = (rHRomVisualizationModelPart.GetRootModelPart().PropertiesEnd()-1)->Id();
    auto p_prop = rHRomVisualizationModelPart.CreateNewProperties(max_prop_id + 1);
    for (auto it_p_geom = skin_geom_prototypes.begin(); it_p_geom != skin_geom_prototypes.end(); ++it_p_geom) {
        // Get condition type from geometry type and create new condition
        const std::string condition_name = AuxiliaryGeometryToConditionMap[(*it_p_geom)->GetGeometryType()];
        rHRomVisualizationModelPart.CreateNewCondition(condition_name, ++max_cond_id, *it_p_geom, p_prop);
    }
}

std::vector<IndexType> RomAuxiliaryUtilities::GetHRomConditionParentsIds(
    const ModelPart& rModelPart,
    std::map<std::string, std::map<IndexType, double>>& rHromWeights)
{
    std::vector<IndexType> parent_ids;
    auto& r_elem_weights = rHromWeights["Elements"];
    auto& r_cond_weights = rHromWeights["Conditions"];

    for (auto it = r_cond_weights.begin(); it != r_cond_weights.end(); ++it) {
        // Get the condition parent
        const auto& r_cond = rModelPart.GetCondition(it->first + 1); //FIXME: FIX THE + 1 --> WE NEED TO WRITE REAL IDS IN THE WEIGHTS!!
        const auto& r_neigh = r_cond.GetValue(NEIGHBOUR_ELEMENTS);
        KRATOS_ERROR_IF(r_neigh.size() == 0) << "Condition "<< r_cond.Id() <<" has no parent element assigned. Check that \'NEIGHBOUR_ELEMENTS\' have been already computed." << std::endl;

        // Add the parent to the HROM weights
        // Note that we check if the condition parent has been already added by the HROM element selection strategy
        if (r_elem_weights.find(r_neigh[0].Id() - 1) == r_elem_weights.end()) { //FIXME: FIX THE + 1 --> WE NEED TO WRITE REAL IDS IN THE WEIGHTS!!
            parent_ids.push_back(r_neigh[0].Id() - 1); //FIXME: FIX THE + 1 --> WE NEED TO WRITE REAL IDS IN THE WEIGHTS!!
        }
    }

    return parent_ids;
}

void RomAuxiliaryUtilities::ProjectRomSolutionIncrementToNodes(
    const std::vector<std::string> &rRomVariableNames,
    ModelPart &rModelPart)
{
    // Create an array with pointers to the ROM variables from the provided names
    // Note that these are assumed to be provided in the same order used to create the basis
    IndexType i_var = 0;
    const SizeType n_rom_vars = rRomVariableNames.size();
    std::vector<const Variable<double>*> rom_var_list(n_rom_vars);
    for (const auto& r_var_name : rRomVariableNames) {
        rom_var_list[i_var++] = &(KratosComponents<Variable<double>>::Get(r_var_name));
    }

    // Project the ROM solution increment onto the nodal basis and append it to the current value
    // Note that the ROM solution increment is retrieved from the root model part
    const auto& r_rom_sol_incr = rModelPart.GetRootModelPart().GetValue(ROM_SOLUTION_INCREMENT);
    block_for_each(rModelPart.Nodes(), [&rom_var_list, &r_rom_sol_incr](NodeType& rNode){
        const auto& r_rom_basis = rNode.GetValue(ROM_BASIS);
        IndexType i_var = 0;
        for (const auto& p_var : rom_var_list) {
            if (!rNode.IsFixed(*p_var)) {
                // It is important to update the values from the old one in buffer position 1
                // Otherwise the update of the nodes shared by this model part and the HROM one
                // would be accumulated to that one performed in the ROM B&S (see ProjectToFineBasis)
                const double& r_old_val = rNode.FastGetSolutionStepValue(*p_var,1);
                rNode.FastGetSolutionStepValue(*p_var) = r_old_val + inner_prod(row(r_rom_basis, i_var++), r_rom_sol_incr);
            }
        }
    });
}

} // namespace Kratos
