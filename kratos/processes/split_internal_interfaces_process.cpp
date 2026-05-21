//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Collaborators:   Miguel Angel Celigueta
//                   Ruben Zorrilla
//

#include "split_internal_interfaces_process.h"

namespace Kratos
{

void SplitInternalInterfacesProcess::ExecuteInitialize()
{
    KRATOS_TRY

    // Find the conditions being affected by the connectivities update
    // To do this we save the parent element of each condition in an auxiliary map
    // Note that we assume that each condition has a unique parent (i.e. there are no interior conditions)
    std::unordered_map<IndexType, IndexType> conditions_parent_map;
    for (auto& r_condition : mrModelPart.Conditions()) {
        auto& r_geom_cond = r_condition.GetGeometry();
        const SizeType n_nodes_cond = r_geom_cond.PointsNumber();
        for (auto& r_element : mrModelPart.Elements()) {
            const auto& r_geom_elem = r_element.GetGeometry();
            const SizeType n_nodes_elem = r_geom_elem.PointsNumber();
            SizeType n_coind = 0;
            for (IndexType i_node = 0; i_node < n_nodes_cond; ++i_node) {
                for (IndexType j_node = 0; j_node < n_nodes_elem; ++j_node) {
                    if (r_geom_cond[i_node].Id() == r_geom_elem[j_node].Id()) {
                        ++n_coind;
                    }
                }
            }
            if (n_coind == n_nodes_cond) {
                conditions_parent_map.insert(std::make_pair(r_element.Id(),r_condition.Id()));
                break;
            }
        }
    }

    std::size_t max_prop_id = 0;
    std::set< std::size_t> property_ids;
    for(auto& rElem : mrModelPart.Elements()) {
        std::size_t elem_prop_id = rElem.GetProperties().Id();
        property_ids.insert(elem_prop_id);
        if (elem_prop_id > max_prop_id) {
            max_prop_id = elem_prop_id;
        }
    }

    if(property_ids.size()) {
        GenericFindElementalNeighboursProcess(mrModelPart).Execute();
        for (auto it=property_ids.begin(); it!=(--property_ids.end()); ++it) {
            std::size_t id = *it;
            KRATOS_INFO("") << "Splitting the interface between the domain identified with property Id "  << id <<" and properties with bigger Ids ..."<< std::endl;
            SplitBoundary(id, ++max_prop_id, conditions_parent_map, mrModelPart);
            KRATOS_INFO("") << "Splitting the interface between the domain identified with property Id "  << id <<" and properties with bigger Ids finished!"<< std::endl;
        }
    }
    KRATOS_CATCH("");
}

void SplitInternalInterfacesProcess::SplitBoundary(
    const std::size_t PropertyIdBeingProcessed,
    const std::size_t InterfaceConditionsPropertyId,
    const std::unordered_map<IndexType, IndexType>& rConditionsParentMap,
    ModelPart& rModelPart)
{
    KRATOS_TRY

    //construct list of faces on the interface
    std::vector< Geometry<Node > > interface_faces;
    std::vector< std::pair< Geometry< Node >::Pointer, Geometry< Node >::Pointer> > neighbouring_elements;

    for(auto& rElem : mrModelPart.Elements()) {
        const auto& neighb = rElem.GetValue(NEIGHBOUR_ELEMENTS);

        for(unsigned int i=0; i<rElem.GetGeometry().size(); ++i) {

            if(rElem.GetProperties().Id() == PropertyIdBeingProcessed && neighb(i).get()!=nullptr && neighb[i].GetProperties().Id() > PropertyIdBeingProcessed) {
                auto boundary_entities = rElem.GetGeometry().GenerateBoundariesEntities();
                interface_faces.push_back(boundary_entities[i]);
                neighbouring_elements.push_back( std::make_pair(rElem.pGetGeometry(), neighb[i].pGetGeometry()) );
            }
        }
    }

    //construct list of nodes on the interface
    std::set< std::size_t > ids_on_interface;
    for(auto& geom : interface_faces) {
        for(auto& rNode : geom)
            ids_on_interface.insert(rNode.Id());
    }

    //create a list with the previous set size to be filled when looping the nodes below
    //this list will be used to add the nodes to the interface submodelpart later on
    std::size_t n_nodes_on_interface = ids_on_interface.size();
    std::vector<std::size_t> ids_on_interface_list(n_nodes_on_interface);

    //create duplicated nodes list
    std::size_t max_node_id = 0;
    for(auto& rNode : mrModelPart.Nodes()) {
        if(rNode.Id() > max_node_id) max_node_id = rNode.Id();
    }
    max_node_id++;

    std::map<std::size_t, Node::Pointer> new_nodes_map;
    std::size_t aux = 0;
    for(auto& id : ids_on_interface) {
        auto& rOrigNode = rModelPart.Nodes()[id];
        auto pNode = (*mpInterfacesSubModelPart).CreateNewNode(max_node_id++, rOrigNode);
        auto& origin_dofs = rOrigNode.GetDofs();
        for (auto it_dof = origin_dofs.begin(); it_dof != origin_dofs.end(); it_dof++) {
            pNode->pAddDof(**it_dof);
        }
        new_nodes_map[id] = pNode;
        ids_on_interface_list[aux++] = id;
    }

    //add the already existing "origin" nodes to the interface submodelpart
    //this is required in order to create the conditions in this submodelpart
    (*mpInterfacesSubModelPart).AddNodes(ids_on_interface_list);

    //now change the nodes to make the split and generate the new conditions
    std::size_t max_cond_id = 0;
    for(auto& rCond : mrModelPart.Conditions()) {
        if(rCond.Id() > max_cond_id) max_cond_id = rCond.Id();
    }
    max_cond_id++;

    Properties::Pointer p_interface_prop = (*mpInterfacesSubModelPart).CreateNewProperties(InterfaceConditionsPropertyId);
    for(std::size_t i=0; i<interface_faces.size(); ++i) {
        //do the split
        auto& pgeom = neighbouring_elements[i].second;
        for(std::size_t k=0; k<pgeom->size(); ++k) {
            auto it = new_nodes_map.find((*pgeom)[k].Id());
            if( it != new_nodes_map.end() )
                (*pgeom)(k) = it->second;
        }
        //create prism(3D) or quadrilateral(2D) as provided in the parameters
        std::vector<std::size_t> interface_condition_ids;
        auto& r_interface_faces_i = interface_faces[i];
        for(std::size_t k=0; k < r_interface_faces_i.size(); ++k)
            interface_condition_ids.push_back(r_interface_faces_i[k].Id());
        for(std::size_t k=0; k<r_interface_faces_i.size(); ++k)
            interface_condition_ids.push_back(new_nodes_map[r_interface_faces_i[k].Id()]->Id());

        (*mpInterfacesSubModelPart).CreateNewCondition(mConditionName, max_cond_id++, interface_condition_ids, p_interface_prop);
    }

    // Loop the remaining elements and conditions to update their connectivity
    // Note that in here we change the connectivity of those elements and conditions that do not have a neighbouring face
    for (auto& r_element : mrModelPart.Elements()) {
        if (r_element.GetProperties().Id() != PropertyIdBeingProcessed) {
            auto& r_geom = r_element.GetGeometry();
            const SizeType n_nodes = r_geom.PointsNumber();
            for (SizeType i_node = 0; i_node < n_nodes; ++i_node) {
                auto it = new_nodes_map.find(r_geom[i_node].Id());
                if (it != new_nodes_map.end()) {
                    r_geom(i_node) = it->second;
                }
            }

            // If the current element has a child condition, update the condition node as well
            auto it_conds = rConditionsParentMap.find(r_element.Id());
            if (it_conds != rConditionsParentMap.end()) {
                auto& r_cond_geom = mrModelPart.GetCondition(it_conds->second).GetGeometry();
                SizeType n_nodes_cond = r_cond_geom.PointsNumber();
                for (SizeType i_cond_node = 0; i_cond_node < n_nodes_cond; ++i_cond_node) {
                    auto it_new_nodes = new_nodes_map.find(r_geom[i_cond_node].Id());
                    if (it_new_nodes != new_nodes_map.end()) {
                        r_cond_geom(i_cond_node) = it_new_nodes->second;
                    }
                }
            }
        }
    }

    KRATOS_CATCH("");
}

}