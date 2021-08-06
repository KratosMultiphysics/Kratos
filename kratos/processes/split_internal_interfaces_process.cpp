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
    std::set< std::size_t> property_ids;

    std::size_t max_prop_id = 0;
    for(auto& rElem : mrModelPart.Elements()) {
        std::size_t elem_prop_id = rElem.GetProperties().Id();
        property_ids.insert(elem_prop_id);
        if (elem_prop_id > max_prop_id) {
            max_prop_id = elem_prop_id;
        }
    }

    if(property_ids.size()) {
        std::size_t domain_size = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension(); //TODO: this may not be a very good solution.
        FindElementalNeighboursProcess(mrModelPart, domain_size).Execute();
        for (auto it=property_ids.begin(); it!=(--property_ids.end()); ++it) {
            std::size_t id = *it;
            KRATOS_INFO("") << "Splitting the interface between the domain identified with property Id "  << id <<" and properties with bigger Ids ..."<< std::endl;
            SplitBoundary(id, ++max_prop_id, mrModelPart);
            KRATOS_INFO("") << "Splitting the interface between the domain identified with property Id "  << id <<" and properties with bigger Ids finished!"<< std::endl;
        }
    }
    KRATOS_CATCH("");
}

void SplitInternalInterfacesProcess::SplitBoundary(
    const std::size_t PropertyIdBeingProcessed,
    const std::size_t InterfaceConditionsPropertyId,
    ModelPart& rModelPart)
{
    KRATOS_TRY

    //construct list of faces on the interface
    std::vector< Geometry<Node<3> > > interface_faces;
    std::vector< std::pair< Geometry< Node<3> >::Pointer, Geometry< Node<3> >::Pointer> > neighbouring_elements;

    for(auto& rElem : mrModelPart.Elements()) {
        const auto& neighb = rElem.GetValue(NEIGHBOUR_ELEMENTS);

        for(unsigned int i=0; i<rElem.GetGeometry().size(); ++i) {

            if(rElem.GetProperties().Id() == PropertyIdBeingProcessed && neighb[i].GetProperties().Id() > PropertyIdBeingProcessed) {
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

    std::map<std::size_t, Node<3>::Pointer> new_nodes_map;
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
    KRATOS_CATCH("");
}

}