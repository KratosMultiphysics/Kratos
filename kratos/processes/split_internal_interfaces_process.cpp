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
//

#include "split_internal_interfaces_process.h"

namespace Kratos
{

void SplitInternalInterfacesProcess::ExecuteInitialize() {
    KRATOS_TRY
    std::set< std::size_t> property_ids;

    for(auto& rElem : mrModelPart.Elements()) {
        property_ids.insert( rElem.GetProperties().Id() );
    }

    if(property_ids.size()) {
        std::size_t domain_size = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension(); //TODO: this may not be a very good solution.
        FindElementalNeighboursProcess(mrModelPart, domain_size).Execute();
        for (auto it=property_ids.begin(); it!=(--property_ids.end()); ++it) {
            std::size_t id = *it;
            KRATOS_INFO("") << "Splitting the interface between the domain identified with property Id "  << id <<" and properties with bigger Ids ..."<< std::endl;
            SplitBoundary(id, mrModelPart);
            KRATOS_INFO("") << "Splitting the interface between the domain identified with property Id "  << id <<" and properties with bigger Ids finished!"<< std::endl;
        }
    }
    KRATOS_CATCH("");
}

void SplitInternalInterfacesProcess::SplitBoundary(const std::size_t PropertyIdBeingProcessed, ModelPart& rModelPart) {
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

    //create duplicated nodes list
    std::size_t max_node_id = 0;
    for(auto& rNode : mrModelPart.Nodes()) {
        if(rNode.Id() > max_node_id) max_node_id = rNode.Id();
    }
    max_node_id++;

    std::map<std::size_t, Node<3>::Pointer> new_nodes_map;
    for(auto& id : ids_on_interface) {
        auto& rOrigNode = rModelPart.Nodes()[id];
        auto pNode = mrModelPart.CreateNewNode(max_node_id++, rOrigNode );
        auto& origin_dofs = rOrigNode.GetDofs();
        for (auto it_dof = origin_dofs.begin(); it_dof != origin_dofs.end(); it_dof++) {
            pNode->pAddDof(**it_dof);
        }
        new_nodes_map[id] = pNode;
    }

    //now change the nodes to make the split and generate the new conditions
    std::size_t max_cond_id = 0;
    for(auto& rCond : mrModelPart.Conditions()) {
        if(rCond.Id() > max_cond_id) max_cond_id = rCond.Id();
    }
    max_cond_id++;

    Properties::Pointer pInterfaceProp = mrModelPart.pGetProperties(1); //TODO: understand if the property 1 is what we want
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
        for(std::size_t k=0; k<interface_faces[i].size(); ++k)
            interface_condition_ids.push_back(interface_faces[i][k].Id());
        for(std::size_t k=0; k<interface_faces[i].size(); ++k)
            interface_condition_ids.push_back(new_nodes_map[interface_faces[i][k].Id()]->Id());

        rModelPart.CreateNewCondition(mConditionName, max_cond_id++, interface_condition_ids, pInterfaceProp ); //TODO: understand if the property 1 is what we want
    }
    KRATOS_CATCH("");
}

}