// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nelson Lafontaine
//                   Jordi Cotela Dalmau
//                   Riccardo Rossi
//    Co-authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "utilities/parallel_utilities.h"
#include "local_refine_geometry_mesh.hpp"

namespace Kratos
{
void LocalRefineGeometryMesh::LocalRefineMesh(
    bool RefineOnReference,
    bool InterpolateInternalVariables
    )
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(RefineOnReference && !(mModelPart.NodesBegin()->SolutionStepsDataHas(DISPLACEMENT))) << "DISPLACEMENT Variable is not in the model part -- needed if refine_on_reference = true" << std::endl;

    compressed_matrix<int> coord;                            // The matrix that stores all the index of the geometry
    std::vector<int> list_of_nodes;                          // The news nodes
    std::vector<array_1d<int, 2>> position_node;             // Edges where are the news nodes

    PointerVector< Element > new_elements;
    new_elements.reserve(20);

    // Initial renumber of nodes and elements and finding refinement level
    // TODO: Using parallel utilities here
    mCurrentRefinementLevel=0;
    std::size_t id = 1;
    for (auto it = mModelPart.NodesBegin(); it != mModelPart.NodesEnd(); ++it) {
        it->SetId(id++);
        int node_refinement_level = 0;
        if (it->Has(REFINEMENT_LEVEL)) {
            node_refinement_level = it->GetValue(REFINEMENT_LEVEL);
        } else {
            it->SetValue(REFINEMENT_LEVEL,0);
        }
        if(node_refinement_level>mCurrentRefinementLevel) {
            mCurrentRefinementLevel=node_refinement_level;
        }
    }
    mCurrentRefinementLevel++;

    IndexPartition<std::size_t>(mModelPart.NumberOfElements()).for_each([&](std::size_t i) {
        auto it = mModelPart.ElementsBegin() + i;
        it->SetId(i + 1);
    });

    IndexPartition<std::size_t>(mModelPart.NumberOfConditions()).for_each([&](std::size_t i) {
        auto it = mModelPart.ConditionsBegin() + i;
        it->SetId(i + 1);
    });

    if (RefineOnReference) {
        block_for_each(mModelPart.Nodes(), [&](Node& rNode) {
            rNode.X() = rNode.X0();
            rNode.Y() = rNode.Y0();
            rNode.Z() = rNode.Z0();
        });
    }

    this->ResetFatherNodes(mModelPart);

    /* Calling all the functions necessaries to refine the mesh */
    CSRRowMatrix(mModelPart, coord);

    SearchEdgeToBeRefined(mModelPart, coord);

    CreateListOfNewNodes(mModelPart, coord, list_of_nodes, position_node);

    CalculateCoordinateAndInsertNewNodes(mModelPart, position_node, list_of_nodes);

    EraseOldElementAndCreateNewElement(mModelPart, coord, new_elements, InterpolateInternalVariables);

    EraseOldConditionsAndCreateNew(mModelPart, coord);

    RenumeringElementsAndNodes(mModelPart, new_elements);

    // Assigning refinement level to newly created nodes
    block_for_each(mModelPart.Nodes(), [&](Node& rNode) {
        if(!rNode.Has(REFINEMENT_LEVEL)){
            rNode.SetValue(REFINEMENT_LEVEL, mCurrentRefinementLevel);
        }
    });

    if (RefineOnReference) {
        block_for_each(mModelPart.Nodes(), [&](Node& rNode) {
            const array_1d<double, 3 >& r_disp = rNode.FastGetSolutionStepValue(DISPLACEMENT);
            rNode.X() = rNode.X0() + r_disp[0];
            rNode.Y() = rNode.Y0() + r_disp[1];
            rNode.Z() = rNode.Z0() + r_disp[2];
        });
    }

    UpdateSubModelPartNodes(mModelPart);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void LocalRefineGeometryMesh::CSRRowMatrix(
    ModelPart& rModelPart,
    compressed_matrix<int>& rCoord
    )
{
    KRATOS_TRY;

    const std::size_t number_of_nodes = rModelPart.NumberOfNodes();
    rCoord.resize(number_of_nodes, number_of_nodes, false);

    for (auto it = rModelPart.NodesBegin(); it != rModelPart.NodesEnd(); ++it) {
        int index_i = it->Id() - 1; // WARNING: MESH MUST BE IN ORDER
        auto& r_neighb_nodes = it->GetValue(NEIGHBOUR_NODES);

        std::vector<unsigned int> aux(r_neighb_nodes.size());
        unsigned int active = 0;
        for (auto inode = r_neighb_nodes.begin(); inode != r_neighb_nodes.end(); inode++) {
            int index_j = inode->Id() - 1;
            if (index_j > index_i) {
                aux[active] = index_j;
                active++;
            }
        }

        std::sort(aux.begin(), aux.begin() + active);
        for (unsigned int k = 0; k < active; k++) {
            rCoord.push_back(index_i, aux[k], -1);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void LocalRefineGeometryMesh::SearchEdgeToBeRefined(
    ModelPart& rModelPart,
    compressed_matrix<int>& rCoord
    )
{
    KRATOS_TRY;

    this->SearchEdgeToBeRefinedGeneric(rModelPart.ElementsBegin(), rModelPart.ElementsEnd(), rCoord);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void LocalRefineGeometryMesh::CreateListOfNewNodes(
    ModelPart& rModelPart,
    compressed_matrix<int>& rCoord,
    std::vector<int>& rListNewNodes,
    std::vector<array_1d<int, 2 > >& rPositionNode
    )
{
    KRATOS_TRY;

    unsigned int number_of_new_nodes = 0;

    NodesArrayType& pNodes = rModelPart.Nodes();

    typedef compressed_matrix<int>::iterator1 i1_t;
    typedef compressed_matrix<int>::iterator2 i2_t;

    for (i1_t i1 = rCoord.begin1(); i1 != rCoord.end1(); ++i1) {
        for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
            if (rCoord(i2.index1(), i2.index2()) == -2) {
                number_of_new_nodes++;
            }
        }
    }

    // New ID of the nodes
    rListNewNodes.resize(number_of_new_nodes, false);
    int total_node = pNodes.size();
    for (unsigned int i = 0; i < number_of_new_nodes; i++) {
        rListNewNodes[i] = total_node + i + 1;
    }

    // Setting edges -2 to the new id of the new node
    rPositionNode.resize(number_of_new_nodes);
    unsigned int index = 0;
    for (i1_t i1 = rCoord.begin1(); i1 != rCoord.end1(); ++i1) {
        for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
            if (rCoord(i2.index1(), i2.index2()) == -2) {
                rCoord(i2.index1(), i2.index2()) = rListNewNodes[index];
                rPositionNode[index][0] = i2.index1() + 1;
                rPositionNode[index][1] = i2.index2() + 1;
                index++;
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void LocalRefineGeometryMesh::CalculateCoordinateAndInsertNewNodes(
    ModelPart& rModelPart,
    const std::vector<array_1d<int, 2>>& rPositionNode,
    const std::vector<int>& rListNewNodes
    )
{
    KRATOS_TRY;

    array_1d<double, 3 > coord_node_1;
    array_1d<double, 3 > coord_node_2;
    std::vector< array_1d<double, 3 > > coordinates_new_nodes; // The coordinate of the new nodes
    coordinates_new_nodes.resize(rPositionNode.size());
    unsigned int step_data_size = rModelPart.GetNodalSolutionStepDataSize();
    auto& r_reference_dofs = (rModelPart.NodesBegin())->GetDofs();

    // Fill an auxiliary array with a pointer to the existing nodes
    std::unordered_map< unsigned int, Node::Pointer > aux_node_list;
    for(auto it = rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); it++)
        aux_node_list[it->Id()] = ( *(it.base()) );

    for (unsigned int i = 0; i < rPositionNode.size(); i++) {
        /* Calculating the coordinate of the news nodes */
        const int node_i = rPositionNode[i][0];
        const int node_j = rPositionNode[i][1];

        auto& it_node1 = aux_node_list[node_i];
        auto& it_node2 = aux_node_list[node_j];

        noalias(coord_node_1) = it_node1->Coordinates();
        noalias(coord_node_2) = it_node2->Coordinates();

        noalias(coordinates_new_nodes[i]) = 0.50 * (coord_node_1 + coord_node_2);

        /* Inserting the news node in the model part */
        auto pnode = Node::Pointer(new Node(rListNewNodes[i], coordinates_new_nodes[i][0], coordinates_new_nodes[i][1], coordinates_new_nodes[i][2]));
        pnode->SetSolutionStepVariablesList( rModelPart.NodesBegin()->pGetVariablesList() );
        pnode->SetBufferSize(rModelPart.NodesBegin()->GetBufferSize());

        pnode->GetValue(FATHER_NODES).resize(0);
        pnode->GetValue(FATHER_NODES).push_back(Node ::WeakPointer( it_node1 ));
        pnode->GetValue(FATHER_NODES).push_back(Node ::WeakPointer( it_node2 ));

        pnode->X0() = 0.5 * (it_node1->X0() + it_node2->X0());
        pnode->Y0() = 0.5 * (it_node1->Y0() + it_node2->Y0());
        pnode->Z0() = 0.5 * (it_node1->Z0() + it_node2->Z0());

        for (auto iii = r_reference_dofs.begin(); iii != r_reference_dofs.end(); iii++) {
            auto& r_dof = **iii;
            auto p_new_dof = pnode->pAddDof(r_dof);

            const auto& r_variable = (r_dof).GetVariable();
            if (it_node1->IsFixed(r_variable) && it_node2->IsFixed(r_variable)) {
                (p_new_dof)->FixDof();
            } else {
                (p_new_dof)->FreeDof();
            }
        }

        // Interpolating the data
        unsigned int buffer_size = pnode->GetBufferSize();
        for (unsigned int step = 0; step < buffer_size; step++) {
            double* new_step_data = pnode->SolutionStepData().Data(step);
            double* step_data1 = it_node1->SolutionStepData().Data(step);
            double* step_data2 = it_node2->SolutionStepData().Data(step);

            // Copying this data in the position of the vector we are interested in
            for (unsigned int j = 0; j < step_data_size; j++) {
                new_step_data[j] = 0.5 * (step_data1[j] + step_data2[j]);
            }
        }

        aux_node_list[pnode->Id()] = pnode;
    }

    auto& r_nodes = rModelPart.Nodes();
    r_nodes.clear();
    r_nodes.reserve(aux_node_list.size());
    for(auto & it : aux_node_list) {
        r_nodes.push_back(it.second);
    }

    r_nodes.Sort();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void LocalRefineGeometryMesh::EraseOldElementAndCreateNewElement(
        ModelPart& rModelPart,
        const compressed_matrix<int>& Coord,
        PointerVector< Element >& New_Elements,
        bool InterpolateInternalVariables
                                        )
{
    KRATOS_TRY;

    KRATOS_ERROR << "Called the virtual function of LocalRefineGeometryMesh for EraseOldElementAndCreateNewElement" << std::endl;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void LocalRefineGeometryMesh::EraseOldConditionsAndCreateNew(
    ModelPart& rModelPart,
    const compressed_matrix<int>& rCoord
    )
{
    KRATOS_TRY;

    KRATOS_ERROR << "Called the virtual function of LocalRefineGeometryMesh for EraseOldConditionsAndCreateNew" << std::endl;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void LocalRefineGeometryMesh::CalculateEdges(
    Element::GeometryType& rGeometry,
    const compressed_matrix<int>& rCoord,
    int* EdgeIds,
    std::vector<int>& rAux
    )
{
    KRATOS_TRY;

    KRATOS_ERROR << "Called the virtual function of LocalRefineGeometryMesh for CalculateEdges" << std::endl;

    KRATOS_CATCH( "" );

}

/***********************************************************************************/
/***********************************************************************************/

void LocalRefineGeometryMesh::RenumeringElementsAndNodes(
    ModelPart& rModelPart,
    PointerVector< Element >& rNewElements
    )
{
    KRATOS_TRY;

    unsigned int id_node = 1;
    unsigned int id_elem = 1;

    for (auto i = rModelPart.NodesBegin(); i != rModelPart.NodesEnd(); ++i) {
        if (i->Id() != id_node) {
            i->SetId(id_node);
        }
        id_node++;
    }

    for (auto it = rModelPart.ElementsBegin(); it != rModelPart.ElementsEnd(); ++it) {
        if (it->Id() != id_elem) {
            it->SetId(id_elem);
        }
        id_elem++;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

inline void LocalRefineGeometryMesh::CreatePartition(
    unsigned int NumberOfThreads,
    const int NumberOfRows,
    vector<unsigned int>& rPartitions
    )
{
    rPartitions.resize(NumberOfThreads + 1, false);
    const unsigned int partition_size = NumberOfRows / NumberOfThreads;
    rPartitions[0] = 0;
    rPartitions[NumberOfThreads] = NumberOfRows;
    for (unsigned int i = 1; i < NumberOfThreads; i++) {
        rPartitions[i] = rPartitions[i - 1] + partition_size;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void LocalRefineGeometryMesh::UpdateSubModelPartNodes(ModelPart &rModelPart)
{
    for (auto iSubModelPart = rModelPart.SubModelPartsBegin(); iSubModelPart != rModelPart.SubModelPartsEnd(); iSubModelPart++) {
        ModelPart::NodesContainerType new_nodes;
        for (auto iNode = rModelPart.Nodes().ptr_begin(); iNode != rModelPart.Nodes().ptr_end(); iNode++) {
            auto& r_father_nodes = (*iNode)->GetValue(FATHER_NODES);
            unsigned int ParentCount = r_father_nodes.size();

            if (ParentCount > 0) {
                unsigned int ParentsInSubModelPart = 0;

                for ( auto iParent = r_father_nodes.begin(); iParent != r_father_nodes.end(); iParent++) {
                    const unsigned int parent_id = iParent->Id();
                    ModelPart::NodeIterator iFound = iSubModelPart->Nodes().find( parent_id );
                    if ( iFound != iSubModelPart->NodesEnd() )
                        ParentsInSubModelPart++;
                }

                if ( ParentCount == ParentsInSubModelPart ) {
                    new_nodes.push_back(*iNode);
                }
            }
        }
        if (new_nodes.size() > 0) {
            ModelPart& r_sub_model_part = *iSubModelPart;
            r_sub_model_part.AddNodes(new_nodes.begin(), new_nodes.end());
            new_nodes.clear();
            UpdateSubModelPartNodes(r_sub_model_part);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void LocalRefineGeometryMesh::ResetFatherNodes(ModelPart &rModelPart)
{
    block_for_each(mModelPart.Nodes(), [&](Node& rNode) {
        if (rNode.Has(FATHER_NODES)) {
            (rNode.GetValue(FATHER_NODES)).clear();
        }
    });
}

} // Namespace Kratos.
