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
//

// System includes
#include <queue>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/model_part_graph_utilities.h"

namespace Kratos
{

Kratos::unique_ptr<SparseContiguousRowGraph<>> ModelPartGraphUtilities::ComputeGraph(const ModelPart& rModelPart)
{
    ModelPartGraphUtilities::IndexType largest_id = 0;
    if(rModelPart.Nodes().size() != 0)
        largest_id = (rModelPart.Nodes().end()-1)->Id()-1;

    auto Agraph = Kratos::make_unique<SparseContiguousRowGraph<>>(largest_id+1);
    std::vector<ModelPartGraphUtilities::IndexType> tls_aux;

    const auto adding_entries_lambda = [&Agraph](const auto& rEntity, auto& aux_list) {
        const auto& r_geom = rEntity.GetGeometry();
        const ModelPartGraphUtilities::IndexType nnodes = r_geom.size();
        if(aux_list.size() != nnodes)
            aux_list.resize(nnodes);

        for(ModelPartGraphUtilities::IndexType i=0; i<nnodes; ++i) {
            aux_list[i] = r_geom[i].Id()-1;
        }
        Agraph->AddEntries(aux_list);
    };

    block_for_each(rModelPart.Elements(), tls_aux, adding_entries_lambda);
    block_for_each(rModelPart.Conditions(), tls_aux, adding_entries_lambda);

    Agraph->Finalize();

    return Agraph;
}

std::pair<DenseVector<ModelPartGraphUtilities::IndexType>, DenseVector<ModelPartGraphUtilities::IndexType>> ModelPartGraphUtilities::ComputeCSRGraph(const ModelPart& rModelPart)
{
    std::pair<DenseVector<ModelPartGraphUtilities::IndexType>, DenseVector<ModelPartGraphUtilities::IndexType>> data;
    auto& rRowIndices = data.first;
    auto& rColIndices = data.second;
    const auto graph = ComputeGraph(rModelPart);
    graph->ExportCSRArrays(rRowIndices,rColIndices);
    return data;
}

std::pair<ModelPartGraphUtilities::IndexType, DenseVector<double>> ModelPartGraphUtilities::ComputeConnectedComponents(
    const ModelPart::NodesContainerType& rNodes,
    const DenseVector<ModelPartGraphUtilities::IndexType>& rRowIndices,
    const DenseVector<ModelPartGraphUtilities::IndexType>& rColIndices
        )
{
    const ModelPartGraphUtilities::IndexType nindices = rNodes.size();

    std::unordered_map<ModelPartGraphUtilities::IndexType, int> visited(nindices);

    //set -1 as a special "not visited" value in visited - nodes will be marked as visited when visited[i]!=-1
    //note that we subtract -1 from the Id of the nodes to have all the indices to start in 0
    for(const auto& r_node : rNodes)
        visited[r_node.Id()-1] = -1;

    ModelPartGraphUtilities::IndexType color=0;
    for(ModelPartGraphUtilities::IndexType i=0; i<nindices; ++i) {
        const ModelPartGraphUtilities::IndexType gid = (rNodes.begin()+i)->Id()-1;
        if(visited[gid] == -1 && (rRowIndices[gid+1]-rRowIndices[gid])!=0) {
            ModelPartGraphUtilities::IndexType root=gid;
            BreadthFirstSearch(root, color, rRowIndices, rColIndices, visited);
            color += 1;
        }
    }

    auto number_of_colors = color;

    //prepare output so that is can be directly mapped to nodes by VariableUtils
    //even in the case in which nodes are not numbered consecutively
    DenseVector<double> colors(nindices);
    ModelPartGraphUtilities::IndexType counter=0;
    for(const auto& r_node : rNodes)
        colors[counter++] = visited[r_node.Id()-1];

    return std::pair{number_of_colors, colors};
}

std::pair<ModelPartGraphUtilities::IndexType, DenseVector<double>> ModelPartGraphUtilities::ComputeConnectedComponentsWithActiveNodesCheck(
    const ModelPart::NodesContainerType& rNodes,
    const DenseVector<ModelPartGraphUtilities::IndexType>& rRowIndices,
    const DenseVector<ModelPartGraphUtilities::IndexType>& rColIndices,
    const std::vector<bool>& active_nodes_list
        )
{
    ModelPartGraphUtilities::IndexType nindices = rNodes.size();

    //here we create a active_nodes map to take into account that there may be gaps in the numbering of rNodes
    std::unordered_map<ModelPartGraphUtilities::IndexType, bool> active_nodes(nindices);
    for(ModelPartGraphUtilities::IndexType i=0; i<nindices; ++i)
        active_nodes[(rNodes.begin()+i)->Id()-1]=active_nodes_list[i];

    std::unordered_map<ModelPartGraphUtilities::IndexType, int> visited(nindices);

    //set -1 as a special "not visited" value in visited - nodes will be marked as visited when visited[i]!=-1
    //note that we subtract -1 from the Id of the nodes to have all the indices to start in 0
    for(const auto& rNode : rNodes)
        visited[rNode.Id()-1] = -1; //here we allocate all of the positions that CAN be visited

    ModelPartGraphUtilities::IndexType color=0;
    for(ModelPartGraphUtilities::IndexType i=0; i<nindices; ++i) {
        const ModelPartGraphUtilities::IndexType gid = (rNodes.begin()+i)->Id()-1;

        if(active_nodes.find(gid)->second &&
                visited.find(gid)->second == -1 &&
                (rRowIndices[gid+1]-rRowIndices[gid])!=0) {
            ModelPartGraphUtilities::IndexType root=gid;
            BreadthFirstSearchWithActiveNodesCheck(root, color, rRowIndices, rColIndices, visited,active_nodes);
            color += 1;
        }
    }

    auto number_of_colors = color;

    //prepare output so that is can be directly mapped to nodes by VariableUtils
    //even in the case in which nodes are not numbered consecutively
    DenseVector<double> colors(nindices);
    ModelPartGraphUtilities::IndexType counter=0;
    for(const auto& r_node : rNodes)
        colors[counter++] = visited[r_node.Id()-1];

    return std::pair{number_of_colors, colors};
}

std::vector<ModelPartGraphUtilities::IndexType> ModelPartGraphUtilities::ApplyMinimalScalarFixity(
    ModelPart::NodesContainerType& rNodes,
    const Variable<double>& rVar,
    const DenseVector<double>& colors,
    const ModelPartGraphUtilities::IndexType ncolors
)
{
    KRATOS_ERROR_IF(rNodes.size() != colors.size()) << "mismatch in the size of rNodes and of colors" << std::endl;
    //count the fixed entries for every color
    DenseVector<ModelPartGraphUtilities::IndexType> v(ncolors);
    v.clear();
    std::vector<ModelPartGraphUtilities::IndexType> fixed_ids;

    //count the fixed nodes
    for(ModelPartGraphUtilities::IndexType i=0; i<rNodes.size(); ++i) {
        if(colors[i]>=0 && (rNodes.begin()+i)->IsFixed(rVar)) {
            v[colors[i]] += 1;
        }
    }

    //now we need to fix one node per each color (unless there are some fixed)
    for(ModelPartGraphUtilities::IndexType i=0; i<rNodes.size(); ++i) {
        int color = colors[i];
        if(color>=0 && v[color]==0) {
            (rNodes.begin()+i)->Fix(rVar);
            fixed_ids.push_back((rNodes.begin()+i)->Id());
            v[color] += 1;
        }
    }

    return fixed_ids;
}

void ModelPartGraphUtilities::BreadthFirstSearch(
    const int startVertex,
    const int color,
    const DenseVector<ModelPartGraphUtilities::IndexType>& rRowIndices,
    const DenseVector<ModelPartGraphUtilities::IndexType>& rColIndices,
    std::unordered_map<ModelPartGraphUtilities::IndexType, int>& visited)
{
    std::queue<int> q;
    visited.find(startVertex)->second = color;
    q.push(startVertex);

    while (!q.empty()) {
        int currVertex = q.front();
        q.pop();
        for (ModelPartGraphUtilities::IndexType i = rRowIndices[currVertex]; i != rRowIndices[currVertex+1]; ++i) {
            int adjVertex = rColIndices[i];
            auto& item = visited.find(adjVertex)->second;
            if (item==-1) {
                item=color;
                q.push(adjVertex);
            }
        }
    }
}

void ModelPartGraphUtilities::BreadthFirstSearchWithActiveNodesCheck(
    const int startVertex,
    const int color,
    const DenseVector<ModelPartGraphUtilities::IndexType>& rRowIndices,
    const DenseVector<ModelPartGraphUtilities::IndexType>& rColIndices,
    std::unordered_map<ModelPartGraphUtilities::IndexType, int>& visited,
    const std::unordered_map<ModelPartGraphUtilities::IndexType, bool>& active_nodes)
{
    std::queue<int> q;

    visited.find(startVertex)->second = color;
    q.push(startVertex);

    while (!q.empty()) {
        int currVertex = q.front();
        q.pop();
        for (ModelPartGraphUtilities::IndexType i = rRowIndices[currVertex]; i != rRowIndices[currVertex+1]; ++i) {
            int adjVertex = rColIndices[i];
            auto& item = visited.find(adjVertex)->second;
            if (item==-1 && active_nodes.find(adjVertex)->second) {
                item = color;
                q.push(adjVertex);
            }
        }
    }
}

}  // namespace Kratos.


