//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:   Riccardo Rossi
//

#if !defined(KRATOS_MODELPART_GRAPH_UTILITIES_H_INCLUDED)
#define  KRATOS_MODELPART_GRAPH_UTILITIES_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <queue>

// External includes


// Project includes
#include "includes/define.h"
#include "containers/sparse_contiguous_row_graph.h"
#include "utilities/parallel_utilities.h"


namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// This file computes the graph representing the connectivity of a modelpart. 
/** Given a modelpart, it returns the csr_representation of its graph.
*/
class ModelPartGraphUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ModelPartGraphUtilities
    KRATOS_CLASS_POINTER_DEFINITION(ModelPartGraphUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ModelPartGraphUtilities(){}

    /// Destructor.
    virtual ~ModelPartGraphUtilities(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**The following function computes the connectivity graph (based on node.Id()-1 so that ids start in 0) for the mesh used as input.
    * @param rModelPart modelpart of which we will compute the graph
    * @retval the Graph being computed
    *
    * NOTE: this function is suboptimal if discontinuous or very large ids are employed 
    (the graph will have empty rows from 0 to the id of the largest node)
    */
    SparseContiguousRowGraph<> ComputeGraph(const ModelPart& rModelPart){

        unsigned int largest_id = 0;
        if(rModelPart.Nodes().size() != 0)
            largest_id = (rModelPart.Nodes().end()-1)->Id()-1;

        SparseContiguousRowGraph<> Agraph(largest_id+1);
        std::vector<unsigned int> tls_aux;
    
        block_for_each(rModelPart.Elements(), tls_aux, [&Agraph](const auto& rElem, auto& aux_list){
            const auto& geom = rElem.GetGeometry();
            const unsigned int nnodes = geom.size();
            if(aux_list.size() != nnodes)
                aux_list.resize(nnodes);

            for(unsigned int i=0; i<nnodes; ++i){
                aux_list[i] = geom[i].Id()-1;
            }
            Agraph.AddEntries(aux_list);
        });

        block_for_each(rModelPart.Conditions(), tls_aux, [&Agraph](const auto& rCond, auto& aux_list){
            const auto& geom = rCond.GetGeometry();
            const unsigned int nnodes = geom.size();
            if(aux_list.size() != nnodes)
                aux_list.resize(nnodes);

            for(unsigned int i=0; i<nnodes; ++i){
                aux_list[i] = geom[i].Id()-1;
            }
            Agraph.AddEntries(aux_list);
        });
        Agraph.Finalize();

        return Agraph;
    }

    /**The following function computes the connectivity graph (based on node.Id()-1 so that ids start in 0)
    * expressed as CSR arrays  for the modelpart used as input.
    * @param rModelPart modelpart of which we will compute the graph
    * @retval the row and col arrays in CSR form
    *
    * NOTE: this function is suboptimal if discontinuous or very large ids are employed (the graph will have empty rows from 0 to the id of the largest node)
    */
    std::pair<DenseVector<unsigned int>, DenseVector<unsigned int>> ComputeCSRGraph(const ModelPart& rModelPart)
    {
        std::pair<DenseVector<unsigned int>, DenseVector<unsigned int>> data;
        auto& rRowIndices = data.first;
        auto& rColIndices = data.second;
        const auto& graph = ComputeGraph(rModelPart);
        graph.ExportCSRArrays(rRowIndices,rColIndices);
        return data;
    }

    /**
    This function computes the Connected Components for a given graph, expressed in terms of its CSR representation
    it returns a pair containing the number of columns identified and the "color" associated to each of the nodes
    the returned "colors" array is such that one can directly employ the function
    VariableUtils.SetSolutionStepValue(model_part.Nodes(), colors) to set the value on the nodes in the modelpart for which it was 
    */
    
    std::pair<unsigned int, DenseVector<double>> ComputeConnectedComponents(
        const ModelPart::NodesContainerType& rNodes,
        const DenseVector<unsigned int>& rRowIndices,
        const DenseVector<unsigned int>& rColIndices
        )
    {
        unsigned int nindices = rNodes.size();

        std::unordered_map<unsigned int, int> visited(nindices);

        //set -1 as a special "not visited" value in visited - nodes will be marked as visited when visited[i]!=-1
        //note that we subtract -1 from the Id of the nodes to have all the indices to start in 0
        for(const auto& rNode : rNodes)
            visited[rNode.Id()-1] = -1;

        unsigned int color=0;
        for(unsigned int i=0; i<nindices; ++i)
        {
            unsigned int gid = (rNodes.begin()+i)->Id()-1;
            if(visited[gid] == -1 && (rRowIndices[gid+1]-rRowIndices[gid])!=0){
                unsigned int root=gid;
                BreadthFirstSearch(root, color, rRowIndices, rColIndices, visited);
                color += 1;
            }
        }

        auto number_of_colors = color;

        //prepare output so that is can be directly mapped to nodes by VariableUtils 
        //even in the case in which nodes are not numbered consecutively
        DenseVector<double> colors(nindices);
        unsigned int counter=0;
        for(const auto& rNode : rNodes)
            colors[counter++] = visited[rNode.Id()-1]; 

        return std::pair{number_of_colors, colors};
    }

    //similar to the previous version, except that a "active_nodes" array needs to be passed
    //the "active_nodes" array needs to be initialized without gaps in the same order in which the nodes are passed in the rNodes array
    std::pair<unsigned int, DenseVector<double>> ComputeConnectedComponents_ActiveNodesCheck(
        const ModelPart::NodesContainerType& rNodes,
        const DenseVector<unsigned int>& rRowIndices,
        const DenseVector<unsigned int>& rColIndices,
        const std::vector<bool>& active_nodes_list
        )
    {
        unsigned int nindices = rNodes.size();

        //here we create a active_nodes map to take into account that there may be gaps in the numbering of rNodes
        std::unordered_map<unsigned int, bool> active_nodes(nindices);
        for(unsigned int i=0; i<nindices; ++i)
            active_nodes[(rNodes.begin()+i)->Id()-1]=active_nodes_list[i];

        std::unordered_map<unsigned int, int> visited(nindices);

        //set -1 as a special "not visited" value in visited - nodes will be marked as visited when visited[i]!=-1
        //note that we subtract -1 from the Id of the nodes to have all the indices to start in 0
        for(const auto& rNode : rNodes)
            visited[rNode.Id()-1] = -1; //here we allocate all of the positions that CAN be visited

        unsigned int color=0;
        for(unsigned int i=0; i<nindices; ++i)
        {
            unsigned int gid = (rNodes.begin()+i)->Id()-1;

            if(active_nodes.find(gid)->second && 
                visited.find(gid)->second == -1 && 
                (rRowIndices[gid+1]-rRowIndices[gid])!=0){
                unsigned int root=gid;
                BreadthFirstSearch_ActiveNodesCheck(root, color, rRowIndices, rColIndices, visited,active_nodes);
                color += 1;
            }
        }

        auto number_of_colors = color;

        //prepare output so that is can be directly mapped to nodes by VariableUtils 
        //even in the case in which nodes are not numbered consecutively
        DenseVector<double> colors(nindices);
        unsigned int counter=0;
        for(const auto& rNode : rNodes)
            colors[counter++] = visited[rNode.Id()-1]; 

        return std::pair{number_of_colors, colors};
    }

    std::vector<unsigned int> ApplyMinimalScalarFixity(
        ModelPart::NodesContainerType& rNodes,
        const Variable<double>& rVar,
        const DenseVector<double>& colors,
        const unsigned int ncolors
        )
    {
        KRATOS_ERROR_IF(rNodes.size() != colors.size()) << "mismatch in the size of rNodes and of colors" << std::endl;
        //count the fixed entries for every color
        DenseVector<unsigned int> v(ncolors);
        v.clear();
        std::vector<unsigned int> fixed_ids;

        //count the fixed nodes
        for(unsigned int i=0; i<rNodes.size(); ++i)
        {
            if(colors[i]>=0 && (rNodes.begin()+i)->IsFixed(rVar)){
                v[colors[i]] += 1;
            }
        }

        //now we need to fix one node per each color (unless there are some fixed)
        for(unsigned int i=0; i<rNodes.size(); ++i)
        {
            int color = colors[i];
            if(color>=0 && v[color]==0){
                (rNodes.begin()+i)->Fix(rVar);
                fixed_ids.push_back((rNodes.begin()+i)->Id());
                v[color] += 1;
            }
        }

        return fixed_ids;
    }


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
    std::stringstream buffer;
    buffer << "ModelPartGraphUtilities" ;
    return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "ModelPartGraphUtilities";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{
    // BFS algorithm for computing connected components - simplest version - only considers connectivity
    void BreadthFirstSearch(
        const int startVertex, 
        const int color, 
        const DenseVector<unsigned int>& rRowIndices,
        const DenseVector<unsigned int>& rColIndices,
        std::unordered_map<unsigned int, int>& visited) 
    {
        std::queue<int> q;
        visited.find(startVertex)->second = color;
        q.push(startVertex);

        while (!q.empty()) {
            int currVertex = q.front();
            q.pop();
            for (unsigned int i = rRowIndices[currVertex]; i != rRowIndices[currVertex+1]; ++i) {
                int adjVertex = rColIndices[i];
                auto& item = visited.find(adjVertex)->second; 
                if (item==-1) {
                    item=color;
                    q.push(adjVertex);
                }
            }
        }
    }

    // BFS algorithm for computing connected components 
    // an edge is active only if both endpoints are marked as active in the "active_nodes" map
    void BreadthFirstSearch_ActiveNodesCheck(
        const int startVertex, 
        const int color, 
        const DenseVector<unsigned int>& rRowIndices,
        const DenseVector<unsigned int>& rColIndices,
        std::unordered_map<unsigned int, int>& visited,
        const std::unordered_map<unsigned int, bool>& active_nodes) 
    {
        std::queue<int> q;

        visited.find(startVertex)->second = color;
        q.push(startVertex);

        while (!q.empty()) {
            int currVertex = q.front();
            q.pop();
            for (unsigned int i = rRowIndices[currVertex]; i != rRowIndices[currVertex+1]; ++i) {
                int adjVertex = rColIndices[i];
                auto& item = visited.find(adjVertex)->second;
                if (item==-1 && active_nodes.find(adjVertex)->second) {
                    item = color;
                    q.push(adjVertex);
                }
            }
        }
    }


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ModelPartGraphUtilities& operator=(ModelPartGraphUtilities const& rOther) = delete;

    /// Copy constructor.
    ModelPartGraphUtilities(ModelPartGraphUtilities const& rOther) = delete;

    ///@}

}; // Class ModelPartGraphUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                ModelPartGraphUtilities& rThis){return rIStream;}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const ModelPartGraphUtilities& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MODELPART_GRAPH_UTILITIES_H_INCLUDED  defined


