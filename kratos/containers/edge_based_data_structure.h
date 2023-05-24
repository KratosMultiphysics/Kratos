//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once


// System includes


// Project includes
#include "containers/sparse_graph.h"
#include "includes/define.h"
#include "includes/global_pointer_variables.h"
#include "includes/model_part.h"

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

class EdgeBasedDataStructure
{
public:
    ///@name Type Definitions
    ///@{

    /// Index type definition
    using IndexType = std::size_t;

    /// Sparse graph definition
    using SparseGraphType = SparseGraph<IndexType>;

    /// Pointer definition of EdgeBasedDataStructure
    KRATOS_CLASS_POINTER_DEFINITION(EdgeBasedDataStructure);

    //TODO: Fake edge data structure to be defined later on
    struct EdgeData final {};

    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Constructor
    EdgeBasedDataStructure() = default;

    /// Destructor.
    ~EdgeBasedDataStructure() = default;

    /// Assignment operator.
    EdgeBasedDataStructure& operator=(const EdgeBasedDataStructure& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void CalculateEdgeDataStructure(const ModelPart& rModelPart)
    {
        // Create the sparse edge container graph
        SparseGraphType edges_graph;
        FillEdgesSparseGraph(rModelPart, edges_graph);
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
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "EdgeBasedDataStructure";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "EdgeBasedDataStructure" << std::endl;
        PrintData(rOStream);
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
    }

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

    void FillEdgesSparseGraph(
        const ModelPart& rModelPart,
        SparseGraphType& rEdgesSparseGraph)
    {
        // Get the edge connectivities from nodal neighbours and set the sparse graph with these
        // Note that we define edge ij such as i is the lower id between i and j
        for (auto& r_node : rModelPart.Nodes()) {
            // i-node values
            const IndexType i_id = r_node.Id();
            std::vector<IndexType> col_ids;
            
            // Loop nodal neighbours
            auto& r_node_neighs = r_node.GetValue(NEIGHBOUR_NODES);
            for (auto& r_neigh : r_node_neighs) {
                // Check ids for adding current edge
                const IndexType j_id = r_neigh.Id();
                if (i_id < j_id) {
                    col_ids.push_back(j_id);
                }
            }
            // Add edges from current node to graph
            rEdgesSparseGraph.AddEntries(i_id, col_ids);
        }

        // Finalize edge graph
        rEdgesSparseGraph.Finalize();
    }


    friend class Serializer;

    void save(Serializer &rSerializer) const
    {
    }

    void load(Serializer &rSerializer)
    {
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}
}; // Class CsrMatrix

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template <class TDataType>
inline std::istream &operator>>(
    std::istream &rIStream,
    EdgeBasedDataStructure &rThis)
{
    return rIStream;
}

/// output stream function
template <class TDataType>
inline std::ostream &operator<<(
    std::ostream &rOStream,
    const EdgeBasedDataStructure &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
///@} addtogroup block

}  // namespace Kratos.

