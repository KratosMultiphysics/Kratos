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

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "includes/define.h"
#include "containers/sparse_contiguous_row_graph.h"


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
class KRATOS_API(KRATOS_CORE) ModelPartGraphUtilities
{
public:
    ///@name Type Definitions
    ///@{
    typedef unsigned int IndexType;

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
    static Kratos::unique_ptr<SparseContiguousRowGraph<>> ComputeGraph(const ModelPart& rModelPart);

    /**The following function computes the connectivity graph (based on node.Id()-1 so that ids start in 0)
    * expressed as CSR arrays  for the modelpart used as input.
    * @param rModelPart modelpart of which we will compute the graph
    * @retval the row and col arrays in CSR form
    *
    * NOTE: this function is suboptimal if discontinuous or very large ids are employed (the graph will have empty rows from 0 to the id of the largest node)
    */
    static std::pair<DenseVector<IndexType>, DenseVector<IndexType>> ComputeCSRGraph(const ModelPart& rModelPart);

    /**
    This function computes the Connected Components for a given graph, expressed in terms of its CSR representation
    it returns a pair containing the number of columns identified and the "color" associated to each of the nodes
    the returned "colors" array is such that one can directly employ the function
    VariableUtils.SetSolutionStepValue(model_part.Nodes(), colors) to set the value on the nodes in the modelpart for which it was
    */

    static std::pair<IndexType, DenseVector<double>> ComputeConnectedComponents(
        const ModelPart::NodesContainerType& rNodes,
        const DenseVector<IndexType>& rRowIndices,
        const DenseVector<IndexType>& rColIndices
        );

    //similar to the previous version, except that a "active_nodes" array needs to be passed
    //the "active_nodes" array needs to be initialized without gaps in the same order in which the nodes are passed in the rNodes array
    static std::pair<IndexType, DenseVector<double>> ComputeConnectedComponents_ActiveNodesCheck(
        const ModelPart::NodesContainerType& rNodes,
        const DenseVector<IndexType>& rRowIndices,
        const DenseVector<IndexType>& rColIndices,
        const std::vector<bool>& active_nodes_list
        );

    static std::vector<IndexType> ApplyMinimalScalarFixity(
        ModelPart::NodesContainerType& rNodes,
        const Variable<double>& rVar,
        const DenseVector<double>& colors,
        const IndexType ncolors
        );

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
    static void BreadthFirstSearch(
        const int startVertex,
        const int color,
        const DenseVector<IndexType>& rRowIndices,
        const DenseVector<IndexType>& rColIndices,
        std::unordered_map<IndexType, int>& rVisited);

    // BFS algorithm for computing connected components
    // an edge is active only if both endpoints are marked as active in the "active_nodes" map
    static void BreadthFirstSearch_ActiveNodesCheck(
        const int startVertex,
        const int color,
        const DenseVector<IndexType>& rRowIndices,
        const DenseVector<IndexType>& rColIndices,
        std::unordered_map<IndexType, int>& rVisited,
        const std::unordered_map<IndexType, bool>& rActiveNodes);

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



