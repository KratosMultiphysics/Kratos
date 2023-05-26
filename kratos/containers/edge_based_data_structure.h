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
#include "utilities/geometry_utilities.h"

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

template<unsigned int TDim>
class EdgeBasedDataStructure
{
public:
    ///@name Type Definitions
    ///@{

    // Number of element nodes (note that simplicial elements are assumed)
    static constexpr std::size_t NumNodes = TDim + 1;

    //TODO: Fake edge data structure to be defined later on
    class EdgeData final
    {
    public:

        EdgeData() = default;

        EdgeData(const EdgeData& rOther) = delete;

        double GetLength() const
        {
            return mLength;
        }

        double GetMassCoefficient() const
        {
            return mMassCoefficient;
        }

        double GetLumpedMassCoefficient() const
        {
            return mLumpedMassCoefficient;
        }

        const array_1d<double,TDim>& GetFirstDerivatives() const
        {
            return mFirstDerivatives;
        }

        void SetLength(const double Length)
        {
            mLength = Length;
        }

        void AddFirstDerivatives(
            const double DomainSize,
            const array_1d<double,TDim>& rDNiDX,
            const array_1d<double,TDim>& rDNjDX)
        {
            const double w_i = DomainSize / NumNodes;
            for (std::size_t d = 0; d < TDim; ++d) {
                mFirstDerivatives[d] = 0.5*w_i*(rDNiDX[d]*0.5 + rDNjDX[d]*0.5);
            }
        }

        void AddMassCoefficients(const double DomainSize)
        {
            const double w_i = DomainSize / NumNodes;
            mMassCoefficient += 0.25 * w_i;
            mLumpedMassCoefficient += 0.5 * w_i;
        }
    
    private:

        double mLength = 0.0;
        double mMassCoefficient = 0.0;
        double mLumpedMassCoefficient = 0.0;
        array_1d<double, TDim> mFirstDerivatives = ZeroVector(TDim);
    };

    /// Index type definition
    using IndexType = std::size_t;

    /// Size type definition
    using SizeType = std::size_t;

    /// Sparse graph definition
    using SparseGraphType = SparseGraph<IndexType>;

    /// CSR storage indices vector definition
    using IndicesVectorType = std::vector<IndexType>;

    /// CSR storage values vector definition
    using ValuesVectorType = std::vector<std::unique_ptr<EdgeData>>;

    /// Pointer definition of EdgeBasedDataStructure
    KRATOS_CLASS_POINTER_DEFINITION(EdgeBasedDataStructure);

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
        // Create the edges container sparse graph
        // Note that this is used to create the CSR storage indices arrays
        SparseGraphType edges_graph;
        FillEdgesSparseGraph(rModelPart, edges_graph);
        edges_graph.ExportCSRArrays(mRowIndices, mColIndices);

        KRATOS_WATCH(mNumEdges)
        KRATOS_WATCH(mRowIndices)
        KRATOS_WATCH(mColIndices)

        // Create the edge data sparse container
        // TODO: Check the element type with the first one
        CalculateEdgeDataValues(rModelPart, edges_graph);
    }

    ///@}
    ///@name Access
    ///@{

    SizeType NumberOfEdges() const
    {
        return mNumEdges;
    }

    const IndicesVectorType& GetRowIndices() const
    {
        return mRowIndices;
    }

    const IndicesVectorType& GetColIndices() const
    {
        return mColIndices;
    }

    const ValuesVectorType& GetValues() const
    {
        return mValues;
    }

    const EdgeData& GetEdgeData(
        IndexType I,
        IndexType J) const
    {
        const IndexType ij_col_vect_index = GetColumVectorIndex(I,J);
        return *(mValues[ij_col_vect_index]);
    }

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
    
    SizeType mNumEdges;
    IndicesVectorType mRowIndices;
    IndicesVectorType mColIndices;
    ValuesVectorType mValues;

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
        mNumEdges = 0;
        for (auto& r_node : rModelPart.Nodes()) {
            // i-node values
            const IndexType i_id = r_node.Id();
            std::vector<IndexType> col_ids;
            
            // Loop nodal neighbours (j-node)
            auto& r_node_neighs = r_node.GetValue(NEIGHBOUR_NODES);
            for (auto& r_neigh : r_node_neighs) {
                // Check ids for adding current edge
                const IndexType j_id = r_neigh.Id();
                if (i_id < j_id) {
                    col_ids.push_back(j_id);
                    mNumEdges++;
                }
            }
            // Add edges from current node to graph
            if (col_ids.size() != 0) {
                KRATOS_WATCH(i_id)
                KRATOS_WATCH(col_ids)
                rEdgesSparseGraph.AddEntries(i_id, col_ids);
            }
        }

        // Finalize edge graph
        rEdgesSparseGraph.Finalize();
    }

    void CalculateEdgeDataValues(
        const ModelPart& rModelPart,
        const SparseGraphType& rEdgesSparseGraph)
    {
        // Resize values vector to the number of edges
        KRATOS_ERROR_IF(mNumEdges == 0) << "No edges in model part '" << rModelPart.FullName() << "'" << std::endl;
        mValues.resize(mNumEdges);

        // Allocate auxiliary arrays for the elementwise calculations
        double domain_size;
        array_1d<double, NumNodes> N;
        BoundedMatrix<double, NumNodes, TDim> DNDX;

        // Loop elements to calculate their corresponding edge contributions
        for (auto& r_element : rModelPart.Elements()) {
            // Getting geometry data of the element
            const auto& r_geom = r_element.GetGeometry();
            GeometryUtils::CalculateGeometryData(r_geom, DNDX, N, domain_size);

            // Loop element edges
            for (IndexType i = 0; i < NumNodes-1; ++i) {
                const IndexType i_id = r_geom[i].Id();
                for (IndexType j = i+1; j < NumNodes; ++j) {
                    const IndexType j_id = r_geom[j].Id();
                    // Check presence of current ij-edge in the sparse graph
                    if (rEdgesSparseGraph.Has(i_id, j_id)) {
                        // Get position in the column indices vector as this is the same one to be used in the values vector
                        const IndexType ij_col_index = GetColumVectorIndex(i_id, j_id);
                        std::cout << "Element " << r_element.Id() << " edge " << i_id << "-" << j_id << " with col_index " << ij_col_index << std::endl;
                        KRATOS_ERROR_IF(ij_col_index < 0) << "Column index cannot be found for ij-edge " << i_id << "-" << j_id << "." << std::endl;

                        // If not created yet, create current edge data
                        auto& rp_edge_data = mValues[ij_col_index];
                        if (rp_edge_data == nullptr) {
                            rp_edge_data = std::move(Kratos::make_unique<EdgeData>());
                        }

                        // Add current edge elementwise contributions
                        rp_edge_data->AddMassCoefficients(domain_size);
                        rp_edge_data->AddFirstDerivatives(domain_size, row(DNDX,i), row(DNDX,j));
                        rp_edge_data->SetLength(norm_2(r_geom[i].Coordinates()-r_geom[j].Coordinates()));
                    }
                }
            }
        }
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

    IndexType GetColumVectorIndex(
        const IndexType I,
        const IndexType J) const
    {
        IndexType j_col_index = -1;
        const IndexType i_row_index = mRowIndices[I];
        for (auto it = mColIndices.begin() + i_row_index; it != mColIndices.end(); ++it) {
            if (*it == J) {
                j_col_index = std::distance(mColIndices.begin(), it);
                break;
            }
        }
        return j_col_index;
    }

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
template <unsigned int TDim>
inline std::istream &operator>>(
    std::istream &rIStream,
    EdgeBasedDataStructure<TDim> &rThis)
{
    return rIStream;
}

/// output stream function
template <unsigned int TDim>
inline std::ostream &operator<<(
    std::ostream &rOStream,
    const EdgeBasedDataStructure<TDim> &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
///@} addtogroup block

}  // namespace Kratos.

