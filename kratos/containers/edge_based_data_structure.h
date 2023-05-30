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

        double GetOffDiagonalConsistentMass() const
        {
            return mMij;
        }

        double GetOffDiagonalLaplacian() const
        {
            return mLij;
        }

        const array_1d<double,TDim>& GetOffDiagonalConvective() const
        {
            return mNiDNj;
        }

        const array_1d<double,TDim>& GetOffDiagonalConvectiveTranspose() const
        {
            return mDNiNj;
        }

        void SetLength(const double Length)
        {
            mLength = Length;
        }

        void AddOffDiagonalValues(
            const double DomainSize,
            const array_1d<double, TDim> &rDNiDX,
            const array_1d<double, TDim> &rDNjDX)
        {
            const double w_i = DomainSize / NumNodes;
            mMij += w_i * 0.25;
            for (std::size_t d = 0; d < TDim; ++d) {
                mLij += w_i * rDNiDX[d] * rDNjDX[d];
                mNiDNj[d] += w_i * 0.5 * rDNjDX[d];
                mDNiNj[d] += w_i * rDNiDX[d] * 0.5;
            }
        }

    private:

        double mLength = 0.0; // l_{ij}
        double mMij = 0.0; // N_{i}*N_{j}
        double mLij = 0.0; // grad(N_{i})^{T}*grad(N_{j})
        array_1d<double, TDim> mNiDNj = ZeroVector(TDim); // N_{i}*grad(N_{j})
        array_1d<double, TDim> mDNiNj = ZeroVector(TDim); // grad(N_{i})*N_{j}
    };

    /// Index type definition
    using IndexType = std::size_t;

    /// Size type definition
    using SizeType = std::size_t;

    /// Sparse graph definition
    using SparseGraphType = SparseGraph<IndexType>;

    /// CSR storage indices vector definition
    using IndicesVectorType = std::vector<IndexType>;

    /// CSR storage off-diagonal values vector definition
    using EdgeDataVectorType = std::vector<std::unique_ptr<EdgeData>>;

    /// CSR storage mass matrices vector definition
    using MassMatrixVectorType = std::vector<double>;

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

        // Create the edge data sparse container
        // TODO: Check the element type with the first one
        CalculateEdgeDataValues(rModelPart, edges_graph);
    }

    void Clear()
    {
        mNumEdges = 0;
        mRowIndices.clear();
        mColIndices.clear();
        mEdgeData.clear();
        mMassMatrixDiagonal.clear();
        mLumpedMassMatrixDiagonal.clear();
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

    const MassMatrixVectorType& GetMassMatrixDiagonal() const
    {
        return mMassMatrixDiagonal;
    }

    double GetMassMatrixDiagonal(IndexType I) const
    {
        return mMassMatrixDiagonal[I-1];
    }

    const MassMatrixVectorType& GetLumpedMassMatrixDiagonal() const
    {
        return mLumpedMassMatrixDiagonal;
    }

    double GetLumpedMassMatrixDiagonal(IndexType I) const
    {
        return mLumpedMassMatrixDiagonal[I-1];
    }

    const EdgeDataVectorType& GetEdgeData() const
    {
        return mEdgeData;
    }

    const EdgeData& GetEdgeData(
        IndexType I,
        IndexType J) const
    {
        const IndexType ij_col_vect_index = GetColumVectorIndex(I,J);
        return *(mEdgeData[ij_col_vect_index]);
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
    MassMatrixVectorType mMassMatrixDiagonal;
    MassMatrixVectorType mLumpedMassMatrixDiagonal;
    EdgeDataVectorType mEdgeData;

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
        mEdgeData.resize(mNumEdges);
        mMassMatrixDiagonal.resize(rModelPart.NumberOfNodes());
        mLumpedMassMatrixDiagonal.resize(rModelPart.NumberOfNodes());

        // Initialize mass matrices diagonal values
        std::fill(mMassMatrixDiagonal.begin(), mMassMatrixDiagonal.end(), 0.0);
        std::fill(mLumpedMassMatrixDiagonal.begin(), mLumpedMassMatrixDiagonal.end(), 0.0);

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

                    // Get ij-edge auxiliary id (edge are stored being i the lowest nodal i and j the largest)
                    const IndexType aux_i_id = std::min(i_id, j_id);
                    const IndexType aux_j_id = std::max(i_id, j_id);

                    // Check presence of current ij-edge in the sparse graph
                    if (rEdgesSparseGraph.Has(aux_i_id, aux_j_id)) {
                        // Add the mass matrices diagonal (II) contributions
                        // Note that in here we take advantage of the fact that there is an integration point at the mid of each edge
                        const IndexType i_row_id = i_id - 1;
                        const double aux_mass = 0.25 * domain_size / NumNodes;
                        mMassMatrixDiagonal[i_row_id] = aux_mass;
                        mLumpedMassMatrixDiagonal[i_row_id] = 2.0 * aux_mass;

                        // Get position in the column indices vector as this is the same one to be used in the values vector
                        const IndexType ij_col_index = GetColumVectorIndex(aux_i_id, aux_j_id);
                        KRATOS_ERROR_IF(ij_col_index < 0) << "Column index cannot be found for ij-edge " << i_id << "-" << j_id << "." << std::endl;

                        // If not created yet, create current edge data container
                        // Note that we also set the length which is the same for all neighbour elements
                        auto& rp_edge_data = mEdgeData[ij_col_index];
                        if (rp_edge_data == nullptr) {
                            rp_edge_data = std::move(Kratos::make_unique<EdgeData>());
                            rp_edge_data->SetLength(norm_2(r_geom[i].Coordinates()-r_geom[j].Coordinates()));
                        }

                        // Add current element off-diagonal (IJ) contributions
                        rp_edge_data->AddOffDiagonalValues(domain_size, row(DNDX,i), row(DNDX,j));
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

