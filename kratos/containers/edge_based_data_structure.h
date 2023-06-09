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

        bool IsBoundary() const
        {
            return mIsBoundary;
        }

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

        const array_1d<double,TDim>& GetConvectiveBoundary() const
        {
            return mNiNjNormal;
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

        void AddConvectiveBoundaryValue(
            const double FaceDomainSize,
            const array_1d<double,TDim>& rUnitNormal)
        {
            // Face weight (length in 2D and 1/3 of the area in 3D)
            double w_i;
            if constexpr (TDim == 2) {
                w_i = FaceDomainSize;
            } else {
                w_i = FaceDomainSize / 3.0;
            }

            // Convective boundary term (we take advantage of the fact that N_i = N_j = 0.5)
            for (std::size_t d = 0; d < TDim; ++d) {
                mNiNjNormal[d] += w_i * 0.25 * rUnitNormal[d];
            }

            // Set the boundary flag to true
            mIsBoundary = true;
        }

    private:

        bool mIsBoundary = false; // flag to indicate if current flag belongs to a boundary
        double mLength = 0.0; // l_{ij}
        double mMij = 0.0; // N_{i}*N_{j}
        double mLij = 0.0; // grad(N_{i})^{T}*grad(N_{j})
        array_1d<double, TDim> mNiDNj = ZeroVector(TDim); // N_{i}*grad(N_{j})
        array_1d<double, TDim> mDNiNj = ZeroVector(TDim); // grad(N_{i})*N_{j}
        array_1d<double, TDim> mNiNjNormal = ZeroVector(TDim); // NiNjn (or NiNin as N value is always 0.5)
    };

    /// Index type definition
    using IndexType = std::size_t;

    /// Size type definition
    using SizeType = std::size_t;

    /// Edge data container pointer type
    using EdgeDataPointerType = std::unique_ptr<EdgeData>;

    /// Sparse graph definition
    using SparseGraphType = SparseGraph<IndexType>;

    /// CSR storage indices vector definition
    using IndicesVectorType = std::vector<IndexType>;

    /// CSR storage off-diagonal values vector definition
    using EdgeDataVectorType = std::vector<EdgeDataPointerType>;

    /// CSR storage mass matrices vector definition
    using MassMatrixVectorType = std::vector<double>;

    /// CSR storage boundary mass matrices vector definition
    using BoundaryMassMatrixVectorType = std::vector<array_1d<double,TDim>>;

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

    const BoundaryMassMatrixVectorType& GetBoundaryMassMatrixDiagonal() const
    {
        return mBoundaryMassMatrixDiagonal;
    }

    const array_1d<double,TDim>& GetBoundaryMassMatrixDiagonal(IndexType I) const
    {
        return mBoundaryMassMatrixDiagonal[I-1];
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
    BoundaryMassMatrixVectorType mBoundaryMassMatrixDiagonal;
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
        KRATOS_ERROR_IF(mNumEdges == 0) << "No edges in model part '" << rModelPart.FullName() << "'." << std::endl;
        mEdgeData.resize(mNumEdges);
        mMassMatrixDiagonal.resize(rModelPart.NumberOfNodes());
        mLumpedMassMatrixDiagonal.resize(rModelPart.NumberOfNodes());
        mBoundaryMassMatrixDiagonal.resize(rModelPart.NumberOfNodes());

        // Initialize mass matrices diagonal values
        std::fill(mMassMatrixDiagonal.begin(), mMassMatrixDiagonal.end(), 0.0);
        std::fill(mLumpedMassMatrixDiagonal.begin(), mLumpedMassMatrixDiagonal.end(), 0.0);
        std::fill(mBoundaryMassMatrixDiagonal.begin(), mBoundaryMassMatrixDiagonal.end(), ZeroVector(TDim));

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

                    // Get ij-edge auxiliary local ids
                    // Note that these are set according to the "i lowest id" storage criterion
                    IndexType aux_i;
                    IndexType aux_j;
                    if (i_id < j_id) {
                        aux_i = i;
                        aux_j = j;
                    } else {
                        aux_i = j;
                        aux_j = i;
                    }

                    // Add the mass matrices diagonal (II) contributions
                    // Note that in here we take advantage of the fact that there is an integration point at the mid of each edge
                    const IndexType i_row_id = i_id - 1;
                    const IndexType j_row_id = j_id - 1;
                    const double aux_mass = 0.25 * domain_size / NumNodes;
                    mMassMatrixDiagonal[i_row_id] += aux_mass;
                    mMassMatrixDiagonal[j_row_id] += aux_mass;
                    mLumpedMassMatrixDiagonal[i_row_id] += 2.0 * aux_mass;
                    mLumpedMassMatrixDiagonal[j_row_id] += 2.0 * aux_mass;

                    // If not created yet, create current edge data container
                    // Note that we also set the length which is the same for all neighbour elements
                    auto& rp_edge_data = pGetEdgeData(i_id, j_id);
                    if (rp_edge_data == nullptr) {
                        rp_edge_data = std::move(Kratos::make_unique<EdgeData>());
                        rp_edge_data->SetLength(norm_2(r_geom[aux_i].Coordinates()-r_geom[aux_j].Coordinates()));
                    }

                    // Add current element off-diagonal (IJ) contributions
                    rp_edge_data->AddOffDiagonalValues(domain_size, row(DNDX,aux_i), row(DNDX,aux_j));
                }
            }
        }

        // Loop the conditions to calculate the normals in the boundary edges
        // Note that in here we are assuming that current model part has conditions in the entire skin
        KRATOS_ERROR_IF(rModelPart.GetCommunicator().GlobalNumberOfConditions() == 0) << "No conditions found in model part '" << rModelPart.FullName() << "'." << std::endl;
        array_1d<double,TDim> unit_normal;
        for (auto& r_condition : rModelPart.Conditions()) {
            // Get condition data
            const auto& r_geom = r_condition.GetGeometry();
            CalculateUnitNormal(r_geom, unit_normal);
            const double domain_size = r_geom.DomainSize();

            // Loop condition edges
            for (IndexType i = 0; i < TDim - 1; ++i) {
                const IndexType i_id = r_geom[i].Id();
                for (IndexType j = i + 1; j < TDim; ++j) {
                    const IndexType j_id = r_geom[j].Id();

                    // Add the boundary mass matrix diagonal (II) contribution
                    // Note that this already includes the unit normal contribution
                    // Also note that in here we take advantage of the fact that there is an integration point at the mid of each edge
                    const IndexType i_row_id = i_id - 1;
                    const IndexType j_row_id = j_id - 1;
                    double aux_mass;
                    if constexpr (TDim == 2) {
                        aux_mass = 0.25 * domain_size;
                    } else {
                        aux_mass = 0.25 * domain_size / 3.0;
                    }
                    noalias(mBoundaryMassMatrixDiagonal[i_row_id]) += aux_mass * unit_normal;
                    noalias(mBoundaryMassMatrixDiagonal[j_row_id]) += aux_mass * unit_normal;

                    // Get current edge data container
                    auto& rp_edge_data = pGetEdgeData(i_id, j_id);
                    rp_edge_data->AddConvectiveBoundaryValue(domain_size, unit_normal);
                }
            }
        }
    }

    void CalculateUnitNormal(
        const Geometry<Node>& rGeometry,
        array_1d<double, TDim>& rUnitNormal) const
    {
        if constexpr (TDim == 2) {
            rUnitNormal[0] = rGeometry[1].Y() - rGeometry[0].Y();
            rUnitNormal[1] = -(rGeometry[1].X() - rGeometry[0].X());
        } else {
            array_1d<double, 3> v1, v2;
            v1[0] = rGeometry[1].X() - rGeometry[0].X();
            v1[1] = rGeometry[1].Y() - rGeometry[0].Y();
            v1[2] = rGeometry[1].Z() - rGeometry[0].Z();

            v2[0] = rGeometry[2].X() - rGeometry[0].X();
            v2[1] = rGeometry[2].Y() - rGeometry[0].Y();
            v2[2] = rGeometry[2].Z() - rGeometry[0].Z();

            MathUtils<double>::CrossProduct(rUnitNormal, v1, v2);
            rUnitNormal *= 0.5;
        }

        const double n_norm = norm_2(rUnitNormal);
        KRATOS_WARNING_IF("EdgeBasedDataStructure", n_norm < 1.0e-12) << "Normal is close to zero." << std::endl;
        rUnitNormal /= n_norm;
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

    EdgeDataPointerType& pGetEdgeData(
        const IndexType GlobalIdI,
        const IndexType GlobalIdJ)
    {
        // Get ij-edge auxiliary ids (edges are stored being i the lowest nodal i and j the largest)
        const IndexType aux_i_id = std::min(GlobalIdI, GlobalIdJ);
        const IndexType aux_j_id = std::max(GlobalIdI, GlobalIdJ);

        // Get position in the column indices vector as this is the same one to be used in the values vector
        const IndexType ij_col_index = GetColumVectorIndex(aux_i_id, aux_j_id);
        KRATOS_ERROR_IF(ij_col_index < 0) << "Column index cannot be found for ij-edge " << GlobalIdI << "-" << GlobalIdJ << "." << std::endl;

        return mEdgeData[ij_col_index];
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

