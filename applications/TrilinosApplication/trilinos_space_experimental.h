//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

/* Trilinos includes */
#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_FEMultiVector.hpp>
#include <Tpetra_FECrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <TpetraExt_MatrixMatrix.hpp>
// #include <TpetraExt_TripleMatrixMultiply.hpp> // NOTE: TripleMatrixMultiply is giving segmentation fault in the version of Trilinos available in Ubuntu 24.04. We cannot use it for the BDB' operation, and instead must do two separate multiplications with an intermediate CrsMatrix. This is less efficient than a single triple-matrix multiply, but it is a temporary workaround until we can upgrade to a newer Trilinos version with a working TpetraExt_TripleMatrixMultiply_def implementation.
// NOTE: MatrixMarket_Tpetra.hpp is intentionally NOT included here — see
// custom_utilities/auxiliary_matrix_market.h for the rationale.

// Project includes
#include "trilinos_application.h"
#include "custom_utilities/auxiliary_matrix_market.h"
#include "includes/ublas_interface.h"
#include "spaces/ublas_space.h"
#include "includes/data_communicator.h"
#include "mpi/includes/mpi_data_communicator.h"

namespace Kratos
{
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

/**
 * @class TrilinosSpaceExperimental
 * @ingroup TrilinosApplication
 * @brief The space adapted for Trilinos vectors and matrices (TPetra)
 * @details This class is experimental and aims to bring support for new TPetra matrices and vectors to TrilinosApplication
 * @author Vicente Mataix Ferrandiz
 * @tparam TMatrixType The matrix type considered
 * @tparam TVectorType the vector type considered
 */
template<class TMatrixType, class TVectorType>
class TrilinosSpaceExperimental
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosSpaceExperimental
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosSpaceExperimental);

    /// Definition of the data type
    using DataType = typename TMatrixType::scalar_type;

    /// Definition of the matrix type
    using MatrixType = TMatrixType;

    /// Definition of the vector type
    using VectorType = TVectorType;

    /// Definition of the index type
    using IndexType = std::size_t;

    /// Definition of the size type
    using SizeType = std::size_t;

    /// Class definition
    using ClassType = TrilinosSpaceExperimental<TMatrixType, TVectorType>;

    /// Tpetra definitions
    // Your scalar type; the type of sparse matrix entries. e.g., double.
    using ST = typename MatrixType::scalar_type;
    // Your local ordinal type; the signed integer type used to store local sparse matrix indices.  e.g., int.
    using LO = typename MatrixType::local_ordinal_type;
    // Your global ordinal type; the signed integer type used to index the matrix globally, over all processes. e.g., int, long, ptrdif_t, int64_t, ...
    using GO = typename MatrixType::global_ordinal_type;
    // The Node type.  e.g., Kokkos::DefaultNode::DefaultNodeType, defined in KokkosCompat_DefaultNode.hpp.
    using NT = typename MatrixType::node_type;

    /// Definition of the CrsMatrix type
    using CrsMatrixType = Tpetra::CrsMatrix<ST, LO, GO, NT>;

    /// Define the import/export types
    using ImportType = Tpetra::Import<LO, GO, NT>;
    using ExportType = Tpetra::Export<LO, GO, NT>;

    /// Define the map type
    using MapType = Tpetra::Map<LO, GO, NT>;
    using MapPointerType = Teuchos::RCP<const MapType>;

    /// Define the graph type
    using GraphType = Tpetra::FECrsGraph<LO, GO, NT>;
    using GraphPointerType = Teuchos::RCP<const GraphType>;

    // Define TPetra communicator
    using CommunicatorType = Teuchos::MpiComm<int>;
    using CommunicatorPointerType = Teuchos::RCP<const CommunicatorType>;

    /// Definition of the pointer types
    using MatrixPointerType = Teuchos::RCP<MatrixType>;
    using VectorPointerType = Teuchos::RCP<VectorType>;

    /// Column indices of row non-zero values type
    using ColumnViewType = typename MatrixType::nonconst_global_inds_host_view_type;

    /// Row non-zero values type
    using ValueViewType = typename MatrixType::nonconst_values_host_view_type;

    /// Some other definitions
    using DofUpdaterType = DofUpdater<ClassType>;
    using DofUpdaterPointerType = typename DofUpdater<ClassType>::UniquePointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TrilinosSpaceExperimental()
    {
    }

    /// Destructor.
    virtual ~TrilinosSpaceExperimental()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method returns the linear algebra library used
     * @return The linear algebra library, TPETRA in this case
     */
    static constexpr TrilinosLinearAlgebraLibrary LinearAlgebraLibrary()
    {
        return TrilinosLinearAlgebraLibrary::TPETRA;
    }

    /**
     * @brief This method returns the rank of the communicator
     * @param rComm The communicator considered
     * @return The rank of the communicator
     */
    inline static int GetRank(const CommunicatorType& rComm)
    {
        return rComm.getRank();
    }

    /**
     * @brief This method returns true if the pointer is null
     * @param pPointer The pointer considered
     * @return True if the pointer is null
     */
    template<class TPointerType>
    inline static bool IsNull(const TPointerType& pPointer)
    {
        return pPointer == Teuchos::null;
    }

    /**
     * @brief This method returns the map of the vector
     * @param rV The vector considered
     * @return The map of the vector
     */
    inline static const MapType& GetMap(const VectorType& rV)
    {
        return *(rV.getMap());
    }

    /**
     * @brief Returns a Tpetra map for this rank's local rows starting at FirstMyId.
     * @param rComm The MPI communicator
     * @param LocalSize The number of locally owned rows on this rank
     * @param FirstMyId The first global row index owned by this rank
     * @return A Tpetra::Map covering the locally owned GIDs [FirstMyId, FirstMyId+LocalSize)
     */
    static MapPointerType GetOrCreateMap(
        CommunicatorType& rComm,
        const IndexType LocalSize,
        const int FirstMyId)
    {
        std::vector<GO> local_ids(LocalSize);
        for (IndexType i = 0; i < LocalSize; ++i) {
            local_ids[i] = static_cast<GO>(FirstMyId + static_cast<int>(i));
        }
        return Teuchos::rcp(new MapType(
            Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
            Teuchos::ArrayView<const GO>(local_ids.data(), static_cast<int>(local_ids.size())),
            0,
            Teuchos::rcp(&rComm, false)));
    }

    /**
     * @brief This method returns the communicator of the vector
     * @param rV The vector considered
     * @return The communicator of the vector
     */
    inline static const CommunicatorType& GetCommunicator(const VectorType& rV)
    {
        return dynamic_cast<const CommunicatorType&>(*(rV.getMap()->getComm()));
    }

    /**
     * @brief This method returns the communicator of the matrix
     * @param rA The matrix considered
     * @return The communicator of the matrix
     */
    inline static const CommunicatorType& GetCommunicator(const MatrixType& rA)
    {
        return dynamic_cast<const CommunicatorType&>(*(rA.getMap()->getComm()));
    }

    /**
     * @brief Global assembly on a Tpetra FECrsMatrix — close matrix if open.
     * @param rA The matrix to assemble
     * @details After AssembleLHS opens the matrix lazily, this finalizes assembly
     * so the matrix is in a closed (fillComplete) state for solver consumption.
     */
    static void GlobalAssemble(MatrixType& rA)
    {
        // endAssembly() is an MPI collective on FECrsMatrix. All ranks must call it,
        // even idle ones that contributed no entries (and therefore never called
        // beginAssembly() via AssembleLHS). Open the matrix first if needed.
        auto p_fe_rA = dynamic_cast<MatrixType*>(&rA);
        if (p_fe_rA) {
            if (!rA.isFillActive()) {
                p_fe_rA->beginAssembly();
            }
            p_fe_rA->endAssembly();
        } else if (rA.isFillActive()) {
            rA.endAssembly();
        }
    }

    /**
     * @brief Global assembly on a Tpetra Vector.
     * @param rV The vector to assemble (no-op)
     * @note No-op: vectors in this space use a null importer (no overlap),
     * so all RHS contributions go directly to locally owned entries via
     * sumIntoGlobalValue. No FE state-machine management or cross-process
     * communication is required.
     */
    static void GlobalAssemble(VectorType& rV)
    {
        // No-op: vectors in this space use a null importer (no overlap),
        // so all RHS contributions go directly to locally owned entries via
        // sumIntoGlobalValue. No FE state-machine management or cross-process communication is required.
    }

    /**
     * @brief Manually finalizes matrix assembly.
     * @param rA The matrix to finalize
     * @details Calls endAssembly() followed by fillComplete() if the matrix is
     * still in fill-active state after endAssembly().
     * @param rA The matrix to finalize
     * @details Calls endAssembly() followed by fillComplete() if the matrix is
     * still in fill-active state after endAssembly().
     */
    static void ManualFinalize(MatrixType& rA)
    {
        rA.endAssembly();
        if (rA.isFillActive()) {
            rA.fillComplete();
        }
    }

    /**
     * @brief This method creates an empty pointer to a map
     * @return The pointer to the map
     */
    inline static MapPointerType CreateEmptyMapPointer()
    {
        return Teuchos::null;
    }

    /**
     * @brief This method creates an empty pointer to a vector
     * @return The pointer to the vector
     */
    inline static VectorPointerType CreateEmptyVectorPointer()
    {
        return VectorPointerType(nullptr);
    }

    /**
     * @brief Creates an empty VectorType from a map.
     * @details Handles both Tpetra::FEMultiVector (needs importer+numVecs) and
     *          plain Vector/MultiVector types.
     * @param pMap The map to be used for the construction of the vector
     * @return The pointer to the created vector
     */
    inline static VectorPointerType CreateVector(const MapPointerType& pMap)
    {
        if constexpr (std::is_same_v<VectorType, Tpetra::FEMultiVector<ST, LO, GO, NT>>) {
            return Teuchos::rcp(new VectorType(pMap, Teuchos::null, 1));
        } else {
            return Teuchos::rcp(new VectorType(pMap));
        }
    }

    /**
     * @brief Creates a copy of srcVector into a new heap-allocated VectorType.
     * Safe for FEMultiVector whose copy ctor is deleted; data is deep-copied via update().
     * @param rVector The source vector
     * @return The pointer to the new vector
     */
    inline static VectorPointerType CreateVectorCopy(const VectorType& rVector)
    {
        auto p = CreateVector(rVector.getMap());
        p->update(1.0, rVector, 0.0);
        return p;
    }

    /**
     * @brief This method creates an empty pointer to a vector using Tpetra communicator
     * @param pComm The Tpetra communicator
     * @return The pointer to the vector
     */
    inline static VectorPointerType CreateEmptyVectorPointer(CommunicatorPointerType pComm)
    {
        const int global_elems = 0;
        MapPointerType map = Teuchos::rcp(new MapType(global_elems, 0, pComm));
        return CreateVector(map);
    }

    /**
     * @brief This method creates an empty pointer to a vector using Tpetra communicator
     * @param rComm The Tpetra communicator
     * @return The pointer to the vector
     */
    inline static VectorPointerType CreateEmptyVectorPointer(CommunicatorType& rComm)
    {
        return CreateEmptyVectorPointer(Teuchos::rcp(&rComm, false));
    }


    /**
     * @brief Creates an empty MatrixType from a graph.
     * @param pGraph The FECrsGraph used to construct the matrix
     * @return The pointer to the created matrix
     */
    inline static MatrixPointerType CreateMatrix(const GraphPointerType& pGraph)
    {
        return Teuchos::rcp(new MatrixType(pGraph));
    }

    /**
     * @brief Creates a deep copy of a matrix.
     * @details Replicates the sparsity pattern of @p rMatrix into a freshly
     *          allocated FECrsMatrix and then copies all values via CopyMatrixValues.
     *          Safe for FECrsMatrix whose copy constructor and assignment operator
     *          are deleted/unreliable.
     * @param rMatrix The source matrix
     * @return The pointer to the new matrix (same graph, same values)
     */
    inline static MatrixPointerType CreateMatrixCopy(const MatrixType& rMatrix)
    {
        // Reproduce the sparsity pattern through a new FECrsGraph
        const auto p_row_map = rMatrix.getRowMap();
        const auto p_col_map = rMatrix.getColMap();
        const LO num_local_rows = static_cast<LO>(rMatrix.getNodeNumRows());

        // Compute max entries per row to size the FE graph allocation.
        // (FECrsGraph accepts a scalar maxNumEntriesPerRow, not a per-row array.)
        std::size_t max_entries_per_row = 0;
        for (LO i = 0; i < num_local_rows; ++i) {
            max_entries_per_row = std::max(max_entries_per_row,
                static_cast<std::size_t>(rMatrix.getNumEntriesInLocalRow(i)));
        }
        Teuchos::RCP<GraphType> p_graph = Teuchos::rcp(new GraphType(
            p_row_map, p_col_map, max_entries_per_row));

        // Insert column indices (global indexing)
        p_graph->beginAssembly();
        for (LO i = 0; i < num_local_rows; ++i) {
            const auto global_row_index = p_row_map->getGlobalElement(i);
            typename MatrixType::local_inds_host_view_type local_cols;
            typename MatrixType::values_host_view_type vals;
            rMatrix.getLocalRowView(i, local_cols, vals);
            if (local_cols.extent(0) > 0) {
                Teuchos::Array<GO> global_cols(local_cols.extent(0));
                for (std::size_t j = 0; j < static_cast<std::size_t>(local_cols.extent(0)); ++j) {
                    global_cols[j] = p_col_map->getGlobalElement(local_cols(j));
                }
                p_graph->insertGlobalIndices(global_row_index,
                    Teuchos::ArrayView<const GO>(global_cols.data(), global_cols.size()));
            }
        }
        p_graph->endAssembly();

        // Construct the new matrix from the closed graph and copy values
        auto p_new_matrix = Teuchos::rcp(new MatrixType(
            Teuchos::rcp_const_cast<const GraphType>(p_graph)));
        CopyMatrixValues(*p_new_matrix, rMatrix);
        return p_new_matrix;
    }

    /**
     * @brief This method creates an empty pointer to a matrix
     * @return The pointer to the matrix
     */
    inline static MatrixPointerType CreateEmptyMatrixPointer()
    {
        return MatrixPointerType(nullptr);
    }

    /**
     * @brief This method creates an empty pointer to a matrix using TPetra communicator
     * @param rComm The Tpetra communicator
     * @return The pointer to the matrix
     */
    inline static MatrixPointerType CreateEmptyMatrixPointer(CommunicatorPointerType pComm)
    {
        const int global_elems = 0;
        MapPointerType map = Teuchos::rcp(new MapType(global_elems, 0, pComm));
        Teuchos::RCP<GraphType> graph = Teuchos::rcp(new GraphType(map, map, 0));
        // Use fillComplete() for the CRS graph (avoids FE state machine issues with endAssembly on a
        // fresh graph whose FE state may already be considered closed)
        if (graph->isFillActive()) {
            graph->fillComplete();
        }
        return Teuchos::rcp(new MatrixType(Teuchos::rcp_const_cast<const GraphType>(graph)));
    }

    /**
     * @brief This method creates an empty pointer to a matrix using TPetra communicator
     * @param rComm The Tpetra communicator
     * @return The pointer to the matrix
     */
    inline static MatrixPointerType CreateEmptyMatrixPointer(CommunicatorType& rComm)
    {
        return CreateEmptyMatrixPointer(Teuchos::rcp(&rComm, false));
    }

    /**
     * @brief Returns size of vector rV
     * @param rV The vector considered
     * @return The size of the vector
     */
    inline static IndexType Size(const VectorType& rV)
    {
        const int size = rV.getGlobalLength();
        return size;
    }

    /**
     * @brief Returns number of rows of rM
     * @param rM The matrix considered
     * @return The number of rows of rM
     */
    inline static IndexType Size1(MatrixType const& rM)
    {
        const int size1 = rM.getGlobalNumRows();
        return size1;
    }

    /**
     * @brief Returns number of columns of rM
     * @param rM The matrix considered
     * @return The number of columns of rM
     */
    inline static IndexType Size2(MatrixType const& rM)
    {
        const int size1 = rM.getGlobalNumCols();
        return size1;
    }

    /**
     * @brief Returns the column of the matrix in the given position
     * @details rXi = rMij
     * @param j The position of the column
     * @param rM The matrix considered
     * @param rX The column considered
     * @todo Implement this method
     */
    inline static void GetColumn(
        const unsigned int j,
        const MatrixType& rM,
        VectorType& rX
        )
    {
        KRATOS_ERROR << "GetColumn method is not currently implemented" << std::endl;
    }

    /**
     * @brief Returns a copy of the matrix rX
     * @details rY = rX
     * @param rX The matrix considered
     * @param rY The copy of the matrix rX
     */
    inline static void Copy(
        const MatrixType& rX,
        MatrixType& rY
        )
    {
        // FECrsMatrix copy assignment is unreliable; use CopyMatrixValues instead
        CopyMatrixValues(rY, rX);
    }

    /**
     * @brief Returns a copy of the vector rX
     * @details rY = rX
     * @param rX The vector considered
     * @param rY The copy of the vector rX
     */
    inline static void Copy(
        const VectorType& rX,
        VectorType& rY
        )
    {
        // FEMultiVector has deleted copy assignment; use update instead
        rY.update(1.0, rX, 0.0);
    }

    /**
     * @brief Returns the product of two vectors
     * @details rX * rY
     * @param rX The first vector considered
     * @param rY The second vector considered
     */
    inline static double Dot(
        const VectorType& rX,
        const VectorType& rY
        )
    {
        // FEMultiVector doesn't have a single-value dot(); compute via local views
        auto localViewX = rX.getLocalViewHost(Tpetra::Access::ReadOnly);
        auto localViewY = rY.getLocalViewHost(Tpetra::Access::ReadOnly);
        ST localDot = Teuchos::ScalarTraits<ST>::zero();
        const auto n = rX.getLocalLength();
        for (std::size_t i = 0; i < n; ++i) {
            localDot += localViewX(i, 0) * localViewY(i, 0);
        }
        ST globalDot = Teuchos::ScalarTraits<ST>::zero();
        Teuchos::reduceAll(*rX.getMap()->getComm(), Teuchos::REDUCE_SUM, localDot, Teuchos::outArg(globalDot));
        return static_cast<double>(globalDot);
    }

    /**
     * @brief Returns the maximum value of the vector rX
     * @param rX The vector considered
     * @return The maximum value of the vector rX
     */
    inline static double Max(const VectorType& rX)
    {
        // Access the local data
        auto localVec = rX.getLocalViewHost(Tpetra::Access::ReadOnly);

        // Find the local maximum
        ST localMax = localVec(0,0);
        auto localLength = rX.getLocalLength();
        for (std::size_t i = 1; i < localLength; ++i) {
            localMax = std::max(localMax, localVec(i,0));
        }

        // Perform a global maximum reduction
        ST globalMax;
        Teuchos::reduceAll(*rX.getMap()->getComm(), Teuchos::REDUCE_MAX, localMax, Teuchos::outArg(globalMax));

        return globalMax;
    }

    /**
     * @brief Returns the minimum value of the vector rX
     * @param rX The vector considered
     * @return The minimum value of the vector rX
     */
    inline static double Min(const VectorType& rX)
    {
        // Access the local data
        auto localVec = rX.getLocalViewHost(Tpetra::Access::ReadOnly);

        // Find the local minimum
        ST localMin = localVec(0,0);
        auto localLength = rX.getLocalLength();
        for (std::size_t i = 1; i < localLength; ++i) {
            localMin = std::min(localMin, localVec(i,0));
        }

        // Perform a global minimum reduction
        ST globalMin;
        Teuchos::reduceAll(*rX.getMap()->getComm(), Teuchos::REDUCE_MIN, localMin, Teuchos::outArg(globalMin));

        return globalMin;
    }

    /**
     * @brief Returns the norm of the vector rX
     * @details ||rX||2
     * @param rX The vector considered
     * @return The norm of the vector rX
     */
    inline static double TwoNorm(const VectorType& rX)
    {
        // FEMultiVector doesn't have a scalar-returning norm2(); compute via local views
        auto localView = rX.getLocalViewHost(Tpetra::Access::ReadOnly);
        ST localSumSq = Teuchos::ScalarTraits<ST>::zero();
        const auto n = rX.getLocalLength();
        for (std::size_t i = 0; i < n; ++i) {
            localSumSq += localView(i, 0) * localView(i, 0);
        }
        ST globalSumSq = Teuchos::ScalarTraits<ST>::zero();
        Teuchos::reduceAll(*rX.getMap()->getComm(), Teuchos::REDUCE_SUM, localSumSq, Teuchos::outArg(globalSumSq));
        return std::sqrt(static_cast<double>(globalSumSq));
    }

    /**
     * @brief Returns the Frobenius norm of the matrix rA
     * @details ||rA||2
     * @param rA The matrix considered
     * @return The Frobenius norm of the matrix rA
     */
    inline static double TwoNorm(const MatrixType& rA)
    {
        // Use the built-in getFrobeniusNorm() method to compute the Frobenius norm
        return rA.getFrobeniusNorm();
    }

    /**
     * @brief Returns the multiplication of a matrix by a vector
     * @details y = A*x
     * @param rA The matrix considered
     * @param rX The vector considered
     * @param rY The result of the multiplication
     */
    inline static void Mult(
        const MatrixType& rA,
        const VectorType& rX,
        VectorType& rY
        )
    {
        rA.apply(rX, rY);
    }

    /**
     * @brief Returns the multiplication matrix-matrix
     * @details C = A*B
     * @param rA The first matrix considered
     * @param rB The second matrix considered
     * @param rC The result of the multiplication
     * @param CallFillCompleteOnResult Whether to call fillComplete on the result matrix
     * @param KeepAllHardZeros If true, keeps all hard zeros in the result matrix
     */
    inline static void Mult(
        const MatrixType& rA,
        const MatrixType& rB,
        MatrixType& rC,
        const bool CallFillCompleteOnResult = true,
        const bool KeepAllHardZeros = false
        )
    {
        KRATOS_TRY

        // Use a temporary CrsMatrix to allow for dynamic sparsity
        Teuchos::RCP<CrsMatrixType> aux_C = Teuchos::rcp(new CrsMatrixType(rC.getRowMap(), 16));
        Tpetra::MatrixMatrix::Multiply(rA, false, rB, false, *aux_C, CallFillCompleteOnResult);

        // Copy values back to rC
        // Inline copy logic from CrsMatrixType to MatrixType
        if (!rC.isFillActive()) rC.resumeFill();
        auto p_fe_rC = dynamic_cast<MatrixType*>(&rC);
        if (p_fe_rC) p_fe_rC->beginAssembly();
        for (LO i = 0; i < static_cast<LO>(aux_C->getNodeNumRows()); ++i) {
            const auto global_row_index = aux_C->getRowMap()->getGlobalElement(i);
            typename MatrixType::local_inds_host_view_type local_cols;
            typename MatrixType::values_host_view_type vals;
            aux_C->getLocalRowView(i, local_cols, vals);
            if (vals.extent(0) > 0) {
                Teuchos::Array<GO> global_cols(local_cols.extent(0));
                for (std::size_t j = 0; j < static_cast<std::size_t>(local_cols.extent(0)); ++j) {
                    global_cols[j] = aux_C->getColMap()->getGlobalElement(local_cols(j));
                }
                rC.sumIntoGlobalValues(global_row_index, static_cast<LO>(global_cols.size()), vals.data(), global_cols.data());
            }
        }
        if (p_fe_rC) p_fe_rC->endAssembly();
        if (rC.isFillActive()) rC.fillComplete();

        KRATOS_CATCH("")
    }

    /**
     * @brief Returns the transpose multiplication of a matrix by a vector
     * @details y = AT*x
     * @param rA The matrix considered
     * @param rX The vector considered
     * @param rY The result of the multiplication
     */
    inline static void TransposeMult(
        const MatrixType& rA,
        const VectorType& rX,
        VectorType& rY
        )
    {
        // Use the apply method with the transpose flag set to true
        rA.apply(rX, rY, Teuchos::TRANS);
    }

    /**
     * @brief Returns the transpose multiplication matrix-matrix
     * @details C = A*B
     * @param rA The first matrix considered
     * @param rB The second matrix considered
     * @param rC The result of the multiplication
     * @param TransposeFlag Flags to transpose the matrices
     * @param CallFillCompleteOnResult	Optional argument, defaults to true. Power users may specify this argument to be false if they DON'T want this function to call C.FillComplete. (It is often useful to allow this function to call C.FillComplete, in cases where one or both of the input matrices are rectangular and it is not trivial to know which maps to use for the domain- and range-maps.)
     * @param KeepAllHardZeros	Optional argument, defaults to false. If true, Multiply, keeps all entries in C corresponding to hard zeros. If false, the following happens by case: A*B^T, A^T*B^T - Does not store entries caused by hard zeros in C. A^T*B (unoptimized) - Hard zeros are always stored (this option has no effect) A*B, A^T*B (optimized) - Hard zeros in corresponding to hard zeros in A are not stored, There are certain cases involving reuse of C, where this can be useful.
     */
    inline static void TransposeMult(
        const MatrixType& rA,
        const MatrixType& rB,
        MatrixType& rC,
        const std::pair<bool, bool> TransposeFlag = {false, false},
        const bool CallFillCompleteOnResult = true,
        const bool KeepAllHardZeros = false
        )
    {
        KRATOS_TRY

        // Use a temporary CrsMatrix to allow for dynamic sparsity
        Teuchos::RCP<CrsMatrixType> aux_C = Teuchos::rcp(new CrsMatrixType(rC.getRowMap(), 16));
        Tpetra::MatrixMatrix::Multiply(rA, TransposeFlag.first, rB, TransposeFlag.second, *aux_C, CallFillCompleteOnResult);

        // Copy values back to rC
        // Inline copy logic from CrsMatrixType to MatrixType
        if (!rC.isFillActive()) rC.resumeFill();
        auto p_fe_rC = dynamic_cast<MatrixType*>(&rC);
        if (p_fe_rC) p_fe_rC->beginAssembly();
        for (LO i = 0; i < static_cast<LO>(aux_C->getNodeNumRows()); ++i) {
            const auto global_row_index = aux_C->getRowMap()->getGlobalElement(i);
            typename MatrixType::local_inds_host_view_type local_cols;
            typename MatrixType::values_host_view_type vals;
            aux_C->getLocalRowView(i, local_cols, vals);
            if (vals.extent(0) > 0) {
                Teuchos::Array<GO> global_cols(local_cols.extent(0));
                for (std::size_t j = 0; j < static_cast<std::size_t>(local_cols.extent(0)); ++j) {
                    global_cols[j] = aux_C->getColMap()->getGlobalElement(local_cols(j));
                }
                rC.sumIntoGlobalValues(global_row_index, static_cast<LO>(global_cols.size()), vals.data(), global_cols.data());
            }
        }
        if (p_fe_rC) p_fe_rC->endAssembly();
        if (rC.isFillActive()) rC.fillComplete();

        KRATOS_CATCH("")
    }

    /**
     * @brief Calculates the product operation B'DB
     * @param rA The resulting matrix
     * @param rD The "center" matrix
     * @param rB The matrices to be transposed
     * @param CallFillCompleteOnResult Reserved for API compatibility with TrilinosSpace; ignored in this experimental implementation.
     * @param KeepAllHardZeros Reserved for API compatibility with TrilinosSpace; ignored in this experimental implementation.
     * @param EnforceInitialGraph Reserved for API compatibility with TrilinosSpace; ignored in this experimental implementation.
     * @todo TripleMatrixMultiply::MultiplyRAP is giving segmentation fault in the version of Trilinos available in Ubuntu 24.04. We cannot use it for the BDB' operation, and instead must do two separate multiplications with an intermediate CrsMatrix. This is less efficient than a single triple-matrix multiply, but it is a temporary workaround until we can upgrade to a newer Trilinos version with a working TpetraExt_TripleMatrixMultiply_def implementation.
     */
    inline static void BtDBProductOperation(
        MatrixType& rA,
        const MatrixType& rD,
        const MatrixType& rB,
        const bool CallFillCompleteOnResult = true,
        const bool KeepAllHardZeros = false,
        const bool EnforceInitialGraph = false
        )
    {
        // Use temporary CrsMatrix for intermediate steps
        Teuchos::RCP<CrsMatrixType> aux_1 = Teuchos::rcp(new CrsMatrixType(rB.getDomainMap(), 16));
        Tpetra::MatrixMatrix::Multiply(rB, true, rD, false, *aux_1);

        Teuchos::RCP<CrsMatrixType> aux_2 = Teuchos::rcp(new CrsMatrixType(aux_1->getRowMap(), 16));
        Tpetra::MatrixMatrix::Multiply(*aux_1, false, rB, false, *aux_2);

        auto build_fe_from_crs = [&](const CrsMatrixType& rSrc) {
            std::size_t max_entries_per_row = 0;
            for (LO i = 0; i < static_cast<LO>(rSrc.getNodeNumRows()); ++i) {
                typename MatrixType::local_inds_host_view_type local_cols;
                typename MatrixType::values_host_view_type vals;
                rSrc.getLocalRowView(i, local_cols, vals);
                max_entries_per_row = std::max(max_entries_per_row, static_cast<std::size_t>(local_cols.extent(0)));
            }

            Teuchos::RCP<GraphType> p_graph = Teuchos::rcp(new GraphType(
                rSrc.getRowMap(), rSrc.getColMap(), max_entries_per_row));
            p_graph->beginAssembly();
            for (LO i = 0; i < static_cast<LO>(rSrc.getNodeNumRows()); ++i) {
                const auto global_row_index = rSrc.getRowMap()->getGlobalElement(i);
                typename MatrixType::local_inds_host_view_type local_cols;
                typename MatrixType::values_host_view_type vals;
                rSrc.getLocalRowView(i, local_cols, vals);
                if (vals.extent(0) == 0) continue;
                Teuchos::Array<GO> global_cols(local_cols.extent(0));
                for (std::size_t j = 0; j < static_cast<std::size_t>(local_cols.extent(0)); ++j) {
                    global_cols[j] = rSrc.getColMap()->getGlobalElement(local_cols(j));
                }
                p_graph->insertGlobalIndices(global_row_index,
                    Teuchos::ArrayView<const GO>(global_cols.data(), static_cast<int>(global_cols.size())));
            }
            p_graph->endAssembly();

            auto p_matrix = Teuchos::rcp(new MatrixType(Teuchos::rcp_const_cast<const GraphType>(p_graph)));
            if (p_matrix->isFillActive()) p_matrix->fillComplete();
            p_matrix->beginAssembly();
            for (LO i = 0; i < static_cast<LO>(rSrc.getNodeNumRows()); ++i) {
                const auto global_row_index = rSrc.getRowMap()->getGlobalElement(i);
                typename MatrixType::local_inds_host_view_type local_cols;
                typename MatrixType::values_host_view_type vals;
                rSrc.getLocalRowView(i, local_cols, vals);
                if (vals.extent(0) == 0) continue;
                Teuchos::Array<GO> global_cols(local_cols.extent(0));
                for (std::size_t j = 0; j < static_cast<std::size_t>(local_cols.extent(0)); ++j) {
                    global_cols[j] = rSrc.getColMap()->getGlobalElement(local_cols(j));
                }
                p_matrix->sumIntoGlobalValues(global_row_index, static_cast<LO>(global_cols.size()), vals.data(), global_cols.data());
            }
            p_matrix->endAssembly();
            if (CallFillCompleteOnResult && p_matrix->isFillActive()) p_matrix->fillComplete();
            return p_matrix;
        };

        auto p_result = build_fe_from_crs(*aux_2);
        if (EnforceInitialGraph) {
            auto p_combined_graph = CombineMatricesGraphs(rA, *p_result);
            auto p_enforced = Teuchos::rcp(new MatrixType(Teuchos::rcp_const_cast<const GraphType>(p_combined_graph)));
            CopyMatrixValues(*p_enforced, *p_result);
            p_result = p_enforced;
        }

        CopyMatrixValues(rA, *p_result);
    }

    /**
     * @brief Calculates the product operation BDB'
     * @param rA The resulting matrix
     * @param rD The "center" matrix
     * @param rB The matrices to be transposed
     * @param CallFillCompleteOnResult Reserved for API compatibility with TrilinosSpace; ignored in this experimental implementation.
     * @param KeepAllHardZeros Reserved for API compatibility with TrilinosSpace; ignored in this experimental implementation.
     * @param EnforceInitialGraph Reserved for API compatibility with TrilinosSpace; ignored in this experimental implementation.
     * @todo TripleMatrixMultiply::MultiplyRAP is giving segmentation fault in the version of Trilinos available in Ubuntu 24.04. We cannot use it for the BDB' operation, and instead must do two separate multiplications with an intermediate CrsMatrix. This is less efficient than a single triple-matrix multiply, but it is a temporary workaround until we can upgrade to a newer Trilinos version with a working TpetraExt_TripleMatrixMultiply_def implementation.
     */
    inline static void BDBtProductOperation(
        MatrixType& rA,
        const MatrixType& rD,
        const MatrixType& rB,
        const bool CallFillCompleteOnResult = true,
        const bool KeepAllHardZeros = false,
        const bool EnforceInitialGraph = false
        )
    {
        // Use temporary CrsMatrix for intermediate steps
        Teuchos::RCP<CrsMatrixType> aux_1 = Teuchos::rcp(new CrsMatrixType(rB.getRowMap(), 16));
        Tpetra::MatrixMatrix::Multiply(rB, false, rD, false, *aux_1);

        Teuchos::RCP<CrsMatrixType> aux_2 = Teuchos::rcp(new CrsMatrixType(aux_1->getRowMap(), 16));
        Tpetra::MatrixMatrix::Multiply(*aux_1, false, rB, true, *aux_2);

        // Copy aux_2 values into rA
        auto p_fe_rA = dynamic_cast<MatrixType*>(&rA);
        if (p_fe_rA) {
            if (rA.isFillActive()) {
                rA.fillComplete();
            }
            p_fe_rA->beginAssembly();
        } else {
            if (!rA.isFillActive()) rA.resumeFill();
        }

        // Preserve existing structure from rD (hard-zero positions) so FE
        // assembly does not shrink the graph to only aux_2 inserted entries.
        for (LO i = 0; i < static_cast<LO>(rD.getNodeNumRows()); ++i) {
            const auto global_row_index = rD.getRowMap()->getGlobalElement(i);
            typename MatrixType::local_inds_host_view_type local_cols_d;
            typename MatrixType::values_host_view_type vals_d;
            rD.getLocalRowView(i, local_cols_d, vals_d);
            if (vals_d.extent(0) > 0) {
                Teuchos::Array<GO> global_cols(local_cols_d.extent(0));
                std::vector<ST> zero_vals(static_cast<std::size_t>(local_cols_d.extent(0)), static_cast<ST>(0));
                for (std::size_t j = 0; j < static_cast<std::size_t>(local_cols_d.extent(0)); ++j) {
                    global_cols[j] = rD.getColMap()->getGlobalElement(local_cols_d(j));
                }
                rA.sumIntoGlobalValues(global_row_index, static_cast<LO>(global_cols.size()), zero_vals.data(), global_cols.data());
            }
        }

        // Zero rA before summing: full replacement, not accumulation.
        rA.setAllToScalar(static_cast<ST>(0));
        for (LO i = 0; i < static_cast<LO>(aux_2->getNodeNumRows()); ++i) {
            const auto global_row_index = aux_2->getRowMap()->getGlobalElement(i);
            typename MatrixType::local_inds_host_view_type local_cols;
            typename MatrixType::values_host_view_type vals;
            aux_2->getLocalRowView(i, local_cols, vals);
            if (vals.extent(0) > 0) {
                Teuchos::Array<GO> global_cols(local_cols.extent(0));
                for (std::size_t j = 0; j < static_cast<std::size_t>(local_cols.extent(0)); ++j) {
                    global_cols[j] = aux_2->getColMap()->getGlobalElement(local_cols(j));
                }
                rA.sumIntoGlobalValues(global_row_index, static_cast<LO>(global_cols.size()), vals.data(), global_cols.data());
            }
        }
        if (p_fe_rA) p_fe_rA->endAssembly();
        if (CallFillCompleteOnResult && rA.isFillActive()) rA.fillComplete();
    }

    /**
     * @brief Returns the multiplication of a vector by a scalar
     * @details y = A*x
     * Checks if a multiplication is needed and tries to do otherwise
     * @param rX The vector considered
     * @param A The scalar considered
     */
    inline static void InplaceMult(
        VectorType& rX,
        const double A
        )
    {
        if (A != 1.00) {
            // Scale the vector: x = A * x
            rX.scale(A);
        }
    }

    /**
     * @brief Returns the multiplication of a vector by a scalar
     * @details x = A*y
     * Checks if a multiplication is needed and tries to do otherwise
     * @note ATTENTION it is assumed no aliasing between rX and rY
     * @param rX The resulting vector considered
     * @param A The scalar considered
     * @param rY The multiplied vector considered
     */
    inline static void Assign(
        VectorType& rX,
        const double A,
        const VectorType& rY
        )
    {
        if (A != 1.00) {
            // Perform the operation x = A * y
            rX.update(A, rY, 0.0);
        } else {
            // FEMultiVector has deleted copy assignment; use update to copy
            rX.update(1.0, rY, 0.0);
        }
    }

    /**
     * @brief Returns the unaliased addition of a vector by a scalar times a vector
     * @details X += A*y;
     * Checks if a multiplication is needed and tries to do otherwise
     * @note ATTENTION it is assumed no aliasing between rX and rY
     * @param rX The resulting vector considered
     * @param A The scalar considered
     * @param rY The multiplied vector considered
     */
    inline static void UnaliasedAdd(
        VectorType& rX,
        const double A,
        const VectorType& rY
        )
    {
        rX.update(A, rY, 1.0);
    }

    /**
     * @brief Returns the unaliased addition of two vectors by a scalar
     * @details rZ = (A * rX) + (B * rY)
     * @param A The scalar considered
     * @param rX The first vector considered
     * @param B The scalar considered
     * @param rY The second vector considered
     * @param rZ The resulting vector considered
     */
    inline static void ScaleAndAdd(
        const double A,
        const VectorType& rX,
        const double B,
        const VectorType& rY,
        VectorType& rZ
        )
    {
        // Compute rZ = A * rX + B * rY
        rZ.update(A, rX, B, rY, 0.0);
    }

    /**
     * @brief Returns the unaliased addition of two vectors by a scalar
     * @details rY = (A * rX) + (B * rY)
     * @param A The scalar considered
     * @param rX The first vector considered
     * @param B The scalar considered
     * @param rY The resulting vector considered
     */
    inline static void ScaleAndAdd(
        const double A,
        const VectorType& rX,
        const double B,
        VectorType& rY
        )
    {
        // Compute rY = A * rX + B * rY
        rY.update(A, rX, B);
    }

    /**
     * @brief Returns the unaliased addition of two matrices by a scalar
     * @details rY = (A * rX) + (B * rY)
     * @param A The scalar considered
     * @param rX The first matrix considered
     * @param B The scalar considered
     * @param rY The resulting matrix considered
     */
    inline static void ScaleAndAdd(
        const double A,
        const MatrixType& rX,
        const double B,
        MatrixType& rY
        )
    {
        // Use a temporary CrsMatrix to allow for dynamic sparsity
        Teuchos::RCP<CrsMatrixType> aux_Y = Teuchos::rcp(new CrsMatrixType(rY.getRowMap(), 100));
        Tpetra::MatrixMatrix::Add(rX, false, A, rY, false, B, aux_Y);
        aux_Y->fillComplete();

        // Inline copy logic from CrsMatrixType to MatrixType
        MatrixType& rDest = rY;
        const CrsMatrixType& rSrc = *aux_Y;
        auto p_fe_Dest = dynamic_cast<MatrixType*>(&rDest);
        if (p_fe_Dest) {
            // For FECrsMatrix: beginAssembly() handles closed->open transition directly.
            // Do NOT call resumeFill() first - it sets isFillActive()=true and then
            // beginAssembly() would crash expecting a closed matrix.
            if (!rDest.isFillActive()) {
                p_fe_Dest->beginAssembly();
            }
        } else if (!rDest.isFillActive()) {
            rDest.resumeFill();
        }
        // Zero out rDest before summing to avoid accumulating onto existing values
        rDest.setAllToScalar(static_cast<ST>(0));
        for (LO i = 0; i < static_cast<LO>(rSrc.getNodeNumRows()); ++i) {
            const auto global_row_index = rSrc.getRowMap()->getGlobalElement(i);
            typename MatrixType::local_inds_host_view_type local_cols;
            typename MatrixType::values_host_view_type vals;
            rSrc.getLocalRowView(i, local_cols, vals);
            if (vals.extent(0) > 0) {
                Teuchos::Array<GO> global_cols(local_cols.extent(0));
                for (std::size_t j = 0; j < static_cast<std::size_t>(local_cols.extent(0)); ++j) {
                    global_cols[j] = rSrc.getColMap()->getGlobalElement(local_cols(j));
                }
                rDest.sumIntoGlobalValues(global_row_index, static_cast<LO>(global_cols.size()), vals.data(), global_cols.data());
            }
        }
        if (p_fe_Dest) {
            p_fe_Dest->endAssembly();
        }
        if (rDest.isFillActive()) rDest.fillComplete();
    }

    /**
     * @brief Sets a value in a vector
     * @param rX The vector considered
     * @param i The index of the value considered
     * @param value The value considered
     */
    inline static void SetValue(
        VectorType& rX,
        IndexType i,
        const double value
        )
    {
        // Get the local index corresponding to the global index `i`
        auto map = rX.getMap();
        IndexType localIndex = map->getLocalElement(i);

        // Check if the index is on this process
        if (localIndex != Tpetra::Details::OrdinalTraits<IndexType>::invalid()) {
            // Set the value at the specified local index
            rX.replaceLocalValue(localIndex, size_t(0), value);
        }
        // If the index `i` is not local, it is ignored (handled by Tpetra's parallel distribution)
    }

    /**
     * @brief Sets a value in a matrix
     * @param rA The matrix considered
     * @param i The row index
     * @param j The column index
     * @param value The value considered
     */
    inline static void SetValue(
        MatrixType& rA,
        IndexType i,
        IndexType j,
        const double value
        )
    {
        const bool was_closed = !rA.isFillActive();
        if (was_closed) {
            rA.beginAssembly();
        }
        const GO globalRow = static_cast<GO>(i);
        const GO globalCol = static_cast<GO>(j);
        const ST val = static_cast<ST>(value);
        rA.replaceGlobalValues(globalRow, 1, &val, &globalCol);
        if (was_closed) {
            rA.endAssembly();
            if (rA.isFillActive()) rA.fillComplete();
        }
    }

    /**
     * @brief assigns a scalar to a vector
     * @details rX = A
     * @param rX The vector considered
     * @param A The scalar considered
     */
    inline static void Set(
        VectorType& rX,
        const DataType A
        )
    {
        rX.putScalar(static_cast<ST>(A));
    }

    /**
     * @brief Resizes a matrix
     * @param rA The matrix to be resized
     * @param m The new number of rows
     * @param n The new number of columns
     */
    inline static void Resize(
        MatrixType& rA,
        const SizeType m,
        const SizeType n
        )
    {
        KRATOS_ERROR << "Resize is not defined for Trilinos Sparse Matrix" << std::endl;
    }

    /**
     * @brief Resizes a vector
     * @param rX The vector to be resized
     * @param n The new size
     */
    inline static void Resize(
        VectorType& rX,
        const SizeType n
        )
    {
        KRATOS_ERROR << "Resize is not defined for a reference to Trilinos Vector - need to use the version passing a Pointer" << std::endl;
    }

    /**
     * @brief Resizes a vector
     * @param pX The pointer to the vector to be resized
     * @param n The new global size
    */
    inline static void Resize(
        VectorPointerType& pX,
        const SizeType n
        )
    {
        auto pComm = pX->getMap()->getComm();
        Tpetra::global_size_t global_elems = static_cast<Tpetra::global_size_t>(n);
        MapPointerType new_map = Teuchos::rcp(new MapType(global_elems, 0, pComm));
        VectorPointerType pNewEmptyX = CreateVector(new_map);
        pX.swap(pNewEmptyX);
    }

    /**
     * @brief Clears a matrix
     * @param pA The pointer to the matrix to be cleared
     */
    inline static void Clear(MatrixPointerType& pA)
    {
        if(pA != Teuchos::null) {
            auto map = Teuchos::rcp(new MapType(0, 0, pA->getMap()->getComm()));
            Teuchos::RCP<GraphType> graph = Teuchos::rcp(new GraphType(map, map, 0));
            // Use fillComplete() instead of endAssembly() to avoid FE state machine issues
            if (graph->isFillActive()) {
                graph->fillComplete();
            }
            MatrixPointerType pNewEmptyA = Teuchos::rcp(new MatrixType(
                Teuchos::rcp_const_cast<const GraphType>(graph)));
            pA.swap(pNewEmptyA);
        }
    }

    /**
     * @brief Clears a vector
     * @param pX The pointer to the vector to be cleared
     */
    inline static void Clear(VectorPointerType& pX)
    {
        if(pX != Teuchos::null) {
            auto map = Teuchos::rcp(new MapType(0, 0, pX->getMap()->getComm()));
            VectorPointerType pNewEmptyX = CreateVector(map);
            pX.swap(pNewEmptyX);
        }
    }

    /**
     * @brief Sets a matrix to zero
     * @param rX The matrix to be set
     */
    inline static void SetToZero(MatrixType& rA)
    {
        // Open for FE assembly if the matrix is currently closed (fillComplete was called).
        // Avoid mixing resumeFill() + beginAssembly() which conflicts in FECrsMatrix.
        const bool was_open = rA.isFillActive();
        if (!was_open) {
            rA.beginAssembly();
        }
        rA.setAllToScalar(static_cast<ST>(0));
        // If we opened it ourselves, close it back so callers see a consistent closed state.
        if (!was_open && rA.isFillActive()) {
            rA.endAssembly();
        }
    }

    /**
     * @brief Sets a vector to zero
     * @param rX The vector to be set
     */
    inline static void SetToZero(VectorType& rX)
    {
        rX.putScalar(static_cast<ST>(0));
    }

    /**
     * @brief Build Tpetra FECrsGraph and create new system matrix + vectors.
     * @details Uses a proper OWNED_PLUS_SHARED map that includes ghost DOFs
     *          (equation IDs referenced by local elements but owned by other
     *          processes). This enables FECrsMatrix to accept cross-partition
     *          element contributions and communicate them in endAssembly().
     * @param rComm The MPI communicator
     * @param LocalSize The number of locally owned DOFs on this rank
     * @param FirstMyId The first global DOF id owned by this rank
     * @param GuessRowSize The estimated number of non-zeros per row
     * @param rAllEquationIds The list of lists of equation ids for each element
     * @param rpA The pointer to the system matrix to be created
     * @param rpb The pointer to the RHS vector to be created
     * @param rpDx The pointer to the solution vector to be created
     * @param rpReactions The pointer to the reactions vector to be created
     * @param equationSystemSize The global size of the equation system
     * @param pMap The map to be used for the construction of the matrix and vectors
     */
    static void BuildSystemStructure(
        CommunicatorType& rComm,
        const IndexType LocalSize,
        const int FirstMyId,
        const int GuessRowSize,
        const std::vector<std::vector<int>>& rAllEquationIds,
        MatrixPointerType& rpA,
        VectorPointerType& rpb,
        VectorPointerType& rpDx,
        VectorPointerType& rpReactions,
        const IndexType equationSystemSize,
        MapPointerType pMap)
    {
        auto comm = pMap->getComm();
        const int nproc = static_cast<int>(comm->getSize());

        const GO firstId = static_cast<GO>(FirstMyId);
        const GO lastId  = firstId + static_cast<GO>(LocalSize);

        // Step 1: collect locally-visible ghost DOF GIDs (referenced by local elements
        // but owned by other processes, i.e. outside [firstId, lastId)).
        std::set<GO> ghost_gids_set;
        for (const auto& eq_ids : rAllEquationIds) {
            for (int id : eq_ids) {
                const GO gid = static_cast<GO>(id);
                if (gid < firstId || gid >= lastId)
                    ghost_gids_set.insert(gid);
            }
        }

        // Step 2: Symmetric ghost expansion.
        //
        // FECrsGraph::endFill() has two code paths controlled solely by whether
        // ownedRowsImporter_ is null.  The importer is set to null whenever the
        // ownedMap and ownedPlusSharedMap are identical (maps_are_the_same==true).
        // The "easy case"  (null importer) skips doExport(); the "hard case" calls it.
        // When some processes take the easy case and others the hard case, their MPI
        // operations diverge → MPI_ERR_TRUNCATE.
        //
        // Fix: every process that will *receive* ghost contributions (i.e. some other
        // process is sending rows that fall in my owned range) must also have at least
        // one ghost DOF so it enters the hard case.  We achieve this with an AllGather
        // of each process's owned-DOF range, followed by a notification round-trip that
        // tells each receiver to add a dummy ghost from the sender's range.

        // Gather every process's [firstId, localSize] so we can find ownership.
        GO myFid = firstId, mySz = static_cast<GO>(LocalSize);
        Teuchos::Array<GO> allFirstIds(nproc), allLocalSizes(nproc);
        Teuchos::gatherAll(*comm, 1, &myFid, nproc, allFirstIds.getRawPtr());
        Teuchos::gatherAll(*comm, 1, &mySz,  nproc, allLocalSizes.getRawPtr());

        // For each ghost GID we own locally, identify its owner process and record
        // that we will send that process a contribution.  We notify the owner by
        // storing our firstId in notify_to[owner].
        Teuchos::Array<GO> notify_to(nproc, static_cast<GO>(-1));
        for (GO g : ghost_gids_set) {
            for (int p = 0; p < nproc; ++p) {
                if (g >= allFirstIds[p] && g < allFirstIds[p] + allLocalSizes[p]) {
                    if (p != static_cast<int>(comm->getRank()))
                        notify_to[p] = firstId;
                    break;
                }
            }
        }

        // AllGather the full nproc×nproc notification matrix.
        // notify_all[q*nproc + k] = value that process q put in notify_to[k].
        Teuchos::Array<GO> notify_all(nproc * nproc, static_cast<GO>(-1));
        Teuchos::gatherAll(*comm, nproc, notify_to.getRawPtr(),
                           nproc * nproc, notify_all.getRawPtr());

        // Each process k checks column k of notify_all.  If process q will send it
        // data (notify_all[q*nproc+k] != -1), add q's firstId as a dummy ghost so
        // process k also enters the hard case.
        const int myrank = static_cast<int>(comm->getRank());
        for (int q = 0; q < nproc; ++q) {
            if (notify_all[q * nproc + myrank] != static_cast<GO>(-1)) {
                const GO q_first = allFirstIds[q];
                if (q_first < firstId || q_first >= lastId)   // not already owned
                    ghost_gids_set.insert(q_first);
            }
        }

        // Force ALL ranks into the FECrsGraph "hard case" if any rank has ghosts.
        // Idle ranks (e.g. LocalSize==0, no elements) would otherwise take the
        // "easy path" in endAssembly(), which skips the MPI collectives that
        // active ranks call, causing a collective mismatch and crash.
        {
            GO has_ghosts_local = ghost_gids_set.empty() ? GO(0) : GO(1);
            GO has_ghosts_global = GO(0);
            Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &has_ghosts_local, &has_ghosts_global);
            if (has_ghosts_global > GO(0) && ghost_gids_set.empty()) {
                // Add a dummy ghost from the first rank that owns any DOFs.
                for (int p = 0; p < nproc; ++p) {
                    const GO q_first = allFirstIds[p];
                    if (allLocalSizes[p] > 0 && (q_first < firstId || q_first >= lastId)) {
                        ghost_gids_set.insert(q_first);
                        break;
                    }
                }
            }
        }

        // Step 3: Build OWNED_PLUS_SHARED map.
        // After symmetric expansion every process with any cross-partition elements
        // has a non-trivially-different OPS map → all enter FECrsGraph's hard case
        // → all call doExport during endAssembly → no MPI divergence.
        MapPointerType owned_plus_shared_map;
        if (ghost_gids_set.empty()) {
            // Purely serial or identical partition — safe easy case for ALL processes.
            owned_plus_shared_map = pMap;
        } else {
            // Owned GIDs first (required by FECrsGraph::setup locality check),
            // then ghost GIDs.
            std::vector<GO> ops_gids;
            ops_gids.reserve(LocalSize + ghost_gids_set.size());
            for (IndexType i = 0; i < LocalSize; ++i)
                ops_gids.push_back(firstId + static_cast<GO>(i));
            for (GO g : ghost_gids_set)
                ops_gids.push_back(g);
            owned_plus_shared_map = Teuchos::rcp(new MapType(
                Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                Teuchos::ArrayView<const GO>(ops_gids.data(),
                                            static_cast<int>(ops_gids.size())),
                0, pMap->getComm()));
        }

        // Step 4: Build FECrsGraph with (ownedMap, ownedPlusSharedMap).
        Teuchos::RCP<GraphType> graph =
            Teuchos::rcp(new GraphType(pMap, owned_plus_shared_map, GuessRowSize));
        graph->beginAssembly();
        std::vector<GO> gids;
        for (const auto& eq_ids : rAllEquationIds) {
            if (eq_ids.empty()) continue;
            gids.resize(eq_ids.size());
            for (std::size_t k = 0; k < eq_ids.size(); ++k)
                gids[k] = static_cast<GO>(eq_ids[k]);
            for (std::size_t row = 0; row < gids.size(); ++row) {
                if (owned_plus_shared_map->getLocalElement(gids[row]) !=
                        Teuchos::OrdinalTraits<LO>::invalid()) {
                    graph->insertGlobalIndices(
                        gids[row],
                        Teuchos::ArrayView<const GO>(gids.data(),
                                                     static_cast<int>(gids.size())));
                }
            }
        }
        // endAssembly: communicates ghost-row structure to owning processes and
        // switches the active graph from OWNED_PLUS_SHARED to OWNED.
        graph->endAssembly();

        rpA = Teuchos::rcp(new MatrixType(
            Teuchos::rcp_const_cast<const GraphType>(graph)));
        // Post-construction normalization: FECrsMatrix ctor leaves isFillActive()=true
        // but fillState_=closed (internal beginFill→resumeFill).  Call fillComplete()
        // to reset so GlobalAssemble does not incorrectly call endAssembly().
        if (rpA->isFillActive()) rpA->fillComplete();

        if (!rpb || Size(*rpb) != equationSystemSize)
            rpb = CreateVector(pMap);
        if (!rpDx || Size(*rpDx) != equationSystemSize)
            rpDx = CreateVector(pMap);
        if (!rpReactions)
            rpReactions = CreateVector(pMap);
    }

    /**
     * @brief Build Tpetra FECrsGraph with separate row and column block lists.
     * @details This overload mirrors the Epetra BuildSystemStructure with separate
     *          rAllRowEquationIds and rAllColEquationIds block lists, enabling
     *          rectangular (non-symmetric) sparsity patterns.
     * @param rComm The MPI communicator
     * @param LocalSize The local size of the system
     * @param FirstMyId The first global id owned by this rank
     * @param GuessRowSize The guess row size for the graph construction
     * @param rAllRowEquationIds The list of lists of row equation ids for each element
     * @param rAllColEquationIds The list of lists of column equation ids for each element
     * @param rpA The pointer to the system matrix to be created
     * @param rpb The pointer to the RHS vector to be created
     * @param rpDx The pointer to the solution vector to be created
     * @param rpReactions The pointer to the reactions vector to be created
     * @param EquationSystemSize The global size of the equation system
     * @param pMap The map to be used for the construction of the matrix and vectors
     */
    static void BuildSystemStructure(
        CommunicatorType& rComm,
        const IndexType LocalSize,
        const int FirstMyId,
        const int GuessRowSize,
        const std::vector<std::vector<int>>& rAllRowEquationIds,
        const std::vector<std::vector<int>>& rAllColEquationIds,
        MatrixPointerType& rpA,
        VectorPointerType& rpb,
        VectorPointerType& rpDx,
        VectorPointerType& rpReactions,
        const IndexType EquationSystemSize,
        MapPointerType pMap)
    {
        KRATOS_ERROR_IF(rAllRowEquationIds.size() != rAllColEquationIds.size())
            << "BuildSystemStructure: row and column block lists must have the same size" << std::endl;

        auto comm = pMap->getComm();
        const int nproc = static_cast<int>(comm->getSize());
        const GO firstId = static_cast<GO>(FirstMyId);
        const GO lastId  = firstId + static_cast<GO>(LocalSize);

        // Collect ghost GIDs from both row and column ids
        std::set<GO> ghost_gids_set;
        for (std::size_t i_block = 0; i_block < rAllRowEquationIds.size(); ++i_block) {
            for (int id : rAllRowEquationIds[i_block]) {
                const GO gid = static_cast<GO>(id);
                if (gid < firstId || gid >= lastId) ghost_gids_set.insert(gid);
            }
            for (int id : rAllColEquationIds[i_block]) {
                const GO gid = static_cast<GO>(id);
                if (gid < firstId || gid >= lastId) ghost_gids_set.insert(gid);
            }
        }

        // Symmetric ghost expansion (same logic as the single-block overload)
        GO myFid = firstId, mySz = static_cast<GO>(LocalSize);
        Teuchos::Array<GO> allFirstIds(nproc), allLocalSizes(nproc);
        Teuchos::gatherAll(*comm, 1, &myFid, nproc, allFirstIds.getRawPtr());
        Teuchos::gatherAll(*comm, 1, &mySz,  nproc, allLocalSizes.getRawPtr());

        Teuchos::Array<GO> notify_to(nproc, static_cast<GO>(-1));
        for (GO g : ghost_gids_set) {
            for (int p = 0; p < nproc; ++p) {
                if (g >= allFirstIds[p] && g < allFirstIds[p] + allLocalSizes[p]) {
                    if (p != static_cast<int>(comm->getRank()))
                        notify_to[p] = firstId;
                    break;
                }
            }
        }

        Teuchos::Array<GO> notify_all(nproc * nproc, static_cast<GO>(-1));
        Teuchos::gatherAll(*comm, nproc, notify_to.getRawPtr(),
                           nproc * nproc, notify_all.getRawPtr());

        const int myrank = static_cast<int>(comm->getRank());
        for (int q = 0; q < nproc; ++q) {
            if (notify_all[q * nproc + myrank] != static_cast<GO>(-1)) {
                const GO q_first = allFirstIds[q];
                if (q_first < firstId || q_first >= lastId)
                    ghost_gids_set.insert(q_first);
            }
        }

        // Force ALL ranks into the FECrsGraph "hard case" if any rank has ghosts.
        // (same logic as the single-block overload above)
        {
            GO has_ghosts_local = ghost_gids_set.empty() ? GO(0) : GO(1);
            GO has_ghosts_global = GO(0);
            Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &has_ghosts_local, &has_ghosts_global);
            if (has_ghosts_global > GO(0) && ghost_gids_set.empty()) {
                for (int p = 0; p < nproc; ++p) {
                    const GO q_first = allFirstIds[p];
                    if (allLocalSizes[p] > 0 && (q_first < firstId || q_first >= lastId)) {
                        ghost_gids_set.insert(q_first);
                        break;
                    }
                }
            }
        }

        // Build OWNED_PLUS_SHARED map
        MapPointerType owned_plus_shared_map;
        if (ghost_gids_set.empty()) {
            owned_plus_shared_map = pMap;
        } else {
            std::vector<GO> ops_gids;
            ops_gids.reserve(LocalSize + ghost_gids_set.size());
            for (IndexType i = 0; i < LocalSize; ++i)
                ops_gids.push_back(firstId + static_cast<GO>(i));
            for (GO g : ghost_gids_set)
                ops_gids.push_back(g);
            owned_plus_shared_map = Teuchos::rcp(new MapType(
                Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                Teuchos::ArrayView<const GO>(ops_gids.data(),
                                            static_cast<int>(ops_gids.size())),
                0, pMap->getComm()));
        }

        // Build FECrsGraph inserting row_ids × col_ids blocks
        Teuchos::RCP<GraphType> graph =
            Teuchos::rcp(new GraphType(pMap, owned_plus_shared_map, GuessRowSize));
        graph->beginAssembly();
        for (std::size_t i_block = 0; i_block < rAllRowEquationIds.size(); ++i_block) {
            const auto& r_row_ids = rAllRowEquationIds[i_block];
            const auto& r_col_ids = rAllColEquationIds[i_block];
            if (r_row_ids.empty() || r_col_ids.empty()) continue;
            std::vector<GO> row_gids(r_row_ids.size());
            std::vector<GO> col_gids(r_col_ids.size());
            for (std::size_t k = 0; k < r_row_ids.size(); ++k)
                row_gids[k] = static_cast<GO>(r_row_ids[k]);
            for (std::size_t k = 0; k < r_col_ids.size(); ++k)
                col_gids[k] = static_cast<GO>(r_col_ids[k]);
            for (const GO row : row_gids) {
                if (owned_plus_shared_map->getLocalElement(row) !=
                        Teuchos::OrdinalTraits<LO>::invalid()) {
                    graph->insertGlobalIndices(
                        row,
                        Teuchos::ArrayView<const GO>(col_gids.data(),
                                                     static_cast<int>(col_gids.size())));
                }
            }
        }
        graph->endAssembly();

        rpA = Teuchos::rcp(new MatrixType(
            Teuchos::rcp_const_cast<const GraphType>(graph)));
        if (rpA->isFillActive()) rpA->fillComplete();

        if (!rpb || Size(*rpb) != EquationSystemSize)
            rpb = CreateVector(pMap);
        if (!rpDx || Size(*rpDx) != EquationSystemSize)
            rpDx = CreateVector(pMap);
        if (!rpReactions)
            rpReactions = CreateVector(pMap);
    }

    /**
     * @brief Build a Tpetra FE constraint graph and create T matrix + constant vector.
     * @param rComm The communicator considered
     * @param LocalSize The local size of the system
     * @param FirstMyId The first global id owned by this rank
     * @param GuessRowSize The guess row size for the graph construction
     * @param rSlaveEquationIds The list of lists of slave equation ids for each local row
     * @param rMasterEquationIds The list of lists of master equation ids for each local row
     * @param rpT The pointer to the T matrix to be created
     * @param rpConstantVector The pointer to the constant vector to be created
     * @param pMap The map to be used for the construction of the matrix and vectors
     */
    static void BuildConstraintsStructure(
        CommunicatorType& rComm,
        const IndexType LocalSize,
        const int FirstMyId,
        const int GuessRowSize,
        const std::vector<std::vector<int>>& rSlaveEquationIds,
        const std::vector<std::vector<int>>& rMasterEquationIds,
        MatrixPointerType& rpT,
        VectorPointerType& rpConstantVector,
        MapPointerType pMap
        )
    {
        Teuchos::RCP<GraphType> graph = Teuchos::rcp(new GraphType(pMap, pMap, GuessRowSize));
        for (IndexType i = 0; i < LocalSize; ++i) {
            const GO gid = static_cast<GO>(FirstMyId + static_cast<int>(i));
            graph->insertGlobalIndices(gid, Teuchos::ArrayView<const GO>(&gid, 1));
        }
        for (std::size_t c = 0; c < rSlaveEquationIds.size(); ++c) {
            const auto& slave_ids = rSlaveEquationIds[c];
            const auto& master_ids = rMasterEquationIds[c];
            if (slave_ids.empty() || master_ids.empty()) continue;
            std::vector<GO> master_gids(master_ids.size());
            for (std::size_t k = 0; k < master_ids.size(); ++k) master_gids[k] = static_cast<GO>(master_ids[k]);
            for (int slave_id : slave_ids) {
                const GO slave_gid = static_cast<GO>(slave_id);
                graph->insertGlobalIndices(slave_gid,
                    Teuchos::ArrayView<const GO>(master_gids.data(), static_cast<int>(master_gids.size())));
            }
        }
        if (graph->isFillActive()) graph->fillComplete();
        rpT = Teuchos::rcp(new MatrixType(Teuchos::rcp_const_cast<const GraphType>(graph)));
        // Post-construction normalization (same as BuildSystemStructure above).
        if (rpT->isFillActive()) rpT->fillComplete();
        rpConstantVector = CreateVector(pMap);
    }

    /**
     * @brief Assembles the LHS of the system
     * @param rA The LHS matrix
     * @param rLHSContribution The contribution to the LHS
     * @param rEquationId The equation ids
     */
    inline static void AssembleLHS(
        MatrixType& rA,
        const Matrix& rLHSContribution,
        const std::vector<std::size_t>& rEquationId
        )
    {
        // Lazily open for FE assembly on first call in a Build() cycle.
        // GlobalAssemble() will close it at the end.
        auto p_fe_rA = dynamic_cast<MatrixType*>(&rA);
        if (p_fe_rA && !rA.isFillActive()) {
            p_fe_rA->beginAssembly();
        }

        const std::size_t system_size = rA.getGlobalNumRows();

        // Count active indices
        std::vector<LO> indices;
        for (std::size_t i = 0; i < rEquationId.size(); ++i) {
            if (rEquationId[i] < system_size) {
                indices.push_back(static_cast<LO>(rEquationId[i]));
            }
        }

        if (!indices.empty()) {
            std::vector<GO> global_indices(indices.size());
            for (std::size_t i = 0; i < indices.size(); ++i) {
                global_indices[i] = static_cast<GO>(indices[i]);
            }

            for (std::size_t i = 0; i < indices.size(); ++i) {
                const GO globalRow = global_indices[i];
                std::vector<ST> row_values(indices.size());
                for (std::size_t j = 0; j < indices.size(); ++j) {
                    row_values[j] = rLHSContribution(i, j);
                }
                const int ierr = rA.sumIntoGlobalValues(globalRow, static_cast<LO>(global_indices.size()), row_values.data(), global_indices.data());
                // Note: sumIntoGlobalValues might return the number of values successfully summed instead of an error code 0 or -1. 
                // Epetra returns 0, Tpetra returns the number of values (indices.size()) if successful.
                KRATOS_ERROR_IF(ierr != static_cast<int>(indices.size())) << "Tpetra failure found" << std::endl;
            }
        }
    }

    /**
     * @brief Assembles the RHS of the system
     * @param rb The RHS vector
     * @param rRHSContribution The RHS contribution
     * @param rEquationId The equation ids
     */
    inline static void AssembleRHS(
        VectorType& rb,
        const Vector& rRHSContribution,
        const std::vector<std::size_t>& rEquationId
        )
    {
        const std::size_t system_size = rb.getGlobalLength();

        // Count active indices
        std::vector<LO> indices;
        for (std::size_t i = 0; i < rEquationId.size(); ++i) {
            if (rEquationId[i] < system_size) {
                indices.push_back(static_cast<LO>(rEquationId[i]));
            }
        }

        if (!indices.empty()) {
            std::vector<GO> global_indices(indices.size());
            std::vector<ST> values(indices.size());
            for (std::size_t i = 0; i < indices.size(); ++i) {
                global_indices[i] = static_cast<GO>(indices[i]);
                values[i] = rRHSContribution[i];
            }
            for (std::size_t i = 0; i < global_indices.size(); ++i) {
                rb.sumIntoGlobalValue(global_indices[i], size_t(0), values[i]);
            }
        }
    }

    /**
     * @brief This function returns if we are in a distributed system
     * @return True if we are in a distributed system, false otherwise (always true in this case)
     */
    inline static constexpr bool IsDistributed()
    {
        return true;
    }

   /**
    * @brief Check if the TrilinosSpaceExperimental is distributed.
    * @details This static member function checks whether the TrilinosSpaceExperimental is distributed or not.
    * @return True if the space is distributed, false otherwise.
    */
    static constexpr bool IsDistributedSpace()
    {
        return true;
    }

    /**
     * @brief Returns a list of the fastest direct solvers.
     * @details This function returns a vector of strings representing the names of the fastest direct solvers. The order of the solvers in the list may need to be updated and reordered depending on the size of the equation system.
     * @return A vector of strings containing the names of the fastest direct solvers.
     */
    inline static std::vector<std::string> FastestDirectSolverList()
    {
        // May need to be updated and reordered. In fact I think it depends of the size of the equation system
        std::vector<std::string> faster_direct_solvers({
            "mumps2",         // Amesos2 (if compiled with MUMPS-support)
            "mumps",          // Amesos (if compiled with MUMPS-support)
            "super_lu_dist2", // Amesos2 SuperLUDist (if compiled with MPI-support)
            "super_lu_dist",  // Amesos SuperLUDist (if compiled with MPI-support)
            "amesos2",        // Amesos2
            "amesos",         // Amesos
            "klu2",           // Amesos2 KLU
            "klu",            // Amesos KLU
            "basker"          // Amesos2 Basker
        });
        return faster_direct_solvers;
    }

    /**
     * @brief Sets a vector entry using a global index and performs assembly.
     * @param rX The vector to be modified.
     * @param i The global index of the entry.
     * @param Value The value to assign.
     * @details This method sets the value at global index i (if locally owned)
     * and then calls GlobalAssemble() so the change is synchronized.
     */
    static void SetGlobalVec(
        VectorType& rX,
        const IndexType i,
        const double Value
        )
    {
        auto map = rX.getMap();
        IndexType localIndex = map->getLocalElement(i);
        if (localIndex != Tpetra::Details::OrdinalTraits<IndexType>::invalid()) {
            rX.replaceLocalValue(localIndex, size_t(0), static_cast<ST>(Value));
        }
        GlobalAssemble(rX);
    }

    /**
     * @brief Sets a vector entry using a global index without assembly.
     * @param rX The vector to be modified.
     * @param i The global index of the entry.
     * @param Value The value to assign.
     * @details This method is intended for batched updates where assembly is
     * deferred and handled separately by the caller.
     */
    static void SetGlobalVecNoAssemble(
        VectorType& rX,
        const IndexType i,
        const double Value
        )
    {
        auto map = rX.getMap();
        IndexType localIndex = map->getLocalElement(i);
        if (localIndex != Tpetra::Details::OrdinalTraits<IndexType>::invalid()) {
            rX.replaceLocalValue(localIndex, size_t(0), static_cast<ST>(Value));
        }
    }

    /**
     * @brief Sets a vector entry using a local index and performs assembly.
     * @param rX The vector to be modified.
     * @param i The local index of the entry on the current rank.
     * @param Value The value to assign.
     * @details This method uses local indexing and then calls GlobalAssemble()
     * to propagate pending contributions.
     */
    static void SetLocalVec(
        VectorType& rX,
        const IndexType i,
        const double Value
        )
    {
        rX.replaceLocalValue(static_cast<LO>(i), size_t(0), static_cast<ST>(Value));
        GlobalAssemble(rX);
    }

    /**
     * @brief Sets a vector entry using a local index without assembly.
     * @param rX The vector to be modified.
     * @param i The local index of the entry on the current rank.
     * @param Value The value to assign.
     * @details This method is useful when several local values are modified
     * and a single later assembly is preferred for efficiency.
     */
    static void SetLocalVecNoAssemble(
        VectorType& rX,
        const IndexType i,
        const double Value
        )
    {
        rX.replaceLocalValue(static_cast<LO>(i), size_t(0), static_cast<ST>(Value));
    }

    /**
     * @brief Sets a matrix entry using global row and column indices and performs assembly.
     * @param rA The matrix to be modified.
     * @param i The global row index.
     * @param j The global column index.
     * @param Value The value to assign.
     * @details The value is inserted with global indexing and the matrix is then
     * globally assembled so the modification becomes consistent in parallel.
     */
    static void SetGlobalMat(
        MatrixType& rA,
        IndexType i,
        IndexType j,
        double Value
        )
    {
        const bool was_closed = !rA.isFillActive();
        if (was_closed) {
            rA.beginAssembly();
        }
        const GO globalRow = static_cast<GO>(i);
        const GO globalCol = static_cast<GO>(j);
        const ST val = static_cast<ST>(Value);
        rA.replaceGlobalValues(globalRow, 1, &val, &globalCol);
        if (was_closed) {
            rA.endAssembly();
            if (rA.isFillActive()) rA.fillComplete();
        }
    }

    /**
     * @brief Sets a matrix entry using global row and column indices without assembly.
     * @param rA The matrix to be modified.
     * @param i The global row index.
     * @param j The global column index.
     * @param Value The value to assign.
     * @details This method is intended for batched matrix insertions where
     * GlobalAssemble() will be invoked separately by the caller.
     */
    static void SetGlobalMatNoAssemble(
        MatrixType& rA,
        const IndexType i,
        const IndexType j,
        const double Value
        )
    {
        if (!rA.isFillActive()) {
            rA.beginAssembly();
        }
        const GO globalRow = static_cast<GO>(i);
        const GO globalCol = static_cast<GO>(j);
        const ST val = static_cast<ST>(Value);
        rA.replaceGlobalValues(globalRow, 1, &val, &globalCol);
        // Assembly must be finalized by the caller via GlobalAssemble()
    }

    /**
     * @brief Sets a matrix entry using local row and column indices and performs assembly.
     * @param rA The matrix to be modified.
     * @param i The local row index on the current rank.
     * @param j The local column index on the current rank.
     * @param Value The value to assign.
     * @details This method uses local indexing and then assembles the FE matrix
     * so pending contributions are synchronized across ranks.
     */
    static void SetLocalMat(
        MatrixType& rA,
        const IndexType i,
        const IndexType j,
        const double Value
        )
    {
        const bool was_closed = !rA.isFillActive();
        if (was_closed) {
            rA.beginAssembly();
        }
        const LO local_row = static_cast<LO>(i);
        const LO local_col = static_cast<LO>(j);
        const ST val = static_cast<ST>(Value);
        rA.replaceLocalValues(local_row, 1, &val, &local_col);
        if (was_closed) {
            rA.endAssembly();
            if (rA.isFillActive()) rA.fillComplete();
        }
    }

    /**
     * @brief Sets a matrix entry using local row and column indices without assembly.
     * @param rA The matrix to be modified.
     * @param i The local row index on the current rank.
     * @param j The local column index on the current rank.
     * @param Value The value to assign.
     * @details This method is useful for repeated local updates followed by a
     * single explicit GlobalAssemble() call.
     */
    static void SetLocalMatNoAssemble(
        MatrixType& rA,
        const IndexType i,
        const IndexType j,
        const double Value
        )
    {
        if (!rA.isFillActive()) {
            rA.beginAssembly();
        }
        const LO local_row = static_cast<LO>(i);
        const LO local_col = static_cast<LO>(j);
        const ST val = static_cast<ST>(Value);
        rA.replaceLocalValues(local_row, 1, &val, &local_col);
        // Assembly must be finalized by the caller via GlobalAssemble()
    }

    /**
     * @brief This function returns a value from a given vector according to a given index
     * @param rX The vector from which values are to be gathered
     * @param I The index of the value to be gathered
     * @return The value of the vector corresponding to the index I
     */
    inline static double GetValue(
        const VectorType& rX,
        const std::size_t I
        )
    {
        // Get the local index corresponding to the global index `I`
        auto map = rX.getMap();
        IndexType localIndex = map->getLocalElement(static_cast<IndexType>(I));

        // Index must be local to this proc
        KRATOS_ERROR_IF(localIndex == Tpetra::Details::OrdinalTraits<IndexType>::invalid()) << " non-local id: " << I << "." << std::endl;

        // Get the value at the specified local index via local view (FEMultiVector-compatible)
        auto localView = rX.getLocalViewHost(Tpetra::Access::ReadOnly);
        return static_cast<double>(localView(localIndex, 0));
    }

    /**
     * @brief This function gathers the values of a given vector according to a given index array
     * @param rX The vector from which values are to be gathered
     * @param IndexArray The array containing the indices of the values to be gathered
     * @param pValues The array containing the gathered values
     */
    inline static void GatherValues(
        const VectorType& rX,
        const std::vector<int>& IndexArray,
        double* pValues
        )
    {
        KRATOS_TRY

        // Get the total size of the index array
        const std::size_t tot_size = IndexArray.size();

        // Convert to the correct GO type for Tpetra
        std::vector<GO> global_indices(tot_size);
        for (std::size_t i = 0; i < tot_size; ++i) {
            global_indices[i] = static_cast<GO>(IndexArray[i]);
        }

        // Create a Map with the desired indices using the correct Tpetra types
        auto pComm = rX.getMap()->getComm();
        auto dof_update_map = Tpetra::createNonContigMapWithNode<LO, GO, NT>(
            Teuchos::ArrayView<const GO>(global_indices.data(), tot_size), pComm);

        // Define importer: source=rX.getMap(), target=dof_update_map
        Tpetra::Import<LO, GO, NT> importer(rX.getMap(), dof_update_map);

        // Use a plain MultiVector (FEMultiVector has no single-map constructor)
        Tpetra::MultiVector<ST, LO, GO, NT> temp(dof_update_map, 1);

        // Import the values from rX into the temp vector
        temp.doImport(rX, importer, Tpetra::INSERT);

        // Extract values via local view
        auto localView = temp.getLocalViewHost(Tpetra::Access::ReadOnly);
        for (std::size_t i = 0; i < tot_size; ++i) {
            pValues[i] = static_cast<double>(localView(i, 0));
        }

        // Synchronize processes
        pComm->barrier();

        KRATOS_CATCH("")
    }

    /**
     * @brief Generates a graph combining the graphs of two matrices
     * @param rA The first matrix
     * @param rB The second matrix
     * @return A pointer to the combined FECrsGraph
     */
    static GraphPointerType CombineMatricesGraphs(
        const MatrixType& rA,
        const MatrixType& rB
        )
    {
        // Row maps must be the same
        KRATOS_ERROR_IF(!rA.getRowMap()->isSameAs(*rB.getRowMap())) << "Row maps are not compatible" << std::endl;

        // Getting the graphs
        auto p_graph_a = rA.getCrsGraph();
        auto p_graph_b = rB.getCrsGraph();

        // Getting the maps
        const auto& r_row_map = rA.getRowMap();

        // New graph with large capacity
        Teuchos::RCP<GraphType> graph = Teuchos::rcp(new GraphType(r_row_map, r_row_map, 100));
        // Open the FE graph for insertion (sets fillState_=open)
        graph->beginAssembly();

        const auto numLocalRows = r_row_map->getNodeNumElements();

        // Combine graphs using global indexing
        for (LO i = 0; i < static_cast<LO>(numLocalRows); ++i) {
            const auto global_row_index = r_row_map->getGlobalElement(i);
            std::set<GO> combined_indices;

            if (p_graph_a->isLocallyIndexed()) {
                typename MatrixType::local_inds_host_view_type cols_a;
                p_graph_a->getLocalRowView(i, cols_a);
                for (std::size_t j = 0; j < static_cast<std::size_t>(cols_a.extent(0)); ++j) combined_indices.insert(p_graph_a->getColMap()->getGlobalElement(cols_a(j)));
            } else {
                typename MatrixType::global_inds_host_view_type cols_a;
                p_graph_a->getGlobalRowView(global_row_index, cols_a);
                for (std::size_t j = 0; j < static_cast<std::size_t>(cols_a.extent(0)); ++j) combined_indices.insert(cols_a(j));
            }

            if (p_graph_b->isLocallyIndexed()) {
                typename MatrixType::local_inds_host_view_type cols_b;
                p_graph_b->getLocalRowView(i, cols_b);
                for (std::size_t j = 0; j < static_cast<std::size_t>(cols_b.extent(0)); ++j) combined_indices.insert(p_graph_b->getColMap()->getGlobalElement(cols_b(j)));
            } else {
                typename MatrixType::global_inds_host_view_type cols_b;
                p_graph_b->getGlobalRowView(global_row_index, cols_b);
                for (std::size_t j = 0; j < static_cast<std::size_t>(cols_b.extent(0)); ++j) combined_indices.insert(cols_b(j));
            }

            std::vector<GO> combined_indices_vector(combined_indices.begin(), combined_indices.end());
            graph->insertGlobalIndices(global_row_index, Teuchos::ArrayView<const GO>(combined_indices_vector));
        }

        // Close via endAssembly which calls endFill() -> activeCrsGraph_ = OWNED -> fillComplete
        graph->endAssembly();
        return graph;
    }

    /**
     * @brief Copy values from one matrix to another
     * @details It is assumed that the sparsity of both matrices is compatible
     * @param rA The matrix where assigning values
     * @param rB The matrix to be copied
     */
    inline static void CopyMatrixValues(
        MatrixType& rA,
        const MatrixType& rB
        )
    {
        // Open for FE assembly if the matrix is closed, then zero all entries.
        // Do NOT call resumeFill() before beginAssembly() on FECrsMatrix: beginAssembly
        // expects a closed matrix and handles the transition itself.
        auto p_fe_rA_copy = dynamic_cast<MatrixType*>(&rA);
        if (p_fe_rA_copy) {
            // Handle the FECrsMatrix post-construction anomaly:
            // constructor leaves isFillActive()=true but fillState_=closed.
            // Normalize with fillComplete() so beginAssembly() can proceed.
            if (rA.isFillActive()) {
                rA.fillComplete();
            }
            rA.beginAssembly();
        } else if (!rA.isFillActive()) {
            rA.resumeFill();
        }
        rA.setAllToScalar(static_cast<ST>(0));

        for (LO i = 0; i < static_cast<LO>(rB.getNodeNumRows()); ++i) {
            const auto global_row_index = rB.getRowMap()->getGlobalElement(i);
            typename MatrixType::local_inds_host_view_type local_cols_b;
            typename MatrixType::values_host_view_type vals;
            rB.getLocalRowView(i, local_cols_b, vals);

            if (vals.extent(0) > 0) {
                Teuchos::Array<GO> global_cols(local_cols_b.extent(0));
                for (std::size_t j = 0; j < static_cast<std::size_t>(local_cols_b.extent(0)); ++j) {
                    global_cols[j] = rB.getColMap()->getGlobalElement(local_cols_b(j));
                }

                // Sum values into global matrix using global row and column indices
                rA.sumIntoGlobalValues(global_row_index, static_cast<LO>(global_cols.size()), vals.data(), global_cols.data());
            }
        }

        // Finalizing the fill process: endAssembly handles fillComplete internally
        rA.endAssembly();
    }

    /**
     * @brief This method checks and corrects the zero diagonal values
     * @details This method returns the scale norm considering scaling the diagonal
     * @param rProcessInfo The problem process info
     * @param rA The LHS matrix
     * @param rb The RHS vector
     * @param ScalingDiagonal The type of scaling diagonal considered
     * @return The scale norm
     */
    inline static double CheckAndCorrectZeroDiagonalValues(
        const ProcessInfo& rProcessInfo,
        MatrixType& rA,
        VectorType& rb,
        const SCALING_DIAGONAL ScalingDiagonal = SCALING_DIAGONAL::NO_SCALING
        )
    {
        KRATOS_TRY

        // Define zero value tolerance
        const double zero_tolerance = std::numeric_limits<double>::epsilon();

        // The diagonal considered
        const double scale_factor = GetScaleNorm(rProcessInfo, rA, ScalingDiagonal);

        auto localMatrix = rA.getLocalMatrixHost();
        auto rowMap = rA.getRowMap();
        auto colMap = rA.getColMap();
        auto localRhs = rb.getLocalViewHost(Tpetra::Access::ReadWrite);

        for (int i = 0; i < localMatrix.numRows(); ++i) {
            auto localRow = localMatrix.row(i);
            const auto row_gid = rowMap->getGlobalElement(i);
            bool empty = true;
            int j;
            for (j = 0; j < localRow.length; ++j) {
                const auto col_gid = colMap->getGlobalElement(localRow.colidx(j));
                // Check diagonal value
                if (col_gid == row_gid) {
                    if (std::abs(localRow.value(j)) > zero_tolerance) {
                        empty = false;
                    }
                    break;
                }
            }

            // If diagonal empty assign scale factor
            if (empty) {
                // replaceValues expects LOCAL column indices, not global IDs.
                // Convert the global diagonal column ID to a local column index.
                // This is safe as long as the colMap contains the diagonal entry,
                // which is guaranteed for owned rows in any properly-structured Tpetra matrix.
                if (!colMap.is_null()) {
                    const LO local_diag_col = static_cast<LO>(colMap->getLocalElement(row_gid));
                    if (local_diag_col != Teuchos::OrdinalTraits<LO>::invalid()) {
                        localMatrix.replaceValues(i, &local_diag_col, 1, &scale_factor, false, true);
                    }
                }
                localRhs(i, 0) = 0.0;
            }
        }

        return scale_factor;

        KRATOS_CATCH("")
    }

    /**
     * @brief This method returns the scale norm considering for scaling the diagonal
     * @param rProcessInfo The problem process info
     * @param rA The LHS matrix
     * @param ScalingDiagonal The type of scaling diagonal considered
     * @return The scale norm
     */
    inline static double GetScaleNorm(
        const ProcessInfo& rProcessInfo,
        const MatrixType& rA,
        const SCALING_DIAGONAL ScalingDiagonal = SCALING_DIAGONAL::NO_SCALING
        )
    {
        KRATOS_TRY

        switch (ScalingDiagonal) {
            case SCALING_DIAGONAL::NO_SCALING:
                return 1.0;
            case SCALING_DIAGONAL::CONSIDER_PRESCRIBED_DIAGONAL: {
                KRATOS_ERROR_IF_NOT(rProcessInfo.Has(BUILD_SCALE_FACTOR)) << "Scale factor not defined at process info" << std::endl;
                return rProcessInfo.GetValue(BUILD_SCALE_FACTOR);
            }
            case SCALING_DIAGONAL::CONSIDER_NORM_DIAGONAL:
                return GetDiagonalNorm(rA)/static_cast<double>(Size1(rA));
            case SCALING_DIAGONAL::CONSIDER_MAX_DIAGONAL:
                return GetMaxDiagonal(rA);
            default:
                return GetMaxDiagonal(rA);
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief This method returns the diagonal norm considering for scaling the diagonal
     * @param rA The LHS matrix
     * @return The diagonal norm
     */
    inline static double GetDiagonalNorm(const MatrixType& rA)
    {
        KRATOS_TRY

        // Create a plain Vector (not FEMultiVector) for diagonal copy — getLocalDiagCopy requires Vector
        Tpetra::Vector<ST, LO, GO, NT> diag(rA.getRowMap());

        // Extract the diagonal entries
        rA.getLocalDiagCopy(diag);

        // Get the local view of the diagonal
        auto diagLocalView = diag.getLocalViewHost(Tpetra::Access::ReadOnly);

        // Compute the local sum of squares of the diagonal
        ST localSumOfSquares = Teuchos::ScalarTraits<ST>::zero();  // Initialize to 0

        auto numLocalEntries = diag.getLocalLength();
        ST value = 0.0;
        for (std::size_t i = 0; i < numLocalEntries; ++i) {
            value = diagLocalView(i, 0);
            localSumOfSquares += value * value;  // Sum of squares
        }

        // Perform a global reduction to sum the squares across all processes
        ST globalSumOfSquares = 0.0;
        Teuchos::reduceAll(*rA.getMap()->getComm(), Teuchos::REDUCE_SUM, localSumOfSquares, Teuchos::outArg(globalSumOfSquares));

        // Compute the two-norm by taking the square root of the global sum of squares
        return std::sqrt(globalSumOfSquares);

        KRATOS_CATCH("");
    }

    /**
     * @brief This method returns the average of the max and min diagonal values
     * @param rA The LHS matrix
     * @return The average diagonal value (0.5 * (max + min))
     */
    inline static double GetAveragevalueDiagonal(const MatrixType& rA)
    {
        KRATOS_TRY

        return 0.5 * (GetMaxDiagonal(rA) + GetMinDiagonal(rA));

        KRATOS_CATCH("");
    }

    /**
     * @brief This method returns the maximum diagonal value
     * @param rA The LHS matrix
     * @return The maximum diagonal value
     */
    inline static double GetMaxDiagonal(const MatrixType& rA)
    {
        KRATOS_TRY

        // Create a plain Vector (not FEMultiVector) for diagonal copy — getLocalDiagCopy requires Vector
        Tpetra::Vector<ST, LO, GO, NT> diag(rA.getRowMap());

        // Extract the diagonal entries
        rA.getLocalDiagCopy(diag);

        // Get the local view of the diagonal
        auto diagLocalView = diag.getLocalViewHost(Tpetra::Access::ReadOnly);

        // Find the local maximum value
        ST localMax = Teuchos::ScalarTraits<ST>::zero();  // Initialize to 0

        auto numLocalEntries = diag.getLocalLength();
        ST value = 0.0;
        for (std::size_t i = 0; i < numLocalEntries; ++i) {
            value = diagLocalView(i, 0);
            if (value > localMax) {
                localMax = value;
            }
        }

        // Perform a global reduction to find the global maximum
        double globalMax = 0.0;
        Teuchos::reduceAll(*rA.getMap()->getComm(), Teuchos::REDUCE_MAX, localMax, Teuchos::outArg(globalMax));

        return globalMax;

        KRATOS_CATCH("");
    }

    /**
     * @brief This method returns the minimum diagonal value
     * @param rA The LHS matrix
     * @return The minimum diagonal value
     */
    inline static double GetMinDiagonal(const MatrixType& rA)
    {
        KRATOS_TRY

        // Create a plain Vector (not FEMultiVector) for diagonal copy — getLocalDiagCopy requires Vector
        Tpetra::Vector<ST, LO, GO, NT> diag(rA.getRowMap());

        // Extract the diagonal entries
        rA.getLocalDiagCopy(diag);

        // Get the local view of the diagonal
        auto diagLocalView = diag.getLocalViewHost(Tpetra::Access::ReadOnly);

        // Find the local minimum value
        ST localMin = Teuchos::ScalarTraits<ST>::rmax();  // Initialize to the max possible value

        auto numLocalEntries = diag.getLocalLength();
        ST value = 0.0;
        for (std::size_t i = 0; i < numLocalEntries; ++i) {
            value = diagLocalView(i, 0);
            if (value < localMin) {
                localMin = value;
            }
        }

        // Perform a global reduction to find the global minimum
        double globalMin = 0.0;
        Teuchos::reduceAll(*rA.getMap()->getComm(), Teuchos::REDUCE_MIN, localMin, Teuchos::outArg(globalMin));

        return globalMin;

        KRATOS_CATCH("");
    }

    /**
     * @brief Read a matrix from a MatrixMarket file
     * @param rFileName The name of the file to read
     * @param rComm The Tpetra MPI communicator
     * @return The FECrsMatrix read from the file
     */
    inline static MatrixPointerType ReadMatrixMarket(
        const std::string& rFileName,
        CommunicatorType& rComm
        )
    {
        return AuxiliaryMatrixMarket<MatrixType, VectorType>::ReadMatrixMarket(rFileName, rComm);
    }

    /**
     * @brief Read a vector from a MatrixMarket file
     * @param rFileName The name of the file to read
     * @param pComm The Tpetra MPI communicator
     * @param N The global size of the vector
     * @return The FEMultiVector read from the file
     */
    inline static VectorPointerType ReadMatrixMarketVector(
        const std::string& rFileName,
        CommunicatorPointerType pComm,
        const int N
        )
    {
        return AuxiliaryMatrixMarket<MatrixType, VectorType>::ReadMatrixMarketVector(rFileName, pComm, N);
    }

    /**
     * @brief Writes a matrix to a file in MatrixMarket format
     * @param pFileName The name of the file to be written
     * @param rA The matrix to be written
     * @param Symmetric If the matrix is symmetric (unused; kept for API compatibility)
     */
    static void WriteMatrixMarketMatrix(
        const char* pFileName,
        const MatrixType& rA,
        const bool Symmetric
        )
    {
        AuxiliaryMatrixMarket<MatrixType, VectorType>::WriteMatrixMarketMatrix(pFileName, rA, Symmetric);
    }

    /**
     * @brief Writes a vector to a file in MatrixMarket format
     * @param pFileName The name of the file to be written
     * @param rV The vector to be written
     */
    static void WriteMatrixMarketVector(
        const char* pFileName,
        const VectorType& rV
        )
    {
        AuxiliaryMatrixMarket<MatrixType, VectorType>::WriteMatrixMarketVector(pFileName, rV);
    }

    /**
     * @brief Creates a new dof updater
     * @return The new dof updater
     */
    inline static DofUpdaterPointerType CreateDofUpdater()
    {
        return DofUpdaterPointerType(new DofUpdater<TrilinosSpaceExperimental<TMatrixType, TVectorType>>());
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

    /**
     * @brief Turn back information as a string.
     * @return Info as a string.
     */
    virtual std::string Info() const
    {
        return "TrilinosSpaceExperimental";
    }

    /**
     * @brief Print information about this object.
     * @param rOStream The output stream to print on.
     */
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TrilinosSpaceExperimental";
    }

    /**
     * @brief Print object's data.
     * @param rOStream The output stream to print on.
     */
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}
private:
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TrilinosSpaceExperimental & operator=(TrilinosSpaceExperimental const& rOther);

    /// Copy constructor.
    TrilinosSpaceExperimental(TrilinosSpaceExperimental const& rOther);

    ///@}
}; // Class TrilinosSpaceExperimental

///@}

} // namespace Kratos.