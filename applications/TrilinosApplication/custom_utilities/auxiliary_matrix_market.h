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
#include <string>
#include <cstddef>

// External includes

/* Trilinos includes — safe headers (no include-guard conflicts) */
#include <Tpetra_FECrsMatrix.hpp>
#include <Tpetra_FECrsGraph.hpp>
#include <Tpetra_FEMultiVector.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultMpiComm.hpp>

/**
 * @file auxiliary_matrix_market.h
 * @brief Declares AuxiliaryMatrixMarket, an isolation wrapper for MatrixMarket I/O.
 *
 * WHY THIS CLASS EXISTS
 * ---------------------
 * MatrixMarket_Tpetra.hpp (pulled in by Tpetra's MatrixMarket reader/writer API)
 * internally includes "mmio_Tpetra.h", which reuses the include guard MM_IO_H.
 * That same guard is also used by Kratos' own "mmio.h" (the classical NIST mmio
 * library), but with a different signature for mm_read_mtx_crd_size (size_t* vs
 * int*).  When a unity-build batch merges an Epetra-based translation unit (which
 * has already pulled in Kratos' mmio.h) with a Tpetra-based translation unit, the
 * include guard causes mmio_Tpetra.h to be silently skipped, leaving the wrong
 * declaration in scope and producing a hard compilation error.
 *
 * By confining the #include <MatrixMarket_Tpetra.hpp> exclusively to
 * auxiliary_matrix_market.cpp — and marking that file SKIP_UNITY_BUILD_INCLUSION
 * in CMakeLists.txt — we ensure it always compiles as its own translation unit
 * with a clean include stack.  All other Tpetra/Epetra files remain free to
 * participate in unity builds without workarounds.
 *
 * The class provides the same four MatrixMarket operations that were previously
 * inlined in TrilinosSpaceExperimental (ReadMatrixMarket, ReadMatrixMarketVector,
 * WriteMatrixMarketMatrix, WriteMatrixMarketVector).  TrilinosSpaceExperimental
 * simply delegates to AuxiliaryMatrixMarket, keeping the public API unchanged.
 */

namespace Kratos
{

/**
 * @class AuxiliaryMatrixMarket
 * @ingroup TrilinosApplication
 * @brief Isolation wrapper for Tpetra MatrixMarket read/write operations.
 * @details See file-level documentation above for the rationale.
 *          The class is templated on the FE matrix and vector types so that it
 *          can be used from the equally-templated TrilinosSpaceExperimental.
 *          Only the explicit specialisation for the default Tpetra scalar/ordinal
 *          types (<Tpetra::FECrsMatrix<>, Tpetra::FEMultiVector<>>) is compiled
 *          and linked via the explicit-instantiation definition in
 *          auxiliary_matrix_market.cpp.
 * @tparam TMatrixType The FE sparse matrix type (e.g. Tpetra::FECrsMatrix<>)
 * @tparam TVectorType The FE multi-vector type (e.g. Tpetra::FEMultiVector<>)
 * @author Vicente Mataix Ferrandiz
 */
template<class TMatrixType, class TVectorType>
class AuxiliaryMatrixMarket
{
public:
    ///@name Type Definitions
    ///@{

    /// Scalar type (e.g. double)
    using ST = typename TMatrixType::scalar_type;
    /// Local ordinal type (e.g. int)
    using LO = typename TMatrixType::local_ordinal_type;
    /// Global ordinal type (e.g. int or long)
    using GO = typename TMatrixType::global_ordinal_type;
    /// Kokkos node type
    using NT = typename TMatrixType::node_type;

    /// FE sparse matrix
    using MatrixType = TMatrixType;
    /// FE multi-vector
    using VectorType = TVectorType;

    /// Standard (non-FE) CRS matrix — required by Tpetra MatrixMarket reader/writer
    using CrsMatrixType = Tpetra::CrsMatrix<ST, LO, GO, NT>;

    /// FE graph for constructing FE matrices
    using GraphType = Tpetra::FECrsGraph<LO, GO, NT>;

    /// Tpetra row-map type
    using MapType = Tpetra::Map<LO, GO, NT>;
    using MapPointerType = Teuchos::RCP<const MapType>;

    /// Managed pointer types
    using MatrixPointerType = Teuchos::RCP<MatrixType>;
    using VectorPointerType = Teuchos::RCP<VectorType>;

    /// MPI communicator
    using CommunicatorType = Teuchos::MpiComm<int>;
    using CommunicatorPointerType = Teuchos::RCP<const CommunicatorType>;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Read a sparse matrix from a MatrixMarket file.
     * @details Reads into a plain CrsMatrix first (required by the Tpetra reader),
     *          then builds an FECrsGraph mirroring the sparsity pattern and
     *          constructs an FECrsMatrix from it.
     * @param rFileName Path to the MatrixMarket (.mm / .mtx) file
     * @param rComm     Tpetra MPI communicator
     * @return Managed pointer to the filled FECrsMatrix
     */
    static MatrixPointerType ReadMatrixMarket(
        const std::string& rFileName,
        CommunicatorType& rComm
    );

    /**
     * @brief Read a dense vector from a MatrixMarket file.
     * @details Reads into a plain MultiVector first, then copies into an
     *          FEMultiVector (which has a deleted copy constructor).
     * @param rFileName Path to the MatrixMarket file
     * @param pComm     Shared pointer to the Tpetra communicator
     * @param N         Global length of the vector
     * @return Managed pointer to the FEMultiVector
     */
    static VectorPointerType ReadMatrixMarketVector(
        const std::string& rFileName,
        CommunicatorPointerType pComm,
        const int N
    );

    /**
     * @brief Write a sparse matrix to a MatrixMarket file.
     * @details FECrsMatrix inherits from CrsMatrix, so the Tpetra Writer can
     *          accept it after a static pointer cast.
     * @param pFileName  Output file path (C string)
     * @param rA         Matrix to write
     * @param Symmetric  Kept for API compatibility; not used by the Tpetra Writer
     */
    static void WriteMatrixMarketMatrix(
        const char* pFileName,
        const MatrixType& rA,
        const bool Symmetric
    );

    /**
     * @brief Write a dense vector to a MatrixMarket file.
     * @details FEMultiVector inherits from MultiVector, so it is cast before
     *          being handed to the Tpetra Writer.
     * @param pFileName Output file path (C string)
     * @param rV        Vector to write
     */
    static void WriteMatrixMarketVector(
        const char* pFileName,
        const VectorType& rV
    );

    ///@}

}; // class AuxiliaryMatrixMarket

// ---------------------------------------------------------------------------
// Explicit-instantiation declaration — the definition lives in
// auxiliary_matrix_market.cpp so that MatrixMarket_Tpetra.hpp is compiled in
// a single dedicated translation unit (SKIP_UNITY_BUILD_INCLUSION).
// ---------------------------------------------------------------------------
extern template class AuxiliaryMatrixMarket<Tpetra::FECrsMatrix<>, Tpetra::FEMultiVector<>>;

} // namespace Kratos
