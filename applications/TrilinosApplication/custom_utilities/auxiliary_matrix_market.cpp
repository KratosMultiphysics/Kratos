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

/**
 * @file auxiliary_matrix_market.cpp
 * @brief Implementation of AuxiliaryMatrixMarket — the *only* file that includes
 *        MatrixMarket_Tpetra.hpp in the TrilinosApplication.
 *
 * Keep this file marked SKIP_UNITY_BUILD_INCLUSION in CMakeLists.txt so that
 * it is always compiled as its own, isolated translation unit.  See
 * auxiliary_matrix_market.h for a full explanation of the include-guard conflict
 * that makes this isolation necessary.
 */

// Project includes — must come before MatrixMarket_Tpetra.hpp so that the safe
// Tpetra headers (pulled in by auxiliary_matrix_market.h) are processed first.
#include "custom_utilities/auxiliary_matrix_market.h"
#include "trilinos_application.h"

// External includes — the conflicting header, isolated here intentionally.
// MatrixMarket_Tpetra.hpp includes mmio_Tpetra.h, which reuses the MM_IO_H
// include guard.  Compiling this file alone guarantees that the guard has not
// already been claimed by Kratos' own mmio.h.
#include <MatrixMarket_Tpetra.hpp>

namespace Kratos
{

///@name AuxiliaryMatrixMarket Implementation
///@{

template<class TMatrixType, class TVectorType>
typename AuxiliaryMatrixMarket<TMatrixType, TVectorType>::MatrixPointerType
AuxiliaryMatrixMarket<TMatrixType, TVectorType>::ReadMatrixMarket(
    const std::string& rFileName,
    CommunicatorType& rComm
    )
{
    KRATOS_TRY

    // Read into a regular (non-FE) CrsMatrix first — the Tpetra reader requires it.
    auto p_comm = Teuchos::rcp(&rComm, false);
    auto p_crs = Tpetra::MatrixMarket::Reader<CrsMatrixType>::readSparseFile(rFileName, p_comm);

    const auto p_row_map = p_crs->getRowMap();
    const auto p_col_map = p_crs->getColMap();
    const LO num_local_rows = static_cast<LO>(p_crs->getNodeNumRows());

    // Compute the per-row maximum entry count (FECrsGraph takes a scalar bound).
    std::size_t max_entries_per_row = 0;
    for (LO i = 0; i < num_local_rows; ++i) {
        max_entries_per_row = std::max(max_entries_per_row,
            static_cast<std::size_t>(p_crs->getNumEntriesInLocalRow(i)));
    }

    // Build an FECrsGraph mirroring the sparsity pattern of the read matrix.
    auto p_graph = Teuchos::rcp(new GraphType(p_row_map, p_col_map, max_entries_per_row));
    p_graph->beginAssembly();
    for (LO i = 0; i < num_local_rows; ++i) {
        const GO global_row = p_row_map->getGlobalElement(i);
        typename CrsMatrixType::local_inds_host_view_type local_cols;
        typename CrsMatrixType::values_host_view_type vals;
        p_crs->getLocalRowView(i, local_cols, vals);
        if (local_cols.extent(0) > 0) {
            Teuchos::Array<GO> global_cols(local_cols.extent(0));
            for (std::size_t j = 0; j < local_cols.extent(0); ++j) {
                global_cols[j] = p_col_map->getGlobalElement(local_cols(j));
            }
            p_graph->insertGlobalIndices(global_row,
                Teuchos::ArrayView<const GO>(global_cols.data(), global_cols.size()));
        }
    }
    p_graph->endAssembly();

    // Construct FECrsMatrix from the closed graph.
    auto p_matrix = Teuchos::rcp(new MatrixType(
        Teuchos::rcp_const_cast<const GraphType>(p_graph)));

    // Normalise fill state, then populate values from the read CrsMatrix.
    if (p_matrix->isFillActive()) {
        p_matrix->fillComplete();
    }
    p_matrix->beginAssembly();
    p_matrix->setAllToScalar(static_cast<ST>(0));
    for (LO i = 0; i < num_local_rows; ++i) {
        const GO global_row = p_row_map->getGlobalElement(i);
        typename CrsMatrixType::local_inds_host_view_type local_cols;
        typename CrsMatrixType::values_host_view_type vals;
        p_crs->getLocalRowView(i, local_cols, vals);
        if (vals.extent(0) > 0) {
            Teuchos::Array<GO> global_cols(local_cols.extent(0));
            for (std::size_t j = 0; j < local_cols.extent(0); ++j) {
                global_cols[j] = p_col_map->getGlobalElement(local_cols(j));
            }
            p_matrix->sumIntoGlobalValues(global_row, static_cast<LO>(global_cols.size()),
                vals.data(), global_cols.data());
        }
    }
    p_matrix->endAssembly();

    return p_matrix;

    KRATOS_CATCH("")
}

template<class TMatrixType, class TVectorType>
typename AuxiliaryMatrixMarket<TMatrixType, TVectorType>::VectorPointerType
AuxiliaryMatrixMarket<TMatrixType, TVectorType>::ReadMatrixMarketVector(
    const std::string& rFileName,
    CommunicatorPointerType pComm,
    const int N
    )
{
    KRATOS_TRY

    // Create a uniform contiguous map over all N global DOFs.
    MapPointerType p_map = Teuchos::rcp(new MapType(static_cast<GO>(N), 0, pComm));

    // Read dense file into a plain (non-FE) MultiVector.
    auto p_mv = Tpetra::MatrixMarket::Reader<CrsMatrixType>::readDenseFile(rFileName, pComm, p_map);

    // Create the FEMultiVector (whose copy constructor is deleted) and copy via update().
    VectorPointerType p_vector;
    if constexpr (std::is_same_v<VectorType, Tpetra::FEMultiVector<ST, LO, GO, NT>>) {
        p_vector = Teuchos::rcp(new VectorType(p_map, Teuchos::null, 1));
    } else {
        p_vector = Teuchos::rcp(new VectorType(p_map));
    }
    p_vector->update(static_cast<ST>(1), *p_mv, static_cast<ST>(0));

    return p_vector;

    KRATOS_CATCH("")
}

template<class TMatrixType, class TVectorType>
void AuxiliaryMatrixMarket<TMatrixType, TVectorType>::WriteMatrixMarketMatrix(
    const char* pFileName,
    const MatrixType& rA,
    const bool Symmetric
    )
{
    KRATOS_TRY

    // FECrsMatrix inherits from CrsMatrix; cast so Writer can accept it.
    auto p_crs = Teuchos::rcp_static_cast<const CrsMatrixType>(Teuchos::rcpFromRef(rA));
    Tpetra::MatrixMarket::Writer<CrsMatrixType>::writeSparseFile(std::string(pFileName), p_crs);

    KRATOS_CATCH("")
}

template<class TMatrixType, class TVectorType>
void AuxiliaryMatrixMarket<TMatrixType, TVectorType>::WriteMatrixMarketVector(
    const char* pFileName,
    const VectorType& rV
    )
{
    KRATOS_TRY

    // FEMultiVector inherits from MultiVector; cast so Writer can accept it.
    using MultiVectorType = Tpetra::MultiVector<ST, LO, GO, NT>;
    auto p_mv = Teuchos::rcp_static_cast<const MultiVectorType>(Teuchos::rcpFromRef(rV));
    Tpetra::MatrixMarket::Writer<CrsMatrixType>::writeDenseFile(std::string(pFileName), p_mv);

    KRATOS_CATCH("")
}

///@}

// ---------------------------------------------------------------------------
// Explicit instantiation — provides the compiled symbols for the only
// concrete type pair used in the TrilinosApplication.  The declaration-side
// "extern template" in the header prevents the compiler from re-generating
// these symbols in every translation unit that includes the header.
// ---------------------------------------------------------------------------
template class AuxiliaryMatrixMarket<Tpetra::FECrsMatrix<>, Tpetra::FEMultiVector<>>;

} // namespace Kratos
