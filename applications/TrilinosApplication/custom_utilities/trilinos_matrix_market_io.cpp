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
 * @file trilinos_matrix_market_io.cpp
 * @brief Implementation of Tpetra MatrixMarket I/O helpers.
 *
 * This is the ONLY file that includes MatrixMarket_Tpetra.hpp, which pulls
 * in mmio_Tpetra.h.  Keep it marked SKIP_UNITY_BUILD_INCLUSION in
 * CMakeLists.txt so that it never shares a TU with Kratos' mmio.h.
 */

// Project includes — must come first so the safe Tpetra headers are seen.
#include "custom_utilities/trilinos_matrix_market_io.h"
#include "includes/define.h"

// External includes — the conflicting header, isolated here.
#include <MatrixMarket_Tpetra.hpp>

namespace Kratos::TpetraMatrixMarketIO
{

// Shorthand type aliases local to this TU.
using ST = FEMatrix::scalar_type;
using LO = FEMatrix::local_ordinal_type;
using GO = FEMatrix::global_ordinal_type;
using NT = FEMatrix::node_type;
using CrsMatrix  = Tpetra::CrsMatrix<ST, LO, GO, NT>;
using FEGraph    = Tpetra::FECrsGraph<LO, GO, NT>;
using MapType    = Tpetra::Map<LO, GO, NT>;

Teuchos::RCP<FEMatrix> ReadMatrix(
    const std::string& rFileName,
    Comm& rComm)
{
    KRATOS_TRY

    auto p_comm = Teuchos::rcp(&rComm, false);
    auto p_crs  = Tpetra::MatrixMarket::Reader<CrsMatrix>::readSparseFile(rFileName, p_comm);

    const auto p_row_map = p_crs->getRowMap();
    const auto p_col_map = p_crs->getColMap();
    const LO num_local_rows = static_cast<LO>(p_crs->getLocalNumRows());

    // Per-row max entry count for the FECrsGraph.
    std::size_t max_entries = 0;
    for (LO i = 0; i < num_local_rows; ++i) {
        max_entries = std::max(max_entries,
            static_cast<std::size_t>(p_crs->getNumEntriesInLocalRow(i)));
    }

    // Build FECrsGraph mirroring the sparsity pattern.
    auto p_graph = Teuchos::rcp(new FEGraph(p_row_map, p_col_map, max_entries));
    p_graph->beginAssembly();
    for (LO i = 0; i < num_local_rows; ++i) {
        const GO global_row = p_row_map->getGlobalElement(i);
        typename CrsMatrix::local_inds_host_view_type local_cols;
        typename CrsMatrix::values_host_view_type vals;
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

    // Construct FECrsMatrix from the graph and fill values.
    auto p_matrix = Teuchos::rcp(new FEMatrix(
        Teuchos::rcp_const_cast<const FEGraph>(p_graph)));

    if (p_matrix->isFillActive()) {
        p_matrix->fillComplete();
    }
    p_matrix->beginAssembly();
    p_matrix->setAllToScalar(static_cast<ST>(0));
    for (LO i = 0; i < num_local_rows; ++i) {
        const GO global_row = p_row_map->getGlobalElement(i);
        typename CrsMatrix::local_inds_host_view_type local_cols;
        typename CrsMatrix::values_host_view_type vals;
        p_crs->getLocalRowView(i, local_cols, vals);
        if (vals.extent(0) > 0) {
            Teuchos::Array<GO> global_cols(local_cols.extent(0));
            for (std::size_t j = 0; j < local_cols.extent(0); ++j) {
                global_cols[j] = p_col_map->getGlobalElement(local_cols(j));
            }
            p_matrix->sumIntoGlobalValues(global_row,
                static_cast<LO>(global_cols.size()),
                vals.data(), global_cols.data());
        }
    }
    p_matrix->endAssembly();

    return p_matrix;

    KRATOS_CATCH("")
}

Teuchos::RCP<FEVector> ReadVector(
    const std::string& rFileName,
    Teuchos::RCP<const Comm> pComm,
    const int N)
{
    KRATOS_TRY

    Teuchos::RCP<const MapType> p_map = Teuchos::rcp(new MapType(static_cast<GO>(N), 0, pComm));
    auto p_mv  = Tpetra::MatrixMarket::Reader<CrsMatrix>::readDenseFile(rFileName, pComm, p_map);

    auto p_vector = Teuchos::rcp(new FEVector(p_map, Teuchos::null, 1));
    p_vector->update(static_cast<ST>(1), *p_mv, static_cast<ST>(0));

    return p_vector;

    KRATOS_CATCH("")
}

void WriteMatrix(
    const char* pFileName,
    const FEMatrix& rA)
{
    KRATOS_TRY

    auto p_crs = Teuchos::rcp_static_cast<const CrsMatrix>(Teuchos::rcpFromRef(rA));
    Tpetra::MatrixMarket::Writer<CrsMatrix>::writeSparseFile(std::string(pFileName), p_crs);

    KRATOS_CATCH("")
}

void WriteVector(
    const char* pFileName,
    const FEVector& rV)
{
    KRATOS_TRY

    using MultiVector = Tpetra::MultiVector<ST, LO, GO, NT>;
    auto p_mv = Teuchos::rcp_static_cast<const MultiVector>(Teuchos::rcpFromRef(rV));
    Tpetra::MatrixMarket::Writer<CrsMatrix>::writeDenseFile(std::string(pFileName), p_mv);

    KRATOS_CATCH("")
}

} // namespace Kratos::TpetraMatrixMarketIO
