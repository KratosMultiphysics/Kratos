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

// External includes
#include <Tpetra_FECrsMatrix.hpp>
#include <Tpetra_FEMultiVector.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultMpiComm.hpp>

// Project includes
#include "includes/define.h"

/**
 * @file trilinos_matrix_market_io.h
 * @brief Free-function wrappers for Tpetra MatrixMarket I/O.
 *
 * MatrixMarket_Tpetra.hpp internally includes mmio_Tpetra.h which declares
 * mm_read_mtx_crd_size(FILE*, size_t*, size_t*, size_t*) with C linkage.
 * Kratos' own mmio.h declares the same function with (int*, int*, int*).
 * Because C linkage forbids overloads, the two headers cannot coexist in
 * the same translation unit.
 *
 * The implementations in trilinos_matrix_market_io.cpp are compiled in their
 * own TU (SKIP_UNITY_BUILD_INCLUSION) so that only mmio_Tpetra.h is visible.
 */

namespace Kratos::TpetraMatrixMarketIO
{

/// Concrete Tpetra types used throughout the TrilinosApplication.
using FEMatrix = Tpetra::FECrsMatrix<>;
using FEVector = Tpetra::FEMultiVector<>;
using Comm     = Teuchos::MpiComm<int>;

/// @brief Read a sparse FECrsMatrix from a MatrixMarket file.
Teuchos::RCP<FEMatrix> KRATOS_API(TRILINOS_APPLICATION) ReadMatrix(
    const std::string& rFileName,
    Comm& rComm);

/// @brief Read a dense FEMultiVector from a MatrixMarket file.
Teuchos::RCP<FEVector> KRATOS_API(TRILINOS_APPLICATION) ReadVector(
    const std::string& rFileName,
    Teuchos::RCP<const Comm> pComm,
    int N);

/// @brief Write a sparse matrix to a MatrixMarket file.
void KRATOS_API(TRILINOS_APPLICATION) WriteMatrix(
    const char* pFileName,
    const FEMatrix& rA);

/// @brief Write a dense vector to a MatrixMarket file.
void KRATOS_API(TRILINOS_APPLICATION) WriteVector(
    const char* pFileName,
    const FEVector& rV);

} // namespace Kratos::TpetraMatrixMarketIO
