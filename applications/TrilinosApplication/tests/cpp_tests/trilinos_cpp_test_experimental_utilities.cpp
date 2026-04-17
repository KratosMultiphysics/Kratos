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

// System includes

// External includes

// Project includes
#include "includes/expect.h"
#include "trilinos_cpp_test_experimental_utilities.h"
#include "trilinos_cpp_test_utilities.h"

namespace Kratos
{

TrilinosCPPTestExperimentalUtilities::TrilinosLocalMatrixType TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalMatrix(
    const int NumGlobalElements,
    const double Offset,
    const bool AddNoDiagonalValues,
    const int InitialCapacity
    )
{
    return TrilinosCPPTestUtilities::GenerateDummyLocalMatrix(NumGlobalElements, Offset, AddNoDiagonalValues);
}

/***********************************************************************************/
/***********************************************************************************/

TrilinosCPPTestExperimentalUtilities::MatrixPointerType TrilinosCPPTestExperimentalUtilities::GenerateDummySparseMatrix(
    const DataCommunicator& rDataCommunicator,
    const int NumGlobalElements,
    const double Offset,
    const bool AddNoDiagonalValues,
    const int InitialCapacity
    )
{
    // Generate Tpetra communicator
    KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rDataCommunicator);
    Teuchos::RCP<const Teuchos::Comm<int>> tpetra_comm = Teuchos::rcp(new CommunicatorType(raw_mpi_comm));

    // Create a map
    Tpetra::global_size_t numGlobalElements = static_cast<Tpetra::global_size_t>(NumGlobalElements);
    Teuchos::RCP<const MapType> map = Teuchos::rcp(new MapType(numGlobalElements, 0, tpetra_comm));
    KRATOS_ERROR_IF_NOT(map->getGlobalNumElements() == static_cast<std::size_t>(NumGlobalElements)) << "Inconsistent number of rows" << std::endl;

    // Create a graph
    const std::size_t maxNodeEntriesPerRow = AddNoDiagonalValues ? 3 : 1;
    Teuchos::RCP<GraphType> graph = Teuchos::rcp(new GraphType(map, map, maxNodeEntriesPerRow));

    // Begin graph assembly
    graph->resumeFill();

    const int NumMyElements = map->getNodeNumElements();
    auto MyGlobalElements = map->getNodeElementList();
    std::vector<GO> nonDiagonalIndices(2); // Define based on your Tpetra types

    for (int i = 0; i < NumMyElements; ++i) {
        // Non-diagonal values
        if (AddNoDiagonalValues) {
            if (MyGlobalElements[i] == 0) {
                nonDiagonalIndices[0] = 1;
                graph->insertGlobalIndices(MyGlobalElements[i], Teuchos::ArrayView<const GO>(&nonDiagonalIndices[0], 1));
            } else if (MyGlobalElements[i] == NumGlobalElements - 1) {
                nonDiagonalIndices[0] = NumGlobalElements - 2;
                graph->insertGlobalIndices(MyGlobalElements[i], Teuchos::ArrayView<const GO>(&nonDiagonalIndices[0], 1));
            } else {
                nonDiagonalIndices[0] = MyGlobalElements[i] - 1;
                nonDiagonalIndices[1] = MyGlobalElements[i] + 1;
                graph->insertGlobalIndices(MyGlobalElements[i], Teuchos::ArrayView<const GO>(nonDiagonalIndices.data(), 2));
            }
        }

        // Insert diagonal entry
        graph->insertGlobalIndices(MyGlobalElements[i], Teuchos::ArrayView<const GO>(&MyGlobalElements[i], 1));
    }

    // End graph assembly
    if (graph->isFillActive()) graph->fillComplete();

    // Create a Tpetra::FECrsMatrix using the graph
    MatrixPointerType A = Teuchos::rcp(new TrilinosSparseMatrixType(graph));

    // Define non-diagonal values and other necessary variables
    std::vector<double> nonDiagonalValues(2, -1.0);
    double value;

    // Begin matrix assembly
    A->beginAssembly();

    for (int i = 0; i < NumMyElements; ++i) {
        // Non-diagonal values
        if (AddNoDiagonalValues) {
            if (MyGlobalElements[i] == 0) {
                nonDiagonalIndices[0] = 1;
                A->replaceGlobalValues(MyGlobalElements[i], Teuchos::ArrayView<const GO>(&nonDiagonalIndices[0], 1), Teuchos::ArrayView<const double>(&nonDiagonalValues[0], 1));
            } else if (MyGlobalElements[i] == NumGlobalElements - 1) {
                nonDiagonalIndices[0] = NumGlobalElements - 2;
                A->replaceGlobalValues(MyGlobalElements[i], Teuchos::ArrayView<const GO>(&nonDiagonalIndices[0], 1), Teuchos::ArrayView<const double>(&nonDiagonalValues[0], 1));
            } else {
                nonDiagonalIndices[0] = MyGlobalElements[i] - 1;
                nonDiagonalIndices[1] = MyGlobalElements[i] + 1;
                A->replaceGlobalValues(MyGlobalElements[i], Teuchos::ArrayView<const GO>(nonDiagonalIndices.data(), 2), Teuchos::ArrayView<const double>(nonDiagonalValues.data(), 2));
            }
        }

        // Insert diagonal entry
        value = Offset + static_cast<double>(MyGlobalElements[i]);
        A->replaceGlobalValues(MyGlobalElements[i], Teuchos::ArrayView<const GO>(&MyGlobalElements[i], 1), Teuchos::ArrayView<const double>(&value, 1));
    }

    // End matrix assembly
    A->endAssembly();

    // Complete the fill process, optimizing data for matrix-vector products
    if (A->isFillActive()) A->fillComplete();

    return A;
}

/***********************************************************************************/
/***********************************************************************************/

TrilinosCPPTestExperimentalUtilities::TrilinosLocalVectorType TrilinosCPPTestExperimentalUtilities::GenerateDummyLocalVector(
    const int NumGlobalElements,
    const double Offset
    )
{
    return TrilinosCPPTestUtilities::GenerateDummyLocalVector(NumGlobalElements, Offset);
}

/***********************************************************************************/
/***********************************************************************************/

TrilinosCPPTestExperimentalUtilities::VectorPointerType TrilinosCPPTestExperimentalUtilities::GenerateDummySparseVector(
    const DataCommunicator& rDataCommunicator,
    const int NumGlobalElements,
    const double Offset
    )
{
    // Generate Tpetra communicator
    KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rDataCommunicator);
    Teuchos::RCP<const Teuchos::Comm<int>> tpetra_comm = Teuchos::rcp(new CommunicatorType(Teuchos::opaqueWrapper(raw_mpi_comm)));

    // Create a map
    Teuchos::RCP<const MapType> map = Teuchos::rcp(new MapType(NumGlobalElements, 0, tpetra_comm));

    // Local number of rows
    const std::size_t NumMyElements = map->getNodeNumElements();

    // Get update list
    auto MyGlobalElements = map->getNodeElementList();

    // Create a Tpetra_Vector
    Teuchos::RCP<Tpetra::FEMultiVector<>> b = Teuchos::rcp(new Tpetra::FEMultiVector<>(map, Teuchos::null, 1));

    // Fill the vector with values
    double value;
    for (std::size_t i = 0; i < NumMyElements; ++i) {
        value = Offset + static_cast<double>(MyGlobalElements[i]);
        b->replaceLocalValue(i, size_t(0), value);
    }

    return b;
}

/***********************************************************************************/
/***********************************************************************************/

void TrilinosCPPTestExperimentalUtilities::CheckSparseVectorFromLocalVector(
    const TrilinosVectorType& rA,
    const TrilinosLocalVectorType& rB,
    const double NegligibleValueThreshold
    )
{
    const std::size_t total_size = rB.size();
    std::vector<int> indexes;
    indexes.reserve(total_size);
    std::vector<double> values;
    values.reserve(total_size);
    for (std::size_t i = 0; i < rB.size(); ++i) {
        indexes.push_back(i);
        values.push_back(rB[i]);
    }
    CheckSparseVector(rA, indexes, values);
}

/***********************************************************************************/
/***********************************************************************************/

void TrilinosCPPTestExperimentalUtilities::CheckSparseVector(
    const TrilinosVectorType& rb,
    const std::vector<int>& rIndexes,
    const std::vector<double>& rValues,
    const double NegligibleValueThreshold
    )
{
    int index;
    double value;
    const auto& r_map = rb.getMap(); // Direct access to the map associated with rb
    auto r_local_map = r_map->getLocalMap(); // Get the local map for local index operations
    for (std::size_t counter = 0; counter < rIndexes.size(); ++counter) {
        index = rIndexes[counter];
        value = rValues[counter];

        if (r_map->isNodeGlobalElement(index)) {
            auto lid = r_map->getLocalElement(index);
            const double ref_value = rb.getData(0)[lid];
            KRATOS_EXPECT_RELATIVE_NEAR(value, ref_value, NegligibleValueThreshold);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(
    const TrilinosSparseMatrixType& rA,
    const TrilinosLocalMatrixType& rB,
    const double NegligibleValueThreshold
    )
{
    const std::size_t total_size = rB.size1() * rB.size2();
    std::vector<int> row_indexes;
    row_indexes.reserve(total_size);
    std::vector<int> column_indexes;
    column_indexes.reserve(total_size);
    std::vector<double> values;
    values.reserve(total_size);
    for (std::size_t i = 0; i < rB.size1(); ++i) {
        for (std::size_t j = 0; j < rB.size2(); ++j) {
            const double value = rB(i, j);
            if (std::abs(value) > std::numeric_limits<double>::epsilon()) {
                row_indexes.push_back(i);
                column_indexes.push_back(j);
                values.push_back(value);
            }
        }
    }
    CheckSparseMatrix(rA, row_indexes, column_indexes, values);
}

/***********************************************************************************/
/***********************************************************************************/

void TrilinosCPPTestExperimentalUtilities::CheckSparseMatrixFromLocalMatrix(
    const TrilinosSparseMatrixType& rA,
    const std::vector<int>& rRowIndexes,
    const std::vector<int>& rColumnIndexes,
    const TrilinosLocalMatrixType& rB,
    const double NegligibleValueThreshold
    )
{

    int row, column;
    std::vector<double> values(rRowIndexes.size());
    for (std::size_t counter = 0; counter < rRowIndexes.size(); ++counter) {
        row = rRowIndexes[counter];
        column = rColumnIndexes[counter];
        values[counter] = rB(row, column);
    }
    CheckSparseMatrix(rA, rRowIndexes, rColumnIndexes, values, NegligibleValueThreshold);
}

/***********************************************************************************/
/***********************************************************************************/

void TrilinosCPPTestExperimentalUtilities::CheckSparseMatrix(
    const TrilinosSparseMatrixType& rA,
    const std::vector<int>& rRowIndexes,
    const std::vector<int>& rColumnIndexes,
    const std::vector<double>& rValues,
    const double NegligibleValueThreshold
    )
{
    int local_validated_values = 0;
    int row, column;
    double value;

    // Get necessary maps and communicator
    auto rowMap = rA.getRowMap();  // rowMap is an RCP<const MapType>
    auto colMap = rA.getColMap();  // colMap is an RCP<const MapType>
    auto comm = rA.getComm();      // comm is an RCP<const Teuchos::Comm<int>>

    for (std::size_t counter = 0; counter < rRowIndexes.size(); ++counter) {
        row = rRowIndexes[counter];
        column = rColumnIndexes[counter];
        value = rValues[counter];

        // Convert global row index to local index
        const auto localRowIndex = rowMap->getLocalElement(row);

        if (localRowIndex != Teuchos::OrdinalTraits<MapType::local_ordinal_type>::invalid()) {
            typename TrilinosSparseMatrixType::values_host_view_type vals;   // Row non-zero values
            typename TrilinosSparseMatrixType::local_inds_host_view_type cols; // Column indices of row non-zero values
            rA.getLocalRowView(localRowIndex, cols, vals);
            const std::size_t numEntries = cols.size(); // Number of non-zero entries
            for (std::size_t j = 0; j < numEntries; j++) {
                // Convert local column index to global index
                const int col_gid = colMap->getGlobalElement(cols[j]);
                if (col_gid == column) {
                    KRATOS_EXPECT_RELATIVE_NEAR(value, vals[j], NegligibleValueThreshold);
                    ++local_validated_values;
                    break;
                }
            }
        }
    }

    // Checking that all the values have been validated
    int global_validated_values;
    comm->barrier();  // Ensure synchronization across processes
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &local_validated_values, &global_validated_values);
    KRATOS_EXPECT_EQ(global_validated_values, static_cast<int>(rValues.size()));
}

/***********************************************************************************/
/***********************************************************************************/

void TrilinosCPPTestExperimentalUtilities::GenerateSparseMatrixIndexAndValuesVectors(
    const TrilinosSparseSpaceType::MatrixType& rA,
    std::vector<int>& rRowIndexes,
    std::vector<int>& rColumnIndexes,
    std::vector<double>& rValues,
    const bool PrintValues,
    const double ThresholdIncludeHardZeros
    )
{
    KRATOS_ERROR_IF_NOT(rA.getComm()->getSize() == 1) << "Debug must be done with one MPI core" << std::endl;

    // If print values
    if (PrintValues) {
        std::cout << "\n        KRATOS_EXPECT_EQ(rA.getGlobalNumRows(), " << rA.getGlobalNumRows() << ");\n";
        std::cout << "        KRATOS_EXPECT_EQ(rA.getGlobalNumCols(), " << rA.getGlobalNumCols() << ");\n";
        std::cout << "        KRATOS_EXPECT_EQ(rA.getGlobalNumEntries(), " << rA.getGlobalNumEntries() << ");\n";
    }

    std::vector<double> values;
    for (std::size_t i = 0; i < rA.getNodeNumRows(); i++) {
        typename TrilinosSparseSpaceType::MatrixType::values_host_view_type vals;   // Row non-zero values
        typename TrilinosSparseSpaceType::MatrixType::local_inds_host_view_type cols;      // Column indices of row non-zero values
        rA.getLocalRowView(i, cols, vals);
        const std::size_t numEntries = cols.size(); // Number of non-zero entries
        const int row_gid = rA.getRowMap()->getGlobalElement(i);
        for (std::size_t j = 0; j < numEntries; j++) {
            const int col_gid = rA.getColMap()->getGlobalElement(cols[j]);
            if (std::abs(vals[j]) > ThresholdIncludeHardZeros) {
                rRowIndexes.push_back(row_gid);
                rColumnIndexes.push_back(col_gid);
                rValues.push_back(vals[j]);
            }
        }
    }

    // If print values
    if (PrintValues) {
        std::cout << "\n        // Values to check\n";
        std::cout << "        std::vector<int> row_indexes = {";
        for (std::size_t i = 0; i < rRowIndexes.size() - 1; ++i) {
            std::cout << rRowIndexes[i] << ", ";
        }
        std::cout << rRowIndexes.back() << "};";
        std::cout << "\n        std::vector<int> column_indexes = {";
        for (std::size_t i = 0; i < rColumnIndexes.size() - 1; ++i) {
            std::cout << rColumnIndexes[i] << ", ";
        }
        std::cout << rColumnIndexes.back() << "};";
        std::cout << "\n        std::vector<double> values = {";
        for (std::size_t i = 0; i < rValues.size() - 1; ++i) {
            std::cout << std::fixed;
            std::cout << std::setprecision(16);
            std::cout << rValues[i] << ", ";
        }
        std::cout << std::fixed;
        std::cout << std::setprecision(16);
        std::cout << rValues.back() << "};" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

TrilinosCPPTestExperimentalUtilities::MatrixPointerType TrilinosCPPTestExperimentalUtilities::GenerateSparseMatrix(
    const DataCommunicator& rDataCommunicator,
    const int NumGlobalElements,
    const std::vector<int>& rRowIndexes,
    const std::vector<int>& rColumnIndexes,
    const std::vector<double>& rValues,
    const MapType* pMap
    )
{
    // Generate Tpetra communicator
    KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
    auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rDataCommunicator);
    auto comm = Teuchos::rcp(new CommunicatorType(raw_mpi_comm));

    // Create a map
    Teuchos::RCP<const MapType> map;
    if (pMap == nullptr) {
        map = Teuchos::rcp(new MapType(NumGlobalElements, 0, comm));
    } else {
        map = Teuchos::rcp(new MapType(*pMap));
    }

    // Local number of rows
    const int NumMyElements = map->getNodeNumElements();

    // Get update list
    auto MyGlobalElements = map->getNodeElementList();

    // Create an integer vector NumNz that is used to build the Tpetra Matrix.
    const int size_global_vector = rRowIndexes.size();
    std::vector<size_t> NumNz(NumMyElements, 0);
    int current_row_index;
    for (current_row_index = 0; current_row_index < size_global_vector; ++current_row_index) {
        if (MyGlobalElements[0] == rRowIndexes[current_row_index]) {
            break;
        }
    }
    int current_id = rRowIndexes[current_row_index];
    std::size_t nnz = 0;
    GO initial_index, end_index;
    std::unordered_map<int, std::pair<int, int>> initial_and_end_index;
    for (int i = 0; i < NumMyElements; i++) {
        if (MyGlobalElements[i] == rRowIndexes[current_row_index]) {
            initial_index = current_row_index;
            const int start_index = current_row_index;
            for (current_row_index = start_index; current_row_index < size_global_vector; ++current_row_index) {
                if (current_id == rRowIndexes[current_row_index]) {
                    ++nnz;
                } else {
                    current_id = rRowIndexes[current_row_index];
                    break;
                }
            }
            end_index = current_row_index;
            initial_and_end_index.insert({i, std::make_pair(initial_index, end_index)});
            NumNz[i] = nnz;
            nnz = 0;
        }
    }

    // Create a Tpetra_CrsMatrix graph
    const int NumNzMax = *std::max_element(NumNz.begin(), NumNz.end());
    Teuchos::RCP<GraphType> graph = Teuchos::rcp(new GraphType(map, map, NumNzMax));

    // Build the graph
    for (int i = 0; i < NumMyElements; ++i) {
        auto it_find = initial_and_end_index.find(i);
        if (it_find != initial_and_end_index.end()) {
            const auto& r_pair = it_find->second;
            initial_index = r_pair.first;
            end_index = r_pair.second;
            std::vector<GO> indexes(rColumnIndexes.begin() + initial_index, rColumnIndexes.begin() + end_index);
            // Insert global indices into the graph
            graph->insertGlobalIndices(MyGlobalElements[i], Teuchos::ArrayView<const GO>(indexes));
        }
    }

    // Finish up the graph
    if (graph->isFillActive()) graph->fillComplete();

    // Create the matrix and set values
    MatrixPointerType A = Teuchos::rcp(new TrilinosSparseMatrixType(graph));
    // Begin matrix assembly
    A->beginAssembly();
    // Set the matrix values
    for (int i = 0; i < NumMyElements; ++i) {
        auto it_find = initial_and_end_index.find(i);
        if (it_find != initial_and_end_index.end()) {
            const auto& r_pair = it_find->second;
            initial_index = r_pair.first;
            end_index = r_pair.second;
            std::vector<GO> indexes(rColumnIndexes.begin() + initial_index, rColumnIndexes.begin() + end_index);
            std::vector<double> values(rValues.begin() + initial_index, rValues.begin() + end_index);
            A->replaceGlobalValues(MyGlobalElements[i], Teuchos::ArrayView<const GO>(indexes), Teuchos::ArrayView<const double>(values));
        }
    }

    // End matrix assembly
    A->endAssembly();

    // Finish up, transforming the matrix entries into local numbering,
    // to optimize data transfer during matrix-vector products
    if (A->isFillActive()) A->fillComplete();

    return A;
}

} /// namespace Kratos