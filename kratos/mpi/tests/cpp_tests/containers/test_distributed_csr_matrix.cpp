//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi

// System includes

// External includes

// Project includes
#include "containers/sparse_contiguous_row_graph.h"
#include "containers/distributed_sparse_graph.h"
#include "containers/distributed_csr_matrix.h"
#include "containers/distributed_system_vector.h"
#include "containers/distributed_vector_importer.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "mpi/testing/mpi_testing.h"
#include "mpi/tests/test_utilities/distributed_sparse_containers_test_utilities.h"
#include "mpi/utilities/amgcl_distributed_csr_conversion_utilities.h"
#include "mpi/utilities/amgcl_distributed_csr_spmm_utilities.h"
#include "utilities/openmp_utils.h"
#include "utilities/builtin_timer.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(DistributedCSRMatrixConstruction, KratosMPICoreFastSuite)
{
    using IndexType = DistributedSparseContainersTestUtilities::IndexType;
    const auto& r_comm = ParallelEnvironment::GetDefaultDataCommunicator();
    int world_size = r_comm.Size();
    int my_rank = r_comm.Rank();

    auto dofs_bounds = DistributedSparseContainersTestUtilities::ComputeBounds<IndexType>(40, world_size, my_rank);
    auto reference_A_map = DistributedSparseContainersTestUtilities::GetReferenceMatrixAsMap(dofs_bounds);

    auto el_bounds = DistributedSparseContainersTestUtilities::ComputeBounds<IndexType>(31, world_size, my_rank);
    const auto connectivities = DistributedSparseContainersTestUtilities::ElementConnectivities(el_bounds);

    DistributedSparseGraph<IndexType> A_graph(dofs_bounds[1]-dofs_bounds[0], r_comm);
    IndexPartition<IndexType>(connectivities.size()).for_each([&](IndexType i) {
        A_graph.AddEntries(connectivities[i]);
    });
    A_graph.Finalize();

    //FEM assembly
    DistributedCsrMatrix<double,IndexType> A(A_graph);
    A.BeginAssemble();
    for(const auto& c : connectivities) {
        Matrix data(c.size(),c.size(),1.0);
        A.Assemble(data,c);
    }
    A.FinalizeAssemble();

    for(const auto& item : reference_A_map) {
        auto global_i = item.first.first;
        auto global_j = item.first.second;
        const double reference_value = item.second;
        const double value = A.GetLocalDataByGlobalId(global_i, global_j);
        KRATOS_EXPECT_NEAR(reference_value, value, 1e-14);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DistributedCSRMatrixSpMV, KratosMPICoreFastSuite)
{
    using IndexType = DistributedSparseContainersTestUtilities::IndexType;
    const auto& r_comm = ParallelEnvironment::GetDefaultDataCommunicator();
    int world_size = r_comm.Size();
    int my_rank = r_comm.Rank();

    auto dofs_bounds = DistributedSparseContainersTestUtilities::ComputeBounds<IndexType>(40, world_size, my_rank);
    auto reference_A_map = DistributedSparseContainersTestUtilities::GetReferenceMatrixAsMap(dofs_bounds);

    auto el_bounds = DistributedSparseContainersTestUtilities::ComputeBounds<IndexType>(31, world_size, my_rank);
    const auto connectivities = DistributedSparseContainersTestUtilities::ElementConnectivities(el_bounds);

    DistributedSparseGraph<IndexType> A_graph(dofs_bounds[1]-dofs_bounds[0], r_comm);
    IndexPartition<IndexType>(connectivities.size()).for_each([&](IndexType i) {
        A_graph.AddEntries(connectivities[i]);
    });
    A_graph.Finalize();

    DistributedCsrMatrix<double, IndexType> A(A_graph);
    A.BeginAssemble();
    for(const auto& c : connectivities) {
        Matrix data(c.size(),c.size(),1.0);
        A.Assemble(data,c);
    }
    A.FinalizeAssemble();

    //Test SPMV
    //here we test SPMV by a vector of 1s
    DistributedSystemVector<> y(A_graph);
    DistributedSystemVector<> b(A_graph);
    y.SetValue(0.0);
    b.SetValue(1.0);

    A.SpMV(b,y); //y+=A*b

    std::vector<double> reference_spmv_res{4,12,8,12,12,12,20,24,16,16,8,16,12,4,24,12,20,12,24,12,0,4,16,4,16,16,24,4,20,8,8,12,4,16,4,20,8,16,4,12};
    for(unsigned int i=0; i<y.LocalSize(); ++i) {
        IndexType global_i = y.GetNumbering().GlobalId(i);
        KRATOS_EXPECT_NEAR(y[i],  reference_spmv_res[global_i], 1e-14 );
    }
}

KRATOS_TEST_CASE_IN_SUITE(DistributedCSRMatrixToAMGCLMatrix, KratosCoreFastSuite)
{
    using IndexType = DistributedSparseContainersTestUtilities::IndexType;
    const auto& r_comm = ParallelEnvironment::GetDefaultDataCommunicator();
    int world_size = r_comm.Size();
    int my_rank = r_comm.Rank();

    auto dofs_bounds = DistributedSparseContainersTestUtilities::ComputeBounds<IndexType>(40, world_size, my_rank);
    auto reference_A_map = DistributedSparseContainersTestUtilities::GetReferenceMatrixAsMap(dofs_bounds);

    auto el_bounds = DistributedSparseContainersTestUtilities::ComputeBounds<IndexType>(31, world_size, my_rank);
    const auto connectivities = DistributedSparseContainersTestUtilities::ElementConnectivities(el_bounds);

    DistributedSparseGraph<IndexType> A_graph(dofs_bounds[1]-dofs_bounds[0], r_comm);
    IndexPartition<IndexType>(connectivities.size()).for_each([&](IndexType i) {
        A_graph.AddEntries(connectivities[i]);
    });
    A_graph.Finalize();

    DistributedCsrMatrix<double, IndexType> A(A_graph);
    A.BeginAssemble();
    for(const auto& c : connectivities) {
        Matrix data(c.size(),c.size(),1.0);
        A.Assemble(data,c);
    }
    A.FinalizeAssemble();

    DistributedSystemVector<> y(A_graph);
    DistributedSystemVector<> b(A_graph);

    //testing AMGCL interface
    bool move_to_backend = true;
    auto offdiag_global_index2 = A.GetOffDiagonalIndex2DataInGlobalNumbering();
    auto pAmgcl = AmgclDistributedCSRConversionUtilities::ConvertToAmgcl<double,IndexType>(A, offdiag_global_index2, move_to_backend);

    y.SetValue(0.0);
    b.SetValue(1.0);
    pAmgcl->mul(1.0,b,1.0,y);

    std::vector<double> reference_spmv_res{4,12,8,12,12,12,20,24,16,16,8,16,12,4,24,12,20,12,24,12,0,4,16,4,16,16,24,4,20,8,8,12,4,16,4,20,8,16,4,12};
    for(unsigned int i=0; i<y.LocalSize(); ++i) {
        IndexType global_i = y.GetNumbering().GlobalId(i);
        KRATOS_EXPECT_NEAR(y[i], reference_spmv_res[global_i], 1e-14 );
    }

    //check round trip krtos_csr->amgcl->kratos_csr
    move_to_backend = false; //note that the "move_to_backend needs to be set to false to be able to do the round trip"
    offdiag_global_index2 = A.GetOffDiagonalIndex2DataInGlobalNumbering();
    pAmgcl = AmgclDistributedCSRConversionUtilities::ConvertToAmgcl<double,IndexType>(A, offdiag_global_index2, move_to_backend);

    //convert back to CSR matrix
    DistributedCsrMatrix<double, IndexType>::Pointer p_A_reconverted =
        AmgclDistributedCSRConversionUtilities::ConvertToCsrMatrix<double,IndexType>(*pAmgcl, r_comm);

    y.SetValue(0.0);
    b.SetValue(1.0);
    p_A_reconverted->SpMV(b,y); //y+=A*b

    for(unsigned int i = 0; i < y.LocalSize(); ++i) {
        IndexType global_i = y.GetNumbering().GlobalId(i);
        KRATOS_EXPECT_NEAR(y[i],  reference_spmv_res[global_i], 1e-14 );
    }
}

KRATOS_TEST_CASE_IN_SUITE(DistributedCSRMatrixRectangularMatrix, KratosMPICoreFastSuite)
{
    using IndexType = DistributedSparseContainersTestUtilities::IndexType;
    const auto& r_comm = ParallelEnvironment::GetDefaultDataCommunicator();
    IndexType col_divider = 3; //ratio of size between columns and row indices

    //*************************************************************************
    //compute reference solution - serial mode
    std::vector<IndexType> all_el_bounds{0,31};
    const auto all_connectivities = DistributedSparseContainersTestUtilities::ElementConnectivities(all_el_bounds);
    SparseContiguousRowGraph<IndexType> A_graph_serial(40);

    IndexPartition<IndexType>(all_connectivities.size()).for_each([&](IndexType i) {
        std::vector<IndexType> row_ids = all_connectivities[i];
        std::vector<IndexType> col_ids{row_ids[0]/col_divider, row_ids[1]/col_divider};
        A_graph_serial.AddEntries(row_ids, col_ids);
    });
    A_graph_serial.Finalize();

    CsrMatrix<double, IndexType> A_serial(A_graph_serial);

    A_serial.BeginAssemble();
    for(const auto& c : all_connectivities) {
        std::vector<IndexType> row_ids = c;
        std::vector<IndexType> col_ids{row_ids[0]/col_divider, row_ids[1]/col_divider};
        Matrix data(row_ids.size(),col_ids.size(),1.0);
        A_serial.Assemble(data,row_ids, col_ids);
    }
    A_serial.FinalizeAssemble();

    auto reference_A_map = A_serial.ToMap();

    //here we test SPMV by a vector of 1s
    SystemVector<> yserial(A_serial.size1()); //destination vector
    yserial.SetValue(0.0);

    SystemVector<> xserial(A_serial.size2()); //origin vector
    xserial.SetValue(1.0);

    A_serial.SpMV(xserial,yserial);

    //*************************************************************************
    int world_size =r_comm.Size();
    int my_rank = r_comm.Rank();

    auto dofs_bounds = DistributedSparseContainersTestUtilities::ComputeBounds<IndexType>(40, world_size, my_rank);
    auto el_bounds = DistributedSparseContainersTestUtilities::ComputeBounds<IndexType>(31, world_size, my_rank);
    const auto connectivities = DistributedSparseContainersTestUtilities::ElementConnectivities(el_bounds);

    DistributedSparseGraph<IndexType> A_graph(dofs_bounds[1]-dofs_bounds[0], r_comm);
    IndexPartition<IndexType>(connectivities.size()).for_each([&](IndexType i) {
        std::vector<IndexType> row_ids = connectivities[i];
        std::vector<IndexType> col_ids{row_ids[0]/col_divider, row_ids[1]/col_divider};
        A_graph.AddEntries(row_ids, col_ids);
    });
    A_graph.Finalize();

    DistributedCsrMatrix<double, IndexType> A(A_graph);

    A.BeginAssemble();
    for(const auto& c : connectivities) {
        std::vector<IndexType> row_ids = c;
        std::vector<IndexType> col_ids{row_ids[0]/col_divider, row_ids[1]/col_divider};
        Matrix data(row_ids.size(),col_ids.size(),1.0);
        A.Assemble(data,row_ids, col_ids);
    }
    A.FinalizeAssemble();

    auto A_map = A.ToMap();
    //check that all "local" values in reference_A_map also appear in A_map
    for (const auto& item : reference_A_map) {
        IndexType i = item.first.first;
        if(A.GetRowNumbering().IsLocal(i)) {
            IndexType j = item.first.second;
            double reference_v = item.second;
            if(A_map.find(item.first) == A_map.end())
                KRATOS_ERROR << "entry " << i << " " << j << "not found in A_map" <<std::endl;

            double v = A_map.find(item.first)->second;
            KRATOS_EXPECT_NEAR(v,reference_v,1e-14);
        }
    }

    //here we test SPMV by a vector of 1s
    DistributedSystemVector<> y(A_graph); //destination vector
    y.SetValue(0.0);

    DistributedSystemVector<> x(A.GetColNumbering()); //origin vector
    x.SetValue(1.0);
    A.SpMV(x,y);

    for(IndexType i_local = 0; i_local < y.LocalSize(); ++i_local) {
        IndexType i_global = y.GetNumbering().GlobalId(i_local);
        KRATOS_EXPECT_NEAR(y[i_local], yserial(i_global), 1e-14);
    }

    //check TransposeSpMV
    x.SetValue(0.0);
    y.SetValue(1.0);
    A.TransposeSpMV(y,x);

    std::vector<double> reference_transpose_spmv_res{20,8,16,20,8,32,28,4,52,16,8,24,12};
    for(unsigned int i=0; i<x.LocalSize(); ++i) {
        IndexType global_i = x.GetNumbering().GlobalId(i);
        KRATOS_EXPECT_NEAR(x[i],  reference_transpose_spmv_res[global_i], 1e-14 );
    }

}

KRATOS_TEST_CASE_IN_SUITE(DistributedCSRMatrix1DLaplacianAMGCLConstruction, KratosMPICoreFastSuite)
{
    using IndexType = DistributedSparseContainersTestUtilities::IndexType;
    const auto& r_comm = ParallelEnvironment::GetDefaultDataCommunicator();
    int world_size = r_comm.Size();
    int my_rank = r_comm.Rank();

    IndexType max_size = 3;
    auto dofs_bounds = DistributedSparseContainersTestUtilities::ComputeBounds<IndexType>(max_size, world_size, my_rank);

    Matrix local_matrix(2,2); //we will assemble a 1D laplacian
    local_matrix(0,0) = 1.0;
    local_matrix(0,1) = -1.0;
    local_matrix(1,0) = -1.0;
    local_matrix(1,1) = 1.0;
    DenseVector<IndexType> connectivities(2);

    DistributedSparseGraph<IndexType> A_graph(dofs_bounds[1]-dofs_bounds[0], r_comm);
    for (IndexType i = dofs_bounds[0]; i < dofs_bounds[1]; ++i) {
        if (i+1<max_size) {
            connectivities[0] = i;
            connectivities[1] = i+1;
            A_graph.AddEntries(connectivities);
        }
    }
    A_graph.Finalize();

    //Test SPMV
    DistributedCsrMatrix<double, IndexType> A(A_graph);

    A.BeginAssemble();
    for (IndexType i = dofs_bounds[0]; i < dofs_bounds[1]; ++i) {
        if (i+1<max_size) {
            connectivities[0] = i;
            connectivities[1] = i+1;
            A.Assemble(local_matrix,connectivities);
        }
    }
    A.FinalizeAssemble();

    bool move_to_backend = false;
    auto offdiag_global_index2 = A.GetOffDiagonalIndex2DataInGlobalNumbering();
    auto pA_amgcl = AmgclDistributedCSRConversionUtilities::ConvertToAmgcl<double,IndexType>(A, offdiag_global_index2, move_to_backend);
    auto A_reconverted = AmgclDistributedCSRConversionUtilities::ConvertToCsrMatrix<double, IndexType>(*pA_amgcl, A.GetComm());
}

KRATOS_TEST_CASE_IN_SUITE(DistributedCSRMatrixSmallRectangularMatrixMultiply, KratosMPICoreFastSuite)
{
    using IndexType = DistributedSparseContainersTestUtilities::IndexType;
    const auto& r_comm = ParallelEnvironment::GetDefaultDataCommunicator();
    int world_size =r_comm.Size();
    int my_rank = r_comm.Rank();

    // matrix A
    //[[1,0,0,7,2],
    // [0,3,0,0,0],
    // [0,0,0,7,7]]
    DistributedSparseContainersTestUtilities::MatrixMapType A_map {
        {{0,0},1.0}, {{0,3},7.0}, {{0,4},2.0},
        {{1,1},3.0},
        {{2,3},7.0}, {{2,4},7.0}
    };

    IndexType A_size1 = 3;
    auto dofs_bounds = DistributedSparseContainersTestUtilities::ComputeBounds<IndexType>(A_size1, world_size, my_rank);
    DistributedSparseGraph<IndexType> A_graph(dofs_bounds[1]-dofs_bounds[0], r_comm);
    for(const auto& item : A_map) {
        IndexType i = item.first.first;
        IndexType j = item.first.second;
        if( A_graph.GetRowNumbering().IsLocal(i)) {
            A_graph.AddEntry(i, j);
        }
    }
    A_graph.Finalize();

    DistributedCsrMatrix<double, IndexType> A(A_graph);
    A.BeginAssemble();
    for(const auto& item : A_map) {
        IndexType i = item.first.first;
        IndexType j = item.first.second;
        double value = item.second;
        if (A.GetRowNumbering().IsLocal(i)) {
            A.AssembleEntry(value, i, j);
        }
    }
    A.FinalizeAssemble();

    // matrix B
    // [[1,0,0],
    // [0,2,3],
    // [0,3,0],
    // [0,0,0],
    // [5,0,6]]
    DistributedSparseContainersTestUtilities::MatrixMapType B_map {
        {{0,0},1.0},
        {{1,1},2.0}, {{1,2},3.0},
        {{2,1},3.0},
        //empty row
        {{4,0},5.0}, {{4,2},6.0},
    };

    IndexType B_size1 = 5;
    dofs_bounds = DistributedSparseContainersTestUtilities::ComputeBounds<IndexType>(B_size1, world_size, my_rank);

    DistributedSparseGraph<IndexType> B_graph(dofs_bounds[1]-dofs_bounds[0], r_comm);
    for (const auto& item : B_map) {
        IndexType i = item.first.first;
        IndexType j = item.first.second;
        if( B_graph.GetRowNumbering().IsLocal(i)) {
            B_graph.AddEntry(i, j);
        }
    }
    B_graph.Finalize();

    DistributedCsrMatrix<double, IndexType> B(B_graph);
    B.BeginAssemble();
    for(const auto& item : B_map) {
        IndexType i = item.first.first;
        IndexType j = item.first.second;
        double value = item.second;
        if(B.GetRowNumbering().IsLocal(i)) {
            B.AssembleEntry(value, i, j);
        }
    }
    B.FinalizeAssemble();

    // //Cref = A@B
    //[[11,  0, 12],
    // [ 0,  6,  9],
    // [35,  0, 42]]
    DistributedSparseContainersTestUtilities::MatrixMapType Cref {
        {{0,0},11.0}, {{0,2},12.0},
        {{1,1},6.0}, {{1,2},9.0},
        {{2,0},35.0}, {{2,2},42.0}
    };

    auto pC = DistributedAmgclCSRSpMMUtilities::SparseMultiply<double,IndexType>(A,B); //C=A*B
    auto& C = *pC;

    for(const auto& item : Cref) {
        IndexType i = item.first.first;
        if( C.GetRowNumbering().IsLocal(i)) {
            IndexType j = item.first.second;
            double ref_value = item.second;
            KRATOS_EXPECT_EQ(ref_value, C.GetLocalDataByGlobalId(i, j));
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(DistributedCSRMatrixDiagonalValues, KratosMPICoreFastSuite)
{
    using IndexType = DistributedSparseContainersTestUtilities::IndexType;
    const auto& r_comm = ParallelEnvironment::GetDefaultDataCommunicator();
    int world_size = r_comm.Size();
    int my_rank = r_comm.Rank();

    auto dofs_bounds = DistributedSparseContainersTestUtilities::ComputeBounds<IndexType>(40, world_size, my_rank);
    auto reference_A_map = DistributedSparseContainersTestUtilities::GetReferenceMatrixAsMap(dofs_bounds);

    auto el_bounds = DistributedSparseContainersTestUtilities::ComputeBounds<IndexType>(31, world_size, my_rank);
    const auto connectivities = DistributedSparseContainersTestUtilities::ElementConnectivities(el_bounds);

    DistributedSparseGraph<IndexType> A_graph(dofs_bounds[1]-dofs_bounds[0], r_comm);
    IndexPartition<IndexType>(connectivities.size()).for_each([&](IndexType i) {
        A_graph.AddEntries(connectivities[i]);
    });
    A_graph.Finalize();

    //FEM assembly
    DistributedCsrMatrix<double,IndexType> A(A_graph);
    A.BeginAssemble();
    for(const auto& c : connectivities) {
        Matrix data(c.size(),c.size(),1.0);
        A.Assemble(data,c);
    }
    A.FinalizeAssemble();

    KRATOS_CHECK_NEAR(A.MaxDiagonal(), 6.0, 1.0e-12);
    KRATOS_CHECK_NEAR(A.MinDiagonal(), 1.0, 1.0e-12);
    KRATOS_CHECK_NEAR(A.NormDiagonal(), 22.135943621179, 1.0e-12);
}

} // namespace Kratos::Testing
