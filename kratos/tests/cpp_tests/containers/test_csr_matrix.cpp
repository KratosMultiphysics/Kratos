//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi

// System includes

// External includes

// Project includes
#include "containers/csr_matrix.h"
#include "containers/sparse_contiguous_row_graph.h"
#include "containers/sparse_graph.h"
#include "containers/system_vector.h"
#include "utilities/amgcl_csr_conversion_utilities.h"
#include "utilities/amgcl_csr_spmm_utilities.h"
#include "tests/test_utilities/sparse_containers_test_utilities.h"
#include "testing/testing.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CSRConstruction, KratosCoreFastSuite)
{
    const auto connectivities = SparseContainersTestUtilities::ElementConnectivities();
    auto reference_A_map = SparseContainersTestUtilities::GetReferenceMatrixAsMap();

    SparseContiguousRowGraph<> Agraph(40);
    for(const auto& c : connectivities)
        Agraph.AddEntries(c);
    Agraph.Finalize();

    CsrMatrix<double> A(Agraph);

    A.BeginAssemble();
    for(const auto& c : connectivities)
    {
        Matrix data(c.size(),c.size(),1.0);
        A.Assemble(data,c);
    }
    A.FinalizeAssemble();

    SparseContainersTestUtilities::CheckCSRMatrix(A, reference_A_map);
}

KRATOS_TEST_CASE_IN_SUITE(SpMV, KratosCoreFastSuite)
{
    const auto connectivities = SparseContainersTestUtilities::ElementConnectivities();
    auto reference_A_map = SparseContainersTestUtilities::GetReferenceMatrixAsMap();

    SparseContiguousRowGraph<> Agraph(40);
    for(const auto& c : connectivities)
        Agraph.AddEntries(c);
    Agraph.Finalize();

    CsrMatrix<double> A(Agraph);

    A.BeginAssemble();
    for(const auto& c : connectivities){
        Matrix data(c.size(),c.size(),1.0);
        A.Assemble(data,c);
    }
    A.FinalizeAssemble();

    SystemVector<> x(Agraph);
    x.SetValue(0.0);

    SystemVector<> y(Agraph);
    y.SetValue(1.0);

    A.SpMV(y,x); //x += A*y

    double sum = 0.0; //for this test the final value is 124
    for(SparseContainersTestUtilities::IndexType i=0; i!=x.size(); ++i){
        sum += x(i);
    }

    double reference_sum = 0.0;
    for(auto& item : reference_A_map)
        reference_sum += item.second;

    KRATOS_EXPECT_EQ(sum, reference_sum);

}

KRATOS_TEST_CASE_IN_SUITE(ToAMGCLMatrix, KratosCoreFastSuite)
{
    const auto connectivities = SparseContainersTestUtilities::ElementConnectivities();
    auto reference_A_map = SparseContainersTestUtilities::GetReferenceMatrixAsMap();

    SparseContiguousRowGraph<> Agraph(40);
    for(const auto& c : connectivities)
        Agraph.AddEntries(c);
    Agraph.Finalize();

    CsrMatrix<double> A(Agraph);

    A.BeginAssemble();
    for(const auto& c : connectivities){
        Matrix data(c.size(),c.size(),1.0);
        A.Assemble(data,c);
    }
    A.FinalizeAssemble();

    auto pAmgcl = AmgclCSRConversionUtilities::ConvertToAmgcl(A);

    std::vector<double> x(A.size1(), 1.0);

    std::vector<double> y(A.size1(), 0.0);

    amgcl::backend::spmv(1.0,*pAmgcl, x, 1.0, y);

    double sum=0;
    for(auto& item : y)
       sum += item;

    KRATOS_EXPECT_EQ(sum,496);

    auto pAconverted = AmgclCSRConversionUtilities::ConvertToCsrMatrix<double,IndexType>(*pAmgcl); //NOTE that A,Aconverted and pAmgcl all have the same data!
    auto reference_map = A.ToMap();
    auto converted_A_map = pAconverted->ToMap();
    for(const auto& item : reference_map)
        KRATOS_EXPECT_EQ(item.second, converted_A_map[item.first]);
    for(const auto& item : converted_A_map)
        KRATOS_EXPECT_EQ(item.second, reference_map[item.first]);

    //matrix matrix multiplication
    CsrMatrix<double>::Pointer pC = AmgclCSRSpMMUtilities::SparseMultiply(A,*pAconverted); //C=A*Aconverted
}

KRATOS_TEST_CASE_IN_SUITE(SmallRectangularMatricMatrixMultiply, KratosCoreFastSuite)
{
    // matrix A
    //[[1,0,0,7,2],
    // [0,3,0,0,0],
    // [0,0,0,7,7]]
    SparseContiguousRowGraph<> Agraph(3);
    Agraph.AddEntry(0,0); Agraph.AddEntry(0,3);  Agraph.AddEntry(0,4);
    Agraph.AddEntry(1,1);
    Agraph.AddEntry(2,3);  Agraph.AddEntry(2,4);
    Agraph.Finalize();

    CsrMatrix<double> A(Agraph);
    A.BeginAssemble();
    A.AssembleEntry(1.0,0,0); A.AssembleEntry(7.0,0,3); A.AssembleEntry(2.0,0,4);
    A.AssembleEntry(3.0,1,1);
    A.AssembleEntry(7.0,2,3); A.AssembleEntry(7.0,2,4);
    A.FinalizeAssemble();

    // matrix B
    //[[1,0,0],
    // [0,2,3],
    // [0,3,0],
    // [0,0,0],
    // [5,0,6]]
    SparseContiguousRowGraph<> Bgraph(5);
    Bgraph.AddEntry(0,0);
    Bgraph.AddEntry(1,1); Bgraph.AddEntry(1,2);
    Bgraph.AddEntry(2,1);
    Bgraph.AddEntry(3,0);
    Bgraph.AddEntry(4,0); Bgraph.AddEntry(4,2);
    Bgraph.Finalize();

    CsrMatrix<double> B(Bgraph);

    B.BeginAssemble();
    B.AssembleEntry(1.0,0,0);
    B.AssembleEntry(2.0,1,1); B.AssembleEntry(3.0,1,2);
    B.AssembleEntry(3.0,2,1);
    B.AssembleEntry(0.0,3,0);
    B.AssembleEntry(5.0,4,0);  B.AssembleEntry(6.0,4,2);
    B.FinalizeAssemble();

    //Cref = A@B
    //[[11,  0, 12],
    // [ 0,  6,  9],
    // [35,  0, 42]]
    SparseContainersTestUtilities::MatrixMapType Cref{
        {{0,0},11.0}, {{0,2},12.0},
        {{1,1},6.0}, {{1,2},9.0},
        {{2,0},35.0}, {{2,2},42.0}
        };

    CsrMatrix<double>::Pointer pC = AmgclCSRSpMMUtilities::SparseMultiply(A,B); //C=A*B
    auto& C = *pC;

    for(const auto& item : Cref)
    {
        IndexType I = item.first.first;
        IndexType J = item.first.second;
        double ref_value = item.second;
        KRATOS_EXPECT_EQ(ref_value,C(I,J));
    }
}

KRATOS_TEST_CASE_IN_SUITE(RectangularMatrixConstruction, KratosCoreFastSuite)
{
    SparseContainersTestUtilities::IndexType col_divider = 3;

    //*************************************************************************
    //compute reference solution - serial mode
    const auto all_connectivities = SparseContainersTestUtilities::ElementConnectivities();
    SparseContiguousRowGraph<SparseContainersTestUtilities::IndexType> Agraph(40);

    IndexPartition<SparseContainersTestUtilities::IndexType>(all_connectivities.size()).for_each([&](SparseContainersTestUtilities::IndexType i)    {
        std::vector<SparseContainersTestUtilities::IndexType> row_ids = all_connectivities[i];
        std::vector<SparseContainersTestUtilities::IndexType> col_ids{row_ids[0]/col_divider, row_ids[1]/col_divider};
        Agraph.AddEntries(row_ids, col_ids);
    });
    Agraph.Finalize();

    CsrMatrix<double, SparseContainersTestUtilities::IndexType> A(Agraph);

    A.BeginAssemble();
    for(const auto& c : all_connectivities){
        std::vector<SparseContainersTestUtilities::IndexType> row_ids = c;
        std::vector<SparseContainersTestUtilities::IndexType> col_ids{row_ids[0]/col_divider, row_ids[1]/col_divider};
        Matrix data(row_ids.size(),col_ids.size(),1.0);
        A.Assemble(data,row_ids, col_ids);
    }
    A.FinalizeAssemble();

    auto reference_A_map = A.ToMap();

    //constructing transpose Map to later verify TransposeSPMV
    SparseContainersTestUtilities::MatrixMapType reference_transpose_A_map;
    for(const auto& item : reference_A_map)
    {
        const SparseContainersTestUtilities::IndexType I = item.first.first;
        const SparseContainersTestUtilities::IndexType J = item.first.second;
        const SparseContainersTestUtilities::IndexType value = item.second;
        reference_transpose_A_map[{J,I}] = value;
    }

    // //here we test SPMV by a vector of 1s
    SystemVector<> y(A.size1()); //destination vector
    y.SetValue(0.0);

    SystemVector<> x(A.size2()); //origin vector
    x.SetValue(1.0);

    A.SpMV(x,y);

    //test TransposeSpMV
    y.SetValue(1.0);
    x.SetValue(0.0);
    A.TransposeSpMV(y,x); //x += Atranspose*y

    SystemVector<> xtranspose_spmv_ref(A.size2()); //origin vector
    xtranspose_spmv_ref.SetValue(0.0);
    for(const auto& item : reference_transpose_A_map)
    {
        const SparseContainersTestUtilities::IndexType I = item.first.first;
        const SparseContainersTestUtilities::IndexType J = item.first.second;
        const SparseContainersTestUtilities::IndexType Atranspose_ij = item.second;
        xtranspose_spmv_ref[I] += Atranspose_ij*y[J];
    }

    for(SparseContainersTestUtilities::IndexType i=0; i<x.size(); ++i)
    {
        KRATOS_EXPECT_NEAR(x[i], xtranspose_spmv_ref[i], 1e-14);
    }
}

KRATOS_TEST_CASE_IN_SUITE(CsrMatrixDiagonalValues, KratosCoreFastSuite)
{
    const auto connectivities = SparseContainersTestUtilities::ElementConnectivities();

    SparseContiguousRowGraph<> Agraph(40);
    for(const auto& c : connectivities)
        Agraph.AddEntries(c);
    Agraph.Finalize();

    CsrMatrix<double> A(Agraph);

    A.BeginAssemble();
    for(const auto& c : connectivities){
        Matrix data(c.size(),c.size(),1.0);
        A.Assemble(data,c);
    }
    A.FinalizeAssemble();

    KRATOS_CHECK_NEAR(A.MaxDiagonal(), 6.0, 1.0e-12);
    KRATOS_CHECK_NEAR(A.MinDiagonal(), 1.0, 1.0e-12);
    KRATOS_CHECK_NEAR(A.NormDiagonal(), 22.135943621179, 1.0e-12);
}

} // namespace Kratos::Testing
