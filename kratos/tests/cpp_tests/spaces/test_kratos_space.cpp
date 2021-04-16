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
//                   Philipp Bucher (https://github.com/philbucher)
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "testing/testing.h"
#include "spaces/kratos_space.h"


namespace Kratos {
namespace Testing {

namespace SparseSpaceTestingInternals {
    typedef std::size_t IndexType;
    typedef KratosSpace<double, CsrMatrix<double,IndexType> , SystemVector<double,IndexType>> SparseSpaceType;
    typedef KratosSpace<double, DenseMatrix<double>, DenseVector<double>> LocalSpaceType;

}


KRATOS_TEST_CASE_IN_SUITE(KratosSpaceNormSparseMatrix, KratosCoreFastSuite)
{
    const SparseSpaceTestingInternals::IndexType size = 10;

    SparseContiguousRowGraph<IndexType> graph(size);
    for (IndexType i=0; i<size; ++i)	{
        if (i>=1) graph.AddEntry(i, i-1);
        graph.AddEntry(i, i);
        if (i+1<size) graph.AddEntry(i, i+1);
    }
    graph.Finalize();

    SparseSpaceTestingInternals::SparseSpaceType::MatrixType mat(graph);
    mat.BeginAssemble();
    for (IndexType i=0; i<mat.size1(); ++i)	{
        if (i>=1) mat.Assemble( -1.123, i, i-1);
        mat.Assemble(4.5,i, i );
        if (i+1<mat.size2()) {mat.Assemble( 2.336, i, i+1);}
    }
    mat.FinalizeAssemble();
    

    KRATOS_CHECK_NEAR(16.216110045260546, SparseSpaceTestingInternals::SparseSpaceType::TwoNorm(mat), 1e-12);
    KRATOS_CHECK_NEAR(31.131, SparseSpaceTestingInternals::SparseSpaceType::JacobiNorm(mat), 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(KratosSpaceNormDenseMatrix, KratosCoreFastSuite)
{
    const SparseSpaceTestingInternals::IndexType size = 10;
    SparseSpaceTestingInternals::LocalSpaceType::MatrixType mat(size, size, 0.0);

    for (std::size_t i=0; i<mat.size1(); ++i)	{
        mat(i,i) = 4.5;

        if (i>=1) {mat(i,i-1) = -1.123;}

        if (i+1<mat.size2()) {mat(i,i+1) = 2.336;}
    }

    KRATOS_CHECK_NEAR(16.216110045260546, SparseSpaceTestingInternals::LocalSpaceType::TwoNorm(mat), 1e-12);
    KRATOS_CHECK_NEAR(31.131, SparseSpaceTestingInternals::LocalSpaceType::JacobiNorm(mat), 1e-12);
}

} // namespace Testing
} // namespace Kratos.
