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
#include "testing/testing.h"
#include "containers/sparse_contiguous_row_graph.h"
#include "containers/sparse_graph.h"
#include "containers/system_vector.h"
#include "tests/test_utilities/sparse_containers_test_utilities.h"
#include "utilities/builtin_timer.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(SystemVectorAssembly, KratosCoreFastSuite)
{
    const auto connectivities = SparseContainersTestUtilities::ElementConnectivities();
    const auto reference_A_map = SparseContainersTestUtilities::GetReferenceMatrixAsMap();

    SparseContiguousRowGraph<> Agraph(40);
    #pragma omp parallel for
    for(int i=0; i<static_cast<int>(connectivities.size()); ++i) //note that this version is threadsafe
        Agraph.AddEntries(connectivities[i]);
    Agraph.Finalize();


    SystemVector<> b(Agraph);

    b.SetValue(0.0);
    b.BeginAssemble();
    for(const auto& c : connectivities)
    {
        Vector vdata(c.size(),1.0);
        b.Assemble(vdata,c);
    }
    b.FinalizeAssemble();

    double reference_sum = 0.0;
    for(const auto& c : connectivities){
        reference_sum += c.size();
    }

    double sum = 0.0; //for this test the final value is 124
    for(SparseContainersTestUtilities::IndexType i=0; i!=b.size(); ++i){
        sum += b(i);
    }

    KRATOS_EXPECT_EQ(sum, reference_sum);
}

KRATOS_TEST_CASE_IN_SUITE(SystemVectorOperations, KratosCoreFastSuite)
{

    SparseContainersTestUtilities::IndexType vector_size = 4;
    SystemVector<double,SparseContainersTestUtilities::IndexType> a(vector_size);
    a.SetValue(5.0);

    SystemVector<double,SparseContainersTestUtilities::IndexType> b(vector_size);
    b.SetValue(3.0);

    SystemVector<double,SparseContainersTestUtilities::IndexType> c(a);
    c += b;
    for(unsigned int i=0; i<c.size(); ++i)
        KRATOS_EXPECT_NEAR(c[i], 8.0,1e-14);

    c -= b;
    for(unsigned int i=0; i<c.size(); ++i)
        KRATOS_EXPECT_NEAR(c[i], 5.0,1e-14);

    c.Add(3.0,a);
    for(unsigned int i=0; i<c.size(); ++i)
        KRATOS_EXPECT_NEAR(c[i], 20.0,1e-14);

    c*=2.0;
    for(unsigned int i=0; i<c.size(); ++i)
        KRATOS_EXPECT_NEAR(c[i], 40.0,1e-14);

    c/=4.0;
    for(unsigned int i=0; i<c.size(); ++i)
        KRATOS_EXPECT_NEAR(c[i], 10.0,1e-14);
}

} // namespace Kratos::Testing
