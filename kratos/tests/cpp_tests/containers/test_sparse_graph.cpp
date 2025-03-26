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
#include "containers/sparse_graph.h"
#include "containers/sparse_contiguous_row_graph.h"
#include "testing/testing.h"
#include "tests/test_utilities/sparse_containers_test_utilities.h"
#include "utilities/builtin_timer.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GraphConstruction, KratosCoreFastSuite)
{
    const auto connectivities = SparseContainersTestUtilities::ElementConnectivities();
    auto reference_A_map = SparseContainersTestUtilities::GetReferenceMatrixAsMap();

    SparseGraph<> Agraph;
    for(const auto& c : connectivities)
        Agraph.AddEntries(c);
    Agraph.Finalize();

    SparseContainersTestUtilities::CheckGraph(Agraph, reference_A_map);

    //check serialization
    StreamSerializer serializer;
    const std::string tag_string("TestString");
    serializer.save(tag_string, Agraph);
    Agraph.Clear();
    serializer.load(tag_string, Agraph);
    SparseContainersTestUtilities::CheckGraph(Agraph, reference_A_map);

    //check exporting
    std::vector<IndexType> row_index;
    std::vector<IndexType> col_index;
    Agraph.ExportCSRArrays(row_index, col_index);

    SparseContainersTestUtilities::CheckCSRGraphArrays(row_index, col_index, reference_A_map);
}

KRATOS_TEST_CASE_IN_SUITE(GraphContiguousRowConstruction, KratosCoreFastSuite)
{
    const auto connectivities = SparseContainersTestUtilities::ElementConnectivities();
    auto reference_A_map = SparseContainersTestUtilities::GetReferenceMatrixAsMap();

    SparseContiguousRowGraph<> Agraph(40);
    for(const auto& c : connectivities)
        Agraph.AddEntries(c);
    Agraph.Finalize();

    SparseContainersTestUtilities::CheckGraph(Agraph, reference_A_map);

    //check serialization
    StreamSerializer serializer;
    const std::string tag_string("TestString");
    serializer.save(tag_string, Agraph);
    Agraph.Clear();
    serializer.load(tag_string, Agraph);
    SparseContainersTestUtilities::CheckGraph(Agraph, reference_A_map);

    // //check exporting
    // vector<SparseContainersTestUtilities::IndexType> row_index, col_index;
    // Agraph.ExportCSRArrays(row_index, col_index);

    // CheckCSRGraphArrays(row_index, col_index, reference_A_map);
}

// Basic Type
KRATOS_TEST_CASE_IN_SUITE(OpenMPGraphContiguousRowConstruction, KratosCoreFastSuite)
{

    const auto connectivities = SparseContainersTestUtilities::ElementConnectivities();
    auto reference_A_map = SparseContainersTestUtilities::GetReferenceMatrixAsMap();

    std::unique_ptr<SparseContiguousRowGraph<>> pAgraph(nullptr);

    #pragma omp parallel
    {
        std::unique_ptr<SparseContiguousRowGraph<>> plocal_graph(new SparseContiguousRowGraph<>(40));
        #pragma omp for
        for(int i=0; i<static_cast<int>(connectivities.size()); ++i){
            plocal_graph->AddEntries(connectivities[i]);
        }

        #pragma omp critical
        {
            if(pAgraph == nullptr )
                pAgraph.swap(plocal_graph);
            else
                pAgraph->AddEntries(*plocal_graph);
        }
    }

    pAgraph->Finalize();

    SparseContainersTestUtilities::CheckGraph(*pAgraph, reference_A_map);
}

KRATOS_TEST_CASE_IN_SUITE(PerformanceBenchmarkSparseGraph, KratosCoreFastSuite)
{
    const SparseContainersTestUtilities::IndexType block_size = 4;
    const SparseContainersTestUtilities::IndexType nodes_in_elem = 4;
    const SparseContainersTestUtilities::IndexType nel = 1e2; //set at least to 1e6 for a realistic benchmark
    const SparseContainersTestUtilities::IndexType ndof = nel/6;
    const SparseContainersTestUtilities::IndexType standard_dev = 100;
    auto connectivities = SparseContainersTestUtilities::RandomElementConnectivities(block_size,nodes_in_elem, 0,nel,ndof,standard_dev);

    const auto timer = BuiltinTimer();

    auto pAgraph = SparseContainersTestUtilities::AssembleGraph<SparseGraph<>>(connectivities, ndof*block_size);

    // std::cout << "SparseGraph time = " << timer.ElapsedSeconds() << std::endl;
}

KRATOS_TEST_CASE_IN_SUITE(PerformanceBenchmarkSparseContiguousRowGraph, KratosCoreFastSuite)
{

    const SparseContainersTestUtilities::IndexType block_size = 4;
    const SparseContainersTestUtilities::IndexType nodes_in_elem = 4;
    const SparseContainersTestUtilities::IndexType nel = 1e2; //set this to at least 1e6 for a realistic benchmark
    const SparseContainersTestUtilities::IndexType ndof = nel/6;
    const SparseContainersTestUtilities::IndexType standard_dev = 100;
    auto connectivities = SparseContainersTestUtilities::RandomElementConnectivities(block_size,nodes_in_elem,0, nel,ndof,standard_dev);

    const auto timer = BuiltinTimer();

    SparseContiguousRowGraph<> Agraph(ndof*block_size);
    #pragma omp parallel for
    for(int i=0; i<static_cast<int>(connectivities.size()); ++i) //note that this version is threadsafe
        Agraph.AddEntries(connectivities[i]);
    Agraph.Finalize();

    // std::cout << "SparseGraphContiguousRow generation time = " << timer.ElapsedSeconds() << std::endl;
}

} // namespace Kratos::Testing
