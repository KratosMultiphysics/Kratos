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
#include <utility>
#include <iostream>
#include <random>

// External includes
#include <span/span.hpp>

// Project includes
#include "testing/testing.h"
#include "containers/sparse_graph.h"
#include "containers/sparse_contiguous_row_graph.h"
#include "containers/csr_matrix.h"
#include "containers/system_vector.h"
#include "includes/key_hash.h"
#include "utilities/builtin_timer.h"
#include "utilities/amgcl_csr_conversion_utilities.h"
#include "utilities/amgcl_csr_spmm_utilities.h"

namespace Kratos {
namespace Testing {
namespace SparseTestingInternals {

typedef std::size_t IndexType;
typedef std::vector<std::vector<std::size_t>> ElementConnectivityType;
typedef std::unordered_map<std::pair<std::size_t, std::size_t>,
                          double,
                          PairHasher<std::size_t, std::size_t>,
                          PairComparor<std::size_t, std::size_t>
                          > MatrixMapType;


ElementConnectivityType RandomElementConnectivities(
    const SparseTestingInternals::IndexType block_size,
    const SparseTestingInternals::IndexType nodes_in_elem,
    const SparseTestingInternals::IndexType index_begin,
    const SparseTestingInternals::IndexType index_end,
    const SparseTestingInternals::IndexType ndof,
    const SparseTestingInternals::IndexType standard_dev
)
{

    std::cout << std::endl;
    std::cout << "beginning generation" << std::endl;
    const auto timer = BuiltinTimer();
    //generating random indices
    ElementConnectivityType connectivities((index_end-index_begin)*block_size);

    #pragma omp parallel for
    for(int i=static_cast<int>(index_begin); i<static_cast<int>(index_end);++i)
    {
        connectivities[i].resize(nodes_in_elem*block_size);
        std::mt19937 gen(i);
        //std::uniform_int_distribution<> dis(0,ndof-1);
        std::normal_distribution<> dis{
            static_cast<double>(ndof/(index_end-index_begin)*i),
            static_cast<double>(standard_dev)
            };

        for(int j = 0; j<static_cast<int>(nodes_in_elem); ++j){
            //SparseTestingInternals::IndexType eq_id = dis(gen)*block_size;
            SparseTestingInternals::IndexType eq_id;
            bool acceptable = false;
            while(!acceptable){
                auto randomid = static_cast<SparseTestingInternals::IndexType>(dis(gen));
                if(static_cast<SparseTestingInternals::IndexType>(randomid) > 0 &&
                   static_cast<SparseTestingInternals::IndexType>(randomid) < ndof-1)
                {
                    acceptable=true;
                    eq_id = randomid * block_size;
                }
            }

            for(SparseTestingInternals::IndexType k = 0; k<block_size; ++k){
                connectivities[i][j*block_size+k] = eq_id+k;
            }
        }
    }
    std::cout << "finishing generation - time = " << timer.ElapsedSeconds() << std::endl;

    return connectivities;
}

ElementConnectivityType ElementConnectivities()
{
    ElementConnectivityType connectivities = {
        {19,11,7,39},
        {33,27,22,9},
        {11,2,3,6},
        {8,26,3,22},
        {0,26,5,31},
        {1,18,35,12},
        {3,36,23,7},
        {16,8,18,15},
        {16,33,10,26},
        {25,2,18,31},
        {33,26,4,6},
        {19,21,22,7},
        {9,37,29,14},
        {18,19,14,39},
        {24,34,37,7},
        {16,9,29,14},
        {17,18,11,4},
        {16,33,28,37},
        {37,26,11,5},
        {8,26,35,14},
        {24,4,30,15},
        {16,17,12,6},
        {32,25,35,28},
        {24,25,14,1},
        {24,35,5,6},
        {28,12,38,15},
        {8,18,35,6},
        {28,31,22,39},
        {1,28,13,7},
        {17,10,36,7},
        {25,14,30,9},
        };

    return connectivities;
}

MatrixMapType GetReferenceMatrixAsMap()
{
    MatrixMapType AMap{{{19,19},3.0},
        {{19,11},1.0},{{19,7},2.0},{{19,39},2.0},{{11,19},1.0},{{11,11},4.0},{{11,7},1.0},{{11,39},1.0},{{7,19},2.0},{{7,11},1.0},{{7,7},6.0},
        {{7,39},1.0},{{39,19},2.0},{{39,11},1.0},{{39,7},1.0},{{39,39},3.0},{{33,33},4.0},{{33,27},1.0},{{33,22},1.0},{{33,9},1.0},{{27,33},1.0},
        {{27,27},1.0},{{27,22},1.0},{{27,9},1.0},{{22,33},1.0},{{22,27},1.0},{{22,22},4.0},{{22,9},1.0},{{9,33},1.0},{{9,27},1.0},{{9,22},1.0},
        {{9,9},4.0},{{11,2},1.0},{{11,3},1.0},{{11,6},1.0},{{2,11},1.0},{{2,2},2.0},{{2,3},1.0},{{2,6},1.0},{{3,11},1.0},{{3,2},1.0},
        {{3,3},3.0},{{3,6},1.0},{{6,11},1.0},{{6,2},1.0},{{6,3},1.0},{{6,6},5.0},{{8,8},4.0},{{8,26},2.0},{{8,3},1.0},{{8,22},1.0},
        {{26,8},2.0},{{26,26},6.0},{{26,3},1.0},{{26,22},1.0},{{3,8},1.0},{{3,26},1.0},{{3,22},1.0},{{22,8},1.0},{{22,26},1.0},{{22,3},1.0},
        {{0,0},1.0},{{0,26},1.0},{{0,5},1.0},{{0,31},1.0},{{26,0},1.0},{{26,5},2.0},{{26,31},1.0},{{5,0},1.0},{{5,26},2.0},{{5,5},3.0},
        {{5,31},1.0},{{31,0},1.0},{{31,26},1.0},{{31,5},1.0},{{31,31},3.0},{{1,1},3.0},{{1,18},1.0},{{1,35},1.0},{{1,12},1.0},{{18,1},1.0},
        {{18,18},6.0},{{18,35},2.0},{{18,12},1.0},{{35,1},1.0},{{35,18},2.0},{{35,35},5.0},{{35,12},1.0},{{12,1},1.0},{{12,18},1.0},{{12,35},1.0},
        {{12,12},3.0},{{3,36},1.0},{{3,23},1.0},{{3,7},1.0},{{36,3},1.0},{{36,36},2.0},{{36,23},1.0},{{36,7},2.0},{{23,3},1.0},{{23,36},1.0},
        {{23,23},1.0},{{23,7},1.0},{{7,3},1.0},{{7,36},2.0},{{7,23},1.0},{{16,16},5.0},{{16,8},1.0},{{16,18},1.0},{{16,15},1.0},{{8,16},1.0},
        {{8,18},2.0},{{8,15},1.0},{{18,16},1.0},{{18,8},2.0},{{18,15},1.0},{{15,16},1.0},{{15,8},1.0},{{15,18},1.0},{{15,15},3.0},{{16,33},2.0},
        {{16,10},1.0},{{16,26},1.0},{{33,16},2.0},{{33,10},1.0},{{33,26},2.0},{{10,16},1.0},{{10,33},1.0},{{10,10},2.0},{{10,26},1.0},{{26,16},1.0},
        {{26,33},2.0},{{26,10},1.0},{{25,25},4.0},{{25,2},1.0},{{25,18},1.0},{{25,31},1.0},{{2,25},1.0},{{2,18},1.0},{{2,31},1.0},{{18,25},1.0},
        {{18,2},1.0},{{18,31},1.0},{{31,25},1.0},{{31,2},1.0},{{31,18},1.0},{{33,4},1.0},{{33,6},1.0},{{26,4},1.0},{{26,6},1.0},{{4,33},1.0},
        {{4,26},1.0},{{4,4},3.0},{{4,6},1.0},{{6,33},1.0},{{6,26},1.0},{{6,4},1.0},{{19,21},1.0},{{19,22},1.0},{{21,19},1.0},{{21,21},1.0},
        {{21,22},1.0},{{21,7},1.0},{{22,19},1.0},{{22,21},1.0},{{22,7},1.0},{{7,21},1.0},{{7,22},1.0},{{9,37},1.0},{{9,29},2.0},{{9,14},3.0},
        {{37,9},1.0},{{37,37},4.0},{{37,29},1.0},{{37,14},1.0},{{29,9},2.0},{{29,37},1.0},{{29,29},2.0},{{29,14},2.0},{{14,9},3.0},{{14,37},1.0},
        {{14,29},2.0},{{14,14},6.0},{{18,19},1.0},{{18,14},1.0},{{18,39},1.0},{{19,18},1.0},{{19,14},1.0},{{14,18},1.0},{{14,19},1.0},{{14,39},1.0},
        {{39,18},1.0},{{39,14},1.0},{{24,24},4.0},{{24,34},1.0},{{24,37},1.0},{{24,7},1.0},{{34,24},1.0},{{34,34},1.0},{{34,37},1.0},{{34,7},1.0},
        {{37,24},1.0},{{37,34},1.0},{{37,7},1.0},{{7,24},1.0},{{7,34},1.0},{{7,37},1.0},{{16,9},1.0},{{16,29},1.0},{{16,14},1.0},{{9,16},1.0},
        {{29,16},1.0},{{14,16},1.0},{{17,17},3.0},{{17,18},1.0},{{17,11},1.0},{{17,4},1.0},{{18,17},1.0},{{18,11},1.0},{{18,4},1.0},{{11,17},1.0},
        {{11,18},1.0},{{11,4},1.0},{{4,17},1.0},{{4,18},1.0},{{4,11},1.0},{{16,28},1.0},{{16,37},1.0},{{33,28},1.0},{{33,37},1.0},{{28,16},1.0},
        {{28,33},1.0},{{28,28},5.0},{{28,37},1.0},{{37,16},1.0},{{37,33},1.0},{{37,28},1.0},{{37,26},1.0},{{37,11},1.0},{{37,5},1.0},{{26,37},1.0},
        {{26,11},1.0},{{11,37},1.0},{{11,26},1.0},{{11,5},1.0},{{5,37},1.0},{{5,11},1.0},{{8,35},2.0},{{8,14},1.0},{{26,35},1.0},{{26,14},1.0},
        {{35,8},2.0},{{35,26},1.0},{{35,14},1.0},{{14,8},1.0},{{14,26},1.0},{{14,35},1.0},{{24,4},1.0},{{24,30},1.0},{{24,15},1.0},{{4,24},1.0},
        {{4,30},1.0},{{4,15},1.0},{{30,24},1.0},{{30,4},1.0},{{30,30},2.0},{{30,15},1.0},{{15,24},1.0},{{15,4},1.0},{{15,30},1.0},{{16,17},1.0},
        {{16,12},1.0},{{16,6},1.0},{{17,16},1.0},{{17,12},1.0},{{17,6},1.0},{{12,16},1.0},{{12,17},1.0},{{12,6},1.0},{{6,16},1.0},{{6,17},1.0},
        {{6,12},1.0},{{32,32},1.0},{{32,25},1.0},{{32,35},1.0},{{32,28},1.0},{{25,32},1.0},{{25,35},1.0},{{25,28},1.0},{{35,32},1.0},{{35,25},1.0},
        {{35,28},1.0},{{28,32},1.0},{{28,25},1.0},{{28,35},1.0},{{24,25},1.0},{{24,14},1.0},{{24,1},1.0},{{25,24},1.0},{{25,14},2.0},{{25,1},1.0},
        {{14,24},1.0},{{14,25},2.0},{{14,1},1.0},{{1,24},1.0},{{1,25},1.0},{{1,14},1.0},{{24,35},1.0},{{24,5},1.0},{{24,6},1.0},{{35,24},1.0},
        {{35,5},1.0},{{35,6},2.0},{{5,24},1.0},{{5,35},1.0},{{5,6},1.0},{{6,24},1.0},{{6,35},2.0},{{6,5},1.0},{{28,12},1.0},{{28,38},1.0},
        {{28,15},1.0},{{12,28},1.0},{{12,38},1.0},{{12,15},1.0},{{38,28},1.0},{{38,12},1.0},{{38,38},1.0},{{38,15},1.0},{{15,28},1.0},{{15,12},1.0},
        {{15,38},1.0},{{8,6},1.0},{{18,6},1.0},{{6,8},1.0},{{6,18},1.0},{{28,31},1.0},{{28,22},1.0},{{28,39},1.0},{{31,28},1.0},{{31,22},1.0},
        {{31,39},1.0},{{22,28},1.0},{{22,31},1.0},{{22,39},1.0},{{39,28},1.0},{{39,31},1.0},{{39,22},1.0},{{1,28},1.0},{{1,13},1.0},{{1,7},1.0},
        {{28,1},1.0},{{28,13},1.0},{{28,7},1.0},{{13,1},1.0},{{13,28},1.0},{{13,13},1.0},{{13,7},1.0},{{7,1},1.0},{{7,28},1.0},{{7,13},1.0},
        {{17,10},1.0},{{17,36},1.0},{{17,7},1.0},{{10,17},1.0},{{10,36},1.0},{{10,7},1.0},{{36,17},1.0},{{36,10},1.0},{{7,17},1.0},{{7,10},1.0},
        {{25,30},1.0},{{25,9},1.0},{{14,30},1.0},{{30,25},1.0},{{30,14},1.0},{{30,9},1.0},{{9,25},1.0},{{9,30},1.0},};

    return AMap;
}

template< class TGraphType>
bool CheckGraph(
        const TGraphType& rAgraph,
        const MatrixMapType& rReferenceGraph)
{
    for(auto it=rAgraph.begin(); it!=rAgraph.end(); ++it)
    {
        const auto I = it.GetRowIndex();
        for(auto J : *it )
        {
            if(rReferenceGraph.find({I,J}) == rReferenceGraph.end()) //implies it is not present
                KRATOS_ERROR << "Entry " << I << "," << J << "not present in A graph"  << std::endl;
        }
    }

    //check that all the entries of reference_A_map are also in Agraph
    for(const auto& item : rReferenceGraph)
    {
        auto I = item.first.first;
        auto J = item.first.second;
        if(!rAgraph.Has(I,J))
            KRATOS_ERROR << "Entry " << I << "," << J << " is in the reference graph but not in Agraph"  << std::endl;
    }

    return true;
}

bool CheckCSRGraphArrays(
    const Kratos::span<SparseTestingInternals::IndexType>& rRowIndices,
    const Kratos::span<SparseTestingInternals::IndexType>& rColIndices,
    const MatrixMapType& rReferenceGraph)
{
    auto N = rRowIndices.size()-1;
    for (SparseTestingInternals::IndexType I = 0; I < N; ++I)
    {
        for (auto k = rRowIndices[I]; k < rRowIndices[I + 1]; ++k)
        {
            SparseTestingInternals::IndexType J = rColIndices[k];
            if (rReferenceGraph.find({I, J}) == rReferenceGraph.end()) //implies it is not present
                KRATOS_ERROR << "Entry " << I << "," << J << "not present in A graph" << std::endl;

            //check that that cols are ordered
            if(k-rRowIndices[I] > 0)
                if(rColIndices[k-1]>rColIndices[k])
                    KRATOS_ERROR << "columns are not ordered in csr" << std::endl;
        }
    }
    return true;
}

bool CheckCSRMatrix(
    const CsrMatrix<>& rA,
    const MatrixMapType& rReferenceGraph)
{
    const auto& rRowIndices = rA.index1_data();
    const auto& rColIndices = rA.index2_data();
    const auto& rValues = rA.value_data();

    auto N = rRowIndices.size()-1;

    KRATOS_EXPECT_EQ(rReferenceGraph.size(), rA.nnz());

    for (SparseTestingInternals::IndexType I = 0; I < N; ++I)
    {
        for (auto k = rRowIndices[I]; k < rRowIndices[I + 1]; ++k)
        {
            SparseTestingInternals::IndexType J = rColIndices[k];
            double value = rValues[k];
            if (rReferenceGraph.find({I, J}) == rReferenceGraph.end()) //implies it is not present
                KRATOS_ERROR << "Entry " << I << "," << J << "not present in A graph" << std::endl;

            KRATOS_EXPECT_NEAR(rReferenceGraph.find({I, J})->second , value, 1e-14);

            //check that that cols are ordered
            if(k-rRowIndices[I] > 0)
                if(rColIndices[k-1]>rColIndices[k])
                    KRATOS_ERROR << "columns are not ordered in csr" << std::endl;
        }
    }
    return true;
}



template<typename TGraphType>
std::unique_ptr<TGraphType> AssembleGraph(
    ElementConnectivityType& rConnectivities,
    SparseTestingInternals::IndexType GraphSize)
{
    std::unique_ptr<TGraphType> pAgraph(nullptr);

    #pragma omp parallel
    {
        std::unique_ptr<TGraphType>
            plocal_graph(new TGraphType(GraphSize));
        #pragma omp for
        for(int i=0; i<static_cast<int>(rConnectivities.size()); ++i){
            plocal_graph->AddEntries(rConnectivities[i]);
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
    return pAgraph;
}

}

// Basic Type
KRATOS_TEST_CASE_IN_SUITE(GraphConstruction, KratosCoreFastSuite)
{
    const auto connectivities = SparseTestingInternals::ElementConnectivities();
    auto reference_A_map = SparseTestingInternals::GetReferenceMatrixAsMap();

    SparseGraph<> Agraph;
    for(const auto& c : connectivities)
        Agraph.AddEntries(c);

    Agraph.Finalize();

    SparseTestingInternals::CheckGraph(Agraph, reference_A_map);

    //check serialization
    StreamSerializer serializer;
    const std::string tag_string("TestString");
    serializer.save(tag_string, Agraph);
    Agraph.Clear();
    serializer.load(tag_string, Agraph);
    SparseTestingInternals::CheckGraph(Agraph, reference_A_map);

    //check exporting
    std::vector<IndexType> row_index;
    std::vector<IndexType> col_index;
    Agraph.ExportCSRArrays(row_index, col_index);

    SparseTestingInternals::CheckCSRGraphArrays(row_index, col_index, reference_A_map);

}

KRATOS_TEST_CASE_IN_SUITE(GraphContiguousRowConstruction, KratosCoreFastSuite)
{
    const auto connectivities = SparseTestingInternals::ElementConnectivities();
    auto reference_A_map = SparseTestingInternals::GetReferenceMatrixAsMap();

    SparseContiguousRowGraph<> Agraph(40);
    for(const auto& c : connectivities)
        Agraph.AddEntries(c);
    Agraph.Finalize();

    SparseTestingInternals::CheckGraph(Agraph, reference_A_map);

    //check serialization
    StreamSerializer serializer;
    const std::string tag_string("TestString");
    serializer.save(tag_string, Agraph);
    Agraph.Clear();
    serializer.load(tag_string, Agraph);
    SparseTestingInternals::CheckGraph(Agraph, reference_A_map);

    // //check exporting
    // vector<SparseTestingInternals::IndexType> row_index, col_index;
    // Agraph.ExportCSRArrays(row_index, col_index);

    // CheckCSRGraphArrays(row_index, col_index, reference_A_map);

}

KRATOS_TEST_CASE_IN_SUITE(CSRConstruction, KratosCoreFastSuite)
{
    
    const auto connectivities = SparseTestingInternals::ElementConnectivities();
    auto reference_A_map = SparseTestingInternals::GetReferenceMatrixAsMap();

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

    SparseTestingInternals::CheckCSRMatrix(A, reference_A_map);



}

// Basic Type
KRATOS_TEST_CASE_IN_SUITE(OpenMPGraphContiguousRowConstruction, KratosCoreFastSuite)
{
    
    const auto connectivities = SparseTestingInternals::ElementConnectivities();
    auto reference_A_map = SparseTestingInternals::GetReferenceMatrixAsMap();

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

    SparseTestingInternals::CheckGraph(*pAgraph, reference_A_map);
}

KRATOS_TEST_CASE_IN_SUITE(PerformanceBenchmarkSparseGraph, KratosCoreFastSuite)
{
    const SparseTestingInternals::IndexType block_size = 4;
    const SparseTestingInternals::IndexType nodes_in_elem = 4;
    const SparseTestingInternals::IndexType nel = 1e2; //set at least to 1e6 for a realistic benchmark
    const SparseTestingInternals::IndexType ndof = nel/6;
    const SparseTestingInternals::IndexType standard_dev = 100;
    auto connectivities = SparseTestingInternals::RandomElementConnectivities(block_size,nodes_in_elem, 0,nel,ndof,standard_dev);

    const auto timer = BuiltinTimer();
    auto pAgraph = SparseTestingInternals::AssembleGraph<SparseGraph<>>(connectivities, ndof*block_size);

    std::cout << "SparseGraph time = " << timer.ElapsedSeconds() << std::endl;
}

KRATOS_TEST_CASE_IN_SUITE(PerformanceBenchmarkSparseContiguousRowGraph, KratosCoreFastSuite)
{
    
    const SparseTestingInternals::IndexType block_size = 4;
    const SparseTestingInternals::IndexType nodes_in_elem = 4;
    const SparseTestingInternals::IndexType nel = 1e2; //set this to at least 1e6 for a realistic benchmark
    const SparseTestingInternals::IndexType ndof = nel/6;
    const SparseTestingInternals::IndexType standard_dev = 100;
    auto connectivities = SparseTestingInternals::RandomElementConnectivities(block_size,nodes_in_elem,0, nel,ndof,standard_dev);

    const auto timer = BuiltinTimer();

    SparseContiguousRowGraph<> Agraph(ndof*block_size);
    #pragma omp parallel for
    for(int i=0; i<static_cast<int>(connectivities.size()); ++i) //note that this version is threadsafe
        Agraph.AddEntries(connectivities[i]);
    Agraph.Finalize();

    std::cout << "SparseGraphContiguousRow generation time = " << timer.ElapsedSeconds() << std::endl;
}

KRATOS_TEST_CASE_IN_SUITE(SystemVectorAssembly, KratosCoreFastSuite)
{
    
    const auto connectivities = SparseTestingInternals::ElementConnectivities();
    const auto reference_A_map = SparseTestingInternals::GetReferenceMatrixAsMap();


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
    for(SparseTestingInternals::IndexType i=0; i!=b.size(); ++i){
        sum += b(i);
    }

    KRATOS_EXPECT_EQ(sum, reference_sum);
}


KRATOS_TEST_CASE_IN_SUITE(SpMV, KratosCoreFastSuite)
{
    
    const auto connectivities = SparseTestingInternals::ElementConnectivities();
    auto reference_A_map = SparseTestingInternals::GetReferenceMatrixAsMap();

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
    for(SparseTestingInternals::IndexType i=0; i!=x.size(); ++i){
        sum += x(i);
    }

    double reference_sum = 0.0;
    for(auto& item : reference_A_map)
        reference_sum += item.second;

    KRATOS_EXPECT_EQ(sum, reference_sum);

}

KRATOS_TEST_CASE_IN_SUITE(SystemVectorOperations, KratosCoreFastSuite)
{
    
    SparseTestingInternals::IndexType vector_size = 4;
    SystemVector<double,SparseTestingInternals::IndexType> a(vector_size);
    a.SetValue(5.0);

    SystemVector<double,SparseTestingInternals::IndexType> b(vector_size);
    b.SetValue(3.0);

    SystemVector<double,SparseTestingInternals::IndexType> c(a);
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

KRATOS_TEST_CASE_IN_SUITE(ToAMGCLMatrix, KratosCoreFastSuite)
{
    const auto connectivities = SparseTestingInternals::ElementConnectivities();
    auto reference_A_map = SparseTestingInternals::GetReferenceMatrixAsMap();

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
    SparseTestingInternals::MatrixMapType Cref{
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
    
    SparseTestingInternals::IndexType col_divider = 3;

    //*************************************************************************
    //compute reference solution - serial mode
    const auto all_connectivities = SparseTestingInternals::ElementConnectivities();
    SparseContiguousRowGraph<SparseTestingInternals::IndexType> Agraph(40);

    IndexPartition<SparseTestingInternals::IndexType>(all_connectivities.size()).for_each([&](SparseTestingInternals::IndexType i)    {
        std::vector<SparseTestingInternals::IndexType> row_ids = all_connectivities[i];
        std::vector<SparseTestingInternals::IndexType> col_ids{row_ids[0]/col_divider, row_ids[1]/col_divider};
        Agraph.AddEntries(row_ids, col_ids);
    });
    Agraph.Finalize();

    CsrMatrix<double, SparseTestingInternals::IndexType> A(Agraph);

    A.BeginAssemble();   
    for(const auto& c : all_connectivities){   
        std::vector<SparseTestingInternals::IndexType> row_ids = c;
        std::vector<SparseTestingInternals::IndexType> col_ids{row_ids[0]/col_divider, row_ids[1]/col_divider};
        Matrix data(row_ids.size(),col_ids.size(),1.0);
        A.Assemble(data,row_ids, col_ids);
    }
    A.FinalizeAssemble();

    auto reference_A_map = A.ToMap();

    //constructing transpose Map to later verify TransposeSPMV
    SparseTestingInternals::MatrixMapType reference_transpose_A_map;
    for(const auto& item : reference_A_map)
    {
        const SparseTestingInternals::IndexType I = item.first.first;
        const SparseTestingInternals::IndexType J = item.first.second;
        const SparseTestingInternals::IndexType value = item.second;
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
        const SparseTestingInternals::IndexType I = item.first.first;
        const SparseTestingInternals::IndexType J = item.first.second;
        const SparseTestingInternals::IndexType Atranspose_ij = item.second;
        xtranspose_spmv_ref[I] += Atranspose_ij*y[J];
    }

    for(SparseTestingInternals::IndexType i=0; i<x.size(); ++i)
    {
        KRATOS_EXPECT_NEAR(x[i], xtranspose_spmv_ref[i], 1e-14);
    }

}



} // namespace Testing
} // namespace Kratos
