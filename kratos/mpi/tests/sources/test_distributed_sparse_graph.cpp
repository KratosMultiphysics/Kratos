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

// Project includes
#include "testing/testing.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "containers/sparse_contiguous_row_graph.h"
#include "containers/distributed_sparse_graph.h"
#include "containers/distributed_csr_matrix.h"
#include "containers/distributed_system_vector.h"
#include "containers/distributed_vector_importer.h"
#include "includes/key_hash.h"
#include "utilities/openmp_utils.h"
#include "utilities/builtin_timer.h"

#include "mpi/utilities/amgcl_distributed_csr_conversion_utilities.h"
#include "mpi/utilities/amgcl_distributed_csr_spmm_utilities.h"

namespace Kratos
{
namespace Testing
{

namespace DistTestingInternals
{

typedef std::size_t IndexType;
typedef std::vector<std::vector<DistTestingInternals::IndexType>> ElementConnectivityType;
typedef std::unordered_map<std::pair<DistTestingInternals::IndexType, DistTestingInternals::IndexType>,
        double,
        PairHasher<DistTestingInternals::IndexType, DistTestingInternals::IndexType>,
        PairComparor<DistTestingInternals::IndexType, DistTestingInternals::IndexType>
        > MatrixMapType;

ElementConnectivityType ElementConnectivities()
{
    ElementConnectivityType connectivities =
    {
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

ElementConnectivityType ElementConnectivities(const std::vector<DistTestingInternals::IndexType> bounds)
{
    KRATOS_TRY
    ElementConnectivityType all_connectivities = ElementConnectivities();

    ElementConnectivityType connectivities;
    for(DistTestingInternals::IndexType i=bounds[0]; i<bounds[1]; ++i)
    {
        connectivities.push_back(all_connectivities[i]);
    }

    return connectivities;
    KRATOS_CATCH("")
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

MatrixMapType GetReferenceMatrixAsMap(const std::vector<DistTestingInternals::IndexType>& bounds)
{
    KRATOS_TRY
    MatrixMapType all_connectivities = GetReferenceMatrixAsMap();
    if(bounds[1]>all_connectivities.size())
        KRATOS_ERROR << "bounds : " << bounds << " exceed the total size : "
                     << all_connectivities.size() << std::endl;
    MatrixMapType connectivities;
    for(auto& item : all_connectivities)
    {
        if(item.first.first >= bounds[0] && item.first.first<bounds[1])
            connectivities.insert(item);
    }

    return connectivities;
    KRATOS_CATCH("")
}


std::map<DistTestingInternals::IndexType, double> GetReferencebAsMap(const std::vector<DistTestingInternals::IndexType>& bounds)
{
    KRATOS_TRY
    MatrixMapType all_connectivities = GetReferenceMatrixAsMap();
    if(bounds[1]>all_connectivities.size())
        KRATOS_ERROR << "bounds : " << bounds << " exceed the total size : "
                     << all_connectivities.size() << std::endl;

    std::vector<double> reference_b_vector{1,3,2,3,3,3,5,6,4,4,2,4,3,1,6,3,5,3,6,3,0,1,4,1,4,4,6,1,5,2,2,3,1,4,1,5,2,4,1,3};

    std::map<DistTestingInternals::IndexType, double> reference_b_map;

    for(unsigned int i=0; i<reference_b_vector.size(); ++i)
    {
        if(i>=bounds[0] && i<bounds[1])
            reference_b_map.insert({i,reference_b_vector[i]});
    }

    return reference_b_map;
    KRATOS_CATCH("")
}



template< class TSparseGraphType>
bool CheckGraph(
    const TSparseGraphType& rAgraph,
    const MatrixMapType& rReferenceGraph)
{
    //check that all entries in Agraph are also in reference_A_map
    for(DistTestingInternals::IndexType local_i = 0; local_i<rAgraph.LocalSize(); ++local_i) //i is the LOCAL index
    {
        auto I = rAgraph.GetRowNumbering().GlobalId(local_i);
        for(auto J : rAgraph[local_i] )
        {
            if(rReferenceGraph.find({I,J}) == rReferenceGraph.end()) //implies it is not present
                KRATOS_ERROR << "Entry " << I << "," << J << "not present in A graph"  << std::endl;
        }
    }

    //check that all the entries of reference_A_map are also in Agraph
    for(auto item : rReferenceGraph)
    {
        auto GlobalI = item.first.first;
        auto GlobalJ = item.first.second;
        if(!rAgraph.Has(GlobalI,GlobalJ))
            KRATOS_ERROR << "Entry " << GlobalI << "," << GlobalJ << " is in the reference graph but not in Agraph"  << std::endl;

    }

    return true;
}

ElementConnectivityType RandomElementConnectivities(
    const DistTestingInternals::IndexType block_size,
    const DistTestingInternals::IndexType nodes_in_elem,
    const DistTestingInternals::IndexType index_begin,
    const DistTestingInternals::IndexType index_end,
    const DistTestingInternals::IndexType ndof,
    const DistTestingInternals::IndexType standard_dev
)
{

    std::cout << std::endl;
    std::cout << "beginning generation" << std::endl;
    const auto timer = BuiltinTimer();
    //generating random indices
    ElementConnectivityType connectivities((index_end-index_begin)*block_size);

    IndexPartition<std::size_t>(connectivities.size()).for_each([&](std::size_t Index){
        connectivities[Index].resize(nodes_in_elem*block_size);
        std::mt19937 gen(Index);
        std::normal_distribution<> dis{
            static_cast<double>(ndof/(index_end-index_begin)* Index),
            static_cast<double>(standard_dev)
        };

        for(int j = 0; j<static_cast<int>(nodes_in_elem); ++j)
        {
            //DistTestingInternals::IndexType eq_id = dis(gen)*block_size;
            DistTestingInternals::IndexType eq_id;
            bool acceptable = false;
            while(!acceptable)
            {
                auto randomid = static_cast<DistTestingInternals::IndexType>(dis(gen));
                if(static_cast<DistTestingInternals::IndexType>(randomid) > 0 &&
                        static_cast<DistTestingInternals::IndexType>(randomid) < ndof-1)
                {
                    acceptable=true;
                    eq_id = randomid * block_size;
                }
            }

            for(DistTestingInternals::IndexType k = 0; k<block_size; ++k){
                connectivities[Index][j*block_size+k] = eq_id+k;
            }
        }
    });

    std::cout << "Finishing generation - time = " << timer.ElapsedSeconds() << std::endl;

    return connectivities;
}

template<class TIndexType>
std::vector<TIndexType> ComputeBounds( TIndexType N,
                                       TIndexType Ndivisions,
                                       TIndexType current_rank
                                     )
{
    std::vector<int> partition;
    OpenMPUtils::DivideInPartitions(N,Ndivisions,partition);

    std::vector<TIndexType> bounds
    {
        static_cast<TIndexType>(partition[current_rank]),
        static_cast<TIndexType>(partition[current_rank+1])
    };
    return bounds;
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedGraphConstructionMPI, KratosMPICoreFastSuite)
{
    DataCommunicator& rComm=ParallelEnvironment::GetDefaultDataCommunicator();
    int world_size =rComm.Size();
    int my_rank = rComm.Rank();

    auto dofs_bounds = DistTestingInternals::ComputeBounds<DistTestingInternals::IndexType>(40, world_size, my_rank);
    auto reference_A_map = DistTestingInternals::GetReferenceMatrixAsMap(dofs_bounds);

    auto el_bounds = DistTestingInternals::ComputeBounds<DistTestingInternals::IndexType>(31, world_size, my_rank);
    const auto connectivities = DistTestingInternals::ElementConnectivities(el_bounds);

    DistributedSparseGraph<DistTestingInternals::IndexType> Agraph(dofs_bounds[1]-dofs_bounds[0], rComm);

    IndexPartition<DistTestingInternals::IndexType>(connectivities.size()).for_each([&](DistTestingInternals::IndexType i)
    {
        Agraph.AddEntries(connectivities[i]);
    });
    Agraph.Finalize();

    DistTestingInternals::CheckGraph(Agraph, reference_A_map);

}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedCSRConstructionMPI, KratosMPICoreFastSuite)
{
    DataCommunicator& rComm=ParallelEnvironment::GetDefaultDataCommunicator();
    int world_size =rComm.Size();
    int my_rank = rComm.Rank();

    auto dofs_bounds = DistTestingInternals::ComputeBounds<DistTestingInternals::IndexType>(40, world_size, my_rank);
    auto reference_A_map = DistTestingInternals::GetReferenceMatrixAsMap(dofs_bounds);

    auto el_bounds = DistTestingInternals::ComputeBounds<DistTestingInternals::IndexType>(31, world_size, my_rank);
    const auto connectivities = DistTestingInternals::ElementConnectivities(el_bounds);

    DistributedSparseGraph<DistTestingInternals::IndexType> Agraph(dofs_bounds[1]-dofs_bounds[0], rComm);

    IndexPartition<DistTestingInternals::IndexType>(connectivities.size()).for_each([&](DistTestingInternals::IndexType i)
    {
        Agraph.AddEntries(connectivities[i]);
    });
    Agraph.Finalize();

    //FEM assembly
    DistributedCsrMatrix<double,DistTestingInternals::IndexType> A(Agraph);
    A.BeginAssemble();
    for(const auto& c : connectivities)
    {
        Matrix data(c.size(),c.size(),1.0);
        A.Assemble(data,c);
    }
    A.FinalizeAssemble();

    for(const auto& item : reference_A_map)
    {
        auto GlobalI = item.first.first;
        auto GlobalJ = item.first.second;
        const double reference_value = item.second;
        const double value = A.GetLocalDataByGlobalId(GlobalI,GlobalJ);
        KRATOS_CHECK_NEAR(reference_value, value, 1e-14);

    }


}


KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(BenchmarkDistributedGraphConstructionMPI, KratosMPICoreFastSuite)
{
    DataCommunicator& rComm=ParallelEnvironment::GetDefaultDataCommunicator();
    int world_size =rComm.Size();
    int my_rank = rComm.Rank();

    const DistTestingInternals::IndexType block_size = 4;
    const DistTestingInternals::IndexType nodes_in_elem = 4;
    const DistTestingInternals::IndexType nel = 1e2; //set to 1e6 or 1e7 for a more realistic benchmark
    const DistTestingInternals::IndexType ndof = nel/6;
    const DistTestingInternals::IndexType standard_dev = 100; //reducing this implies using more optimized numbering

    auto el_bounds = DistTestingInternals::ComputeBounds<DistTestingInternals::IndexType>(nel, world_size, my_rank);
    auto dofs_bounds = DistTestingInternals::ComputeBounds<DistTestingInternals::IndexType>(ndof*block_size, world_size, my_rank);

    auto connectivities = DistTestingInternals::RandomElementConnectivities(block_size,
                          nodes_in_elem,
                          el_bounds[0],
                          el_bounds[1],
                          ndof,
                          standard_dev);
    rComm.Barrier(); //to ensure fair timings
    auto timer =  BuiltinTimer();
    DistributedSparseGraph<DistTestingInternals::IndexType> Agraph(dofs_bounds[1]-dofs_bounds[0], rComm);

    IndexPartition<DistTestingInternals::IndexType>(connectivities.size()).for_each([&](DistTestingInternals::IndexType i)
    {
        Agraph.AddEntries(connectivities[i]);
    });
    Agraph.Finalize();
    rComm.Barrier(); //to ensure fair timings
    std::cout << "graph - time = " << timer.ElapsedSeconds() << std::endl;

}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedSystemVectorConstructionMPI, KratosMPICoreFastSuite)
{
    DataCommunicator& rComm=ParallelEnvironment::GetDefaultDataCommunicator();
    int world_size =rComm.Size();
    int my_rank = rComm.Rank();

    auto dofs_bounds = DistTestingInternals::ComputeBounds<DistTestingInternals::IndexType>(40, world_size, my_rank);
    auto reference_A_map = DistTestingInternals::GetReferenceMatrixAsMap(dofs_bounds);
    auto reference_b_map = DistTestingInternals::GetReferencebAsMap(dofs_bounds);

    auto el_bounds = DistTestingInternals::ComputeBounds<DistTestingInternals::IndexType>(31, world_size, my_rank);
    const auto connectivities = DistTestingInternals::ElementConnectivities(el_bounds);

    DistributedSparseGraph<DistTestingInternals::IndexType> Agraph(dofs_bounds[1]-dofs_bounds[0], rComm);

    IndexPartition<DistTestingInternals::IndexType>(connectivities.size()).for_each([&](DistTestingInternals::IndexType i)
    {
        Agraph.AddEntries(connectivities[i]);
    });
    Agraph.Finalize();

    //FEM assembly
    DistributedSystemVector<> b(Agraph);

    b.SetValue(0.0);
    b.BeginAssemble();

    for(const auto& c : connectivities)
    {
        Vector vdata(c.size(),1.0);
        b.Assemble(vdata,c);
    }

    b.FinalizeAssemble();
    DistTestingInternals::IndexType local_size = b.LocalSize();
    for(unsigned int i=0; i<local_size; ++i)
    {
        auto global_i = b.GetNumbering().GlobalId(i);
        auto it = reference_b_map.find(global_i);
        const auto& ref_value = it->second;
        KRATOS_CHECK_NEAR(b(i),  ref_value, 1e-14 );
    }
    //test importing
    std::vector<DistTestingInternals::IndexType> to_import{39, 0,37,2};
    DistributedVectorImporter<double,DistTestingInternals::IndexType> importer(b.GetComm(),to_import,b.GetNumbering());   //this operation is expensive since it requires mounting the communication plan
    auto x = importer.ImportData(b);
    KRATOS_CHECK_NEAR(x[0],  3, 1e-14 );
    KRATOS_CHECK_NEAR(x[1],  1, 1e-14 );
    KRATOS_CHECK_NEAR(x[2],  4, 1e-14 );
    KRATOS_CHECK_NEAR(x[3],  2, 1e-14 );

    //Test SPMV
    DistributedCsrMatrix<double, DistTestingInternals::IndexType> A(Agraph);
    A.BeginAssemble();
    for(const auto& c : connectivities)
    {
        Matrix data(c.size(),c.size(),1.0);
        A.Assemble(data,c);
    }
    A.FinalizeAssemble();

    //here we test SPMV by a vector of 1s
    DistributedSystemVector<> y(Agraph);
    y.SetValue(0.0);
    b.SetValue(1.0);

    A.SpMV(b,y); //y+=A*b

    std::vector<double> reference_spmv_res{4,12,8,12,12,12,20,24,16,16,8,16,12,4,24,12,20,12,24,12,0,4,16,4,16,16,24,4,20,8,8,12,4,16,4,20,8,16,4,12};
    for(unsigned int i=0; i<y.LocalSize(); ++i)
    {
        DistTestingInternals::IndexType global_i = y.GetNumbering().GlobalId(i);
        KRATOS_CHECK_NEAR(y[i],  reference_spmv_res[global_i], 1e-14 );
    }

    //testing AMGCL interface
    bool move_to_backend = true;
    auto offdiag_global_index2 = A.GetOffDiagonalIndex2DataInGlobalNumbering();
    auto pAmgcl = AmgclDistributedCSRConversionUtilities::ConvertToAmgcl<double,IndexType>(A,offdiag_global_index2,move_to_backend);

    y.SetValue(0.0);
    b.SetValue(1.0);
    pAmgcl->mul(1.0,b,1.0,y);

    for(unsigned int i=0; i<y.LocalSize(); ++i)
    {
        IndexType global_i = y.GetNumbering().GlobalId(i);
        KRATOS_CHECK_NEAR(y[i],  reference_spmv_res[global_i], 1e-14 );
    }

    //check round trip krtos_csr->amgcl->kratos_csr
    move_to_backend = false; //note that the "move_to_backend needs to be set to false to be able to do the round trip"
    offdiag_global_index2 = A.GetOffDiagonalIndex2DataInGlobalNumbering();
    pAmgcl = AmgclDistributedCSRConversionUtilities::ConvertToAmgcl<double,IndexType>(A,offdiag_global_index2,move_to_backend);

    //convert back to CSR matrix
    DistributedCsrMatrix<double, DistTestingInternals::IndexType>::Pointer pAreconverted = 
        AmgclDistributedCSRConversionUtilities::ConvertToCsrMatrix<double,IndexType>(*pAmgcl, rComm);

    y.SetValue(0.0);
    b.SetValue(1.0);

    pAreconverted->SpMV(b,y); //y+=A*b
    for(unsigned int i=0; i<y.LocalSize(); ++i)
    {
        DistTestingInternals::IndexType global_i = y.GetNumbering().GlobalId(i);
        KRATOS_CHECK_NEAR(y[i],  reference_spmv_res[global_i], 1e-14 );
    }
}


KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(RectangularMatrixConstructionMPI, KratosMPICoreFastSuite)
{
    DataCommunicator& rComm=ParallelEnvironment::GetDefaultDataCommunicator();
    DistTestingInternals::IndexType col_divider = 3; //ratio of size between columns and row indices

    //*************************************************************************
    //compute reference solution - serial mode
    std::vector<DistTestingInternals::IndexType> all_el_bounds{0,31};
    const auto all_connectivities = DistTestingInternals::ElementConnectivities(all_el_bounds);
    SparseContiguousRowGraph<DistTestingInternals::IndexType> Agraph_serial(40);

    IndexPartition<DistTestingInternals::IndexType>(all_connectivities.size()).for_each([&](DistTestingInternals::IndexType i)
    {
        std::vector<DistTestingInternals::IndexType> row_ids = all_connectivities[i];
        std::vector<DistTestingInternals::IndexType> col_ids{row_ids[0]/col_divider, row_ids[1]/col_divider};

        Agraph_serial.AddEntries(row_ids, col_ids);
    });
    Agraph_serial.Finalize();

    CsrMatrix<double, DistTestingInternals::IndexType> Aserial(Agraph_serial);

    Aserial.BeginAssemble();
    for(const auto& c : all_connectivities)
    {
        std::vector<DistTestingInternals::IndexType> row_ids = c;
        std::vector<DistTestingInternals::IndexType> col_ids{row_ids[0]/col_divider, row_ids[1]/col_divider};
        Matrix data(row_ids.size(),col_ids.size(),1.0);
        Aserial.Assemble(data,row_ids, col_ids);
    }
    Aserial.FinalizeAssemble();

    auto reference_A_map = Aserial.ToMap();

    // //here we test SPMV by a vector of 1s
    SystemVector<> yserial(Aserial.size1()); //destination vector
    yserial.SetValue(0.0);

    SystemVector<> xserial(Aserial.size2()); //origin vector
    xserial.SetValue(1.0);

    Aserial.SpMV(xserial,yserial);

    //*************************************************************************
    int world_size =rComm.Size();
    int my_rank = rComm.Rank();

    auto dofs_bounds = DistTestingInternals::ComputeBounds<DistTestingInternals::IndexType>(40, world_size, my_rank);
    auto el_bounds = DistTestingInternals::ComputeBounds<DistTestingInternals::IndexType>(31, world_size, my_rank);
    const auto connectivities = DistTestingInternals::ElementConnectivities(el_bounds);

    DistributedSparseGraph<DistTestingInternals::IndexType> Agraph(dofs_bounds[1]-dofs_bounds[0], rComm);

    IndexPartition<DistTestingInternals::IndexType>(connectivities.size()).for_each([&](DistTestingInternals::IndexType i)
    {
        std::vector<DistTestingInternals::IndexType> row_ids = connectivities[i];
        std::vector<DistTestingInternals::IndexType> col_ids{row_ids[0]/col_divider, row_ids[1]/col_divider};
        Agraph.AddEntries(row_ids, col_ids);
    });
    Agraph.Finalize();

    DistributedCsrMatrix<double, DistTestingInternals::IndexType> A(Agraph);

    A.BeginAssemble();
    for(const auto& c : connectivities)
    {
        std::vector<DistTestingInternals::IndexType> row_ids = c;
        std::vector<DistTestingInternals::IndexType> col_ids{row_ids[0]/col_divider, row_ids[1]/col_divider};
        Matrix data(row_ids.size(),col_ids.size(),1.0);
        A.Assemble(data,row_ids, col_ids);
    }
    A.FinalizeAssemble();

    auto Amap = A.ToMap();
    //check that all "local" values in reference_A_map also appear in A_map
    for(const auto& item : reference_A_map)
    {
        DistTestingInternals::IndexType i = item.first.first;
        if(A.GetRowNumbering().IsLocal(i))
        {
            DistTestingInternals::IndexType j = item.first.second;
            double reference_v = item.second;
            if(Amap.find(item.first) == Amap.end())
                KRATOS_ERROR << "entry " << i << " " << j << "not found in A_map" <<std::endl;

            double v = Amap.find(item.first)->second;
            KRATOS_CHECK_NEAR(v,reference_v,1e-14);
        }
    }

    // //here we test SPMV by a vector of 1s
    DistributedSystemVector<> y(Agraph); //destination vector
    y.SetValue(0.0);

    DistributedSystemVector<> x(A.GetColNumbering()); //origin vector
    x.SetValue(1.0);
    A.SpMV(x,y);

    for(DistTestingInternals::IndexType i_local=0; i_local<y.LocalSize(); ++i_local)
    {
        DistTestingInternals::IndexType i_global = y.GetNumbering().GlobalId(i_local);
        KRATOS_CHECK_NEAR(y[i_local], yserial(i_global), 1e-14);
    }

    //check TransposeSpMV
    x.SetValue(0.0);
    y.SetValue(1.0);
    A.TransposeSpMV(y,x);

    std::vector<double> reference_transpose_spmv_res{20,8,16,20,8,32,28,4,52,16,8,24,12};
    for(unsigned int i=0; i<x.LocalSize(); ++i)
    {
        DistTestingInternals::IndexType global_i = x.GetNumbering().GlobalId(i);
        KRATOS_CHECK_NEAR(x[i],  reference_transpose_spmv_res[global_i], 1e-14 );
    }

}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedSystemVectorOperationsMPI, KratosMPICoreFastSuite)
{
    DataCommunicator& rComm=ParallelEnvironment::GetDefaultDataCommunicator();

    DistTestingInternals::IndexType local_size = 4;
    DistributedNumbering<DistTestingInternals::IndexType> numbering(rComm,local_size);
    DistributedSystemVector<double,DistTestingInternals::IndexType> a(numbering);
    KRATOS_CHECK_EQUAL(a.LocalSize(), local_size);
    a.SetValue(5.0);
    for(unsigned int i=0; i<a.LocalSize(); ++i)
        KRATOS_CHECK_NEAR(a[i], 5.0,1e-14);

    DistributedSystemVector<double,DistTestingInternals::IndexType> b(numbering);
    b.SetValue(3.0);
    KRATOS_CHECK_EQUAL(b.LocalSize(), local_size);
    for(unsigned int i=0; i<b.LocalSize(); ++i)
        KRATOS_CHECK_NEAR(b[i], 3.0,1e-14);

    DistributedSystemVector<double,DistTestingInternals::IndexType> c(a);
    KRATOS_CHECK_EQUAL(c.LocalSize(), local_size);
    KRATOS_CHECK_EQUAL(c.Size(), local_size*rComm.Size());
    for(unsigned int i=0; i<c.LocalSize(); ++i)
        KRATOS_CHECK_NEAR(c[i], 5.0,1e-14);

    c += b;
    for(unsigned int i=0; i<c.LocalSize(); ++i)
        KRATOS_CHECK_NEAR(c[i], 8.0,1e-14);

    c -= b;
    for(unsigned int i=0; i<c.LocalSize(); ++i)
        KRATOS_CHECK_NEAR(c[i], 5.0,1e-14);

    c.Add(3.0,a);
    for(unsigned int i=0; i<c.LocalSize(); ++i)
        KRATOS_CHECK_NEAR(c[i], 20.0,1e-14);

    c*=2.0;
    for(unsigned int i=0; i<c.LocalSize(); ++i)
        KRATOS_CHECK_NEAR(c[i], 40.0,1e-14);

    c/=4.0;
    for(unsigned int i=0; i<c.LocalSize(); ++i)
        KRATOS_CHECK_NEAR(c[i], 10.0,1e-14);
}


KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(Small1dLaplacianAmgclConstruction, KratosMPICoreFastSuite)
{
    typedef std::size_t IndexType;
    DataCommunicator& rComm=ParallelEnvironment::GetDefaultDataCommunicator();
    int world_size =rComm.Size();
    int my_rank = rComm.Rank();

    IndexType max_size = 3;
    auto dofs_bounds = DistTestingInternals::ComputeBounds<IndexType>(max_size, world_size, my_rank);

    Matrix local_matrix(2,2); //we will assemble a 1D laplacian
    local_matrix(0,0) = 1.0;
    local_matrix(0,1) = -1.0;
    local_matrix(1,0) = -1.0;
    local_matrix(1,1) = 1.0;
    DenseVector<IndexType> connectivities(2);

    DistributedSparseGraph<IndexType> Agraph(dofs_bounds[1]-dofs_bounds[0], rComm);

    for(IndexType i=dofs_bounds[0]; i<dofs_bounds[1]; ++i)
    {
        if(i+1<max_size)
        {
            connectivities[0] = i;
            connectivities[1] = i+1;
            Agraph.AddEntries(connectivities);
        }
    }

    Agraph.Finalize();

    //Test SPMV
    DistributedCsrMatrix<double, IndexType> A(Agraph);

    A.BeginAssemble();

    for(IndexType i=dofs_bounds[0]; i<dofs_bounds[1]; ++i)
    {
        if(i+1<max_size)
        {
            connectivities[0] = i;
            connectivities[1] = i+1;
            A.Assemble(local_matrix,connectivities);
        }
    }
    A.FinalizeAssemble();

    auto offdiag_global_index2 = A.GetOffDiagonalIndex2DataInGlobalNumbering();

    bool move_to_backend=false;
    auto pAmgcl = AmgclDistributedCSRConversionUtilities::ConvertToAmgcl<double,IndexType>(A,offdiag_global_index2,move_to_backend);

    auto Areconverted 
        = AmgclDistributedCSRConversionUtilities::ConvertToCsrMatrix<double,IndexType>(*pAmgcl, A.GetComm());

}





KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(SmallRectangularDistributedMatrixMatrixMultiply, KratosMPICoreFastSuite)
{
    typedef std::size_t IndexType;
    DataCommunicator& rComm=ParallelEnvironment::GetDefaultDataCommunicator();
    int world_size =rComm.Size();
    int my_rank = rComm.Rank();

    // matrix A
    //[[1,0,0,7,2],
    // [0,3,0,0,0],
    // [0,0,0,7,7]]
    DistTestingInternals::MatrixMapType Amap
    {
        {{0,0},1.0}, {{0,3},7.0}, {{0,4},2.0},
        {{1,1},3.0},
        {{2,3},7.0}, {{2,4},7.0}
    };

    IndexType Asize1 = 3;
    auto dofs_bounds = DistTestingInternals::ComputeBounds<IndexType>(Asize1, world_size, my_rank);

    DistributedSparseGraph<IndexType> Agraph(dofs_bounds[1]-dofs_bounds[0], rComm);
    for(const auto& item : Amap)
    {
        IndexType I = item.first.first;
        IndexType J = item.first.second;
        if( Agraph.GetRowNumbering().IsLocal(I))
            Agraph.AddEntry(I,J);
    }
    Agraph.Finalize();

    DistributedCsrMatrix<double, IndexType> A(Agraph);
    A.BeginAssemble();
    for(const auto& item : Amap)
    {
        IndexType I = item.first.first;
        IndexType J = item.first.second;
        double value = item.second;
        if( A.GetRowNumbering().IsLocal(I))
            A.AssembleEntry(value,I,J);
    }
    A.FinalizeAssemble();

    // matrix B
    // [[1,0,0],
    // [0,2,3],
    // [0,3,0],
    // [0,0,0],
    // [5,0,6]]
    DistTestingInternals::MatrixMapType Bmap
    {
        {{0,0},1.0},
        {{1,1},2.0}, {{1,2},3.0},
        {{2,1},3.0},
        //empty row
        {{4,0},5.0}, {{4,2},6.0},
    };

    IndexType Bsize1 = 5;
    dofs_bounds = DistTestingInternals::ComputeBounds<IndexType>(Bsize1, world_size, my_rank);

    DistributedSparseGraph<IndexType> Bgraph(dofs_bounds[1]-dofs_bounds[0], rComm);
    for(const auto& item : Bmap)
    {
        IndexType I = item.first.first;
        IndexType J = item.first.second;
        if( Bgraph.GetRowNumbering().IsLocal(I))
            Bgraph.AddEntry(I,J);
    }
    Bgraph.Finalize();

    DistributedCsrMatrix<double, IndexType> B(Bgraph);
    B.BeginAssemble();
    for(const auto& item : Bmap)
    {
        IndexType I = item.first.first;
        IndexType J = item.first.second;
        double value = item.second;
        if( B.GetRowNumbering().IsLocal(I))
            B.AssembleEntry(value,I,J);
    }
    B.FinalizeAssemble();

    // //Cref = A@B
    //[[11,  0, 12],
    // [ 0,  6,  9],
    // [35,  0, 42]]
    DistTestingInternals::MatrixMapType Cref
    {
        {{0,0},11.0}, {{0,2},12.0},
        {{1,1},6.0}, {{1,2},9.0},
        {{2,0},35.0}, {{2,2},42.0}
    };

    auto pC = DistributedAmgclCSRSpMMUtilities::SparseMultiply<double,IndexType>(A,B); //C=A*B
    auto& C = *pC;

    for(const auto& item : Cref)
    {
        IndexType I = item.first.first;
        if( C.GetRowNumbering().IsLocal(I))
        {
            IndexType J = item.first.second;
            double ref_value = item.second;
            std::cout << I << " " << J << " " << ref_value << std::endl;
            KRATOS_CHECK_EQUAL(ref_value,C.GetLocalDataByGlobalId(I,J));
        }
    }
}




} // namespace Testing
} // namespace Kratos
