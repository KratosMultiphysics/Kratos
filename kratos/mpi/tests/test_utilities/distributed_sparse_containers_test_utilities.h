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
//                   Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "utilities/openmp_utils.h"
#include "tests/test_utilities/sparse_containers_test_utilities.h"

namespace Kratos::Testing
{

class KRATOS_API(KRATOS_MPI_CORE) DistributedSparseContainersTestUtilities
{
public:

    using IndexType = typename SparseContainersTestUtilities::IndexType;

    using ElementConnectivityType = typename SparseContainersTestUtilities::ElementConnectivityType;

    using MatrixMapType = typename SparseContainersTestUtilities::MatrixMapType;

    static ElementConnectivityType ElementConnectivities(const std::vector<IndexType>& rBounds);

    static MatrixMapType GetReferenceMatrixAsMap(const std::vector<IndexType>& rBounds);

    static std::map<IndexType, double> GetReferenceVectorAsMap(const std::vector<IndexType>& rBounds);

    template <class TIndexType>
    static std::vector<TIndexType> ComputeBounds(
        TIndexType N,
        TIndexType Ndivisions,
        TIndexType CurrentRank)
    {
        std::vector<int> partition;
        OpenMPUtils::DivideInPartitions(N, Ndivisions, partition);
        std::vector<TIndexType> bounds{
            static_cast<TIndexType>(partition[CurrentRank]),
            static_cast<TIndexType>(partition[CurrentRank + 1])};

        return bounds;
    }

    template< class TSparseGraphType>
    static bool CheckGraph(
        const TSparseGraphType& rAgraph,
        const MatrixMapType& rReferenceGraph)
    {
        //check that all entries in Agraph are also in reference_A_map
        for(IndexType local_i = 0; local_i < rAgraph.LocalSize(); ++local_i) //i is the LOCAL index
        {
            auto i = rAgraph.GetRowNumbering().GlobalId(local_i);
            for(auto j : rAgraph[local_i] ) {
                KRATOS_ERROR_IF(rReferenceGraph.find({i, j}) == rReferenceGraph.end()) << "Entry " << i << "," << j << "not present in A graph" << std::endl;
            }
        }

        //check that all the entries of reference_A_map are also in Agraph
        for(auto item : rReferenceGraph) {
            auto global_i = item.first.first;
            auto global_j = item.first.second;
            KRATOS_ERROR_IF_NOT(rAgraph.Has(global_i, global_j)) << "Entry " << global_i << "," << global_j << " is in the reference graph but not in Agraph" << std::endl;
        }

        return true;
    }

};

} // namespace Kratos::Testing
