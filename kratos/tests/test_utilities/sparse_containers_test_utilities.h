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
#include <utility>
#include <iostream>
#include <random>

// External includes

// Project includes
#include "containers/csr_matrix.h"
#include "containers/sparse_graph.h"
#include "containers/sparse_contiguous_row_graph.h"
#include "containers/system_vector.h"
#include "includes/key_hash.h"

namespace Kratos
{

class SparseContainersTestUtilities
{
public:

    using IndexType = std::size_t;

    using ElementConnectivityType = std::vector<std::vector<IndexType>>;

    using MatrixMapType = std::unordered_map<std::pair<IndexType, IndexType>, double, PairHasher<IndexType, IndexType>, PairComparor<IndexType, IndexType>>;

    static ElementConnectivityType RandomElementConnectivities(
        const IndexType BlockSize,
        const IndexType NodesInElem,
        const IndexType IndexBegin,
        const IndexType IndexEnd,
        const IndexType NDof,
        const IndexType StandardDev);

    static ElementConnectivityType ElementConnectivities();

    static MatrixMapType GetReferenceMatrixAsMap();

    template< class TGraphType>
    static bool CheckGraph(
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

    static bool CheckCSRGraphArrays(
        const Kratos::span<IndexType>& rRowIndices,
        const Kratos::span<IndexType>& rColIndices,
        const MatrixMapType& rReferenceGraph);

    static bool CheckCSRMatrix(
        const CsrMatrix<>& rA,
        const MatrixMapType& rReferenceGraph);

    template<typename TGraphType>
    static std::unique_ptr<TGraphType> AssembleGraph(
        ElementConnectivityType& rConnectivities,
        IndexType GraphSize)
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

};

} // namespace Kratos
