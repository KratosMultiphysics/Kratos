//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <vector>
#include <cstddef>

// External includes
#include "kaHIP_interface.h"

// Project includes
#include "includes/define.h"
#include "includes/io.h"

namespace Kratos
{

///@addtogroup KaHIPApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class KaHIPCSRConverter
 * @ingroup KaHIPApplication
 * @brief Utilities for converting Kratos mesh connectivity to the CSR graph format
 *        expected by the KaHIP partitioning API.
 * @details KaHIP (like METIS) requires the adjacency graph of the mesh in
 *          Compressed Sparse Row (CSR) format:
 *
 *          - `n`       : number of graph vertices
 *          - `xadj`    : size n+1; xadj[i] is the start of vertex i's neighbor list in adjncy
 *          - `adjncy`  : the concatenated neighbor lists, 0-indexed
 *          - `vwgt`    : optional vertex weights (uniform if nullptr)
 *          - `adjcwgt` : optional edge weights (uniform if nullptr)
 *
 *          Each undirected edge (u, v) must appear in both u's list and v's list with
 *          equal weight. This class enforces that invariant.
 *
 *          Node IDs in Kratos are arbitrary positive integers starting from 1; CSR
 *          indices are 0-based. The conversion subtracts 1 from each ID.
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KAHIP_APPLICATION) KaHIPCSRConverter
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KaHIPCSRConverter
    KRATOS_CLASS_POINTER_DEFINITION(KaHIPCSRConverter);

    /// Alias for Kratos IO size type
    using SizeType = IO::SizeType;

    /// Alias for Kratos connectivity container (vector of vectors of node IDs, 1-indexed)
    using ConnectivitiesContainerType = IO::ConnectivitiesContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Not instantiable — all methods are static utilities.
    KaHIPCSRConverter() = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Convert Kratos-format nodal connectivity to KaHIP CSR format.
     * @details The input is the adjacency list returned by @c IO::ReadNodalGraph, where
     *          each entry @c KratosFormatNodeConnectivities[i] lists the 1-indexed
     *          neighbour IDs of node i+1. This function converts that representation to
     *          the 0-indexed CSR arrays expected by @c kaffpa().
     *
     *          Memory layout:
     *          - @p rXAdj   has size n+1; @p rXAdj[i] is the start of vertex i's list
     *          - @p rAdjncy has size @p rXAdj[n] = 2*m; each undirected edge appears twice
     *
     * @param rKratosFormatConnectivities  [in]  Nodal adjacency in Kratos 1-indexed format
     * @param rXAdj                        [out] CSR row offsets (size n+1)
     * @param rAdjncy                      [out] CSR column indices (0-indexed)
     */
    static void ConvertToCSRFormat(
        const ConnectivitiesContainerType& rKratosFormatConnectivities,
        std::vector<kahip_idx>& rXAdj,
        std::vector<kahip_idx>& rAdjncy);

    /**
     * @brief Convert Kratos-format nodal connectivity to KaHIP CSR format with uniform weights.
     * @details Same as ConvertToCSRFormat() but also fills uniform vertex weights (all 1)
     *          and uniform edge weights (all 1) so the caller does not need to allocate them.
     *
     * @param rKratosFormatConnectivities  [in]  Nodal adjacency in Kratos 1-indexed format
     * @param rXAdj                        [out] CSR row offsets (size n+1)
     * @param rAdjncy                      [out] CSR column indices (0-indexed)
     * @param rVWgt                        [out] Vertex weights (all 1, size n)
     * @param rAdjcWgt                     [out] Edge weights (all 1, size 2*m)
     */
    static void ConvertToCSRFormatWithWeights(
        const ConnectivitiesContainerType& rKratosFormatConnectivities,
        std::vector<kahip_idx>& rXAdj,
        std::vector<kahip_idx>& rAdjncy,
        std::vector<int>& rVWgt,
        std::vector<kahip_idx>& rAdjcWgt);

    /**
     * @brief Validate a CSR graph for KaHIP correctness.
     * @details Checks:
     *          1. xadj is monotone non-decreasing and xadj[0] == 0
     *          2. adjncy entries are in [0, n)
     *          3. The graph is symmetric (if (u,v) is present then (v,u) is present)
     *          4. No self-loops
     * @param n        Number of vertices
     * @param rXAdj    CSR row offsets
     * @param rAdjncy  CSR column indices
     * @return @c true if the graph is valid, @c false otherwise (with KRATOS_WARNING)
     */
    static bool ValidateCSRGraph(
        int n,
        const std::vector<kahip_idx>& rXAdj,
        const std::vector<kahip_idx>& rAdjncy);

    ///@}

}; // class KaHIPCSRConverter

///@}

} // namespace Kratos
