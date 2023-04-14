//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#pragma once

// System includes


// External includes
#include "metis.h"


// Project includes
#include "includes/io.h"


namespace Kratos {

///@addtogroup MetisApplication
///@{

///@name Kratos Classes
///@{

/// This class contains legacy versions of utilities used by the metis partitioners.
/** The new impleementation of these functionalities is part of the MetisPartitioningUtilities.
*/
class KRATOS_API(METIS_APPLICATION) LegacyPartitioningUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LegacyPartitioningUtilities
    KRATOS_CLASS_POINTER_DEFINITION(LegacyPartitioningUtilities);

    using idxtype = idx_t; // from metis
    using PartitionIndicesType = std::vector<idxtype>;
    using SizeType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LegacyPartitioningUtilities() = delete;

    /// Copy constructor.
    LegacyPartitioningUtilities(LegacyPartitioningUtilities const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    LegacyPartitioningUtilities& operator=(LegacyPartitioningUtilities const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    static void CalculateDomainsGraph(
        IO::GraphType& rDomainsGraph,
        SizeType NumberOfElements,
        IO::ConnectivitiesContainerType& ElementsConnectivities,
        PartitionIndicesType const& NPart,
        PartitionIndicesType const&  EPart);

    static void DividingNodes(
        IO::PartitionIndicesContainerType& rNodesAllPartitions,
        IO::ConnectivitiesContainerType& ElementsConnectivities,
        IO::ConnectivitiesContainerType& ConditionsConnectivities,
        PartitionIndicesType const& NodesPartitions,
        PartitionIndicesType const& ElementsPartitions,
        PartitionIndicesType const& ConditionsPartitions);

    static void DividingElements(
        IO::PartitionIndicesContainerType& rElementsAllPartitions,
        PartitionIndicesType const& ElementsPartitions);

    static void DividingConditions(
        IO::PartitionIndicesContainerType& rConditionsAllPartitions,
        PartitionIndicesType const& ConditionsPartitions);

    static void ConvertKratosToCSRFormat(
        IO::ConnectivitiesContainerType& KratosFormatNodeConnectivities,
        idxtype** NodeIndices,
        idxtype** NodeConnectivities);

    ///@}

}; // Class LegacyPartitioningUtilities

///@}

///@} addtogroup block

} // namespace Kratos
