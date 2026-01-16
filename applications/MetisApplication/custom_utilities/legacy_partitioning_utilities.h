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
//                   Vicente Mataix Ferrandiz
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

///@}
///@name Kratos Classes
///@{

/**
 * @class LegacyPartitioningUtilities
 * @ingroup MetisApplication
 * @brief A utility class for legacy partitioning operations in the MetisApplication.
 * @details This class provides static methods for performing various partitioning-related operations,
 * such as calculating domain graphs, dividing nodes, elements, and conditions among partitions,
 * and converting Kratos format connectivities to CSR format.
 * @note The new implementation of these functionalities is part of the MetisPartitioningUtilities.
 * @author Philipp Bucher
 */
class KRATOS_API(METIS_APPLICATION) LegacyPartitioningUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LegacyPartitioningUtilities
    KRATOS_CLASS_POINTER_DEFINITION(LegacyPartitioningUtilities);

    /// Type alias for Metis index type.
    using idxtype = idx_t;

    /// Type alias for partition indices.
    using PartitionIndicesType = std::vector<idxtype>;

    /// Type alias for index type.
    using IndexType = std::size_t;

    /// Type alias for size type.
    using SizeType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor is deleted to prevent instantiation.
     */
    LegacyPartitioningUtilities() = delete;

    /**
     * @brief Copy constructor is deleted to prevent copying.
     */
    LegacyPartitioningUtilities(LegacyPartitioningUtilities const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Assignment operator is deleted to prevent assignment.
     */
    LegacyPartitioningUtilities& operator=(LegacyPartitioningUtilities const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Calculates the domains graph for partitioning.
     * @param rDomainsGraph The resulting domains graph.
     * @param NumberOfElements The number of elements in the graph.
     * @param ElementsConnectivities The connectivities of the elements.
     * @param NPart The node partitions.
     * @param EPart The element partitions.
     */
    static void CalculateDomainsGraph(
        IO::GraphType& rDomainsGraph,
        SizeType NumberOfElements,
        IO::ConnectivitiesContainerType& ElementsConnectivities,
        PartitionIndicesType const& NPart,
        PartitionIndicesType const& EPart
        );

    /**
     * @brief Divides nodes among partitions.
     * @details Legacy version not including geometries and constraints.
     * @param rNodesAllPartitions The resulting nodes for all partitions.
     * @param rElementsConnectivities The connectivities of the elements.
     * @param rConditionsConnectivities The connectivities of the conditions.
     * @param rNodesPartitions The node partitions.
     * @param rElementsPartitions The element partitions.
     * @param rConditionsPartitions The condition partitions.
     */
    static void DividingNodes(
        IO::PartitionIndicesContainerType& rNodesAllPartitions,
        IO::ConnectivitiesContainerType& rElementsConnectivities,
        IO::ConnectivitiesContainerType& rConditionsConnectivities,
        const PartitionIndicesType& rNodesPartitions,
        const PartitionIndicesType& rElementsPartitions,
        const PartitionIndicesType& rConditionsPartitions
        );

    /**
     * @brief Divides nodes among partitions including geometries and constraints.
     * @param rNodesAllPartitions The resulting nodes for all partitions.
     * @param rGeometriesAllPartitions The resulting geometries for all partitions.
     * @param rElementsConnectivities The connectivities of the elements.
     * @param rConditionsConnectivities The connectivities of the conditions.
     * @param rConstraintsAllPartitions The resulting constraints for all partitions.
     * @param rNodesPartitions The node partitions.
     * @param rGeometriesPartitions The geometry partitions.
     * @param rElementsPartitions The element partitions.
     * @param rConditionsPartitions The condition partitions.
     * @param rConstraintsPartitions The constraints partitions.
     */
    static void DividingNodes(
        IO::PartitionIndicesContainerType& rNodesAllPartitions,
        IO::PartitionIndicesContainerType& rGeometriesAllPartitions,
        IO::ConnectivitiesContainerType& rElementsConnectivities,
        IO::ConnectivitiesContainerType& rConditionsConnectivities,
        IO::PartitionIndicesContainerType& rConstraintsAllPartitions,
        const PartitionIndicesType& rNodesPartitions,
        const PartitionIndicesType& rGeometriesPartitions,
        const PartitionIndicesType& rElementsPartitions,
        const PartitionIndicesType& rConditionsPartitions,
        const PartitionIndicesType& rConstraintsPartitions
        );

    /**
     * @brief Divides elements among partitions.
     * @param rGeometriesAllPartitions The resulting geometries for all partitions.
     * @param rGeometriesPartitions The geometry partitions.
     */
    static void DividingGeometries(
        IO::PartitionIndicesContainerType& rGeometriesAllPartitions,
        const PartitionIndicesType& rGeometriesPartitions
        );

    /**
     * @brief Divides elements among partitions.
     * @param rElementsAllPartitions The resulting elements for all partitions.
     * @param ElementsPartitions The element partitions.
     */
    static void DividingElements(
        IO::PartitionIndicesContainerType& rElementsAllPartitions,
        PartitionIndicesType const& ElementsPartitions
        );

    /**
     * @brief Divides conditions among partitions.
     * @param rConditionsAllPartitions The resulting conditions for all partitions.
     * @param ConditionsPartitions The condition partitions.
     */
    static void DividingConditions(
        IO::PartitionIndicesContainerType& rConditionsAllPartitions,
        PartitionIndicesType const& ConditionsPartitions
        );

    /**
     * @brief Divides constraints among partitions.
     * @param rConstraintsAllPartitions The resulting constraints for all partitions.
     * @param rConstraintsPartitions The constraint partitions.
     */
    static void DividingConstraints(
        IO::PartitionIndicesContainerType& rConstraintsAllPartitions,
        const PartitionIndicesType& rConstraintsPartitions
        );

    /**
     * @brief Converts Kratos format connectivities to CSR format.
     * @param KratosFormatNodeConnectivities The Kratos format node connectivities.
     * @param NodeIndices The resulting CSR format node indices.
     * @param NodeConnectivities The resulting CSR format node connectivities.
     */
    static void ConvertKratosToCSRFormat(
        IO::ConnectivitiesContainerType& KratosFormatNodeConnectivities,
        idxtype** NodeIndices,
        idxtype** NodeConnectivities
        );

    ///@}
};

///@}

///@} addtogroup block

} // namespace Kratos
