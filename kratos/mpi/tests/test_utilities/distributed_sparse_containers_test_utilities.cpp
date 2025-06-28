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
#include "distributed_sparse_containers_test_utilities.h"

namespace Kratos::Testing
{

DistributedSparseContainersTestUtilities::ElementConnectivityType DistributedSparseContainersTestUtilities::ElementConnectivities(const std::vector<IndexType> &rBounds)
{
    KRATOS_TRY

    const auto all_connectivities = SparseContainersTestUtilities::ElementConnectivities();

    ElementConnectivityType connectivities;
    for (IndexType i = rBounds[0]; i < rBounds[1]; ++i) {
        connectivities.push_back(all_connectivities[i]);
    }

    return connectivities;

    KRATOS_CATCH("")
}

DistributedSparseContainersTestUtilities::MatrixMapType DistributedSparseContainersTestUtilities::GetReferenceMatrixAsMap(const std::vector<IndexType> &rBounds)
{
    KRATOS_TRY

    const auto all_connectivities = SparseContainersTestUtilities::GetReferenceMatrixAsMap();

    KRATOS_ERROR_IF(rBounds[1] > all_connectivities.size())
        << "Bounds : " << rBounds << " exceed the total size: " << all_connectivities.size() << std::endl;

    MatrixMapType connectivities;
    for (auto& item : all_connectivities) {
        if (item.first.first >= rBounds[0] && item.first.first < rBounds[1]) {
            connectivities.insert(item);
        }
    }

    return connectivities;

    KRATOS_CATCH("")
}

std::map<IndexType, double> DistributedSparseContainersTestUtilities::GetReferenceVectorAsMap(const std::vector<IndexType> &rBounds)
{
    KRATOS_TRY

    const auto all_connectivities = SparseContainersTestUtilities::GetReferenceMatrixAsMap();

    KRATOS_ERROR_IF(rBounds[1] > all_connectivities.size())
        << "Bounds : " << rBounds << " exceed the total size: " << all_connectivities.size() << std::endl;

    std::vector<double> reference_vector{1,3,2,3,3,3,5,6,4,4,2,4,3,1,6,3,5,3,6,3,0,1,4,1,4,4,6,1,5,2,2,3,1,4,1,5,2,4,1,3};

    std::map<IndexType, double> reference_vector_map;
    for (IndexType i = 0; i < reference_vector.size(); ++i) {
        if (i >= rBounds[0] && i < rBounds[1]) {
            reference_vector_map.insert({i, reference_vector[i]});
        }
    }

    return reference_vector_map;

    KRATOS_CATCH("")
}

} // namespace Kratos::Testing
