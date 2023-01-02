//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "nearest_neighbor_mapper.h"
#include "mapping_application_variables.h"

namespace Kratos {

void NearestNeighborInterfaceInfo::ProcessSearchResult(const InterfaceObject& rInterfaceObject)
{
    SetLocalSearchWasSuccessful();

    const double neighbor_distance = MapperUtilities::ComputeDistance(this->Coordinates(), rInterfaceObject.Coordinates());

    if (neighbor_distance < mNearestNeighborDistance) {
        mNearestNeighborDistance = neighbor_distance;
        mNearestNeighborId = rInterfaceObject.pGetBaseNode()->GetValue(INTERFACE_EQUATION_ID);
    }
}

void NearestNeighborLocalSystem::CalculateAll(MatrixType& rLocalMappingMatrix,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds,
                    MapperLocalSystem::PairingStatus& rPairingStatus) const
{
    if (mInterfaceInfos.size() > 0) {
        rPairingStatus = MapperLocalSystem::PairingStatus::InterfaceInfoFound;

        if (rOriginIds.size()      != 1) rOriginIds.resize(1);
        if (rDestinationIds.size() != 1) rDestinationIds.resize(1);

        int nearest_neighbor_id;
        double nearest_neighbor_distance;
        mInterfaceInfos[0]->GetValue(nearest_neighbor_id, MapperInterfaceInfo::InfoType::Dummy);
        mInterfaceInfos[0]->GetValue(nearest_neighbor_distance, MapperInterfaceInfo::InfoType::Dummy);

        for (std::size_t i=1; i<mInterfaceInfos.size(); ++i) {
            // no check if this InterfaceInfo is an approximation is necessary
            // bcs this does not exist for NearestNeighbor
            double distance;
            mInterfaceInfos[i]->GetValue(distance, MapperInterfaceInfo::InfoType::Dummy);

            if (distance < nearest_neighbor_distance) {
                nearest_neighbor_distance = distance;
                mInterfaceInfos[i]->GetValue(nearest_neighbor_id, MapperInterfaceInfo::InfoType::Dummy);
                rOriginIds.resize(1);
                rOriginIds[0] = nearest_neighbor_id;
            } else if (distance == nearest_neighbor_distance) {
                mInterfaceInfos[i]->GetValue(nearest_neighbor_id, MapperInterfaceInfo::InfoType::Dummy);
                rOriginIds.push_back(nearest_neighbor_id);
            }
        }

        if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != rOriginIds.size()) {
            rLocalMappingMatrix.resize(1, rOriginIds.size(), false);
        }
        const double mapping_weight = 1.0/rOriginIds.size();
        for (IndexType i=0; i<rOriginIds.size(); ++i) {
            rLocalMappingMatrix(0,i) = mapping_weight;
        }

        KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;
        rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);
    }
    else ResizeToZero(rLocalMappingMatrix, rOriginIds, rDestinationIds, rPairingStatus);
}

void NearestNeighborLocalSystem::PairingInfo(std::ostream& rOStream, const int EchoLevel) const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

    rOStream << "NearestNeighborLocalSystem based on " << mpNode->Info();
    if (EchoLevel > 3) {
        rOStream << " at Coordinates " << Coordinates()[0] << " | " << Coordinates()[1] << " | " << Coordinates()[2];
    }
}

void NearestNeighborLocalSystem::SetPairingStatusForPrinting()
{
    if (mPairingStatus == MapperLocalSystem::PairingStatus::Approximation) {
        mpNode->SetValue(PAIRING_STATUS, 0);
    } else {
        mpNode->SetValue(PAIRING_STATUS, -1);
    }
}

}  // namespace Kratos.
