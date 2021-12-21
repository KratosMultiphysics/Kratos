//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes

// External includes

// Project includes
#include "barycentric_mapper.h"
#include "mapping_application_variables.h"

namespace Kratos {

BarycentricInterfaceInfo::BarycentricInterfaceInfo(const std::size_t NumInterpolationNodes)
{
    Initialize(NumInterpolationNodes);
}

BarycentricInterfaceInfo::BarycentricInterfaceInfo(const CoordinatesArrayType& rCoordinates,
                                    const IndexType SourceLocalSystemIndex,
                                    const IndexType SourceRank,
                                    const std::size_t NumInterpolationNodes)
    : MapperInterfaceInfo(rCoordinates, SourceLocalSystemIndex, SourceRank)
{
    Initialize(NumInterpolationNodes);
}

MapperInterfaceInfo::Pointer BarycentricInterfaceInfo::Create() const
{
    return Kratos::make_shared<BarycentricInterfaceInfo>(mNodeIds.size());
}

MapperInterfaceInfo::Pointer BarycentricInterfaceInfo::Create(const CoordinatesArrayType& rCoordinates,
                                    const IndexType SourceLocalSystemIndex,
                                    const IndexType SourceRank) const
{
    return Kratos::make_shared<BarycentricInterfaceInfo>(
        rCoordinates,
        SourceLocalSystemIndex,
        SourceRank,
        mNodeIds.size());
}

InterfaceObject::ConstructionType BarycentricInterfaceInfo::GetInterfaceObjectType() const
{
    return InterfaceObject::ConstructionType::Node_Coords;
}

void BarycentricInterfaceInfo::GetValue(std::vector<int>& rValue,
                const InfoType ValueType) const
{
    rValue = mNodeIds;
}

void BarycentricInterfaceInfo::GetValue(std::vector<double>& rValue,
                const InfoType ValueType) const
{
    rValue = mNeighborCoordinates;
}

void BarycentricInterfaceInfo::ProcessSearchResult(const InterfaceObject& rInterfaceObject)
{
    SetLocalSearchWasSuccessful();

    const double neighbor_distance = MapperUtilities::ComputeDistance(this->Coordinates(), rInterfaceObject.Coordinates());

    // if (neighbor_distance < mNearestNeighborDistance) {
    //     mNearestNeighborDistance = neighbor_distance;
    //     mNearestNeighborId = rInterfaceObject.pGetBaseNode()->GetValue(INTERFACE_EQUATION_ID);
    // }
}

void BarycentricInterfaceInfo::Initialize(const std::size_t NumInterpolationNodes)
{
    mNodeIds.resize(NumInterpolationNodes);
    mNeighborCoordinates.resize(3*NumInterpolationNodes);
    std::fill(mNodeIds.begin(), mNodeIds.end(), -1);
    std::fill(mNeighborCoordinates.begin(), mNeighborCoordinates.end(), std::numeric_limits<double>::max());
}

void BarycentricInterfaceInfo::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MapperInterfaceInfo );
    rSerializer.save("NodeIds", mNodeIds);
    rSerializer.save("NeighborCoords", mNeighborCoordinates);
}

void BarycentricInterfaceInfo::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MapperInterfaceInfo );
    rSerializer.load("NodeIds", mNodeIds);
    rSerializer.load("NeighborCoords", mNeighborCoordinates);
}

void BarycentricLocalSystem::CalculateAll(MatrixType& rLocalMappingMatrix,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds,
                    MapperLocalSystem::PairingStatus& rPairingStatus) const
{
    if (mInterfaceInfos.size() > 0) {
        rPairingStatus = MapperLocalSystem::PairingStatus::InterfaceInfoFound;

        if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != 1) {
            rLocalMappingMatrix.resize(1, 1, false);
        }
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
            }
        }

        rLocalMappingMatrix(0,0) = 1.0;
        rOriginIds[0] = nearest_neighbor_id;
        KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;
        rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);
    }
    else ResizeToZero(rLocalMappingMatrix, rOriginIds, rDestinationIds, rPairingStatus);
}

void BarycentricLocalSystem::PairingInfo(std::ostream& rOStream, const int EchoLevel) const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

    rOStream << "BarycentricLocalSystem based on " << mpNode->Info();
    if (EchoLevel > 3) {
        rOStream << " at Coodinates " << Coordinates()[0] << " | " << Coordinates()[1] << " | " << Coordinates()[2];
    }
}

void BarycentricLocalSystem::SetPairingStatusForPrinting()
{
    if (mPairingStatus == MapperLocalSystem::PairingStatus::Approximation) {
        mpNode->SetValue(PAIRING_STATUS, 0);
    } else {
        mpNode->SetValue(PAIRING_STATUS, -1);
    }
}

}  // namespace Kratos.
