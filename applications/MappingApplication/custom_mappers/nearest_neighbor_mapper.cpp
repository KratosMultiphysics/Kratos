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
#include "custom_utilities/mapper_typedefs.h"

namespace Kratos
{
/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
void NearestNeighborInterfaceInfo::ProcessSearchResult(const InterfaceObject::Pointer& rpInterfaceObject,
                                                      const double NeighborDistance)
{
    SetLocalSearchWasSuccessful();

    if (NeighborDistance < mNearestNeighborDistance) {
        mNearestNeighborDistance = NeighborDistance;
        mNearestNeighborId = rpInterfaceObject->pGetBaseNode()->GetValue(INTERFACE_EQUATION_ID);
    }
}

void NearestNeighborLocalSystem::CalculateAll(MatrixType& rLocalMappingMatrix,
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
        mInterfaceInfos[0]->GetValue(nearest_neighbor_id);
        mInterfaceInfos[0]->GetValue(nearest_neighbor_distance);

        for (SizeType i=1; i<mInterfaceInfos.size(); ++i) {
            // no check if this InterfaceInfo is an approximation is necessary
            // bcs this does not exist for NearestNeighbor
            double distance;
            mInterfaceInfos[i]->GetValue(distance);

            if (distance < nearest_neighbor_distance) {
                nearest_neighbor_distance = distance;
                mInterfaceInfos[i]->GetValue(nearest_neighbor_id);
            }
        }

        rLocalMappingMatrix(0,0) = 1.0;
        rOriginIds[0] = nearest_neighbor_id;
        KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;
        rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);
    }
    else ResizeToZero(rLocalMappingMatrix, rOriginIds, rDestinationIds, rPairingStatus);
}

std::string NearestNeighborLocalSystem::PairingInfo(const int EchoLevel, const int CommRank) const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

    std::stringstream buffer;
    buffer << "NearestNeighborLocalSystem based on " << mpNode->Info();
    if (EchoLevel > 1) // TODO leave here?
        buffer << " at Coodinates " << Coordinates()[0] << " | " << Coordinates()[1] << " | " << Coordinates()[2];
    buffer << " in rank " << CommRank;
    return buffer.str();
}

/***********************************************************************************/
/* PROTECTED Methods */
/***********************************************************************************/


/***********************************************************************************/
/* PRIVATE Methods */
/***********************************************************************************/

// /// input stream function
// inline std::istream & operator >> (std::istream& rIStream, NearestNeighborMapper& rThis);

// /// output stream function
// inline std::ostream & operator << (std::ostream& rOStream, const NearestNeighborMapper& rThis) {
//   rThis.PrintInfo(rOStream);
//   rOStream << " : " << std::endl;
//   rThis.PrintData(rOStream);
//   return rOStream;
// }

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class NearestNeighborMapper< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
template class NearestNeighborMapper< MapperDefinitions::MPISparseSpaceType, MapperDefinitions::DenseSpaceType >;
#endif

}  // namespace Kratos.
