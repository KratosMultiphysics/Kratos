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
void NearestNeigborInterfaceInfo::ProcessSearchResult(const InterfaceObject::Pointer& rpInterfaceObject,
                                                      const double NeighborDistance)
{
    SetLocalSearchWasSuccessful();

    if (NeighborDistance < mNearestNeighborDistance)
    {
        mNearestNeighborDistance = NeighborDistance;
        mNearestNeighborId = rpInterfaceObject->pGetBaseNode()->GetValue(INTERFACE_EQUATION_ID);
    }
}

void NearestNeighborLocalSystem::CalculateAll(MappingWeightsVector& rMappingWeights,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds) const
{
    if (rMappingWeights.size() != 1) rMappingWeights.resize(1);
    if (rOriginIds.size()      != 1) rOriginIds.resize(1);
    if (rDestinationIds.size() != 1) rDestinationIds.resize(1);

    if (mInterfaceInfos.size() > 0)
    {
        int nearest_neighbor_id;
        double nearest_neighbor_distance;
        mInterfaceInfos[0]->GetValue(nearest_neighbor_id);
        mInterfaceInfos[0]->GetValue(nearest_neighbor_distance);

        for (SizeType i=1; i<mInterfaceInfos.size(); ++i)
        {
            double distance;
            mInterfaceInfos[i]->GetValue(distance);

            if (distance < nearest_neighbor_distance)
            {
                nearest_neighbor_distance = distance;
                mInterfaceInfos[i]->GetValue(nearest_neighbor_id);
            }
        }

        rMappingWeights[0] = 1.0;
        rOriginIds[0] = nearest_neighbor_id;
        rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);
    }
    else
    {
        KRATOS_WARNING("NearestNeighborMapper")
            << "MapperLocalSystem No xxx" << "xxx" << " has not found a neighbor" << std::endl;

        // TODO is this ok? => I guess it would be better to do this in a more general way in the baseclass...
        // TODO resize to zero, then it wont be assembled! (might be a bit slower though...)
        rMappingWeights[0] = 0.0;
        rOriginIds[0]      = 0;
        rDestinationIds[0] = 0;
    }

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
