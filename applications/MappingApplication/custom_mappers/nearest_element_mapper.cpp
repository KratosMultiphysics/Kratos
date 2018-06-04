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
#include "nearest_element_mapper.h"
#include "custom_utilities/mapper_typedefs.h"
#include "utilities/mortar_utilities.h" // for projections

namespace Kratos
{
/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
void NearestElementInterfaceInfo::ProcessSearchResult(const InterfaceObject::Pointer& rpInterfaceObject,
                                                      const double NeighborDistance)
{
    const auto& p_geom = rpInterfaceObject->pGetBaseGeometry();

    // trying to project to the geometry
    // ...

        // checking whether the projection is inside or outside
        // ...

            // SetLocalSearchWasSuccessful();
            // checking whether this geometry gives a closer projection

                // if it is closer, then we update the members to make this geometry the closest projection
                // Store SF-values
                // Store INTERFACE_EQUATION_ID s of the nodes of the geometry
}



void NearestElementInterfaceInfo::ProcessSearchResultForApproximation(const InterfaceObject::Pointer& rpInterfaceObject,
                                                                      const double NeighborDistance)
{
    const auto& p_geom = rpInterfaceObject->pGetBaseGeometry();

    // looping the points of the geometry and finding the nearest neighbor

    // SetIsApproximation()
}

void NearestElementLocalSystem::CalculateAll(MappingWeightsVector& rMappingWeights,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds) const
{
    // if (rMappingWeights.size() != 1) rMappingWeights.resize(1);
    // if (rOriginIds.size() != 1)      rOriginIds.resize(1);
    // if (rDestinationIds.size() != 1) rDestinationIds.resize(1);

    // if (this->mInterfaceInfos.size() > 0)
    // {
    //     // IndexType nearest_neighbor_id = this->mInterfaceInfos[0]->GetInterfaceData().GetNearestNeighborId();
    //     // double nearest_neighbor_distance = this->mInterfaceInfos[0]->GetInterfaceData().GetNearestNeighborDistance();

    //     // for (SizeType i=1; i<this->mInterfaceInfos.size(); ++i)
    //     // {
    //     //     TDataHolder data = this->mInterfaceInfos[i]->GetInterfaceData();
    //     //     double distance = data.GetNearestNeighborDistance();

    //     //     if (distance < nearest_neighbor_distance)
    //     //     {
    //     //         nearest_neighbor_distance = distance;
    //     //         nearest_neighbor_id = data.GetNearestNeighborId();
    //     //     }
    //     // }

    //     rMappingWeights[0] = 1.0;
    //     // rOriginIds[0] = nearest_neighbor_id;
    //     // rDestinationIds[0] = this->mrNode.GetValue(INTERFACE_EQUATION_ID); //TODO
    // }
    // else
    // {
    //     KRATOS_WARNING_IF("NearestNeighborMapper", this->mInterfaceInfos.size() == 0)
    //         << "MapperLocalSystem No xxx" << "xxx" << " has not found a neighbor" << std::endl;

    //     // TODO is this ok? => I guess it would be better to do this in a mor general way in the baseclass...
    //     // TODO resize to zero, then it wont be assembled! (might be a bit slower though...)
    //     rMappingWeights[0] = 0.0;
    //     rOriginIds[0]      = 0;
    //     rDestinationIds[0] = 0;
    // }

}

/***********************************************************************************/
/* PROTECTED Methods */
/***********************************************************************************/


/***********************************************************************************/
/* PRIVATE Methods */
/***********************************************************************************/


// /// input stream function
// inline std::istream & operator >> (std::istream& rIStream, NearestElementMapper& rThis);

// /// output stream function
// inline std::ostream & operator << (std::ostream& rOStream, const NearestElementMapper& rThis) {
//   rThis.PrintInfo(rOStream);
//   rOStream << " : " << std::endl;
//   rThis.PrintData(rOStream);
//   return rOStream;
// }

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class NearestElementMapper< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
template class NearestElementMapper< MapperDefinitions::MPISparseSpaceType, MapperDefinitions::DenseSpaceType >;
#endif

}  // namespace Kratos.
