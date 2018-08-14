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
#include "custom_utilities/mapper_utilities.h"
#include "utilities/geometrical_projection_utilities.h"

namespace Kratos
{
/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
void NearestElementInterfaceInfo::ProcessSearchResult(const InterfaceObject::Pointer& rpInterfaceObject,
                                                      const double NeighborDistance)
{
    const auto& p_geom = rpInterfaceObject->pGetBaseGeometry();

    Point proj_point;
    Point point_to_proj(this->Coordinates());

    typedef typename Geometry<Node<3>>::CoordinatesArrayType CoordinatesArrayType;

    CoordinatesArrayType local_coords_init;
    CoordinatesArrayType local_coords;

    p_geom->PointLocalCoordinates(local_coords_init, (*p_geom)[0]);

    // // trying to project to the geometry
    const double proj_dist = GeometricalProjectionUtilities::FastProjectDirection(
        *p_geom, point_to_proj, proj_point, p_geom->UnitNormal(local_coords_init), p_geom->UnitNormal(local_coords_init));

    const bool is_inside = p_geom->IsInside(proj_point, local_coords);

    // if it is closer, then we update the members to make this geometry the closest projection
    if (is_inside && proj_dist < mClosestProjectionDistance)
    {
        SetLocalSearchWasSuccessful();
        mClosestProjectionDistance = proj_dist;
        mShapeFunctionValues.clear();
        mNodeIds.clear();
        mShapeFunctionValues.resize(local_coords.size());
        mNodeIds.resize(local_coords.size());
        for (IndexType i=0; i<local_coords.size(); ++i)
        {
            mShapeFunctionValues[i] = local_coords[i];
            mNodeIds[i] = (*p_geom)[i].GetValue(INTERFACE_EQUATION_ID);
        }
    }
}


void NearestElementInterfaceInfo::ProcessSearchResultForApproximation(const InterfaceObject::Pointer& rpInterfaceObject,
                                                                      const double NeighborDistance)
{
    const auto& p_geom = rpInterfaceObject->pGetBaseGeometry();

    // looping the points of the geometry and finding the nearest neighbor
    for (const auto& r_point : p_geom->Points())
    {
        const double dist = MapperUtilities::ComputeDistance(this->Coordinates(), r_point.Coordinates());

        if (dist < mClosestProjectionDistance)
        {
            mClosestProjectionDistance = dist;
            if (mNodeIds.size() != 1) mNodeIds.resize(1);
            if (mShapeFunctionValues.size() != 1) mShapeFunctionValues.resize(1);
            mNodeIds[0] = r_point.GetValue(INTERFACE_EQUATION_ID);
            mShapeFunctionValues[0] = 1.0; // Approximation is nearest node
        }
    }

    SetIsApproximation();
}

void NearestElementLocalSystem::CalculateAll(MappingWeightsVector& rMappingWeights,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds,
                    MapperLocalSystem::PairingStatus& rPairingStatus) const
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
