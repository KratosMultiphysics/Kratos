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
#include "geometries/line_2d_2.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "barycentric_mapper.h"
#include "mapping_application_variables.h"

namespace Kratos {
typedef Node<3> NodeType;
typedef Geometry<NodeType> GeometryType;

namespace { // anonymous namespace

// array_1d<double, 3> CreateArrayFromVector(const std::vector<double>& rVector,
//                                           const std::size_t StartPosition)
// {
//     array_1d<double,3> the_array;
//     the_array[0] = rVector[StartPosition];
//     the_array[1] = rVector[StartPosition+1];
//     the_array[2] = rVector[StartPosition+2];
//     return the_array;
// }

// void ComputeNeighborDistances(const array_1d<double,3>& rCoords,
//                               const std::vector<double>& rNeighborCoords,
//                               std::vector<double>& rDistances)
// {
//     for (std::size_t i=0; i<rDistances.size(); ++i) {
//         const auto neighbor_coords = CreateArrayFromVector(rNeighborCoords, i*3);
//         rDistances[i] = MapperUtilities::ComputeDistance(rCoords, neighbor_coords);
//     }
// }

// void InsertIfBetter(const array_1d<double,3>& rRefCoords,
//                     const int CandidateEquationId,
//                     const array_1d<double,3>& rCandidateCoords,
//                     std::vector<int>& rNeighborIds,
//                     std::vector<double>& rNeighborCoods)
// {
//     const std::size_t num_interpolation_nodes = rNeighborIds.size();

//     std::vector<double> neighbor_distances(num_interpolation_nodes, std::numeric_limits<double>::max());

//     // compute the distances to the currently closest nodes to the candidate to check for coinciding nodes
//     ComputeNeighborDistances(rCandidateCoords, rNeighborCoods, neighbor_distances);
//     for (const double dist : neighbor_distances) {
//         if (dist < 1e-12) return; // the candidate is coinciding with an already used point, hence we don't use it
//     }

//     // compute the distances to the currently closest nodes based on the coordinates
//     ComputeNeighborDistances(rRefCoords, rNeighborCoods, neighbor_distances);

//     // compute distance to candidate
//     const double candidate_distance = MapperUtilities::ComputeDistance(rRefCoords, rCandidateCoords);

//     // check if the candidate is closer than the previously found nodes and save it if it is closer
//     std::vector<int> temp_neighbor_ids(rNeighborIds);
//     std::vector<double> temp_neighbor_coords(rNeighborCoods);

//     for (std::size_t i=0; i<num_interpolation_nodes; ++i) {
//         if (candidate_distance < neighbor_distances[i]) {
//             temp_neighbor_ids.insert(temp_neighbor_ids.begin()+i, CandidateEquationId);
//             temp_neighbor_coords.insert(temp_neighbor_coords.begin()+(i*3), std::begin(rCandidateCoords), std::end(rCandidateCoords));
//             break;
//         }
//     }

//     // check the quality of the resulting geometry
//     if (num_interpolation_nodes == 3) {
//         if (std::count(temp_neighbor_ids.begin(), temp_neighbor_ids.begin()+num_interpolation_nodes, -1) == 0) { // make sure enough points were found
//             GeometryType::PointsArrayType geom_points;
//             for (std::size_t i=0; i<num_interpolation_nodes; ++i) {
//                 geom_points.push_back(Kratos::make_intrusive<NodeType>(0, temp_neighbor_coords[i*3], temp_neighbor_coords[i*3+1], temp_neighbor_coords[i*3+2]));
//             }
//             Kratos::unique_ptr<GeometryType> p_geom = Kratos::make_unique<Triangle3D3<NodeType>>(geom_points);

//             const double quality(p_geom->Quality(GeometryType::QualityCriteria::INRADIUS_TO_CIRCUMRADIUS));
//             if (quality < 0.05) { // TODO check if this is suitable ...
//                 KRATOS_INFO("Bad Quality") << quality << std::endl;
//                 // removing the 3rd closest point
//                 temp_neighbor_ids.erase(temp_neighbor_ids.begin()+2);
//                 temp_neighbor_coords.erase(temp_neighbor_coords.begin()+6, temp_neighbor_coords.begin()+9);
//             }
//         }
//     } else if (num_interpolation_nodes == 4) {

//         if (std::count(temp_neighbor_ids.begin(), temp_neighbor_ids.begin()+num_interpolation_nodes, -1) == 0) { // make sure enough points were found
//             GeometryType::PointsArrayType geom_points;
//             for (std::size_t i=0; i<num_interpolation_nodes; ++i) {
//                 geom_points.push_back(Kratos::make_intrusive<NodeType>(0, temp_neighbor_coords[i*3], temp_neighbor_coords[i*3+1], temp_neighbor_coords[i*3+2]));
//             }
//             Kratos::unique_ptr<GeometryType> p_geom = Kratos::make_unique<Tetrahedra3D4<NodeType>>(geom_points);
//             std::cout << "\n  >>> HHHEEERRREEE <<<" << std::endl;

//             // TODO I think I have to check first the Tri (only when filling up with points, after that should be enough to only check the tetra(...?)) and then the Tetra ...?
//         }
//     }

//     // copy back the final results
//     for (std::size_t i=0; i<num_interpolation_nodes; ++i) {
//         rNeighborIds[i] = temp_neighbor_ids[i];
//         rNeighborCoods[(i*3)]   = temp_neighbor_coords[(i*3)];
//         rNeighborCoods[(i*3)+1] = temp_neighbor_coords[(i*3)+1];
//         rNeighborCoods[(i*3)+2] = temp_neighbor_coords[(i*3)+2];
//     }
// }

bool BarycentricInterpolateInEntity(const array_1d<double,3>& rRefCoords,
                                    const std::vector<double>& rCoordinates,
                                    Vector& rShapeFunctionValues,
                                    std::vector<int>& rEquationIds,
                                    ProjectionUtilities::PairingIndex& rPairingIndex)
{
    // Check how many "proper" results were found
    const std::size_t num_interpolation_nodes = rEquationIds.size() - std::count(rEquationIds.begin(), rEquationIds.end(), -1);

    const bool is_full_projection = num_interpolation_nodes == rEquationIds.size();

    KRATOS_DEBUG_ERROR_IF(num_interpolation_nodes < 2 || num_interpolation_nodes > 4) << "Wrong number of interpolation nodes" << std::endl;
    KRATOS_DEBUG_ERROR_IF(rCoordinates.size() < num_interpolation_nodes*3) << "Not enough coords" << std::endl;
    KRATOS_DEBUG_ERROR_IF(rCoordinates.size()%3 != 0) << "Coords have wrong size" << std::endl;

    GeometryType::PointsArrayType geom_points;
    for (std::size_t i=0; i<num_interpolation_nodes; ++i) {
        geom_points.push_back(Kratos::make_intrusive<NodeType>(0, rCoordinates[i*3], rCoordinates[i*3+1], rCoordinates[i*3+2]));
        geom_points[i].SetValue(INTERFACE_EQUATION_ID, rEquationIds[i]);
    }

    Kratos::unique_ptr<GeometryType> p_geom;
    if      (num_interpolation_nodes == 2) p_geom = Kratos::make_unique<Line2D2<NodeType>>(geom_points);
    else if (num_interpolation_nodes == 3) p_geom = Kratos::make_unique<Triangle3D3<NodeType>>(geom_points);
    else if (num_interpolation_nodes == 4) p_geom = Kratos::make_unique<Tetrahedra3D4<NodeType>>(geom_points);

    double dummy_dist;

    return is_full_projection && ProjectionUtilities::ComputeProjection(*p_geom, Point(rRefCoords), 0.25, rShapeFunctionValues, rEquationIds, dummy_dist, rPairingIndex, true);
}

int GetNumPointsApprox(const BarycentricInterpolationType InterpolationType)
{
    switch(InterpolationType)
    {
    case BarycentricInterpolationType::LINE:
        return 2;
    case BarycentricInterpolationType::TRIANGLE:
        return 10;
    case BarycentricInterpolationType::TETRAHEDRA:
        return 15;
    };
}

Kratos::unique_ptr<GeometryType> ReconstructGeometry(
    const ClosestPointsContainer& rClosestPoints,
    const BarycentricInterpolationType InterpolationType)
{
    KRATOS_TRY

    switch(InterpolationType)
    {
    case BarycentricInterpolationType::LINE:
        KRATOS_ERROR_IF(rClosestPoints.GetPoints().size() != 2) << "Wrong size!" << std::endl;
        GeometryType::PointsArrayType geom_points;
        for (const auto& r_point : rClosestPoints.GetPoints()) {
            auto new_node = Kratos::make_intrusive<NodeType>(0, r_point[0], r_point[1], r_point[2]);
            new_node->SetValue(INTERFACE_EQUATION_ID, r_point.GetId());
            geom_points.push_back(new_node);
        }
        return Kratos::make_unique<Line2D2<NodeType>>(geom_points);

    // case BarycentricInterpolationType::TRIANGLE:
    //     return 10;

    // case BarycentricInterpolationType::TETRAHEDRA:
    //     return 15;
    };

    KRATOS_CATCH("")
}

} // anonymous namespace


BarycentricInterfaceInfo::BarycentricInterfaceInfo(const BarycentricInterpolationType InterpolationType)
    : mInterpolationType(InterpolationType), mClosestPoints(GetNumPointsApprox(InterpolationType)) {}

BarycentricInterfaceInfo::BarycentricInterfaceInfo(const CoordinatesArrayType& rCoordinates,
                                    const IndexType SourceLocalSystemIndex,
                                    const IndexType SourceRank,
                                    const BarycentricInterpolationType InterpolationType)
    : MapperInterfaceInfo(rCoordinates, SourceLocalSystemIndex, SourceRank),
      mInterpolationType(InterpolationType),
      mClosestPoints(GetNumPointsApprox(InterpolationType)) {}

MapperInterfaceInfo::Pointer BarycentricInterfaceInfo::Create() const
{
    return Kratos::make_shared<BarycentricInterfaceInfo>(mInterpolationType);
}

MapperInterfaceInfo::Pointer BarycentricInterfaceInfo::Create(const CoordinatesArrayType& rCoordinates,
                                    const IndexType SourceLocalSystemIndex,
                                    const IndexType SourceRank) const
{
    return Kratos::make_shared<BarycentricInterfaceInfo>(
        rCoordinates,
        SourceLocalSystemIndex,
        SourceRank,
        mInterpolationType);
}

InterfaceObject::ConstructionType BarycentricInterfaceInfo::GetInterfaceObjectType() const
{
    return InterfaceObject::ConstructionType::Node_Coords;
}

void BarycentricInterfaceInfo::ProcessSearchResult(const InterfaceObject& rInterfaceObject)
{
    KRATOS_TRY

    SetLocalSearchWasSuccessful();

    const Node<3>& r_node = *rInterfaceObject.pGetBaseNode();

    PointWithId point(
        r_node.GetValue(INTERFACE_EQUATION_ID),
        r_node.Coordinates(),
        MapperUtilities::ComputeDistance(this->Coordinates(), r_node.Coordinates()));

    mClosestPoints.Add(point);

    KRATOS_CATCH("")
}

void BarycentricInterfaceInfo::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MapperInterfaceInfo );
    rSerializer.save("InterpolationType", static_cast<int>(mInterpolationType));
    rSerializer.save("ClosestPoints", mClosestPoints);
}

void BarycentricInterfaceInfo::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MapperInterfaceInfo );
    int temp;
    rSerializer.load("InterpolationType", temp);
    mInterpolationType = static_cast<BarycentricInterpolationType>(temp);
    rSerializer.load("ClosestPoints", mClosestPoints);
}


void BarycentricLocalSystem::CalculateAll(MatrixType& rLocalMappingMatrix,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds,
                    MapperLocalSystem::PairingStatus& rPairingStatus) const
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

    if (mInterfaceInfos.size() < 1) {
        ResizeToZero(rLocalMappingMatrix, rOriginIds, rDestinationIds, rPairingStatus);
        return;
    }

    const BarycentricInterfaceInfo& r_first_info = static_cast<const BarycentricInterfaceInfo&>(*mInterfaceInfos[0]);
    ClosestPointsContainer closest_points(GetNumPointsApprox(r_first_info.GetInterpolationType()));
    closest_points.Merge(r_first_info.GetClosestPoints());

    for (std::size_t i=0; i<mInterfaceInfos.size(); ++i) {
        // in MPI there might be results also from other ranks
        const BarycentricInterfaceInfo& r_info = static_cast<const BarycentricInterfaceInfo&>(*mInterfaceInfos[i]);
        closest_points.Merge(r_info.GetClosestPoints());
    }

    KRATOS_ERROR_IF(closest_points.GetPoints().size() == 0) << "Not even an approximation is found, this should not happen!" << std::endl; // TODO should this be an error?

    if (rDestinationIds.size() != 1) rDestinationIds.resize(1);
    rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);

    if (closest_points.GetPoints().size() == 1) { // only one node was found => using nearest-neighbor
        rPairingStatus = MapperLocalSystem::PairingStatus::Approximation;
        mPairingIndex = ProjectionUtilities::PairingIndex::Closest_Point;
        if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != 1) rLocalMappingMatrix.resize(1, 1, false);
        rLocalMappingMatrix(0,0) = 1.0;
        if (rOriginIds.size() != 1) rOriginIds.resize(1);
        rOriginIds[0] = closest_points.GetPoints().begin()->GetId();
        return;
    }

    Kratos::unique_ptr<GeometryType> p_reconstr_geom = ReconstructGeometry(closest_points, r_first_info.GetInterpolationType());

    // bool ComputeProjection(const GeometryType& rGeometry,
    //                    const Point& rPointToProject,
    //                    const double LocalCoordTol,
    //                    Vector& rShapeFunctionValues,
    //                    std::vector<int>& rEquationIds,
    //                    double& rProjectionDistance,
    //                    PairingIndex& rPairingIndex,
    //                    const bool ComputeApproximation)

    Vector sf_values;
    double proj_dist;

    const bool is_approx = ProjectionUtilities::ComputeProjection(
        *p_reconstr_geom,
        Point(Coordinates()),
        0.1,
        sf_values,
        rOriginIds,
        proj_dist,
        mPairingIndex,
        true);

    if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != sf_values.size()) {
        rLocalMappingMatrix.resize(1, sf_values.size(), false);
    }
    for (IndexType i=0; i<sf_values.size(); ++i) {
        rLocalMappingMatrix(0,i) = sf_values[i];
    }


    // std::vector<int> node_ids;
    // std::vector<double> neighbor_coods;

    // // allocate final vectors, using the max possible size (in case of a volume interpolation)
    // std::vector<int> final_node_ids(4);
    // std::vector<double> final_neighbor_coords(12);
    // std::fill(final_node_ids.begin(), final_node_ids.end(), -1);
    // std::fill(final_neighbor_coords.begin(), final_neighbor_coords.end(), std::numeric_limits<double>::max());

    // for (std::size_t i=0; i<mInterfaceInfos.size(); ++i) {
    //     mInterfaceInfos[i]->GetValue(node_ids, MapperInterfaceInfo::InfoType::Dummy);
    //     mInterfaceInfos[i]->GetValue(neighbor_coods, MapperInterfaceInfo::InfoType::Dummy);

    //     const std::size_t num_interpolation_nodes = node_ids.size();

    //     if (final_node_ids.size() != num_interpolation_nodes) {
    //         final_node_ids.resize(num_interpolation_nodes);
    //     }

    //     for (std::size_t j=0; j<num_interpolation_nodes; ++j) {
    //         InsertIfBetter(
    //             Coordinates(),
    //             node_ids[j],
    //             CreateArrayFromVector(neighbor_coods, j*3),
    //             final_node_ids,
    //             final_neighbor_coords);
    //     }
    // }



    // const std::size_t num_interpolation_nodes = final_node_ids.size();

    // KRATOS_DEBUG_ERROR_IF(num_interpolation_nodes < 2 || num_interpolation_nodes > 4) << "Wrong number of interpolation nodes" << std::endl;



    // if (final_node_ids[1] == -1) { // only one node was found => using nearest-neighbor
    //     rPairingStatus = MapperLocalSystem::PairingStatus::Approximation;
    //     mPairingIndex = ProjectionUtilities::PairingIndex::Closest_Point;
    //     if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != 1) rLocalMappingMatrix.resize(1, 1, false);
    //     rLocalMappingMatrix(0,0) = 1.0;
    //     if (rOriginIds.size() != 1) rOriginIds.resize(1);
    //     rOriginIds[0] = final_node_ids[0];
    // } else {
    //     Vector shape_function_values;
    //     const bool is_full_projection = BarycentricInterpolateInEntity(Coordinates(), final_neighbor_coords, shape_function_values, final_node_ids, mPairingIndex);
    //     rPairingStatus = is_full_projection ? MapperLocalSystem::PairingStatus::InterfaceInfoFound : MapperLocalSystem::PairingStatus::Approximation;

    //     if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != shape_function_values.size()) {
    //         rLocalMappingMatrix.resize(1, shape_function_values.size(), false);
    //     }
    //     for (std::size_t i=0; i<shape_function_values.size(); ++i) {
    //         rLocalMappingMatrix(0,i) = shape_function_values[i];
    //     }

    //     if (rOriginIds.size() != final_node_ids.size()) rOriginIds.resize(final_node_ids.size());
    //     for (std::size_t i=0; i<final_node_ids.size(); ++i) {
    //         rOriginIds[i] = final_node_ids[i];
    //     }
    // }

    KRATOS_CATCH("")
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
        mpNode->SetValue(PAIRING_STATUS, (int)mPairingIndex);
    }
}

}  // namespace Kratos.
