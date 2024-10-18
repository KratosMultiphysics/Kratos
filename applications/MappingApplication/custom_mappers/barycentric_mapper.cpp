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

// System includes

// External includes

// Project includes
#include "geometries/line_3d_2.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "barycentric_mapper.h"
#include "mapping_application_variables.h"

namespace Kratos {
typedef Node NodeType;
typedef Geometry<NodeType> GeometryType;

namespace { // anonymous namespace

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
    default:
        KRATOS_ERROR << "Wrong interpolation type!" << std::endl;
    };
}

Kratos::unique_ptr<GeometryType> ReconstructLine(const ClosestPointsContainer& rClosestPoints)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rClosestPoints.GetPoints().size() != 2) << "Wrong size! Expected 2 points but got " << rClosestPoints.GetPoints().size() << std::endl;
    GeometryType::PointsArrayType geom_points;
    for (const auto& r_point : rClosestPoints.GetPoints()) {
        auto new_node = Kratos::make_intrusive<NodeType>(0, r_point[0], r_point[1], r_point[2]);
        new_node->SetValue(INTERFACE_EQUATION_ID, r_point.GetId());
        geom_points.push_back(new_node);
    }
    return Kratos::make_unique<Line3D2<NodeType>>(geom_points);

    KRATOS_CATCH("")
}

Kratos::unique_ptr<GeometryType> ReconstructTriangle(const ClosestPointsContainer& rClosestPoints)
{
    KRATOS_TRY

    ClosestPointsContainer points_copy(rClosestPoints);

    while(true) {
        if (points_copy.GetPoints().size() < 3) {
            return ReconstructLine(points_copy);
        }

        GeometryType::PointsArrayType geom_points;
        for (const auto& r_point : points_copy.GetPoints()) {
            auto new_node = Kratos::make_intrusive<NodeType>(0, r_point[0], r_point[1], r_point[2]);
            new_node->SetValue(INTERFACE_EQUATION_ID, r_point.GetId());
            geom_points.push_back(new_node);
            if (geom_points.size() == 3) break;
        }

        // construct triangle from closest points
        auto p_geom = Kratos::make_unique<Triangle3D3<NodeType>>(geom_points);

        // check the quality of the constructed triangle
        const double quality = p_geom->Quality(GeometryType::QualityCriteria::INRADIUS_TO_CIRCUMRADIUS);

        if (quality < 0.05) {
            const double d1 = geom_points[0].Distance(geom_points[1]);
            const double d2 = geom_points[0].Distance(geom_points[2]);
            auto it_point_to_remove = points_copy.GetPoints().begin();
            if (d1 * 10 < d2) {
                // this means that the first and the second point are much closer to each other than the
                // first and the third point
                // hence the second point is removed as the triangle could never be of good condition
                std::advance(it_point_to_remove, 1);
            } else {
                // else remove third point
                std::advance(it_point_to_remove, 2);
            }
            points_copy.GetPoints().erase(it_point_to_remove);
        } else {
            return std::move(p_geom);
        }
    }

    KRATOS_CATCH("")
}

Kratos::unique_ptr<GeometryType> ReconstructTetrahedra(const ClosestPointsContainer& rClosestPoints)
{
    KRATOS_TRY

    ClosestPointsContainer points_copy(rClosestPoints);

    while(true) {
        if (points_copy.GetPoints().size() < 4) {
            return ReconstructTriangle(points_copy);
        }

        GeometryType::PointsArrayType geom_points;
        for (const auto& r_point : points_copy.GetPoints()) {
            auto new_node = Kratos::make_intrusive<NodeType>(0, r_point[0], r_point[1], r_point[2]);
            new_node->SetValue(INTERFACE_EQUATION_ID, r_point.GetId());
            geom_points.push_back(new_node);
            if (geom_points.size() == 4) break;
        }

        // construct triangle from closest points
        auto p_geom = Kratos::make_unique<Tetrahedra3D4<NodeType>>(geom_points);

        // check the quality of the constructed triangle
        const double quality = p_geom->Quality(GeometryType::QualityCriteria::INRADIUS_TO_CIRCUMRADIUS);

        if (quality < 0.05) {
            const double d1 = geom_points[0].Distance(geom_points[1]);
            const double d2 = geom_points[0].Distance(geom_points[2]);
            auto it_point_to_remove = points_copy.GetPoints().begin();
            if (d1 * 10 < d2) {
                // this means that the first and the second point are much closer to each other than the
                // first and the third point
                // hence the second point is removed as the triangle could never be of good condition
                std::advance(it_point_to_remove, 1);
            } else if (MapperUtilities::PointsAreCollinear(geom_points[0], geom_points[1], geom_points[2])) {
                std::advance(it_point_to_remove, 2);
            } else {
                // else remove fourth point
                std::advance(it_point_to_remove, 3);
            }
            points_copy.GetPoints().erase(it_point_to_remove);
        } else {
            return std::move(p_geom);
        }
    }

    KRATOS_CATCH("")
}

Kratos::unique_ptr<GeometryType> ReconstructGeometry(
    const ClosestPointsContainer& rClosestPoints,
    const BarycentricInterpolationType InterpolationType)
{
    KRATOS_TRY

    switch(InterpolationType)
    {
    case BarycentricInterpolationType::LINE:
        return ReconstructLine(rClosestPoints);
    case BarycentricInterpolationType::TRIANGLE:
        return ReconstructTriangle(rClosestPoints);
    case BarycentricInterpolationType::TETRAHEDRA:
        return ReconstructTetrahedra(rClosestPoints);
    default:
        KRATOS_ERROR << "Wrong interpolation type!" << std::endl;
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

    mNumSearchResults++;

    const Node& r_node = *rInterfaceObject.pGetBaseNode();

    PointWithId point(
        r_node.GetValue(INTERFACE_EQUATION_ID),
        r_node.Coordinates(),
        MapperUtilities::ComputeDistance(this->Coordinates(), r_node.Coordinates())
        );

    mClosestPoints.Add(point);

    const int num_found_points = mClosestPoints.GetPoints().size();
    if (num_found_points >= GetNumPointsApprox(mInterpolationType)) {
        SetLocalSearchWasSuccessful();
    } else if (num_found_points > 0) {
        SetIsApproximation();
    }

    KRATOS_CATCH("")
}

void BarycentricInterfaceInfo::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MapperInterfaceInfo );
    rSerializer.save("InterpolationType", static_cast<int>(mInterpolationType));
    rSerializer.save("ClosestPoints", mClosestPoints);
    rSerializer.save("NumSearchResults", mNumSearchResults);
}

void BarycentricInterfaceInfo::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MapperInterfaceInfo );
    int temp;
    rSerializer.load("InterpolationType", temp);
    mInterpolationType = static_cast<BarycentricInterpolationType>(temp);
    rSerializer.load("ClosestPoints", mClosestPoints);
    rSerializer.load("NumSearchResults", mNumSearchResults);
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

    for (std::size_t i=1; i<mInterfaceInfos.size(); ++i) {
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

    Vector sf_values;
    double proj_dist;

    const bool is_full_projection = ProjectionUtilities::ComputeProjection(
        *p_reconstr_geom,
        Point(Coordinates()),
        0.25,
        sf_values,
        rOriginIds,
        proj_dist,
        mPairingIndex,
        true);

    rPairingStatus = is_full_projection ? MapperLocalSystem::PairingStatus::InterfaceInfoFound : MapperLocalSystem::PairingStatus::Approximation;

    if (rPairingStatus == MapperLocalSystem::PairingStatus::InterfaceInfoFound) {
        const auto interpol_type = r_first_info.GetInterpolationType();
        if ((interpol_type == BarycentricInterpolationType::LINE       && p_reconstr_geom->PointsNumber() != 2) ||
            (interpol_type == BarycentricInterpolationType::TRIANGLE   && p_reconstr_geom->PointsNumber() != 3) ||
            (interpol_type == BarycentricInterpolationType::TETRAHEDRA && p_reconstr_geom->PointsNumber() != 4)) {
            rPairingStatus = MapperLocalSystem::PairingStatus::Approximation;
        }
    }

    if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != sf_values.size()) {
        rLocalMappingMatrix.resize(1, sf_values.size(), false);
    }
    for (IndexType i=0; i<sf_values.size(); ++i) {
        rLocalMappingMatrix(0,i) = sf_values[i];
    }

    KRATOS_CATCH("")
}

void BarycentricLocalSystem::PairingInfo(std::ostream& rOStream, const int EchoLevel) const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

    rOStream << "BarycentricLocalSystem based on " << mpNode->Info();
    if (EchoLevel > 3) {
        rOStream << " at Coordinates " << Coordinates()[0] << " | " << Coordinates()[1] << " | " << Coordinates()[2];
    }
}

void BarycentricLocalSystem::SetPairingStatusForPrinting()
{
    if (mPairingStatus == MapperLocalSystem::PairingStatus::Approximation) {
        mpNode->SetValue(PAIRING_STATUS, (int)mPairingIndex);
    }
}

bool BarycentricLocalSystem::IsDoneSearching() const
{
    if (HasInterfaceInfoThatIsNotAnApproximation()) {return true;};

    if (mInterfaceInfos.empty()) {return false;}

    // collect the results from all partitions and check if enough points were found
    const BarycentricInterfaceInfo& r_first_info = static_cast<const BarycentricInterfaceInfo&>(*mInterfaceInfos[0]);
    const std::size_t num_approx = GetNumPointsApprox(r_first_info.GetInterpolationType());
    std::size_t sum_search_results = r_first_info.GetNumSearchResults();

    return sum_search_results > num_approx*2;
}

}  // namespace Kratos.
