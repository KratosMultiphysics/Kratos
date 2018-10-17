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

typedef std::size_t IndexType;
typedef std::size_t SizeType;

typedef Node<3> NodeType;
typedef Geometry<NodeType> GeometryType;
typedef Kratos::shared_ptr<GeometryType> GeometryPointerType;
typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;


bool ProjectTo1D2D(const GeometryPointerType pGeometry,
                   const Point& rPointToProject,
                   CoordinatesArrayType& rLocalCoords,
                   double& rDistance)
{
    CoordinatesArrayType local_coords_init;

    Point projected_point;

    // using the center as trial for the projection
    pGeometry->PointLocalCoordinates(local_coords_init, pGeometry->Center());

    // trying to project to the geometry
    rDistance = std::abs(GeometricalProjectionUtilities::FastProjectDirection(
        *pGeometry,
        rPointToProject,
        projected_point,
        pGeometry->UnitNormal(local_coords_init),
        pGeometry->UnitNormal(local_coords_init)));

    bool is_inside = pGeometry->IsInside(projected_point, rLocalCoords);

    return is_inside;
}

bool ProjectIntoVolume(const GeometryPointerType pGeometry,
                       const Point& rPointToProject,
                       CoordinatesArrayType& rLocalCoords,
                       double& rDistance)
{
    bool is_inside = pGeometry->IsInside(rPointToProject, rLocalCoords);

    if (is_inside) {
        // Calculate Distance
        rDistance = MapperUtilities::ComputeDistance(rPointToProject, pGeometry->Center());
        rDistance /= pGeometry->Volume(); // Normalize Distance by Volume
    }

    return is_inside;
}

/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
void NearestElementInterfaceInfo::ProcessSearchResult(const InterfaceObject::Pointer& rpInterfaceObject,
                                                      const double NeighborDistance)
{
    const auto& p_geom = rpInterfaceObject->pGetBaseGeometry();
    const SizeType num_nodes = p_geom->PointsNumber();

    const auto geom_family = p_geom->GetGeometryFamily();

    bool is_inside;
    double proj_dist;
    CoordinatesArrayType local_coords;

    const Point point_to_proj(this->Coordinates());

    // select projection depending on type of geometry
    if ((geom_family == GeometryData::Kratos_Linear        && num_nodes == 2) || // linear line
        (geom_family == GeometryData::Kratos_Triangle      && num_nodes == 3) || // linear triangle
        (geom_family == GeometryData::Kratos_Quadrilateral && num_nodes == 4)) { // linear quad
        is_inside = ProjectTo1D2D(p_geom, point_to_proj, local_coords, proj_dist);
    }
    else if (geom_family == GeometryData::Kratos_Tetrahedra ||
             geom_family == GeometryData::Kratos_Prism ||
             geom_family == GeometryData::Kratos_Hexahedra) { // Volume projection
        is_inside = ProjectIntoVolume(p_geom, point_to_proj, local_coords, proj_dist);
    }
    else {
        is_inside = false;
        KRATOS_WARNING_ONCE("NearestElementMapper") << "Unsupported type of geometry,"
            << "trying to use an approximation (Nearest Neighbor)" << std::endl;
    }

    Vector shape_function_values;

    // if it is closer, then we update the members to make this geometry the closest projection
    if (is_inside && proj_dist < mClosestProjectionDistance) {
        SetLocalSearchWasSuccessful();
        mClosestProjectionDistance = proj_dist;
        mShapeFunctionValues.clear();
        mNodeIds.clear();

        p_geom->ShapeFunctionsValues(shape_function_values, local_coords);
        KRATOS_DEBUG_ERROR_IF_NOT(shape_function_values.size() == num_nodes)
            << "Number of SFs is different from number of nodes!" << std::endl;

        if (mShapeFunctionValues.size() != num_nodes) mShapeFunctionValues.resize(num_nodes);
        if (mNodeIds.size() != num_nodes)             mNodeIds.resize(num_nodes);
        for (IndexType i=0; i<num_nodes; ++i) {
            mShapeFunctionValues[i] = shape_function_values[i];
            KRATOS_DEBUG_ERROR_IF_NOT((*p_geom)[i].Has(INTERFACE_EQUATION_ID))
                << "Node #" << (*p_geom)[i].Id() << " does not have an Interface Id!\n" << (*p_geom)<< "\n"
                <<  (*p_geom)[i] << std::endl;
            mNodeIds[i] = (*p_geom)[i].GetValue(INTERFACE_EQUATION_ID);
        }
    }
}

void NearestElementInterfaceInfo::ProcessSearchResultForApproximation(const InterfaceObject::Pointer& rpInterfaceObject,
                                                                      const double NeighborDistance)
{
    const auto& p_geom = rpInterfaceObject->pGetBaseGeometry();

    // looping the points of the geometry and finding the nearest neighbor
    for (const auto& r_point : p_geom->Points()) {
        const double dist = MapperUtilities::ComputeDistance(this->Coordinates(), r_point.Coordinates());

        // in case of an approximation this is the actual distance,
        // not the projected one bcs no valid projection could be found!
        if (dist < mClosestProjectionDistance) {
            mClosestProjectionDistance = dist;
            if (mNodeIds.size() != 1) mNodeIds.resize(1);
            if (mShapeFunctionValues.size() != 1) mShapeFunctionValues.resize(1);
            KRATOS_DEBUG_ERROR_IF_NOT(r_point.Has(INTERFACE_EQUATION_ID));
            mNodeIds[0] = r_point.GetValue(INTERFACE_EQUATION_ID);
            mShapeFunctionValues[0] = 1.0; // Approximation is nearest node
        }
    }

    SetIsApproximation();
}

void NearestElementLocalSystem::CalculateAll(MatrixType& rLocalMappingMatrix,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds,
                    MapperLocalSystem::PairingStatus& rPairingStatus) const
{
    if (mInterfaceInfos.size() > 0) {
        double distance;
        double min_distance = std::numeric_limits<double>::max();
        int found_idx = -1;
        for (IndexType i=0; i<mInterfaceInfos.size(); ++i) {
            // the approximations will be processed in the next step if necessary
            if (!mInterfaceInfos[i]->GetIsApproximation()) {
                mInterfaceInfos[i]->GetValue(distance);
                if (distance < min_distance) {
                    min_distance = distance;
                    found_idx = static_cast<int>(i); // TODO explicit conversion needed?
                    rPairingStatus = MapperLocalSystem::PairingStatus::InterfaceInfoFound;
                }
            }
        }

        if (found_idx == -1) { // no valid projection exists => using an approximation
            for (IndexType i=0; i<mInterfaceInfos.size(); ++i) {
                // now the approximations are being checked
                if (mInterfaceInfos[i]->GetIsApproximation()) {
                    mInterfaceInfos[i]->GetValue(distance);
                    if (distance < min_distance) {
                        min_distance = distance;
                        found_idx = static_cast<int>(i); // TODO explicit conversion needed?
                        rPairingStatus = MapperLocalSystem::PairingStatus::Approximation;
                    }
                }
            }
        }

        KRATOS_ERROR_IF(found_idx == -1) << "Not even an approximation is found, this should not happen!"
            << std::endl; // TODO should thi sbe an error?

        std::vector<double> sf_values;

        mInterfaceInfos[found_idx]->GetValue(sf_values);

        if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != sf_values.size()) {
            rLocalMappingMatrix.resize(1, sf_values.size(), false);
        }
        for (IndexType i=0; i<sf_values.size(); ++i) {
            rLocalMappingMatrix(0,i) = sf_values[i];
        }

        mInterfaceInfos[found_idx]->GetValue(rOriginIds);

        KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

        if (rDestinationIds.size() != 1) rDestinationIds.resize(1);
        rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);
    }
    else ResizeToZero(rLocalMappingMatrix, rOriginIds, rDestinationIds, rPairingStatus);
}

std::string NearestElementLocalSystem::PairingInfo(const int EchoLevel, const int CommRank) const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

    std::stringstream buffer;
    buffer << "NearestElementLocalSystem based on " << mpNode->Info();
    if (EchoLevel > 1) // TODO leave here?
        buffer << " at Coodinates " << Coordinates()[0] << " | " << Coordinates()[1] << " | " << Coordinates()[2];
    buffer << " in rank " << CommRank;
    return buffer.str();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class NearestElementMapper< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
template class NearestElementMapper< MapperDefinitions::MPISparseSpaceType, MapperDefinitions::DenseSpaceType >;
#endif

}  // namespace Kratos.
