//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#ifndef KRATOS_MPM_SEARCH_ELEMENT_UTILITY
#define KRATOS_MPM_SEARCH_ELEMENT_UTILITY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/quadrature_points_utility.h"

#include "particle_mechanics_application_variables.h"
#include "geometries/geometry_shape_function_container.h"
#include "custom_geometries/quadrature_point_partitioned_geometry.h"

#include "geometries/geometry.h"
#include "includes/model_part.h"

#include "boost/geometry/geometry.hpp"
#include "boost/geometry/geometries/register/point.hpp"
#include "boost/geometry/geometries/register/ring.hpp"


namespace Kratos
{
namespace MPMSearchElementUtility
{

    // Standard types
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef Node<3> NodeType;
    typedef typename ModelPart::GeometryType GeometryType;

    // Boost Polygon
    typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> Boost2DPointType;
    typedef boost::geometry::model::polygon<Boost2DPointType> Boost2DPolygonType;

    // Container types
    typedef typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::IntegrationPointsContainerType IntegrationPointsContainerType;
    typedef typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;
    typedef typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;


    inline double CrossProductDet2D(array_1d<double, 3> VectorA, array_1d<double, 3> VectorB)
    {
        return (VectorA[0] * VectorB[1] - VectorB[0] * VectorA[1]);
    }

    inline bool CheckIsInside(const GeometryType& rGeom, array_1d<double, 3>& LocalCoords, const array_1d<double, 3>& Coords, const double Tolerance, const bool IsCalcLocalCoords = true)
    {
        bool is_inside = true;
        if (rGeom.Dimension() == 2)
        {
            is_inside = true;
            // Do walk around method
            Vector cross_products(rGeom.PointsNumber());
            for (size_t i = 0; i < rGeom.PointsNumber(); ++i)
            {
                if (rGeom.Points()[i].Coordinates()[2] != 0.0) {
                    return rGeom.IsInside(Coords, LocalCoords, Tolerance);
                    break;
                }
                cross_products[i] = CrossProductDet2D(Coords - rGeom.Points()[i].Coordinates(),
                    rGeom.Points()[(i + 1) % rGeom.PointsNumber()].Coordinates() - rGeom.Points()[i].Coordinates());
            }
            for (size_t i = 1; i < cross_products.size(); ++i)
            {
                if (cross_products[i] * cross_products[0] < -std::abs(Tolerance))
                {
                    is_inside = false;
                    break;
                }
            }

        }
        else return rGeom.IsInside(Coords, LocalCoords, Tolerance);

        if (is_inside) {
            if (IsCalcLocalCoords) return rGeom.IsInside(Coords, LocalCoords, Tolerance);
            else return true;
        }


        return false;
    }

    inline void ConstructNeighbourRelations(GeometryType& rGeom, const ModelPart& rBackgroundGridModelPart)
    {
        std::vector<typename Geometry<Node<3>>::Pointer> geometry_neighbours;
        for (IndexType j = 0; j < rBackgroundGridModelPart.NumberOfElements(); j++)
        {
            auto p_geometry_neighbour = (rBackgroundGridModelPart.ElementsBegin() + j)->pGetGeometry();
            if (p_geometry_neighbour->Id() != rGeom.Id()) // dont add the parent as its own neighbour
            {
                for (IndexType n = 0; n < p_geometry_neighbour->size(); n++)
                {
                    for (IndexType k = 0; k < rGeom.size(); k++)
                    {
                        if (rGeom[k].Id() == (*p_geometry_neighbour)[n].Id()) {
                            // Prevent duplicate additions
                            bool add_entry = true;
                            for (size_t i = 0; i < geometry_neighbours.size(); i++)
                            {
                                if (geometry_neighbours[i]->Id() == p_geometry_neighbour->Id())
                                {
                                    add_entry = false;
                                    break;
                                }
                            }
                            if (add_entry)
                            {
                                geometry_neighbours.push_back(p_geometry_neighbour);
                            }
                            break;
                        }
                    }
                }
            }

        }
        #pragma omp critical
        rGeom.SetValue(GEOMETRY_NEIGHBOURS, geometry_neighbours);
    }

    inline void CreateBoundingBoxPoints(std::vector<array_1d<double, 3>>& rPointVector,
        const array_1d<double, 3>& rCenter, const double SideHalfLength, const SizeType WorkingDim)
    {
        KRATOS_TRY

            if (WorkingDim == 2)
            {
                if (rPointVector.size() != 4) rPointVector.resize(4);
                for (size_t i = 0; i < 4; ++i) {
                    rPointVector[i].clear();
                    rPointVector[i] += rCenter;
                }
                rPointVector[0][0] -= SideHalfLength;
                rPointVector[1][0] += SideHalfLength;
                rPointVector[2][0] += SideHalfLength;
                rPointVector[3][0] -= SideHalfLength;

                rPointVector[0][1] -= SideHalfLength;
                rPointVector[1][1] -= SideHalfLength;
                rPointVector[2][1] += SideHalfLength;
                rPointVector[3][1] += SideHalfLength;
            }
            else
            {
                if (rPointVector.size() != 8) rPointVector.resize(8);
                for (size_t i = 0; i < 8; ++i) {
                    rPointVector[i].clear();
                    rPointVector[i] += rCenter;
                }
                rPointVector[0][0] -= SideHalfLength;
                rPointVector[1][0] += SideHalfLength;
                rPointVector[2][0] += SideHalfLength;
                rPointVector[3][0] -= SideHalfLength;
                rPointVector[4][0] -= SideHalfLength;
                rPointVector[5][0] += SideHalfLength;
                rPointVector[6][0] += SideHalfLength;
                rPointVector[7][0] -= SideHalfLength;

                rPointVector[0][1] -= SideHalfLength;
                rPointVector[1][1] -= SideHalfLength;
                rPointVector[2][1] += SideHalfLength;
                rPointVector[3][1] += SideHalfLength;
                rPointVector[4][1] -= SideHalfLength;
                rPointVector[5][1] -= SideHalfLength;
                rPointVector[6][1] += SideHalfLength;
                rPointVector[7][1] += SideHalfLength;

                rPointVector[0][2] -= SideHalfLength;
                rPointVector[1][2] -= SideHalfLength;
                rPointVector[2][2] -= SideHalfLength;
                rPointVector[3][2] -= SideHalfLength;
                rPointVector[4][2] += SideHalfLength;
                rPointVector[5][2] += SideHalfLength;
                rPointVector[6][2] += SideHalfLength;
                rPointVector[7][2] += SideHalfLength;
            }

        KRATOS_CATCH("")
    }

    inline void CheckPQMPM(IntegrationPointsArrayType& rIntergrationSubPoints, const double Tolerance, const Matrix& rN, const DenseVector<Matrix>& rDN_De)
    {
        KRATOS_TRY

            if (rIntergrationSubPoints.size() != rN.size1()) {
                KRATOS_INFO("MPMSearchElementUtility::Check - ") << "Shape function rows must equal number of sub-points!";
                KRATOS_ERROR << "ERROR";
            }

        for (size_t i = 0; i < rIntergrationSubPoints.size(); ++i) {
            if (rIntergrationSubPoints[i].Weight() < Tolerance) {
                KRATOS_INFO("MPMSearchElementUtility::Check - ") << "Volume fraction of sub-points is too small!";
                KRATOS_ERROR << "ERROR";
            }

            if (rIntergrationSubPoints[i].Weight() > 1.0) {
                KRATOS_INFO("MPMSearchElementUtility::Check - ") << "Volume fraction of sub-points is too large!";
                KRATOS_ERROR << "ERROR";
            }

        }

        for (size_t j = 0; j < rN.size2(); ++j) {
            SizeType nonzero_entries = 0;
            for (size_t i = 0; i < rIntergrationSubPoints.size(); i++) if (rN(i, j) > 0.0) nonzero_entries += 1;
            if (nonzero_entries != 1) {
                KRATOS_INFO("MPMSearchElementUtility::Check - ") << "There must be only one nonzero entry per shape function column!"
                    << "\nrN = " << rN;
                KRATOS_ERROR << "ERROR";
            }
        }

        KRATOS_CATCH("")
    }

    inline bool CheckAllPointsAreInGeom(
        const std::vector<array_1d<double, 3>>& rPoints,
        const GeometryType& rReferenceGeom,
        const double Tolerance)
    {
        KRATOS_TRY

        array_1d<double, 3> dummy_local_coords;
        for (size_t i = 0; i < rPoints.size(); ++i) {
            if (!CheckIsInside(rReferenceGeom, dummy_local_coords, rPoints[i], Tolerance, false)) {
                return false;
            }
        }
        return true;

        KRATOS_CATCH("")
    }

    inline void Check3DBackGroundMeshIsCubicAxisAligned(const std::vector<GeometryType*>& rIntersectedGeometries)
    {
        KRATOS_TRY

            NodeType point_low, point_high;
        for (size_t i = 0; i < rIntersectedGeometries.size(); ++i) {
            if (rIntersectedGeometries[i]->GetGeometryType() != GeometryData::Kratos_Hexahedra3D8) {
                #pragma omp single
                KRATOS_ERROR << "MPMSearchElementUtility::Check3DBackGroundMeshIsCubicAxisAligned - "
                    << "3D PQMPM CAN ONLY BE USED FOR AXIS-ALIGNED RECTANGULAR-PRISM BACKGROUND GRIDS" << std::endl;
            }
            rIntersectedGeometries[i]->BoundingBox(point_low, point_high);
            for (size_t j = 0; j < rIntersectedGeometries[i]->PointsNumber(); ++j) {
                for (size_t k = 0; k < 3; ++k) {
                    if (rIntersectedGeometries[i]->GetPoint(j).Coordinates()[k] != point_low[k]) {
                        if (rIntersectedGeometries[i]->GetPoint(j).Coordinates()[k] != point_high[k]) {
                            #pragma omp single
                            KRATOS_ERROR << "MPMSearchElementUtility::Check3DBackGroundMeshIsCubicAxisAligned - "
                                << "3D PQMPM CAN ONLY BE USED FOR AXIS-ALIGNED RECTANGULAR-PRISM BACKGROUND GRIDS" << std::endl;
                        }
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }


    inline Boost2DPolygonType Create2DPolygonBoundingSquareFromPointsFast(const std::vector<array_1d<double, 3>>& rPoints,
        const bool XActive = true, const bool YActive = true, const bool ZActive = false)
    {
        KRATOS_TRY

            Boost2DPolygonType rPolygon;
        std::vector<Boost2DPointType> rPolygonPoints(5);

        if (!XActive || !YActive || ZActive)  if (rPoints.size() != 8) {
            KRATOS_INFO("MPMSearchElementUtility::Create2DPolygonBoundingSquareFromPointsFast - ")
                << "ALL BOUNDING SQUARES SHOULD BE CONSTRUCTED IN XY SPACE EXCEPT FOR HEX BACKGROUND GRID\n";
            KRATOS_ERROR << "ERROR";
        }

        if (XActive && YActive && !ZActive)
        {
            for (size_t i = 0; i < 4; ++i) {
                rPolygonPoints[i] = Boost2DPointType(rPoints[i][0], rPoints[i][1]);
            }
        }
        else if (!XActive && YActive && ZActive) // 3D case only!
        {
            rPolygonPoints[0] = Boost2DPointType(rPoints[0][1], rPoints[0][2]);
            rPolygonPoints[1] = Boost2DPointType(rPoints[4][1], rPoints[4][2]);
            rPolygonPoints[2] = Boost2DPointType(rPoints[7][1], rPoints[7][2]);
            rPolygonPoints[3] = Boost2DPointType(rPoints[3][1], rPoints[3][2]); // as per Hexahedra3D8 node ordering
        }
        else if (XActive && !YActive && ZActive)
        {
            rPolygonPoints[0] = Boost2DPointType(rPoints[0][0], rPoints[0][2]);
            rPolygonPoints[1] = Boost2DPointType(rPoints[1][0], rPoints[1][2]);
            rPolygonPoints[2] = Boost2DPointType(rPoints[5][0], rPoints[5][2]);
            rPolygonPoints[3] = Boost2DPointType(rPoints[4][0], rPoints[4][2]);
        }
        else
        {
            KRATOS_INFO("MPMSearchElementUtility::Create2DPolygonBoundingSquareFromPointsFast - ")
                << "INVALID PLANE TO MAKE 2D POLYGON IN\n";
            KRATOS_ERROR << "ERROR";
        }
        rPolygonPoints[4] = rPolygonPoints[0];
        rPolygon.outer().assign(rPolygonPoints.begin(), rPolygonPoints.end());
        boost::geometry::correct(rPolygon); // to close the polygon

        return rPolygon;

        KRATOS_CATCH("")
    }

    inline Boost2DPolygonType Create2DPolygonFromGeometryFast(const GeometryType& rGeom,
        const bool XActive = true, const bool YActive = true, const bool ZActive = false)
    {
        KRATOS_TRY

            Boost2DPolygonType rPolygon;

        if (rGeom.WorkingSpaceDimension() == 3)
        {
            std::vector<Boost2DPointType> rPolygonPoints(5);
            NodeType point_low, point_high;
            rGeom.BoundingBox(point_low, point_high);

            if (XActive && YActive && !ZActive)
            {
                rPolygonPoints[0] = Boost2DPointType(point_low[0], point_low[1]);
                rPolygonPoints[1] = Boost2DPointType(point_high[0], point_low[1]);
                rPolygonPoints[2] = Boost2DPointType(point_high[0], point_high[1]);
                rPolygonPoints[3] = Boost2DPointType(point_low[0], point_high[1]);
            }
            else if (XActive && !YActive && ZActive)
            {
                rPolygonPoints[0] = Boost2DPointType(point_low[0], point_low[2]);
                rPolygonPoints[1] = Boost2DPointType(point_high[0], point_low[2]);
                rPolygonPoints[2] = Boost2DPointType(point_high[0], point_high[2]);
                rPolygonPoints[3] = Boost2DPointType(point_low[0], point_high[2]);
            }
            else if (!XActive && YActive && ZActive)
            {
                rPolygonPoints[0] = Boost2DPointType(point_low[1], point_low[2]);
                rPolygonPoints[1] = Boost2DPointType(point_high[1], point_low[2]);
                rPolygonPoints[2] = Boost2DPointType(point_high[1], point_high[2]);
                rPolygonPoints[3] = Boost2DPointType(point_low[1], point_high[2]);
            }
            else
            {
                KRATOS_INFO("MPMSearchElementUtility::Create2DPolygonFromGeometryFast - ")
                    << "INVALID PLANE TO MAKE 2D POLYGON IN\n";
                KRATOS_ERROR << "ERROR";
            }
            rPolygonPoints[4] = rPolygonPoints[0];
            rPolygon.outer().assign(rPolygonPoints.begin(), rPolygonPoints.end());
        }
        else
        {
            std::vector<Boost2DPointType> rPolygonPoints(rGeom.PointsNumber() + 1);
            for (size_t i = 0; i < rGeom.PointsNumber(); ++i) {
                rPolygonPoints[i] = Boost2DPointType(rGeom.GetPoint(i).X(), rGeom.GetPoint(i).Y());
            }
            rPolygonPoints[rGeom.PointsNumber()] = rPolygonPoints[0];
            rPolygon.outer().assign(rPolygonPoints.begin(), rPolygonPoints.end());
        }
        boost::geometry::correct(rPolygon); // to close the polygon
        return rPolygon;

        KRATOS_CATCH("")
    }

    inline IntegrationPoint<3> CreateSubPoint(const array_1d<double, 3>& rGlobalCoords, const double rVolumeFraction,
        const GeometryType& rBackgroundGridElementGeom, Vector& rN, Matrix& rDN_De)
    {
        KRATOS_TRY

            array_1d<double, 3> local_coordinates;
        rBackgroundGridElementGeom.PointLocalCoordinates(local_coordinates, rGlobalCoords);
        rBackgroundGridElementGeom.ShapeFunctionsValues(rN, local_coordinates);
        rBackgroundGridElementGeom.ShapeFunctionsLocalGradients(rDN_De, local_coordinates);

        return IntegrationPoint<3>(local_coordinates, rVolumeFraction);

        KRATOS_CATCH("")
    }

    inline void Determine2DSubPoint(const GeometryType& rGridElement, const std::vector<array_1d<double, 3>>& rMasterDomainPoints,
        array_1d<double, 3>& rSubPointCoord, double& rSubPointVolume)
    {
        KRATOS_TRY

            // make boost polygon of current background element geometry
            Boost2DPolygonType polygon_grid = Create2DPolygonFromGeometryFast(rGridElement);

        // make boost polygon of bounding box
        Boost2DPolygonType polygon_box = Create2DPolygonBoundingSquareFromPointsFast(rMasterDomainPoints);


        // make boost polygon result container
        std::vector<Boost2DPolygonType> polygon_result_container;

        // reset accumulated quantities
        rSubPointVolume = 0.0;
        rSubPointCoord.clear();
        Boost2DPointType centroid_result;

        // accumulate result over intersected sub-polygons
        if (boost::geometry::intersection(polygon_grid, polygon_box, polygon_result_container)) {
            for (auto& polygon_result : polygon_result_container) {
                rSubPointVolume += boost::geometry::area(polygon_result);
                boost::geometry::centroid(polygon_result, centroid_result);
                rSubPointCoord[0] += centroid_result.get<0>();
                rSubPointCoord[1] += centroid_result.get<1>();
            }
        }
        else
        {
            KRATOS_INFO("MPMSearchElementUtility::Determine2DSubPoint - ")
                << "BOOST INTERSECTION FAILED ALTHOUGH KRATOS INTERSECTION WORKED\n";
            KRATOS_ERROR << "ERROR";
        }

        rSubPointCoord /= double(polygon_result_container.size());

        KRATOS_CATCH("")
    }

    inline void Determine3DSubPoint(const GeometryType& rGridElement, const std::vector<array_1d<double, 3>>& rMasterDomainPoints,
        array_1d<double, 3>& rSubPointCoord, double& rSubPointVolume)
    {
        KRATOS_TRY

            // NOTE: THIS FUNCTION ASSUMES THE BACKGROUND GRID ELEMENT IS PERFECTLY CUBIC
            //    AND THE RESULTING INTERSECTION VOLUME IS A RECTANGULAR PRISM

            // make boost xy polygon of current background element geometry
            Boost2DPolygonType polygon_grid_xy = Create2DPolygonFromGeometryFast(rGridElement);

        // make boost yz polygon of current background element geometry
        Boost2DPolygonType polygon_grid_yz = Create2DPolygonFromGeometryFast(rGridElement, false, true, true);

        // make boost xy polygon of bounding box
        Boost2DPolygonType polygon_box_xy = Create2DPolygonBoundingSquareFromPointsFast(rMasterDomainPoints);

        // make boost yz polygon of bounding box
        Boost2DPolygonType polygon_box_yz = Create2DPolygonBoundingSquareFromPointsFast(rMasterDomainPoints, false, true, true);

        // make boost polygon result container
        std::vector<Boost2DPolygonType> polygon_xy_result_container;
        std::vector<Boost2DPolygonType> polygon_yz_result_container;

        // reset accumulated quantities
        double sub_volume_area = 0.0;

        rSubPointCoord.clear();
        Boost2DPointType centroid_result;

        // Determine area and x y coordinates from xy polygons
        if (boost::geometry::intersection(polygon_grid_xy, polygon_box_xy, polygon_xy_result_container)) {
            for (auto& polygon_result : polygon_xy_result_container) {
                sub_volume_area += boost::geometry::area(polygon_result);
                boost::geometry::centroid(polygon_result, centroid_result);
                rSubPointCoord[0] += centroid_result.get<0>();
                rSubPointCoord[1] += centroid_result.get<1>();
            }
        }
        else
        {
            KRATOS_INFO("MPMSearchElementUtility::Determine3DSubPoint - ")
                << "BOOST INTERSECTION FAILED ALTHOUGH KRATOS INTERSECTION WORKED\n";
            KRATOS_ERROR << "ERROR";
        }

        rSubPointCoord /= double(polygon_xy_result_container.size()); // at the moment this is just the xy coords!

        // Perform yz polygon intersection to determine depth and z-position of sub-point
        // local x = global y
        // local y = global z
        array_1d<double, 2> sub_point_z_coord = ZeroVector(2);
        bool is_initialized = false;
        double min_z = 0.0;
        double max_z = 0.0;
        if (boost::geometry::intersection(polygon_grid_yz, polygon_box_yz, polygon_yz_result_container)) {
            for (auto& polygon_result : polygon_yz_result_container) {
                for (auto& result_point : polygon_result.outer()) {
                    if (!is_initialized) {
                        min_z = result_point.get<1>();
                        max_z = result_point.get<1>();
                        is_initialized = true;
                    }
                    else if (result_point.get<1>() < min_z) min_z = result_point.get<1>();
                    else if (result_point.get<1>() > max_z) max_z = result_point.get<1>();
                }
            }
        }
        else
        {
            KRATOS_INFO("MPMSearchElementUtility::Determine3DSubPoint - ")
                << "BOOST INTERSECTION FAILED ALTHOUGH KRATOS INTERSECTION WORKED\n";
            KRATOS_ERROR << "ERROR";
        }

        rSubPointCoord[2] = 0.5 * (min_z + max_z);
        rSubPointVolume = sub_volume_area * (max_z - min_z);

        KRATOS_CATCH("")
    }

    inline typename Geometry<Node<3>>::Pointer CreateCustomQuadraturePoint(
        SizeType WorkingSpaceDimension,
        SizeType LocalSpaceDimension,
        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rShapeFunctionContainer,
        typename Geometry<Node<3>>::PointsArrayType rPoints,
        GeometryType* pGeometryParent)
    {
        KRATOS_TRY

            if (WorkingSpaceDimension == 1 && LocalSpaceDimension == 1)
                return Kratos::make_shared<
                QuadraturePointPartitionedGeometry<Node<3>, 1>>(
                    rPoints, rShapeFunctionContainer, pGeometryParent);
            else if (WorkingSpaceDimension == 2 && LocalSpaceDimension == 1)
                return Kratos::make_shared<
                QuadraturePointPartitionedGeometry<Node<3>, 2, 1>>(
                    rPoints, rShapeFunctionContainer, pGeometryParent);
            else if (WorkingSpaceDimension == 2 && LocalSpaceDimension == 2)
                return Kratos::make_shared<
                QuadraturePointPartitionedGeometry<Node<3>, 2>>(
                    rPoints, rShapeFunctionContainer, pGeometryParent);
            else if (WorkingSpaceDimension == 3 && LocalSpaceDimension == 2)
                return Kratos::make_shared<
                QuadraturePointPartitionedGeometry<Node<3>, 3, 2>>(
                    rPoints, rShapeFunctionContainer, pGeometryParent);
            else if (WorkingSpaceDimension == 3 && LocalSpaceDimension == 3)
                return Kratos::make_shared<
                QuadraturePointPartitionedGeometry<Node<3>, 3>>(
                    rPoints, rShapeFunctionContainer, pGeometryParent);
            else {
                KRATOS_ERROR << "Working/Local space dimension combinations are "
                    << "not provided for QuadraturePointGeometry. WorkingSpaceDimension: "
                    << WorkingSpaceDimension << ", LocalSpaceDimension: " << LocalSpaceDimension
                    << std::endl;
            }
        KRATOS_CATCH("")
    }

    inline bool DetermineIfDomainOverlapsBoundaryConditions(std::vector<GeometryType*>& IntersectedGeometries,
        const array_1d<double, 3>& rCoordinates, const double SideHalfLength)
    {
        const double reach = (IntersectedGeometries[0]->WorkingSpaceDimension() == 3)
            ? 1.7321 * SideHalfLength : 1.414214 * SideHalfLength;
        for (size_t i = 0; i < IntersectedGeometries.size(); ++i)
        {
            for (size_t j = 0; j < IntersectedGeometries[i]->PointsNumber(); ++j) {
                auto node_it = IntersectedGeometries[i]->pGetPoint(j);
                bool is_fixed = false;
                if (node_it->IsFixed(DISPLACEMENT_X)) is_fixed = true;
                else if (node_it->IsFixed(DISPLACEMENT_Y)) is_fixed = true;
                else if (node_it->HasDofFor(DISPLACEMENT_Z))  if (node_it->IsFixed(DISPLACEMENT_Z)) is_fixed = true;
                if (is_fixed) {
                    const double fixed_point_to_cog = norm_2(node_it->Coordinates() - rCoordinates);
                    if (fixed_point_to_cog <= reach) return true;
                }
            }
        }

        return false;
    }

    inline bool IntersectionCheckWithBoundingBox(const GeometryType& rGeom, const array_1d<double, 3>& rCoord, const double SideHalfLength)
    {
        const double z_coord = (rGeom.WorkingSpaceDimension() == 3) ? SideHalfLength : 0.0;
        const Point point_low(rCoord[0] - SideHalfLength, rCoord[1] - SideHalfLength, rCoord[2] - z_coord);
        const Point point_high(rCoord[0] + SideHalfLength, rCoord[1] + SideHalfLength, rCoord[2] + z_coord);
        NodeType ele_point_low, ele_point_high;

        const double dimension_45_degree_factor = (rGeom.WorkingSpaceDimension() == 3) ? 1.7321 : 1.414214;
        double center_to_center = norm_2(rGeom.Center() - rCoord);
        rGeom.BoundingBox(ele_point_low, ele_point_high);
        double maximum_contact_range = dimension_45_degree_factor * SideHalfLength +
            norm_2(ele_point_high - ele_point_low);
        if (center_to_center <= maximum_contact_range) return true;
        return false;
    }


    inline void RecursivePQMPMNeighbourSearch(const ModelPart& rBackgroundGridModelPart,
        std::vector<GeometryType*>& rIntersectedGeometries,
        const Point& rPointLow,
        const Point& rPointHigh,
        IndexType& RecursionCount,
        const array_1d<double, 3>& rCoordinates,
        const double SideHalfLength,
        const SizeType MaxRecursions = 100)
    {
        RecursionCount += 1;
        if (RecursionCount < MaxRecursions)
        {
            if (!rIntersectedGeometries.back()->Has(GEOMETRY_NEIGHBOURS))
                ConstructNeighbourRelations(*rIntersectedGeometries.back(), rBackgroundGridModelPart);
            auto& geometry_neighbours = rIntersectedGeometries.back()->GetValue(GEOMETRY_NEIGHBOURS);
            bool check_geom;
            for (size_t i = 0; i < geometry_neighbours.size(); ++i) {
                // dont check elements we have already intersected with
                check_geom = true;
                for (size_t j = 0; j < rIntersectedGeometries.size(); ++j) {
                    if (geometry_neighbours[i]->Id() == rIntersectedGeometries[j]->Id()) {
                        check_geom = false;
                        break;
                    }
                }

                if (check_geom) {
                    // check if this background grid and the MP domain overlap
                    if (IntersectionCheckWithBoundingBox(*geometry_neighbours[i], rCoordinates, SideHalfLength)) {
                        if (geometry_neighbours[i]->HasIntersection(rPointLow, rPointHigh)) {
                            // add to container and then search its neighbours
                            rIntersectedGeometries.push_back(geometry_neighbours[i].get());

                                RecursivePQMPMNeighbourSearch(rBackgroundGridModelPart,
                                    rIntersectedGeometries, rPointLow, rPointHigh,
                                    RecursionCount,
                                    rCoordinates, SideHalfLength,
                                    MaxRecursions);
                        }
                    }
                }
            }
        }
        else KRATOS_INFO("RecursivePQMPMNeighbourSearch:: ") <<
                "Recursion count of " << MaxRecursions << " exceeded\n" << std::endl;
    }


    inline bool CheckFixedNodesWithinBoundingBox(const PointerVector<Node<3>>& rNodesList,
        const Point& rPointHigh, const Point& rPointLow, const SizeType WorkingDim)
    {
        for (auto& node_it : rNodesList)
        {
            bool is_fixed = false;

            if (node_it.IsFixed(DISPLACEMENT_X)) is_fixed = true;
            else if (node_it.IsFixed(DISPLACEMENT_Y)) is_fixed = true;
            else if (WorkingDim == 3) if (node_it.IsFixed(DISPLACEMENT_Z)) is_fixed = true;

            if (is_fixed)
            {
                bool is_inside = true;
                for (size_t i = 0; i < WorkingDim; i++)
                {
                    if (rPointLow.Coordinates()[i] <= node_it.Coordinates()[i] &&
                        node_it.Coordinates()[i] <= rPointHigh.Coordinates()[i])
                    {
                        // we are inside for this dimension
                    }
                    else
                    {
                        is_inside = false;
                        break;
                    }
                }
                if (is_inside) return true;
            }
        }

        return false;
    }


    inline void PartitionMasterMaterialPointsIntoSubPoints(const ModelPart& rBackgroundGridModelPart,
        const array_1d<double, 3>& rCoordinates,
        const array_1d<double, 3>& rLocalCoords,
        Element& rMasterMaterialPoint,
        typename GeometryType::Pointer pQuadraturePointGeometry,
        const double Tolerance)
    {
        KRATOS_TRY;

        GeometryType& rParentGeom = pQuadraturePointGeometry->GetGeometryParent(0);

        // If axisymmetric make normal MP
        if (rBackgroundGridModelPart.GetProcessInfo().Has(IS_AXISYMMETRIC)) {
            if (rBackgroundGridModelPart.GetProcessInfo().GetValue(IS_AXISYMMETRIC)) {
                CreateQuadraturePointsUtility<Node<3>>::UpdateFromLocalCoordinates(
                    pQuadraturePointGeometry, rLocalCoords,
                    rMasterMaterialPoint.GetGeometry().IntegrationPoints()[0].Weight(),
                    rParentGeom);
                return;
            }
        }

        const SizeType working_dim = rParentGeom.WorkingSpaceDimension();
        const double pqmpm_min_fraction = (rBackgroundGridModelPart.GetProcessInfo().Has(PQMPM_SUBPOINT_MIN_VOLUME_FRACTION))
            ? std::max(rBackgroundGridModelPart.GetProcessInfo()[PQMPM_SUBPOINT_MIN_VOLUME_FRACTION], std::numeric_limits<double>::epsilon())
            : std::numeric_limits<double>::epsilon();

        // Get volume and set up master domain bounding points
        std::vector<double> mp_volume_vec;
        rMasterMaterialPoint.CalculateOnIntegrationPoints(MP_VOLUME, mp_volume_vec, rBackgroundGridModelPart.GetProcessInfo());
        if (rBackgroundGridModelPart.GetProcessInfo()[DOMAIN_SIZE] == 2 && rMasterMaterialPoint.GetProperties().Has(THICKNESS))
            mp_volume_vec[0] /= rMasterMaterialPoint.GetProperties()[THICKNESS];
        const double side_half_length = std::pow(mp_volume_vec[0], 1.0 / double(working_dim)) / 2.0;
        std::vector<array_1d<double, 3>> master_domain_points(std::pow(2.0, working_dim));
        CreateBoundingBoxPoints(master_domain_points, rCoordinates, side_half_length, working_dim);

        // Initially check if the bounding box volume scalar is less than the element volume scalar
        if (mp_volume_vec[0] <= rParentGeom.DomainSize()) {
            if (CheckAllPointsAreInGeom(master_domain_points, rParentGeom, Tolerance)) {
                CreateQuadraturePointsUtility<Node<3>>::UpdateFromLocalCoordinates(
                    pQuadraturePointGeometry, rLocalCoords,
                    rMasterMaterialPoint.GetGeometry().IntegrationPoints()[0].Weight(),
                    rParentGeom);
                return;
            }
        }

        // we need to do splitting. Initially determine all grid elements we intersect with
        Point point_low(rCoordinates[0] - side_half_length, rCoordinates[1] - side_half_length, rCoordinates[2]);
        if (working_dim == 3) point_low[2] -= side_half_length;
        Point point_high(rCoordinates[0] + side_half_length, rCoordinates[1] + side_half_length, rCoordinates[2]);
        if (working_dim == 3) point_high[2] += side_half_length;
        SizeType number_of_nodes = 0;
        std::vector<GeometryType*> intersected_geometries;

        // Do neighbour searching to determine the intersected geometries
        IndexType recursion_count = 0;
        intersected_geometries.push_back(&rParentGeom);
        RecursivePQMPMNeighbourSearch(rBackgroundGridModelPart,
            intersected_geometries, point_low , point_high , recursion_count,
            rCoordinates, side_half_length);
        for (size_t i = 0; i < intersected_geometries.size(); ++i)
            number_of_nodes += intersected_geometries[i]->PointsNumber();

        // If we are 3D, check background mesh are axis-aligned perfect rectangular prisms
        if (working_dim == 3) Check3DBackGroundMeshIsCubicAxisAligned(intersected_geometries);

        // Prepare containers to hold all sub-points
        const SizeType number_of_sub_material_points = intersected_geometries.size();
        PointerVector<Node<3>> nodes_list(number_of_nodes);
        IntegrationPointsArrayType ips(number_of_sub_material_points);
        Matrix N_matrix(number_of_sub_material_points, number_of_nodes, -1.0);
        DenseVector<Matrix> DN_De_vector(number_of_sub_material_points);

        // Temporary local containers
        double sub_point_volume;
        array_1d<double, 3> sub_point_position;
        IndexType active_node_index = 0;
        IndexType active_subpoint_index = 0;

        // Loop over all intersected grid elements and make subpoints in each
        for (size_t i = 0; i < number_of_sub_material_points; ++i) {
            Matrix DN_De(intersected_geometries[i]->PointsNumber(), working_dim);
            Vector N(intersected_geometries[i]->PointsNumber());
            sub_point_position.clear();
            sub_point_volume = 0.0;
            IntegrationPoint<3> trial_subpoint;

            if (working_dim == 2) {
                Determine2DSubPoint(*intersected_geometries[i], master_domain_points, sub_point_position, sub_point_volume);
                sub_point_position[2] = rCoordinates[2]; // set z coord of sub point to that of the master
            }
            else Determine3DSubPoint(*intersected_geometries[i], master_domain_points, sub_point_position, sub_point_volume);

            trial_subpoint = CreateSubPoint(sub_point_position, sub_point_volume / mp_volume_vec[0],
                *intersected_geometries[i], N, DN_De);

            // Transfer local data to containers
            if (trial_subpoint.Weight() > pqmpm_min_fraction) {
                ips[active_subpoint_index] = trial_subpoint;
                DN_De_vector[active_subpoint_index] = DN_De;
                for (size_t j = 0; j < N.size(); ++j) {
                    N_matrix(active_subpoint_index, active_node_index) = N[j];
                    nodes_list(active_node_index) = intersected_geometries[i]->pGetPoint(j);

                    active_node_index += 1;
                }
                active_subpoint_index += 1;
            }
        }
        if (active_subpoint_index == 1) {
            CreateQuadraturePointsUtility<Node<3>>::UpdateFromLocalCoordinates(
                pQuadraturePointGeometry, rLocalCoords, rMasterMaterialPoint.GetGeometry().IntegrationPoints()[0].Weight(),
                rParentGeom);
            return;
        }

        IntegrationPointsArrayType ips_active(active_subpoint_index);
        PointerVector<Node<3>> nodes_list_active(active_node_index);
        if (ips_active.size() == ips.size()) {
            ips_active = ips;
            nodes_list_active = nodes_list;
        }
        else {
            N_matrix.resize(active_subpoint_index, active_node_index, true);
            DN_De_vector.resize(active_subpoint_index, true);
            for (size_t i = 0; i < active_subpoint_index; ++i) ips_active[i] = ips[i];
            for (size_t i = 0; i < active_node_index; ++i) nodes_list_active(i) = nodes_list(i);
        }

        // check if there are any fixed nodes within the bounding box
        if (CheckFixedNodesWithinBoundingBox(nodes_list_active, point_high, point_low, working_dim))
        {
            CreateQuadraturePointsUtility<Node<3>>::UpdateFromLocalCoordinates(
                pQuadraturePointGeometry, rLocalCoords, rMasterMaterialPoint.GetGeometry().IntegrationPoints()[0].Weight(),
                rParentGeom);
            return;
        }

        // Check volume fractions sum to unity
        double vol_sum = 0.0;
        for (size_t i = 0; i < ips_active.size(); ++i) vol_sum += ips_active[i].Weight();
        if (std::abs(vol_sum - 1.0) > Tolerance) {
            const bool is_pqmpm_fallback = (rBackgroundGridModelPart.GetProcessInfo().Has(IS_MAKE_NORMAL_MP_IF_PQMPM_FAILS))
                ? rBackgroundGridModelPart.GetProcessInfo().GetValue(IS_MAKE_NORMAL_MP_IF_PQMPM_FAILS) : false;
            if (is_pqmpm_fallback) {
                CreateQuadraturePointsUtility<Node<3>>::UpdateFromLocalCoordinates(
                    pQuadraturePointGeometry, rLocalCoords, rMasterMaterialPoint.GetGeometry().IntegrationPoints()[0].Weight(),
                    rParentGeom);
                return;
            }
            else {
                #pragma omp critical
                KRATOS_INFO("MPMSearchElementUtility::Check")
                    << "Volume fraction of sub-points does not approximately sum to 1.0."
                    << " This probably means the background grid is not big enough."
                    << "\nPosition = " << rCoordinates
                    << "\nNumber of active sub points = " << ips_active.size()
                    << "\nNumber of trial sub points = " << ips.size()
                    << "\nMaterial point volume = " << mp_volume_vec[0]
                    << "\nTotal volume fraction = " << vol_sum << "\nIndividual volume fractions:\n";
                for (size_t i = 0; i < ips_active.size(); ++i) std::cout << "\t" << ips_active[i].Weight()
                    << "\t\t" << ips_active[i].Coordinates() << std::endl;
                KRATOS_ERROR << "ERROR";
            }
        } else CheckPQMPM(ips_active, std::numeric_limits<double>::epsilon(), N_matrix, DN_De_vector);

        // Transfer data over
        GeometryData::IntegrationMethod ThisDefaultMethod = pQuadraturePointGeometry->GetDefaultIntegrationMethod();
        IntegrationPointsContainerType ips_container;
        ips_container[ThisDefaultMethod] = ips_active;
        ShapeFunctionsValuesContainerType shape_function_container;
        shape_function_container[ThisDefaultMethod] = N_matrix;
        ShapeFunctionsLocalGradientsContainerType shape_function_derivatives_container;
        shape_function_derivatives_container[ThisDefaultMethod] = DN_De_vector;
        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container( ThisDefaultMethod,
            ips_container, shape_function_container, shape_function_derivatives_container);

        for (size_t i = 0; i < nodes_list_active.size(); ++i) nodes_list_active[i].Set(ACTIVE);
        if (pQuadraturePointGeometry->IntegrationPointsNumber() == 1) {
            pQuadraturePointGeometry = CreateCustomQuadraturePoint(working_dim, pQuadraturePointGeometry->LocalSpaceDimension(),
                data_container, nodes_list_active, &rParentGeom);
            rMasterMaterialPoint.SetGeometry(pQuadraturePointGeometry);
        }

        pQuadraturePointGeometry->SetGeometryShapeFunctionContainer(data_container);
        pQuadraturePointGeometry->Points() = nodes_list_active;

        KRATOS_CATCH("");
    }

    inline bool IsExplicitAndNeedsCorrection(GeometryType::Pointer pQuadraturePoint, const ProcessInfo& rProcessInfo)
    {
        if (rProcessInfo.Has(IS_FIX_EXPLICIT_MP_ON_GRID_EDGE)) {
            if (rProcessInfo.GetValue(IS_FIX_EXPLICIT_MP_ON_GRID_EDGE)) {
                if (pQuadraturePoint->IntegrationPointsNumber() == 1)
                {
                    for (size_t i = 0; i < pQuadraturePoint->ShapeFunctionsValues().size2(); ++i)
                    {
                        if (pQuadraturePoint->ShapeFunctionsValues()(0, i) < std::numeric_limits<double>::epsilon()) return true;
                    }
                }
            }
        }

        return false;
    }

    inline GeometryType& FindGridGeom(GeometryType& rParentGeom,
        const ModelPart& rBackgroundGridModelPart,
        const double Tolerance,
        const array_1d<double, 3>& xg,
        array_1d<double, 3>& rLocalCoords,
        const ProcessInfo& rProcessInfo,
        bool& IsFound)
    {
        IsFound = false;

        if (CheckIsInside(rParentGeom, rLocalCoords, xg, Tolerance)) {
            IsFound = true;
            return rParentGeom;
        }
        else
        {
            if (!rParentGeom.Has(GEOMETRY_NEIGHBOURS))
                ConstructNeighbourRelations(rParentGeom, rBackgroundGridModelPart);

            auto& geometry_neighbours = rParentGeom.GetValue(GEOMETRY_NEIGHBOURS);
            for (IndexType k = 0; k < geometry_neighbours.size(); ++k) {
                if (CheckIsInside(*geometry_neighbours[k], rLocalCoords, xg, Tolerance)) {
                    IsFound = true;
                    return *(geometry_neighbours[k].get());
                }
            }
        }

        return rParentGeom;
    }


    inline void UpdatePartitionedQuadraturePoint(const ModelPart& rBackgroundGridModelPart,
        const array_1d<double, 3>& rCoordinates,
        Element& rMasterMaterialPoint,
        typename GeometryType::Pointer pQuadraturePointGeometry,
        const double Tolerance)
    {
        array_1d<double, 3> local_coords;
        pQuadraturePointGeometry->IsInside(rCoordinates, local_coords, Tolerance);
        PartitionMasterMaterialPointsIntoSubPoints(rBackgroundGridModelPart, rCoordinates,
            local_coords, rMasterMaterialPoint, pQuadraturePointGeometry, Tolerance);
    }


    inline void NeighbourSearchElements(const ModelPart& rMPMModelPart,
        const ModelPart& rBackgroundGridModelPart,
        std::vector<typename Element::Pointer>& rMissingElements,
        const double Tolerance)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rMPMModelPart.Elements().size()); ++i) {
            auto element_itr = (rMPMModelPart.ElementsBegin() + i);
            array_1d<double, 3> local_coordinates;
            bool is_found = false;
            std::vector<array_1d<double, 3>> xg;
            element_itr->CalculateOnIntegrationPoints(MP_COORD, xg, rBackgroundGridModelPart.GetProcessInfo());

            GeometryType& r_found_geom = FindGridGeom(element_itr->GetGeometry().GetGeometryParent(0),
                rBackgroundGridModelPart, Tolerance, xg[0], local_coordinates,
                rMPMModelPart.GetProcessInfo(), is_found);

            if (is_found)
            {
                const bool is_pqmpm = (rBackgroundGridModelPart.GetProcessInfo().Has(IS_PQMPM))
                    ? rBackgroundGridModelPart.GetProcessInfo().GetValue(IS_PQMPM) : false;
                if (is_pqmpm)
                {
                    // Updates the quadrature point geometry.
                    PartitionMasterMaterialPointsIntoSubPoints(rBackgroundGridModelPart, xg[0],
                        local_coordinates, *element_itr, element_itr->pGetGeometry(), Tolerance);
                }
                else
                {
                    CreateQuadraturePointsUtility<Node<3>>::UpdateFromLocalCoordinates(
                        element_itr->pGetGeometry(), local_coordinates,
                        element_itr->GetGeometry().IntegrationPoints()[0].Weight(), r_found_geom);
                }

                if (IsExplicitAndNeedsCorrection(element_itr->pGetGeometry(), rBackgroundGridModelPart.GetProcessInfo()))
                    is_found = false;
                else {
                    for (IndexType j = 0; j < r_found_geom.PointsNumber(); ++j)
                        r_found_geom.Points()[j].Set(ACTIVE);
                }
            }
            if(!is_found)
            {
                #pragma omp critical
                rMissingElements.push_back(&*element_itr);
            }
        }
    }


    //


    inline void NeighbourSearchConditions(const ModelPart& rMPMModelPart,
        const ModelPart& rBackgroundGridModelPart,
        std::vector<typename Condition::Pointer>& rMissingConditions,
        const double Tolerance)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rMPMModelPart.Conditions().size()); ++i) {
            auto condition_itr = rMPMModelPart.Conditions().begin() + i;

            std::vector<array_1d<double, 3>> xg;
            condition_itr->CalculateOnIntegrationPoints(MPC_COORD, xg, rMPMModelPart.GetProcessInfo());

            if (xg.size() > 0 && condition_itr->Is(BOUNDARY))
            {
                array_1d<double, 3> local_coordinates;
                bool is_found = false;

                GeometryType& r_found_geom = FindGridGeom(condition_itr->GetGeometry(),
                    rBackgroundGridModelPart, Tolerance, xg[0], local_coordinates,
                    rMPMModelPart.GetProcessInfo(), is_found);

                if (is_found)
                {
                    condition_itr->GetGeometry() = r_found_geom;
                    for (IndexType j = 0; j < r_found_geom.PointsNumber(); ++j)
                        r_found_geom[j].Set(ACTIVE);
                }
                else
                {
                    #pragma omp critical
                    rMissingConditions.push_back(&*condition_itr);
                }
            }
        }
    }


    inline bool IsFixExplicitAndOnElementEdge(const Vector& N, const ProcessInfo& rProcessInfo)
    {
        if (rProcessInfo.Has(IS_FIX_EXPLICIT_MP_ON_GRID_EDGE)) {
            if (rProcessInfo.GetValue(IS_FIX_EXPLICIT_MP_ON_GRID_EDGE)) {
                // check if MP is exactly on the edge of the element, this gives spurious strains in explicit
                for (SizeType i = 0; i < N.size(); ++i) {
                    if (std::abs(N[i]) < std::numeric_limits<double>::epsilon()) {
                        return true;
                    }
                }
            }
        }

        return false;
    }


    template <std::size_t TDimension>
    void BinBasedSearchElementsAndConditions(ModelPart& rMPMModelPart,
        ModelPart& rBackgroundGridModelPart,
        std::vector<typename Element::Pointer>& rMissingElements,
        std::vector<typename Condition::Pointer>& rMissingConditions,
        const std::size_t MaxNumberOfResults, const double Tolerance)
    {
        const ProcessInfo& r_process_info = rBackgroundGridModelPart.GetProcessInfo();
        bool is_pqmpm = (r_process_info.Has(IS_PQMPM))
            ? r_process_info.GetValue(IS_PQMPM) : false;

        // Search background grid and make element active
        Vector N;
        const int max_result = 1000;

        #pragma omp parallel
        {
            BinBasedFastPointLocator<TDimension> SearchStructure(rBackgroundGridModelPart);
            SearchStructure.UpdateSearchDatabase();
            typename BinBasedFastPointLocator<TDimension>::ResultContainerType results(max_result);

            // Element search and assign background grid
            #pragma omp for
            for (int i = 0; i < static_cast<int>(rMissingElements.size()); ++i) {
                auto element_itr = *(rMissingElements.begin() + i);
                std::vector<array_1d<double, 3>> xg;
                element_itr->CalculateOnIntegrationPoints(MP_COORD, xg, rMPMModelPart.GetProcessInfo());
                typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();
                Element::Pointer pelem;

                // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
                bool is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);

                if (is_found == true) {
                    if (IsFixExplicitAndOnElementEdge(N, r_process_info) && !is_pqmpm) {
                        // MP is exactly on the edge. Now we give it a little 'nudge'
                        array_1d<double, 3> xg_nudged = array_1d<double, 3>(xg[0]);
                        std::vector<array_1d<double, 3>> mp_vel;
                        element_itr->CalculateOnIntegrationPoints(MP_VELOCITY, mp_vel, rMPMModelPart.GetProcessInfo());
                        xg_nudged += r_process_info[DELTA_TIME] / 1000.0 * mp_vel[0];
                        if (SearchStructure.FindPointOnMesh(xg_nudged, N, pelem, result_begin, MaxNumberOfResults, Tolerance)) {
                            element_itr->SetValuesOnIntegrationPoints(MP_COORD, { xg_nudged }, rMPMModelPart.GetProcessInfo());
                            KRATOS_INFO("MPMSearchElementUtility") << "WARNING: To prevent spurious explicit stresses, Material Point "
                                << element_itr->Id() << " was nudged." << std::endl;
                        } else {
                            is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);
                            KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Material Point " << element_itr->Id()
                                << " lies exactly on an element edge and may give spurious results." << std::endl;
                        }
                    }
                    pelem->Set(ACTIVE);

                    const bool is_pqmpm = (rBackgroundGridModelPart.GetProcessInfo().Has(IS_PQMPM))
                        ? rBackgroundGridModelPart.GetProcessInfo().GetValue(IS_PQMPM) : false;
                    if (is_pqmpm)
                    {
                        // Updates the quadrature point geometry.
                        UpdatePartitionedQuadraturePoint(rBackgroundGridModelPart, xg[0],
                            *element_itr, pelem->pGetGeometry(), Tolerance);
                    }
                    else
                    {
                        auto p_quadrature_point_geometry = element_itr->pGetGeometry();
                        array_1d<double, 3> local_coordinates;
                        p_quadrature_point_geometry->PointLocalCoordinates(local_coordinates, xg[0]);
                        CreateQuadraturePointsUtility<Node<3>>::UpdateFromLocalCoordinates(
                            p_quadrature_point_geometry, local_coordinates,
                            p_quadrature_point_geometry->IntegrationPoints()[0].Weight(), pelem->GetGeometry());
                    }

                    auto& r_geometry = element_itr->GetGeometry();
                    for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j)
                        r_geometry[j].Set(ACTIVE);
                }
                else {
                    KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Search Element for Material Point: " << element_itr->Id()
                        << " is failed. Geometry is cleared." << std::endl;

                    element_itr->GetGeometry().clear();
                    element_itr->Reset(ACTIVE);
                    element_itr->Set(TO_ERASE);
                }
            }

            // Condition search and assign background grid
            #pragma omp for
            for (int i = 0; i < static_cast<int>(rMissingConditions.size()); ++i) {
                auto condition_itr = *(rMissingConditions.begin() + i);
                std::vector<array_1d<double, 3>> xg;
                condition_itr->CalculateOnIntegrationPoints(MPC_COORD, xg, rMPMModelPart.GetProcessInfo());

                if (xg.size() > 0) {
                    // Only search for particle based BCs!
                    // Grid BCs are still applied on MP_model_part but we don't want to search for them.
                    typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();
                    Element::Pointer pelem;

                    // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
                    bool is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);

                    if (is_found == true) {
                        condition_itr->GetGeometry() = pelem->GetGeometry();
                        auto& r_geometry = condition_itr->GetGeometry();

                        for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j)
                            r_geometry[j].Set(ACTIVE);
                    } else {
                        KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Search Element for Material Point Condition: " << condition_itr->Id()
                            << " is failed. Geometry is cleared." << std::endl;

                        condition_itr->GetGeometry().clear();
                        condition_itr->Reset(ACTIVE);
                        condition_itr->Set(TO_ERASE);
                    }
                }
            }
        }
    }

    inline void ResetElementsAndNodes(ModelPart& rBackgroundGridModelPart)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rBackgroundGridModelPart.Elements().size()); ++i) {
            auto element_itr = rBackgroundGridModelPart.Elements().begin() + i;
            auto& r_geometry = element_itr->GetGeometry();
            element_itr->Reset(ACTIVE);

            for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j)
                r_geometry[j].Reset(ACTIVE);

        }
    }

    /**
     * @brief Search element connectivity for each particle
     * @details A search is performed to know in which grid element the material point falls.
     * If one or more material points fall in the grid element, the grid element is
     * set to be active and its connectivity is associated to the material point
     * element.
     * STEPS:
     * 1) All the elements are set to be INACTIVE
     * 2) A searching is performed and the grid elements which contain at least a MP are set to be ACTIVE
     *
     */
    template<std::size_t TDimension>
    void SearchElement(ModelPart& rBackgroundGridModelPart, ModelPart& rMPMModelPart, const std::size_t MaxNumberOfResults,
        const double Tolerance)
    {
        ResetElementsAndNodes(rBackgroundGridModelPart);

        std::vector<typename Element::Pointer> missing_elements;
        std::vector<typename Condition::Pointer> missing_conditions;

        NeighbourSearchElements(rMPMModelPart, rBackgroundGridModelPart, missing_elements, Tolerance);
        NeighbourSearchConditions(rMPMModelPart, rBackgroundGridModelPart, missing_conditions, Tolerance);

        if (missing_conditions.size() > 0 || missing_elements.size() > 0)
            BinBasedSearchElementsAndConditions<TDimension>(rMPMModelPart,
                rBackgroundGridModelPart, missing_elements, missing_conditions,
                MaxNumberOfResults, Tolerance);
    }
} // end namespace MPMSearchElementUtility

} // end namespace Kratos

#endif // KRATOS_MPM_SEARCH_ELEMENT_UTILITY

