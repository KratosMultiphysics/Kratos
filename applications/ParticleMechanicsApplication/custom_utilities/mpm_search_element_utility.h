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
    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef Geometry<Node<3>> GeometryType;

    typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> Boost2DPointType;

    typedef boost::geometry::model::polygon<Boost2DPointType> Boost2DPolygonType;

    typedef typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::IntegrationPointsArrayType IntegrationPointsArrayType;

    typedef Node<3> NodeType;

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


    inline void Check(IntegrationPointsArrayType& rIntergrationSubPoints, const double Tolerance, const Matrix& rN, const DenseVector<Matrix>& rDN_De)
    {
        KRATOS_TRY

            double vol_frac_accum = 0.0;

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

            vol_frac_accum += rIntergrationSubPoints[i].Weight();
        }

        if (std::abs(vol_frac_accum - 1.0) < rIntergrationSubPoints.size() * Tolerance)
        {
            KRATOS_INFO("MPMSearchElementUtility::Check - ")
                << "Volume fraction of sub-points does not approximately sum to 1.0."
                << " This probably means the background grid is not big enough or that the PQMPM search factor is too small";
            KRATOS_ERROR << "ERROR";
        }

        for (size_t j = 0; j < rN.size2(); ++j) {
            SizeType nonzero_entries = 0;
            for (size_t i = 0; i < rIntergrationSubPoints.size(); i++) if (rN(i, j) != 0.0) nonzero_entries += 1;
            if (nonzero_entries != 1) {
                KRATOS_INFO("MPMSearchElementUtility::Check - ") << "There must be only one nonzero entry per shape function column!";
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
        bool is_coincident;
        for (size_t i = 0; i < rPoints.size(); ++i) {
            if (!rReferenceGeom.IsInside(rPoints[i], dummy_local_coords, Tolerance)) {
                // the test point may directly lie on one of the ref geom nodes - test this
                is_coincident = false;
                for (size_t j = 0; j < rReferenceGeom.PointsNumber(); ++j) {
                    if (norm_2(rPoints[i] - rReferenceGeom.GetPoint(j).Coordinates()) < Tolerance)
                    {
                        is_coincident = true;
                        break;
                    }
                }
                if (!is_coincident) return false;
            }
        }
        return true;

        KRATOS_CATCH("")
    }


    inline void Check3DBackGroundMeshIsCubicAxisAligned(const std::vector<typename GeometryType::Pointer> rIntersectedGeometries)
    {
        KRATOS_TRY

            NodeType point_low, point_high;
        for (size_t i = 0; i < rIntersectedGeometries.size(); ++i) {
            if (rIntersectedGeometries[i]->GetGeometryType() != GeometryData::Kratos_Hexahedra3D8) {
                KRATOS_INFO("MPMSearchElementUtility::Check3DBackGroundMeshIsCubicAxisAligned - ")
                    << "3D PQMPM CAN ONLY BE USED FOR AXIS-ALIGNED RECTANGULAR-PRISM BACKGROUND GRIDS";
                KRATOS_ERROR << "ERROR";
            }
            rIntersectedGeometries[i]->BoundingBox(point_low, point_high);
            for (size_t j = 0; j < rIntersectedGeometries[i]->PointsNumber(); ++j) {
                for (size_t k = 0; k < 3; ++k) {
                    if (rIntersectedGeometries[i]->GetPoint(j).Coordinates()[k] != point_low[k]) {
                        if (rIntersectedGeometries[i]->GetPoint(j).Coordinates()[k] != point_high[k]) {
                            KRATOS_INFO("MPMSearchElementUtility::Check3DBackGroundMeshIsCubicAxisAligned - ")
                                << "3D PQMPM CAN ONLY BE USED FOR AXIS-ALIGNED RECTANGULAR-PRISM BACKGROUND GRIDS";
                            KRATOS_ERROR << "ERROR";
                        }
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }


    inline bool CheckNoPointsAreInGeom(
        const std::vector<array_1d<double, 3>>& rPoints,
        const GeometryType& rReferenceGeom,
        const double Tolerance)
    {
        KRATOS_TRY

            array_1d<double, 3> dummy_local_coords;
        for (size_t i = 0; i < rPoints.size(); ++i) {
            if (rReferenceGeom.IsInside(rPoints[i], dummy_local_coords), Tolerance) return false;
        }
        return true;

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
        typename Geometry<Node<3>>::PointsArrayType rPoints)
    {
        KRATOS_TRY

            if (WorkingSpaceDimension == 1 && LocalSpaceDimension == 1)
                return Kratos::make_shared<
                QuadraturePointPartitionedGeometry<Node<3>, 1>>(
                    rPoints, rShapeFunctionContainer);
            else if (WorkingSpaceDimension == 2 && LocalSpaceDimension == 1)
                return Kratos::make_shared<
                QuadraturePointPartitionedGeometry<Node<3>, 2, 1>>(
                    rPoints, rShapeFunctionContainer);
            else if (WorkingSpaceDimension == 2 && LocalSpaceDimension == 2)
                return Kratos::make_shared<
                QuadraturePointPartitionedGeometry<Node<3>, 2>>(
                    rPoints, rShapeFunctionContainer);
            else if (WorkingSpaceDimension == 3 && LocalSpaceDimension == 2)
                return Kratos::make_shared<
                QuadraturePointPartitionedGeometry<Node<3>, 3, 2>>(
                    rPoints, rShapeFunctionContainer);
            else if (WorkingSpaceDimension == 3 && LocalSpaceDimension == 3)
                return Kratos::make_shared<
                QuadraturePointPartitionedGeometry<Node<3>, 3>>(
                    rPoints, rShapeFunctionContainer);
            else {
                KRATOS_ERROR << "Working/Local space dimension combinations are "
                    << "not provided for QuadraturePointGeometry. WorkingSpaceDimension: "
                    << WorkingSpaceDimension << ", LocalSpaceDimension: " << LocalSpaceDimension
                    << std::endl;
            }
        KRATOS_CATCH("")
    }


    inline bool DetermineIfDomainOverlapsBoundaryConditions(std::vector<typename GeometryType::Pointer>& IntersectedGeometries,
        const array_1d<double, 3>& rCoordinates, const double RangeFactor, const double SideHalfLength)
    {
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
                    if (fixed_point_to_cog < RangeFactor * SideHalfLength) return true;
                }
            }
        }

        return false;
    }


    inline typename Geometry<Node<3>>::Pointer PartitionMasterMaterialPointsIntoSubPoints(const ModelPart& rBackgroundGridModelPart,
        const array_1d<double, 3>& rCoordinates,
        Element& rMasterMaterialPoint,
        const typename Geometry<Node<3>>::Pointer pGeometry,
        const double Tolerance)
    {
        KRATOS_TRY;

        const SizeType working_dim = pGeometry->WorkingSpaceDimension();
        const bool is_axisymmetric = (rBackgroundGridModelPart.GetProcessInfo().Has(IS_AXISYMMETRIC))
            ? rBackgroundGridModelPart.GetProcessInfo().GetValue(IS_AXISYMMETRIC) : false;
        const double pqmpm_search_factor = (rBackgroundGridModelPart.GetProcessInfo().Has(PQMPM_SEARCH_FACTOR))
            ? rBackgroundGridModelPart.GetProcessInfo().GetValue(PQMPM_SEARCH_FACTOR) : 0.0;
        if (pqmpm_search_factor < 0.0)
        {
            KRATOS_INFO("MPMSearchElementUtility::PartitionMasterMaterialPointsIntoSubPoints - ")
                << "The PQMPM search factor enables fast filtering of background grid intersections and must be positive."
                << " It should be 2-5x larger the maximum aspect ratio of the most distorted background grid element."
                << " If problems presist, disable the fast filtering by setting pqmpm_search_factor = 0.0\n";
            KRATOS_ERROR << "ERROR";
        }


        // Get volume and set up master domain bounding points
        std::vector<double> mp_volume_vec;
        rMasterMaterialPoint.CalculateOnIntegrationPoints(MP_VOLUME, mp_volume_vec, rBackgroundGridModelPart.GetProcessInfo());
        const double side_half_length = std::pow(mp_volume_vec[0], 1.0 / double(working_dim)) / 2.0;
        const SizeType n_bounding_box_vertices = std::pow(2.0, working_dim);
        std::vector<array_1d<double, 3>> master_domain_points(n_bounding_box_vertices);
        CreateBoundingBoxPoints(master_domain_points, rCoordinates, side_half_length, working_dim);

        // If axisymmetric, we can't make a sub-point with x<0.
        if (is_axisymmetric) {
            if ((rCoordinates[0] - side_half_length) < std::numeric_limits<double>::epsilon()) {
                return CreateQuadraturePointsUtility<Node<3>>::CreateFromCoordinates(
                    pGeometry, rCoordinates, rMasterMaterialPoint.GetGeometry().IntegrationPoints()[0].Weight());
            }
        }

        // Initially check if the bounding box volume scalar is less than the element volume scalar
        if (mp_volume_vec[0] <= pGeometry->DomainSize()) {
            if (CheckAllPointsAreInGeom(master_domain_points, *pGeometry, Tolerance)) {
                return CreateQuadraturePointsUtility<Node<3>>::CreateFromCoordinates(
                    pGeometry, rCoordinates, rMasterMaterialPoint.GetGeometry().IntegrationPoints()[0].Weight());
            }
        }

        // we need to do splitting. Initially determine all grid elements we intersect with
        const double z_mod = (working_dim == 3) ? 1.0 : 0.0;
        const Point point_low(rCoordinates[0] - side_half_length, rCoordinates[1] - side_half_length, rCoordinates[2] - z_mod * side_half_length);
        const Point point_high(rCoordinates[0] + side_half_length, rCoordinates[1] + side_half_length, rCoordinates[2] + z_mod * side_half_length);
        const double range_factor = (working_dim == 3) ? 2.0 : 1.414214; // 45 deg for each dim
        double center_to_center, maximum_contact_range, char_length;
        NodeType ele_point_low, ele_point_high;

        SizeType number_of_nodes = 0;
        std::vector<typename GeometryType::Pointer> intersected_geometries;
        std::vector<typename Element::Pointer> intersected_elements;

        for (auto ele_it : rBackgroundGridModelPart.Elements())
        {
            char_length = std::pow(ele_it.GetGeometry().DomainSize(), 1.0 / double(working_dim)) * pqmpm_search_factor + side_half_length;
            if (pqmpm_search_factor == 0.0 || std::abs(ele_it.GetGeometry().Center().X() - rCoordinates[0]) < char_length) {
                if (pqmpm_search_factor == 0.0 || std::abs(ele_it.GetGeometry().Center().Y() - rCoordinates[1]) < char_length) {
                    if (pqmpm_search_factor == 0.0 || working_dim == 2 || std::abs(ele_it.GetGeometry().Center().Z() - rCoordinates[2]) < char_length) {
                        center_to_center = norm_2(ele_it.GetGeometry().Center() - rCoordinates);
                        ele_it.GetGeometry().BoundingBox(ele_point_low, ele_point_high);
                        maximum_contact_range = range_factor * side_half_length + norm_2(ele_point_high - ele_point_low) / 2.0;
                        if (center_to_center <= maximum_contact_range)
                        {
                            if (ele_it.GetGeometry().HasIntersection(point_low, point_high)) {
                                number_of_nodes += ele_it.GetGeometry().PointsNumber();
                                intersected_geometries.push_back(ele_it.pGetGeometry());
                                intersected_elements.push_back(&ele_it);
                            }
                        }
                    }
                }
            }
        }

        // Prevent splitting particles over fixed nodes
        if (DetermineIfDomainOverlapsBoundaryConditions(intersected_geometries, rCoordinates, range_factor, side_half_length))
            return CreateQuadraturePointsUtility<Node<3>>::CreateFromCoordinates(
                pGeometry, rCoordinates, rMasterMaterialPoint.GetGeometry().IntegrationPoints()[0].Weight());

        // If we are 3D, check background mesh are axis-aligned perfect cubes
        if (working_dim == 3)  Check3DBackGroundMeshIsCubicAxisAligned(intersected_geometries);

        // Prepare containers to hold all sub-points
        const SizeType number_of_sub_material_points = intersected_geometries.size();
        PointerVector<Node<3>> nodes_list(number_of_nodes);
        IntegrationPointsArrayType ips(number_of_sub_material_points);
        Matrix N_matrix(number_of_sub_material_points, number_of_nodes, 0.0);
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

            if (CheckNoPointsAreInGeom(master_domain_points, *intersected_geometries[i], Tolerance)) {
                // whole element is completely inside bounding box

                trial_subpoint = CreateSubPoint(intersected_geometries[i]->Center(),
                    intersected_geometries[i]->DomainSize() / mp_volume_vec[0],
                    *intersected_geometries[i], N, DN_De);
            }
            else {
                // only some of the background element is within the bounding box - most expensive check

                if (working_dim == 2) {
                    Determine2DSubPoint(*intersected_geometries[i], master_domain_points, sub_point_position, sub_point_volume);
                    sub_point_position[2] = rCoordinates[2]; // set z coord of sub point to that of the master
                }
                else
                    Determine3DSubPoint(*intersected_geometries[i], master_domain_points, sub_point_position, sub_point_volume);
                trial_subpoint = CreateSubPoint(sub_point_position, sub_point_volume / mp_volume_vec[0],
                    *intersected_geometries[i], N, DN_De);
            }

            // Transfer local data to containers
            if (trial_subpoint.Weight() > std::numeric_limits<double>::epsilon())
            {
                intersected_elements[i]->Set(ACTIVE);
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

        if (active_subpoint_index == 1) return CreateQuadraturePointsUtility<Node<3>>::CreateFromCoordinates(
            pGeometry, rCoordinates, rMasterMaterialPoint.GetGeometry().IntegrationPoints()[0].Weight());

        IntegrationPointsArrayType ips_active(active_subpoint_index);
        PointerVector<Node<3>> nodes_list_active(active_node_index);
        if (ips_active.size() == ips.size())
        {
            ips_active = ips;
            nodes_list_active = nodes_list;
        }
        else
        {
            N_matrix.resize(active_subpoint_index, active_node_index, true);
            DN_De_vector.resize(active_subpoint_index, true);
            for (size_t i = 0; i < active_subpoint_index; i++) ips_active[i] = ips[i];
            for (size_t i = 0; i < active_node_index; i++) nodes_list_active(i) = nodes_list(i);
        }

        Check(ips_active, std::numeric_limits<double>::epsilon(), N_matrix, DN_De_vector);

        GeometryData::IntegrationMethod ThisDefaultMethod = pGeometry->GetDefaultIntegrationMethod();
        typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::IntegrationPointsContainerType ips_container;
        ips_container[ThisDefaultMethod] = ips_active;
        typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::ShapeFunctionsValuesContainerType shape_function_container;
        shape_function_container[ThisDefaultMethod] = N_matrix;
        typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::ShapeFunctionsLocalGradientsContainerType shape_function_derivatives_container;
        shape_function_derivatives_container[ThisDefaultMethod] = DN_De_vector;

        GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
            ThisDefaultMethod,
            ips_container,
            shape_function_container,
            shape_function_derivatives_container);

        return CreateCustomQuadraturePoint(working_dim, pGeometry->LocalSpaceDimension(), data_container, nodes_list_active);

        KRATOS_CATCH("");
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
        const ProcessInfo& r_process_info = rBackgroundGridModelPart.GetProcessInfo();
        const bool is_explicit = (r_process_info.Has(IS_EXPLICIT))
            ? r_process_info.GetValue(IS_EXPLICIT)
            : false;
        const bool is_pqmpm = (r_process_info.Has(IS_PQMPM))
            ? r_process_info.GetValue(IS_PQMPM)
            : false;

        // Reset elements to inactive
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rBackgroundGridModelPart.Elements().size()); ++i) {
            auto element_itr = rBackgroundGridModelPart.Elements().begin() + i;
            auto& r_geometry = element_itr->GetGeometry();
            element_itr->Reset(ACTIVE);

            for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j)
                r_geometry[j].Reset(ACTIVE);

        }
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
            for (int i = 0; i < static_cast<int>(rMPMModelPart.Elements().size()); ++i) {
                auto element_itr = rMPMModelPart.Elements().begin() + i;

                std::vector<array_1d<double, 3>> xg;
                element_itr->CalculateOnIntegrationPoints(MP_COORD, xg, rMPMModelPart.GetProcessInfo());
                typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();

                Element::Pointer pelem;

                // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
                bool is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);


                if (is_found && is_explicit && !is_pqmpm) {
                    // check if MP is exactly on the edge of the element, this gives spurious strains in explicit
                    bool isOnEdge = false;
                    for (SizeType i = 0; i < N.size(); ++i) {
                        if (std::abs(N[i]) < std::numeric_limits<double>::epsilon()) {
                            isOnEdge = true;
                            break;
                        }
                    }
                    if (isOnEdge) {
                        // MP is exactly on the edge. Now we give it a little 'nudge'
                        array_1d<double, 3> xg_nudged = array_1d<double, 3>(xg[0]);
                        const double& delta_time = r_process_info[DELTA_TIME];
                        std::vector<array_1d<double, 3>> mp_vel;
                        element_itr->CalculateOnIntegrationPoints(MP_VELOCITY, mp_vel, rMPMModelPart.GetProcessInfo());
                        array_1d<double, 3> nudge_displacement = delta_time / 1000.0 * mp_vel[0];
                        xg_nudged += nudge_displacement;
                        is_found = SearchStructure.FindPointOnMesh(xg_nudged, N, pelem, result_begin, MaxNumberOfResults, Tolerance);
                        // check if the nudged point is found...
                        if (is_found){
                            // store the nudged MP position
                            element_itr->SetValuesOnIntegrationPoints(MP_COORD, { xg_nudged }, rMPMModelPart.GetProcessInfo());
                            KRATOS_INFO("MPMSearchElementUtility") << "WARNING: To prevent spurious explicit stresses, Material Point " << element_itr->Id()
                                << " was nudged by " << nudge_displacement << std::endl;
                        }
                        else {
                            // find the un-nudged MP again
                            is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);
                            KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Material Point " << element_itr->Id()
                                << " lies exactly on an element edge and may give spurious results."<< std::endl;
                        }
                    }
                }

                if (is_found == true) {
                    pelem->Set(ACTIVE);

                    typename GeometryType::Pointer p_new_geometry = (is_pqmpm)
                        ? PartitionMasterMaterialPointsIntoSubPoints(rBackgroundGridModelPart,
                            xg[0], *element_itr, pelem->pGetGeometry(), Tolerance)
                        : CreateQuadraturePointsUtility<Node<3>>::CreateFromCoordinates(
                            pelem->pGetGeometry(), xg[0],
                            element_itr->GetGeometry().IntegrationPoints()[0].Weight());

                    // Update geometry of particle element
                    element_itr->SetGeometry(p_new_geometry);

                    for (IndexType j = 0; j < p_new_geometry->PointsNumber(); ++j) (*p_new_geometry)[j].Set(ACTIVE);
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
            for (int i = 0; i < static_cast<int>(rMPMModelPart.Conditions().size()); ++i) {

                auto condition_itr = rMPMModelPart.Conditions().begin() + i;
                std::vector<array_1d<double, 3>> xg;
                condition_itr->CalculateOnIntegrationPoints(MPC_COORD, xg, rMPMModelPart.GetProcessInfo());

                if (xg.size() == 1) {
                    // Only search for particle based BCs!
                    // Grid BCs are still applied on MP_model_part but we don't want to search for them.
                    typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();

                    Element::Pointer pelem;

                    // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
                    bool is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);

                    if (is_found == true) {
                        pelem->Set(ACTIVE);
                        condition_itr->GetGeometry() = pelem->GetGeometry();
                        auto& r_geometry = condition_itr->GetGeometry();

                        for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j)
                            r_geometry[j].Set(ACTIVE);
                    }
                    else {
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



} // end namespace MPMSearchElementUtility
} // end namespace Kratos
#endif // KRATOS_MPM_SEARCH_ELEMENT_UTILITY

