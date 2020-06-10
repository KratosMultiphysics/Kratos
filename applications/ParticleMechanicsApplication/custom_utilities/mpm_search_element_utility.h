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

    typedef typename ModelPart::GeometryType GeometryType;

    typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> Boost2DPointType;

    typedef boost::geometry::model::polygon<Boost2DPointType> Boost2DPolygonType;

    typedef typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::IntegrationPointsArrayType IntegrationPointsArrayType;

    typedef Node<3> NodeType;

    inline double CrossProductDet2D(array_1d<double, 3> VectorA, array_1d<double, 3> VectorB)
    {
        return (VectorA[0] * VectorB[1] - VectorB[0] * VectorA[1]);
    }

    inline bool CheckIsInside(const GeometryType& rGeom, array_1d<double, 3>& LocalCoords, const array_1d<double, 3>& Coords, const double Tolerance)
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
                    rGeom.Points()[(i+1)% rGeom.PointsNumber()].Coordinates()- rGeom.Points()[i].Coordinates());
            }
            for (size_t i = 1; i < cross_products.size(); ++i)
            {
                if (cross_products[i] * cross_products[0] < 0.0)
                {
                    is_inside = false;
                    break;
                }
            }

        }

        if (is_inside) return rGeom.IsInside(Coords, LocalCoords, Tolerance);

        return false;
    }

    inline void ConstructNeighbourRelations(GeometryType& rGeom, const ModelPart& rBackgroundGridModelPart)
    {
        std::vector<typename Geometry<Node<3>>::Pointer> geometry_neighbours;
        for (IndexType j = 0; j < rBackgroundGridModelPart.NumberOfElements(); j++)
        {
            auto p_geometry_neighbour = (rBackgroundGridModelPart.ElementsBegin() + j)->pGetGeometry();
            for (IndexType n = 0; n < p_geometry_neighbour->size(); n++)
            {
                for (IndexType k = 0; k < rGeom.size(); k++)
                {
                    if (rGeom[k].Id() == (*p_geometry_neighbour)[n].Id()) {
                        geometry_neighbours.push_back(p_geometry_neighbour);
                        break;
                    }
                }
            }
        }
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


    inline void Check(IntegrationPointsArrayType& rIntergrationSubPoints, const double Tolerance, const Matrix& rN, const DenseVector<Matrix>& rDN_De)
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
                    << "3D PQMPM CAN ONLY BE USED FOR AXIS-ALIGNED RECTANGULAR-PRISM BACKGROUND GRIDS" << std::endl;
                KRATOS_ERROR << "ERROR" << std::endl;
            }
            rIntersectedGeometries[i]->BoundingBox(point_low, point_high);
            for (size_t j = 0; j < rIntersectedGeometries[i]->PointsNumber(); ++j) {
                for (size_t k = 0; k < 3; ++k) {
                    if (rIntersectedGeometries[i]->GetPoint(j).Coordinates()[k] != point_low[k]) {
                        if (rIntersectedGeometries[i]->GetPoint(j).Coordinates()[k] != point_high[k]) {
                            KRATOS_INFO("MPMSearchElementUtility::Check3DBackGroundMeshIsCubicAxisAligned - ")
                                << "3D PQMPM CAN ONLY BE USED FOR AXIS-ALIGNED RECTANGULAR-PRISM BACKGROUND GRIDS" << std::endl;
                            KRATOS_ERROR << "ERROR" << std::endl;
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

    inline void NeighbourSearchElements(const ModelPart& rMPMModelPart,
        const ModelPart& rBackgroundGridModelPart,
        std::vector<typename Element::Pointer>& rMissingElements,
        const double Tolerance,
        const ProcessInfo& rProcessInfo)
    {
        #pragma omp for
        for (int i = 0; i < static_cast<int>(rMPMModelPart.Elements().size()); ++i) {
            auto element_itr = rMPMModelPart.Elements().begin() + i;

            std::vector<array_1d<double, 3>> xg;
            element_itr->CalculateOnIntegrationPoints(MP_COORD, xg, rMPMModelPart.GetProcessInfo());

            GeometryType& r_parent_geometry = element_itr->GetGeometry().GetGeometryParent(0);

            array_1d<double, 3> local_coordinates;
            bool is_found = CheckIsInside(r_parent_geometry,local_coordinates,xg[0],Tolerance);
            if (!is_found)
            {
                if (!r_parent_geometry.Has(GEOMETRY_NEIGHBOURS))
                    ConstructNeighbourRelations(r_parent_geometry, rBackgroundGridModelPart);

                auto& geometry_neighbours = r_parent_geometry.GetValue(GEOMETRY_NEIGHBOURS);
                for (IndexType k = 0; k < geometry_neighbours.size(); k++)
                {
                    if (CheckIsInside(*geometry_neighbours[k], local_coordinates, xg[0], Tolerance))
                    {
                        auto p_new_geometry = CreateQuadraturePointsUtility<Node<3>>::CreateFromLocalCoordinates(
                            *(geometry_neighbours[k].get()), local_coordinates,
                            element_itr->GetGeometry().IntegrationPoints()[0].Weight());

                        if (!IsExplicitAndNeedsCorrection(p_new_geometry, rProcessInfo)) {
                            is_found = true;
                            // Update geometry of particle element
                            element_itr->SetGeometry(p_new_geometry);

                            for (IndexType j = 0; j < geometry_neighbours[k]->PointsNumber(); ++j)
                                geometry_neighbours[k]->Points()[j].Set(ACTIVE);
                            break;
                        }
                    }
                }
            }
            else {
                //pelem->Set(ACTIVE);

                auto p_new_geometry = CreateQuadraturePointsUtility<Node<3>>::CreateFromLocalCoordinates(
                    r_parent_geometry, local_coordinates,
                    element_itr->GetGeometry().IntegrationPoints()[0].Weight());

                if (IsExplicitAndNeedsCorrection(p_new_geometry, rProcessInfo)) is_found = false;
                else
                {
                    // Update geometry of particle element
                    element_itr->SetGeometry(p_new_geometry);

                    for (IndexType j = 0; j < r_parent_geometry.PointsNumber(); ++j)
                        r_parent_geometry[j].Set(ACTIVE);
                }
            }
            if (!is_found) {
                #pragma omp critical
                rMissingElements.push_back(&*element_itr);
            }
        }
    }

    inline void NeighbourSearchConditions(const ModelPart& rMPMModelPart,
        const ModelPart& rBackgroundGridModelPart,
        std::vector<typename Condition::Pointer>& rMissingConditions,
        const double Tolerance,
        const ProcessInfo& rProcessInfo)
    {
        #pragma omp for
        for (int i = 0; i < static_cast<int>(rMPMModelPart.Conditions().size()); ++i) {
            auto condition_itr = rMPMModelPart.Conditions().begin() + i;

            std::vector<array_1d<double, 3>> xg;
            condition_itr->CalculateOnIntegrationPoints(MPC_COORD, xg, rMPMModelPart.GetProcessInfo());
            if (xg.size() > 0)
            {
                GeometryType& r_parent_geometry = condition_itr->GetGeometry();

                array_1d<double, 3> local_coordinates;
                //bool is_found = r_parent_geometry.IsInside(xg[0], local_coordinates, Tolerance);
                bool is_found = CheckIsInside(r_parent_geometry, local_coordinates, xg[0], Tolerance);
                if (!is_found)
                {
                    if (!r_parent_geometry.Has(GEOMETRY_NEIGHBOURS))
                        ConstructNeighbourRelations(r_parent_geometry, rBackgroundGridModelPart);
                    auto& geometry_neighbours = r_parent_geometry.GetValue(GEOMETRY_NEIGHBOURS);
                    for (IndexType k = 0; k < geometry_neighbours.size(); k++)
                    {
                        //if (geometry_neighbours[k]->IsInside(xg[0], local_coordinates, Tolerance))
                        if (CheckIsInside(*geometry_neighbours[k], local_coordinates, xg[0], Tolerance))
                        {
                            is_found = true;

                            condition_itr->GetGeometry() = *(geometry_neighbours[k].get());

                            for (IndexType j = 0; j < geometry_neighbours[k]->PointsNumber(); ++j)
                                geometry_neighbours[k]->Points()[j].Set(ACTIVE);
                            break;
                        }
                    }
                }
                else {
                    //pelem->Set(ACTIVE);

                    condition_itr->GetGeometry() = r_parent_geometry;

                    for (IndexType j = 0; j < r_parent_geometry.PointsNumber(); ++j)
                        r_parent_geometry[j].Set(ACTIVE);
                }
                if (!is_found) {
                    #pragma omp critical
                    rMissingConditions.push_back(&*condition_itr);
                }
            }
        }
    }

    template <std::size_t TDimension>
    void BinBasedSearchElementsAndConditions(ModelPart& rMPMModelPart,
        ModelPart& rBackgroundGridModelPart,
        std::vector<typename Element::Pointer>& rMissingElements,
        std::vector<typename Condition::Pointer>& rMissingConditions,
        const std::size_t MaxNumberOfResults, const double Tolerance)
    {
        const ProcessInfo& r_process_info = rBackgroundGridModelPart.GetProcessInfo();
        const bool is_fix_explicit_mp_on_grid_edge = (r_process_info.Has(IS_FIX_EXPLICIT_MP_ON_GRID_EDGE))
            ? r_process_info.GetValue(IS_FIX_EXPLICIT_MP_ON_GRID_EDGE)
            : false;

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


                if (is_found && is_fix_explicit_mp_on_grid_edge) {
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
                        if (is_found) {
                            // store the nudged MP position
                            element_itr->SetValuesOnIntegrationPoints(MP_COORD, { xg_nudged }, rMPMModelPart.GetProcessInfo());
                            KRATOS_INFO("MPMSearchElementUtility") << "WARNING: To prevent spurious explicit stresses, Material Point " << element_itr->Id()
                                << " was nudged by " << nudge_displacement << std::endl;
                        }
                        else {
                            // find the un-nudged MP again
                            is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);
                            KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Material Point " << element_itr->Id()
                                << " lies exactly on an element edge and may give spurious results." << std::endl;
                        }
                    }
                }


                if (is_found == true) {
                    pelem->Set(ACTIVE);

                    auto p_new_geometry = CreateQuadraturePointsUtility<Node<3>>::CreateFromCoordinates(
                        pelem->pGetGeometry(), xg[0],
                        element_itr->GetGeometry().IntegrationPoints()[0].Weight());

                    // Update geometry of particle element
                    element_itr->SetGeometry(p_new_geometry);

                    for (IndexType j = 0; j < p_new_geometry->PointsNumber(); ++j)
                        (*p_new_geometry)[j].Set(ACTIVE);
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
                        //pelem->Set(ACTIVE);
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
        const ProcessInfo& r_process_info = rBackgroundGridModelPart.GetProcessInfo();

        ResetElementsAndNodes(rBackgroundGridModelPart);

        const bool use_neighbour_search = true; //TODO delete - for testing only

        std::vector<typename Element::Pointer> missing_elements;
        std::vector<typename Condition::Pointer> missing_conditions;

        if (use_neighbour_search)
        {
            NeighbourSearchElements(rMPMModelPart, rBackgroundGridModelPart, missing_elements, Tolerance, r_process_info);
            NeighbourSearchConditions(rMPMModelPart, rBackgroundGridModelPart, missing_conditions, Tolerance, r_process_info);
        }
        else
        {
            missing_elements.resize(rMPMModelPart.NumberOfElements());
            for (size_t i = 0; i < rMPMModelPart.NumberOfElements(); ++i)
            {
                missing_elements[i] = &*(rMPMModelPart.Elements().begin() + i);
            }

            missing_conditions.resize(rMPMModelPart.NumberOfConditions());
            for (size_t i = 0; i < rMPMModelPart.NumberOfConditions(); ++i)
            {
                missing_conditions[i] = &*(rMPMModelPart.Conditions().begin() + i);
            }
        }

        if (missing_conditions.size() > 0 || missing_elements.size() > 0)
            BinBasedSearchElementsAndConditions<TDimension>(rMPMModelPart,
                rBackgroundGridModelPart, missing_elements, missing_conditions,
                MaxNumberOfResults, Tolerance);

    }
} // end namespace MPMSearchElementUtility

} // end namespace Kratos

#endif // KRATOS_MPM_SEARCH_ELEMENT_UTILITY

