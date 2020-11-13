//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//

#ifndef KRATOS_PQMPM_PARTITION_UTILITIES
#define KRATOS_PQMPM_PARTITION_UTILITIES

// Project includes



#include "includes/define.h"
#include "utilities/quadrature_points_utility.h"

#include "particle_mechanics_application_variables.h"
#include "geometries/geometry_shape_function_container.h"
#include "custom_geometries/quadrature_point_partitioned_geometry.h"

#include "geometries/geometry.h"
#include "includes/model_part.h"

#include "boost/geometry/geometry.hpp"
#include "boost/geometry/geometries/register/point.hpp"
#include "boost/geometry/geometries/register/ring.hpp"


#include "boost/geometry/geometry.hpp"
#include "boost/geometry/geometries/register/point.hpp"
#include "boost/geometry/geometries/register/ring.hpp"


namespace Kratos
{
    ///@name Kratos Classes
    ///@{
    /**
     * @namespace PQMPMPartitionUtilities
     * @ingroup ParticleMechanicsApplication
     * @brief This namespace includes several utilities necessaries for the PQMPM
     * @author Peter Wilson
     */
    namespace PQMPMPartitionUtilities
    {
        // Standard types
        typedef std::size_t IndexType;
        typedef std::size_t SizeType;
        typedef Node<3> NodeType;
        typedef typename ModelPart::GeometryType GeometryType;
        typedef typename Geometry<Node<3>>::Pointer GeometryNodePointerType;

        // Boost Polygon
        typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> Boost2DPointType;
        typedef boost::geometry::model::polygon<Boost2DPointType> Boost2DPolygonType;

        // Container types
        typedef typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::IntegrationPointsArrayType IntegrationPointsArrayType;
        typedef typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::IntegrationPointsContainerType IntegrationPointsContainerType;
        typedef typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;
        typedef typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

        void KRATOS_API(PARTICLE_MECHANICS_APPLICATION) PartitionMasterMaterialPointsIntoSubPoints(const ModelPart& rBackgroundGridModelPart,
            const array_1d<double, 3>& rCoordinates,
            const array_1d<double, 3>& rLocalCoords,
            Element& rMasterMaterialPoint,
            typename GeometryType::Pointer pQuadraturePointGeometry,
            const double Tolerance);

        void KRATOS_API(PARTICLE_MECHANICS_APPLICATION) CreateBoundingBoxPoints(std::vector<array_1d<double, 3>>& rPointVector,
            const array_1d<double, 3>& rCenter, const double SideHalfLength, const SizeType WorkingDim);

        bool KRATOS_API(PARTICLE_MECHANICS_APPLICATION) CheckAllPointsAreInGeom(
            const std::vector<array_1d<double, 3>>& rPoints,
            const GeometryType& rReferenceGeom,
            const double Tolerance);

        void KRATOS_API(PARTICLE_MECHANICS_APPLICATION) Check3DBackGroundMeshIsCubicAxisAligned(const std::vector<GeometryType*>& rIntersectedGeometries);

        Boost2DPolygonType KRATOS_API(PARTICLE_MECHANICS_APPLICATION) Create2DPolygonBoundingSquareFromPointsFast(const std::vector<array_1d<double, 3>>& rPoints,
            const bool XActive = true, const bool YActive = true, const bool ZActive = false);

        Boost2DPolygonType KRATOS_API(PARTICLE_MECHANICS_APPLICATION) Create2DPolygonFromGeometryFast(const GeometryType& rGeom,
            const bool XActive = true, const bool YActive = true, const bool ZActive = false);

        IntegrationPoint<3> KRATOS_API(PARTICLE_MECHANICS_APPLICATION) CreateSubPoint(const array_1d<double, 3>& rGlobalCoords, const double rVolumeFraction,
            const GeometryType& rBackgroundGridElementGeom, Vector& rN, Matrix& rDN_De);


        void KRATOS_API(PARTICLE_MECHANICS_APPLICATION) Determine2DSubPoint(const GeometryType& rGridElement, const std::vector<array_1d<double, 3>>& rMasterDomainPoints,
            array_1d<double, 3>& rSubPointCoord, double& rSubPointVolume);

        void KRATOS_API(PARTICLE_MECHANICS_APPLICATION) Determine3DSubPoint(const GeometryType& rGridElement, const std::vector<array_1d<double, 3>>& rMasterDomainPoints,
            array_1d<double, 3>& rSubPointCoord, double& rSubPointVolume);

        void KRATOS_API(PARTICLE_MECHANICS_APPLICATION) CheckPQMPM(IntegrationPointsArrayType& rIntergrationSubPoints,
            const double Tolerance, Matrix& rN, const DenseVector<Matrix>& rDN_De);

        bool KRATOS_API(PARTICLE_MECHANICS_APPLICATION) DetermineIfDomainOverlapsBoundaryConditions(std::vector<GeometryType*>& IntersectedGeometries,
            const array_1d<double, 3>& rCoordinates, const double SideHalfLength);

        bool KRATOS_API(PARTICLE_MECHANICS_APPLICATION) IntersectionCheckWithBoundingBox(const GeometryType& rGeom,
            const array_1d<double, 3>& rCoord, const double SideHalfLength);

        void KRATOS_API(PARTICLE_MECHANICS_APPLICATION) RecursivePQMPMNeighbourSearch(const ModelPart& rBackgroundGridModelPart,
            std::vector<GeometryType*>& rIntersectedGeometries,
            const Point& rPointLow,
            const Point& rPointHigh,
            IndexType& RecursionCount,
            const array_1d<double, 3>& rCoordinates,
            const double SideHalfLength,
            const SizeType MaxRecursions = 100);

        bool KRATOS_API(PARTICLE_MECHANICS_APPLICATION) CheckFixedNodesWithinBoundingBox(const PointerVector<Node<3>>& rNodesList,
            const Point& rPointHigh, const Point& rPointLow, const SizeType WorkingDim);

        GeometryNodePointerType KRATOS_API(PARTICLE_MECHANICS_APPLICATION) CreateCustomQuadraturePoint(
            SizeType WorkingSpaceDimension,
            SizeType LocalSpaceDimension,
            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rShapeFunctionContainer,
            typename Geometry<Node<3>>::PointsArrayType rPoints,
            GeometryType* pGeometryParent);
    }; // namespace PQMPMPartitionUtilities
}  // namespace Kratos
#endif /* KRATOS_PQMPM_PARTITION_UTILITIES defined */