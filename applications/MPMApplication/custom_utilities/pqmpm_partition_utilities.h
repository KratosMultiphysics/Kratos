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
#include "mpm_application_variables.h"
#include "geometries/geometry_shape_function_container.h"
#include "custom_geometries/quadrature_point_partitioned_geometry.h"

#include "geometries/geometry.h"
#include "includes/model_part.h"

#include "boost/geometry/geometry.hpp"

namespace Kratos
{
    ///@name Kratos Classes
    ///@{
    /**
     * @class PQMPMPartitionUtilities
     * @ingroup MPMApplication
     * @brief This class includes several utilities necessaries for the PQMPM
     * @author Peter Wilson
     */
    class KRATOS_API(MPM_APPLICATION) PQMPMPartitionUtilities
    {
    public:
        // Standard types
        typedef std::size_t IndexType;
        typedef std::size_t SizeType;
        typedef Node NodeType;
        typedef typename ModelPart::GeometryType GeometryType;
        typedef typename Geometry<Node>::Pointer GeometryNodePointerType;

        // Boost Polygon
        typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> Boost2DPointType;
        typedef boost::geometry::model::polygon<Boost2DPointType> Boost2DPolygonType;

        // Container types
        typedef typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::IntegrationPointsArrayType IntegrationPointsArrayType;
        typedef typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::IntegrationPointsContainerType IntegrationPointsContainerType;
        typedef typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;
        typedef typename GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

        /// Constructor
        PQMPMPartitionUtilities() = default;

        /// Destructor.
        virtual ~PQMPMPartitionUtilities() = default;

        static void PartitionMasterMaterialPointsIntoSubPoints(const ModelPart& rBackgroundGridModelPart,
            const array_1d<double, 3>& rCoordinates,
            const array_1d<double, 3>& rLocalCoords,
            Element& rMasterMaterialPoint,
            typename GeometryType::Pointer pQuadraturePointGeometry,
            const double Tolerance);

    private:
        static void CreateBoundingBoxPoints(std::vector<array_1d<double, 3>>& rPointVector,
            const array_1d<double, 3>& rCenter, const double SideHalfLength, const SizeType WorkingDim);

        static bool CheckAllPointsAreInGeom(
            const std::vector<array_1d<double, 3>>& rPoints,
            const GeometryType& rReferenceGeom,
            const double Tolerance);

        static void Check3DBackGroundMeshIsCubicAxisAligned(const std::vector<GeometryType*>& rIntersectedGeometries);

        static Boost2DPolygonType Create2DPolygonBoundingSquareFromPointsFast(const std::vector<array_1d<double, 3>>& rPoints,
            const bool XActive = true, const bool YActive = true, const bool ZActive = false);

        static Boost2DPolygonType Create2DPolygonFromGeometryFast(const GeometryType& rGeom,
            const bool XActive = true, const bool YActive = true, const bool ZActive = false);

        static IntegrationPoint<3> CreateSubPoint(const array_1d<double, 3>& rGlobalCoords, const double rVolumeFraction,
            const GeometryType& rBackgroundGridElementGeom, Vector& rN, Matrix& rDN_De);

        static void Determine2DSubPoint(const GeometryType& rGridElement, const std::vector<array_1d<double, 3>>& rMasterDomainPoints,
            array_1d<double, 3>& rSubPointCoord, double& rSubPointVolume);

        static void Determine3DSubPoint(const GeometryType& rGridElement, const std::vector<array_1d<double, 3>>& rMasterDomainPoints,
            array_1d<double, 3>& rSubPointCoord, double& rSubPointVolume);

        static void CheckPQMPM(IntegrationPointsArrayType& rIntergrationSubPoints,
            const double Tolerance, Matrix& rN, const DenseVector<Matrix>& rDN_De);

        static bool DetermineIfDomainOverlapsBoundaryConditions(std::vector<GeometryType*>& IntersectedGeometries,
            const array_1d<double, 3>& rCoordinates, const double SideHalfLength);

        static bool IntersectionCheckWithBoundingBox(const GeometryType& rGeom,
            const array_1d<double, 3>& rCoord, const double SideHalfLength);

        static void RecursivePQMPMNeighbourSearch(const ModelPart& rBackgroundGridModelPart,
            std::vector<GeometryType*>& rIntersectedGeometries,
            const Point& rPointLow,
            const Point& rPointHigh,
            IndexType& RecursionCount,
            const array_1d<double, 3>& rCoordinates,
            const double SideHalfLength,
            const SizeType MaxRecursions = 100);

        static bool CheckFixedNodesWithinBoundingBox(const PointerVector<Node>& rNodesList,
            const Point& rPointHigh, const Point& rPointLow, const SizeType WorkingDim);

        static GeometryNodePointerType CreateCustomQuadraturePoint(
            SizeType WorkingSpaceDimension,
            SizeType LocalSpaceDimension,
            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rShapeFunctionContainer,
            typename Geometry<Node>::PointsArrayType rPoints,
            GeometryType* pGeometryParent);
    }; // class PQMPMPartitionUtilities
}  // namespace Kratos
#endif /* KRATOS_PQMPM_PARTITION_UTILITIES defined */