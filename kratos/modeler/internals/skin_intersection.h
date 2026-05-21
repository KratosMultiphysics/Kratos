//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/ublas_interface.h"
#include "utilities/divide_triangle_3d_3.h"

namespace Kratos::Internals {

class SkinIntersection
{
public:
    using SplitGeometryType = DivideTriangle3D3<Geometry<Point>>::IndexedPointGeometryType;
    using SplitGeometryPointerType = SplitGeometryType::Pointer;
    using SplitGeometryPointType = DivideTriangle3D3<Geometry<Point>>::IndexedPointType;
    using SplitGeometryPointsArrayType = Kratos::PointerVector<SplitGeometryPointType>;

    SkinIntersection(const array_1d<std::size_t,3>& NumDivisions):
        mNumDivisions(NumDivisions),
        mCutGeometries(NumDivisions[0]*NumDivisions[1]*NumDivisions[2])
    {}

    ~SkinIntersection() = default;


    template<class GeometryType>
    static const SplitGeometryPointerType ConvertToWorkGeometry(const GeometryType& rGeometry) {
        SplitGeometryPointsArrayType points;
        points.reserve(rGeometry.size());
        std::size_t index = 0;
        for (const auto& point: rGeometry) {
            points.push_back(Kratos::make_shared<SplitGeometryPointType>(point.Coordinates(), ++index));
        }
        return Kratos::make_shared<Triangle3D3<SplitGeometryPointType>>(points);
    }


    void ComputeCutInsideCell(
        SplitGeometryPointerType pGeometry,
        const Point& rLowPoint,
        const Point& rHighPoint,
        std::size_t I,
        std::size_t J,
        std::size_t K)
    {
        std::vector<SplitGeometryPointerType> geometries;
        geometries.push_back(pGeometry);

        for (std::size_t i = 0; i < 3; i++) {
            geometries = TrimLower<SplitGeometryType>(geometries, rLowPoint, i);
            geometries = TrimHigher<SplitGeometryType>(geometries, rHighPoint, i);
        }

        auto& r_cell_geometries = CutGeometriesInCell(I, J, K);
        r_cell_geometries.insert(r_cell_geometries.end(), geometries.begin(), geometries.end());
    }

    const std::vector<SplitGeometryPointerType>& CutGeometriesInCell(std::size_t I, std::size_t J, std::size_t K) const
    {
        return mCutGeometries[K + mNumDivisions[2]*(J + mNumDivisions[1]*I)];
    }


    double GetIntegrationPointData(
        std::vector<Vector>& rIntegrationPointData,
        const Point& rLowPoint,
        const Point& rHighPoint,
        std::size_t I,
        std::size_t J,
        std::size_t K,
        const GeometryData::IntegrationMethod Method,
        double Thickness) const
    {
        rIntegrationPointData.resize(0);
        auto& r_cut_geometries = CutGeometriesInCell(I, J, K);

        for (auto& p_geom: r_cut_geometries) {
            AddTriangleIntegrationPointData(rIntegrationPointData, *p_geom, Method);
        }

        // Assuming that the element is a hexahedra with plane faces and local coordinates in [-1,1]^3
        array_1d<double,3> len = rHighPoint - rLowPoint;
        double vol = len[0]*len[1]*len[2];
        double weight_factor = Thickness*8.0/vol; // Integration point weights add to 8 for hexas
        double effective_area = 0.0;

        for (auto& r_data: rIntegrationPointData) {
            for (std::size_t d = 0; d < 3; d++) {
                r_data[d] = 2*(r_data[d] - rLowPoint[d])/len[d] - 1;
            }
            effective_area += r_data[3];
            r_data[3] *= weight_factor;
            //TODO: PointIsOnSharedFace should account for same-color neighbors only
            //if (PointIsOnSharedFace(r_data, I, J, K)) {
            //    r_data[3] *= 0.5*weight_factor;
            //}
            //else {
            //    r_data[3] *= weight_factor;
            //}
        }

        const double full_fraction = Thickness*effective_area / vol;
        return full_fraction;
    }

private:
    array_1d<std::size_t, 3> mNumDivisions;
    std::vector<std::vector<SplitGeometryPointerType>> mCutGeometries = {};

    std::vector<SplitGeometryPointerType>& CutGeometriesInCell(std::size_t I, std::size_t J, std::size_t K)
    {
        return mCutGeometries[K + mNumDivisions[2]*(J + mNumDivisions[1]*I)];
    }

    template<class GeometryType>
    static std::size_t SignedDistanceUnder(
        Vector& rNodalDistances,
        const GeometryType& rGeometry,
        const Point& rReferencePoint,
        const std::size_t Component)
    {
        auto num_nodes = rGeometry.PointsNumber();
        if (rNodalDistances.size() != num_nodes) {
            rNodalDistances.resize(num_nodes);
            noalias(rNodalDistances) = ZeroVector(num_nodes);
        }

        auto idistance = rNodalDistances.begin();
        const double reference_coord = rReferencePoint[Component];
        std::size_t points_under = 0;

        for (const auto& point: rGeometry) {
            const double d = reference_coord - point[Component];
            if (d >= 0) points_under++;
            *(idistance++) = d;
        }

        return points_under;
    }



    template<class GeometryType>
    static std::size_t SignedDistanceOver(
        Vector& rNodalDistances,
        const GeometryType& rGeometry,
        const Point& rReferencePoint,
        const std::size_t Component)
    {
        auto num_nodes = rGeometry.PointsNumber();
        if (rNodalDistances.size() != num_nodes) {
            rNodalDistances.resize(num_nodes);
            noalias(rNodalDistances) = ZeroVector(num_nodes);
        }

        auto idistance = rNodalDistances.begin();
        const double reference_coord = rReferencePoint[Component];
        std::size_t points_over = 0;

        for (const auto& point: rGeometry) {
            const double d = point[Component] - reference_coord;
            if (d >= 0) points_over++;
            *(idistance++) = d;
        }

        return points_over;
    }

    template<class GeometryType, class GeometryPointerType = typename GeometryType::Pointer>
    static std::vector<GeometryPointerType> TrimHigher(
        const std::vector<GeometryPointerType>& rGeometries,
        const Point& rLowPoint,
        std::size_t Component)
    {
        using PointType = typename GeometryType::PointType;
        std::vector<GeometryPointerType> cut, trimmed_geometries;
        Vector nodal_distances;

        for (auto iter = rGeometries.begin(); iter < rGeometries.end();  ++iter) {
            auto& r_geom = **iter;
            cut.resize(0);

            std::size_t points_under = SignedDistanceUnder(nodal_distances, r_geom, rLowPoint, Component);
            if (points_under == 3) {
                trimmed_geometries.push_back(*iter);
            }
            else if (points_under > 0) {
                auto splitter = DivideTriangle3D3<PointType>(r_geom, nodal_distances);
                splitter.GenerateDivision();
                cut = splitter.GetPositiveSubdivisions();
                trimmed_geometries.insert(trimmed_geometries.end(), cut.begin(), cut.end());
            }
        }

        return trimmed_geometries;
    }

    template<class GeometryType, class GeometryPointerType = typename GeometryType::Pointer>
    static std::vector<GeometryPointerType> TrimLower(
        const std::vector<GeometryPointerType>& rGeometries,
        const Point& rHighPoint,
        std::size_t Component)
    {
        using PointType = typename GeometryType::PointType;
        std::vector<GeometryPointerType> cut, trimmed_geometries;
        Vector nodal_distances;

        for (auto iter = rGeometries.begin(); iter < rGeometries.end();  ++iter) {
            auto& r_geom = **iter;
            cut.resize(0);

            std::size_t points_over = SignedDistanceOver(nodal_distances, r_geom, rHighPoint, Component);
            if (points_over == 3) {
                trimmed_geometries.push_back(*iter);
            }
            else if (points_over > 0) {
                auto splitter = DivideTriangle3D3<PointType>(r_geom, nodal_distances);
                splitter.GenerateDivision();
                cut = splitter.GetPositiveSubdivisions();
                trimmed_geometries.insert(trimmed_geometries.end(), cut.begin(), cut.end());
            }
        }

        return trimmed_geometries;
    }

    // Each component of rIntegrationPointData is a 4-component vector
    // Components 0-2 are the GLOBAL coordinates of the integration point
    // Component 3 is the integration point weight in terms of area (DetJ*Gauss_W)
    void AddTriangleIntegrationPointData(
        std::vector<Vector>& rIntegrationPointData,
        const SplitGeometryType& rGeom,
        const GeometryData::IntegrationMethod Method) const
    {
        const std::size_t num_nodes = rGeom.PointsNumber();
        const auto& r_integration_points = rGeom.IntegrationPoints(Method);
        const std::size_t num_gauss = r_integration_points.size();
        std::size_t new_point_offset = rIntegrationPointData.size();

        Matrix shape_functions = rGeom.ShapeFunctionsValues(Method);
        Vector det_j;
        rGeom.DeterminantOfJacobian(det_j, Method);

        for (std::size_t g = 0; g < num_gauss; g++) {
            Vector data = ZeroVector(4);
            data[3] = r_integration_points[g].Weight() * det_j[g];
            rIntegrationPointData.push_back(data);
        }

        for (std::size_t n = 0; n < num_nodes; n++) {
            const auto& r_coordinates = rGeom[n].Coordinates();
            for (std::size_t g = 0; g < num_gauss; g++) {
                for (std::size_t d = 0; d < 3; d++) {
                    rIntegrationPointData[new_point_offset+g][d] += shape_functions(g, n)*r_coordinates[d];
                }
            }
        }
    }

    bool PointIsOnSharedFace(Vector& rPointData, std::size_t i, std::size_t j, std::size_t k) const
    {
        constexpr double coordinate_tol = 1e-9;
        if (std::abs(rPointData[0] + 1.0) < coordinate_tol && i != 0) return true;
        if (std::abs(rPointData[0] - 1.0) < coordinate_tol && i != mNumDivisions[0]) return true;
        if (std::abs(rPointData[1] + 1.0) < coordinate_tol && j != 0) return true;
        if (std::abs(rPointData[1] - 1.0) < coordinate_tol && j != mNumDivisions[1]) return true;
        if (std::abs(rPointData[2] + 1.0) < coordinate_tol && k != 0) return true;
        if (std::abs(rPointData[2] - 1.0) < coordinate_tol && k != mNumDivisions[1]) return true;
        return false;
    }

};

}