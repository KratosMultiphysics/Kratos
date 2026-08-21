//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_volume_geometry.h"
#include "geometries/brep_surface.h"
#include "geometries/brep_surface_on_volume.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    typedef Node NodeType;

    NurbsVolumeGeometry<PointerVector<NodeType>>::Pointer GenerateReferenceVolumeForBrep() {

        Geometry<NodeType>::PointsArrayType points;

        // Create 3x3x3 grid of points for the volume
        int id = 1;
        for (std::size_t k = 0; k < 3; ++k) {
            for (std::size_t j = 0; j < 3; ++j) {
                for (std::size_t i = 0; i < 3; ++i) {
                    points.push_back(NodeType::Pointer(new NodeType(id++, i, j, k)));
                }
            }
        }

        // Define degree in each direction
        int p = 1;
        int q = 1;
        int r = 1;

        // KnotsURefined: [3](0, 1, 2)
        Vector knot_u = ZeroVector(3);
        knot_u[0] = 0.0;
        knot_u[1] = 1.0;
        knot_u[2] = 2.0;

        // KnotsVRefined: [3](0, 1, 2)
        Vector knot_v = ZeroVector(3);
        knot_v[0] = 0.0;
        knot_v[1] = 1.0;
        knot_v[2] = 2.0;

        // KnotsWRefined: [3](0, 1, 2)
        Vector knot_w = ZeroVector(3);
        knot_w[0] = 0.0;
        knot_w[1] = 1.0;
        knot_w[2] = 2.0;

        return Kratos::make_shared <NurbsVolumeGeometry<PointerVector<NodeType>>>(points, p, q, r, knot_u, knot_v, knot_w);
    }


    ///// Tests
    KRATOS_TEST_CASE_IN_SUITE(BrepSurfaceOnVolume, KratosCoreGeometriesFastSuite) {

        auto p_volume = GenerateReferenceVolumeForBrep();

        Vector knot_vector_u = ZeroVector(2); 
        Vector knot_vector_v = ZeroVector(2); 
        Vector normal = ZeroVector(3);
        const SizeType p = 1;
        const SizeType q = 1;
        Point rCoordsA(0.0, 0.0, 0.0);
        Point rCoordsB(2.0, 2.0, 2.0);
        IndexType rLastGeometryId = 1;

        // lower
        Geometry<NodeType>::PointsArrayType points_lower;
        knot_vector_u[0] = 0.0; knot_vector_u[1] = std::abs(rCoordsA[0]-rCoordsB[0]);  
        knot_vector_v[0] = 0.0; knot_vector_v[1] = std::abs(rCoordsA[1]-rCoordsB[1]); 
        points_lower.push_back(NodeType::Pointer(new NodeType(1, rCoordsA[0], rCoordsA[1], rCoordsA[2])));
        points_lower.push_back(NodeType::Pointer(new NodeType(2, rCoordsA[0], rCoordsB[1], rCoordsA[2])));
        points_lower.push_back(NodeType::Pointer(new NodeType(3, rCoordsB[0], rCoordsA[1], rCoordsA[2])));
        points_lower.push_back(NodeType::Pointer(new NodeType(4, rCoordsB[0], rCoordsB[1], rCoordsA[2])));
        normal = ZeroVector(3);
        normal[2] = -1;
        auto p_surface_1 = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(points_lower, p, q, knot_vector_u, knot_vector_v);
        auto brep_p_surface_1 = Kratos::make_shared<BrepSurface<PointerVector<NodeType>,  false, PointerVector<Point>>>(p_surface_1);
        auto p_brep_surface_on_volume_1 = Kratos::make_shared<BrepSurfaceOnVolume< PointerVector<NodeType>, false, PointerVector<NodeType>>>(p_volume, p_surface_1);
        p_brep_surface_on_volume_1->SetId(rLastGeometryId);
        p_brep_surface_on_volume_1->SetNormalSbm(normal);

        // back
        Geometry<NodeType>::PointsArrayType points_back;
        knot_vector_u[0] = 0.0; knot_vector_u[1] = std::abs(rCoordsA[0]-rCoordsB[0]);  
        knot_vector_v[0] = 0.0; knot_vector_v[1] = std::abs(rCoordsA[2]-rCoordsB[2]); 
        points_back.push_back(NodeType::Pointer(new NodeType(1, rCoordsA[0], rCoordsB[1], rCoordsA[2])));
        points_back.push_back(NodeType::Pointer(new NodeType(2, rCoordsA[0], rCoordsB[1], rCoordsB[2])));
        points_back.push_back(NodeType::Pointer(new NodeType(3, rCoordsB[0], rCoordsB[1], rCoordsA[2])));
        points_back.push_back(NodeType::Pointer(new NodeType(4, rCoordsB[0], rCoordsB[1], rCoordsB[2])));
        normal = ZeroVector(3);
        normal[1] = 1;
        auto p_surface_4 = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(points_back, p, q, knot_vector_u, knot_vector_v);
        auto brep_p_surface_4 = Kratos::make_shared<BrepSurface<PointerVector<NodeType>,  false, PointerVector<Point>>>(p_surface_4);
        auto p_brep_surface_on_volume_4 = Kratos::make_shared<BrepSurfaceOnVolume< PointerVector<NodeType>, false, PointerVector<NodeType>>>(p_volume, p_surface_4);
        p_brep_surface_on_volume_4->SetId(++rLastGeometryId);
        p_brep_surface_on_volume_4->SetNormalSbm(normal);

        auto p_nurbs_surface_on_volume_1 = p_brep_surface_on_volume_1->pGetSurfaceOnVolume();
        auto p_nurbs_surface_on_volume_4 = p_brep_surface_on_volume_4->pGetSurfaceOnVolume();

        // Check general information, input to ouput
        KRATOS_EXPECT_EQ(p_brep_surface_on_volume_1->WorkingSpaceDimension(), 3);
        KRATOS_EXPECT_EQ(p_brep_surface_on_volume_1->LocalSpaceDimension(), 1);
        KRATOS_EXPECT_EQ(p_brep_surface_on_volume_1->Center().X(), 1.0);
        KRATOS_EXPECT_EQ(p_brep_surface_on_volume_1->Center().Y(), 1.0);
        KRATOS_EXPECT_EQ(p_brep_surface_on_volume_1->Center().Z(), 0.0);
        KRATOS_EXPECT_EQ(p_nurbs_surface_on_volume_1->Center().X(), 1.0);
        KRATOS_EXPECT_EQ(p_nurbs_surface_on_volume_1->Center().Y(), 1.0);
        KRATOS_EXPECT_EQ(p_nurbs_surface_on_volume_1->Center().Z(), 0.0);

        KRATOS_EXPECT_EQ(p_brep_surface_on_volume_4->LocalSpaceDimension(), 1);
        KRATOS_EXPECT_EQ(p_brep_surface_on_volume_4->Center().X(), 1.0);
        KRATOS_EXPECT_EQ(p_brep_surface_on_volume_4->Center().Y(), 2.0);
        KRATOS_EXPECT_EQ(p_brep_surface_on_volume_4->Center().Z(), 1.0);
        KRATOS_EXPECT_EQ(p_nurbs_surface_on_volume_4->Center().X(), 1.0);
        KRATOS_EXPECT_EQ(p_nurbs_surface_on_volume_4->Center().Y(), 2.0);
        KRATOS_EXPECT_EQ(p_nurbs_surface_on_volume_4->Center().Z(), 1.0);

        const auto geometry_family = GeometryData::KratosGeometryFamily::Kratos_Brep;
        const auto geometry_type = GeometryData::KratosGeometryType::Kratos_Brep_Surface_On_Volume;
        KRATOS_EXPECT_EQ(p_brep_surface_on_volume_1->GetGeometryFamily(), geometry_family);
        KRATOS_EXPECT_EQ(p_brep_surface_on_volume_1->GetGeometryType(), geometry_type);
    }

} // namespace Testing.
} // namespace Kratos.
