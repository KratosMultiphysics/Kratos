//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_volume_geometry.h"
#include "geometries/brep_surface_on_volume.h"
#include "geometries/brep_volume.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    typedef Node NodeType;

    NurbsVolumeGeometry<PointerVector<NodeType>>::Pointer GenerateReferenceVolume() {
        // Geometry<NodeType>::PointsArrayType points;

        // points.push_back(NodeType::Pointer(new NodeType(1, 10, 0, 0)));
        // points.push_back(NodeType::Pointer(new NodeType(2, 10, 0, 10)));
        // points.push_back(NodeType::Pointer(new NodeType(3, 0, 0, 10)));
        // points.push_back(NodeType::Pointer(new NodeType(4, -10, 0, 10)));
        // points.push_back(NodeType::Pointer(new NodeType(5, -10, 0, 0)));
        // points.push_back(NodeType::Pointer(new NodeType(6, 10, 10, 0)));
        // points.push_back(NodeType::Pointer(new NodeType(7, 10, 10, 10)));
        // points.push_back(NodeType::Pointer(new NodeType(8, 0, 10, 10)));
        // points.push_back(NodeType::Pointer(new NodeType(9, -10, 10, 10)));
        // points.push_back(NodeType::Pointer(new NodeType(10, -10, 10, 0)));

        // Vector knot_u = ZeroVector(6);
        // knot_u[0] = 0.0;
        // knot_u[1] = 0.0;
        // knot_u[2] = 15.707963267948966;
        // knot_u[3] = 15.707963267948966;
        // knot_u[4] = 31.415926535897931;
        // knot_u[5] = 31.415926535897931;
        // Vector knot_v = ZeroVector(2);
        // knot_v[0] = 0.0;
        // knot_v[1] = 10.0;

        // int p = 2;
        // int q = 1;

        // Vector weights = ZeroVector(10);
        // weights[0] = 1.0;
        // weights[1] = 0.70710678118654757;
        // weights[2] = 1.0;
        // weights[3] = 0.70710678118654757;
        // weights[4] = 1.0;
        // weights[5] = 1.0;
        // weights[6] = 0.70710678118654757;
        // weights[7] = 1.0;
        // weights[8] = 0.70710678118654757;
        // weights[9] = 1.0;

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

        // Define knot vectors (uniform, open)
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

    // NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateReference1Curve2dPointer()
    // {
    //     NurbsCurveGeometry<2, PointerVector<Point>>::PointsArrayType points;

    //     points.push_back(Point::Pointer(new Point(0, 0)));
    //     points.push_back(Point::Pointer(new Point(31.415926535897931, 0)));

    //     Vector knot_vector = ZeroVector(2);
    //     knot_vector[0] = 0.0;
    //     knot_vector[1] = 31.415926535897931;

    //     int p = 1;

    //     return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    // }

    // NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateReference2Curve2dPointer()
    // {
    //     NurbsCurveGeometry<2, PointerVector<Point>>::PointsArrayType points;

    //     points.push_back(Point::Pointer(new Point(31.415926535897931, 0)));
    //     points.push_back(Point::Pointer(new Point(31.415926535897931, 10)));

    //     Vector knot_vector = ZeroVector(2);
    //     knot_vector[0] = 0.0;
    //     knot_vector[1] = 10.0;

    //     int p = 1;

    //     return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    // }

    // NurbsCurveGeometry<2, PointerVector<Point>>::Pointer GenerateReference3Curve2dPointer()
    // {
    //     NurbsCurveGeometry<2, PointerVector<Point>>::PointsArrayType points;

    //     points.push_back(Point::Pointer(new Point(31.415926535897931, 10)));
    //     points.push_back(Point::Pointer(new Point(28.500619950938574, 10)));
    //     points.push_back(Point::Pointer(new Point(26.033141879879658, 9.48077148945451)));
    //     points.push_back(Point::Pointer(new Point(23.323541783172974, 8.4455316972070165)));
    //     points.push_back(Point::Pointer(new Point(20.805481541235128, 7.4834723273026622)));
    //     points.push_back(Point::Pointer(new Point(18.157204808069945, 6.1025460603903392)));
    //     points.push_back(Point::Pointer(new Point(15.707963267948966, 5)));
    //     points.push_back(Point::Pointer(new Point(13.258721727828151, 3.8974539396097354)));
    //     points.push_back(Point::Pointer(new Point(10.610444994662638, 2.5165276726972694)));
    //     points.push_back(Point::Pointer(new Point(8.09238475272496,   1.5544683027929818)));
    //     points.push_back(Point::Pointer(new Point(5.382784656018293,  0.51922851054549679)));
    //     points.push_back(Point::Pointer(new Point(2.91530658495934,   0)));
    //     points.push_back(Point::Pointer(new Point(0, 0)));

    //     Vector knot_vector = ZeroVector(15);
    //     knot_vector[0] =  -34.16313447529194;
    //     knot_vector[1] =  -34.16313447529194;
    //     knot_vector[2] =  -34.16313447529194;
    //     knot_vector[3] =  -25.309419191442988;
    //     knot_vector[4] =  -25.309419191442988;
    //     knot_vector[5] =  -25.309419191442988;
    //     knot_vector[6] =  -17.081567237645665;
    //     knot_vector[7] =  -17.081567237645665;
    //     knot_vector[8] =  -17.081567237645665;
    //     knot_vector[9] =  -8.8537152838488975;
    //     knot_vector[10] = -8.8537152838488975;
    //     knot_vector[11] = -8.8537152838488975;
    //     knot_vector[12] = 0;
    //     knot_vector[13] = 0;
    //     knot_vector[14] = 0;

    //     int p = 3;

    //     return Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(points, p, knot_vector);
    // }

    KRATOS_TEST_CASE_IN_SUITE(NurbsBrepVolumeSurrogate, KratosCoreGeometriesFastSuite)
    {
        auto p_volume = GenerateReferenceVolume();

        using BrepVolumeType = BrepVolume<PointerVector<Node>, true, PointerVector<Point>>;
        using BrepSurfaceOnVolumeLoopArrayType = typename BrepVolumeType::BrepSurfaceOnVolumeLoopArrayType;

        using GeometrySurrogateArrayType = typename BrepVolumeType::GeometrySurrogateArrayType;
        

        BrepSurfaceOnVolumeLoopArrayType outer_loops, inner_loops;
        
        Model model;
        ModelPart& rSurrogateModelPartOuter = model.CreateModelPart("surrogate_model_part_outer");
        rSurrogateModelPartOuter.CreateNewProperties(0);
        rSurrogateModelPartOuter.CreateNewNode(1, 0.0, 0.0, 0.0);
        rSurrogateModelPartOuter.CreateNewNode(2, 2.0, 0.0, 0.0);
        rSurrogateModelPartOuter.CreateNewNode(3, 2.0, 2.0, 0.0);
        rSurrogateModelPartOuter.CreateNewNode(4, 0.0, 2.0, 0.0);

        Properties::Pointer p_prop = rSurrogateModelPartOuter.pGetProperties(0);
        rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", 1, {{1, 2}}, p_prop);
        rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", 2, {{2, 3}}, p_prop);
        rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", 3, {{3, 4}}, p_prop);
        rSurrogateModelPartOuter.CreateNewCondition("LineCondition2D2N", 4, {{4, 1}}, p_prop);

        GeometrySurrogateArrayType surrogate_outer_loop_geometries(rSurrogateModelPartOuter.NumberOfConditions());
        GeometrySurrogateArrayType surrogate_inner_loop_geometries(rSurrogateModelPartOuter.NumberOfConditions());

        int count = 0;
        for (auto i_cond : rSurrogateModelPartOuter.Conditions())
        {
            surrogate_outer_loop_geometries[count] = i_cond.pGetGeometry();
            count++;
        }

        count = 0;
        for (auto i_cond : rSurrogateModelPartOuter.Conditions())
        {
            surrogate_inner_loop_geometries[count] = i_cond.pGetGeometry();
            count++;
        }

        auto p_brep_volume =
            Kratos::make_shared<BrepVolumeType>(
                p_volume, 
                outer_loops,
                inner_loops);

        auto p_outer = Kratos::make_shared<GeometrySurrogateArrayType>(surrogate_outer_loop_geometries);
        auto p_inner = Kratos::make_shared<GeometrySurrogateArrayType>(surrogate_inner_loop_geometries);

        p_brep_volume->SetSurrogateOuterLoopGeometries(p_outer);
        p_brep_volume->SetSurrogateInnerLoopGeometries(p_inner);

        KRATOS_EXPECT_EQ(p_brep_volume->GetSurrogateInnerLoopGeometries().size(), 4);
        KRATOS_EXPECT_EQ(p_brep_volume->GetSurrogateOuterLoopGeometries().size(), 4);

    }
} // namespace Testing.
} // namespace Kratos.
