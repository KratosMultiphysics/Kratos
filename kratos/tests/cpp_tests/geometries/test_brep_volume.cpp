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
        Properties::Pointer p_prop = rSurrogateModelPartOuter.pGetProperties(0);

        rSurrogateModelPartOuter.CreateNewNode(1, 0.0, 0.0, 0.0);
        rSurrogateModelPartOuter.CreateNewNode(2, 2.0, 0.0, 0.0);
        rSurrogateModelPartOuter.CreateNewNode(3, 2.0, 2.0, 0.0);
        rSurrogateModelPartOuter.CreateNewNode(4, 0.0, 2.0, 0.0);
        rSurrogateModelPartOuter.CreateNewNode(5, 0.0, 0.0, 2.0);
        rSurrogateModelPartOuter.CreateNewNode(6, 2.0, 0.0, 2.0);
        rSurrogateModelPartOuter.CreateNewNode(7, 2.0, 2.0, 2.0);
        rSurrogateModelPartOuter.CreateNewNode(8, 0.0, 2.0, 2.0);

        // (XY, z=0)
        rSurrogateModelPartOuter.CreateNewCondition("SurfaceCondition3D4N", 1, {{1, 2, 3, 4}}, p_prop);

        // (XY, z=2)
        rSurrogateModelPartOuter.CreateNewCondition("SurfaceCondition3D4N", 2, {{5, 6, 7, 8}}, p_prop);

        // (YZ, x=0)
        rSurrogateModelPartOuter.CreateNewCondition("SurfaceCondition3D4N", 3, {{1, 4, 8, 5}}, p_prop);

        // (YZ, x=2)
        rSurrogateModelPartOuter.CreateNewCondition("SurfaceCondition3D4N", 4, {{2, 6, 7, 3}}, p_prop);

        // (XZ, y=0)
        rSurrogateModelPartOuter.CreateNewCondition("SurfaceCondition3D4N", 5, {{1, 2, 6, 5}}, p_prop);

        // (XZ, y=2)
        rSurrogateModelPartOuter.CreateNewCondition("SurfaceCondition3D4N", 6, {{4, 3, 7, 8}}, p_prop);


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
        
        KRATOS_EXPECT_EQ(p_brep_volume->GetSurrogateInnerLoopGeometries().size(), 6);
        KRATOS_EXPECT_EQ(p_brep_volume->GetSurrogateOuterLoopGeometries().size(), 6);

    }
} // namespace Testing.
} // namespace Kratos.
