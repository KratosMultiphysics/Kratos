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
//                   Pooyan Dadvand
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/line_3d_2.h"

#include "containers/geometry_container.h"

namespace Kratos::Testing {

namespace
{

Line3D2<Node>::Pointer GenerateLineGeometry() {
    Geometry<Node>::PointsArrayType nodes_vector;
    nodes_vector.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes_vector.push_back(Kratos::make_intrusive<Node>(2, 1.0, 1.0, 1.0));
    return Kratos::make_shared<Line3D2<Node>>(nodes_vector);
}

void FillGeometryContainer(
    Geometry<Node>::PointsArrayType& rAuxNodeList,
    GeometryContainer<Geometry<Node>>& rGeometryContainer)
{
    rAuxNodeList.push_back(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    rAuxNodeList.push_back(Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    rAuxNodeList.push_back(Kratos::make_intrusive<Node>(3, 0.0, 1.0, 0.0));
    rAuxNodeList.push_back(Kratos::make_intrusive<Node>(4, 1.0, 1.0, 0.0));

    Geometry<Node>::PointsArrayType nodes_vector_1;
    nodes_vector_1.push_back(rAuxNodeList(0));
    nodes_vector_1.push_back(rAuxNodeList(1));
    nodes_vector_1.push_back(rAuxNodeList(3));
    auto p_geom_1 = Kratos::make_shared<Triangle2D3<Node>>(1, nodes_vector_1);

    Geometry<Node>::PointsArrayType nodes_vector_2;
    nodes_vector_2.push_back(rAuxNodeList(0));
    nodes_vector_2.push_back(rAuxNodeList(3));
    nodes_vector_2.push_back(rAuxNodeList(2));
    auto p_geom_2 = Kratos::make_shared<Triangle2D3<Node>>(2, nodes_vector_2);

    rGeometryContainer.AddGeometry(p_geom_1);
    rGeometryContainer.AddGeometry(p_geom_2);
}

}

    ///// Test Geometry Container
    KRATOS_TEST_CASE_IN_SUITE(TestGeometryContainer, KratosCoreGeometryContainerFastSuite) {
        auto geometry_container = GeometryContainer<Geometry<Node>>();

        auto p_line_1 = GenerateLineGeometry();
        p_line_1->SetId(1);

        geometry_container.AddGeometry(p_line_1);
        KRATOS_EXPECT_EQ(geometry_container.NumberOfGeometries(), 1);
        geometry_container.AddGeometry(p_line_1); // adding same geomerty does not fail
        KRATOS_EXPECT_EQ(geometry_container.NumberOfGeometries(), 1);

        auto p_line_2 = GenerateLineGeometry();
        p_line_2->SetId(1);

        // check that we do nothing if the same geometry (same id and connectivities) is to be added multiple times
        KRATOS_EXPECT_EQ(geometry_container.NumberOfGeometries(), 1);

        p_line_2->SetId(2);
        geometry_container.AddGeometry(p_line_2);

        // check correct number of geometries
        KRATOS_EXPECT_EQ(geometry_container.NumberOfGeometries(), 2);

        // check adding with string
        auto p_line_3 = GenerateLineGeometry();
        p_line_3->SetId("GeometryLine1");
        geometry_container.AddGeometry(p_line_3);

        KRATOS_EXPECT_EQ(geometry_container.NumberOfGeometries(), 3);

        // check if correct element is returned
        KRATOS_EXPECT_EQ(geometry_container.GetGeometry(1).Id(), 1);
        KRATOS_EXPECT_EQ(geometry_container.pGetGeometry(1)->Id(), 1);

        // check remove functions
        geometry_container.RemoveGeometry("GeometryLine1");
        geometry_container.RemoveGeometry(1);
        KRATOS_EXPECT_EQ(geometry_container.NumberOfGeometries(), 1);
        KRATOS_EXPECT_FALSE(geometry_container.HasGeometry("GeometryLine1"));
    }

    KRATOS_TEST_CASE_IN_SUITE(TestGeometryContainerWithRepeatedGeometries1, KratosCoreGeometryContainerFastSuite)
    {
        // Prepare geometry container
        Geometry<Node>::PointsArrayType aux_node_pt_list;
        GeometryContainer<Geometry<Node>> geometry_container;
        FillGeometryContainer(aux_node_pt_list, geometry_container);

        // Create and add a geometry with the same connectivities but different Id (valid)
        Geometry<Node>::PointsArrayType nodes_vector_3;
        nodes_vector_3.push_back(aux_node_pt_list(0));
        nodes_vector_3.push_back(aux_node_pt_list(1));
        nodes_vector_3.push_back(aux_node_pt_list(3));
        auto p_geom_3 = Kratos::make_shared<Triangle2D3<Node>>(3, nodes_vector_3);
        geometry_container.AddGeometry(p_geom_3);

        // Check results
        KRATOS_EXPECT_EQ(geometry_container.NumberOfGeometries(), 3);
    }

    KRATOS_TEST_CASE_IN_SUITE(TestGeometryContainerWithRepeatedGeometries2, KratosCoreGeometryContainerFastSuite)
    {
        // Prepare geometry container
        Geometry<Node>::PointsArrayType aux_node_pt_list;
        GeometryContainer<Geometry<Node>> geometry_container;
        FillGeometryContainer(aux_node_pt_list, geometry_container);

        // Keep the existing when trying to add an already existing geometry (same Id and same connectivities) (valid)
        Geometry<Node>::PointsArrayType nodes_vector_2;
        nodes_vector_2.push_back(aux_node_pt_list(0));
        nodes_vector_2.push_back(aux_node_pt_list(3));
        nodes_vector_2.push_back(aux_node_pt_list(2));
        auto p_geom_2 = Kratos::make_shared<Triangle2D3<Node>>(2, nodes_vector_2);
        geometry_container.AddGeometry(p_geom_2);

        // Check results
        KRATOS_EXPECT_EQ(geometry_container.NumberOfGeometries(), 2);
    }

    KRATOS_TEST_CASE_IN_SUITE(TestGeometryContainerWithRepeatedGeometries3, KratosCoreGeometryContainerFastSuite)
    {
        // Prepare geometry container
        Geometry<Node>::PointsArrayType aux_node_pt_list;
        GeometryContainer<Geometry<Node>> geometry_container;
        FillGeometryContainer(aux_node_pt_list, geometry_container);

        // Throw an error when trying to add a different geometry with an existing Id (not defined)
        Geometry<Node>::PointsArrayType nodes_vector_2;
        nodes_vector_2.push_back(aux_node_pt_list(0));
        nodes_vector_2.push_back(aux_node_pt_list(1));
        nodes_vector_2.push_back(aux_node_pt_list(3));
        nodes_vector_2.push_back(aux_node_pt_list(2));
        auto p_geom_2 = Kratos::make_shared<Quadrilateral2D4<Node>>(2, nodes_vector_2);
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry_container.AddGeometry(p_geom_2), "Attempting to add geometry with Id: 2. A different geometry with the same Id already exists.")
    }

    KRATOS_TEST_CASE_IN_SUITE(TestGeometryContainerWithRepeatedGeometries4, KratosCoreGeometryContainerFastSuite)
    {
        // Prepare geometry container
        Geometry<Node>::PointsArrayType aux_node_pt_list;
        GeometryContainer<Geometry<Node>> geometry_container;
        FillGeometryContainer(aux_node_pt_list, geometry_container);

        // Throw an error when trying to add a same type geometry but with different connectivities and existing Id (not defined)
        Geometry<Node>::PointsArrayType nodes_vector_2;
        nodes_vector_2.push_back(aux_node_pt_list(2));
        nodes_vector_2.push_back(aux_node_pt_list(3));
        nodes_vector_2.push_back(aux_node_pt_list(0));
        auto p_geom_2 = Kratos::make_shared<Triangle2D3<Node>>(2, nodes_vector_2);
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geometry_container.AddGeometry(p_geom_2), "Attempting to add a new geometry with Id :2. A same type geometry with same Id but different connectivities already exists.")
    }

} // namespace Kratos::Testing.
