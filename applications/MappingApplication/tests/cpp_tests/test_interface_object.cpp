//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// Project includes
#include "geometries/quadrilateral_2d_4.h"

// Application includes
#include "custom_searching/interface_object.h"
#include "tests/cpp_tests/mapping_fast_suite.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometryObject, KratosMappingApplicationSerialTestSuite)
{
    array_1d<double, 3> init_coords;
    init_coords[0] = 0.0;
    init_coords[1] = 0.0;
    init_coords[2] = 0.0;

    InterfaceObject interface_obj(init_coords);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(interface_obj.pGetBaseNode(),
        "Error: Base class function called!");

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(interface_obj.pGetBaseGeometry(),
        "Error: Base class function called!");

    for (std::size_t i=0; i<3; ++i)
        KRATOS_EXPECT_DOUBLE_EQ(interface_obj.Coordinates()[i], 0.0);

    array_1d<double, 3> new_coords;
    new_coords[0] = -5.3;
    new_coords[1] = 18.7993;
    new_coords[2] = -547.1;

    interface_obj.Coordinates() = new_coords;

    for (std::size_t i=0; i<3; ++i)
        KRATOS_EXPECT_DOUBLE_EQ(interface_obj.Coordinates()[i], new_coords[i]);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceNode, KratosMappingApplicationSerialTestSuite)
{
    Point coords(1.0, 2.5, 30.0);

    const auto node_1(Kratos::make_shared<Node>(1, coords));

    Kratos::unique_ptr<InterfaceObject> p_interface_obj(Kratos::make_unique<InterfaceNode>(node_1.get()));

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_interface_obj->pGetBaseGeometry(),
        "Error: Base class function called!");

    for (std::size_t i=0; i<3; ++i)
        KRATOS_EXPECT_DOUBLE_EQ(p_interface_obj->Coordinates()[i], coords[i]);

    KRATOS_EXPECT_EQ(*(p_interface_obj->pGetBaseNode()), *node_1);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceObject, KratosMappingApplicationSerialTestSuite)
{
    Node::Pointer p_point1(new Node(1, 0.00, 0.00, 0.00));
    Node::Pointer p_point2(new Node(2, 0.00, 10.00, 0.00));
    Node::Pointer p_point3(new Node(3, 10.00, 10.00, 0.00));
    Node::Pointer p_point4(new Node(4, 10.00, 0.00, 0.00));

    const Kratos::shared_ptr<Geometry<Node>> p_quad(
        Kratos::make_shared<Quadrilateral2D4<Node>>(
            p_point1, p_point2, p_point3, p_point4));

    Kratos::unique_ptr<InterfaceObject> p_interface_obj(Kratos::make_unique<InterfaceGeometryObject>(p_quad.get()));

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_interface_obj->pGetBaseNode(),
        "Error: Base class function called!");

    Point center_coords(5.0, 5.0, 0.0);

    for (std::size_t i=0; i<3; ++i)
        KRATOS_EXPECT_DOUBLE_EQ(p_interface_obj->Coordinates()[i], center_coords[i]);

    // KRATOS_EXPECT_EQ(*(p_interface_obj->pGetBaseGeometry()), *p_quad); // does not compile ...
}

}  // namespace Kratos::Testing