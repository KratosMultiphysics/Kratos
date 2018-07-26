//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// Project includes
#include "testing/testing.h"
#include "geometries/quadrilateral_2d_4.h"
#include "custom_searching/interface_object.h"

namespace Kratos {
namespace Testing {


KRATOS_TEST_CASE_IN_SUITE(InterfaceGeometryObject, KratosMappingApplicationSerialTestSuite)
{
    InterfaceObject interface_obj(array_1d<double, 3>(0.0));

    KRATOS_CHECK_EXCEPTION_IS_THROWN(interface_obj.UpdateCoordinates(),
        "Error: Base class function called!");

    KRATOS_CHECK_EXCEPTION_IS_THROWN(interface_obj.pGetBaseNode(),
        "Error: Base class function called!");

    KRATOS_CHECK_EXCEPTION_IS_THROWN(interface_obj.pGetBaseGeometry(),
        "Error: Base class function called!");

    for (std::size_t i=0; i<3; ++i)
        KRATOS_CHECK_DOUBLE_EQUAL(interface_obj.Coordinates()[i], 0.0);

    array_1d<double, 3> new_coords;
    new_coords[0] = -5.3;
    new_coords[1] = 18.7993;
    new_coords[2] = -547.1;

    interface_obj.UpdateCoordinates(new_coords);

    for (std::size_t i=0; i<3; ++i)
        KRATOS_CHECK_DOUBLE_EQUAL(interface_obj.Coordinates()[i], new_coords[i]);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceNode, KratosMappingApplicationSerialTestSuite)
{
    Point coords(1.0, 2.5, 30.0);

    const auto node_1(Kratos::make_shared<Node<3>>(1, coords));

    Kratos::unique_ptr<InterfaceObject> p_interface_obj(Kratos::make_unique<InterfaceNode>(node_1));

    KRATOS_CHECK_EXCEPTION_IS_THROWN(p_interface_obj->pGetBaseGeometry(),
        "Error: Base class function called!");

    for (std::size_t i=0; i<3; ++i)
        KRATOS_CHECK_DOUBLE_EQUAL(p_interface_obj->Coordinates()[i], coords[i]);

    Point new_coords(18.3, -89.123, 125.7);

    noalias(node_1->Coordinates()) = new_coords;

    p_interface_obj->UpdateCoordinates();

    for (std::size_t i=0; i<3; ++i)
        KRATOS_CHECK_DOUBLE_EQUAL(p_interface_obj->Coordinates()[i], new_coords[i]);

    KRATOS_CHECK_EQUAL(*(p_interface_obj->pGetBaseNode()), *node_1);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceObject, KratosMappingApplicationSerialTestSuite)
{
    Node<3>::Pointer p_point1(new Node<3>(1, 0.00, 0.00, 0.00));
    Node<3>::Pointer p_point2(new Node<3>(2, 0.00, 10.00, 0.00));
    Node<3>::Pointer p_point3(new Node<3>(3, 10.00, 10.00, 0.00));
    Node<3>::Pointer p_point4(new Node<3>(4, 10.00, 0.00, 0.00));

    const Kratos::shared_ptr<Geometry<Node<3>>> p_quad(
        Kratos::make_shared<Quadrilateral2D4<Node<3>>>(
            p_point1, p_point2, p_point3, p_point4));

    Kratos::unique_ptr<InterfaceObject> p_interface_obj(Kratos::make_unique<InterfaceGeometryObject>(p_quad));

    KRATOS_CHECK_EXCEPTION_IS_THROWN(p_interface_obj->pGetBaseNode(),
        "Error: Base class function called!");

    Point center_coords(5.0, 5.0, 0.0);

    for (std::size_t i=0; i<3; ++i)
        KRATOS_CHECK_DOUBLE_EQUAL(p_interface_obj->Coordinates()[i], center_coords[i]);

    p_point1->X() -= 2.5;
    p_point2->X() -= 2.5;
    p_point3->X() -= 2.5;
    p_point4->X() -= 2.5;

    p_point1->Y() += 2.0;
    p_point4->Y() += 2.0;

    p_point1->Z() = 2.0;
    p_point2->Z() = 2.0;
    p_point3->Z() = 2.0;
    p_point4->Z() = 2.0;

    Point new_center_coords(2.5, 6.0, 2.0);

    p_interface_obj->UpdateCoordinates();

    for (std::size_t i=0; i<3; ++i)
        KRATOS_CHECK_DOUBLE_EQUAL(p_interface_obj->Coordinates()[i], new_center_coords[i]);

    // KRATOS_CHECK_EQUAL(*(p_interface_obj->pGetBaseGeometry()), *p_quad); // does not compile ...
}

}  // namespace Testing
}  // namespace Kratos