//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "utilities/compare_elements_and_conditions_utility.h"

#include "includes/element.h"
#include "includes/condition.h"
#include "includes/kratos_components.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/quadrilateral_2d_4.h"

namespace Kratos {
namespace Testing {

class CustomElement : public Element
{
    public:
    CustomElement(Element::IndexType NewId, Element::GeometryType::Pointer pGeometry) :
        Element(NewId, pGeometry)
    {}
};

class CustomCondition : public Condition
{
    public:
    CustomCondition(Condition::IndexType NewId, Condition::GeometryType::Pointer pGeometry) :
        Condition(NewId, pGeometry)
    {}
};

Element::Pointer CreateCustomElement2D3N(Element::IndexType Id)
{
    Element::GeometryType::Pointer p_geom(new Triangle2D3<Element::NodeType>(Element::GeometryType::PointsArrayType(3)));
    return Element::Pointer(new CustomElement(Id, p_geom));
}

Element::Pointer CreateCustomElement2D4N(Element::IndexType Id)
{
    Element::GeometryType::Pointer p_geom(new Quadrilateral2D4<Element::NodeType>(Element::GeometryType::PointsArrayType(4)));
    return Element::Pointer(new CustomElement(Id, p_geom));
}

Condition::Pointer CreateCustomLineCondition2D3N(Condition::IndexType Id)
{
    Condition::GeometryType::Pointer p_geom(new Triangle2D3<Condition::NodeType>(Condition::GeometryType::PointsArrayType(3)));
    return Condition::Pointer(new CustomCondition(Id, p_geom));
}

Condition::Pointer CreateCustomCondition2D4N(Condition::IndexType Id)
{
    Condition::GeometryType::Pointer p_geom(new Quadrilateral2D4<Condition::NodeType>(Condition::GeometryType::PointsArrayType(4)));
    return Condition::Pointer(new CustomCondition(Id, p_geom));
}

KRATOS_TEST_CASE_IN_SUITE(IsSame_Element_Element, KratosCoreFastSuite)
{
    auto p_custom_element_2d_3n_1 = CreateCustomElement2D3N(1);
    const Element& r_custom_element_2d_3n_1 = *p_custom_element_2d_3n_1;
    auto p_custom_element_2d_3n_2 = CreateCustomElement2D3N(2);
    // same element type, same geometry type.
    KRATOS_CHECK(GeometricalObject::IsSame(r_custom_element_2d_3n_1, *p_custom_element_2d_3n_2));
    KRATOS_CHECK(GeometricalObject::IsSame(&r_custom_element_2d_3n_1, &(*p_custom_element_2d_3n_2)));
    KRATOS_CHECK(GeometricalObject::IsSame(*p_custom_element_2d_3n_2, r_custom_element_2d_3n_1));
    KRATOS_CHECK(GeometricalObject::IsSame(&(*p_custom_element_2d_3n_2), &r_custom_element_2d_3n_1));
    // different element type, same geometry type.
    const Element& r_element_2d_3n = KratosComponents<Element>::Get("Element2D3N");
    KRATOS_CHECK_IS_FALSE(GeometricalObject::IsSame(r_element_2d_3n, *p_custom_element_2d_3n_1));
    KRATOS_CHECK_IS_FALSE(GeometricalObject::IsSame(&r_element_2d_3n, &(*p_custom_element_2d_3n_1)));
    KRATOS_CHECK_IS_FALSE(GeometricalObject::IsSame(*p_custom_element_2d_3n_1, r_element_2d_3n));
    KRATOS_CHECK_IS_FALSE(GeometricalObject::IsSame(&(*p_custom_element_2d_3n_1), &r_element_2d_3n));
    // same element type, different geometry type.
    auto p_custom_element_2d_4n_3 = CreateCustomElement2D4N(3);
    const Element& r_custom_element_2d_4n_3 = *p_custom_element_2d_4n_3;
    KRATOS_CHECK_IS_FALSE(GeometricalObject::IsSame(r_custom_element_2d_3n_1, r_custom_element_2d_4n_3));
    KRATOS_CHECK_IS_FALSE(GeometricalObject::IsSame(&r_custom_element_2d_3n_1, &r_custom_element_2d_4n_3));
    KRATOS_CHECK_IS_FALSE(GeometricalObject::IsSame(r_custom_element_2d_4n_3, r_custom_element_2d_3n_1));
    KRATOS_CHECK_IS_FALSE(GeometricalObject::IsSame(&r_custom_element_2d_4n_3, &r_custom_element_2d_3n_1));
}

KRATOS_TEST_CASE_IN_SUITE(IsSame_Condition_Condition, KratosCoreFastSuite)
{
    auto p_custom_condition_2d_3n_1 = CreateCustomLineCondition2D3N(1);
    const Condition& r_custom_condition_2d_3n_1 = *p_custom_condition_2d_3n_1;
    auto p_custom_condition_2d_3n_2 = CreateCustomLineCondition2D3N(2);
    // same condition type, same geometry type.
    KRATOS_CHECK(GeometricalObject::IsSame(r_custom_condition_2d_3n_1, *p_custom_condition_2d_3n_2));
    KRATOS_CHECK(GeometricalObject::IsSame(&r_custom_condition_2d_3n_1, &(*p_custom_condition_2d_3n_2)));
    KRATOS_CHECK(GeometricalObject::IsSame(*p_custom_condition_2d_3n_2, r_custom_condition_2d_3n_1));
    KRATOS_CHECK(GeometricalObject::IsSame(&(*p_custom_condition_2d_3n_2), &r_custom_condition_2d_3n_1));
    // different condition type, same geometry type.
    const Condition& r_condition_2d_3n = KratosComponents<Condition>::Get("LineCondition2D3N");
    KRATOS_CHECK_IS_FALSE(GeometricalObject::IsSame(r_condition_2d_3n, *p_custom_condition_2d_3n_1));
    KRATOS_CHECK_IS_FALSE(GeometricalObject::IsSame(&r_condition_2d_3n, &(*p_custom_condition_2d_3n_1)));
    KRATOS_CHECK_IS_FALSE(GeometricalObject::IsSame(*p_custom_condition_2d_3n_1, r_condition_2d_3n));
    KRATOS_CHECK_IS_FALSE(GeometricalObject::IsSame(&(*p_custom_condition_2d_3n_1), &r_condition_2d_3n));
    // same condition type, different geometry type.
    auto p_custom_condition_2d_4n_3 = CreateCustomCondition2D4N(3);
    const Condition& r_custom_condition_2d_4n_3 = *p_custom_condition_2d_4n_3;
    KRATOS_CHECK_IS_FALSE(GeometricalObject::IsSame(r_custom_condition_2d_3n_1, r_custom_condition_2d_4n_3));
    KRATOS_CHECK_IS_FALSE(GeometricalObject::IsSame(&r_custom_condition_2d_3n_1, &r_custom_condition_2d_4n_3));
    KRATOS_CHECK_IS_FALSE(GeometricalObject::IsSame(r_custom_condition_2d_4n_3, r_custom_condition_2d_3n_1));
    KRATOS_CHECK_IS_FALSE(GeometricalObject::IsSame(&r_custom_condition_2d_4n_3, &r_custom_condition_2d_3n_1));
}

KRATOS_TEST_CASE_IN_SUITE(GetRegisteredNameElement, KratosCoreFastSuite)
{
    // Element name
    const Element& r_element_2d_3n = KratosComponents<Element>::Get("Element2D3N");
    std::string component_name;
    CompareElementsAndConditionsUtility::GetRegisteredName(r_element_2d_3n, component_name);
    KRATOS_CHECK_EQUAL(component_name, "Element2D3N");
    CompareElementsAndConditionsUtility::GetRegisteredName(&r_element_2d_3n, component_name);
    KRATOS_CHECK_EQUAL(component_name, "Element2D3N");
    CompareElementsAndConditionsUtility::GetRegisteredName(r_element_2d_3n, component_name);
    KRATOS_CHECK_NOT_EQUAL(component_name, "Element");
    CompareElementsAndConditionsUtility::GetRegisteredName(&r_element_2d_3n, component_name);
    KRATOS_CHECK_NOT_EQUAL(component_name, "Element");
}

KRATOS_TEST_CASE_IN_SUITE(GetRegisteredNameCondition, KratosCoreFastSuite)
{
    // Condition name
    const Condition& r_condition_2d_3n = KratosComponents<Condition>::Get("LineCondition2D3N");
    std::string component_name;
    CompareElementsAndConditionsUtility::GetRegisteredName(r_condition_2d_3n, component_name);
    KRATOS_CHECK_EQUAL(component_name, "LineCondition2D3N");
    CompareElementsAndConditionsUtility::GetRegisteredName(&r_condition_2d_3n, component_name);
    KRATOS_CHECK_EQUAL(component_name, "LineCondition2D3N");
    CompareElementsAndConditionsUtility::GetRegisteredName(r_condition_2d_3n, component_name);
    KRATOS_CHECK_NOT_EQUAL(component_name, "Condition");
    CompareElementsAndConditionsUtility::GetRegisteredName(&r_condition_2d_3n, component_name);
    KRATOS_CHECK_NOT_EQUAL(component_name, "Condition");
}

}  // namespace Testing.
}  // namespace Kratos.
