//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application

// System includes

// External includes

// Project includes
#include "iga_application.h"
#include "iga_application_variables.h"

namespace Kratos {

KratosIgaApplication::KratosIgaApplication()
    : KratosApplication("IgaApplication")
    , mShell3pElement(0, Element::GeometryType::Pointer(
        new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1))))
    , mIgaTrussElement(0, Element::GeometryType::Pointer(
        new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1))))
    , mShellKLDiscreteElement(0, Element::GeometryType::Pointer(
        new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1))))
    , mLoadCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node<3>>(Condition::GeometryType::PointsArrayType(1))))
    , mPenaltyCouplingCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node<3>>(Condition::GeometryType::PointsArrayType(1))))
{
}

void KratosIgaApplication::Register() {
    KratosApplication::Register();
    KRATOS_INFO("") << "Initializing KratosIgaApplication..." << std::endl;

    // ELEMENTS
    KRATOS_REGISTER_ELEMENT("Shell3pElement", mShell3pElement)
    KRATOS_REGISTER_ELEMENT("IgaTrussElement", mIgaTrussElement)
    KRATOS_REGISTER_ELEMENT("ShellKLDiscreteElement", mShellKLDiscreteElement)

    // CONDITIONS
    KRATOS_REGISTER_CONDITION("LoadCondition", mLoadCondition)
    KRATOS_REGISTER_CONDITION("PenaltyCouplingCondition", mPenaltyCouplingCondition)

    // VARIABLES
    KRATOS_REGISTER_VARIABLE(NURBS_CONTROL_POINT_WEIGHT)

    KRATOS_REGISTER_VARIABLE(COORDINATES)
    KRATOS_REGISTER_VARIABLE(TANGENTS)

    KRATOS_REGISTER_VARIABLE(CROSS_AREA)
    KRATOS_REGISTER_VARIABLE(PRESTRESS_CAUCHY)

    KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_VALUES)
    KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_LOCAL_DERIVATIVES)
    KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES)

    KRATOS_REGISTER_VARIABLE(RAYLEIGH_ALPHA)
    KRATOS_REGISTER_VARIABLE(RAYLEIGH_BETA)

    KRATOS_REGISTER_VARIABLE(POINT_LOAD)
    KRATOS_REGISTER_VARIABLE(LINE_LOAD)
    KRATOS_REGISTER_VARIABLE(SURFACE_LOAD)

    KRATOS_REGISTER_VARIABLE(PENALTY_FACTOR)
}

}  // namespace Kratos
