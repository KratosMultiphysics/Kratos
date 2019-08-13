/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

// System includes

// External includes

// Project includes
#include "iga_application.h"
#include "iga_application_variables.h"

namespace Kratos {

KratosIgaApplication::KratosIgaApplication()
    : KratosApplication("IgaApplication")
    , mIgaTrussElement(0, Element::GeometryType::Pointer(
        new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1))))
    , mShellKLDiscreteElement(0, Element::GeometryType::Pointer(
        new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1))))
    , mIgaShell3pElement(0, Element::GeometryType::Pointer(
        new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1))))
    , mIgaShell5pElement(0, Element::GeometryType::Pointer(
        new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1))))
    , mIgaShell5pElementPreInt(0, Element::GeometryType::Pointer(
        new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1))))
    , mIgaShell5pElementStuttgart(0, Element::GeometryType::Pointer(
        new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1))))
    , mIgaShell7pElement(0, Element::GeometryType::Pointer(
        new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1))))
        
    , mCouplingPenaltyDiscreteCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node<3>>(Condition::GeometryType::PointsArrayType(1))))
    , mSupportPenaltyCurveDiscreteCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node<3>>(Condition::GeometryType::PointsArrayType(1))))
    , mLoadSurfaceDiscreteCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node<3>>(Condition::GeometryType::PointsArrayType(1))))
    , mLoadCurveDiscreteCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node<3>>(Condition::GeometryType::PointsArrayType(1))))
{
}

void KratosIgaApplication::Register() {
    KratosApplication::Register();
    KRATOS_INFO("") << "Initializing KratosIgaApplication..." << std::endl;

    // ELEMENTS
    KRATOS_REGISTER_ELEMENT("IgaTrussElement", mIgaTrussElement)
    KRATOS_REGISTER_ELEMENT("ShellKLDiscreteElement", mShellKLDiscreteElement)
    KRATOS_REGISTER_ELEMENT("IgaShell3pElement", mIgaShell3pElement)
    KRATOS_REGISTER_ELEMENT("IgaShell5pElement", mIgaShell5pElement)
    KRATOS_REGISTER_ELEMENT("IgaShell5pElementPreInt", mIgaShell5pElementPreInt)
    KRATOS_REGISTER_ELEMENT("IgaShell5pElementStuttgart", mIgaShell5pElementStuttgart)
    KRATOS_REGISTER_ELEMENT("IgaShell7pElement", mIgaShell7pElement)

    KRATOS_REGISTER_CONDITION("CouplingPenaltyDiscreteCondition", mCouplingPenaltyDiscreteCondition)
    KRATOS_REGISTER_CONDITION("SupportPenaltyCurveDiscreteCondition", mSupportPenaltyCurveDiscreteCondition)
    KRATOS_REGISTER_CONDITION("LoadSurfaceDiscreteCondition", mLoadSurfaceDiscreteCondition)
    KRATOS_REGISTER_CONDITION("LoadCurveDiscreteCondition", mLoadCurveDiscreteCondition)


    // Variables
    KRATOS_REGISTER_VARIABLE(BREP_ID)

    KRATOS_REGISTER_VARIABLE(NURBS_CONTROL_POINT_WEIGHT)

    KRATOS_REGISTER_VARIABLE(COORDINATES)
    KRATOS_REGISTER_VARIABLE(LOCAL_COORDINATES)
    KRATOS_REGISTER_VARIABLE(TANGENTS)
    KRATOS_REGISTER_VARIABLE(TANGENTS_SLAVE)

    KRATOS_REGISTER_VARIABLE(CROSS_AREA)
    KRATOS_REGISTER_VARIABLE(PRESTRESS_CAUCHY)

    KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_VALUES)
    KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_LOCAL_DERIVATIVES)
    KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES)
    KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_LOCAL_THIRD_DERIVATIVES)

    KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_VALUES_SLAVE)
    KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE)
    KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES_SLAVE)

    KRATOS_REGISTER_VARIABLE(RAYLEIGH_ALPHA)
    KRATOS_REGISTER_VARIABLE(RAYLEIGH_BETA)

    KRATOS_REGISTER_VARIABLE(POINT_LOAD)
    KRATOS_REGISTER_VARIABLE(LINE_LOAD)
    KRATOS_REGISTER_VARIABLE(SURFACE_LOAD)

    KRATOS_REGISTER_VARIABLE(PENALTY_FACTOR)

    KRATOS_REGISTER_VARIABLE(STRESS_CAUCHY_11)
    KRATOS_REGISTER_VARIABLE(STRESS_CAUCHY_12)
    KRATOS_REGISTER_VARIABLE(STRESS_CAUCHY_22)
    KRATOS_REGISTER_VARIABLE(STRESS_CAUCHY_23)
    KRATOS_REGISTER_VARIABLE(STRESS_CAUCHY_13)
    KRATOS_REGISTER_VARIABLE(STRESS_CAUCHY_TOP_11)
    KRATOS_REGISTER_VARIABLE(STRESS_CAUCHY_TOP_12)
    KRATOS_REGISTER_VARIABLE(STRESS_CAUCHY_TOP_22)
    KRATOS_REGISTER_VARIABLE(STRESS_CAUCHY_BOTTOM_11)
    KRATOS_REGISTER_VARIABLE(STRESS_CAUCHY_BOTTOM_12)
    KRATOS_REGISTER_VARIABLE(STRESS_CAUCHY_BOTTOM_22)
    KRATOS_REGISTER_VARIABLE(INTERNAL_FORCE_11)
    KRATOS_REGISTER_VARIABLE(INTERNAL_FORCE_12)
    KRATOS_REGISTER_VARIABLE(INTERNAL_FORCE_22)
    KRATOS_REGISTER_VARIABLE(INTERNAL_MOMENT_11)
    KRATOS_REGISTER_VARIABLE(INTERNAL_MOMENT_12)
    KRATOS_REGISTER_VARIABLE(INTERNAL_MOMENT_22)
    KRATOS_REGISTER_VARIABLE(SHEAR_FORCE_1)
    KRATOS_REGISTER_VARIABLE(SHEAR_FORCE_2)

    // DOF
    KRATOS_REGISTER_VARIABLE(W_BAR)
    KRATOS_REGISTER_VARIABLE(W_BAR_REACTION)
}

}  // namespace Kratos
