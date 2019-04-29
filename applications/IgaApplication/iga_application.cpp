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
    //, mIgaTrussElement(0, Element::GeometryType::Pointer(
    //    new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1))))
    , mShellKLDiscreteElement(0, Element::GeometryType::Pointer(
        new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1))))
    , mIgaShell3pElement(0, Element::GeometryType::Pointer(
        new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1))))
    , mIgaCheckCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node<3>>(Condition::GeometryType::PointsArrayType(1))))
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
    //KRATOS_REGISTER_ELEMENT("IgaTrussElement", mIgaTrussElement)
    KRATOS_REGISTER_ELEMENT("ShellKLDiscreteElement", mShellKLDiscreteElement)
    KRATOS_REGISTER_ELEMENT("IgaShell3pElement", mShellKLDiscreteElement)

    KRATOS_REGISTER_CONDITION("IgaCheckCondition", mIgaCheckCondition)

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

    KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_VALUES_SLAVE)
    KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE)
    KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES_SLAVE)

    KRATOS_REGISTER_VARIABLE(RAYLEIGH_ALPHA)
    KRATOS_REGISTER_VARIABLE(RAYLEIGH_BETA)
    KRATOS_REGISTER_VARIABLE(NODAL_INERTIA)

        //Postprocessing variables
        KRATOS_REGISTER_VARIABLE(STRESS_RESULTANT_FORCE)
        KRATOS_REGISTER_VARIABLE(STRESS_RESULTANT_MOMENT)

    KRATOS_REGISTER_VARIABLE(POINT_LOAD)
    KRATOS_REGISTER_VARIABLE(LINE_LOAD)
    KRATOS_REGISTER_VARIABLE(SURFACE_LOAD)

    KRATOS_REGISTER_VARIABLE(PENALTY_FACTOR)
}

}  // namespace Kratos
