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
    , mIgaMembraneElement(0, Element::GeometryType::Pointer(
        new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1))))
    , mIgaEdgeCableElement(0, Element::GeometryType::Pointer(
        new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1))))
    , mOutputCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node<3>>(Condition::GeometryType::PointsArrayType(1))))
    , mLoadCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node<3>>(Condition::GeometryType::PointsArrayType(1))))
    , mPenaltyCouplingCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node<3>>(Condition::GeometryType::PointsArrayType(1))))
    , mLagrangeCouplingCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node<3>>(Condition::GeometryType::PointsArrayType(1))))
    , mNitscheCouplingCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node<3>>(Condition::GeometryType::PointsArrayType(1))))
    , mSupportPenaltyCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node<3>>(Condition::GeometryType::PointsArrayType(1))))
    , mSupportLagrangeCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node<3>>(Condition::GeometryType::PointsArrayType(1))))
    , mSupportNitscheCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node<3>>(Condition::GeometryType::PointsArrayType(1))))
{
}

void KratosIgaApplication::Register() {

KRATOS_INFO("") << "    KRATOS  _____ _____\n"
                << "           |_   _/ ____|   /\\\n"
                << "             | || |  __   /  \\\n"
                << "             | || | |_ | / /\\ \\\n"
                << "            _| || |__| |/ ____ \\\n"
                << "           |_____\\_____/_/    \\_\\\n"
                << "Initializing KratosIgaApplication..." << std::endl;

    // ELEMENTS
    KRATOS_REGISTER_ELEMENT("Shell3pElement", mShell3pElement)
    KRATOS_REGISTER_ELEMENT("IgaMembraneElement", mIgaMembraneElement)
    KRATOS_REGISTER_ELEMENT("IgaEdgeCableElement", mIgaEdgeCableElement)

    // CONDITIONS
    KRATOS_REGISTER_CONDITION("OutputCondition", mOutputCondition)
    KRATOS_REGISTER_CONDITION("LoadCondition", mLoadCondition)
    KRATOS_REGISTER_CONDITION("PenaltyCouplingCondition", mPenaltyCouplingCondition)
    KRATOS_REGISTER_CONDITION("LagrangeCouplingCondition", mLagrangeCouplingCondition)
    KRATOS_REGISTER_CONDITION("NitscheCouplingCondition", mNitscheCouplingCondition)
    KRATOS_REGISTER_CONDITION("SupportPenaltyCondition", mSupportPenaltyCondition)
    KRATOS_REGISTER_CONDITION("SupportLagrangeCondition", mSupportLagrangeCondition)
    KRATOS_REGISTER_CONDITION("SupportNitscheCondition", mSupportNitscheCondition)

    KRATOS_REGISTER_MODELER("IgaModeler", mIgaModeler);

    // VARIABLES
    KRATOS_REGISTER_VARIABLE(NURBS_CONTROL_POINT_WEIGHT)

    KRATOS_REGISTER_VARIABLE(COORDINATES)
    KRATOS_REGISTER_VARIABLE(TANGENTS)
    KRATOS_REGISTER_VARIABLE(LOCAL_ELEMENT_ORIENTATION)
    KRATOS_REGISTER_VARIABLE(LOCAL_PRESTRESS_AXIS_1)
    KRATOS_REGISTER_VARIABLE(LOCAL_PRESTRESS_AXIS_2)

    KRATOS_REGISTER_VARIABLE(CABLE_FORCE)
    KRATOS_REGISTER_VARIABLE(CROSS_AREA)
    KRATOS_REGISTER_VARIABLE(PRESTRESS_CAUCHY)
    KRATOS_REGISTER_VARIABLE(PRESTRESS)
    KRATOS_REGISTER_VARIABLE(PRINCIPAL_STRESS_1)
    KRATOS_REGISTER_VARIABLE(PRINCIPAL_STRESS_2)

    KRATOS_REGISTER_VARIABLE(RAYLEIGH_ALPHA)
    KRATOS_REGISTER_VARIABLE(RAYLEIGH_BETA)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(POINT_LOAD)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LINE_LOAD)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SURFACE_LOAD)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DEAD_LOAD)

    KRATOS_REGISTER_VARIABLE(PENALTY_FACTOR)
    KRATOS_REGISTER_VARIABLE(NITSCHE_STABILIZATION_PARAMETER)

    KRATOS_REGISTER_VARIABLE(CHARACTERISTIC_LENGTH)
    KRATOS_REGISTER_VARIABLE(LOCAL_TANGENTS_MASTER)
    KRATOS_REGISTER_VARIABLE(LOCAL_TANGENTS_SLAVE)

    // Generalized eigenvalue problem
    KRATOS_REGISTER_VARIABLE(BUILD_LEVEL)
    KRATOS_REGISTER_VARIABLE(EIGENVALUE_VECTOR)
    KRATOS_REGISTER_VARIABLE(EIGENVECTOR_MATRIX)

    //Lagrange Multiplier Variables
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VECTOR_LAGRANGE_MULTIPLIER_VELOCITY)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VECTOR_LAGRANGE_MULTIPLIER_ACCELERATION)
}

}  // namespace Kratos
