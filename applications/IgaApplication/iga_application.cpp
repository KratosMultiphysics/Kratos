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
    , mOutputCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node<3>>(Condition::GeometryType::PointsArrayType(1))))
    , mLoadCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node<3>>(Condition::GeometryType::PointsArrayType(1))))
    , mPenaltyCouplingCondition(0, Condition::GeometryType::Pointer(
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

    // CONDITIONS
    KRATOS_REGISTER_CONDITION("OutputCondition", mOutputCondition)
    KRATOS_REGISTER_CONDITION("LoadCondition", mLoadCondition)
    KRATOS_REGISTER_CONDITION("PenaltyCouplingCondition", mPenaltyCouplingCondition)

    KRATOS_REGISTER_MODELER("IgaModeler", mIgaModeler);

    // VARIABLES
    KRATOS_REGISTER_VARIABLE(CROSS_AREA)
    KRATOS_REGISTER_VARIABLE(PRESTRESS_CAUCHY)

    KRATOS_REGISTER_VARIABLE(RAYLEIGH_ALPHA)
    KRATOS_REGISTER_VARIABLE(RAYLEIGH_BETA)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(POINT_LOAD)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LINE_LOAD)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SURFACE_LOAD)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DEAD_LOAD)

    KRATOS_REGISTER_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(PK2_STRESS)
    KRATOS_REGISTER_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(CAUCHY_STRESS)
    KRATOS_REGISTER_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(CAUCHY_STRESS_TOP)
    KRATOS_REGISTER_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(CAUCHY_STRESS_BOTTOM)
    KRATOS_REGISTER_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(MEMBRANE_FORCE)
    KRATOS_REGISTER_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(INTERNAL_MOMENT)
        
    KRATOS_REGISTER_VARIABLE(SHEAR_FORCE_1)
    KRATOS_REGISTER_VARIABLE(SHEAR_FORCE_2)

    KRATOS_REGISTER_VARIABLE(PENALTY_FACTOR)
}

}  // namespace Kratos
