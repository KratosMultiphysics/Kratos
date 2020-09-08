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

<<<<<<< HEAD
    // Variables
    KRATOS_REGISTER_VARIABLE(BREP_ID)
    // CONDITIONS
    KRATOS_REGISTER_CONDITION("LoadCondition", mLoadCondition)
    KRATOS_REGISTER_CONDITION("PenaltyCouplingCondition", mPenaltyCouplingCondition)

    KRATOS_REGISTER_MODELER("IgaModeler", mIgaModeler);

    // VARIABLES
    KRATOS_REGISTER_VARIABLE(CROSS_AREA)
    KRATOS_REGISTER_VARIABLE(PRESTRESS_CAUCHY)

>>>>>>> origin/master
    KRATOS_REGISTER_VARIABLE(RAYLEIGH_ALPHA)
    KRATOS_REGISTER_VARIABLE(RAYLEIGH_BETA)

        //Postprocessing variables
        KRATOS_REGISTER_VARIABLE(STRESS_RESULTANT_FORCE)
        KRATOS_REGISTER_VARIABLE(STRESS_RESULTANT_MOMENT)

    KRATOS_REGISTER_VARIABLE(POINT_LOAD)
    KRATOS_REGISTER_VARIABLE(LINE_LOAD)
    KRATOS_REGISTER_VARIABLE(SURFACE_LOAD)

    KRATOS_REGISTER_VARIABLE(PENALTY_FACTOR)
}

}  // namespace Kratos
