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
