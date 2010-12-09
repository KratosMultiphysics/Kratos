#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
    KRATOS_CREATE_VARIABLE(double,TAUONE)
    KRATOS_CREATE_VARIABLE(double,TAUTWO)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VORTICITY)
}
