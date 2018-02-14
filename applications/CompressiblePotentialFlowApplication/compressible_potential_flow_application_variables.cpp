#include "compressible_potential_flow_application_variables.h"

namespace Kratos
{
KRATOS_CREATE_VARIABLE( bool, UPPER_SURFACE )
KRATOS_CREATE_VARIABLE( bool, LOWER_SURFACE )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WAKE_NORMAL )
}
