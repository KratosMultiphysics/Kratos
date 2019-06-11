#include "chimera_application_variables.h"

namespace Kratos
{
// Flag for distinguishing b/w velocity and pressure constraints. Used in fractional-step approach
KRATOS_CREATE_FLAG(FS_CHIMERA_VEL_CONSTRAINT, 100);
KRATOS_CREATE_FLAG(FS_CHIMERA_PRE_CONSTRAINT, 100);

}
