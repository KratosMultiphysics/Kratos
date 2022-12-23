// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#include "convection_diffusion_application_variables.h"

namespace Kratos
{

KRATOS_CREATE_VARIABLE(double, AUX_FLUX)
KRATOS_CREATE_VARIABLE(double, AUX_TEMPERATURE)
KRATOS_CREATE_VARIABLE(double, BFECC_ERROR)
KRATOS_CREATE_VARIABLE(double, BFECC_ERROR_1)
KRATOS_CREATE_VARIABLE(double, DELTA_SCALAR1)
KRATOS_CREATE_VARIABLE(double, MEAN_SIZE)
KRATOS_CREATE_VARIABLE(double, MEAN_VEL_OVER_ELEM_SIZE)
KRATOS_CREATE_VARIABLE(double, MELT_TEMPERATURE_1)
KRATOS_CREATE_VARIABLE(double, MELT_TEMPERATURE_2)
KRATOS_CREATE_VARIABLE(double, EMBEDDED_SCALAR)
KRATOS_CREATE_VARIABLE(double, PROJECTED_SCALAR1)
KRATOS_CREATE_VARIABLE(double, TRANSFER_COEFFICIENT)
KRATOS_CREATE_VARIABLE(double, ADJOINT_HEAT_TRANSFER)
KRATOS_CREATE_VARIABLE(double, SCALAR_PROJECTION)

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)

}
