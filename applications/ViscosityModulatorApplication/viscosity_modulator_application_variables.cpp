// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
// | | | (_)        |  \/  |         | |
// | | | |_ ___  ___| .  . | ___   __| |
// | | | | / __|/ __| |\/| |/ _ \ / _` |
// \ \_/ / \__ \ (__| |  | | (_) | (_| |
//  \___/|_|___/\___\_|  |_/\___/ \__,_|  APPLICATION
//                                      
//
//  License:         BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aniol Sala
//

#include "viscosity_modulator_application_variables.h"

namespace Kratos
{

KRATOS_CREATE_VARIABLE(double, AUX_SCALAR)
KRATOS_CREATE_VARIABLE(double, SHOCK_CAPTURING_INTENSITY)
KRATOS_CREATE_VARIABLE(bool, USE_ANISOTROPIC_DISC_CAPTURING)

}
