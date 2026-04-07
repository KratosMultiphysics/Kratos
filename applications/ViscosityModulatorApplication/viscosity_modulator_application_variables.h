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

#if !defined(KRATOS_VISCOSITY_MODULATOR_APPLICATION_VARIABLES_H_INCLUDED)
#define KRATOS_VISCOSITY_MODULATOR_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/global_pointer_variables.h"

namespace Kratos
{

// Variables definition
KRATOS_DEFINE_APPLICATION_VARIABLE( VISCOSITY_MODULATOR_APPLICATION, double, AUX_SCALAR)
KRATOS_DEFINE_APPLICATION_VARIABLE( VISCOSITY_MODULATOR_APPLICATION, double, SHOCK_CAPTURING_INTENSITY)
KRATOS_DEFINE_APPLICATION_VARIABLE( VISCOSITY_MODULATOR_APPLICATION, bool, USE_ANISOTROPIC_DISC_CAPTURING)

}

#endif /* KRATOS_VISCOSITY_MODULATOR_APPLICATION_VARIABLES_H_INCLUDED */
