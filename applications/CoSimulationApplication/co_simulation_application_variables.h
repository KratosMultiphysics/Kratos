// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//					 license: CoSimulationApplication/license.txt
//
//  Main authors:    Aditya Ghantasala
//                   Philipp Bucher
//

#if !defined(KRATOS_CO_SIMULATION_APPLICATION_VARIABLES_H_INCLUDED)
#define KRATOS_CO_SIMULATION_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable.h"

namespace Kratos
{

KRATOS_DEFINE_APPLICATION_VARIABLE( CO_SIMULATION_APPLICATION, double, SCALAR_DISPLACEMENT );
KRATOS_DEFINE_APPLICATION_VARIABLE( CO_SIMULATION_APPLICATION, double, SCALAR_ROOT_POINT_DISPLACEMENT );
KRATOS_DEFINE_APPLICATION_VARIABLE( CO_SIMULATION_APPLICATION, double, SCALAR_REACTION );
KRATOS_DEFINE_APPLICATION_VARIABLE( CO_SIMULATION_APPLICATION, double, SCALAR_FORCE );
KRATOS_DEFINE_APPLICATION_VARIABLE( CO_SIMULATION_APPLICATION, double, SCALAR_VOLUME_ACCELERATION );

}

#endif /* KRATOS_CO_SIMULATION_APPLICATION_VARIABLES_H_INCLUDED */