// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Julio Marti
//

#if !defined(KRATOS_PFEM_MELTING_APPLICATION_VARIABLES_H_INCLUDED)
#define KRATOS_PFEM_MELTING_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/global_pointer_variables.h"
//
////////
namespace Kratos
{

// Variables definition

KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_MELTING_APPLICATION, double, ACTIVATION_ENERGY)
KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_MELTING_APPLICATION, double, ARRHENIUS_COEFFICIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_MELTING_APPLICATION, double, RADIOUS)

KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( PFEM_MELTING_APPLICATION, INITIAL_POSITION)

}

#endif /* KRATOS_CONVECTION_DIFUSSION_APPLICATION_VARIABLES_H_INCLUDED */
