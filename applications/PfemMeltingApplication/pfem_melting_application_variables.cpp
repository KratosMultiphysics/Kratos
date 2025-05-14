// KRATOS 
// _____   __               __  __      _ _   _             
//|  __ \ / _|             |  \/  |    | | | (_)            
//| |__) | |_ ___ _ __ ___ | \  / | ___| | |_ _ _ __   __ _ 
//|  ___/|  _/ _ \ '_ ` _ \| |\/| |/ _ \ | __| | '_ \ / _` |
//| |    | ||  __/ | | | | | |  | |  __/ | |_| | | | | (_| |
//|_|    |_| \___|_| |_| |_|_|  |_|\___|_|\__|_|_| |_|\__, |
//                                                     __/ |
//                                                    |___/ APPLICATION
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Julio Marti
//


#include "pfem_melting_application_variables.h"

namespace Kratos
{

KRATOS_CREATE_VARIABLE(double, ACTIVATION_ENERGY)
KRATOS_CREATE_VARIABLE(double, ARRHENIUS_COEFFICIENT)
KRATOS_CREATE_VARIABLE(double, RADIOUS)
KRATOS_CREATE_VARIABLE(double, HEAT_OF_VAPORIZATION)
KRATOS_CREATE_VARIABLE(double, ARRHENIUS_VALUE)

KRATOS_CREATE_VARIABLE(double, NODAL_VOLUME)
KRATOS_CREATE_VARIABLE(bool, IS_SOLID)

    KRATOS_CREATE_VARIABLE(double, DELTA_SIGMA_XX)
    KRATOS_CREATE_VARIABLE(double, DELTA_SIGMA_XY)
    KRATOS_CREATE_VARIABLE(double, DELTA_SIGMA_XZ)
    
    KRATOS_CREATE_VARIABLE(double, DELTA_SIGMA_YX)
    KRATOS_CREATE_VARIABLE(double, DELTA_SIGMA_YY)
    KRATOS_CREATE_VARIABLE(double, DELTA_SIGMA_YZ)
    
    KRATOS_CREATE_VARIABLE(double, DELTA_SIGMA_ZX)
    KRATOS_CREATE_VARIABLE(double, DELTA_SIGMA_ZY)
    KRATOS_CREATE_VARIABLE(double, DELTA_SIGMA_ZZ)
    
    
    KRATOS_CREATE_VARIABLE(double, HISTORICAL_SIGMA_XX)
    KRATOS_CREATE_VARIABLE(double, HISTORICAL_SIGMA_XY)
    KRATOS_CREATE_VARIABLE(double, HISTORICAL_SIGMA_XZ)
    
    KRATOS_CREATE_VARIABLE(double, HISTORICAL_SIGMA_YX)
    KRATOS_CREATE_VARIABLE(double, HISTORICAL_SIGMA_YY)
    KRATOS_CREATE_VARIABLE(double, HISTORICAL_SIGMA_YZ)

    KRATOS_CREATE_VARIABLE(double, HISTORICAL_SIGMA_ZX)
    KRATOS_CREATE_VARIABLE(double, HISTORICAL_SIGMA_ZY)
    KRATOS_CREATE_VARIABLE(double, HISTORICAL_SIGMA_ZZ)
    KRATOS_CREATE_VARIABLE(double, PRESSUREAUX)
    
    //element
    KRATOS_CREATE_VARIABLE(Matrix, TOTAL_CAUCHY_STRESS )
    
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(INITIAL_POSITION)



}
