//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//










#if !defined(KRATOS_CFD_VARIABLES_H_INCLUDED )
#define  KRATOS_CFD_VARIABLES_H_INCLUDED



// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "includes/define.h"
#include "containers/variable.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"
#include "includes/kratos_components.h"
#include "includes/ublas_interface.h"
#include "containers/array_1d.h"
#include "containers/weak_pointer_vector.h"
#include "containers/periodic_variables_container.h"

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_API

//TODO: move to the FluidDynamics application
namespace Kratos
{

    // Useful variables

    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( ADVPROJ )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( CONV_PROJ )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PRESS_PROJ )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( MATERIAL_ACCELERATION )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( VORTICITY )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( RELAXED_ACCELERATION )

    KRATOS_DEFINE_VARIABLE( double, DIVPROJ )
    KRATOS_DEFINE_VARIABLE( double, PRESSURE_OLD_IT )
    KRATOS_DEFINE_VARIABLE( double, C_SMAGORINSKY )
    KRATOS_DEFINE_VARIABLE( double, CFL_NUMBER )
    KRATOS_DEFINE_VARIABLE( double, MOLECULAR_VISCOSITY )
    KRATOS_DEFINE_VARIABLE( double, TURBULENT_VISCOSITY )
    KRATOS_DEFINE_VARIABLE( double, Y_WALL)
    KRATOS_DEFINE_VARIABLE( double, PRESSURE_COEFFICIENT)
    KRATOS_DEFINE_VARIABLE( int, FRACTIONAL_STEP )
    KRATOS_DEFINE_VARIABLE( int, OSS_SWITCH )

    // Legacy variables
    KRATOS_DEFINE_VARIABLE( double, DYNAMIC_TAU )
    KRATOS_DEFINE_VARIABLE( double, DYNAMIC_VISCOSITY)
    KRATOS_DEFINE_VARIABLE( double, EFFECTIVE_VISCOSITY )
    KRATOS_DEFINE_VARIABLE( double, KINEMATIC_VISCOSITY)
    KRATOS_DEFINE_VARIABLE( double, THAWONE )
    KRATOS_DEFINE_VARIABLE( double, THAWTWO )
    KRATOS_DEFINE_VARIABLE( double, M )

    KRATOS_DEFINE_VARIABLE( double, CROSS_WIND_STABILIZATION_FACTOR )

} // namespace Kratos

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT

#endif // KRATOS_CFD_VARIABLES_H_INCLUDED  defined
