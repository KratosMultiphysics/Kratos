//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//











#if !defined(KRATOS_C2C_VARIABLES_H_INCLUDED )
#define  KRATOS_C2C_VARIABLES_H_INCLUDED



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

//TODO: move to the KratosC2C application or eventually to the FluidDynamicsAsNeeded
namespace Kratos
{

    KRATOS_DEFINE_VARIABLE( double, LATENT_HEAT )
    KRATOS_DEFINE_VARIABLE( double, SOLID_TEMPERATURE )
    KRATOS_DEFINE_VARIABLE( double, FLUID_TEMPERATURE )
    KRATOS_DEFINE_VARIABLE( double, AVERAGE_TEMPERATURE )
    KRATOS_DEFINE_VARIABLE( double, INLET_TEMPERATURE)
//     KRATOS_DEFINE_VARIABLE( double, AMBIENT_TEMPERATURE )
    KRATOS_DEFINE_VARIABLE( double, Y_WALL)
    KRATOS_DEFINE_VARIABLE( double, COUNTER )
    KRATOS_DEFINE_VARIABLE( double, DISTANCE_CORRECTION )
    KRATOS_DEFINE_VARIABLE( double, COMPUTED_DISTANCE )
    KRATOS_DEFINE_VARIABLE( double, MATERIAL )
    KRATOS_DEFINE_VARIABLE( double, MOULD_AVERAGE_TEMPERATURE)

    KRATOS_DEFINE_VARIABLE( double,  LAST_AIR)
    KRATOS_DEFINE_VARIABLE( double,  PRESSURES)
    KRATOS_DEFINE_VARIABLE( double,  TEMPERATURES)
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( VELOCITIES )

    KRATOS_DEFINE_VARIABLE( double, MOULD_DENSITY )
    KRATOS_DEFINE_VARIABLE( double, MOULD_SPECIFIC_HEAT )
    KRATOS_DEFINE_VARIABLE( double, MOULD_THICKNESS )
    KRATOS_DEFINE_VARIABLE( double, MOULD_SFACT )
    KRATOS_DEFINE_VARIABLE( double, MOULD_VFACT )
    KRATOS_DEFINE_VARIABLE( double, MOULD_CONDUCTIVITY )
    KRATOS_DEFINE_VARIABLE( double, MOULD_HTC_ENVIRONMENT )
    KRATOS_DEFINE_VARIABLE( double, MOULD_TEMPERATURE )
    KRATOS_DEFINE_VARIABLE( double, MOULD_INNER_TEMPERATURE )

    KRATOS_DEFINE_VARIABLE( int, NODE_PROPERTY_ID )
    KRATOS_DEFINE_VARIABLE( double, HTC )
    KRATOS_DEFINE_VARIABLE( int, REF_ID )
    KRATOS_DEFINE_VARIABLE( double, PARTICLE_RADIUS )
    KRATOS_DEFINE_VARIABLE( double, POSETIVE_DISTANCE )
    KRATOS_DEFINE_VARIABLE( double, NAGATIVE_DISTANCE )
    KRATOS_DEFINE_VARIABLE( bool, IS_ESCAPED )
    KRATOS_DEFINE_VARIABLE( int, IS_SOLIDIFIED )
    KRATOS_DEFINE_VARIABLE( double, SOLIDFRACTION )
	KRATOS_DEFINE_VARIABLE(double, SOLIDFRACTION_RATE)
    KRATOS_DEFINE_VARIABLE( double, SOLIDIF_TIME )
    KRATOS_DEFINE_VARIABLE( double, SOLIDIF_MODULUS )
    KRATOS_DEFINE_VARIABLE( double, FILLTIME )
    KRATOS_DEFINE_VARIABLE( double, MACRO_POROSITY )
    KRATOS_DEFINE_VARIABLE( double, SHRINKAGE_POROSITY )
    KRATOS_DEFINE_VARIABLE( double, MAX_VEL )
    KRATOS_DEFINE_VARIABLE( int, IS_GRAVITY_FILLING )
    KRATOS_DEFINE_VARIABLE( double, VOLUME_FRACTION )
    KRATOS_DEFINE_VARIABLE( double, KAPPA )
    KRATOS_DEFINE_VARIABLE( double, EPSILON )
    KRATOS_DEFINE_VARIABLE( double, SHRINKAGE_POROSITY_US )
    KRATOS_DEFINE_VARIABLE( double, SOLIDIF_MODULUS_US )
    KRATOS_DEFINE_VARIABLE( double, TEMPERATURES_US )
    KRATOS_DEFINE_VARIABLE( double, FRONT_MEETING )
	KRATOS_DEFINE_VARIABLE( double, LIQUID_TIME)
	KRATOS_DEFINE_VARIABLE( double, FLOW_LENGTH )
	KRATOS_DEFINE_VARIABLE( double, FLOW_LENGTH2 )
	KRATOS_DEFINE_VARIABLE( double, COOLING_RATE)
	KRATOS_DEFINE_VARIABLE (double, LIQUID_TO_SOLID_TIME)
	KRATOS_DEFINE_VARIABLE (double, TIME_CRT)
	KRATOS_DEFINE_VARIABLE(double, SINKMARK)

}  // namespace Kratos.

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT

#endif // KRATOS_C2C_VARIABLES_H_INCLUDED  defined

