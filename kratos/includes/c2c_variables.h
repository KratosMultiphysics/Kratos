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
    KRATOS_DEFINE_VARIABLE( double, MOLD_AVERAGE_TEMPERATURE)

    KRATOS_DEFINE_VARIABLE( double,  LAST_AIR)
    KRATOS_DEFINE_VARIABLE( double,  PRESSURES)

    KRATOS_DEFINE_VARIABLE( double, MOLD_DENSITY )
    KRATOS_DEFINE_VARIABLE( double, MOLD_SPECIFIC_HEAT )
    KRATOS_DEFINE_VARIABLE( double, MOLD_THICKNESS )
    KRATOS_DEFINE_VARIABLE( double, MOLD_SFACT )
    KRATOS_DEFINE_VARIABLE( double, MOLD_VFACT )
    KRATOS_DEFINE_VARIABLE( double, MOLD_CONDUCTIVITY )
    KRATOS_DEFINE_VARIABLE( double, MOLD_HTC_ENVIRONMENT )
    KRATOS_DEFINE_VARIABLE( double, MOLD_TEMPERATURE )
    KRATOS_DEFINE_VARIABLE( double, MOLD_INNER_TEMPERATURE )

    KRATOS_DEFINE_VARIABLE( double, HTC )
    KRATOS_DEFINE_VARIABLE( double, SOLIDFRACTION )
	KRATOS_DEFINE_VARIABLE(double, SOLIDFRACTION_RATE)
    KRATOS_DEFINE_VARIABLE( double, SOLIDIF_TIME )
    KRATOS_DEFINE_VARIABLE( double, SOLIDIF_MODULUS )
    KRATOS_DEFINE_VARIABLE( double, FILLTIME )
    KRATOS_DEFINE_VARIABLE( double, MACRO_POROSITY )
    KRATOS_DEFINE_VARIABLE( double, SHRINKAGE_POROSITY )
    KRATOS_DEFINE_VARIABLE( double, MAX_VEL )
    KRATOS_DEFINE_VARIABLE( double, FRONT_MEETING )
	KRATOS_DEFINE_VARIABLE( double, LIQUID_TIME)
	KRATOS_DEFINE_VARIABLE( double, FLOW_LENGTH )
	KRATOS_DEFINE_VARIABLE( double, FLOW_LENGTH2 )
	KRATOS_DEFINE_VARIABLE( double, COOLING_RATE)
	KRATOS_DEFINE_VARIABLE (double, LIQUID_TO_SOLID_TIME)
	KRATOS_DEFINE_VARIABLE (double, TIME_CRT)
	KRATOS_DEFINE_VARIABLE(double, SINKMARK)
	KRATOS_DEFINE_VARIABLE(double, COLD_SHUTS)
    KRATOS_DEFINE_VARIABLE(double, NIYAMA)
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( TEMPERATURE_GRADIENT )
    KRATOS_DEFINE_VARIABLE(double, INITIAL_TEMPERATURE)
}  // namespace Kratos.

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT

#endif // KRATOS_C2C_VARIABLES_H_INCLUDED  defined

