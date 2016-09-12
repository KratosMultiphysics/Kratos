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












// This define must be HERE
#define DKRATOS_EXPORT_INTERFACE_2 1

// System includes
#include <string>
#include <iostream>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/c2c_variables.h"
#include "includes/kernel.h"
#include "includes/node.h"
// #include "includes/element.h"
// #include "includes/condition.h"
// #include "includes/constitutive_law.h"
// #include "includes/geometrical_object.h"

// #include "geometries/line_2d.h"
// #include "geometries/line_2d_2.h"
// #include "geometries/line_2d_3.h"
// #include "geometries/line_3d_2.h"
// #include "geometries/line_3d_3.h"
// #include "geometries/point.h"
// #include "geometries/point_2d.h"
// #include "geometries/point_3d.h"
// #include "geometries/sphere_3d_1.h"
// #include "geometries/triangle_2d_3.h"
// #include "geometries/triangle_2d_6.h"
// #include "geometries/triangle_3d_3.h"
// #include "geometries/triangle_3d_6.h"
// #include "geometries/quadrilateral_2d_4.h"
// #include "geometries/quadrilateral_2d_8.h"
// #include "geometries/quadrilateral_2d_9.h"
// #include "geometries/quadrilateral_3d_4.h"
// #include "geometries/quadrilateral_3d_8.h"
// #include "geometries/quadrilateral_3d_9.h"
// #include "geometries/tetrahedra_3d_4.h"
// #include "geometries/tetrahedra_3d_10.h"
// #include "geometries/prism_3d_6.h"
// #include "geometries/prism_3d_15.h"
// #include "geometries/hexahedra_3d_8.h"
// #include "geometries/hexahedra_3d_20.h"
// #include "geometries/hexahedra_3d_27.h"

// #include "python/add_c2c_variables_to_python.h"

// #include "includes/convection_diffusion_settings.h"
// #include "includes/radiation_settings.h"

#include "includes/kratos_flags.h"

namespace Kratos
{

    KRATOS_CREATE_VARIABLE( double, SOLID_TEMPERATURE )
    KRATOS_CREATE_VARIABLE( double, FLUID_TEMPERATURE )
    KRATOS_CREATE_VARIABLE( double, AVERAGE_TEMPERATURE )
    KRATOS_CREATE_VARIABLE( double, INLET_TEMPERATURE)
//     KRATOS_CREATE_VARIABLE( double, AMBIENT_TEMPERATURE )
//     KRATOS_CREATE_VARIABLE( double, Y_WALL)
    KRATOS_CREATE_VARIABLE( double, COUNTER )
    KRATOS_CREATE_VARIABLE( double, DISTANCE_CORRECTION )
    KRATOS_CREATE_VARIABLE( double, COMPUTED_DISTANCE )
    KRATOS_CREATE_VARIABLE( double, MATERIAL )

    Kratos::Variable<double> LAST_AIR( "LAST AIR" );
    Kratos::Variable<double> PRESSURES( "PRESSURES" );
    Kratos::Variable<Kratos::array_1d<double, 3> > VELOCITIES( "VELOCITIES", Kratos::zero_vector<double>( 3 ) );
    Kratos::Variable<double> TEMPERATURES( "TEMPERATURES" );
    /*const*/
    Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > >
    VELOCITIES_X( "X-VELOCITIES", Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> >( VELOCITIES, 0 ) );

    /*const*/
    Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > >
    VELOCITIES_Y( "Y-VELOCITIES)", Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> >( VELOCITIES, 1 ) );

    /*const*/
    Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > >
    VELOCITIES_Z( "Z-VELOCITIES", Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> >( VELOCITIES, 2 ) );

    // for Vulcan application virtual mould properties
    KRATOS_CREATE_VARIABLE(double,  MOULD_DENSITY)
    KRATOS_CREATE_VARIABLE(double,  MOULD_SPECIFIC_HEAT)
    KRATOS_CREATE_VARIABLE(double,  MOULD_THICKNESS)
    KRATOS_CREATE_VARIABLE(double,  MOULD_SFACT)
    KRATOS_CREATE_VARIABLE(double,  MOULD_VFACT)
    KRATOS_CREATE_VARIABLE(double,  MOULD_CONDUCTIVITY)
    KRATOS_CREATE_VARIABLE(double,  MOULD_HTC_ENVIRONMENT)
    KRATOS_CREATE_VARIABLE(double,  MOULD_TEMPERATURE)
    KRATOS_CREATE_VARIABLE(double,  MOULD_INNER_TEMPERATURE)
    // for Click2Cast Application
    KRATOS_CREATE_VARIABLE(int, NODE_PROPERTY_ID)
    KRATOS_CREATE_VARIABLE(double,  HTC)
    KRATOS_CREATE_VARIABLE(int, REF_ID)
    KRATOS_CREATE_VARIABLE(double, PARTICLE_RADIUS)
    KRATOS_CREATE_VARIABLE(double, POSETIVE_DISTANCE)
    KRATOS_CREATE_VARIABLE(double, NAGATIVE_DISTANCE)
    KRATOS_CREATE_VARIABLE(bool, IS_ESCAPED)
    KRATOS_CREATE_VARIABLE(int, IS_SOLIDIFIED)
    Kratos::Variable<double> SOLIDFRACTION( "SOLID_FRACTION" );
	KRATOS_CREATE_VARIABLE(double, SOLIDFRACTION_RATE)
    Kratos::Variable<double> SOLIDIF_TIME( "SOLIDIF TIME" );
    Kratos::Variable<double> SOLIDIF_MODULUS( "SOLIDIF MODULUS" );
    Kratos::Variable<double> FILLTIME( "FILLTIME" );
    KRATOS_CREATE_VARIABLE(double, MACRO_POROSITY )
    Kratos::Variable<double> SHRINKAGE_POROSITY( "SHRINKAGE_POROSITY" );
    Kratos::Variable<double> MAX_VEL( "MAX VEL" );
    KRATOS_CREATE_VARIABLE(int, IS_GRAVITY_FILLING)
    KRATOS_CREATE_VARIABLE(double, VOLUME_FRACTION )
    KRATOS_CREATE_VARIABLE(double, KAPPA )
    KRATOS_CREATE_VARIABLE(double, EPSILON )
    Kratos::Variable<double> SHRINKAGE_POROSITY_US( "SHRINKAGE_POROSITY" );
    Kratos::Variable<double> SOLIDIF_MODULUS_US( "SOLIDIF MODULUS" );
    Kratos::Variable<double> TEMPERATURES_US( "TEMPERATURES" );
    KRATOS_CREATE_VARIABLE(double,FRONT_MEETING)
    KRATOS_CREATE_VARIABLE( double, MOULD_AVERAGE_TEMPERATURE )  
	Kratos::Variable<double> LIQUID_TIME("TIME TO START SOLIDIFICATION");
	Kratos::Variable<double>  FLOW_LENGTH ("FLOW LENGTH ESTIMATION 1");
	Kratos::Variable<double>  FLOW_LENGTH2("FLOW LENGTH ESTIMATION 2" );
	Kratos::Variable<double> COOLING_RATE("COOLING RATE");
	Kratos::Variable<double> LIQUID_TO_SOLID_TIME("TIME TO COMPLETE SOLIDIF.");
	//KRATOS_CREATE_VARIABLE(double, SOLID_FRACTION)
	//KRATOS_CREATE_VARIABLE(double, SOLID_FRACTION_RATE)
	KRATOS_CREATE_VARIABLE( double, TIME_CRT )
	KRATOS_CREATE_VARIABLE(double, SINKMARK)

  void KratosApplication::RegisterC2CVariables()
  {
        KRATOS_REGISTER_VARIABLE(  SOLID_TEMPERATURE )
        KRATOS_REGISTER_VARIABLE(  FLUID_TEMPERATURE )
        KRATOS_REGISTER_VARIABLE(  AVERAGE_TEMPERATURE )
        KRATOS_REGISTER_VARIABLE(  INLET_TEMPERATURE)
//         KRATOS_REGISTER_VARIABLE(  AMBIENT_TEMPERATURE )
//         KRATOS_REGISTER_VARIABLE(  Y_WALL)
        KRATOS_REGISTER_VARIABLE(  COUNTER )
        KRATOS_REGISTER_VARIABLE(  DISTANCE_CORRECTION )
        KRATOS_REGISTER_VARIABLE(  COMPUTED_DISTANCE )
        KRATOS_REGISTER_VARIABLE(  MATERIAL )
        KRATOS_REGISTER_VARIABLE(  MOULD_AVERAGE_TEMPERATURE)

        KRATOS_REGISTER_VARIABLE(   LAST_AIR)
        KRATOS_REGISTER_VARIABLE(   PRESSURES)
        KRATOS_REGISTER_VARIABLE(   TEMPERATURES)
        // KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(VELOCITIES)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VELOCITIES)

        KRATOS_REGISTER_VARIABLE(   MOULD_DENSITY)
        KRATOS_REGISTER_VARIABLE(   MOULD_SPECIFIC_HEAT)
        KRATOS_REGISTER_VARIABLE(   MOULD_THICKNESS)
        KRATOS_REGISTER_VARIABLE(   MOULD_SFACT)
        KRATOS_REGISTER_VARIABLE(   MOULD_VFACT)
        KRATOS_REGISTER_VARIABLE(   MOULD_CONDUCTIVITY)
        KRATOS_REGISTER_VARIABLE(   MOULD_HTC_ENVIRONMENT)
        KRATOS_REGISTER_VARIABLE(   MOULD_TEMPERATURE)
        KRATOS_REGISTER_VARIABLE(   MOULD_INNER_TEMPERATURE)

        KRATOS_REGISTER_VARIABLE(  NODE_PROPERTY_ID)
        KRATOS_REGISTER_VARIABLE(   HTC)
        KRATOS_REGISTER_VARIABLE(  REF_ID)
        KRATOS_REGISTER_VARIABLE(  PARTICLE_RADIUS)
        KRATOS_REGISTER_VARIABLE(  POSETIVE_DISTANCE)
        KRATOS_REGISTER_VARIABLE(  NAGATIVE_DISTANCE)
        KRATOS_REGISTER_VARIABLE( IS_ESCAPED)
        KRATOS_REGISTER_VARIABLE(  IS_SOLIDIFIED)  
        KRATOS_REGISTER_VARIABLE(  SOLIDFRACTION )
		KRATOS_REGISTER_VARIABLE( SOLIDFRACTION_RATE)
        KRATOS_REGISTER_VARIABLE(  SOLIDIF_TIME  )
        KRATOS_REGISTER_VARIABLE(  SOLIDIF_MODULUS  )
        KRATOS_REGISTER_VARIABLE(  FILLTIME  )
        KRATOS_REGISTER_VARIABLE(  MACRO_POROSITY  )
        KRATOS_REGISTER_VARIABLE(  SHRINKAGE_POROSITY  )
        KRATOS_REGISTER_VARIABLE(  MAX_VEL  )
        KRATOS_REGISTER_VARIABLE(  IS_GRAVITY_FILLING)
        KRATOS_REGISTER_VARIABLE( VOLUME_FRACTION )

        KRATOS_REGISTER_VARIABLE( KAPPA ) 
        KRATOS_REGISTER_VARIABLE( EPSILON )
        KRATOS_REGISTER_VARIABLE( SHRINKAGE_POROSITY_US)
        KRATOS_REGISTER_VARIABLE( SOLIDIF_MODULUS_US)
        KRATOS_REGISTER_VARIABLE( TEMPERATURES_US)
        KRATOS_REGISTER_VARIABLE( FRONT_MEETING)
		KRATOS_REGISTER_VARIABLE( LIQUID_TIME)
		KRATOS_REGISTER_VARIABLE( FLOW_LENGTH )
		KRATOS_REGISTER_VARIABLE( FLOW_LENGTH2 )
		KRATOS_REGISTER_VARIABLE( COOLING_RATE)
		KRATOS_REGISTER_VARIABLE (LIQUID_TO_SOLID_TIME)
		KRATOS_REGISTER_VARIABLE (TIME_CRT)
		KRATOS_REGISTER_VARIABLE(SINKMARK)
  }


}  // namespace Kratos.

// This define must be HERE
#undef DKRATOS_EXPORT_INTERFACE_2



