//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.3 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d.h"
#include "geometries/line_2d_3.h"
#include "geometries/point_3d.h"
#include "pfem_fluid_dynamic_application.h"
#include "includes/variables.h"


namespace Kratos
{
	//Example
// 	KRATOS_CREATE_VARIABLE(double, AUX_MESH_VAR)
//	KRATOS_CREATE_VARIABLE(double, IS_INTERFACE);
//	KRATOS_CREATE_VARIABLE(double, NODAL_AREA);
//
//    KRATOS_CREATE_VARIABLE( double, NODAL_AREA )
//   KRATOS_CREATE_VARIABLE( double, IS_STRUCTURE )
//    KRATOS_CREATE_VARIABLE( double, IS_FLUID )
//    KRATOS_CREATE_VARIABLE( double, IS_BOUNDARY )
//    KRATOS_CREATE_VARIABLE( double, IS_FREE_SURFACE )

//    KRATOS_CREATE_VARIABLE( double, NODAL_MASS )
//    KRATOS_CREATE_VARIABLE( double, AUX_INDEX )
//    KRATOS_CREATE_VARIABLE( double, EXTERNAL_PRESSURE )
//    KRATOS_CREATE_VARIABLE( double, PRESSURE_OLD_IT )
//    KRATOS_CREATE_VARIABLE( double, IMPOSED_VELOCITY_X_VALUE )
//    KRATOS_CREATE_VARIABLE( double, IMPOSED_VELOCITY_Y_VALUE )
//    KRATOS_CREATE_VARIABLE( double, IMPOSED_VELOCITY_Z_VALUE )

//    KRATOS_CREATE_VARIABLE( double, DIVPROJ )
//    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( ADVPROJ )
//    KRATOS_CREATE_VARIABLE( std::string, IDENTIFIER )

    KratosPFEMFluidDynamicApplication::KratosPFEMFluidDynamicApplication()
 	{}
 	
 	void KratosPFEMFluidDynamicApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosPFEMFluidDynamicApplication... " << std::endl;
 
// 		KRATOS_REGISTER_VARIABLE( AUX_MESH_VAR )
// 		KRATOS_REGISTER_VARIABLE(IS_INTERFACE);
// 		KRATOS_REGISTER_VARIABLE(NODAL_AREA);

//        KRATOS_REGISTER_VARIABLE( NODAL_AREA )
//        KRATOS_REGISTER_VARIABLE( IS_STRUCTURE )
//        KRATOS_REGISTER_VARIABLE( IS_FLUID )
//        KRATOS_REGISTER_VARIABLE( IS_BOUNDARY )
//        KRATOS_REGISTER_VARIABLE( IS_FREE_SURFACE )

//        KRATOS_REGISTER_VARIABLE( NODAL_MASS )
//        KRATOS_REGISTER_VARIABLE( AUX_INDEX )
//        KRATOS_REGISTER_VARIABLE( EXTERNAL_PRESSURE )
//        KRATOS_REGISTER_VARIABLE( PRESSURE_OLD_IT )
//        KRATOS_REGISTER_VARIABLE( IMPOSED_VELOCITY_X_VALUE )
//        KRATOS_REGISTER_VARIABLE( IMPOSED_VELOCITY_Y_VALUE )
//        KRATOS_REGISTER_VARIABLE( IMPOSED_VELOCITY_Z_VALUE )

//        KRATOS_REGISTER_VARIABLE( DIVPROJ )
//        KRATOS_REGISTER_VARIABLE( ADVPROJ )
//        KRATOS_REGISTER_VARIABLE( IDENTIFIER )
//        KRATOS_REGISTER_VARIABLE(TAUONE);
//        KRATOS_REGISTER_VARIABLE(TAUTWO);
//        KRATOS_REGISTER_VARIABLE(SUBSCALE_VELOCITY);
//        KRATOS_REGISTER_VARIABLE(VORTICITY);



        std::cout << "Initializing KratosPFEMFluidDynamicApplication...elements succesfully registered " << std::endl;
 	}

}  // namespace Kratos.


