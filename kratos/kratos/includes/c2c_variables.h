// Kratos Multi-Physics
// 
// Copyright (c) 2015, Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
// 
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement: 
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// 	
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY 
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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
#include "includes/kratos_export_dll.h"

#undef KRATOS_DEFINE_VARIABLE
#undef KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS
#define KRATOS_DEFINE_VARIABLE KRATOS_DEFINE_VARIABLE_DLL
#define KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS_DLL

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
    KRATOS_DEFINE_VARIABLE(double, MOULD_AVERAGE_TEMPERATURE)

    KRATOS_DEFINE_VARIABLE(double,  LAST_AIR)
    KRATOS_DEFINE_VARIABLE(double,  PRESSURES)
    KRATOS_DEFINE_VARIABLE(double,  TEMPERATURES)
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( VELOCITIES )

    KRATOS_DEFINE_VARIABLE(double,  MOULD_DENSITY)
    KRATOS_DEFINE_VARIABLE(double,  MOULD_SPECIFIC_HEAT)
    KRATOS_DEFINE_VARIABLE(double,  MOULD_THICKNESS)
    KRATOS_DEFINE_VARIABLE(double,  MOULD_SFACT)
    KRATOS_DEFINE_VARIABLE(double,  MOULD_VFACT)
    KRATOS_DEFINE_VARIABLE(double,  MOULD_CONDUCTIVITY)
    KRATOS_DEFINE_VARIABLE(double,  MOULD_HTC_ENVIRONMENT)
    KRATOS_DEFINE_VARIABLE(double,  MOULD_TEMPERATURE)
    KRATOS_DEFINE_VARIABLE(double,  MOULD_INNER_TEMPERATURE)

    KRATOS_DEFINE_VARIABLE(int, NODE_PROPERTY_ID)
    KRATOS_DEFINE_VARIABLE(double,  HTC)
    KRATOS_DEFINE_VARIABLE(int, REF_ID)
    KRATOS_DEFINE_VARIABLE(double, PARTICLE_RADIUS)
    KRATOS_DEFINE_VARIABLE(double, POSETIVE_DISTANCE)
    KRATOS_DEFINE_VARIABLE(double, NAGATIVE_DISTANCE)
    KRATOS_DEFINE_VARIABLE(bool, IS_ESCAPED)
    KRATOS_DEFINE_VARIABLE(int, IS_SOLIDIFIED)  
    KRATOS_DEFINE_VARIABLE(double, SOLIDFRACTION ) 
    KRATOS_DEFINE_VARIABLE(double, SOLIDIF_TIME  )
    KRATOS_DEFINE_VARIABLE(double, SOLIDIF_MODULUS  )
    KRATOS_DEFINE_VARIABLE(double, FILLTIME  )
    KRATOS_DEFINE_VARIABLE(double, MACRO_POROSITY  )
    KRATOS_DEFINE_VARIABLE(double, SHRINKAGE_POROSITY  )
    KRATOS_DEFINE_VARIABLE(double, MAX_VEL  )
    KRATOS_DEFINE_VARIABLE(int, IS_GRAVITY_FILLING)
    KRATOS_DEFINE_VARIABLE(double,VOLUME_FRACTION ) 
    KRATOS_DEFINE_VARIABLE(double,KAPPA ) 
    KRATOS_DEFINE_VARIABLE(double,EPSILON )
    KRATOS_DEFINE_VARIABLE(double,SHRINKAGE_POROSITY_US)
    KRATOS_DEFINE_VARIABLE(double,SOLIDIF_MODULUS_US)
    KRATOS_DEFINE_VARIABLE(double,TEMPERATURES_US)
    KRATOS_DEFINE_VARIABLE(double,FRONT_MEETING)



}  // namespace Kratos.



// Resotre the default defines
#undef KRATOS_DEFINE_VARIABLE
#undef KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS
#define KRATOS_DEFINE_VARIABLE KRATOS_DEFINE_VARIABLE_NO_DLL
#define KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS_NO_DLL

#endif // KRATOS_C2C_VARIABLES_H_INCLUDED  defined 
