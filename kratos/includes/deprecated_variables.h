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


#if !defined(KRATOS_DEPRECATED_VARIABLES_H_INCLUDED )
#define  KRATOS_DEPRECATED_VARIABLES_H_INCLUDED



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

namespace Kratos
{

    //Define Variables by type:

    //bools

    //for Structural application:
    KRATOS_DEFINE_VARIABLE( bool, IS_INACTIVE )

    //for Level Set application:
    KRATOS_DEFINE_VARIABLE( bool, IS_DUPLICATED )
    KRATOS_DEFINE_VARIABLE( bool, SPLIT_ELEMENT )
    KRATOS_DEFINE_VARIABLE( bool, SPLIT_NODAL )


    //for PFEM fluids application:
    KRATOS_DEFINE_VARIABLE( int, IS_JACK_LINK )
    KRATOS_DEFINE_VARIABLE( int, IMPOSED_PRESSURE )
    KRATOS_DEFINE_VARIABLE( int, IMPOSED_VELOCITY_X )
    KRATOS_DEFINE_VARIABLE( int, IMPOSED_VELOCITY_Y )
    KRATOS_DEFINE_VARIABLE( int, IMPOSED_VELOCITY_Z )
    KRATOS_DEFINE_VARIABLE( int, IMPOSED_ANGULAR_VELOCITY_X )
    KRATOS_DEFINE_VARIABLE( int, IMPOSED_ANGULAR_VELOCITY_Y )
    KRATOS_DEFINE_VARIABLE( int, IMPOSED_ANGULAR_VELOCITY_Z )

    // For the DEM Application:
    KRATOS_DEFINE_VARIABLE( double, IMPOSED_VELOCITY_X_VALUE )
    KRATOS_DEFINE_VARIABLE( double, IMPOSED_VELOCITY_Y_VALUE )
    KRATOS_DEFINE_VARIABLE( double, IMPOSED_VELOCITY_Z_VALUE )
    KRATOS_DEFINE_VARIABLE( double, IMPOSED_ANGULAR_VELOCITY_X_VALUE )
    KRATOS_DEFINE_VARIABLE( double, IMPOSED_ANGULAR_VELOCITY_Y_VALUE )
    KRATOS_DEFINE_VARIABLE( double, IMPOSED_ANGULAR_VELOCITY_Z_VALUE )

    KRATOS_DEFINE_VARIABLE( double, IS_INLET )
    KRATOS_DEFINE_VARIABLE( double, IS_INTERFACE )
    KRATOS_DEFINE_VARIABLE( double, IS_VISITED )
    KRATOS_DEFINE_VARIABLE( double, IS_EROSIONABLE )

    KRATOS_DEFINE_VARIABLE( double, IS_STRUCTURE )
    KRATOS_DEFINE_VARIABLE( double, IS_POROUS )
    KRATOS_DEFINE_VARIABLE( double, IS_WATER )
    KRATOS_DEFINE_VARIABLE( double, IS_FLUID )
    KRATOS_DEFINE_VARIABLE( double, IS_BOUNDARY )
    KRATOS_DEFINE_VARIABLE( double, IS_FREE_SURFACE )
    KRATOS_DEFINE_VARIABLE( double, IS_AIR_EXIT )
    KRATOS_DEFINE_VARIABLE( double, IS_LAGRANGIAN_INLET )
    KRATOS_DEFINE_VARIABLE( double, IS_WATER_ELEMENT )


    KRATOS_DEFINE_VARIABLE( double, IS_BURN )
    KRATOS_DEFINE_VARIABLE( double, IS_DRIPPING )
    KRATOS_DEFINE_VARIABLE( double, IS_PERMANENT )
    KRATOS_DEFINE_VARIABLE( double, IS_WALL )

    KRATOS_DEFINE_VARIABLE( double, Ypr ) //var name does not follow standard
    KRATOS_DEFINE_VARIABLE( double, Yox )
    KRATOS_DEFINE_VARIABLE( double, Yfuel )
    KRATOS_DEFINE_VARIABLE( double, Hfuel )
    KRATOS_DEFINE_VARIABLE( double, Hpr )
    KRATOS_DEFINE_VARIABLE( double, Hpr1 )
    KRATOS_DEFINE_VARIABLE( double, Hox )

    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_1 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_2 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_3 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_4 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_5 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_6 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_7 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_8 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_9 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_10 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_11 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_12 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_13 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_14 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_15 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_16 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_17 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_18 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_19 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_20 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_21 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_22 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_23 )
    KRATOS_DEFINE_VARIABLE( double, RADIATIVE_INTENSITY_24 )

    KRATOS_DEFINE_VARIABLE( double, rhoD )
    KRATOS_DEFINE_VARIABLE( double, xi )
    KRATOS_DEFINE_VARIABLE( double, a )
    KRATOS_DEFINE_VARIABLE( double, b )


    KRATOS_DEFINE_VARIABLE( double, IS_SLIP )

    //for Level Set application:
    KRATOS_DEFINE_VARIABLE( double, IS_DIVIDED )

    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( xi_c )

}  // namespace Kratos.

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT

#endif // KRATOS_DEPRECATED_VARIABLES_H_INCLUDED  defined
