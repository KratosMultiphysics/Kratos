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



// This define must be HERE
#define DKRATOS_EXPORT_INTERFACE_2 1

// System includes
#include <string>
#include <iostream>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/cfd_variables.h"
#include "includes/kernel.h"

namespace Kratos
{
// Useful variables
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( ADVPROJ );
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( CONV_PROJ );
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( PRESS_PROJ );
KRATOS_CREATE_VARIABLE( double, DIVPROJ );
KRATOS_CREATE_VARIABLE( double, PRESSURE_OLD_IT );
KRATOS_CREATE_VARIABLE( double, C_SMAGORINSKY );
KRATOS_CREATE_VARIABLE( double, MOLECULAR_VISCOSITY );
KRATOS_CREATE_VARIABLE( double, TURBULENT_VISCOSITY );
KRATOS_CREATE_VARIABLE( double, Y_WALL);
KRATOS_CREATE_VARIABLE( int, OSS_SWITCH );

// Legacy variables
KRATOS_CREATE_VARIABLE( double, DYNAMIC_TAU );
KRATOS_CREATE_VARIABLE( double, DYNAMIC_VISCOSITY);
KRATOS_CREATE_VARIABLE( double, EFFECTIVE_VISCOSITY );
KRATOS_CREATE_VARIABLE( double, KINEMATIC_VISCOSITY);
KRATOS_CREATE_VARIABLE( double, THAWONE );
KRATOS_CREATE_VARIABLE( double, THAWTWO );
KRATOS_CREATE_VARIABLE( double, M );

void KratosApplication::RegisterCFDVariables()
{
    // Useful variables
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( ADVPROJ );
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( CONV_PROJ );
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( PRESS_PROJ );
    KRATOS_REGISTER_VARIABLE( DIVPROJ );
    KRATOS_REGISTER_VARIABLE( PRESSURE_OLD_IT );
    KRATOS_REGISTER_VARIABLE( C_SMAGORINSKY );
    KRATOS_REGISTER_VARIABLE( MOLECULAR_VISCOSITY );
    KRATOS_REGISTER_VARIABLE( TURBULENT_VISCOSITY );
    KRATOS_REGISTER_VARIABLE( Y_WALL);
    KRATOS_REGISTER_VARIABLE( OSS_SWITCH );

    // Legacy variables
    KRATOS_REGISTER_VARIABLE( DYNAMIC_TAU );
    KRATOS_REGISTER_VARIABLE( DYNAMIC_VISCOSITY);
    KRATOS_REGISTER_VARIABLE( EFFECTIVE_VISCOSITY );
    KRATOS_REGISTER_VARIABLE( KINEMATIC_VISCOSITY);
    KRATOS_REGISTER_VARIABLE( THAWONE );
    KRATOS_REGISTER_VARIABLE( THAWTWO );
    KRATOS_REGISTER_VARIABLE( M );
}


}  // namespace Kratos.

// This define must be HERE
#undef DKRATOS_EXPORT_INTERFACE_2



