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
#include "includes/kratos_export_dll.h"

#undef KRATOS_DEFINE_VARIABLE
#undef KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS
#define KRATOS_DEFINE_VARIABLE KRATOS_DEFINE_VARIABLE_DLL
#define KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS_DLL

//TODO: move to the FluidDynamics application
namespace Kratos
{

    // Useful variables
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( ADVPROJ )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( CONV_PROJ )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PRESS_PROJ )
    KRATOS_DEFINE_VARIABLE( double, DIVPROJ )
    KRATOS_DEFINE_VARIABLE( double, PRESSURE_OLD_IT )
    KRATOS_DEFINE_VARIABLE( double, C_SMAGORINSKY )
    KRATOS_DEFINE_VARIABLE( double, MOLECULAR_VISCOSITY )
    KRATOS_DEFINE_VARIABLE( double, TURBULENT_VISCOSITY )
    KRATOS_DEFINE_VARIABLE( double, Y_WALL)
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

} // namespace Kratos


// Resotre the default defines
#undef KRATOS_DEFINE_VARIABLE
#undef KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS
#define KRATOS_DEFINE_VARIABLE KRATOS_DEFINE_VARIABLE_NO_DLL
#define KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS_NO_DLL

#endif // KRATOS_CFD_VARIABLES_H_INCLUDED  defined
