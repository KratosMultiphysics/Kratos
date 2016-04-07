// KratosFluidDynamicsApplication
// 
// Copyright (c) 2015, Pooyan Dadvand, Riccardo Rossi,Jordi Cotela CIMNE (International Center for Numerical Methods in Engineering)
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


//
//   Project Name:        Kratos
//   Last Modified by:    $Author: jcotela $
//   Date:                $Date: 2010-11-11 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_FLUID_DYNAMICS_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_FLUID_DYNAMICS_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/dem_variables.h"  //TODO: must be removed eventually
#include "includes/legacy_structural_app_vars.h"  //TODO: must be removed eventually

namespace Kratos
{
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, int,PATCH_INDEX)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,TAUONE)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,TAUTWO)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,PRESSURE_MASSMATRIX_COEFFICIENT)
//KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,Y_WALL)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,SUBSCALE_PRESSURE)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,C_DES)
//    KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,C_SMAGORINSKY)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( FLUID_DYNAMICS_APPLICATION, SUBSCALE_VELOCITY)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( FLUID_DYNAMICS_APPLICATION, VORTICITY)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( FLUID_DYNAMICS_APPLICATION, COARSE_VELOCITY)

// Non-Newtonian constitutive relations
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, REGULARIZATION_COEFFICIENT)

// To be removed
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, BINGHAM_SMOOTHER)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double, GEL_STRENGTH)

// Q-Criterion (for vortex visualization)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,Q_VALUE)

// Vorticity
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, double,VORTICITY_MAGNITUDE)
KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(FLUID_DYNAMICS_APPLICATION, RECOVERED_PRESSURE_GRADIENT)
KRATOS_DEFINE_APPLICATION_VARIABLE( FLUID_DYNAMICS_APPLICATION, Vector,NODAL_WEIGHTS)
}

#endif	/* KRATOS_FLUID_DYNAMICS_APPLICATION_VARIABLES_H_INCLUDED */
