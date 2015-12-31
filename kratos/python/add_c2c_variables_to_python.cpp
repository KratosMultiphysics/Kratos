// Kratos Multi-Physics
//
// Copyright (c) 2016, Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
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


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "containers/data_value_container.h"
//#include "containers/hash_data_value_container.h"
// #include "containers/variables_list_data_value_container.h"
// #include "containers/fix_data_value_container.h"
// #include "containers/vector_component_adaptor.h"
#include "containers/flags.h"
//#include "containers/all_variables_data_value_container.h"
// #include "includes/kratos_flags.h"
#include "includes/c2c_variables.h"
// #include "includes/constitutive_law.h"
#include "python/add_c2c_variables_to_python.h"
// #include "python/vector_python_interface.h"
// #include "python/vector_scalar_operator_python.h"
// #include "python/vector_vector_operator_python.h"
// #include "python/bounded_vector_python_interface.h"

// #include "includes/convection_diffusion_settings.h"
// #include "includes/radiation_settings.h"
#include "utilities/timer.h"



#ifdef KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION
#undef KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION
#endif
#define KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION(flag) \
 scope().attr(#flag) = boost::ref(flag)      \
 
#ifdef KRATOS_REGISTER_IN_PYTHON_FLAG
#undef KRATOS_REGISTER_IN_PYTHON_FLAG
#endif
#define KRATOS_REGISTER_IN_PYTHON_FLAG(flag) \
    KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION(flag);   \
    KRATOS_REGISTER_IN_PYTHON_FLAG_IMPLEMENTATION(NOT_##flag)

#ifdef KRATOS_REGISTER_IN_PYTHON_VARIABLE
#undef KRATOS_REGISTER_IN_PYTHON_VARIABLE
#endif
#define KRATOS_REGISTER_IN_PYTHON_VARIABLE(variable) \
    scope().attr(#variable) = boost::ref(variable);


#ifdef KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS
#undef KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS
#endif
#define KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(name) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(name) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(name##_X) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(name##_Y) \
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(name##_Z)

namespace Kratos
{
//KRATOS_CREATE_FLAG(STRUCTURE,   63);

namespace Python
{
using namespace boost::python;

void  AddC2CVariablesToPython()
{
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  SOLID_TEMPERATURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  FLUID_TEMPERATURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  AVERAGE_TEMPERATURE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  INLET_TEMPERATURE)
//         KRATOS_REGISTER_IN_PYTHON_VARIABLE(  AMBIENT_TEMPERATURE )
//         KRATOS_REGISTER_IN_PYTHON_VARIABLE(  Y_WALL)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  COUNTER )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  DISTANCE_CORRECTION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  COMPUTED_DISTANCE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  MATERIAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  MOULD_AVERAGE_TEMPERATURE)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(   LAST_AIR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(   PRESSURES)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(   TEMPERATURES)
    // KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(VELOCITIES)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(VELOCITIES)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(   MOULD_DENSITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(   MOULD_SPECIFIC_HEAT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(   MOULD_THICKNESS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(   MOULD_SFACT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(   MOULD_VFACT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(   MOULD_CONDUCTIVITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(   MOULD_HTC_ENVIRONMENT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(   MOULD_TEMPERATURE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(   MOULD_INNER_TEMPERATURE)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  NODE_PROPERTY_ID)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(   HTC)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  REF_ID)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  PARTICLE_RADIUS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  POSETIVE_DISTANCE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  NAGATIVE_DISTANCE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_ESCAPED)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  IS_SOLIDIFIED)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  SOLIDFRACTION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  SOLIDIF_TIME  )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  SOLIDIF_MODULUS  )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  FILLTIME  )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  MACRO_POROSITY  )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  SHRINKAGE_POROSITY  )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  MAX_VEL  )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(  IS_GRAVITY_FILLING)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( VOLUME_FRACTION )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( KAPPA )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( EPSILON )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SHRINKAGE_POROSITY_US)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SOLIDIF_MODULUS_US)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( TEMPERATURES_US)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( FRONT_MEETING)



}
}  // namespace Python.
} // Namespace Kratos

#undef KRATOS_REGISTER_IN_PYTHON_FLAG
#undef KRATOS_REGISTER_IN_PYTHON_VARIABLE
#undef KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS
