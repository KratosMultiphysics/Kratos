// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
//         -        Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//         -        Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
//                  in the documentation and/or other materials provided with the distribution.
//         -        All advertising materials mentioning features or use of this software must display the following acknowledgement:
//                         This product includes Kratos Multi-Physics technology.
//         -        Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "iga_structural_mechanics_application.h"
#include "iga_structural_mechanics_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos
{

namespace Python
{

  using namespace boost::python;



  BOOST_PYTHON_MODULE(KratosIGAStructuralMechanicsApplication)
  {

	  class_<KratosIGAStructuralMechanicsApplication,
			  KratosIGAStructuralMechanicsApplication::Pointer,
			  bases<KratosApplication>, boost::noncopyable >("KratosIGAStructuralMechanicsApplication")
			;

	AddCustomStrategiesToPython();
	AddCustomUtilitiesToPython();

	//registering variables in python

	KRATOS_REGISTER_IN_PYTHON_VARIABLE(INTEGRATION_WEIGHT)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(SHAPE_FUNCTION_VALUES)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(SHAPE_FUNCTION_LOCAL_DERIVATIVES)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(SHAPE_FUNCTION_LOCAL_DERIVATIVES_MASTER)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE)

	KRATOS_REGISTER_IN_PYTHON_VARIABLE(TANGENTS)

	KRATOS_REGISTER_IN_PYTHON_VARIABLE(PENALTY_FACTOR)
	
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(DISPLACEMENT_ROTATION_FIX)

	KRATOS_REGISTER_IN_PYTHON_VARIABLE(LOAD_TYPE)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(DISTRIBUTED_LOAD_FACTOR)

  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
