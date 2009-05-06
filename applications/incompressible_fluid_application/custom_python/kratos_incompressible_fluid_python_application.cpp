/*
==============================================================================
KratosIncompressibleFluidApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: kazem $
//   Date:                $Date: 2009-01-15 18:46:07 $
//   Revision:            $Revision: 1.15 $
//
//

// System includes 

#if defined(KRATOS_PYTHON)
// External includes 
#include <boost/python.hpp>


// Project includes 
#include "includes/define.h"
#include "incompressible_fluid_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_io_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_processes_to_python.h"
 
namespace Kratos
{

namespace Python
{

  using namespace boost::python;


  
  BOOST_PYTHON_MODULE(KratosIncompressibleFluidApplication)
  {

	  class_<KratosIncompressibleFluidApplication, 
			  KratosIncompressibleFluidApplication::Pointer, 
			  bases<KratosApplication>, boost::noncopyable >("KratosIncompressibleFluidApplication")
			;
		AddCustomStrategiesToPython();
		AddCustomUtilitiesToPython();
		AddCustomIOToPython();
		AddProcessesToPython();

		//registering variables in python
		KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(FRACT_VEL)
	//	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(PRESS_PROJ)
		KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(CONV_PROJ)
		KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(ADVPROJ)
		KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(REACTION)


		KRATOS_REGISTER_IN_PYTHON_VARIABLE(  MACH_NUMBER )
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(REACTION_WATER_PRESSURE)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(  PRESSURE_COEFFICIENT )
		KRATOS_REGISTER_IN_PYTHON_VARIABLE( PRESSURE_OLD_IT )
		KRATOS_REGISTER_IN_PYTHON_VARIABLE( BDF_COEFFICIENTS );
		KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_MASS)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE( AUX_INDEX)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE( EXTERNAL_PRESSURE)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE( DIAMETER)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE( PERMEABILITY_INV)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(  WATER_PRESSURE )
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(  AIR_PRESSURE )

		KRATOS_REGISTER_IN_PYTHON_VARIABLE(  AIR_SOUND_VELOCITY )
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(  SOUND_VELOCITY )

		KRATOS_REGISTER_IN_PYTHON_VARIABLE( DIVPROJ)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE( THAWONE)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE( THAWTWO)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE( OSS_SWITCH)

		KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(RHS_VECTOR)
		KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(AUX_VECTOR)
		 
  }
  
  
}  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
