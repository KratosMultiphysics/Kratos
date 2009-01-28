/*
==============================================================================
KratosULFApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Pawel Ryzhakov
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
//   Last modified by:    $Author: anonymous $
//   Date:                $Date: 2008-10-23 12:50:01 $
//   Revision:            $Revision: 1.6 $
//
//

// System includes 

#if defined(KRATOS_PYTHON)
// External includes 
#include <boost/python.hpp>


// Project includes 
#include "includes/define.h"
#include "ULF_application.h"
#include "custom_python/add_custom_io_to_python.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_processes_to_python.h"

 
namespace Kratos
{

namespace Python
{

  using namespace boost::python;


  
  BOOST_PYTHON_MODULE(KratosULFApplication)
  {

	  class_<KratosULFApplication, 
			  KratosULFApplication::Pointer, 
			  bases<KratosApplication>, boost::noncopyable >("KratosULFApplication")
			;
		AddCustomUtilitiesToPython();
		AddCustomIOToPython();
		AddCustomStrategiesToPython();
		AddProcessesToPython();

		//registering variables in python
/*		KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_H)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(IS_STRUCTURE)
		KRATOS_REGISTER_VKRATOS_REGISTER_IN_PYTHON_VARIABLEARIABLE(IS_FLUID)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(IS_BOUNDARY)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(IS_FREE_SURFACE)
*/		
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(IS_LAGRANGIAN_INLET)

		KRATOS_REGISTER_IN_PYTHON_VARIABLE(PRESSURE_FORCE)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(FRACT_VEL)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(VAUX)
		
  }
  
  
}  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
