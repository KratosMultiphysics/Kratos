/*
==============================================================================
KratosR1ElectrostaticApplication 
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
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2007-10-25 10:14:31 $
//   Revision:            $Revision: 1.4 $
//
//

// System includes 

#if defined(KRATOS_PYTHON)
// External includes 
#include <boost/python.hpp>


// Project includes 
#include "includes/define.h"
#include "kElectrostatic.h"
//#include "custom_python/add_custom_strategies_to_python.h"
//#include "custom_python/add_custom_utilities_to_python.h"
//#include "custom_python/add_custom_io_to_python.h"

 
namespace Kratos
{

namespace Python
{

  using namespace boost::python;


  
  BOOST_PYTHON_MODULE(KratosR1ElectrostaticApplication)
  {

	  class_<KratosR1ElectrostaticApplication, 
			  KratosR1ElectrostaticApplication::Pointer, 
			  bases<KratosApplication>, boost::noncopyable >("KratosR1ElectrostaticApplication")
			;

		//AddCustomStrategiesToPython();
		//AddCustomUtilitiesToPython();
//	    AddCustomIOToPython();

		//registering variables in python

		KRATOS_REGISTER_IN_PYTHON_VARIABLE( BDF_COEFFICIENTS );
		//KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_AREA)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE( AUX_INDEX)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE( TEMP_CONV_PROJ)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONDUCTIVITY)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(SPECIFIC_HEAT)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(HEAT_FLUX)	

		KRATOS_REGISTER_IN_PYTHON_VARIABLE(AMBIENT_TEMPERATURE)	
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(EMISSIVITY)	
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONVECTION_COEFFICIENT)	
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(FACE_HEAT_FLUX)	
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(ELECTRIC_POTENTIAL)

		// for electromagnetic applications
		// for kElectrostatic application
		KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(ELECTRICAL_PERMITTIVITY)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(ELECTROSTATIC_POTENTIAL)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(ELECTROSTATIC_POINT_CHARGE)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(ELECTROSTATIC_SURFACE_CHARGE)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(ELECTROSTATIC_VOLUME_CHARGE)
		KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(ELECTRIC_FIELD)
		KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(ELECTRIC_DISPLACEMENT_FIELD)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(INFINIT_COEFFICIENT)

		KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)
			 
  }
  
  
}  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
