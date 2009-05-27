/*
==============================================================================
KratosConvectionDiffusionApplication 
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
//   Last modified by:    $Author: anonymous $
//   Date:                $Date: 2008-12-15 16:03:54 $
//   Revision:            $Revision: 1.5 $
//
//

// System includes 

#if defined(KRATOS_PYTHON)
// External includes 
#include <boost/python.hpp>


// Project includes 
#include "includes/define.h"
#include "convection_diffusion_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"

 
namespace Kratos
{

namespace Python
{

  using namespace boost::python;


  
  BOOST_PYTHON_MODULE(KratosConvectionDiffusionApplication)
  {

	  class_<KratosConvectionDiffusionApplication, 
			  KratosConvectionDiffusionApplication::Pointer, 
			  bases<KratosApplication>, boost::noncopyable >("KratosConvectionDiffusionApplication")
			;
		AddCustomStrategiesToPython();
		AddCustomUtilitiesToPython();

		//registering variables in python

//		KRATOS_REGISTER_IN_PYTHON_VARIABLE( BDF_COEFFICIENTS );
		//KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_AREA)
//		KRATOS_REGISTER_IN_PYTHON_VARIABLE( AUX_INDEX)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE( TEMP_CONV_PROJ)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONDUCTIVITY)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(SPECIFIC_HEAT)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(HEAT_FLUX)	

		KRATOS_REGISTER_IN_PYTHON_VARIABLE(AMBIENT_TEMPERATURE)	
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(EMISSIVITY)	
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONVECTION_COEFFICIENT)	
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(FACE_HEAT_FLUX)	

		KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)
		//ADded by Annelie and Pavel
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(ENTHALPY)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(LATENT_HEAT)	
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(MELT_TEMPERATURE_1)
		KRATOS_REGISTER_IN_PYTHON_VARIABLE(MELT_TEMPERATURE_2)
			 
  }
  
  
}  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
