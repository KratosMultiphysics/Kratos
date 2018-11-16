/*
==============================================================================
KratosPFEM2Application
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
//   Last modified by:    $Author:  $
//   Date:                $Date: $
//   Revision:            $Revision: 1.3 $
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
//#include <boost/python.hpp>


// Project includes
//#include "includes/define.h"
#include "includes/define_python.h"
#include "pfem_2_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos
{

namespace Python
{

  namespace py = pybind11;


  PYBIND11_MODULE(KratosPFEM2Application,m)
  {

    py::class_<KratosPFEM2Application, KratosPFEM2Application::Pointer,KratosApplication>(m,"KratosPFEM2Application").def(py::init<>());

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,G_VALUE)
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,PREVIOUS_ITERATION_PRESSURE)
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,FIRST_ITERATION_PRESSURE)
      KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,PRESS_PROJ_NO_RO)
      KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,DELTA_VELOCITY)
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,MEAN_SIZE)
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SPECIFIC_HEAT_CAPACITY_WATER)
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SPECIFIC_HEAT_CAPACITY_AIR)
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,DELTA_TEMPERATURE)
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,USEFUL_ELEMENT_FOR_COMBUSTION)
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,AVAILABLE_AIR_VOLUME)
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,AVAILABLE_UNBURNED_AIR_VOLUME)
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,OXYGEN_FRACTION)
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,CORRECTED_DISTANCE)
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SOLID_PRESSURE)
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SOLID_YP)
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,WATER_DISTANCE)
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,VOLUMETRIC_STRAIN)
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,ELASTIC_PRESSURE)
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,USE_PRESS_PROJ)
      KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,WATER_VELOCITY)
      KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,WATER_MESH_VELOCITY)

      KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,PROJECTED_VELOCITY)
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,VOLUME_CORRECTION )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,INLET_VELOCITY)



	//KRATOS_REGISTER_IN_PYTHON_VARIABLE(NEIGHBOUR_ELEMENTS)
	//KRATOS_REGISTER_IN_PYTHON_VARIABLE(IS_AIR);

  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
