/*
==============================================================================
KratosSwimmingDEMApplication
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
//   Last modified by:    $Author:  G. Casas$
//   Date:                $Date: $
//   Revision:            $Revision: 1.3 $
//
//
#if defined(KRATOS_PYTHON)

// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/define_python.h"
#include "swimming_DEM_application.h"
#include "swimming_dem_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;

PYBIND11_MODULE(KratosSwimmingDEMApplication, m)
{
    py::class_<KratosSwimmingDEMApplication, KratosSwimmingDEMApplication::Pointer, KratosApplication>(m, "KratosSwimmingDEMApplication")
        .def(py::init<>())
        ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);
    AddCustomConstitutiveLawsToPython(m);

    //registering variables in python

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DISPLACEMENT_OLD)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, AVERAGED_FLUID_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VECTORIAL_ERROR)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EXACT_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EXACT_PRESSURE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_ERROR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ERROR_X )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ERROR_Y )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ERROR_Z )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ERROR_P )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PERMEABILITY_1_DAY )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SDEM_HYDRODYNAMIC_INTERACTION_LAW_NAME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SDEM_BUOYANCY_LAW_NAME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SDEM_DRAG_LAW_NAME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SDEM_INVISCID_FORCE_LAW_NAME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SDEM_HISTORY_FORCE_LAW_NAME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SDEM_VORTICITY_LIFT_LAW_NAME )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SDEM_STEADY_VISCOUS_TORQUE_LAW_NAME )

}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
