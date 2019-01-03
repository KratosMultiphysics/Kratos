//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

// System includes

#if defined(KRATOS_PYTHON)

// External includes
#include "pybind11/pybind11.h"


// Project includes
#include "includes/define_python.h"
#include "stabilized_cfd_application.h"
#include "stabilized_cfd_application_variables.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos
{

namespace Python
{

PYBIND11_MODULE(KratosStabilizedCFDApplication,m)
{
  namespace py = pybind11;
  py::class_<KratosStabilizedCFDApplication,
    KratosStabilizedCFDApplication::Pointer,
    KratosApplication >(m,"KratosStabilizedCFDApplication")
    .def(py::init<>())
    ;

	AddCustomUtilitiesToPython(m);

	//registering variables in python
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FIC_BETA )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DIRECTIONAL_BETA )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RECORDED_STEPS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MEAN_KINETIC_ENERGY )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MEAN_VELOCITY )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MEAN_PRESSURE )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VELOCITY_COVARIANCES )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TURBULENCE_STATISTICS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TRACE_XI )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DIV_XI )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MOMENTUM_PROJECTION )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MOMENTUM_PROJECTION_RHS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MASS_PROJECTION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MASS_PROJECTION_RHS )

}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
