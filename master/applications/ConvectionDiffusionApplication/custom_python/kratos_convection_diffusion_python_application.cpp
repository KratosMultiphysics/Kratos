// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes

// Project includes
#include "includes/define_python.h"
#include "convection_diffusion_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_response_functions_to_python.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;

PYBIND11_MODULE(KratosConvectionDiffusionApplication,m)
{

    py::class_<KratosConvectionDiffusionApplication,
           KratosConvectionDiffusionApplication::Pointer,
           KratosApplication >(m,"KratosConvectionDiffusionApplication")
           .def(py::init<>())
           ;
    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomResponseFunctionsToPython(m);

    // Registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AMBIENT_TEMPERATURE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AUX_FLUX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AUX_TEMPERATURE)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BFECC_ERROR)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BFECC_ERROR_1)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONVECTION_COEFFICIENT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DELTA_SCALAR1)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MEAN_SIZE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MELT_TEMPERATURE_1)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MELT_TEMPERATURE_2)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EMBEDDED_SCALAR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PROJECTED_SCALAR1)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TRANSFER_COEFFICIENT)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ADJOINT_HEAT_TRANSFER)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_PROJECTION)

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,CONVECTION_VELOCITY)
}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
