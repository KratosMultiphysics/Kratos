// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
// | | | (_)        |  \/  |         | |
// | | | |_ ___  ___| .  . | ___   __| |
// | | | | / __|/ __| |\/| |/ _ \ / _` |
// \ \_/ / \__ \ (__| |  | | (_) | (_| |
//  \___/|_|___/\___\_|  |_/\___/ \__,_|  APPLICATION
//                                      
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aniol Sala
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "viscosity_modulator_application.h"
#include "viscosity_modulator_application_variables.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

PYBIND11_MODULE(KratosViscosityModulatorApplication,m)
{

    class_<KratosViscosityModulatorApplication,
            KratosViscosityModulatorApplication::Pointer,
            KratosApplication>(m,"KratosViscosityModulatorApplication")
            .def(init<>())
            ;

    AddCustomProcessesToPython(m);
    AddCustomUtilitiesToPython(m);

    //registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AUX_SCALAR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHOCK_CAPTURING_INTENSITY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, USE_ANISOTROPIC_DISC_CAPTURING)


}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
