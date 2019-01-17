//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"

#include "custom_python/add_trilinos_space_to_python.h"
#include "custom_python/add_trilinos_convergence_criterias_to_python.h"
#include "custom_python/add_trilinos_schemes_to_python.h"
#include "custom_python/add_trilinos_linear_solvers_to_python.h"
#include "custom_python/add_trilinos_processes_to_python.h"
#include "custom_python/add_trilinos_strategies_to_python.h"
#include "custom_python/add_custom_io_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_zoltan_processes_to_python.h"

// Project includes
#include "trilinos_application.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

PYBIND11_MODULE(KratosTrilinosApplication,m)
{

    py::class_<KratosTrilinosApplication,
        KratosTrilinosApplication::Pointer,
        KratosApplication > (m,"KratosTrilinosApplication")
        .def(py::init<>())
        ;

    AddBasicOperations(m);
    AddConvergenceCriterias(m);
    AddSchemes(m);
    AddLinearSolvers(m);
    AddProcesses(m);
    AddStrategies(m);
    AddCustomIOToPython(m);
    AddCustomUtilitiesToPython(m);
    AddZoltanProcessesToPython(m);

    //registering variables in python

}


} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
