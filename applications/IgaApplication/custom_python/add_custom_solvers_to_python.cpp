//
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Manuel Messmer
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_solvers_to_python.h"
#include "custom_solvers/additive_schwarz_preconditioner.h"

#include "iga_application_variables.h"


namespace Kratos {
namespace Python {

void AddCustomSolversToPython(
    pybind11::module& m)
{
    namespace py = pybind11;

}

} // namespace Python
} // Namespace Kratos
