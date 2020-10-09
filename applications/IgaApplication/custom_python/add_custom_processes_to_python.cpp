//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_processes_to_python.h"

#include "custom_processes/output_quadrature_domain_process.h"

#include "iga_application_variables.h"


namespace Kratos {
namespace Python {

void AddCustomProcessesToPython(
    pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<OutputQuadratureDomainProcess, OutputQuadratureDomainProcess::Pointer, Process>(m, "OutputQuadratureDomainProcess")
        .def(py::init<Model&, Parameters >())
        ;
}

} // namespace Python
} // Namespace Kratos
