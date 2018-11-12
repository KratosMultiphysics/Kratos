//
// Author: Salva Latorre latorre@cimne.upc.edu
//

#if defined(KRATOS_PYTHON)

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"
#include "NvidiaFlex_application.h"
#include "custom_python/add_custom_utilities_to_python.h"

namespace Kratos {

    namespace Python {

        namespace py = pybind11;

        PYBIND11_MODULE(KratosNvidiaFlexApplication,m) {

            py::class_<KratosNvidiaFlexApplication, KratosNvidiaFlexApplication::Pointer, KratosApplication>(m, "KratosNvidiaFlexApplication")
                .def(py::init<>());
                ;

            AddCustomUtilitiesToPython(m);

            //KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MASS_FLOW)
            //KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, LINEAR_VELOCITY)

        }

    }  // namespace Python

}  // namespace Kratos

#endif // KRATOS_PYTHON defined
