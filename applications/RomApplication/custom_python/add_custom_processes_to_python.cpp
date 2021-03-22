// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Altug Emiroglu, http://github.com/emiroglu
//

// System includes

// External includes


// Project includes
#include "custom_python/add_custom_processes_to_python.h"
#include "rom_application_variables.h"

//Processes
#include "custom_processes/postprocess_rom_basis_process.h"

namespace Kratos {
namespace Python {

void  AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    /// Processes
    py::class_<PostprocessRomBasisProcess, PostprocessRomBasisProcess::Pointer, Process>(m,"PostprocessRomBasisProcess")
        .def(py::init<ModelPart&, Parameters>());

}

}  // namespace Python.
} // Namespace Kratos

