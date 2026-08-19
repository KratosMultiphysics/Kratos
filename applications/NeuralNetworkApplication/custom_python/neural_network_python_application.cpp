//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Rishith Ellath Meethal (https://github.com/rishithellathmeethal)
//                   Daniel Andrés Arcones https://github.com/danielandresarcones
//

// System includes

#if defined(KRATOS_PYTHON)

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"
#include "neural_network_application.h"
#include "neural_network_application_variables.h"

namespace Kratos::Python {

PYBIND11_MODULE(KratosNeuralNetworkApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosNeuralNetworkApplication,
        KratosNeuralNetworkApplication::Pointer,
        KratosApplication>(m, "KratosNeuralNetworkApplication")
        .def(py::init<>())
        ;
}

} // namespace Kratos::Python.

#endif // KRATOS_PYTHON defined
