//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//                   Riccardo Rossi
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "constitutive_laws_application.h"
#include "constitutive_laws_application_variables.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosConstitutiveLawsApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosConstitutiveLawsApplication,
        KratosConstitutiveLawsApplication::Pointer,
        KratosApplication>(m, "KratosConstitutiveLawsApplication")
        .def(py::init<>())
        ;

}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
