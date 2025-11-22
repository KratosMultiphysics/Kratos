//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//


// System includes

// External includes
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_utilities_to_python.h"

// Application includes
#include "custom_utilities/method_utilities.h"
#include "custom_utilities/norms.h"

namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto norms = m.def_submodule("Norms");

    // adding norms
    py::class_<Norms::L2, Norms::L2::Pointer>(norms, "L2")
        .def(py::init<>())
        .def("Evaluate", &Norms::L2::Evaluate<int>)
        .def("Evaluate", &Norms::L2::Evaluate<double>)
        .def("Evaluate", &Norms::L2::Evaluate<array_1d<double, 3>>)
        .def("Evaluate", &Norms::L2::Evaluate<array_1d<double, 4>>)
        .def("Evaluate", &Norms::L2::Evaluate<array_1d<double, 6>>)
        .def("Evaluate", &Norms::L2::Evaluate<array_1d<double, 9>>)
        .def("Evaluate", &Norms::L2::Evaluate<Vector>)
        .def("Evaluate", &Norms::L2::Evaluate<Matrix>)
        .def("__str__", &Norms::L2::Info)
        ;
    py::class_<Norms::Infinity, Norms::Infinity::Pointer>(norms, "Infinity")
        .def(py::init<>())
        .def("Evaluate", &Norms::Infinity::Evaluate<int>)
        .def("Evaluate", &Norms::Infinity::Evaluate<double>)
        .def("Evaluate", &Norms::Infinity::Evaluate<array_1d<double, 3>>)
        .def("Evaluate", &Norms::Infinity::Evaluate<array_1d<double, 4>>)
        .def("Evaluate", &Norms::Infinity::Evaluate<array_1d<double, 6>>)
        .def("Evaluate", &Norms::Infinity::Evaluate<array_1d<double, 9>>)
        .def("Evaluate", &Norms::Infinity::Evaluate<Vector>)
        .def("Evaluate", &Norms::Infinity::Evaluate<Matrix>)
        .def("__str__", &Norms::Infinity::Info)
        ;

    py::class_<Norms::P, Norms::P::Pointer>(norms, "P")
        .def(py::init<const double>(), py::arg("p_coefficient"))
        .def("Evaluate", &Norms::P::Evaluate<array_1d<double, 3>>)
        .def("Evaluate", &Norms::P::Evaluate<array_1d<double, 4>>)
        .def("Evaluate", &Norms::P::Evaluate<array_1d<double, 6>>)
        .def("Evaluate", &Norms::P::Evaluate<array_1d<double, 9>>)
        .def("Evaluate", &Norms::P::Evaluate<Vector>)
        .def("Evaluate", &Norms::P::Evaluate<Matrix>)
        .def("__str__", &Norms::P::Info)
        ;
    py::class_<Norms::Trace, Norms::Trace::Pointer>(norms, "Trace")
        .def(py::init<>())
        .def("Evaluate", &Norms::Trace::Evaluate)
        .def("__str__", &Norms::Trace::Info)
        ;
    py::class_<Norms::LPQ, Norms::LPQ::Pointer>(norms, "LPQ")
        .def(py::init<const double, const double>(), py::arg("p_coefficient"), py::arg("q_coefficient"))
        .def("Evaluate", &Norms::LPQ::Evaluate)
        .def("__str__", &Norms::LPQ::Info)
        ;

    m.def_submodule("MethodUtilities")
        .def("GetNormMethod", &MethodUtilities::GetNormMethod<int>)
        .def("GetNormMethod", &MethodUtilities::GetNormMethod<double>)
        .def("GetNormMethod", &MethodUtilities::GetNormMethod<array_1d<double, 3>>)
        .def("GetNormMethod", &MethodUtilities::GetNormMethod<Vector>)
        .def("GetNormMethod", &MethodUtilities::GetNormMethod<Matrix>)
        .def("RaiseToPower", &MethodUtilities::RaiseToPower<int>)
        .def("RaiseToPower", &MethodUtilities::RaiseToPower<double>)
        .def("RaiseToPower", &MethodUtilities::RaiseToPower<array_1d<double, 3>>)
        .def("RaiseToPower", &MethodUtilities::RaiseToPower<Vector>)
        .def("RaiseToPower", &MethodUtilities::RaiseToPower<Matrix>)
        ;

}

} // namespace Python.
} // Namespace Kratos
