//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

// --- External Includes ---
#include <pybind11/stl.h>
#include <pybind11/operators.h>

// --- Core Includes ---
#include "add_adjoint_strategies_to_python.hpp"
#include "adjoint/adjoint_scheme.hpp"
#include "adjoint/static_adjoint_scheme.hpp"
#include "spaces/ublas_space.h"


namespace Kratos::Python {


template <class TS, class TD>
void AddAdjointSchemesToPython(pybind11::module_& rModule) {
    pybind11::class_<
        AdjointScheme<TS,TD>,
        Scheme<TS,TD>,
        typename AdjointScheme<TS,TD>::Pointer
    >(rModule, "AdjointScheme")
        .def(pybind11::init<>())
        .def(pybind11::init<Parameters,ResponseFunction::Pointer>())
        .def(
            "SetResponseFunction",
            &AdjointScheme<TS,TD>::SetResponseFunction,
            pybind11::arg("pResponseFunction"))
        .def(
            "GetResponseFunction",
            &AdjointScheme<TS,TD>::GetResponseFunction)
        ;

    pybind11::class_<
        StaticAdjointScheme<TS,TD>,
        AdjointScheme<TS,TD>,
        typename StaticAdjointScheme<TS,TD>::Pointer
    >(rModule, "StaticAdjointScheme")
        .def(pybind11::init<>())
        .def(pybind11::init<Parameters,ResponseFunction::Pointer>())
        ;
} // void AddAdjointSchemesToPython


void AddAdjointStrategiesToPython(pybind11::module_& rModule) {
    AddAdjointSchemesToPython<
        TUblasSparseSpace<double>,
        TUblasDenseSpace<double>
    >(rModule);
}


} // namespace Kratos::Python
