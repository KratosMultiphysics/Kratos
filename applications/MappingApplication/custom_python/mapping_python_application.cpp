//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

#if defined(KRATOS_PYTHON)
// External includes

// Project includes
#include "includes/define_python.h"
#include "mapping_application.h"
#include "custom_python/add_custom_mappers_to_python.h"


namespace Kratos {
namespace Python {


PYBIND11_MODULE(KratosMappingApplication, m)
{
    namespace py = pybind11;

    py::class_<KratosMappingApplication,
        KratosMappingApplication::Pointer,
        KratosApplication >(m,"KratosMappingApplication")
        .def(py::init<>())
        ;

    AddCustomMappersToPython(m);
}

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
