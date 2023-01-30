//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi, Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_mappers_to_python.h"
#include "custom_utilities/shared_points_mapper.h"

namespace Kratos::Python
{
void AddMappersToPython(pybind11::module &m)
{
    namespace py = pybind11;

    py::class_<SharedPointsMapper>(m,"SharedPointsMapper")
        .def(py::init<const ModelPart::NodesContainerType &, const ModelPart::NodesContainerType &, double>())
        .def("ScalarMap", &SharedPointsMapper::ScalarMap)
        .def("InverseScalarMap", &SharedPointsMapper::InverseScalarMap)
        .def("VectorMap", &SharedPointsMapper::VectorMap)
        .def("InverseVectorMap", &SharedPointsMapper::InverseVectorMap);
}

} // Namespace Kratos::Python.
