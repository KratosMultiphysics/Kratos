
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "mappers/mapper_define.h"
#include "custom_mappers/beam_mapper.h"

// Include base h
#include "add_custom_mappers_to_python.h"


namespace Kratos {
namespace Python {

void  AddCustomMappersToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using BeamMapperType = BeamMapper<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>;
    const std::string beam_mapper_name("BeamMapper");
    py::class_< BeamMapperType, typename BeamMapperType::Pointer, Mapper<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>>(m, beam_mapper_name.c_str())
        .def("Map", py::overload_cast<const Variable<array_1d<double, 3>>&, const Variable<array_1d<double, 3>>&, const Variable<array_1d<double, 3>>&, Flags>(&BeamMapperType::Map), py::arg("origin_displacement_variable"), py::arg("origin_rotation_variable"), py::arg("destination_displacement_variable"), py::arg("mapping_options") = Kratos::Flags())
        .def("InverseMap", py::overload_cast<const Variable<array_1d<double, 3>>&, const Variable<array_1d<double, 3>>&, const Variable<array_1d<double, 3>>&, Flags>(&BeamMapperType::InverseMap), py::arg("origin_force_variable"), py::arg("origin_moment_variable"), py::arg("destination_force_variable"), py::arg("mapping_options") = Kratos::Flags())
        ;
}

}  // namespace Python.
} // Namespace Kratos
