//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/model_part.h"
#include "mpi/python/add_mpi_search_strategies_to_python.h"
#include "mpi/spatial_containers/specialized_spatial_search_mpi.h"
#include "mpi/spatial_containers/specialized_spatial_search_mpi_factory.h"

namespace Kratos::Python
{

void AddMPISearchStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<SpecializedSpatialSearchMPIFactory, SpecializedSpatialSearchMPIFactory::Pointer, SpatialSearch>(m, "SpecializedSpatialSearch")
    .def(py::init< >())
    .def(py::init<Parameters>())
    ;

    py::class_<SpecializedSpatialSearchMPI<SpatialContainer::KDTree>, SpecializedSpatialSearchMPI<SpatialContainer::KDTree>::Pointer, SpecializedSpatialSearch<SpatialContainer::KDTree>>(m, "SpatialSearchKDTree")
    .def(py::init< >())
    .def(py::init<Parameters>())
    ;

    py::class_<SpecializedSpatialSearchMPI<SpatialContainer::Octree>, SpecializedSpatialSearchMPI<SpatialContainer::Octree>::Pointer, SpecializedSpatialSearch<SpatialContainer::Octree>>(m, "SpatialSearchOctree")
    .def(py::init< >())
    .def(py::init<Parameters>())
    ;

    py::class_<SpecializedSpatialSearchMPI<SpatialContainer::BinsStatic>, SpecializedSpatialSearchMPI<SpatialContainer::BinsStatic>::Pointer, SpecializedSpatialSearch<SpatialContainer::BinsStatic>>(m, "SpatialSearchBinsStatic")
    .def(py::init< >())
    .def(py::init<Parameters>())
    ;

    py::class_<SpecializedSpatialSearchMPI<SpatialContainer::BinsDynamic>, SpecializedSpatialSearchMPI<SpatialContainer::BinsDynamic>::Pointer, SpecializedSpatialSearch<SpatialContainer::BinsDynamic>>(m, "SpatialSearchBinsDynamic")
    .def(py::init< >())
    .def(py::init<Parameters>())
    ;

}

}  // namespace Kratos::Python