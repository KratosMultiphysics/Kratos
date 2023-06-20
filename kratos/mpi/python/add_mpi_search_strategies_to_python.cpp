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

template <SpatialContainer TSearchBackend>
void DefineSpecializedSpatialSearchMPI(pybind11::module& m, const std::string& rClassName)
{
    using SpatialSearchType = SpecializedSpatialSearchMPI<TSearchBackend>;
    using SpatialSearchPointerType = typename SpecializedSpatialSearchMPI<TSearchBackend>::Pointer;
    using BaseSpatialSearchType = SpecializedSpatialSearch<TSearchBackend>;

    pybind11::class_<SpatialSearchType, SpatialSearchPointerType, BaseSpatialSearchType>(m, rClassName.c_str())
    .def(pybind11::init<>())
    .def(pybind11::init<Parameters>())
    ;
}

void AddMPISearchStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // The factory of the MPI search strategies
    py::class_<SpecializedSpatialSearchMPIFactory, SpecializedSpatialSearchMPIFactory::Pointer, SpatialSearch>(m, "SpecializedSpatialSearchMPI")
    .def(py::init< >())
    .def(py::init<Parameters>())
    ;

    // Register the specializations
    DefineSpecializedSpatialSearchMPI<SpatialContainer::KDTree>(m, "SpatialSearchKDTreeMPI");
    DefineSpecializedSpatialSearchMPI<SpatialContainer::Octree>(m, "SpatialSearchOctreeMPI");
    DefineSpecializedSpatialSearchMPI<SpatialContainer::BinsStatic>(m, "SpatialSearchBinsStaticMPI");
    DefineSpecializedSpatialSearchMPI<SpatialContainer::BinsDynamic>(m, "SpatialSearchBinsDynamicMPI");
}

}  // namespace Kratos::Python