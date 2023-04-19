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
#include "mpi/spatial_containers/geometrical_objects_bins_mpi.h"

namespace Kratos::Python
{

void AddMPISearchStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using ElementsContainerType = ModelPart::ElementsContainerType;
    using ConditionsContainerType = ModelPart::ConditionsContainerType;

    py::class_<GeometricalObjectsBinsMPI, GeometricalObjectsBinsMPI::Pointer, GeometricalObjectsBins>(m, "GeometricalObjectsBinsMPI")
    .def(py::init<ElementsContainerType&>())
    .def(py::init<ConditionsContainerType&>())
    ;
}

}