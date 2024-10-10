//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes

// External includes

// Project includes
#include "python/add_globals_to_python.h"
#include "includes/define_python.h"
#include "includes/global_variables.h"

namespace Kratos::Python {

void AddGlobalsToPython(pybind11::module& m)
{
    namespace py = pybind11;
    using DataLocation = Kratos::Globals::DataLocation;

    auto m_globals = m.def_submodule("Globals");

    py::enum_<DataLocation>(m_globals,"DataLocation")
        .value("NodeHistorical",    DataLocation::NodeHistorical)
        .value("NodeNonHistorical", DataLocation::NodeNonHistorical)
        .value("Element",           DataLocation::Element)
        .value("Condition",         DataLocation::Condition)
        .value("ProcessInfo",       DataLocation::ProcessInfo)
        .value("ModelPart",         DataLocation::ModelPart)
        ;
}

} // namespace Kratos::Python
