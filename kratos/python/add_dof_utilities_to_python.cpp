//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/model_part.h"
#include "python/add_dof_utilities_to_python.h"
#include "utilities/dof_utilities/dof_array_utilities.h"
#include "utilities/dof_utilities/block_build_dof_array_utility.h"
#include "utilities/dof_utilities/elimination_build_dof_array_utility.h"

PYBIND11_MAKE_OPAQUE(Kratos::DofArrayUtilities::SlaveToMasterDofsMap); // prevents pybind11 from trying to convert the unordered_map into a dict

namespace Kratos::Python
{

void AddDofUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<BlockBuildDofArrayUtility>(m, "BlockBuildDofArrayUtility")
        .def_static("SetUpDofArray", [](const ModelPart& rModelPart, BlockBuildDofArrayUtility::DofsArrayType& rDofArray){BlockBuildDofArrayUtility::SetUpDofArray(rModelPart, rDofArray);})
        .def_static("SetUpDofArray", [](const ModelPart& rModelPart, BlockBuildDofArrayUtility::DofsArrayType& rDofArray, const unsigned int EchoLevel){BlockBuildDofArrayUtility::SetUpDofArray(rModelPart, rDofArray, EchoLevel);})
        .def_static("SetUpDofArray", [](const ModelPart& rModelPart, BlockBuildDofArrayUtility::DofsArrayType& rDofArray, const unsigned int EchoLevel, const bool CheckReactionDofs){BlockBuildDofArrayUtility::SetUpDofArray(rModelPart, rDofArray, EchoLevel, CheckReactionDofs);})
    ;

    py::class_<EliminationBuildDofArrayUtility>(m, "EliminationBuildDofArrayUtility")
        .def_static("SetUpDofArray", [](const ModelPart& rModelPart, EliminationBuildDofArrayUtility::DofsArrayType& rDofArray){EliminationBuildDofArrayUtility::SetUpDofArray(rModelPart, rDofArray);})
        .def_static("SetUpDofArray", [](const ModelPart& rModelPart, EliminationBuildDofArrayUtility::DofsArrayType& rDofArray, const unsigned int EchoLevel){EliminationBuildDofArrayUtility::SetUpDofArray(rModelPart, rDofArray, EchoLevel);})
        .def_static("SetUpDofArray", [](const ModelPart& rModelPart, EliminationBuildDofArrayUtility::DofsArrayType& rDofArray, const unsigned int EchoLevel, const bool CheckReactionDofs){EliminationBuildDofArrayUtility::SetUpDofArray(rModelPart, rDofArray, EchoLevel, CheckReactionDofs);})
    ;

    py::class_<DofArrayUtilities::SlaveToMasterDofsMap>(m, "SlaveToMasterDofsMap")
        .def(py::init<>())
        .def("__len__", [](const DofArrayUtilities::SlaveToMasterDofsMap &rSelf)
            { return rSelf.size(); })
        .def("__contains__", [](const DofArrayUtilities::SlaveToMasterDofsMap &rSelf, Node::DofType::Pointer Key)
            { return rSelf.find(Key) != rSelf.end(); })
        .def("__getitem__", [](const DofArrayUtilities::SlaveToMasterDofsMap &rSelf, Node::DofType::Pointer Key)
            { return rSelf.at(Key); }) // returns by value (copy)
        .def("__iter__", [](const DofArrayUtilities::SlaveToMasterDofsMap &rSelf)
            { return py::make_key_iterator(rSelf.begin(), rSelf.end()); }, py::keep_alive<0, 1>()) // iteration over keys
        .def("items", [](const DofArrayUtilities::SlaveToMasterDofsMap &rSelf)
            { return py::make_iterator(rSelf.begin(), rSelf.end()); }, py::keep_alive<0, 1>()); // items() for key/value iteration

    py::class_<DofArrayUtilities>(m, "DofArrayUtilities")
        .def_static("SetUpDofArray", &DofArrayUtilities::SetUpDofArray, py::arg("model_part"), py::arg("dof_array"), py::arg("echo_level") = 0)
        .def_static("SetUpEffectiveDofArray", &DofArrayUtilities::SetUpEffectiveDofArray, py::arg("model_part"), py::arg("dof_array"), py::arg("effective_dof_array"), py::arg("slave_to_master_dofs_map"), py::arg("echo_level") = 0)
        .def_static("SetDofEquationIds", &DofArrayUtilities::SetDofEquationIds, py::arg("dof_array"))
        .def_static("SetEffectiveDofEquationIds", &DofArrayUtilities::SetEffectiveDofEquationIds, py::arg("dof_array"), py::arg("effective_dof_array"))
    ;
}

}  // namespace Kratos::Python.
