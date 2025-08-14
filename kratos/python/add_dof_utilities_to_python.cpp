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

    py::class_<DofArrayUtilities>(m, "DofArrayUtilities")
        .def_static("SetUpDofArray", [](const ModelPart& rModelPart, DofArrayUtilities::DofsArrayType& rDofArray){DofArrayUtilities::SetUpDofArray(rModelPart, rDofArray);})
        .def_static("SetUpDofArray", [](const ModelPart& rModelPart, DofArrayUtilities::DofsArrayType& rDofArray, const unsigned int EchoLevel){DofArrayUtilities::SetUpDofArray(rModelPart, rDofArray, EchoLevel);})
    ;
}

}  // namespace Kratos::Python.
