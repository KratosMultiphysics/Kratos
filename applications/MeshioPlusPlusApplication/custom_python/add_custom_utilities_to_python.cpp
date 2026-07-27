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
#include <vector>

// External includes
#include <pybind11/stl.h>

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/meshioplusplus_mesh_operations.h"
#include "includes/model_part.h"

namespace Kratos::Python
{

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<MeshioPlusPlusMeshOperations, MeshioPlusPlusMeshOperations::Pointer>(
        m, "MeshioPlusPlusMeshOperations",
        "The meshio++ mesh and data operations (clean, transform, split, refine, decimate, "
        "smooth, reorder, partition, crop, slice, isosurface, quality, stats, ...) exposed as "
        "Kratos utilities. Every operation is driven by Parameters keyed by an \"operation\" "
        "name mirroring the meshio++ command line verbs.")
        .def_static("GetSupportedOperations", &MeshioPlusPlusMeshOperations::GetSupportedOperations,
            "The names accepted by the \"operation\" setting of Execute.")
        .def_static("GetDefaultParameters", &MeshioPlusPlusMeshOperations::GetDefaultParameters,
            "The default settings shared by every operation.")
        .def_static("Execute", &MeshioPlusPlusMeshOperations::Execute,
            py::arg("source_model_part"), py::arg("settings"), py::arg("destination_model_part"),
            "Applies an operation, filling the destination model part (or, for \"split\" and "
            "\"partition\", one sub model part per piece). Returns the operation's report.")
        .def_static("Merge",
            [](const py::list& rSources, Parameters Settings, ModelPart& rDestination) {
                std::vector<const ModelPart*> sources;
                sources.reserve(rSources.size());
                for (const auto& r_item : rSources) {
                    sources.push_back(r_item.cast<const ModelPart*>());
                }
                return MeshioPlusPlusMeshOperations::Merge(sources, Settings, rDestination);
            },
            py::arg("source_model_parts"), py::arg("settings"), py::arg("destination_model_part"),
            "Merges several model parts into one, optionally welding coincident nodes.")
        .def_static("Diff", &MeshioPlusPlusMeshOperations::Diff,
            py::arg("first_model_part"), py::arg("second_model_part"),
            py::arg("settings") = Parameters(R"({})"),
            "Compares two model parts with absolute and relative tolerances.")
        .def_static("ComputeStatistics", &MeshioPlusPlusMeshOperations::ComputeStatistics,
            py::arg("model_part"),
            "Bounding box, per-cell-type counts, area, signed/unsigned volume, inverted cells.")
        .def_static("ComputeQuality", &MeshioPlusPlusMeshOperations::ComputeQuality,
            py::arg("model_part"),
            "Per-cell quality metrics summarized per metric, plus inverted/degenerate counts.")
        .def_static("ComputeBandwidth", &MeshioPlusPlusMeshOperations::ComputeBandwidth,
            py::arg("model_part"),
            "The bandwidth of the node adjacency graph (what \"reorder\" reduces).")
        ;
}

} // namespace Kratos::Python
