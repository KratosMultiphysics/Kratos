//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░                        Application
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
        "smooth, reorder, partition, crop, slice, isosurface, quality, stats, data_calc, "
        "data_condition, data_manage, data_info, point_data_to_cell_data, "
        "cell_data_to_point_data, ...) exposed as Kratos utilities. Every operation is "
        "driven by Parameters keyed by an \"operation\" name mirroring the meshio++ command "
        "line verbs. Field data (nodal/elemental/conditional variables and flags) can be "
        "carried through with the same settings MeshioPlusPlusIO uses; see GetDefaultParameters.")
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
        .def_static("Interpolate", &MeshioPlusPlusMeshOperations::Interpolate,
            py::arg("source_model_part"), py::arg("target_model_part"), py::arg("settings"),
            py::arg("destination_model_part"),
            "Samples source_model_part's field data onto target_model_part's geometry, "
            "filling destination_model_part with the target's topology plus the interpolated "
            "arrays. Not reachable through Execute: unlike every other operation this needs "
            "two independent meshes.")
        .def_static("ConservativeInterpolate", &MeshioPlusPlusMeshOperations::ConservativeInterpolate,
            py::arg("source_model_part"), py::arg("target_model_part"), py::arg("settings"),
            py::arg("destination_model_part"),
            "The conservative sibling of Interpolate: redistributes by intersected cell measure "
            "so the summed quantity is preserved, rather than sampling point values. The right "
            "choice for an extensive field (a mass, a heat load), the wrong one for an intensive "
            "one (a temperature). Also needs two independent meshes, so also not reachable "
            "through Execute.")
        .def_static("Grid", &MeshioPlusPlusMeshOperations::Grid,
            py::arg("settings"), py::arg("destination_model_part"),
            "Builds a regular hexahedron lattice from \"dims\"/\"origin\"/\"spacing\". Not "
            "reachable through Execute: it takes no source model part, being the only "
            "generator rather than a transformation.")
        .def_static("DistanceToSurface", &MeshioPlusPlusMeshOperations::DistanceToSurface,
            py::arg("query_model_part"), py::arg("surface_model_part"), py::arg("settings"),
            py::arg("destination_model_part"),
            "Attaches the signed distance from query_model_part to surface_model_part. Set "
            "\"output\" to a registered Variable name (DISTANCE, say) to get the result back: "
            "meshio++ names it \"sdf:distance\", which no Kratos Variable can be.")
        .def_static("CheckSurfaceWatertight", &MeshioPlusPlusMeshOperations::CheckSurfaceWatertight,
            py::arg("surface_model_part"),
            "Reports whether a surface is watertight, with the boundary-edge, non-manifold, "
            "inconsistent-winding and degenerate-triangle counts behind the verdict.")
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
