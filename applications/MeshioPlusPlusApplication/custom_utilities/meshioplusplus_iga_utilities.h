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

#pragma once

// System includes
#include <string>

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class MeshioPlusPlusIgaUtilities
 * @ingroup MeshioPlusPlusApplication
 * @brief The isogeometric layer of meshio++ (v16.43.0): describe, sample, refine and evaluate
 * NURBS models.
 * @details A NURBS model - the trimmed B-rep of a Kratos CAD JSON, STEP or IGES file, or the
 * patches of an MFEM NURBS mesh - is not a mesh: it has no ModelPart representation (a ModelPart
 * holds linear and Lagrange entities only). The operations of this class therefore take the
 * model from a *file* and hand back either a sampled ModelPart (@ref Sample), another file
 * (@ref Refine) or plain values (@ref Describe, @ref Evaluate), so the exact geometry never
 * round-trips through a lossy mesh. The format is resolved from the file name as
 * @ref MeshioPlusPlusIO does; "cad_json", "step" and "iges" give a B-rep (STEP and IGES need a
 * meshio++ built with OpenCASCADE), "mfem" a NURBS mesh.
 *
 * Refinement is exact (knot insertion, then degree elevation applied before it, so k-refinement
 * keeps the highest continuity): the geometry and its parametrisation do not change, so a trimmed
 * face's trimming curves stay exact. A refined B-rep is written back as Kratos CAD JSON and a
 * refined MFEM model as an MFEM NURBS mesh.
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_MESHIOPLUSPLUS_APPLICATION) MeshioPlusPlusIgaUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MeshioPlusPlusIgaUtilities
    KRATOS_CLASS_POINTER_DEFINITION(MeshioPlusPlusIgaUtilities);

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief The default settings of the operations of this class.
     * @details "format" ("auto" resolves it from the file name), "cad_deflection" and
     * "cad_angular_deflection" (the triangulation of a STEP/IGES file, 0 deflection takes 1e-3 of
     * the bounding-box diagonal), "subdivisions" and "elevate" (@ref Refine: per direction),
     * "sample_subdivisions", "trim_samples", "trim" and "edges" (@ref Sample), "index",
     * "parameters", "derivatives" and "trim_samples" (@ref Evaluate).
     */
    static Parameters GetDefaultParameters();

    /**
     * @brief Summarizes a NURBS model.
     * @param rFileName The model file.
     * @param Settings "format", "cad_deflection", "cad_angular_deflection".
     * @return `{"kind" : "brep"|"nurbs", "number_of_patches", "number_of_breps", "number_of_faces",
     * "number_of_trimmed_faces", "number_of_edges", "number_of_vertices", "patches" : [{"degrees",
     * "number_of_control_points", "rational"}]}`; the B-rep counts are 0 for a NURBS mesh.
     */
    static Parameters Describe(
        const std::string& rFileName,
        Parameters Settings = Parameters(R"({})")
        );

    /**
     * @brief Samples a NURBS model into the destination model part.
     * @details Faces become surface cells (conditions), edges line cells, a volume patch
     * elements, all as Lagrange cells of the highest degree the model has - only the ones Kratos
     * has a geometry for can become entities.
     * @param rFileName The model file.
     * @param Settings "sample_subdivisions" (cells per knot span and direction, 0 automatic),
     * "trim_samples", "trim" (keep only the cells inside the trimming loops) and "edges".
     * @param rDestination The (empty) model part to fill.
     * @return `{"number_of_nodes", "number_of_elements", "number_of_conditions"}`
     */
    static Parameters Sample(
        const std::string& rFileName,
        Parameters Settings,
        ModelPart& rDestination
        );

    /**
     * @brief Refines a NURBS model and writes it to a file.
     * @param rInputFileName The model file.
     * @param rOutputFileName Where to write: Kratos CAD JSON for a B-rep, an MFEM NURBS mesh for a
     * NURBS mesh.
     * @param Settings "subdivisions" (h-refinement: equal parts every knot span is split into, 3
     * integers, 1 for none) and "elevate" (p-refinement: degree elevation, 3 integers, 0 for none).
     * @return `{"kind", "number_of_patches", "patches" : [...]}` of the refined model, as @ref Describe.
     */
    static Parameters Refine(
        const std::string& rInputFileName,
        const std::string& rOutputFileName,
        Parameters Settings
        );

    /**
     * @brief Evaluates one patch (NURBS mesh) or face (B-rep) at parametric coordinates.
     * @param rFileName The model file.
     * @param Settings "index" (the patch or face, 0-based), "parameters" (an N x parametric
     * dimension matrix), "derivatives" (0 or 1) and "trim_samples".
     * @return `{"points" : N x space dimension matrix, and with "derivatives" = 1 also
     * "derivatives" (N x parametric dimension * space dimension, row-major per point) and, for a
     * surface in 3-D, "normals" (N x 3); for a face "inside" (N-vector, 1 inside the trimming
     * loops); for a NURBS mesh "field_<name>" matrices}`
     */
    static Parameters Evaluate(
        const std::string& rFileName,
        Parameters Settings
        );

    ///@}
}; // Class MeshioPlusPlusIgaUtilities

///@}

} // namespace Kratos
