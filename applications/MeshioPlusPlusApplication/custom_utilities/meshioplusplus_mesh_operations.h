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
#include <vector>

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class MeshioPlusPlusMeshOperations
 * @ingroup MeshioPlusPlusApplication
 * @brief Exposes the meshio++ mesh and data operations to Kratos.
 * @details Every operation is data-driven through @ref Parameters, keyed by an "operation"
 * name that mirrors the meshio++ command line verbs, so adding an upstream operation is a
 * table entry rather than a new class. The conversion to and from meshio++ is shared with
 * @ref MeshioPlusPlusIO (see meshioplusplus_conversion_utilities.h), so no operation
 * re-implements the staging.
 *
 * Operations fall into two groups:
 *  - mesh-producing (clean, transform, refine, ...): the result is written into the
 *    destination model part, and the returned @ref Parameters carry the operation's own
 *    report (counts, tolerances actually applied, ...);
 *  - report-only (stats, quality, diff, ...): the destination is left untouched and the
 *    whole result is in the returned @ref Parameters.
 *
 * @note meshio++ is serial: these operations do not support distributed model parts. The
 * intended distributed workflow is "partition" with ghost layers feeding an MPI assembly.
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_MESHIOPLUSPLUS_APPLICATION) MeshioPlusPlusMeshOperations
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MeshioPlusPlusMeshOperations
    KRATOS_CLASS_POINTER_DEFINITION(MeshioPlusPlusMeshOperations);

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief The names of every operation @ref Execute accepts.
     * @return The supported operation names, sorted.
     */
    static std::vector<std::string> GetSupportedOperations();

    /**
     * @brief The default settings shared by every operation.
     * @return The default parameters.
     */
    static Parameters GetDefaultParameters();

    /**
     * @brief Applies a meshio++ operation to a model part.
     * @details Mesh-producing operations fill @p rDestination. The multi-output ones
     * ("split", "partition") instead create one model part per piece, named
     * `<destination>_<operation>_<index>` and registered as siblings in the same @ref Model;
     * their names are listed in the returned report. They cannot be sub model parts, since a
     * Kratos sub model part is a view into its root's entity containers while each piece is
     * independently renumbered from 1. Report-only operations leave @p rDestination untouched.
     * @param rSource The model part to operate on.
     * @param Settings The operation name plus its own settings (see @ref GetDefaultParameters).
     * @param rDestination The model part receiving the result (expected empty).
     * @return The operation's report.
     */
    static Parameters Execute(
        const ModelPart& rSource,
        Parameters Settings,
        ModelPart& rDestination
        );

    /**
     * @brief Merges several model parts into one.
     * @param rSources The model parts to merge, in order.
     * @param Settings Merge settings ("weld", "tolerance", "source_tag", "data_policy",
     *                 "drop_duplicate_cells").
     * @param rDestination The model part receiving the merged mesh (expected empty).
     * @return The merge report.
     */
    static Parameters Merge(
        const std::vector<const ModelPart*>& rSources,
        Parameters Settings,
        ModelPart& rDestination
        );

    /**
     * @brief Compares two model parts with tolerances.
     * @param rFirst The first model part.
     * @param rSecond The second model part.
     * @param Settings Diff settings ("absolute_tolerance", "relative_tolerance", "unordered").
     * @return The diff report, including the per-section verdicts.
     */
    static Parameters Diff(
        const ModelPart& rFirst,
        const ModelPart& rSecond,
        Parameters Settings = Parameters(R"({})")
        );

    /**
     * @brief Aggregate mesh statistics: bounding box, centroid, per-cell-type counts,
     * total area, signed and unsigned volume, inverted-cell count.
     * @param rSource The model part to measure.
     * @return The statistics report.
     */
    static Parameters ComputeStatistics(const ModelPart& rSource);

    /**
     * @brief Per-cell quality metrics summarized per metric (scaled Jacobian, aspect ratio,
     * skewness, ...) plus the inverted and degenerate cell counts.
     * @param rSource The model part to measure.
     * @return The quality report.
     */
    static Parameters ComputeQuality(const ModelPart& rSource);

    /**
     * @brief The bandwidth of the model part's node adjacency graph.
     * @details Useful to quantify what the "reorder" operation achieves.
     * @param rSource The model part to measure.
     * @return The bandwidth.
     */
    static std::size_t ComputeBandwidth(const ModelPart& rSource);

    ///@}
}; // Class MeshioPlusPlusMeshOperations

///@}

} // namespace Kratos
