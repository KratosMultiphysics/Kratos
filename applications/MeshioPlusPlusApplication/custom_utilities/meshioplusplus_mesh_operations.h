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
 *  - mesh-producing (clean, transform, refine, remesh, ...): the result is written into the
 *    destination model part, and the returned @ref Parameters carry the operation's own
 *    report (counts, tolerances actually applied, ...);
 *  - report-only (stats, quality, diff, data_info, data_integrate): the destination is left
 *    untouched and the whole result is in the returned @ref Parameters.
 *
 * Four operations are *not* reachable through @ref Execute because they do not have its
 * one-mesh-in/one-mesh-out shape: @ref Interpolate and @ref ConservativeInterpolate and
 * @ref UndoGreen each need two independent meshes, and @ref Grid needs none at all.
 *
 * @ref Execute optionally carries field data through the operation, using the same
 * "nodal_solution_step_data_variables" / "nodal_data_value_variables" / "nodal_flags" /
 * "element_data_value_variables" / "element_flags" / "condition_data_value_variables" /
 * "condition_flags" / "gauss_point_variables_in_elements" / "write_ids" settings
 * @ref MeshioPlusPlusIO uses (all empty/false by default, so an operation run without them
 * behaves exactly as before this was added). A resulting array is written back onto the
 * destination as non-historical data only when its name matches a registered @ref Variable
 * whose component count agrees - an operation's own invented array names never carry
 * through (see @ref Internals::MeshToModelPart), which is the one thing to know before
 * relying on "data_calc"'s "output" setting or similar: point it at an existing variable
 * name to get the result back into Kratos.
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
     * @brief Samples one model part's field data onto another's geometry.
     * @details Not reachable through @ref Execute: unlike every other operation this needs
     * two independent meshes rather than one. The shared field data settings (see the
     * class-level details) apply to *both* @p rSource and @p rTarget, using the same
     * variable/flag names for each - the usual case is that @p rSource carries the named
     * data and @p rTarget does not, but a name present on both is exactly what
     * "on_conflict" resolves. The target's own topology is what is written into
     * @p rDestination, with the interpolated arrays added.
     * @param rSource The model part carrying the field data to sample.
     * @param rTarget The model part whose geometry the data is sampled onto.
     * @param Settings Interpolation settings ("method": "nearest"/"barycentric"; "names":
     *                 array names to interpolate, empty = every array @p rSource's field
     *                 data settings collect; "extrapolate"; "default_value"; "on_conflict":
     *                 "error"/"overwrite"/"suffix") plus the shared field data settings.
     * @param rDestination The model part receiving the target's geometry plus the
     *                     interpolated data (expected empty).
     * @return The interpolation report.
     */
    static Parameters Interpolate(
        const ModelPart& rSource,
        const ModelPart& rTarget,
        Parameters Settings,
        ModelPart& rDestination
        );

    /**
     * @brief Samples one model part's field data onto another's geometry, conserving the
     * integral rather than the point values.
     * @details The conservative sibling of @ref Interpolate, and not reachable through
     * @ref Execute for the same reason: it needs two independent meshes. Where @ref Interpolate
     * samples (nearest or barycentric), this redistributes by intersected cell measure, so the
     * summed quantity is preserved across the transfer - the right choice for an extensive
     * field (a mass, a heat load) and the wrong one for an intensive one (a temperature).
     *
     * The shared field data settings apply to *both* @p rSource and @p rTarget, exactly as in
     * @ref Interpolate.
     * @param rSource The model part carrying the field data to transfer.
     * @param rTarget The model part whose geometry the data is transferred onto.
     * @param Settings "names" (array names to transfer, empty = every array @p rSource's field
     *                 data settings collect), "default_value", "on_conflict"
     *                 ("error"/"overwrite"/"suffix") plus the shared field data settings.
     * @param rDestination The model part receiving the target's geometry plus the transferred
     *                     data (expected empty).
     * @return The transfer report.
     */
    static Parameters ConservativeInterpolate(
        const ModelPart& rSource,
        const ModelPart& rTarget,
        Parameters Settings,
        ModelPart& rDestination
        );

    /**
     * @brief Removes the green closure cells a previous @ref Execute "refine" added.
     * @details Not reachable through @ref Execute: it needs both the coarse mesh and the
     * refined one. A selective `refine(closure="redgreen")` adds green cells purely to keep the
     * mesh conforming; they carry no refinement information of their own and degrade element
     * quality, so a solver that has finished with a level undoes them before refining again.
     * Takes no settings - the coarse mesh is what identifies the green groups.
     * @param rCoarse The mesh as it was before the refinement.
     * @param rFine The refined mesh whose green closure is to be undone.
     * @param rDestination The model part receiving the result (expected empty).
     * @return The report: the undone-group and removed-cell counts plus the entity counts.
     */
    static Parameters UndoGreen(
        const ModelPart& rCoarse,
        const ModelPart& rFine,
        ModelPart& rDestination
        );

    /**
     * @brief Builds a regular hexahedron lattice from nothing.
     * @details Not reachable through @ref Execute: it is meshio++'s only *generator*, taking
     * no source model part at all where every other operation transforms one.
     * @param Settings Lattice settings ("dims": cell counts per axis; "origin": the lattice's
     *                 low corner; "spacing": the cell size per axis; "max_cells").
     * @param rDestination The model part receiving the grid (expected empty).
     * @return The lattice report ("dims", "origin", "spacing" and the entity counts).
     */
    static Parameters Grid(
        Parameters Settings,
        ModelPart& rDestination
        );

    /**
     * @brief Attaches the signed distance from a query model part to a surface.
     * @details Not reachable through @ref Execute: like @ref Interpolate this needs two
     * independent meshes. @p rQuery's geometry is copied unchanged and annotated.
     *
     * meshio++ names the result `"sdf:distance"`, which no Kratos `Variable` is (they are
     * never colon-namespaced), so it could never be read back off @p rDestination. Set
     * "output" to a registered variable name - `DISTANCE`, for a level-set initialization -
     * and the array is renamed before the write-back, which is what makes the result usable
     * rather than merely computed.
     * @param rQuery The model part to annotate.
     * @param rSurface The surface to measure against.
     * @param Settings Surface-distance settings ("sdf_sign", "sdf_weight", "sdf_location",
     *                 "band", "record_closest_cell", "record_inside", "watertight_check",
     *                 "grid_cell_size", "max_winding_work") plus "output".
     * @param rDestination The model part receiving the annotated query mesh (expected empty).
     * @return The report: the surface verdict and the banded-query count.
     */
    static Parameters DistanceToSurface(
        const ModelPart& rQuery,
        const ModelPart& rSurface,
        Parameters Settings,
        ModelPart& rDestination
        );

    /**
     * @brief Reports what is wrong with a surface, in numbers rather than a bare flag.
     * @details Report-only, so it takes no destination - the @ref ComputeQuality shape.
     * A signed distance needs a closed, consistently wound surface; this is how to find out
     * whether one is, before trusting a sign.
     * @param rSurface The surface to check.
     * @return The verdict: "watertight" plus the boundary-edge, non-manifold-edge,
     *         inconsistent-pair and degenerate-triangle counts.
     */
    static Parameters CheckSurfaceWatertight(const ModelPart& rSurface);

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
