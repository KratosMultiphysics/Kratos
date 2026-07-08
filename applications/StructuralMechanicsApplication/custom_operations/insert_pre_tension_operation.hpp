// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// --- Core Includes ---
#include "operations/operation.h" // Operation
#include "containers/model.h" // Model
#include "modeler/modeler.h" // Modeler

// --- STL Includes ---
#include <unordered_set>
#include <unordered_map>


namespace Kratos {


/// @defgroup pre_tensioning Pre-Tensioning
/// @brief Apply pre-tensioning to 1, 2, or 3D structural parts.
/// @details Pre-tensioning approximately models connector parts (e.g.: bolts, screws, etc.) subject to
///          initial loading within an assembly. It is meant to be a simplified alternative to high-fidelity
///          analyses of such assemblies, with less focus around the connector part but similar effects on
///          the rest of the structure.
///
/// @subsection pre_tensioning_surface Pre-Tensioning Surface
///          The pre-tensioning system acts on a surface defined by a set of @ref Geometry "geometries"
///          within a @ref ModelPart "model part". The following requirements apply to this set of geometries:
///          - they @b must be defined on element boundaries (faces of polyhedra, edges of polygons, vertices of
///            lines),
///          - they @b must partition the connector part's elements into 2 distinct sets,
///          - they @b should lie on the same plane as much as possible, and
///          - they @b must be stored in a separate sub model part of the structure that only contains
///            these geometries.
///
///          A <em>pre-tensioning plane</em> is computed by averaging the <em>pre-tensioning surface</em>.
///          The normal of this plane will be the basis of displacement constraints introduced later in this
///          process.
///
/// @subsection pre_tensioning_node_duplication Node Duplication
///          The <em>pre-tensioning surface</em> partitions adjacent elements into two distinct groups:
///          the @a positive side and the @a negative side. Nodes lying on the surface are duplicated, and
///          negative side elements' nodes are replaced with them, while elements on the positive side
///          are left unchanged. This effectively means that the connector part gets cut in half along
///          the <em>pre-tensioning surface</em>. Duplicated nodes are inserted into the sub model part
///          containing the geometries that define the surface.
///
/// @subsection pre_tensioning_in_plane_constraints In-Plane Constraints
///          After duplicating nodes, some @ref MultifreedomConstraint "constraints" must be inserted
///          that restrict displacements (and rotations). First of all, if the nodes have DoFs for rotations,
///          the rotations of the duplicates must match their original counterparts.
///          @f[
///             \phi_o^i - \phi_d^i = 0 \quad \forall \quad 1 \leq i \leq m
///          @f]
///          where
///          - @f$\phi_o^i@f$ is the rotation of the original node @f$i@f$,
///          - @f$\phi_d^i@f$ is the rotation of node @f$i@f$'s duplicate, and
///          - @f$m@f$ is the number of original nodes on the <em>pre-tensioning surface</em>.
///
///          Also, the relative
///          in-plane displacement component of original-duplicate node pairs must vanish.
///          @f[
///             u_o^i - u_d^i - \left< u_o^i - u_d^i, n \right> n = 0 \quad \forall \quad 1 \leq i \leq m
///          @f]
///          where
///          - @f$u_o^i@f$ is the displacement of the original node @f$i@f$,
///          - @f$u_d^i@f$ is the displacement of node @f$i@f$'s duplicate,
///          - @f$n@f$ is the unit normal of the <em>pre-tensioning plane</em>, and
///          - @f$\left< \cdot , \cdot \right>@f$ denotes an inner product.
///
/// @subsection pre_tensioning_out_of_plane_constraints Out-of-Plane Constraints
///          To apply a single pre-tensioning value (load or prescribed displacement) to the surface,
///          the out-of-plane displacement components on each set of nodes (original or duplicated) are
///          averaged and tied to a newly inserted virtual DoF.
///          @f[
///                 u_v - \frac{1}{m} \sum_{i=1}^m{ \left< u_o^i, n \right> - \left< u_d^i, n \right> } = 0
///          @f]
///          where
///          - @f$u_v^o@f$ is the newly inserted virtual DoF.
///
///          Once the virtual degree-of-freedom @f$u_v@f$ is defined, it can either
///          be loaded (Neumann-type) or fixed (Dirichlet-type).
///
/// @subsection pre_tensioning_implementation_notes Implementation Notes
///          Since in-plane displacement constraints reuse DoFs and are dense, they cannot be imposed
///          via master-slave elimination unless DoFs are rotated to line up with the
///          <em>pre-tensioning plane</em>. Since there is currently no robust way of doing this in
///          Kratos, constraints must be imposed via the method of augmented lagrange multipliers,
///          which is implemented in @ref p_multigrid "PMultigridBuilderAndSolver".
///          @note Any analysis involving pre-tensioning must use @ref PMultigridBuilderAndSolver.



/// @brief Base class for @ref InsertDirichletPreTensionOperation and @ref InsertNeumannPreTensionOperation.
/// @details Cut the mesh along the provided surface and apply fix the average
///          out-of-plane displacement to a newly inserted DoF. This DoF can later
///          be loaded or constrained by applying subsequent processes to the analysis.
/// @see @ref pre_tensioning "Pre-Tensioning"
/// @ingroup pre_tensioning
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) InsertPreTensionOperation final : public Operation {
public:
    KRATOS_CLASS_POINTER_DEFINITION(InsertPreTensionOperation);

    InsertPreTensionOperation() noexcept = default;

    InsertPreTensionOperation(Model& rModel, Parameters Settings);

    InsertPreTensionOperation(InsertPreTensionOperation&&) noexcept;

    ~InsertPreTensionOperation();

    void Execute() override;

    //void Apply(double Magnitude);

    const Parameters GetDefaultParameters() const override;

    virtual std::string Info() const;

protected:
    /// @details Inserts a constraint that ties the average out-of-plane
    ///          relative displacement of duplicated nodes to a single DoF.
    void InsertControlNodeConstraints(
        ModelPart& rModelPart,
        array_1d<double,3> SurfaceNormal,
        const std::unordered_map<Node*,Node::Pointer> rDuplicateNodeMap,
        Node::Pointer pControlNode,
        const std::unordered_set<const Dof<double>*> rPositiveSideDofs);

    void AddControlDoFs(Node& rNode);

    struct Impl;
    std::unique_ptr<Impl> mpImpl;

private:
    InsertPreTensionOperation(const InsertPreTensionOperation&) = delete;

    InsertPreTensionOperation& operator=(const InsertPreTensionOperation&) = delete;

    Node::Pointer mpControlNode;

    std::size_t mPreTensionSurfaceSize;
}; // class InsertPreTensionOperation


/// @see @ref pre_tensioning "Pre-Tensioning"
/// @ingroup pre_tensioning
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) PreTensionModeler final : public Modeler {
public:
    KRATOS_CLASS_POINTER_DEFINITION(PreTensionModeler);

    PreTensionModeler() noexcept = default;

    PreTensionModeler(
        Model& rModel,
        Parameters Settings);

    Modeler::Pointer Create(
        Model& rModel,
        Parameters Settings) const override;

    void SetupModelPart() override;

    const Parameters GetDefaultParameters() const override;

private:
    Model* mpModel;
}; // class PreTensionModeler


} // namespace Kratos
