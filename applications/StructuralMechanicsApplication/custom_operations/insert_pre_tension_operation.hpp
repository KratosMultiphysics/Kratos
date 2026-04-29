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
///          averaged and tied to newly inserted virtual DoFs.
///          @f[
///                 u_v^o - \frac{1}{m} \sum_{i=1}^m{ \left< u_o^i, n \right> } = 0
///          @f]
///          @f[
///                 u_v^d - \frac{1}{m} \sum_{i=1}^m{ \left< u_d^i, n \right> } = 0
///          @f]
///          where
///          - @f$u_v^o@f$ is the newly inserted virtual DoF belonging to the original set of nodes, and
///          - @f$u_v^d@f$ is the newly inserted virtual DoF belonging to the duplicate set of nodes.
///
///          Once the virtual degrees-of-freedom @f$u_v^o@f$ and @f$u_v^d@f$ are defined, they can either
///          be loaded (Neumann-type) or fixed (Dirichlet-type).
///
/// @subsection pre_tensioning_neumann_type Neumann-Type Pre-Tensioning
///          Given a pre-tensioning force @f$f@f$, fix the reaction of @f$u_v^o@f$ to @f$f@f$ and
///          the reaction of @f$u_v^d@f$ to @f$-f@f$.
///          @f[
///                 r_v^o - f = 0
///          @f]
///          @f[
///                 r_v^d + f = 0
///          @f]
///          where
///          - @f$r_v^0@f$ is the reaction of @f$u_v^o@f$, and
///          - @f$r_v^d@f$ is the reaction of @f$u_v^d@f$.
///
///          @see @ref NeumannPreTensionProcess
///          @see @ref InsertNeumannPreTensionOperation
///
/// @subsection pre_tensioning_dirichlet_type Dirichlet-Type Pre-Tensioning
///          Fix the relative average out-of-plane displacement to a prescribed value @f$\alpha@f$.
///          @f[
///             u_v^o - u_v^d - \alpha = 0
///          @f]
///
///          @see @ref DirichletPreTensionProcess
///          @see @ref InsertDirichletPreTensionOperation
///
/// @subsection pre_tensioning_implementation_notes Implementation Notes
///          Since in-plane displacement constraints reuse DoFs and are dense, they cannot be imposed
///          via master-slave elimination unless DoFs are rotated to line up with the
///          <em>pre-tensioning plane</em>. Since there is currently no robust way of doing this in
///          Kratos, constraints must be imposed via the method of augmented lagrange multipliers,
///          which is implemented in @ref p_multigrid "PMultigridBuilderAndSolver".
///          @note Any analysis involving pre-tensioning must use @ref PMultigridBuilderAndSolver.



/// @brief Base class for @ref InsertDirichletPreTensionOperation and @ref InsertNeumannPreTensionOperation.
/// @see @ref InsertDirichletPreTensionOperation
/// @see @ref InsertNeumannPreTensionOperation
/// @see @ref pre_tensioning "Pre-Tensioning"
/// @ingroup pre_tensioning
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) InsertPreTensionOperation : public Operation {
public:
    KRATOS_CLASS_POINTER_DEFINITION(InsertPreTensionOperation);

    InsertPreTensionOperation(Model& rModel, Parameters Settings);

    InsertPreTensionOperation(InsertPreTensionOperation&&) noexcept;

    ~InsertPreTensionOperation();

    void Execute() override;

    virtual void Apply(double Magnitude) = 0;

    const Parameters GetDefaultParameters() const override;

    virtual std::string Info() const;

protected:
    /// @brief Handle the average out-of-plane displacement.
    virtual void InsertControlNodeConstraints(
        ModelPart& rModelPart,
        array_1d<double,3> SurfaceNormal,
        const std::unordered_map<Node*,Node::Pointer> rDuplicateNodeMap,
        Node::Pointer pControlNode,
        const std::unordered_set<const Dof<double>*> rPositiveSideDofs) = 0;

    struct Impl;
    std::unique_ptr<Impl> mpImpl;

private:
    InsertPreTensionOperation(const InsertPreTensionOperation&) = delete;

    InsertPreTensionOperation& operator=(const InsertPreTensionOperation&) = delete;
}; // class InsertPreTensionOperation


/// @brief Pre-tensioning defined by a surface and a prescribed displacement.
/// @details Cut the mesh along the provided surface and apply fix the average
///             out-of-plane displacement to the provided value while forbidding relative
///             in-plane displacements.
/// @see @ref pre_tensioning "Pre-Tensioning"
/// @ingroup pre_tensioning
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) InsertDirichletPreTensionOperation final : public InsertPreTensionOperation {
public:
    KRATOS_CLASS_POINTER_DEFINITION(InsertDirichletPreTensionOperation);

    using InsertPreTensionOperation::InsertPreTensionOperation;

    void Apply(double Magnitude) override;

    std::string Info() const override;

protected:
    /// @details Inserts a constraint that ties the average out-of-plane
    ///          relative displacement of duplicated nodes to a prescribed value.
    void InsertControlNodeConstraints(
        ModelPart& rModelPart,
        array_1d<double,3> SurfaceNormal,
        const std::unordered_map<Node*,Node::Pointer> rDuplicateNodeMap,
        Node::Pointer pControlNode,
        const std::unordered_set<const Dof<double>*> rPositiveSideDofs) override;

private:
    Node::Pointer mpControlNode;

    std::size_t mPreTensionSurfaceSize;
}; // class InsertDirichletPreTensionOperation


/// @brief Pre-tensioning defined by a surface and a force.
/// @details Cut the mesh along the provided surface and apply a force to
///             out-of-plane displacement components while forbidding relative
///             in-plane displacements.
/// @see @ref pre_tensioning "Pre-Tensioning"
/// @ingroup pre_tensioning
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) InsertNeumannPreTensionOperation final : public InsertPreTensionOperation {
public:
    KRATOS_CLASS_POINTER_DEFINITION(InsertNeumannPreTensionOperation);

    using InsertPreTensionOperation::InsertPreTensionOperation;

    void Apply(double Magnitude) override;

    std::string Info() const override;

protected:
    /// @details Loads the average out-of-plane displacement with a prescribed
    ///          value on both sides of the pre-tension surface in opposite
    ///          directions.
    void InsertControlNodeConstraints(
        ModelPart& rModelPart,
        array_1d<double,3> SurfaceNormal,
        const std::unordered_map<Node*,Node::Pointer> rDuplicateNodeMap,
        Node::Pointer pControlNode,
        const std::unordered_set<const Dof<double>*> rPositiveSideDofs) override;

private:
    Condition::Pointer mpPositiveSideLoad, mpNegativeSideLoad;
}; // class InsertNeumannPreTensionOperation


} // namespace Kratos
