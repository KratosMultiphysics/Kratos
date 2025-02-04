//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// Project includes
#include "geometries/geometry_data.h" // GeometryData
#include "geometries/geometry.h" // Geometry
#include "spaces/ublas_space.h" // TUblasSparseSpace
#include "includes/lock_object.h" // LockObject

// System includes
#include <tuple> // std::tuple
#include <vector> // std::vector
#include <unordered_map> // std::unordered_map
#include <limits> // std::numeric_limits
#include <cstdint> // std::uint8_t
#include <algorithm> // std::remove_if


namespace Kratos {


namespace detail {
template <class TTable>
struct UnionReduction
{
    using value_type = TTable;

    using return_type = TTable;

    return_type mValue;

    return_type&& GetValue()
    {
        return std::move(mValue);
    }

    void LocalReduce(const value_type& rValue)
    {
        mValue.insert(rValue.begin(), rValue.end());
    }

    void ThreadSafeReduce(UnionReduction& rOther)
    {
        KRATOS_CRITICAL_SECTION
        {
            for (auto& r_item : rOther) {
                mValue.insert(std::move(r_item));
            }
        }
    }
}; // class UnionReduction
} // namespace detail


/// @brief Compute the p-restriction operator of a single geometry.
/// @tparam OrderReduction Defines how many polynomial order to reduce the incoming
///         geometry, at the maximum. If 0, the function returns a correctly sized
///         identity matrix. Reduction to linear geometries will always happen if
///         @p OrderReduction is greater or equal than the geometry's order.
/// @tparam TValue Value type of entries.
/// @tparam TIndex
/// @tparam TNode Node type (usually @ref Node).
/// @tparam TOutputIterator Output iterator of triplets. Nonzero entries in the restriction operator
///         are provided as @p std::tuple<unsigned,unsigned,TValue> triplets of {coarse node index, fine node index, value}.
/// @param rGeometry Geometry to construct a restriction operator for.
/// @param itOutput Output iterator where nonzero entries in the restriction operator are written to. Entries are
///                 represented as {coarse node index, fine node index, value} triplets (COO matrix).
/// @details This function computes the (local) restriction operator @f$ R @f$ for a geometry that
///          transforms the input (higher order @f$ q @f$) geometry's shape functions to an equivalent
///          lower order (@f$ p @f$). Due to the linear relationship between shape functions and
///          degrees-of-freedom, this mapping holds for DoFs as well:
///          @f[
///             u^p_i = R_{ij} u^q_j
///          @f]
///          The reduction of the input geometry's order is defined by `OrderReduction` and is limited
///          to reducing to a linear geometry. Any `OrderReduction` greater than that will only reduce
///          to linear geometries as well.
template <unsigned OrderReduction,
          class TValue,
          class TIndex,
          class TNode,
          class TOutputIterator>
void MakePRestrictionOperator(const Geometry<TNode>& rGeometry,
                              TOutputIterator itOutput)
{
    using G = GeometryData::KratosGeometryType;
    const G geometry_type = rGeometry.GetGeometryType();

    using Triplet = std::tuple<TIndex,TIndex,TValue>;

    switch(geometry_type) {
        /// - @ref GeometryData::KratosGeometryType::Kratos_Point2D "Kratos_Point2D" and @ref GeometryData::KratosGeometryType::Kratos_Point3D "Kratos_Point3D"
        ///   can't be higher order, so they just map back to themselves.
        ///   \f[\begin{bmatrix}
        ///     1
        ///   \end{bmatrix}\f]
        case G::Kratos_Point2D:
        case G::Kratos_Point3D: {
            *itOutput = Triplet(0, 0, 1);
            break;
        } // Kratos_Point2D || Kratos_Point3D

        /// - @ref GeometryData::KratosGeometryType::Kratos_Line2D2 "Kratos_Line2D2" and @ref GeometryData::KratosGeometryType::Kratos_Line3D2 "Kratos_Line3D2"
        ///   linear segment maps to itself.
        ///   \f[\begin{bmatrix}
        ///     1 &   \\
        ///       & 1
        ///   \end{bmatrix}\f]
        case G::Kratos_Line3D2:
        case G::Kratos_Line2D2: {
            *itOutput++ = Triplet(0, 0, 1);
            *itOutput   = Triplet(1, 1, 1);
            break;
        } // Kratos_Line2D2 || Kratos_Line3D2

        /// - @ref GeometryData::KratosGeometryType::Kratos_Line2D3 "Kratos_Line2D3" and @ref GeometryData::KratosGeometryType::Kratos_Line3D3 "Kratos_Line3D3"
        ///   quadratic line segments reduces to linear line segments
        ///   \ref GeometryData::KratosGeometryType::Kratos_Line2D2 "Kratos_Line2D2" and @ref GeometryData::KratosGeometryType::Kratos_Line3D2 "Kratos_Line3D2".
        ///   \f[\begin{bmatrix}
        ///     1 &   & \frac{1}{2} \\
        ///       & 1 & \frac{1}{2}
        ///   \end{bmatrix}\f]
        case G::Kratos_Line3D3:
        case G::Kratos_Line2D3: {
            if constexpr (OrderReduction == 0) {
                *itOutput++ = Triplet(0, 0, 1.0);
                *itOutput++ = Triplet(1, 1, 1.0);
                *itOutput   = Triplet(2, 2, 1.0);
            } /*if OrderReduction == 0*/ else {
                *itOutput++ = Triplet(0, 0, 1.0);
                *itOutput++ = Triplet(0, 2, 0.5);
                *itOutput++ = Triplet(1, 1, 1.0);
                *itOutput   = Triplet(1, 2, 0.5);
            } // if OrderReduction != 0
            break;
        } // Kratos_Line2D3 || Kratos_Line3D3

        /// - @ref GeometryData::KratosGeometryType::Kratos_Triangle2D3 "Kratos_Triangle2D3" and @ref GeometryData::KratosGeometryType::Kratos_Triangle3D3 "Kratos_Triangle3D3"
        ///   linear triangles map to themselves.
        ///   \f[\begin{bmatrix}
        ///     1 &   &   \\
        ///       & 1 &   \\
        ///       &   & 1
        ///   \end{bmatrix}\f]
        case G::Kratos_Triangle2D3:
        case G::Kratos_Triangle3D3: {
            *itOutput++ = Triplet(0, 0, 1);
            *itOutput++ = Triplet(1, 1, 1);
            *itOutput   = Triplet(2, 2, 1);
            break;
        } // Kratos_Triangle2D3 || Kratos_Triangle3D3

        /// - @ref GeometryData::KratosGeometryType::Kratos_Triangle2D6 "Kratos_Triangle2D6" and @ref GeometryData::KratosGeometryType::Kratos_Triangle3D6 "Kratos_Triangle3D6"
        ///   quadratic triangles map to linear triangles
        ///   @ref GeometryData::KratosGeometryType::Kratos_Triangle2D3 "Kratos_Triangle2D3" and @ref GeometryData::KratosGeometryType::Kratos_Triangle3D3 "Kratos_Triangle3D3".
        ///   \f[\begin{bmatrix}
        ///     1 &   &   & \frac{1}{2} &   & \frac{1}{2} \\
        ///       & 1 &   & \frac{1}{2} & \frac{1}{2} &   \\
        ///       &   & 1 &   & \frac{1}{2} & \frac{1}{2}
        ///   \end{bmatrix}\f]
        case G::Kratos_Triangle2D6:
        case G::Kratos_Triangle3D6: {
            if constexpr (OrderReduction == 0) {
                *itOutput++ = Triplet(0, 0, 1.0);
                *itOutput++ = Triplet(1, 1, 1.0);
                *itOutput++ = Triplet(2, 2, 1.0);
                *itOutput++ = Triplet(3, 3, 1.0);
                *itOutput++ = Triplet(4, 4, 1.0);
                *itOutput   = Triplet(5, 5, 1.0);
            } /*if OrderReduction == 0*/ else {
                *itOutput++ = Triplet(0, 0, 1.0);
                *itOutput++ = Triplet(0, 3, 0.5);
                *itOutput++ = Triplet(0, 5, 0.5);
                *itOutput++ = Triplet(1, 1, 1.0);
                *itOutput++ = Triplet(1, 3, 0.5);
                *itOutput++ = Triplet(1, 4, 0.5);
                *itOutput++ = Triplet(2, 2, 1.0);
                *itOutput++ = Triplet(2, 4, 0.5);
                *itOutput   = Triplet(2, 5, 0.5);
            } // if OrderReduction != 0
            break;
        } // Kratos_Triangle2D6 || Kratos_Triangle3D6

        /// - @ref GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4 "Kratos_Tetrahedra3D4"
        ///   linear tetrahedron maps to itself.
        ///   \f[\begin{bmatrix}
        ///     1 &   &   &   \\
        ///       & 1 &   &   \\
        ///       &   & 1 &   \\
        ///       &   &   & 1
        ///   \end{bmatrix}\f]
        case G::Kratos_Tetrahedra3D4: {
            *itOutput++ = Triplet(0, 0, 1);
            *itOutput++ = Triplet(1, 1, 1);
            *itOutput++ = Triplet(2, 2, 1);
            *itOutput   = Triplet(3, 3, 1);
            break;
        } // Kratos_Tetrahedra3D4

        /// - @ref GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10 "Kratos_Tetrahedra3D10"
        ///   quadratic tetrahedron maps to linear tetrahedron @ref GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4 "Kratos_Tetrahedra3D4".
        ///   \f[\begin{bmatrix}
        ///     1 &   &   &   & \frac{1}{2} &   & \frac{1}{2} & \frac{1}{2} &   &   \\
        ///       & 1 &   &   & \frac{1}{2} & \frac{1}{2} &   &   & \frac{1}{2} &   \\
        ///       &   & 1 &   &   & \frac{1}{2} & \frac{1}{2} &   &   & \frac{1}{2} \\
        ///       &   &   & 1 &   &   &   & \frac{1}{2} & \frac{1}{2} & \frac{1}{2}
        ///   \end{bmatrix}\f]
        case G::Kratos_Tetrahedra3D10: {
            if constexpr (OrderReduction == 0) {
                *itOutput++ = Triplet(0, 0, 1.0);
                *itOutput++ = Triplet(1, 1, 1.0);
                *itOutput++ = Triplet(2, 2, 1.0);
                *itOutput++ = Triplet(3, 3, 1.0);
                *itOutput++ = Triplet(4, 4, 1.0);
                *itOutput++ = Triplet(5, 5, 1.0);
                *itOutput++ = Triplet(6, 6, 1.0);
                *itOutput++ = Triplet(7, 7, 1.0);
                *itOutput++ = Triplet(8, 8, 1.0);
                *itOutput   = Triplet(9, 9, 1.0);
            } /*if OrderReduction == 0*/ else {
                *itOutput++ = Triplet(0, 0, 1.0);
                *itOutput++ = Triplet(0, 4, 0.5);
                *itOutput++ = Triplet(0, 6, 0.5);
                *itOutput++ = Triplet(0, 7, 0.5);
                *itOutput++ = Triplet(1, 1, 1.0);
                *itOutput++ = Triplet(1, 4, 0.5);
                *itOutput++ = Triplet(1, 5, 0.5);
                *itOutput++ = Triplet(1, 8, 0.5);
                *itOutput++ = Triplet(2, 2, 1.0);
                *itOutput++ = Triplet(2, 5, 0.5);
                *itOutput++ = Triplet(2, 6, 0.5);
                *itOutput++ = Triplet(2, 9, 0.5);
                *itOutput++ = Triplet(3, 3, 1.0);
                *itOutput++ = Triplet(3, 7, 0.5);
                *itOutput++ = Triplet(3, 8, 0.5);
                *itOutput   = Triplet(3, 9, 0.5);
            }
            break;
        } // Kratos_Tetrahedra3D10

        default: KRATOS_ERROR << "unsupported geometry: " << rGeometry.Name() << "\n";
    } // switch geometry_type
}


/// @brief Compute the p-multigrid restriction operator for the provided mesh.
/// @tparam TValue Number type of stored values in the sparse matrix.
/// @tparam OrderReduction
/// @param rModelPart @ref ModelPart representing the system to be solved.
/// @param FineSystemSize Number of rows in the fine system's left hand side matrix.
/// @param rRestrictionOperator Output matrix to which the restriction operator will be written.
/// @param rDofSet @ref Output vector containing the @ref Dof "DoFs" of the coarse system.
/// @warning This function assumes that elements use every @ref Dof of their nodes. This may
///          not always be true (e.g.: coupled analyses on overlapping domains and shared
///          @ref Node "nodes" but different Dofs).
template <unsigned OrderReduction,
          class TValue>
void MakePRestrictionOperator(ModelPart& rModelPart,
                              const std::size_t FineSystemSize,
                              typename TUblasSparseSpace<TValue>::MatrixType& rRestrictionOperator,
                              std::vector<Dof<double>>& rDofSet,
                              PointerVectorSet<Dof<double>>& rIndirectDofSet)
{
    static_assert(OrderReduction == std::numeric_limits<unsigned>::max(),
                  "The current implementation requires geometries to be always reduced to their linear equivalents in a single step.");

    KRATOS_TRY

    using GlobalIndex = typename TUblasSparseSpace<TValue>::MatrixType::size_type;
    using LocalIndex = std::uint8_t;

    // Construct an extended representation of the restriction operator.
    // Each entry in the array represents a row in the matrix, with the
    // stored map defining its nonzero entries. Empty rows are ignored and
    // won't be in the final operator.
    // Since the final coarse DoF index assignment can only happen once we
    // know which rows have entries and which don't, we can't assign global
    // column indices while constructing the restriction operator. Instead,
    // fine row indices are stored, from which coarse row indices can later
    // be computed.
    std::vector<std::pair<
        std::unordered_map<
            GlobalIndex,    // <== column index
            TValue          // <== value
        >,
        Dof<double>*  // <== pointer to the related DoF on the fine grid
    >> rows(FineSystemSize);

    // Construct a COO representation of the restriction operator.
    {
        std::vector<LockObject> locks(FineSystemSize);
        using Triplet = std::tuple<LocalIndex,LocalIndex,TValue>;

        // A vector of bools indicating whether the node at the same position
        // is a hanging node (i.e.: not part of any element/condition) or not.
        // In case you're wondering, I'm not using an std::vector<std::atomic<bool>>
        // because I don't want someone to change it into a vector of naked
        // bools in the future, which would lead to catastrophe. Bools and
        // std::uint8_t s are memory and performance-wise identical anyway.
        std::vector<std::atomic<std::uint8_t>> hanging_nodes(rModelPart.Nodes().size());
        block_for_each(hanging_nodes, [](auto& r_flag) {r_flag = 1;});

        // Collect terms from elements.
        KRATOS_TRY
        block_for_each(rModelPart.Elements(),
                       std::vector<Triplet>(),
                       [&rows, &locks, &hanging_nodes, &rModelPart](const Element& r_element, std::vector<Triplet>& r_tls) {
            if (r_element.IsActive()) {
                r_tls.clear();
                const auto& r_geometry = r_element.GetGeometry();

                // Fetch the local restriction operator of the element.
                MakePRestrictionOperator<OrderReduction,TValue,LocalIndex>(r_geometry, std::back_inserter(r_tls));

                for (const auto& [i_row, i_column, value] : r_tls) {
                    auto& r_row_dofs = r_geometry[i_row].GetDofs();
                    const auto& r_column_dofs = r_geometry[i_column].GetDofs();
                    KRATOS_ERROR_IF_NOT(r_row_dofs.size() == r_column_dofs.size());

                    if (!r_row_dofs.empty()) {
                        // No need acquire the mutexes of all row DoFs, since all
                        // nodes are assumed to have an identical set of DoFs.
                        // ==> set the lock related to the first DoF in the set.
                        [[maybe_unused]] std::scoped_lock<LockObject> row_lock(locks[r_row_dofs.front()->EquationId()]);

                        const unsigned component_count = r_row_dofs.size();
                        for (unsigned i_component=0u; i_component<component_count; ++i_component) {
                            Dof<double>* p_fine_row_dof = r_row_dofs[i_component].get();
                            const std::size_t i_row_dof = p_fine_row_dof->EquationId();
                            rows[i_row_dof].second = p_fine_row_dof;

                            const std::size_t i_column_dof = r_column_dofs[i_component]->EquationId();
                            [[maybe_unused]] const auto [it_emplace, inserted] = rows[i_row_dof].first.emplace(i_column_dof, value);

                            // Check whether the insertion was successful, and if not,
                            // make sure that the existing restriction component is
                            // identical to what we tried to insert.
                            KRATOS_DEBUG_ERROR_IF_NOT(!inserted || value == it_emplace->second);
                        } // for i_component in range(component_count)
                    } // if !r_row_dofs.empty()
                } // for triplet in r_tls

            } // if r_element.IsActive()

            // Mark nodes as not hanging.
            for (const Node& r_node : r_element.GetGeometry()) {
                const std::size_t i_node = std::distance(
                    rModelPart.Nodes().begin(),
                    std::lower_bound(rModelPart.Nodes().begin(),
                                     rModelPart.Nodes().end(),
                                     r_node,
                                     [](const Node& r_left, const Node& r_right){
                                        return r_left.Id() < r_right.Id();
                                     }));
                KRATOS_DEBUG_ERROR_IF_NOT(i_node < hanging_nodes.size());
                hanging_nodes[i_node] = 0u;
            } // for r_node in r_element.GetGeometry()
        }); // for r_element in rModelPart.Elements()
        KRATOS_CATCH("")

        // Flag conditions' nodes as not hanging.
        KRATOS_TRY
        block_for_each(rModelPart.Conditions(), [&hanging_nodes, &rModelPart](const Condition& r_condition) {
            for (const Node& r_node : r_condition.GetGeometry()) {
                const std::size_t i_node = std::distance(
                    rModelPart.Nodes().begin(),
                    std::lower_bound(rModelPart.Nodes().begin(),
                                     rModelPart.Nodes().end(),
                                     r_node,
                                     [](const Node& r_left, const Node& r_right){
                                        return r_left.Id() < r_right.Id();
                                     }));
                hanging_nodes[i_node] = 0u;
            } // for r_node in r_condition.GetGeometry()
        }); // for r_condition in rModelPart.Conditions()
        KRATOS_CATCH("")

        // Include hanging nodes' DoFs.
        KRATOS_TRY
        IndexPartition<std::size_t>(hanging_nodes.size()).for_each([&hanging_nodes, &rows, &locks, &rModelPart](const std::size_t i_node){
            if (hanging_nodes[i_node]) {
                Node& r_node = rModelPart.Nodes().begin()[i_node];
                for (auto& rp_dof : r_node.GetDofs()) {
                    const std::size_t i_dof = rp_dof->EquationId();
                    rows[i_dof].second = rp_dof.get();
                    [[maybe_unused]] std::scoped_lock<LockObject> lock(locks[i_dof]);
                    if (rows[i_dof].first.empty()) {
                        rows[i_dof].first.emplace(i_dof, 1.0);
                    } // if rows[i_dof].first.empty()
                } // for rp_dof in r_node.GetDofs()
            } // if hanging_nodes[i_node]
        }); // for i_node in range(hanging_nodes.size())
        KRATOS_CATCH("")
    } // destroy locks

    // Construct a fine DoF index => coarse DoF index map.
    std::vector<std::size_t> dof_map(FineSystemSize);

    {
        KRATOS_TRY
        std::size_t i_coarse_dof = 0ul;
        std::size_t entry_count = 0ul;
        for (std::size_t i_fine_dof=0ul; i_fine_dof<FineSystemSize; ++i_fine_dof) {
            dof_map[i_fine_dof] = rows[i_fine_dof].first.empty() ? i_coarse_dof : i_coarse_dof++;
            entry_count += rows[i_fine_dof].first.size();
        } // for i_fine_dof in range(FineSystemSize)

        // No need to keep the higher order rows (i.e.: the empty rows) anymore => erase them.
        rows.erase(std::remove_if(rows.begin(),
                                  rows.end(),
                                  [](const auto& r_pair) {return r_pair.first.empty();}),
                   rows.end());
        KRATOS_ERROR_IF_NOT(rows.size() == i_coarse_dof);

        // Resize the restriction operator.
        rRestrictionOperator = typename TUblasSparseSpace<TValue>::MatrixType(i_coarse_dof, FineSystemSize, entry_count);
        KRATOS_CATCH("")
    }

    // todo
    // Normally, the output DofSet would be filled here but since constructing
    // Dofs is a pain in the ass and the current implementation of constraint assemblers
    // don't need it, I'll skip this task for now. I'm leaving the DofSet in the arguments though
    // because it will eventually become necessary when imposition with lagrange multipliers
    // (not augmented ones) is implemented.
    // Filling the Dof pointers is relatively easy with the current implementation though.
    // This function only allows a direct restriction to a linear mesh, and linear geometries'
    // nodes (i.e.: corner nodes) are always a strict subset of their high order counterparts.
    // This means that fine Dofs can be reused for the coarse system.
    rDofSet.clear();
    rIndirectDofSet.clear();
    rIndirectDofSet.reserve(rows.size());
    for (const auto& [r_row, rp_dof] : rows) {
        rIndirectDofSet.insert(rIndirectDofSet.end(), rp_dof);
    }

    // Fill restriction operator row extents.
    // Note: it's a cumulative sum: can't do it in parallel and the standard
    //       doesn't have an algorithm flexible enough for this purpose.
    KRATOS_TRY
    rRestrictionOperator.index1_data()[0] = 0;
    for (std::size_t i_coarse_dof=0ul; i_coarse_dof<rows.size(); ++i_coarse_dof) {
        rRestrictionOperator.index1_data()[i_coarse_dof + 1] = rRestrictionOperator.index1_data()[i_coarse_dof] + rows[i_coarse_dof].first.size();
    } // for i_coarse_dof in range(rows.size())
    KRATOS_CATCH("")

    // Fill CSR from COO representation.
    KRATOS_TRY
    IndexPartition<std::size_t>(rows.size()).for_each([&rows, &rRestrictionOperator, &dof_map](const std::size_t i_coarse_dof) {
        const std::size_t i_entry_begin = rRestrictionOperator.index1_data()[i_coarse_dof];
        [[maybe_unused]] const std::size_t i_entry_end = rRestrictionOperator.index1_data()[i_coarse_dof + 1];

        KRATOS_DEBUG_ERROR_IF(i_entry_end - i_entry_begin != rows[i_coarse_dof].first.size());

        const auto& r_row = rows[i_coarse_dof].first;
        std::size_t i_entry = i_entry_begin;
        for (const auto [i_fine_column, entry] : r_row) {
            KRATOS_DEBUG_ERROR_IF_NOT(i_entry < rRestrictionOperator.index2_data().size()) << i_entry;
            KRATOS_DEBUG_ERROR_IF_NOT(i_entry < rRestrictionOperator.value_data().size()) << i_entry;
            rRestrictionOperator.index2_data()[i_entry] = i_fine_column;
            rRestrictionOperator.value_data()[i_entry] = entry;
            ++i_entry;
        } // for i_fine_column, entry in r_row
    }); // for i_coarse_dof in range(rows.size())
    KRATOS_CATCH("")

    KRATOS_TRY
    rRestrictionOperator.set_filled(rRestrictionOperator.index1_data().size(),
                                    rRestrictionOperator.index2_data().size());
    KRATOS_CATCH("")

    KRATOS_CATCH("")
}


} // namespace Kratos
