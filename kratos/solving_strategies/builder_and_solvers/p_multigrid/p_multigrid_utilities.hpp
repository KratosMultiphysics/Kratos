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

        /** - @ref GeometryData::KratosGeometryType::Kratos_Line2D2 "Kratos_Line2D2" and @ref GeometryData::KratosGeometryType::Kratos_Line3D2 "Kratos_Line3D2"
         *    linear segment maps to itself.
         *    \f[\begin{bmatrix}
         *      1 &   \\
         *        & 1
         *    \end{bmatrix}\f]
         */
        case G::Kratos_Line3D2:
        case G::Kratos_Line2D2: {
            *itOutput++ = Triplet(0, 0, 1);
            *itOutput   = Triplet(1, 1, 1);
            break;
        } // Kratos_Line2D2 || Kratos_Line3D2

        /** - @ref GeometryData::KratosGeometryType::Kratos_Line2D3 "Kratos_Line2D3" and @ref GeometryData::KratosGeometryType::Kratos_Line3D3 "Kratos_Line3D3"
         *    quadratic line segments reduces to linear line segments
         *    \ref GeometryData::KratosGeometryType::Kratos_Line2D2 "Kratos_Line2D2" and @ref GeometryData::KratosGeometryType::Kratos_Line3D2 "Kratos_Line3D2".
         *    \f[\begin{bmatrix}
         *      1 &   & \frac{1}{2} \\
         *        & 1 & \frac{1}{2}
         *    \end{bmatrix}\f]
         */
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

        /** - @ref GeometryData::KratosGeometryType::Kratos_Triangle2D3 "Kratos_Triangle2D3" and @ref GeometryData::KratosGeometryType::Kratos_Triangle3D3 "Kratos_Triangle3D3"
         *    linear triangles map to themselves.
         *    \f[\begin{bmatrix}
         *      1 &   &   \\
         *        & 1 &   \\
         *        &   & 1
         *    \end{bmatrix}\f]
         */
        case G::Kratos_Triangle2D3:
        case G::Kratos_Triangle3D3: {
            *itOutput++ = Triplet(0, 0, 1);
            *itOutput++ = Triplet(1, 1, 1);
            *itOutput   = Triplet(2, 2, 1);
            break;
        } // Kratos_Triangle2D3 || Kratos_Triangle3D3

        /** - @ref GeometryData::KratosGeometryType::Kratos_Triangle2D6 "Kratos_Triangle2D6" and @ref GeometryData::KratosGeometryType::Kratos_Triangle3D6 "Kratos_Triangle3D6"
         *    quadratic triangles map to linear triangles
         *    @ref GeometryData::KratosGeometryType::Kratos_Triangle2D3 "Kratos_Triangle2D3" and @ref GeometryData::KratosGeometryType::Kratos_Triangle3D3 "Kratos_Triangle3D3".
         *    \f[\begin{bmatrix}
         *      1 &   &   & \frac{1}{2} &   & \frac{1}{2} \\
         *        & 1 &   & \frac{1}{2} & \frac{1}{2} &   \\
         *        &   & 1 &   & \frac{1}{2} & \frac{1}{2}
         *    \end{bmatrix}\f]
         */
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

        /** - @ref GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4 "Kratos_Tetrahedra3D4"
         *    linear tetrahedron maps to itself.
         *    \f[\begin{bmatrix}
         *      1 &   &   &   \\
         *        & 1 &   &   \\
         *        &   & 1 &   \\
         *        &   &   & 1
         *    \end{bmatrix}\f]
         */
        case G::Kratos_Tetrahedra3D4: {
            *itOutput++ = Triplet(0, 0, 1);
            *itOutput++ = Triplet(1, 1, 1);
            *itOutput++ = Triplet(2, 2, 1);
            *itOutput   = Triplet(3, 3, 1);
            break;
        } // Kratos_Tetrahedra3D4

        /** - @ref GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10 "Kratos_Tetrahedra3D10"
         *   quadratic tetrahedron maps to linear tetrahedron @ref GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4 "Kratos_Tetrahedra3D4".
         *   \f[\begin{bmatrix}
         *     1 &   &   &   & \frac{1}{2} &   & \frac{1}{2} & \frac{1}{2} &   &   \\
         *       & 1 &   &   & \frac{1}{2} & \frac{1}{2} &   &   & \frac{1}{2} &   \\
         *       &   & 1 &   &   & \frac{1}{2} & \frac{1}{2} &   &   & \frac{1}{2} \\
         *       &   &   & 1 &   &   &   & \frac{1}{2} & \frac{1}{2} & \frac{1}{2}
         *   \end{bmatrix}\f]
         */
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
/// @param rIndirectDofSet An array of pointers pointing to @ref Dof "dofs" in @p rDofSet.
/// @param rDofMap An index map relating fine DoFs to their coarse counterparts.
/// @warning This function assumes that elements use every @ref Dof of their nodes. This may
///          not always be true (e.g.: coupled analyses on overlapping domains and shared
///          @ref Node "nodes" but different Dofs).
template <unsigned OrderReduction, class TValue>
void MakePRestrictionOperator(ModelPart& rModelPart,
                              const std::size_t FineSystemSize,
                              const PointerVectorSet<Dof<double>>& rParentIndirectDofSet,
                              typename TUblasSparseSpace<TValue>::MatrixType& rRestrictionOperator,
                              const VariablesList::Pointer& rpVariableList,
                              std::vector<std::pair<NodalData,Dof<double>>>& rDofSet,
                              PointerVectorSet<Dof<double>>& rIndirectDofSet,
                              std::vector<std::size_t>& rDofMap)
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

        const auto find_node_index = [&rModelPart](const Node& r_node) -> std::size_t {
            return std::distance(
                rModelPart.Nodes().begin(),
                std::lower_bound(rModelPart.Nodes().begin(),
                                 rModelPart.Nodes().end(),
                                 r_node,
                                 [](const Node& r_left, const Node& r_right){
                                    return r_left.Id() < r_right.Id();
                                 })
            );
        }; // def find_node_index

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
                       [&rows, &locks, &hanging_nodes, &find_node_index, &rParentIndirectDofSet](const Element& r_element, std::vector<Triplet>& r_tls) {
            if (r_element.IsActive()) {
                const auto& r_geometry = r_element.GetGeometry();

                // Fetch the local restriction operator of the element.
                r_tls.clear();
                MakePRestrictionOperator<OrderReduction,TValue,LocalIndex>(r_geometry, std::back_inserter(r_tls));

                for (const auto& [i_row, i_column, value] : r_tls) {
                    auto& r_row_dofs = r_geometry[i_row].GetDofs();
                    const auto& r_column_dofs = r_geometry[i_column].GetDofs();
                    KRATOS_DEBUG_ERROR_IF_NOT(r_row_dofs.size() == r_column_dofs.size());

                    if (!r_row_dofs.empty()) {
                        // Find the first DoF in the parent's set.
                        auto it_first_dof = std::lower_bound(rParentIndirectDofSet.begin(),
                                                             rParentIndirectDofSet.end(),
                                                             **r_row_dofs.begin(),
                                                             [](const Dof<double>& r_parent_dof, const Dof<double>& r_dof){
                                                                return r_parent_dof.Id() < r_dof.Id();
                                                             });
                        if (it_first_dof == rParentIndirectDofSet.end()) continue;

                        // Find which DoFs are included from the node.
                        std::vector<bool> included_dofs(r_row_dofs.size(), false);
                        {
                            const auto i_row_node = r_row_dofs.front()->Id();
                            auto it_parent_dof = it_first_dof;
                            for (std::size_t i_component=0ul; i_component<r_row_dofs.size(); ++i_component) {
                                const Dof<double>& r_row_dof = *r_row_dofs[i_component];
                                if (it_parent_dof->Id() == i_row_node) {
                                    if (it_parent_dof->GetVariable().Key() == r_row_dof.GetVariable().Key()) {
                                        included_dofs[i_component] = true;
                                        ++it_parent_dof;
                                    }
                                } else break;
                            } // for i_component in range(r_row_dofs.size())
                        }

                        // No need acquire the mutexes of all row DoFs, since all
                        // nodes are assumed to have an identical set of DoFs.
                        // ==> set the lock related to the first DoF in the set.
                        [[maybe_unused]] std::scoped_lock<LockObject> row_lock(locks[r_row_dofs.front()->EquationId()]);

                        const unsigned component_count = r_row_dofs.size();
                        for (unsigned i_component=0u; i_component<component_count; ++i_component) {
                            if (!included_dofs[i_component]) continue;
                            Dof<double>* p_fine_row_dof = r_row_dofs[i_component].get();

                            const std::size_t i_row_dof = p_fine_row_dof->EquationId();
                            rows[i_row_dof].second = p_fine_row_dof;

                            const std::size_t i_column_dof = r_column_dofs[i_component]->EquationId();
                            [[maybe_unused]] const auto [it_emplace, inserted] = rows[i_row_dof].first.emplace(i_column_dof, value);

                            // Check whether the insertion was successful, and if not,
                            // make sure that the existing restriction component is
                            // identical to what we tried to insert.
                            KRATOS_DEBUG_ERROR_IF_NOT(inserted || value == it_emplace->second);
                        } // for i_component in range(component_count)
                    } // if !r_row_dofs.empty()
                } // for triplet in r_tls

                // Mark nodes as not hanging.
                for (Node& r_node : r_element.GetGeometry()) {
                    const std::size_t i_node = find_node_index(r_node);
                    KRATOS_DEBUG_ERROR_IF_NOT(i_node < hanging_nodes.size());
                    hanging_nodes[i_node] = 0u;
                } // for r_node in r_element.GetGeometry()
            } // if r_element.IsActive()
        }); // for r_element in rModelPart.Elements()
        KRATOS_CATCH("")

        // Here comes the tricky part.
        // Some models have nodes that are not connected to any element or condition,
        // but have multifreedom constraints to other nodes that are. These connected
        // external DoFs must also be included in the coarse set of DoFs.
        // To find these DoFs, loop through the input DoF set and check whether their
        // nodes are hanging. If they are, they fall into the above mentioned category.
        KRATOS_TRY
        block_for_each(rParentIndirectDofSet, [&hanging_nodes, &locks, &rows, &rModelPart](Dof<double>& r_dof){
            const std::size_t node_id = r_dof.Id();
            const std::size_t i_node = std::distance(
                rModelPart.Nodes().begin(),
                std::lower_bound(rModelPart.Nodes().begin(),
                                 rModelPart.Nodes().end(),
                                 node_id,
                                 [](const Node& r_left, const std::size_t i_node){
                                    return r_left.Id() < i_node;
                                 }));

            // Carry the DoFs of such nodes to the coarse system.
            if (hanging_nodes[i_node]) {
                const std::size_t i_dof = r_dof.EquationId();
                [[maybe_unused]] std::scoped_lock<LockObject> lock(locks[i_dof]);
                rows[i_dof].second = &r_dof;
                rows[i_dof].first.emplace(i_dof, 1.0);
            } // if hanging_nodes[i_node]
        });
        KRATOS_CATCH("")
    } // destroy locks

    // Construct a fine DoF index => coarse DoF index map.
    rDofMap.clear();
    rDofMap.resize(FineSystemSize, std::numeric_limits<std::size_t>::max());

    KRATOS_TRY
    std::size_t i_coarse_dof = 0ul;
    std::size_t entry_count = 0ul;
    for (std::size_t i_fine_dof=0ul; i_fine_dof<FineSystemSize; ++i_fine_dof) {
        entry_count += rows[i_fine_dof].first.size();
        if (!rows[i_fine_dof].first.empty()) {
            rDofMap[i_coarse_dof] = i_fine_dof;
            ++i_coarse_dof;
        } // if not ros[i_fine_dof].first.empty()
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
    KRATOS_TRY
    // Collect solution step variables and store them in a hash map.
    std::unordered_map<typename Variable<double>::KeyType,const Variable<double>*> solution_step_variable_map;

    for (const auto& r_variable_data : rModelPart.GetNodalSolutionStepVariablesList()) {
        const auto key = r_variable_data.Key();
        const std::string& r_name = r_variable_data.Name();
        KRATOS_DEBUG_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(r_name));
        const Variable<double>& r_variable = KratosComponents<Variable<double>>::Get(r_name);
        solution_step_variable_map.emplace(key, &r_variable);
    } // for r_variable_data in rModelPart.GetNodalSolutionStepVariablesList

    rDofSet.clear();
    rIndirectDofSet.clear();

    rDofSet.reserve(rows.size());
    rIndirectDofSet.reserve(rows.size());

    Dof<double>::IndexType i_dof = 0;
    for ([[maybe_unused]] auto& [r_row, rp_dof] : rows) {
        {
            // The CI's rocky linux is using an ancient version of GCC (8.5) that
            // has some bug with non-copyable classes in std::pair, so the following
            // line has to be rewritten in a manner that this fossil understands:
            //rDofSet.emplace_back(NodalData(rp_dof->Id(), rpVariableList), Dof<double>());
            std::pair<NodalData,Dof<double>> entry(1ul, Dof<double>());
            entry.first = NodalData(rp_dof->Id(), rpVariableList);
            rDofSet.emplace_back(std::move(entry));
        }
        rDofSet.back().first.SetSolutionStepData(*rp_dof->GetSolutionStepsData());

        const auto& r_variable_data = rp_dof->GetVariable();
        auto it_variable = solution_step_variable_map.find(r_variable_data.Key());

        if (it_variable == solution_step_variable_map.end()) {
            const auto key = r_variable_data.Key();
            const std::string& r_name = r_variable_data.Name();
            KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(r_name));
            const Variable<double>& r_variable = KratosComponents<Variable<double>>::Get(r_name);
            it_variable = solution_step_variable_map.emplace(key, &r_variable).first;
        }

        const auto& r_reaction_data = rp_dof->GetReaction();
        auto it_reaction = solution_step_variable_map.find(r_reaction_data.Key());

        if (it_reaction == solution_step_variable_map.end()) {
            const auto key = r_reaction_data.Key();
            const std::string& r_name = r_reaction_data.Name();
            KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(r_name));
            const Variable<double>& r_variable = KratosComponents<Variable<double>>::Get(r_name);
            it_reaction = solution_step_variable_map.emplace(key, &r_variable).first;
        }

        rDofSet.back().second = Dof<double>(&rDofSet.back().first,
                                            *it_variable->second,
                                            *it_reaction->second);
        rDofSet.back().second.SetEquationId(i_dof++);

        if (rp_dof->IsFixed()) {
            rDofSet.back().second.FixDof();
        }

        rIndirectDofSet.insert(rIndirectDofSet.end(), &rDofSet.back().second);
    }
    KRATOS_CATCH("")

    // Fill restriction operator row extents.
    // Note: it's a cumulative sum: can't do it in parallel and the standard
    //       doesn't have an algorithm flexible enough for this purpose.
    KRATOS_TRY
    rRestrictionOperator.index1_data()[0] = 0;
    for (std::size_t i_coarse_dof=0ul; i_coarse_dof<rows.size(); ++i_coarse_dof) {
        rRestrictionOperator.index1_data()[i_coarse_dof + 1] = rRestrictionOperator.index1_data()[i_coarse_dof] + rows[i_coarse_dof].first.size();
    } // for i_coarse_dof in range(rows.size())
    KRATOS_CATCH("")

    // Fill the CSR matrix from COO representation.
    KRATOS_TRY
    using TLS = std::vector<std::pair<std::size_t,TValue>>;
    IndexPartition<std::size_t>(rows.size()).for_each(TLS(),
                                                      [&rows, &rRestrictionOperator](const std::size_t i_coarse_dof,
                                                                                     TLS& r_tls) {
        const std::size_t i_entry_begin = rRestrictionOperator.index1_data()[i_coarse_dof];
        [[maybe_unused]] const std::size_t i_entry_end = rRestrictionOperator.index1_data()[i_coarse_dof + 1];

        KRATOS_DEBUG_ERROR_IF(i_entry_end - i_entry_begin != rows[i_coarse_dof].first.size());

        // Copy and sort the current row.
        const auto& r_row = rows[i_coarse_dof].first;
        r_tls.resize(r_row.size());
        std::copy(r_row.begin(), r_row.end(), r_tls.begin());
        std::sort(r_tls.begin(),
                  r_tls.end(),
                  [](const auto& r_left, const auto& r_right){
                    return r_left.first < r_right.first;
                  });

        // Copy row into the CSR matrix.
        std::size_t i_entry = i_entry_begin;
        for (const auto& [i_fine_column, entry] : r_tls) {
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
