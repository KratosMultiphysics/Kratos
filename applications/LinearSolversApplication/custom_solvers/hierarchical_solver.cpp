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

// Project includes
#include "hierarchical_solver.h"
#include "includes/code_location.h"
#include "includes/global_variables.h"
#include "spaces/ublas_space.h"
#include "utilities/profiler.h"
#include "factories/linear_solver_factory.h"
#include "utilities/proxies.h"
#include "utilities/sparse_matrix_multiplication_utility.h"
#include "utilities/builtin_timer.h"

// STL includes
#include <type_traits> // remove_reference_t, remove_cv_t
#include <sstream> // stringstream


namespace Kratos {


namespace detail {


class ParametersWithDefaults
{
public:
    ParametersWithDefaults(Parameters Settings,
                           Parameters Defaults)
        : mSettings(Settings),
          mDefaults(Defaults)
    {}

    ParametersWithDefaults operator[](const std::string& rEntryName)
    {
        if (mSettings.Has(rEntryName)) {
            if (mDefaults.Has(rEntryName)) {
                return ParametersWithDefaults(mSettings[rEntryName], mDefaults[rEntryName]);
            } else {
                KRATOS_ERROR << "Input parameters has entry '" << rEntryName
                             << "' but the defaults do not. Input parameters:\n"
                             << mSettings << "\nDefault parameters:\n"
                             << mDefaults;
            }
        } else {
            if (mDefaults.Has(rEntryName)) {
                Parameters defaults = mDefaults[rEntryName];
                mSettings.AddValue(rEntryName, defaults);
                return ParametersWithDefaults(mSettings[rEntryName], defaults);
            } else {
                KRATOS_ERROR << "Neither input parameters, nor defaults have an entry for "
                             << "'" << rEntryName << "'. Input parameters:\n"
                             << mSettings << "\nDefault parameters:\n"
                             << mDefaults;
            } // else (mDefaults.Has(rEntryName))
        } // else (mSettings.Has(rEntryName))
    }

    template <class T>
    auto Get()
    {
        using TValue = std::remove_cv_t<std::remove_reference_t<T>>;

        if (!mSettings.IsNull()) {
            return mSettings.Get<TValue>();
        }

        if (!mDefaults.IsNull()) {
            auto output = mDefaults.Get<TValue>();
            mSettings.Set<TValue>(output);
            return output;
        }

        KRATOS_ERROR << "Both the input parameters and the defaults are empty.\n";
        return TValue();
    }

    Parameters GetInput()
    {
        return mSettings;
    }

    Parameters GetDefaults()
    {
        return mDefaults;
    }


private:
    Parameters mSettings;

    Parameters mDefaults;
};


} // namespace detail



template <class TSparseSpace,
          class TDenseSpace,
          class TReorderer>
struct HierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::Impl
{
    double mTolerance;

    int mVerbosity;

    std::size_t mMaxIterations;

    typename LinearSolver<TSparseSpace,TDenseSpace/*reorderer is omitted on purpose*/>::Pointer mpFineSolver;

    typename TSparseSpace::MatrixType mRestrictionOperator;

    typename TSparseSpace::MatrixType mInterpolationOperator;

    struct CoarseSystem {
        typename TSparseSpace::MatrixType mA;

        typename LinearSolver<TSparseSpace,TDenseSpace/*reorderer is omitted on purpose*/>::Pointer mpSolver;
    } mCoarseSystem;
}; // struct HierarchicalSolver::Impl



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
HierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::HierarchicalSolver()
    : HierarchicalSolver(HierarchicalSolver::GetDefaultParameters())
{
}



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
HierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::HierarchicalSolver(Parameters parameters)
    : mpImpl(new Impl)
{
    KRATOS_TRY
    Parameters default_parameters = HierarchicalSolver::GetDefaultParameters();
    parameters.ValidateAndAssignDefaults(default_parameters);

    KRATOS_ERROR_IF_NOT(parameters["solver_type"].GetString() == "hierarchical")
        << "Requested a(n) '" << parameters["solver_type"].GetString() << "' solver,"
        << " but constructing a HierarchicalSolver";

    detail::ParametersWithDefaults safe_parameters(parameters, default_parameters);
    mpImpl->mTolerance = safe_parameters["tolerance"].Get<double>();
    mpImpl->mVerbosity = safe_parameters["verbosity"].Get<int>();
    mpImpl->mMaxIterations = safe_parameters["max_iterations"].Get<int>();

    // Construct nested solvers
    KRATOS_TRY
    LinearSolverFactory<TSparseSpace,TDenseSpace> solver_factory;
    mpImpl->mCoarseSystem.mpSolver = solver_factory.Create(safe_parameters["coarse_solver_settings"].GetInput());
    mpImpl->mpFineSolver = solver_factory.Create(safe_parameters["fine_solver_settings"].GetInput());
    KRATOS_CATCH("")

    KRATOS_CATCH("")
}



// Necessary for PIMPL
template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
HierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::~HierarchicalSolver()
{
}



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
bool HierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::Solve(SparseMatrix& rA,
                                                                    Vector& rX,
                                                                    Vector& rB)
{
    KRATOS_TRY
    KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);

    // Declarations
    const std::size_t coarse_size = TSparseSpace::Size1(mpImpl->mRestrictionOperator);
    const std::size_t fine_size = TSparseSpace::Size2(mpImpl->mRestrictionOperator);
    Vector fine_delta(fine_size);
    Vector fine_residual(fine_size);
    Vector fine_tmp(fine_size);
    Vector coarse_delta(coarse_size);
    Vector coarse_residual(coarse_size);

    const double rhs_norm = TSparseSpace::TwoNorm(rB);

    // Initializations
    // r = b - A * x
    TSparseSpace::Mult(rA, rX, fine_tmp);
    noalias(fine_residual) = rB - fine_tmp;

    double residual_norm = 1.0;
    const std::size_t max_iterations = mpImpl->mMaxIterations;
    std::size_t i_iteration = 0ul;
    while (true) {
        // Restrict the residual
        // r_coarse = R * r_fine
        TSparseSpace::Mult(mpImpl->mRestrictionOperator,
                           fine_residual,
                           coarse_residual);

        // Solve coarse system
        // A_coarse * d_coarse = r_coarse
        KRATOS_TRY
        TSparseSpace::SetToZero(coarse_delta);
        mpImpl->mCoarseSystem.mpSolver->Solve(mpImpl->mCoarseSystem.mA,
                                              coarse_delta,
                                              coarse_residual);
        TSparseSpace::Mult(mpImpl->mInterpolationOperator, coarse_delta, fine_delta);
        TSparseSpace::UnaliasedAdd(rX, 1.0, fine_delta);
        //KRATOS_INFO("after preconditioning: ") << rX << "\n";
        KRATOS_CATCH("")

        // Update the fine residual
        // r_fine -= A_fine * d_fine
        TSparseSpace::Mult(rA, fine_delta, fine_tmp);
        //fine_residual -= fine_tmp;
        TSparseSpace::UnaliasedAdd(fine_residual, -1.0, fine_tmp);

        KRATOS_INFO_IF("HierarchicalSolver", 3 <= mpImpl->mVerbosity)
            << "iteration " << i_iteration
            << " residual after coarse preconditioning "
            << TSparseSpace::TwoNorm(fine_residual) / rhs_norm
            << "\n";

        // Perform a few smoothing passes on the fine system
        // A_fine * d_fine = r_fine
        KRATOS_TRY
        TSparseSpace::SetToZero(fine_delta);
        mpImpl->mpFineSolver->Solve(rA, fine_delta, fine_residual);
        TSparseSpace::UnaliasedAdd(rX, 1.0, fine_delta);
        //KRATOS_INFO("after smoothing: ") << rX << "\n";
        KRATOS_CATCH("")

        // Update the fine residual
        // r_fine = b_fine - A_fine * x_fine
        TSparseSpace::Mult(rA, fine_delta, fine_tmp);
        //fine_residual = rB - fine_tmp;
        TSparseSpace::UnaliasedAdd(fine_residual, -1.0, fine_tmp);

        // Update status
        residual_norm = TSparseSpace::TwoNorm(fine_residual) / rhs_norm;
        KRATOS_INFO_IF("HierarchicalSolver", 1 <= mpImpl->mVerbosity)
            << "iteration " << i_iteration
            << " residual " << residual_norm << "\n";

        // Early exit if converged
        if (residual_norm < mpImpl->mTolerance || max_iterations <= ++i_iteration) {
            break;
        }
    } // while (true)

    if (4 <= mpImpl->mVerbosity) {
        KRATOS_INFO("HierarchicalSolver") << "writing solution.mm\n";
        TSparseSpace::WriteMatrixMarketVector("solution.mm", rX);
    }

    return residual_norm < mpImpl->mTolerance;
    KRATOS_CATCH("")
}


std::size_t GetLinearGeometryVertexCount(GeometryData::KratosGeometryFamily GeometryFamily)
{
    using Family = GeometryData::KratosGeometryFamily;
    switch(GeometryFamily) {
        case Family::Kratos_Point:          return 1;
        case Family::Kratos_Linear:         return 2;
        case Family::Kratos_Triangle:       return 3;
        case Family::Kratos_Quadrilateral:  return 4;
        case Family::Kratos_Tetrahedra:     return 4;
        case Family::Kratos_Pyramid:        return 5;
        case Family::Kratos_Hexahedra:      return 6;
        default: {
            KRATOS_ERROR << "unsupported geometry family";
        }
    } // switch (GeometryFamily)
}


template <Globals::DataLocation TEntityLocation, class TSparseSpace>
void GetLowerOrderDofs(typename TSparseSpace::MatrixType& rA,
                       const ModelPart& rModelPart,
                       const ModelPart::DofsArrayType& rDofs,
                       std::vector<bool>& rCoarseMask)
{
    KRATOS_TRY
    // Loop through elements and collect the DoFs' IDs
    // that are related to corner vertices.
    // @todo the parallel version of this loop relies on undefined behaviour,
    //       because all threads write to the same coarse mask at once => implement
    //       an ANY reduction. @matekelemen
    for (auto proxy : MakeProxy<TEntityLocation>(rModelPart)) {
        const auto& r_entity = proxy.GetEntity();
        const Geometry<Node>& r_geometry = r_entity.GetGeometry();

        if (r_entity.IsActive() && r_geometry.PointsNumber()) {
            // Get the number of vertices on the equivalent lower order element
            const std::size_t coarse_node_count = GetLinearGeometryVertexCount(r_geometry.GetGeometryFamily());
            KRATOS_DEBUG_ERROR_IF_NOT(coarse_node_count <= r_geometry.PointsNumber())
                << "coarse node count " << coarse_node_count
                << " is higher than the number of fine vertices " << r_geometry.PointsNumber();

            // Mark each DoF of all nodes on the equivalent lower order element
            // (unless they are constrained as slave dofs in multi-point constraints)
            for (std::size_t i_node=0ul; i_node<coarse_node_count; ++i_node) {
                if (!r_geometry[i_node].Is(SLAVE)) {
                    const auto& r_node_dofs = r_geometry[i_node].GetDofs();
                    for (const auto& rp_dof : r_node_dofs) {
                        const std::size_t i_fine_dof = rp_dof->EquationId();
                        KRATOS_DEBUG_ERROR_IF_NOT(i_fine_dof < TSparseSpace::Size1(rA));

                        // Flag the DoF as part of the coarse system
                        if (rp_dof->IsFree()) rCoarseMask[i_fine_dof] = true;
                    }
                } // if the node is not flagged SLAVE
            }
        } // if (r_entity.IsActive())
    } // for entity in rModelPart

    // Loop through the active DoFs and remove any
    // slave or dirichlet DoFs from the coarse system.
    for (const Dof<double>& r_dof : rDofs) {
        if (r_dof.IsFixed()) {
            rCoarseMask[r_dof.EquationId()] = false;
        } else {
            const auto node_id = r_dof.Id();
            const auto it_node = rModelPart.Nodes().find(node_id);
            KRATOS_ERROR_IF(it_node == rModelPart.Nodes().end())
                << "DoF " << r_dof.EquationId() << " has no corresponding node in model part "
                << rModelPart.Name() << std::endl;

            if (it_node->Is(SLAVE)) {
                rCoarseMask[r_dof.EquationId()] = false;
            }
        }
    }
    KRATOS_CATCH("")
}


namespace Detail {


struct IndexPairTraits
{
    struct IsLess {
        bool operator()(std::pair<std::size_t,std::size_t> Left,
                        std::pair<std::size_t,std::size_t> Right) const noexcept
        {
            if (Left.first < Right.first) {
                return true;
            } else if (Left.first == Right.first) {
                return Left.second < Right.second;
            } else {
                return false;
            }
        }
    }; // struct IsLess
}; // struct IndexPairTraits


} // namespace Detail



/// @brief Get the relationship between lower and higher order DoFs for the provided @ref Geometry.
/// @tparam TOutputIterator output iterator with the following value type:
///         @code
///         std::vector<std::pair<
///             std::size_t,        // <== index of the quadratic basis function
///             double              // <== coefficient of the quadratic basis function
///         >>
///         @endcode
/// @param Geometry type of geometry to get the relationships for.
/// @param itOutput output iterator for the DoF relationships.
/// @details each output item consists of the index of a higher order DoF and a related coefficient.
///          The weighted sum of the output items equals the lower order basis function at the index
///          of the output item. Note that higher order DoFs whose coefficient vanish are omitted.
///          More formally, the output expresses the lower order basis functions (\f@ \psi \f@) as
///          a linear combination of the higher order basis functions (\f@ \phi \f@).
///          \f[
///             \begin{align}
///                 \psi_i(\bold x) &= \phi_i(\bold x) + \frac{1}{2} \sum_{j \in S_i}{\phi_j(\bold x)}
///                 S_i &= \{ j | \psi_i(\bold x_j) \neq 0 \}
///             \end{align}
///          \f]
/// @note The basis functions related to DoFs are assumed to be Lagrangian.
template <class TOutputIterator>
void GetNonzeroQuadraticCoefficients(GeometryData::KratosGeometryType Geometry,
                                     TOutputIterator itOutput)
{
    using G = GeometryData::KratosGeometryType;
    using Pair = std::pair<
        std::size_t,    // index of the quadratic basis function
        double          // coefficient of the quadratic basis function
    >;
    using PairVector = std::vector<Pair>;
    switch (Geometry) {
        case G::Kratos_Point2D: {
            *itOutput = PairVector {{0ul, 1.0}};
            break;
        }
        case G::Kratos_Point3D: {
            *itOutput = PairVector {{0ul, 1.0}};
            break;
        }
        case G::Kratos_Line2D2: {
            *itOutput++ = PairVector {{0ul, 1.0}};
            *itOutput = PairVector   {{1ul, 1.0}};
            break;
        }
        case G::Kratos_Line2D3: {
            *itOutput++ = PairVector {{0ul, 1.0},
                                      {2ul, 0.5}};
            *itOutput = PairVector   {{1ul, 1.0},
                                      {2ul, 0.5}};
            break;
        }
        case G::Kratos_Line3D2: {
            *itOutput++ = PairVector {{0ul, 1.0}};
            *itOutput = PairVector   {{1ul, 1.0}};
            break;
        }
        case G::Kratos_Line3D3: {
            *itOutput++ = PairVector {{0ul, 1.0},
                                      {2ul, 0.5}};
            *itOutput = PairVector   {{1ul, 1.0},
                                      {2ul, 0.5}};
            break;
        }
        case G::Kratos_Triangle2D3: {
            *itOutput++ = PairVector {{0ul, 1.0}};
            *itOutput++ = PairVector {{1ul, 1.0}};
            *itOutput = PairVector   {{2ul, 1.0}};
            break;
        }
        case G::Kratos_Triangle2D6: {
            *itOutput++ = PairVector {{0ul, 1.0},
                                      {3ul, 0.5},
                                      {5ul, 0.5}};
            *itOutput++ = PairVector {{1ul, 1.0},
                                      {3ul, 0.5},
                                      {4ul, 0.5}};
            *itOutput = PairVector   {{2ul, 1.0},
                                      {4ul, 0.5},
                                      {5ul, 0.5}};
            break;
        }
        case G::Kratos_Triangle3D6: {
            *itOutput++ = PairVector {{0ul, 1.0},
                                      {3ul, 0.5},
                                      {5ul, 0.5}};
            *itOutput++ = PairVector {{1ul, 1.0},
                                      {3ul, 0.5},
                                      {4ul, 0.5}};
            *itOutput = PairVector   {{2ul, 1.0},
                                      {4ul, 0.5},
                                      {5ul, 0.5}};
            break;
        }
        case G::Kratos_Tetrahedra3D10: {
            *itOutput++ = PairVector {{0ul, 1.0},
                                      {4ul, 0.5},
                                      {6ul, 0.5},
                                      {7ul, 0.5}};
            *itOutput++ = PairVector {{1ul, 1.0},
                                      {4ul, 0.5},
                                      {5ul, 0.5},
                                      {8ul, 0.5}};
            *itOutput++ = PairVector {{2ul, 1.0},
                                      {5ul, 0.5},
                                      {6ul, 0.5},
                                      {9ul, 0.5}};
            *itOutput = PairVector   {{3ul, 1.0},
                                      {7ul, 0.5},
                                      {8ul, 0.5},
                                      {9ul, 0.5}};
            break;
        }
        default: KRATOS_ERROR << "unsupported geometry type " << (int) Geometry << "\n";
    }
}



namespace Detail {

using SlaveMap = std::unordered_map<
    std::size_t,    // <== equation ID of the slave DoF
    std::pair<
        std::size_t,            // <== equation ID of the master DoF
        std::vector<std::tuple<
            const Element*,     // <== element that contains the MPC DoF
            std::size_t,        // <== index of the coarse node within the element that contains the MPC DoF
            std::size_t         // <== index of the MPC DoF within the coarse node
        >>
    >
>;

using MasterMap = std::unordered_map<
    std::size_t,                    // <== equation ID of the master DoF
    std::vector<std::pair<
        SlaveMap::const_iterator,   // <== iterator pointing to the slave's info
        double                      // <== coefficient of the MPC relationship
    >>
>;

} // namespace Detail



template <Globals::DataLocation TEntityLocation>
void MapHigherToLowerOrder(const ModelPart& rModelPart,
                           const Detail::MasterMap& rMasterMap,
                           Detail::SlaveMap& rSlaveMap,
                           const std::vector<bool>& rCoarseDofMask,
                           const std::vector<std::size_t>& rCoarseDofMap,
                           std::map<
                                std::pair<std::size_t,std::size_t>,
                                double,
                                Detail::IndexPairTraits::IsLess
                           >& rRestrictionMap)
{
    KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);
    KRATOS_TRY

    // Assume all nodes have the same number of DoFs
    const std::size_t dofs_per_node = rModelPart.Nodes().empty() ?
        0 : rModelPart.Nodes().front().GetDofs().size();

    // The entries in rSlaveMap are extended, but no entry is removed or added
    // to the map itself, so its sentinel can be safely used throughout the function.
    const auto it_slave_end = rSlaveMap.end();

    for (auto proxy : MakeProxy<TEntityLocation>(rModelPart)) {
        const auto& r_entity = proxy.GetEntity();

        // Get the geometry of the current element, and its
        // associated lower order geometry.
        const Geometry<Node>& r_fine_geometry = r_entity.GetGeometry();
        const std::size_t fine_vertex_count = r_fine_geometry.PointsNumber();

        if (r_entity.IsActive() && fine_vertex_count) {
            std::vector<
                std::vector<std::pair<
                    std::size_t,    // <== fine DoF index
                    double          // <== fine DoF coefficient
                >>
            > restriction_coefficients;
            restriction_coefficients.reserve(4ul); // <== reserve assuming a tetrahedron (linear vertex count)
            GetNonzeroQuadraticCoefficients(r_fine_geometry.GetGeometryType(),
                                            std::back_inserter(restriction_coefficients));

            for (std::size_t i_coarse_vertex=0ul; i_coarse_vertex<restriction_coefficients.size(); ++i_coarse_vertex) {
                const auto& r_restriction_terms = restriction_coefficients[i_coarse_vertex];
                const auto& r_coarse_node_dofs = r_fine_geometry[i_coarse_vertex].GetDofs();

                for (std::size_t i_node_dof=0ul; i_node_dof<dofs_per_node; ++i_node_dof) {
                    const std::size_t i_coarse_dof_in_fine_system = r_coarse_node_dofs[i_node_dof]->EquationId();

                    // If the coarse DoF is constrained by a Dirichlet condition,
                    // no fine vertices should be added to the restiction operator.
                    const std::size_t i_term_end = r_coarse_node_dofs[i_node_dof]->IsFixed() || r_fine_geometry[i_coarse_vertex].Is(SLAVE)
                                                   ? 1 : r_restriction_terms.size();

                    for (std::size_t i_term=0ul; i_term<i_term_end; ++i_term) {
                        const auto& r_fine_coefficient_pair = r_restriction_terms[i_term];
                        const std::size_t i_fine_vertex = r_fine_coefficient_pair.first;
                        const auto& r_fine_node_dofs = r_fine_geometry[i_fine_vertex].GetDofs();
                        KRATOS_ERROR_IF_NOT(r_fine_node_dofs.size() == r_coarse_node_dofs.size());

                        if (rCoarseDofMask[i_coarse_dof_in_fine_system]) { // <== ignore DoF if not in the coarse set
                            const std::size_t i_coarse_dof = rCoarseDofMap[i_coarse_dof_in_fine_system];
                            const std::size_t i_fine_dof = r_fine_node_dofs[i_node_dof]->EquationId();
                            const auto it_slave = rSlaveMap.find(i_fine_dof);
                            const bool is_slave = it_slave != rSlaveMap.end();
                            if (!is_slave) {
                                // The DoF is not a slave, so it can be added to the restriction map
                                rRestrictionMap.emplace(std::make_pair(i_coarse_dof, i_fine_dof), r_fine_coefficient_pair.second);
                            } else {
                                // The DoF is a slave, so its master should be added to the restriction map instead
                                rRestrictionMap.emplace(std::make_pair(i_coarse_dof, it_slave->second.first), r_fine_coefficient_pair.second);
                            }
                        } // if rCoarseDofMask[i_coarse_dof_in_fine_system]

                        // Check whether the coarse DoF participates in an MPC, and record
                        // data about its connectivities if it does.
                        const auto it_slave_dof = rSlaveMap.find(i_coarse_dof_in_fine_system);
                        if (it_slave_dof != it_slave_end) {
                            it_slave_dof->second.second.emplace_back(&r_entity,
                                                                     i_coarse_vertex,
                                                                     i_node_dof);
                        } // if it_slave_dof != rMasterMap.end
                    } // for i_node_dof in range(dofs_per_node)
                } // for coefficient_pair in restriction_terms
            } // for i_coarse_vertex in restriction_coefficients.size()

            // Loop over the rest of the fine vertices to fill the slave map
            if (!rSlaveMap.empty()) {
                for (std::size_t i_fine_vertex=restriction_coefficients.size(); i_fine_vertex<r_fine_geometry.size(); ++i_fine_vertex) {
                    const Node& r_fine_node = r_fine_geometry[i_fine_vertex];
                    const std::size_t dofs_per_node = r_fine_node.GetDofs().size();
                    for (std::size_t i_node_dof=0ul; i_node_dof<dofs_per_node; ++i_node_dof) {
                        const std::size_t i_dof = r_fine_node.GetDofs()[i_node_dof]->EquationId();
                        const auto it_slave_info = rSlaveMap.find(i_dof);
                        if (it_slave_info != it_slave_end) {
                            it_slave_info->second.second.emplace_back(&r_entity,
                                                                      i_fine_vertex,
                                                                      i_node_dof);
                        } // if it_slave_info != rMasterMap.end
                    } // for i_node_dof in range(dofs_per_node)
                } // for i_fine_vertex in range(coarse_vertex_count, fine_vertex_count)
            } // if MasterSlaveDofMap.has_value()

        } // if (r_entity.IsActive())
    } // for entity in rModelPart

    // Loop over master DoFs and register their contributions
    for (const auto& r_master_info : rMasterMap) {
        const std::size_t i_coarse_dof_in_fine_system = r_master_info.first;

        // Skip the DoF if it is not part of the coarse system
        if (!rCoarseDofMask[i_coarse_dof_in_fine_system]) {
            continue;
        }

        const std::size_t i_coarse_dof = rCoarseDofMap[i_coarse_dof_in_fine_system];

        for (auto slave_pair : r_master_info.second) {
            const auto it_slave_info = slave_pair.first;
            const std::size_t i_slave = it_slave_info->first;
            const auto mpc_coefficient = slave_pair.second;

            KRATOS_ERROR_IF(it_slave_info == it_slave_end);
            const auto& r_slave_info = it_slave_info->second;

            for (auto entry : r_slave_info.second) {
                const Geometry<Node>& r_fine_geometry = std::get<0>(entry)->GetGeometry();
                const std::size_t i_vertex = std::get<1>(entry);
                const std::size_t i_node_dof = std::get<2>(entry);

                std::vector<
                    std::vector<std::pair<
                        std::size_t,    // <== fine DoF index
                        double          // <== fine DoF coefficient
                    >>
                > restriction_coefficients;
                restriction_coefficients.reserve(4ul); // <== reserve assuming a tetrahedron (linear vertex count)
                GetNonzeroQuadraticCoefficients(r_fine_geometry.GetGeometryType(),
                                                std::back_inserter(restriction_coefficients));
                if (i_vertex < restriction_coefficients.size()) {
                    KRATOS_ERROR_IF_NOT(i_vertex < restriction_coefficients.size());

                    const auto& r_restriction_terms = restriction_coefficients[i_vertex];

                    // The first term always corresponds to the linear equivalent of the vertex,
                    // which would not have an impact on the restricted stiffness under normal
                    // circumstances, but Kratos takes a dangerous approach to enforcing MPCs.
                    // Instead of removing rows and columns related to slave DoFs from the system,
                    // an arbitrary nonzero value is placed on the main diagonal. If the first term
                    // of slave DoFs is included in the restriction operator, this arbitrary nonzero
                    // value will completely destroy the rows of the coarse system matrix related to
                    // master DoFs.
                    // ==> these first terms must be omitted from the restriction operator.
                    std::size_t i_term = 1ul;

                    for (; i_term<r_restriction_terms.size(); ++i_term) {
                        const auto& r_fine_coefficient_pair = r_restriction_terms[i_term];
                        const std::size_t i_fine_vertex = r_fine_coefficient_pair.first;
                        const std::size_t i_fine_dof = r_fine_geometry[i_fine_vertex].GetDofs()[i_node_dof]->EquationId();
                        const auto it_slave = rSlaveMap.find(i_fine_dof);
                        const bool is_slave = it_slave != it_slave_end;
                        if (!is_slave) {
                            const auto restriction_coefficient = r_fine_coefficient_pair.second * mpc_coefficient;
                            rRestrictionMap.emplace(std::make_pair(i_coarse_dof, i_fine_dof), restriction_coefficient);
                            //rRestrictionMap.emplace(std::make_pair(i_coarse_dof, i_fine_dof), 0.0).first->second += restriction_coefficient;
                        } else {
                            const auto restriction_coefficient = r_fine_coefficient_pair.second * mpc_coefficient;
                            rRestrictionMap.emplace(std::make_pair(i_coarse_dof, it_slave->second.first), restriction_coefficient);
                        }
                    }
                } /* if (i_vertex < restriction_coefficients.size()) */ else {
                    const auto restriction_coefficient = 0.5 * mpc_coefficient;
                    rRestrictionMap.emplace(std::make_pair(i_coarse_dof, i_slave), restriction_coefficient);
                }
            } // for entry in r_slave_info.mContainingElements
        } // for slave_pair : r_master_info.second.mSlaves
    } // for r_master_info in rMasterSlaveDofMap.value()
    KRATOS_CATCH("")
}



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
void HierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::ProvideAdditionalData(SparseMatrix& rA,
                                                                                    Vector& rX,
                                                                                    Vector& rB,
                                                                                    ModelPart::DofsArrayType& rDofs,
                                                                                    ModelPart& rModelPart)
{
    KRATOS_TRY
    KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);
    const auto timer = BuiltinTimer();

    KRATOS_INFO_IF("HierarchicalSolver", 2 <= mpImpl->mVerbosity)
        << "construct coarse system\n";

    if (4 <= mpImpl->mVerbosity) {
        KRATOS_INFO("HierarchicalSolver") << "write system_matrix.mm\n";
        TSparseSpace::WriteMatrixMarketMatrix("system_matrix.mm", rA, false);

        KRATOS_INFO("HierarchicalSolver") << "write rhs.mm\n";
        TSparseSpace::WriteMatrixMarketVector("rhs.mm", rB);
    }

    const std::size_t system_size = TSparseSpace::Size1(rA);
    std::vector<bool> coarse_mask(system_size, false); // <== yes, an intentional vector<bool>
    GetLowerOrderDofs<Globals::DataLocation::Element,TSparseSpace>(rA,
                                                                   rModelPart,
                                                                   rDofs,
                                                                   coarse_mask);

    // Conditions shouldn't introduce coarse DoFs, but I'll leave this
    // here for the time being in case I'm missing some exocotic boundary
    // condition implementations @matekelemen.
    //GetLowerOrderDofs<Globals::DataLocation::Condition,TSparseSpace>(rA,
    //                                                                 rModelPart,
    //                                                                 rDofs,
    //                                                                 coarse_mask);

    // Compute an index map associating fine DoF indices with coarse ones
    // and collect coarse DoFs.
    // Indices at fine DoFs should not be accessed
    // (doing so will result in a segfault instead of silent miscalculations)
    std::vector<std::size_t> coarse_dof_indices(system_size,
                                                std::numeric_limits<std::size_t>::max());

    // If the coarse solver's AdditionalPhysicalDataIsNeeded is true, the coarse
    // DoF list must be constructed, and it must map equation IDs to the coarse
    // system's index space. This means constructing an entirely new list of DoFs,
    // so it should be avoided if the coarse solver doesn't need it
    // (hence the std::optional).
    std::optional<std::vector<ModelPart::DofType>> coarse_dofs;
    std::size_t masked_count = 0ul;
    if (mpImpl->mCoarseSystem.mpSolver->AdditionalPhysicalDataIsNeeded()) {
        // Reserve space for coarse DoFs assuming the system consists exclusively
        // of quadratic tetrahedral elements (10 nodes per element) that get reduced
        // to linear tetrahedra (4 nodes per element)
        coarse_dofs.emplace();
        auto& r_coarse_dofs = coarse_dofs.value();
        r_coarse_dofs.reserve(std::size_t(0.4 * system_size));

        // Fill the fine=>coarse DoF index map and set coarse DoFs
        for (std::size_t i_mask=0; i_mask<system_size; ++i_mask) {
            if (coarse_mask[i_mask]) {
                r_coarse_dofs.push_back(**(rDofs.ptr_begin() + i_mask));
                r_coarse_dofs.back().SetEquationId(masked_count);
                coarse_dof_indices[i_mask] = masked_count++;
            }
        }
    } /* if CoarseSolver.AdditionalPhysicalDataIsNeeded */ else {
        // Fill the fine => coarse DoF index map
        for (std::size_t i_mask=0; i_mask<system_size; ++i_mask) {
            if (coarse_mask[i_mask]) {
                coarse_dof_indices[i_mask] = masked_count++;
            }
        }
    }

    KRATOS_INFO_IF("HierarchicalSolver", 1 <= mpImpl->mVerbosity)
        << "coarse system size: " << masked_count << "\n";

    // Compute the restriction operator and the interpolation operator (transpose of the restriction)
    typename TSparseSpace::MatrixType& r_restriction_operator = mpImpl->mRestrictionOperator;
    typename TSparseSpace::MatrixType& r_interpolation_operator = mpImpl->mInterpolationOperator;

    {
        // The restriction and interpolation operators are linear, so they can be stored in
        // matrix format. In this case they are stored as sparse matrices, which complicates
        // their construction somewhat. Entries are calculated by looping over elements and
        // their nodes, which does not guarantee that DoFs are iterated in sorted order of
        // their equation IDs (row indices of the restriction operator), and each DoF may
        // be iterated over more than once. This poses a problem for the efficiency of
        // filling the sparse matrices, so the construction is broken down into two stages:
        // 1) collect the sparse matrix entries in a set, sorted in row-major order
        // 2) insert entries to the sparse matrix from the set
        std::map<
            std::pair<std::size_t,std::size_t>,
            double,
            Detail::IndexPairTraits::IsLess
        > restriction_map;

        // Construct a map relating slave DoFs to their masters, which effectively couples the
        // corresponding basis functions' support.
        // @todo properly explain what these maps represent and why they're necessary @matekelemen
        Detail::MasterMap master_map;
        Detail::SlaveMap slave_map;

        // First loop through the constraints:
        // ==> collect slaves
        for (const auto& r_constraint : rModelPart.MasterSlaveConstraints()) {
            for (const auto* p_slave : r_constraint.GetSlaveDofsVector()) {
                slave_map.emplace(p_slave->EquationId(), Detail::SlaveMap::mapped_type {r_constraint.GetMasterDofsVector()[0]->EquationId(), {}});
            } // for p_slave in r_constraint.GetSlaveDofsVector
        } // for r_constraint in rModelPart.MasterSlaveConstraints

        // Second loop through the constraints:
        // ==> collect slaves and coefficients for each master
        for (const auto& r_constraint : rModelPart.MasterSlaveConstraints()) {
            MasterSlaveConstraint::MatrixType coefficients;
            MasterSlaveConstraint::VectorType constants;
            r_constraint.GetLocalSystem(coefficients, constants, rModelPart.GetProcessInfo());

            for (std::size_t i_master=0ul; i_master<r_constraint.GetMasterDofsVector().size(); ++i_master) {
                const auto* p_master = r_constraint.GetMasterDofsVector()[i_master];
                auto it_master = master_map.emplace(p_master->EquationId(), Detail::MasterMap::mapped_type {}).first;

                for (std::size_t i_slave=0ul; i_slave<r_constraint.GetSlaveDofsVector().size(); ++i_slave) {
                    const auto coefficient = coefficients(i_slave, i_master);
                    if (coefficient) {
                        const auto* p_slave = r_constraint.GetSlaveDofsVector()[i_slave];
                        const auto it_slave = slave_map.find(p_slave->EquationId());
                        KRATOS_ERROR_IF(it_slave == slave_map.end());
                        it_master->second.emplace_back(it_slave, coefficients(i_slave, i_master));
                    } // if coefficient != 0
                } // for p_slave_dof in r_constraint.GetSlaveDof
            } // for p_master_dof in r_constraint.GetMasterDofsVector
        } // for r_constraint in rModelPart.MasterSlaveConstraints

        // Collect entries of the restriction operator's sparse matrix
        // from elements
        MapHigherToLowerOrder<Globals::DataLocation::Element>(rModelPart,
                                                              master_map,
                                                              slave_map,
                                                              coarse_mask,
                                                              coarse_dof_indices,
                                                              restriction_map);

        // Conditions shouldn't introduce coarse DoFs, but I'll leave this
        // here for the time being in case I'm missing some exotic boundary
        // condition implementations @matekelemen.
        // Collect entries of the restriction operator's sparse matrix
        // from conditions
        //MapHigherToLowerOrder<Globals::DataLocation::Condition>(rModelPart,
        //                                                        master_slave_dof_map,
        //                                                        coarse_dof_indices,
        //                                                        restriction_map);

        {
            KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);

            // Construct the restriction matrix from the collected entries
            TSparseSpace::Resize(r_restriction_operator,
                                 masked_count,
                                 system_size);
            for (auto pair : restriction_map) {
                const std::size_t i_row = pair.first.first;
                const std::size_t i_column = pair.first.second;
                KRATOS_ERROR_IF_NOT(i_row < masked_count) << i_row << " is out of range " << masked_count << "\n";
                KRATOS_ERROR_IF_NOT(i_column < system_size) << i_column << " is out of range " << system_size << "\n";
                const double value = pair.second;
                r_restriction_operator.insert_element(i_row, i_column, value);
            }
            r_restriction_operator.complete_index1_data();

            // The (approximate) interpolation operator is the transpose of the restriction operator
            SparseMatrixMultiplicationUtility::TransposeMatrix(r_interpolation_operator, r_restriction_operator, 1.0);
        } // end profiling
    } // destroy restriction_map

    if (4 <= mpImpl->mVerbosity) {
        KRATOS_INFO("HierarchicalSolver") << "write restriction_operator.mm\n";
        TSparseSpace::WriteMatrixMarketMatrix("restriction_operator.mm", mpImpl->mRestrictionOperator, false);

        KRATOS_INFO("HierarchicalSolver") << "write interpolation_operator.mm\n";
        TSparseSpace::WriteMatrixMarketMatrix("interpolation_operator.mm", mpImpl->mInterpolationOperator, false);
    }

    // Compute coarse system components
    // coarse_system = restriction_operator * fine_system * restriction_operator^T
    {
        KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);
        typename TSparseSpace::MatrixType tmp;
        TSparseSpace::Resize(tmp, masked_count, system_size);
        TSparseSpace::Resize(mpImpl->mCoarseSystem.mA, masked_count, masked_count);

        SparseMatrixMultiplicationUtility::MatrixMultiplication(mpImpl->mRestrictionOperator,
                                                                rA,
                                                                tmp);
        SparseMatrixMultiplicationUtility::MatrixMultiplication(tmp,
                                                                mpImpl->mInterpolationOperator,
                                                                mpImpl->mCoarseSystem.mA);

        if (4 <= mpImpl->mVerbosity) {
            KRATOS_INFO("HierarchicalSolver") << "write coarse_system_matrix.mm\n";
            TSparseSpace::WriteMatrixMarketMatrix("coarse_system_matrix.mm", mpImpl->mCoarseSystem.mA, false);

            KRATOS_INFO("HierarchicalSolver") << "write coarse_rhs.mm\n";
            typename TSparseSpace::VectorType coarse_rhs(masked_count);
            TSparseSpace::Mult(mpImpl->mRestrictionOperator, rB, coarse_rhs);
            TSparseSpace::WriteMatrixMarketVector("coarse_rhs.mm", coarse_rhs);
        }
    } // end profiling


    if (mpImpl->mCoarseSystem.mpSolver->AdditionalPhysicalDataIsNeeded()) {
        // Create a vector of pointers to coarse system DoFs
        KRATOS_ERROR_IF_NOT(coarse_dofs.has_value());
        auto& r_coarse_dofs = coarse_dofs.value();
        ModelPart::DofsArrayType coarse_dof_pointers;

        // Cannot use standard algorithms on PointerVectorSet because
        // its value and iterator types are ill-formed
        coarse_dof_pointers.reserve(r_coarse_dofs.size());
        for (auto& rDof : r_coarse_dofs) {
            coarse_dof_pointers.push_back(&rDof);
        }
        coarse_dof_pointers.Sort();

        // Restrict solution and RHS vectors
        Vector coarse_x(masked_count), coarse_b(masked_count);
        TSparseSpace::Mult(mpImpl->mRestrictionOperator,
                           rX,
                           coarse_x);
        TSparseSpace::Mult(mpImpl->mRestrictionOperator,
                           rB,
                           coarse_b);

        // Provide additional physical data for the coarse system
        mpImpl->mCoarseSystem.mpSolver->ProvideAdditionalData(mpImpl->mCoarseSystem.mA,
                                                              coarse_x,
                                                              coarse_b,
                                                              coarse_dof_pointers,
                                                              rModelPart);
    } // destroy coarse system containers (excluding the system matrix)

    if (mpImpl->mpFineSolver->AdditionalPhysicalDataIsNeeded()) {
        // Provide additional physical data for the fine system
        mpImpl->mpFineSolver->ProvideAdditionalData(rA,
                                                    rX,
                                                    rB,
                                                    rDofs,
                                                    rModelPart);
    }

    KRATOS_INFO_IF("HierarchicalSolver", 2 <= mpImpl->mVerbosity)
        << "constructing coarse system took " << timer << "\n";
    KRATOS_CATCH("")
}



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
Parameters
HierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::GetDefaultParameters()
{
    return Parameters(R"(
{
    "solver_type" : "hierarchical",
    "verbosity" : 0,
    "max_iterations" : 50,
    "tolerance" : 1e-6,
    "coarse_solver_settings" : {
        "solver_type" : "cg",
        "max_iteration" : 50,
        "tolerance" : 2e-1
    },
    "fine_solver_settings" : {
        "solver_type" : "gauss_seidel",
        "max_iterations" : 3
    }
}
    )");
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
bool HierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::Solve(SparseMatrix& rA,
                                                                    DenseMatrix& rX,
                                                                    DenseMatrix& rB)
{
    return false;
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
void HierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::PrintInfo(std::ostream& rStream) const
{
    rStream << "HierarchicalSolver";
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
void HierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::PrintData(std::ostream& rStream) const
{
    rStream
        << "tolerance     : " << mpImpl->mTolerance << "\n"
        << "max iterations: " << mpImpl->mMaxIterations << "\n"
        << "verbosity     : " << mpImpl->mVerbosity << "\n"
        ;
}



template
class HierarchicalSolver<
    TUblasSparseSpace<double>,
    TUblasDenseSpace<double>,
    Reorderer<
        TUblasSparseSpace<double>,
        TUblasDenseSpace<double>
    >
>;


} // namespace Kratos
