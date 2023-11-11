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

#ifndef AMGCL_PARAM_UNKNOWN
    #include "input_output/logger.h"
    #define AMGCL_PARAM_UNKNOWN(NAME)                       \
        KRATOS_ERROR                                        \
            << KRATOS_CODE_LOCATION                         \
            << "Unknown parameter " << NAME << std::endl
#endif

// Project includes
#include "amgcl_hierarchical_solver.h"
#include "spaces/ublas_space.h"
#include "utilities/profiler.h"

// External includes
#include "amgcl/adapter/ublas.hpp"
#include "amgcl/adapter/zero_copy.hpp"
#include "amgcl/backend/builtin.hpp"
#include "amgcl/value_type/static_matrix.hpp"
#include "amgcl/make_solver.hpp"
#include "amgcl/make_block_solver.hpp"
#include "amgcl/solver/runtime.hpp"
#include "amgcl/preconditioner/runtime.hpp"
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"

// STL includes
#include <type_traits> // remove_reference_t, remove_cv_t
#include <sstream> // stringstream
#include <optional> // optional
#include <bitset> // bitset


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


/// Composite preconditioner traits
struct CPTraits
{
    using Mask = std::vector<char>;

    struct AMGCL
    {
        using ScalarBackend = amgcl::backend::builtin<double>;
    }; // struct amgcl
}; // struct CPTraits


/** @brief Generate four submatrices from a square matrix.
 *
 *  @details Let @f(A@f) be the square input matrix, and
 *           @f(A_{ll}@f), @f(A_{lq}@f), @f(A_{ql}@f), @f(A_{qq}@f) the
 *           output matrices in this exact order. Then this function resizes
 *           and populates the output matrices such that
 *           @f[
 *              A = \begin{bmatrix}
 *                  A_{ll} & A_{lq} \\
 *                  A_{ql} & A_{qq}
 *              \end{bmatrix}
 *           @f]
 *
 *          @f(A_{ll}@f) is defined by a boolean mask, represented by a @c char
 *          array consisting of 1s at row indices where @f(A_{ll}@f) has
 *          components, and 0s everywhere else.
 *
 *  @throws if the mask @a rMask contains items other than 0 or 1.
 */
void MakeSubblocks(const TUblasSparseSpace<double>::MatrixType& rRootMatrix,
                   const CPTraits::Mask& rMask,
                   const std::array<TUblasSparseSpace<double>::MatrixType*,4>& rOutput)
{
    KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);
    KRATOS_ERROR_IF_NOT(rMask.size() == rRootMatrix.size1())
        << "Mask size mismatch: mask of size " << rMask.size()
        << " provided for matrix of size " << rRootMatrix.size1() << "x" << rRootMatrix.size2();

    KRATOS_ERROR_IF_NOT(rRootMatrix.size1() == rRootMatrix.size2())
        << "Expecting a square matrix, but got " << rRootMatrix.size1() << "x" << rRootMatrix.size2();

    KRATOS_ERROR_IF_NOT(std::count(rOutput.begin(), rOutput.end(), nullptr) == 0)
        << "Output pointers must point to existing matrices";

    const std::size_t system_size = rMask.size();

    // Compute the index of each component local to the submatrix it belongs to.
    std::vector<std::size_t> block_local_indices(system_size);
    std::size_t masked_count = 0;
    std::size_t unmasked_count = 0;
    for (std::size_t i_mask=0; i_mask<system_size; ++i_mask) {
        const auto mask_value = rMask[i_mask];
        if (mask_value == 1) {
            block_local_indices[i_mask] = masked_count++;
        } else if (mask_value == 0) {
            block_local_indices[i_mask] = unmasked_count++;
        } else {
            KRATOS_ERROR << "Invalid mask value '" << static_cast<int>(mask_value)
                         << "' at index " << i_mask << ". The mask must only consist of "
                         << " 0s or 1s.";
        }
    }

    // Resize submatrices.
    KRATOS_TRY
    rOutput[0]->resize(  masked_count,   masked_count, false);
    rOutput[1]->resize(  masked_count, unmasked_count, false);
    rOutput[2]->resize(unmasked_count,   masked_count, false);
    rOutput[3]->resize(unmasked_count, unmasked_count, false);
    KRATOS_CATCH("")

    // Fill submatrices
    // Note: the matrices are filled in proper row-major order,
    //       so each insertion should be O(1) (unless a reallocation
    //       is triggered).
    KRATOS_TRY
    for (auto it_row=rRootMatrix.begin1(); it_row!=rRootMatrix.end1(); ++it_row) {
        const std::size_t i_row = it_row.index1();
        const char row_mask = rMask[i_row];
        for (auto it_entry=it_row.begin(); it_entry!=it_row.end(); ++it_entry) {
            const std::size_t i_column = it_entry.index2();
            const char column_mask = rMask[i_column];
            const char i_block = (~(column_mask | (row_mask << 1))) & 0b11; // <== flattened index of the block
            rOutput[i_block]->insert_element(block_local_indices[i_row],
                                             block_local_indices[i_column],
                                             *it_entry);
        }
    }
    KRATOS_CATCH("")
}


} // namespace detail



template <class TSparseSpace,
          class TDenseSpace,
          class TReorderer>
struct AMGCLHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::Impl
{
    double mTolerance;

    int mVerbosity;

    std::size_t mDoFCount;

    std::size_t mMaxIterations;

    // Indicates whether the variable at the given index is lower order.
    std::optional<std::vector<char>> mLowerDofMask;

    boost::property_tree::ptree mLowerAMGCLSettings;

    boost::property_tree::ptree mFullAMGCLSettings;
}; // struct AMGCLHierarchicalSolver::Impl



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
AMGCLHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::AMGCLHierarchicalSolver(Parameters parameters)
    : mpImpl(new Impl)
{
    KRATOS_TRY
    Parameters default_parameters = this->GetDefaultParameters();
    parameters.ValidateAndAssignDefaults(default_parameters);

    KRATOS_ERROR_IF_NOT(parameters["solver_type"].GetString() == "amgcl_hierarchical")
        << "Requested a(n) '" << parameters["solver_type"].GetString() << "' solver,"
        << " but constructing an AMGCLHierarchicalSolver";

    detail::ParametersWithDefaults safe_parameters(parameters, default_parameters);
    mpImpl->mTolerance = safe_parameters["amgcl_settings"]["coarse"]["solver"]["tol"].Get<double>();
    mpImpl->mVerbosity = safe_parameters["verbosity"].Get<int>();
    mpImpl->mMaxIterations = safe_parameters["max_iterations"].Get<int>();
    mpImpl->mDoFCount = 0ul;

    // Convert parameters to AMGCL settings
    std::stringstream json_stream;

    json_stream << safe_parameters["amgcl_settings"]["coarse"].GetInput().PrettyPrintJsonString();
    boost::property_tree::read_json(
        json_stream,
        mpImpl->mLowerAMGCLSettings
    );

    json_stream.clear();
    json_stream << safe_parameters["amgcl_settings"]["full"].GetInput().PrettyPrintJsonString();
    boost::property_tree::read_json(
        json_stream,
        mpImpl->mFullAMGCLSettings
    );
    KRATOS_CATCH("")
}



// Necessary for PIMPL
template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
AMGCLHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::~AMGCLHierarchicalSolver()
{
}



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
bool AMGCLHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::Solve(SparseMatrix& rA,
                                                                         Vector& rX,
                                                                         Vector& rB)
{
    KRATOS_TRY
    KRATOS_ERROR_IF_NOT(mpImpl->mLowerDofMask.has_value())
        << "AMGCLHierarchicalSolver::Solve called before AMGCLHierarchicalSolver::ProvideAdditionalData";
    auto& r_lower_dof_mask = mpImpl->mLowerDofMask.value();
    KRATOS_ERROR_IF_NOT(r_lower_dof_mask.size() == TSparseSpace::Size1(rA))
        << "DoF mask size " << r_lower_dof_mask.size()
        << " does not match system size " << TSparseSpace::Size1(rA);

    mpImpl->mLowerAMGCLSettings.put("solver.verbose", 1 < mpImpl->mVerbosity);

    std::tuple<std::size_t,double> solver_results {
        std::numeric_limits<std::size_t>::max(),    // iteration count
        std::numeric_limits<double>::max()          // residual
    };

    switch (mpImpl->mDoFCount) {
        case 3:
            solver_results = this->SolveImpl<2>(rA, rX, rB);
            break;
        case 4:
            solver_results = this->SolveImpl<3>(rA, rX, rB);
            break;
        default:
            solver_results = this->SolveImpl<1>(rA, rX, rB);
    }

    const auto [iteration_count, residual] = solver_results;

    KRATOS_WARNING_IF("AMGCLHierarchicalSolver", mpImpl->mTolerance <= residual)
        << "Failed to converge. Residual: " << residual << "\n";

    if(1 < mpImpl->mVerbosity) {
        std::cout << "Iterations: " << iteration_count << "\n"
                  << "Error: " << residual << "\n"
                  << "\n";
    }

    return residual < mpImpl->mTolerance ? true : false;
    KRATOS_CATCH("")
}



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
void AMGCLHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::ProvideAdditionalData(SparseMatrix& rA,
                                                                                         Vector& rX,
                                                                                         Vector& rB,
                                                                                         ModelPart::DofsArrayType& rDofs,
                                                                                         ModelPart& rModelPart)
{
    KRATOS_TRY
    KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);

    std::optional<std::size_t> old_dof_count;
    std::size_t dof_count = 0;

    // Compute block size
    if (!rModelPart.IsDistributed()) {
        std::size_t old_dof_id = rDofs.empty() ? 0ul : rDofs.begin()->Id();
        for (const auto& rDof : rDofs) {
            if (rDof.EquationId() < TSparseSpace::Size1(rA)) {
                const auto dof_id = rDof.Id();
                if(dof_id != old_dof_id) {
                    old_dof_id = dof_id;
                    if (!old_dof_count.has_value()) {
                        old_dof_count = dof_count;
                    } else if (old_dof_count != dof_count) { //if it is different than the block size is 1
                        old_dof_count = -1;
                        break;
                    }
                    dof_count = 1;
                } else { // else (dof_id != old_dof_id)
                    ++dof_count;
                }
            }
        }
        mpImpl->mDoFCount = old_dof_count.has_value() ? 1 : dof_count;
    } else { //distributed
        const std::size_t system_size = TSparseSpace::Size1(rA);
        int current_rank = rModelPart.GetCommunicator().GetDataCommunicator().Rank();
        unsigned int old_node_id = rDofs.size() ? rDofs.begin()->Id() : 0;
        for (auto it = rDofs.begin(); it!=rDofs.end(); it++) {
            if(it->EquationId() < system_size  && it->GetSolutionStepValue(PARTITION_INDEX) == current_rank) {
                IndexType id = it->Id();
                if(id != old_node_id) {
                    old_node_id = id;
                    if(!old_dof_count.has_value()) {
                        old_dof_count = dof_count;
                    } else if (old_dof_count.value() != dof_count) { //if it is different than the block size is 1
                        old_dof_count = -1;
                        break;
                    }

                    dof_count=1;
                } else {
                    dof_count++;
                }
            }
        }

        if(old_dof_count.has_value()) {
            mpImpl->mDoFCount = dof_count;
        }

        const std::size_t max_block_size = rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(mpImpl->mDoFCount);

        if(!old_dof_count.has_value()) {
            mpImpl->mDoFCount = max_block_size;
        }

        KRATOS_ERROR_IF(mpImpl->mDoFCount != max_block_size) << "Block size is not consistent. Local: " << mpImpl->mDoFCount  << " Max: " << max_block_size << std::endl;
    }

    KRATOS_INFO_IF("AMGCLHierarchicalSolver", 1 < mpImpl->mVerbosity)
        << "Number of DoFs: " << mpImpl->mDoFCount << "\n";

    // Construct mask
    if (!mpImpl->mLowerDofMask.has_value()) {
        mpImpl->mLowerDofMask.emplace();
    }
    auto& r_lower_dof_mask = mpImpl->mLowerDofMask.value();
    r_lower_dof_mask.resize(TSparseSpace::Size1(rA), false);

    // Construct lower order DoF mask
    // Loop through elements and collect the DoFs' IDs
    // that are related to corner vertices.
    Element::DofsVectorType element_dofs;
    for (const Element& r_element : rModelPart.Elements()) {
        const auto& r_geometry = r_element.GetGeometry();
        element_dofs.clear();
        r_element.GetDofList(element_dofs, rModelPart.GetProcessInfo());
        switch (r_geometry.GetGeometryType()) {
            case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10: {
                for (unsigned i_dof=0; i_dof<4; ++i_dof) {
                    const std::size_t equation_id = element_dofs[i_dof]->EquationId();
                    r_lower_dof_mask[equation_id] = true;
                }
                break;
            }
            default: {
                KRATOS_ERROR << "Unsupported element geometry: " << r_geometry;
            }
        } // switch GeometryType
    }
    KRATOS_CATCH("")
}



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
Parameters
AMGCLHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::GetDefaultParameters()
{
    return Parameters(R"(
{
    "solver_type" : "amgcl_hierarchical",
    "verbosity" : 1,
    "max_iterations" : 50,
    "amgcl_settings" : {
        "coarse" : {
            "precond" : {
                "class" : "amg",
                "relax" : {
                    "type" : "ilu0"
                },
                "coarsening" : {
                    "type" : "aggregation",
                    "aggr" : {
                        "eps_strong" : 0.8,
                        "block_size" : 1
                    }
                },
                "coarse_enough" : 333,
                "npre" : 1,
                "npost" : 1
            },
            "solver" : {
                "type" : "cg",
                "maxiter" : 5,
                "tol" : 1e-6
            }
        },
        "full" : {
            "precond" : {
                "class" : "relaxation",
                "type" : "gauss_seidel"
            },
            "solver" : {
                "type" : "preonly"
            }
        }
    }
}
    )");
}



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
template <unsigned BlockSize>
std::tuple<std::size_t,double>
AMGCLHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::SolveImpl(SparseMatrix& rA,
                                                                        Vector& rX,
                                                                        Vector& rB) const
{
    KRATOS_TRY

    // Separate components related to lower and higher order DoFs.
    SparseMatrix A_ll, A_lh, A_hl, A_hh;
    detail::MakeSubblocks(rA, mpImpl->mLowerDofMask.value(), {&A_ll, &A_lh, &A_hl, &A_hh});

    // Dump the generated matrices to disk f requested.
    if (4 <= mpImpl->mVerbosity) {
        KRATOS_INFO("AMGCLHierarchicalSolver")
            << "writing system matrices A.mm b.mm\n";
        TSparseSpace::WriteMatrixMarketMatrix("A.mm", rA, false);
        TSparseSpace::WriteMatrixMarketVector("b.mm", rB);
        KRATOS_INFO("AMGCLHierarchicalSolver")
            << "writing submatrices: A_ll.mm A_lh.mm A_hl.mm A_hh.mm\n";
        TSparseSpace::WriteMatrixMarketMatrix("A_ll.mm", A_ll, false);
        TSparseSpace::WriteMatrixMarketMatrix("A_lh.mm", A_lh, false);
        TSparseSpace::WriteMatrixMarketMatrix("A_hl.mm", A_hl, false);
        TSparseSpace::WriteMatrixMarketMatrix("A_hh.mm", A_hh, false);
        KRATOS_ERROR << "verbosity >=4 prints system matrices and terminates\n";
    }

    // Construct adapters for AMGCL
    auto p_a = amgcl::adapter::zero_copy(TSparseSpace::Size1(rA),
                                         rA.index1_data().begin(),
                                         rA.index2_data().begin(),
                                         rA.value_data().begin());

    auto p_a_ll = amgcl::adapter::zero_copy(TSparseSpace::Size1(A_ll),
                                            A_ll.index1_data().begin(),
                                            A_ll.index2_data().begin(),
                                            A_ll.value_data().begin());

    // Define solver type
    using ScalarBackend = amgcl::backend::builtin<double>;
    using ScalarSolver = amgcl::make_solver<
        amgcl::runtime::preconditioner<ScalarBackend>,
        amgcl::runtime::solver::wrapper<ScalarBackend>
    >;

    using Block = amgcl::static_matrix<double,BlockSize,BlockSize>;
    using BlockBackend = amgcl::backend::builtin<Block>;
    using BlockSolver = amgcl::make_block_solver<
        amgcl::runtime::preconditioner<BlockBackend>,
        amgcl::runtime::solver::wrapper<BlockBackend>
    >;

    using Solver = std::conditional_t<
        BlockSize == 1,
        ScalarSolver,
        BlockSolver
    >;

    // Construct solvers
    auto lower_solver = Solver(*p_a_ll, mpImpl->mLowerAMGCLSettings);
    auto full_solver = ScalarSolver(*p_a, mpImpl->mFullAMGCLSettings);
    KRATOS_INFO_IF("AMGCLHierarchicalSolver", 1 < mpImpl->mVerbosity)
        << "Solver memory usage: " << amgcl::human_readable_memory(amgcl::backend::bytes(lower_solver)
                                                                   + amgcl::backend::bytes(full_solver))
        << "\n";

    // Construct block maps for vectors
    const auto& r_mask = mpImpl->mLowerDofMask.value();
    const std::size_t full_size = r_mask.size();

    const std::size_t unmasked_count = std::count(r_mask.begin(),
                                                  r_mask.end(),
                                                  0);
    const std::size_t masked_count = full_size - unmasked_count;

    // Define a map for restriction, that collects the indices
    // of components local to the block they belong to.
    std::vector<std::size_t> block_local_indices(full_size);

    // Fill the block map
    std::array<std::size_t,2> block_counters {0ul, 0ul}; // {l, h}
    for (std::size_t i_global=0ul; i_global<full_size; ++i_global) {
        const char mask_value = r_mask[i_global];
        const std::size_t i_block = !mask_value;
        block_local_indices[i_global] = block_counters[i_block]++;
    }

    const auto coarsen = [&block_local_indices, &r_mask](const Vector& rRoot,
                                                         std::array<Vector,2>& rSubs) {
        for (std::size_t i_root=0ul; i_root<rRoot.size(); ++i_root) {
            const char mask_value = r_mask[i_root];
            const std::size_t i_block = !mask_value;
            const std::size_t i_sub = block_local_indices[i_root];
            rSubs[i_block][i_sub] = rRoot[i_root];
        }
    };

    const auto aggregate = [&block_local_indices, &r_mask](const std::array<Vector,2>& rSubs,
                                                           Vector& rRoot) {
        for (std::size_t i_root=0ul; i_root<rRoot.size(); ++i_root) {
            const char mask_value = r_mask[i_root];
            const std::size_t i_block = !mask_value;
            const std::size_t i_sub = block_local_indices[i_root];
            rRoot[i_root] += rSubs[i_block][i_sub];
        }
    };

    // Solve
    std::fill(rX.begin(),
              rX.end(),
              0.0);
    Vector residual = rB;
    const double rhs_norm = norm_2(rB);
    KRATOS_ERROR_IF_NOT(rhs_norm);

    std::array<Vector,2> x, r, d; // <== {x_l, x_h}, {r_l, r_h}, {d_l, d_h}
    x[0].resize(  masked_count);
    x[1].resize(unmasked_count);
    r[0].resize(  masked_count);
    r[1].resize(unmasked_count);
    d[0].resize(  masked_count);
    d[1].resize(unmasked_count);

    Vector delta(full_size, 0.0);
    Vector lower_residual(masked_count, 0.0);
    Vector lower_delta(masked_count, 0.0);

    std::tuple<std::size_t,double> return_pair {0ul, 1.0};
    const std::size_t max_iterations = mpImpl->mMaxIterations;
    for (std::size_t i_iteration=0ul; i_iteration<max_iterations; ++i_iteration) {
        // Perform a few relaxation passes on the full system
        KRATOS_TRY
        full_solver(*p_a, residual, delta);
        KRATOS_CATCH("")
        coarsen(delta, d);

        // Compute restricted residual
        coarsen(residual, r);
        lower_residual = r[0] - prod(A_ll, d[0]) - prod(A_lh, d[1]);

        // Solve lower
        KRATOS_TRY
        lower_solver(*p_a_ll, lower_residual, lower_delta);
        KRATOS_CATCH("")
        d[0] += lower_delta;

        // Update solution
        aggregate(d, rX);

        // Compute the residual
        noalias(residual) = rB - prod(rA, rX);

        std::get<0>(return_pair)++;
        std::get<1>(return_pair) = norm_2(residual) / rhs_norm;
        KRATOS_INFO_IF("AMGCLHierarchicalSolver", 1 <= mpImpl->mVerbosity)
            << "iteration " << std::get<0>(return_pair)
            << " residual " << std::get<1>(return_pair) << "\n";
    }
    return return_pair;
    KRATOS_CATCH("")
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
bool AMGCLHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::Solve(SparseMatrix& rA,
                                                                         DenseMatrix& rX,
                                                                         DenseMatrix& rB)
{
    return false;
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
void AMGCLHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::PrintInfo(std::ostream& rStream) const
{
    rStream << "AMGCLHierarchicalSolver";
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
void AMGCLHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::PrintData(std::ostream& rStream) const
{
    rStream
        << (mpImpl->mLowerDofMask.has_value() ? "initialized" : "uninitialized") << "\n"
        << "tolerance     : " << mpImpl->mTolerance << "\n"
        << "DoF size      : " << mpImpl->mDoFCount << "\n"
        << "verbosity     : " << mpImpl->mVerbosity << "\n"
        << "AMGCL settings: "
        ;
    boost::property_tree::json_parser::write_json(rStream, mpImpl->mLowerAMGCLSettings);
}



template
class AMGCLHierarchicalSolver<
    TUblasSparseSpace<double>,
    TUblasDenseSpace<double>,
    Reorderer<
        TUblasSparseSpace<double>,
        TUblasDenseSpace<double>
    >
>;


} // namespace Kratos
