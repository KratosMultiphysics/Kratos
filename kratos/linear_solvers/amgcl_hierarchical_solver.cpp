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

// External includes
#include "amgcl/adapter/crs_tuple.hpp"
#include "amgcl/adapter/ublas.hpp"
#include "amgcl/adapter/zero_copy.hpp"
#include "amgcl/backend/builtin.hpp"
#include "amgcl/value_type/static_matrix.hpp"
#include "amgcl/make_solver.hpp"
#include "amgcl/make_block_solver.hpp"
#include "amgcl/solver/runtime.hpp"
#include "amgcl/preconditioner/schur_pressure_correction.hpp"
#include "amgcl/preconditioner/runtime.hpp"
#include "amgcl/io/mm.hpp"
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"

// Project includes
#include "amgcl_hierarchical_solver.h"
#include "spaces/ublas_space.h"
#include "utilities/profiler.h"

// STL includes
#include <amgcl/util.hpp>
#include <type_traits> // remove_reference_t, remove_cv_t
#include <sstream> // stringstream
#include <optional> // optional
#include <bitset> // bitset


namespace Kratos {


namespace detail {


template <unsigned BlockSize>
struct CompositePreconditionedSolver
{
private:
    static_assert(0 < BlockSize, "Invalid block size");

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

public:
    using Type = amgcl::make_solver<
        amgcl::preconditioner::schur_pressure_correction<BlockSolver,ScalarSolver>,
        amgcl::runtime::solver::wrapper<ScalarBackend>
    >;
}; // struct CompositePreconditionedSolver



template <>
struct CompositePreconditionedSolver<1>
{
private:
    using Backend = amgcl::backend::builtin<double>;
    using Solver = amgcl::make_solver<
        amgcl::runtime::preconditioner<Backend>,
        amgcl::runtime::solver::wrapper<Backend>
    >;

public:
    using Type = amgcl::make_solver<
        amgcl::preconditioner::schur_pressure_correction<Solver,Solver>,
        amgcl::runtime::solver::wrapper<Backend>
    >;
};



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


/** @brief Scatter components of a vector into two other vectors.
 *  @warning Segfaults if the mask @a rMask contains items other than 0 or 1.
 */
void MakeSubblocks(const TUblasSparseSpace<double>::VectorType& rRootVector,
                   const CPTraits::Mask& rMask,
                   const std::array<TUblasSparseSpace<double>::VectorType*,2>& rOutput)
{
    KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);
    KRATOS_ERROR_IF_NOT(rRootVector.size() == rMask.size())
        << "input vector size " << rRootVector.size()
        << " does not match the mask size " << rMask.size() << "\n";

    KRATOS_TRY

    // Make sure that all input/output vectors are distinct
    {
        std::set<const TUblasSparseSpace<double>::VectorType*> vectors {
            &rRootVector,
            rOutput[0],
            rOutput[1]
        };
        KRATOS_ERROR_IF_NOT(vectors.size() == 3) << "all input/output vectors must be distinct";
    }

    // Resize output vectors
    rOutput[0]->clear();
    rOutput[1]->clear();
    const std::size_t lower_count = std::count_if(rMask.begin(),
                                                  rMask.end(),
                                                  [](char v) -> bool {
                                                    return v;
                                                  });
    const std::size_t higher_count = rMask.size() - lower_count;
    rOutput[0]->resize(lower_count);
    rOutput[1]->resize(higher_count);

    // Copy components to the output vectors
    std::array<std::size_t,2> i_subs {0ul, 0ul}; // <== {i_lower, i_higher}
    for (std::size_t i_root=0ul; i_root<rRootVector.size(); ++i_root) {
        const char is_higher = !rMask[i_root];
        rOutput[is_higher]->data()[i_subs[is_higher]] = rRootVector[i_root];
        ++i_subs[is_higher];
    }

    KRATOS_CATCH("")
}


} // namespace detail



template <unsigned BlockSize>
using MakeCompositePreconditionedSolver = typename detail::CompositePreconditionedSolver<BlockSize>::Type;



template <class TSparseSpace,
          class TDenseSpace,
          class TReorderer>
struct AMGCLHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::Impl
{
    double mTolerance;

    int mVerbosity;

    std::size_t mDoFCount;

    // Indicates whether the variable at the given index is lower order.
    std::optional<std::vector<char>> mLowerDofMask;

    boost::property_tree::ptree mAMGCLSettings;
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
    mpImpl->mTolerance = safe_parameters["amgcl_settings"]["solver"]["tol"].Get<double>();
    mpImpl->mVerbosity = safe_parameters["verbosity"].Get<int>();
    mpImpl->mDoFCount = 0ul;

    std::stringstream json_stream;
    json_stream << safe_parameters["amgcl_settings"].GetInput().PrettyPrintJsonString();
    boost::property_tree::read_json(
        json_stream,
        mpImpl->mAMGCLSettings
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

    mpImpl->mAMGCLSettings.put("precond.pmask_size", r_lower_dof_mask.size());
    mpImpl->mAMGCLSettings.put("precond.pmask", static_cast<void*>(r_lower_dof_mask.data()));
    mpImpl->mAMGCLSettings.put("solver.verbose", 1 < mpImpl->mVerbosity);

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
    "scaling": false,
    "schur_variable" : "PRESSURE",
    "amgcl_settings" : {
        "solver": {
            "type": "lgmres",
            "M": 50,
            "maxiter": 1000,
            "tol": 1e-8,
            "verbose": true
        },
        "precond": {
            "pmask_size": -1,
            "adjust_p": 0,
            "type": 2,
            "usolver": {
                "solver": {
                    "type": "preonly"
                },
                "precond": {
                    "relax": {
                        "type": "ilup"
                    },
                    "coarsening": {
                        "type": "aggregation",
                        "aggr": {
                            "eps_strong": 0
                        }
                    }
                }
            },
            "psolver": {
                "solver": {
                    "type": "preonly"
                }
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
    Vector b_l, b_h;
    detail::MakeSubblocks(rA, mpImpl->mLowerDofMask.value(), {&A_ll, &A_lh, &A_hl, &A_hh});
    detail::MakeSubblocks(rB, mpImpl->mLowerDofMask.value(), {&b_l, &b_h});

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

    // Construct solver
    auto solver = MakeCompositePreconditionedSolver<BlockSize>(*p_a, mpImpl->mAMGCLSettings);
    KRATOS_INFO_IF("AMGCLHierarchicalSolver", 1 < mpImpl->mVerbosity)
        << "Solver memory usage: " << amgcl::human_readable_memory(amgcl::backend::bytes(solver))
        << "\n";

    // Solve
    return solver(*p_a, rB, rX);
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
        << "tolerance: " << mpImpl->mTolerance << "\n"
        << "DoF size: " << mpImpl->mDoFCount << "\n"
        << "verbosity: " << mpImpl->mVerbosity << "\n"
        << "AMGCL settings: "
        ;
    boost::property_tree::json_parser::write_json(rStream, mpImpl->mAMGCLSettings);
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
