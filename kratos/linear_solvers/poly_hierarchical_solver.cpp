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
#include "poly_hierarchical_solver.h"
#include "spaces/ublas_space.h"
#include "utilities/profiler.h"
#include "factories/linear_solver_factory.h"

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
                   const std::vector<char>& rMask,
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
struct PolyHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::Impl
{
    double mTolerance;

    int mVerbosity;

    std::size_t mMaxIterations;

    // Indicates whether the variable at the given index is in the coarse space.
    std::optional<std::vector<char>> mCoarseMask;

    typename LinearSolver<TSparseSpace,TDenseSpace/*reorderer is omitted on purpose*/>::Pointer mpCoarseSolver;

    typename LinearSolver<TSparseSpace,TDenseSpace/*reorderer is omitted on purpose*/>::Pointer mpFineSolver;

    typename TSparseSpace::MatrixType mCoarseMatrix;

    typename TSparseSpace::MatrixType mOffDiagonal;
}; // struct PolyHierarchicalSolver::Impl



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
PolyHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::PolyHierarchicalSolver(Parameters parameters)
    : mpImpl(new Impl)
{
    KRATOS_TRY
    Parameters default_parameters = this->GetDefaultParameters();
    parameters.ValidateAndAssignDefaults(default_parameters);

    KRATOS_ERROR_IF_NOT(parameters["solver_type"].GetString() == "poly_hierarchical")
        << "Requested a(n) '" << parameters["solver_type"].GetString() << "' solver,"
        << " but constructing an PolyHierarchicalSolver";

    detail::ParametersWithDefaults safe_parameters(parameters, default_parameters);
    mpImpl->mTolerance = safe_parameters["tolerance"].Get<double>();
    mpImpl->mVerbosity = safe_parameters["verbosity"].Get<int>();
    mpImpl->mMaxIterations = safe_parameters["max_iterations"].Get<int>();

    // Construct subsolvers
    KRATOS_TRY
    LinearSolverFactory<TSparseSpace,TDenseSpace> solver_factory;
    mpImpl->mpCoarseSolver = solver_factory.Create(safe_parameters["coarse_settings"].GetInput());
    mpImpl->mpFineSolver = solver_factory.Create(safe_parameters["fine_settings"].GetInput());
    KRATOS_CATCH("")

    KRATOS_CATCH("")
}



// Necessary for PIMPL
template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
PolyHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::~PolyHierarchicalSolver()
{
}



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
bool PolyHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::Solve(SparseMatrix& rA,
                                                                         Vector& rX,
                                                                         Vector& rB)
{
    KRATOS_ERROR_IF_NOT(mpImpl->mCoarseMask.has_value())
        << "PolyHierarchicalSolver::Solve called before PolyHierarchicalSolver::ProvideAdditionalData";
    auto& r_coarse_mask = mpImpl->mCoarseMask.value();
    KRATOS_ERROR_IF_NOT(r_coarse_mask.size() == TSparseSpace::Size1(rA))
        << "DoF mask size " << r_coarse_mask.size()
        << " does not match system size " << TSparseSpace::Size1(rA);

    KRATOS_TRY

    // Dump the generated matrices to disk f requested.
    if (4 <= mpImpl->mVerbosity) {
        KRATOS_ERROR << "verbosity >=4 prints system matrices and terminates\n";
    }

    // Construct block maps for vectors
    const auto& r_mask = mpImpl->mCoarseMask.value();
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
        //for (std::size_t i_root=0ul; i_root<rRoot.size(); ++i_root) {
        IndexPartition<std::size_t>(rRoot.size()).for_each([&](std::size_t i_root){
            const char mask_value = r_mask[i_root];
            const std::size_t i_block = !mask_value;
            const std::size_t i_sub = block_local_indices[i_root];
            rSubs[i_block][i_sub] = rRoot[i_root];
        });
        //}
    };

    const auto aggregate = [&block_local_indices, &r_mask](const std::array<Vector,2>& rSubs,
                                                           Vector& rRoot) {
        //for (std::size_t i_root=0ul; i_root<rRoot.size(); ++i_root) {
        IndexPartition<std::size_t>(rRoot.size()).for_each([&](std::size_t i_root){
            const char mask_value = r_mask[i_root];
            const std::size_t i_block = !mask_value;
            const std::size_t i_sub = block_local_indices[i_root];
            rRoot[i_root] += rSubs[i_block][i_sub];
        });
        //}
    };

    // Solve
    std::array<Vector,2> r, d; // <== {r_l, r_h}, {d_l, d_h}
    r[0].resize(  masked_count);
    r[1].resize(unmasked_count);
    d[0].resize(  masked_count);
    d[1].resize(unmasked_count);

    Vector delta(full_size, 0.0);
    Vector coarse_residual(masked_count, 0.0);
    Vector coarse_delta(masked_count, 0.0);
    Vector residual = rB;

    KRATOS_INFO_IF("PolyHierarchicalSolver", 1 <= mpImpl->mVerbosity)
        << "coarse size: " << masked_count
        << " fine size: " << unmasked_count << "\n";

    const double rhs_norm = norm_2(rB);
    KRATOS_ERROR_IF_NOT(rhs_norm);

    std::tuple<std::size_t,double> results {0ul, 1.0};
    const std::size_t max_iterations = mpImpl->mMaxIterations;
    for (std::size_t i_iteration=0ul; i_iteration<max_iterations; ++i_iteration) {
        // Perform a few relaxation passes on the full system
        KRATOS_TRY
        mpImpl->mpFineSolver->Solve(rA, delta, residual);
        KRATOS_CATCH("")
        coarsen(delta, d);

        // Compute restricted residual
        coarsen(residual, r);
        noalias(coarse_residual) = r[0] - prod(mpImpl->mCoarseMatrix, d[0]) - prod(mpImpl->mOffDiagonal, d[1]);

        // Solve coarse
        KRATOS_TRY
        mpImpl->mpCoarseSolver->Solve(mpImpl->mCoarseMatrix, coarse_delta, coarse_residual);
        KRATOS_CATCH("")

        if (2 <= mpImpl->mVerbosity) {
            const double diff_norm = norm_2(coarse_delta) / norm_2(d[0]);
            KRATOS_INFO("PolyHierarchicalSolver")
                << "norm of coarse order correction: " << diff_norm << "\n";
        }

        // Update solution
        noalias(d[0]) += coarse_delta;
        aggregate(d, rX);

        // Compute the residual
        noalias(residual) = rB - prod(rA, rX);

        // Update status
        std::get<0>(results)++;
        std::get<1>(results) = norm_2(residual) / rhs_norm;
        KRATOS_INFO_IF("PolyHierarchicalSolver", 1 <= mpImpl->mVerbosity)
            << "iteration " << std::get<0>(results)
            << " residual " << std::get<1>(results) << "\n";

        // Early exit if converged
        if (std::get<1>(results) < mpImpl->mTolerance) {
            break;
        }
    }

    return std::get<1>(results) < mpImpl->mTolerance;
    KRATOS_CATCH("")
}



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
void PolyHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::ProvideAdditionalData(SparseMatrix& rA,
                                                                                        Vector& rX,
                                                                                        Vector& rB,
                                                                                        ModelPart::DofsArrayType& rDofs,
                                                                                        ModelPart& rModelPart)
{
    KRATOS_TRY
    KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);

    // Construct coarse mask
    if (!mpImpl->mCoarseMask.has_value()) {
        mpImpl->mCoarseMask.emplace();
    }
    auto& r_coarse_mask = mpImpl->mCoarseMask.value();
    r_coarse_mask.resize(TSparseSpace::Size1(rA), false);

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
                    r_coarse_mask[equation_id] = true;
                }
                break;
            }
            default: {
                KRATOS_ERROR << "Unsupported element geometry: " << r_geometry;
            }
        } // switch GeometryType
    }

    // Separate components related to coarse and higher order DoFs.
    SparseMatrix A_hl, A_hh;
    detail::MakeSubblocks(rA,
                          mpImpl->mCoarseMask.value(),
                          {&mpImpl->mCoarseMatrix, &mpImpl->mOffDiagonal, &A_hl, &A_hh});

    // Dump the generated matrices to disk f requested.
    if (4 <= mpImpl->mVerbosity) {
        KRATOS_INFO("PolyHierarchicalSolver")
            << "writing system matrices A.mm b.mm\n";
        TSparseSpace::WriteMatrixMarketMatrix("A.mm", rA, false);
        TSparseSpace::WriteMatrixMarketVector("b.mm", rB);
        KRATOS_INFO("PolyHierarchicalSolver")
            << "writing submatrices: A_ll.mm A_lh.mm A_hl.mm A_hh.mm\n";
        TSparseSpace::WriteMatrixMarketMatrix("A_ll.mm", mpImpl->mCoarseMatrix, false);
        TSparseSpace::WriteMatrixMarketMatrix("A_lh.mm", mpImpl->mOffDiagonal, false);
        TSparseSpace::WriteMatrixMarketMatrix("A_hl.mm", A_hl, false);
        TSparseSpace::WriteMatrixMarketMatrix("A_hh.mm", A_hh, false);
    }

    if (mpImpl->mpCoarseSolver->AdditionalPhysicalDataIsNeeded()) {
        mpImpl->mpCoarseSolver->ProvideAdditionalData(mpImpl->mCoarseMatrix, rX, rB, rDofs, rModelPart);
    }

    if (mpImpl->mpFineSolver->AdditionalPhysicalDataIsNeeded()) {
        mpImpl->mpFineSolver->ProvideAdditionalData(rA, rX, rB, rDofs, rModelPart);
    }

    KRATOS_CATCH("")
}



template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
Parameters
PolyHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::GetDefaultParameters()
{
    return Parameters(R"(
{
    "solver_type" : "poly_hierarchical",
    "verbosity" : 0,
    "max_iterations" : 50,
    "tolerance" : 1e-6,
    "coarse_settings" : {},
    "fine_settings" : {}
}
    )");
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
bool PolyHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::Solve(SparseMatrix& rA,
                                                                        DenseMatrix& rX,
                                                                        DenseMatrix& rB)
{
    return false;
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
void PolyHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::PrintInfo(std::ostream& rStream) const
{
    rStream << "PolyHierarchicalSolver";
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
void PolyHierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>::PrintData(std::ostream& rStream) const
{
    rStream
        << (mpImpl->mCoarseMask.has_value() ? "initialized" : "uninitialized") << "\n"
        << "tolerance     : " << mpImpl->mTolerance << "\n"
        << "max iterations: " << mpImpl->mMaxIterations << "\n"
        << "verbosity     : " << mpImpl->mVerbosity << "\n"
        ;
}



template
class PolyHierarchicalSolver<
    TUblasSparseSpace<double>,
    TUblasDenseSpace<double>,
    Reorderer<
        TUblasSparseSpace<double>,
        TUblasDenseSpace<double>
    >
>;


} // namespace Kratos
