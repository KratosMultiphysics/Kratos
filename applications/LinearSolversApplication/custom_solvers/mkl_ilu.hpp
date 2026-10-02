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
#include "custom_solvers/mkl_smoother_base.hpp" // MKLSmootherBase

// STL includes
#include <memory> // std::unique_ptr
#include <vector> // std::vector
#include <array> // std::array


namespace Kratos {


/// @brief Base class for @ref MKLILU0Smoother and @ref MKLILUTSmoother
/// @internal
template <class TSparse, class TDense>
class KRATOS_API(LINEARSOLVERS_APPLICATION) MKLILUSmootherBase
    : public MKLSmootherBase<TSparse,TDense>
{
private:
    using Base = MKLSmootherBase<TSparse,TDense>;

    using BaseLinearSolver = LinearSolver<TSparse,TDense>;

    using BaseLinearSolver::Solve;

public:
    MKLILUSmootherBase();

    MKLILUSmootherBase(Parameters Settings);

    ~MKLILUSmootherBase();

    void InitializeSolutionStep(typename Base::SparseMatrix& rLhs,
                                typename Base::Vector& rSolution,
                                typename Base::Vector& rRhs) final override;

    static Parameters GetDefaultParameters();

    void Clear() final override;

protected:
    bool Solve(typename Base::CSRView LhsView,
               typename Base::template VectorView</*IsMutable=*/true> SolutionView,
               typename Base::template VectorView</*IsMutable=*/false> RhsView) final override;

    virtual void Factorize(std::vector<int>& rRowExtents,
                           std::vector<int>& rColumnIndices,
                           std::vector<typename TSparse::DataType>& rEntries,
                           typename Base::CSRView LhsView,
                           const std::array<int,128>& rIntegerSettings,
                           const std::array<typename TSparse::DataType,128>& rNumericSettings) = 0;

    typename TSparse::DataType GetRelaxation() const noexcept;

    int GetIterations() const noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
}; // class MKLILUSmootherBase


/** @brief ILU0 smoother/preconditioner/solver from MKL.
 *  @details This solver is intended to be used as a smoother, performing a fixed number of iterations.
 *           Default parameters:
 *           @code
 *           {
 *              "solver_type" : "mkl_ilu0",
 *              "iterations"  : 1,
 *              "relaxation"  : 1.0
 *           }
 *           @endcode
 *           - @p "solver_type" Name of the @ref LinearSolver "solver" to refer to from the JSON layer.
 *           - @p "iterations" Perform this many sweeps.
 *           - @p "relaxation" Damping to apply on the right hand side.
 *  @see https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2025-1/dcsrilu0.html
 *  @see https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2025-1/mkl-sparse-trsv.html
 */
template <class TSparse, class TDense>
class KRATOS_API(LINEARSOLVERS_APPLICATION) MKLILU0Smoother final
    : public MKLILUSmootherBase<TSparse,TDense>
{
private:
    using Base = MKLILUSmootherBase<TSparse,TDense>;

public:
    KRATOS_CLASS_POINTER_DEFINITION(MKLILU0Smoother);

    MKLILU0Smoother();

    MKLILU0Smoother(Parameters Settings);

    static Parameters GetDefaultParameters();

protected:
    void Factorize(std::vector<int>& rRowExtents,
                   std::vector<int>& rColumnIndices,
                   std::vector<typename TSparse::DataType>& rEntries,
                   typename Base::CSRView LhsView,
                   const std::array<int,128>& rIntegerSettings,
                   const std::array<typename TSparse::DataType,128>& rNumericSettings) override;
}; // class MKLILU0Smoother


/** @brief ILUT smoother/preconditioner/solver from MKL.
 *  @details This solver is intended to be used as a smoother, performing a fixed number of iterations.
 *           Default parameters:
 *           @code
 *           {
 *              "solver_type"               : "mkl_ilut",
 *              "iterations"                : 1,
 *              "relaxation"                : 1.0,
 *              "fill_factor"               : 2,
 *              "factorization_tolerance"   : 1e-2
 *           }
 *           @endcode
 *           - @p "solver_type" Name of the @ref LinearSolver "solver" to refer to from the JSON layer.
 *           - @p "iterations" Perform this many sweeps.
 *           - @p "relaxation" Damping to apply on the right hand side.
 *           - @p "fill_factor" Maximum fill-in (half the bandwidth). Larger values improve stability
 *                              and accuracy but degrade performance and increase memory consumption.
 *           - @p "factorization_tolerance" Tolerance for the fill-in threshold during factorization.
 *                                          Larger values improve performance but degrade accuracy and
 *                                          stability.
 *  @see https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2025-1/dcsrilut.html
 *  @see https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2025-1/mkl-sparse-trsv.html
 */
template <class TSparse, class TDense>
class KRATOS_API(LINEARSOLVERS_APPLICATION) MKLILUTSmoother final
    : public MKLILUSmootherBase<TSparse,TDense>
{
private:
    using Base = MKLILUSmootherBase<TSparse,TDense>;

public:
    KRATOS_CLASS_POINTER_DEFINITION(MKLILUTSmoother);

    MKLILUTSmoother();

    MKLILUTSmoother(Parameters Settings);

    static Parameters GetDefaultParameters();

protected:
    void Factorize(std::vector<int>& rRowExtents,
                   std::vector<int>& rColumnIndices,
                   std::vector<typename TSparse::DataType>& rEntries,
                   typename Base::CSRView LhsView,
                   const std::array<int,128>& rIntegerSettings,
                   const std::array<typename TSparse::DataType,128>& rNumericSettings) override;

private:
    typename TSparse::DataType mFactorizationTolerance;

    int mFillFactor;
}; // class MKLILUTSmoother


} // namespace Kratos
