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
#include "linear_solvers/linear_solver.h" // LinearSolver
#include "includes/kratos_parameters.h" // Parameters
#include "includes/code_location.h" // KRATOS_CODE_LOCATION

// STL includes
#include <memory> // std::unique_ptr


namespace Kratos {


/// @brief Base class for MKL smoothers.
/// @internal
template <class TSparse, class TDense>
class KRATOS_API(LINEARSOLVERS_APPLICATION) MKLSmootherBase : public LinearSolver<TSparse,TDense>
{
public:
    /// @brief View over a CSR matrix.
    struct CSRView
    {
        using Value = typename TSparse::DataType;

        /// @brief Sparse index type that must be identical @p MKL_INT.
        using Index = int;

        Index row_count;
        Index column_count;
        Index entry_count;
        const Index* it_row_begin;
        const Index* it_column_begin;
        const Value* it_entry_begin;
    }; // CSRView

    /// @brief View over a contiguous dense array.
    template <bool IsMutable>
    struct VectorView
    {
        using Value = std::conditional_t<
            IsMutable,
            typename TSparse::DataType,
            const typename TSparse::DataType>;

        int size;
        Value* it_begin;
    }; // struct VectorView

    using Base = LinearSolver<TSparse,TDense>;

    using Base::Solve;

    using SparseMatrix = typename TSparse::MatrixType;

    using Vector = typename TSparse::VectorType;

    using DenseMatrix = typename TDense::MatrixType;

    MKLSmootherBase();

    ~MKLSmootherBase();

    /// @copydoc Base::InitializeSolutionStep
    void InitializeSolutionStep(SparseMatrix& rLhs, Vector& rSolution, Vector& rRhs) override;

    /// @copydoc Base::PerformSolutionStep
    bool PerformSolutionStep(SparseMatrix& rLhs, Vector& rSolution, Vector& rRhs) final override;

    /// @copydoc Base::Solve(SparseMatrix&,Vector&,Vector&)
    bool Solve(SparseMatrix& rLhs, Vector& rSolution, Vector& rRhs) final override;

    /// @copydoc Base::Solve(SparseMatrix&,DenseMatrix&,DenseMatrix&)
    bool Solve(SparseMatrix& rLhs, DenseMatrix& rSolutions, DenseMatrix& rRhsVectors) final override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not supported";}

    void FinalizeSolutionStep(SparseMatrix& rLhs,
                              Vector& rSolution,
                              Vector& rRhs) override;

    void Clear() override;

protected:
    std::tuple<
        CSRView,
        VectorView</*IsMutable=*/true>,
        VectorView</*IsMutable=*/false>>
    MakeSystemView(const SparseMatrix& rLhs,
                   typename TSparse::DataType* itSolutionBegin,
                   typename TSparse::DataType* itSolutionEnd,
                   const typename TSparse::DataType* itRhsBegin,
                   const typename TSparse::DataType* itRhsEnd) const;

    virtual bool Solve(CSRView Lhs,
                       VectorView</*IsMutable=*/true> Solution,
                       VectorView</*IsMutable=*/false> Rhs) = 0;

private:
    MKLSmootherBase(MKLSmootherBase&&) = delete;

    MKLSmootherBase(const MKLSmootherBase&) = delete;

    MKLSmootherBase& operator=(MKLSmootherBase&&) = delete;

    MKLSmootherBase& operator=(const MKLSmootherBase&) = delete;

    struct Impl;
    class std::unique_ptr<Impl> mpImpl;
}; // class MKLSmootherBase


} // namespace Kratos
