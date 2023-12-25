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

// External includes

// Project includes
#include "includes/ublas_interface.h"
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"


namespace Kratos {


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer = Reorderer<TSparseSpace, TDenseSpace> >
class HierarchicalSolver final : public LinearSolver<TSparseSpace,
                                                     TDenseSpace,
                                                     TReorderer>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(HierarchicalSolver);

    using Base =  LinearSolver<TSparseSpace, TDenseSpace, TReorderer>;

    using SparseMatrix = typename TSparseSpace::MatrixType;

    using Vector = typename TSparseSpace::VectorType;

    using DenseMatrix = typename TDenseSpace::MatrixType;

    HierarchicalSolver();

    HierarchicalSolver(Parameters rParameters);

    HierarchicalSolver(HierarchicalSolver&&) noexcept = default;

    ~HierarchicalSolver() override;

    /// @copydoc LinearSolver::Solve
    bool Solve(SparseMatrix& rA, Vector& rX, Vector& rB) override;

    /// @copydoc LinearSolver::Solve
    bool Solve(SparseMatrix& rA, DenseMatrix& rX, DenseMatrix& rB) override;

    /// @copydoc LinearSolver::PrintInfo
    void PrintInfo(std::ostream& rOStream) const override;

    /// @copydoc LinearSolver::PrintData
    void PrintData(std::ostream& rOStream) const override;

    /// @copydoc LinearSolver::AdditionalPhysicalDataIsNeeded
    bool AdditionalPhysicalDataIsNeeded() override
    {
        return true;
    }

    /// @copydoc LinearSolver::ProvideAdditionalData
    void ProvideAdditionalData (
        SparseMatrix& rA,
        Vector& rX,
        Vector& rB,
        typename ModelPart::DofsArrayType& rDofSet,
        ModelPart& rModelPart
    ) override;

    static Parameters GetDefaultParameters();

private:
    HierarchicalSolver(const HierarchicalSolver& Other) = delete;

    struct Impl;
    std::unique_ptr<Impl> mpImpl;
}; // class HierarchicalSolver


template<class TSparseSpace, class TDenseSpace,class TReorderer>
std::istream& operator >> (std::istream& rStream,
                           HierarchicalSolver< TSparseSpace,TDenseSpace,TReorderer>& rSolver)
{
    KRATOS_ERROR << "deserializing a HierarchicalSolver from a string stream is not supported";
    return rStream;
}

/**
 * output stream function
 */
template<class TSparseSpace, class TDenseSpace, class TReorderer>
std::ostream& operator << (std::ostream& rStream,
                           const HierarchicalSolver<TSparseSpace,TDenseSpace,TReorderer>& rSolver)
{
    rSolver.PrintInfo(rStream);
    rStream << "\n";
    rSolver.PrintData(rStream);
    return rStream;
}

}  // namespace Kratos
