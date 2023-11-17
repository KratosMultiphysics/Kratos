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
class PolyHierarchicalSolver final : public LinearSolver<TSparseSpace,
                                                       TDenseSpace,
                                                       TReorderer>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PolyHierarchicalSolver);

    using Base =  LinearSolver<TSparseSpace, TDenseSpace, TReorderer>;

    using SparseMatrix = typename TSparseSpace::MatrixType;

    using Vector = typename TSparseSpace::VectorType;

    using DenseMatrix = typename TDenseSpace::MatrixType;

    PolyHierarchicalSolver(Parameters rParameters);

    PolyHierarchicalSolver(PolyHierarchicalSolver&&) noexcept = default;

    ~PolyHierarchicalSolver() override;

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
    PolyHierarchicalSolver(const PolyHierarchicalSolver& Other) = delete;

    struct Impl;
    std::unique_ptr<Impl> mpImpl;
}; // class PolyHierarchicalSolver


template<class TSparseSpace, class TDenseSpace,class TReorderer>
inline std::istream& operator >> (std::istream& rIStream, PolyHierarchicalSolver< TSparseSpace,
                                  TDenseSpace, TReorderer>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpace, class TDenseSpace, class TReorderer>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const PolyHierarchicalSolver<TSparseSpace,
                                  TDenseSpace, TReorderer>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos
