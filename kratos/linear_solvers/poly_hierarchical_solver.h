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

    /**
     * Multi solve method for solving a set of linear systems with same coefficient matrix.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrix& rA, DenseMatrix& rX, DenseMatrix& rB) override;

    void  PrintInfo(std::ostream& rOStream) const override;

    void  PrintData(std::ostream& rOStream) const override;

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function tells if the solver requires such data
     */
    bool AdditionalPhysicalDataIsNeeded() override
    {
        return true;
    }

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function is the place to eventually provide such data
     */
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
