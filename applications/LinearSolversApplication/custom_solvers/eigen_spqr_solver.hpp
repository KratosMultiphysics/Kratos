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
#ifdef KRATOS_USE_EIGEN_SUITESPARSE

// External includes
#include <Eigen/Sparse>
#include <Eigen/SPQRSupport>

// LinearSolversApplication includes
#include "linear_solvers_define.h" // EigenSparseMatrix, EigenDynamicVector

// Core includes
#include "includes/kratos_parameters.h" // Parameters


namespace Kratos {


template <typename TScalar>
class EigenSPQRSolver
{
public:
    using Scalar = TScalar;
    using SparseMatrix = EigenSparseMatrix<Scalar>;
    using Vector = EigenDynamicVector<Scalar>;

    void Initialize(Parameters settings)
    {
    }

    bool Compute(Eigen::Map<const SparseMatrix> rLhs)
    {
        KRATOS_TRY

        mDecomposition.compute(rLhs);
        const auto report = mDecomposition.info();

        switch (report) {
            case Eigen::ComputationInfo::Success:
                break;
            case Eigen::ComputationInfo::NumericalIssue:
                KRATOS_ERROR << "SPQR factorization ran into a numerical issue";
            case Eigen::ComputationInfo::NoConvergence:
                KRATOS_ERROR << "SPQR factorization failed to converge";
            case Eigen::ComputationInfo::InvalidInput:
                KRATOS_ERROR << "invalid input was provided for SPQR";
        } // switch report
        return true;

        KRATOS_CATCH("")
    }

    bool Solve(Eigen::Ref<const Vector> rRhs, Eigen::Ref<Vector> rSolution) const
    {
        rSolution = mDecomposition.solve(rRhs);
        const auto report = mDecomposition.info();

        switch (report) {
            case Eigen::ComputationInfo::Success:
                break;
            case Eigen::ComputationInfo::NumericalIssue:
                KRATOS_ERROR << "SPQR solver ran into a numerical issue";
            case Eigen::ComputationInfo::NoConvergence:
                KRATOS_ERROR << "SPQR solver failed to converge";
            case Eigen::ComputationInfo::InvalidInput:
                KRATOS_ERROR << "invalid input was provided for SPQR";
        } // switch report
        return true;
    }

    void PrintInfo(std::ostream &rOStream) const
    {
    }

    std::string GetSolverErrorMessages() const
    {
        return "";
    }

    static std::string Name()
    {
        return "spqr";
    }

private:
    Eigen::SPQR<SparseMatrix> mDecomposition;
};


} // namespace Kratos

#endif
