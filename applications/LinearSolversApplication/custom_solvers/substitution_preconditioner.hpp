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

// --- Kratos Core Includes ---
#include "linear_solvers/preconditioner.h"

// --- STL Includes ---
#include <memory> // std::unique_ptr, std::shared_ptr


namespace Kratos {


template <class TSparse, class TDense>
class SubstitutionPreconditioner final : public Preconditioner<TSparse,TDense> {
public:
    KRATOS_CLASS_POINTER_DEFINITION(SubstitutionPreconditioner);

    SubstitutionPreconditioner();

    SubstitutionPreconditioner(std::shared_ptr<typename TSparse::MatrixType> pTriangle);

    SubstitutionPreconditioner(
        std::shared_ptr<typename TSparse::MatrixType> pLowerTriangle,
        std::shared_ptr<typename TSparse::MatrixType> pUpperTriangle);

    SubstitutionPreconditioner(SubstitutionPreconditioner&& rRhs) noexcept;

    SubstitutionPreconditioner(const SubstitutionPreconditioner&) = delete;

    ~SubstitutionPreconditioner();

    SubstitutionPreconditioner& operator=(SubstitutionPreconditioner&& rRhs) noexcept;

    SubstitutionPreconditioner& operator=(const SubstitutionPreconditioner) = delete;

    void Mult(
        typename TSparse::MatrixType& rLhs,
        typename TSparse::VectorType& rSolution,
        typename TSparse::VectorType& rRhs) override;

    void TransposeMult(
        typename TSparse::MatrixType& rLhs,
        typename TSparse::VectorType& rSolution,
        typename TSparse::VectorType& rRhs) override;

    typename TSparse::VectorType& ApplyLeft(typename TSparse::VectorType& rSolution) override;

    typename TSparse::VectorType& ApplyTransposeLeft(typename TSparse::VectorType& rSolution) override;

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
}; // class SubstitutionPreconditioner


} // namespace Kratos
