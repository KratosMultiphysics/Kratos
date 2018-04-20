//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//
//

// System includes
#include <iostream>
#include <complex>
#include <vector>
#include <algorithm>

// External includes
#include <memory>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <amgcl/value_type/complex.hpp>
#include <amgcl/solver/skyline_lu.hpp>

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/direct_solver.h"

#if !defined(KRATOS_SKYLINE_LU_CUSTOM_SCALAR_SOLVER_H_INCLUDED)
#define  KRATOS_SKYLINE_LU_CUSTOM_SCALAR_SOLVER_H_INCLUDED

namespace Kratos {

///@name Kratos Classes
///@{

template <class TSparseSpaceType, class TDenseSpaceType, class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class SkylineLUCustomScalarSolver
    : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SkylineLUCustomScalarSolver);

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename TSparseSpaceType::DataType DataType;

    typedef typename amgcl::backend::builtin<DataType>::matrix BuiltinMatrixType;

    typedef amgcl::solver::skyline_lu<DataType> SolverType;

    SkylineLUCustomScalarSolver()
        : DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>()
    {
    }

    SkylineLUCustomScalarSolver(Parameters& rParam)
        : DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>(rParam)
    {
    }

    ~SkylineLUCustomScalarSolver() override
    {
        Clear();
    }

    void InitializeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        Clear();

        pBuiltinMatrix = amgcl::adapter::zero_copy(
                rA.size1(),
                rA.index1_data().begin(),
                rA.index2_data().begin(),
                rA.value_data().begin());

        pSolver = Kratos::make_shared<SolverType>(*pBuiltinMatrix);
    }

    void PerformSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        std::vector<DataType> x(rX.size());
        std::vector<DataType> b(rB.size());

        std::copy(std::begin(rB), std::end(rB), std::begin(b));

        (*pSolver)(b, x);

        std::copy(std::begin(x), std::end(x), std::begin(rX));
    }

    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        InitializeSolutionStep(rA, rX, rB);
        PerformSolutionStep(rA, rX, rB);
        FinalizeSolutionStep(rA, rX, rB);

        return true;
    }

    void FinalizeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        Clear();
    }

    void Clear() override
    {
        pSolver.reset();
        pBuiltinMatrix.reset();
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Skyline LU custom scalar solver";
    }

private:

    Kratos::shared_ptr<BuiltinMatrixType> pBuiltinMatrix;

    Kratos::shared_ptr<SolverType> pSolver;
};

///@}

///@name Input and output
///@{

/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >>(std::istream& rIStream,
        SkylineLUCustomScalarSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    return rIStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator <<(std::ostream& rOStream,
        const SkylineLUCustomScalarSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}// namespace Kratos.

#endif // KRATOS_SKYLINE_LU_CUSTOM_SCALAR_SOLVER_H_INCLUDED  defined
