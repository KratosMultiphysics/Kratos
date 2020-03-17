/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Authors: Thomas Oberbichler
*/

#if !defined(KRATOS_MAKE_CONSISTENT_H_INCLUDED)
#define KRATOS_MAKE_CONSISTENT_H_INCLUDED

// External includes
#include <Eigen/Core>
#include <Eigen/Sparse>

// Project includes
#include "includes/define.h"
#include "custom_utilities/ublas_wrapper.h"
#include "custom_solvers/eigen_direct_solver.h"
#include "custom_solvers/eigen_pardiso_lu_solver.h"
#include "custom_solvers/eigen_sparse_cg_solver.h"

namespace Kratos
{

typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
typedef SparseSpaceType::MatrixType SparseMatrixType;
typedef SparseSpaceType::VectorType VectorType;

void MakeConsistent(SparseMatrixType& matrix, VectorType& x){

    x.resize(matrix.size2());

    UblasWrapper<double> matrix_wrapper(matrix);
    const auto& m = matrix_wrapper.matrix(); // todo make csc

    SparseMatrixType LHS = matrix;

    std::vector<int> m_index1(LHS.index1_data().begin(), LHS.index1_data().end());
    std::vector<int> m_index2(LHS.index2_data().begin(), LHS.index2_data().end());
    Eigen::Map<Eigen::SparseMatrix<double, Eigen::RowMajor, int>> lhs(
        LHS.size1(),
        LHS.size2(),
        LHS.nnz(),
        m_index1.data(),
        m_index2.data(),
        LHS.value_data().begin()
    );

    Eigen::VectorXd col_sums = m.transpose() * Eigen::VectorXd::Ones(m.rows());

    int nz_count = 0;
    // KRATOS_WATCH(m.outerSize())
// #pragma omp parallel for // TODO nz_count
    for (int i=0; i<m.outerSize(); ++i) { // row
        // if (!(i%100)) KRATOS_WATCH(i);
        for (int inner = m.outerIndexPtr()[i]; inner<m.outerIndexPtr()[i+1]; ++inner) {
            const int j = m.innerIndexPtr()[inner];
            const double& m_ij = m.valuePtr()[nz_count];
            const double v = col_sums[j] * m_ij;
            lhs.valuePtr()[nz_count] = v;
            nz_count++;
        }
    } // TODO check if this works

    // very slow loop over all values including zeros
    // for (int i =0; i < int(matrix.size1()); ++i) {
    //     if (!i%100) KRATOS_WATCH(i);
    //     for (int j =0; j < int(matrix.size2()); ++j) {
    //         const double& m_ij = matrix(i,j);
    //         if (m_ij == 0.0) continue;
    //         const double v = m_ij*2;//(m_ij * m.col(j)).sum();
    //         if (v == 0) continue;
    //         LHS(i,j) = v;
    //     }
    // }

    VectorType ones;
    ones.resize(x.size(), 1.0);
    for (int i=0; i<ones.size(); ++i){
        ones[i] = 1.0;
    }

    KRATOS_WATCH("Start solving")
    // EigenDirectSolver<EigenSparseCGSolver<double>> solver; //4356x4356: - does not converge...
    // solver.Solve(LHS, x, ones);
    EigenDirectSolver<EigenPardisoLUSolver<double>> solver; //4356x4356: 12s
    solver.Solve(LHS, x, ones);

    Eigen::Map<Eigen::VectorXd> x_vec(x.data().begin(), x.size());
    x_vec = x_vec.cwiseSqrt();
    Eigen::SparseMatrix<double, Eigen::RowMajor, int> diag(x_vec.asDiagonal());
    Eigen::SparseMatrix<double, Eigen::RowMajor, int> _m = m * diag;

    for (int i=0; i<m.nonZeros(); ++i){
        matrix.value_data()[i] = _m.valuePtr()[i];
    }
}


} // namespace Kratos

#endif // defined(KRATOS_MAKE_CONSISTENT_H_INCLUDED)