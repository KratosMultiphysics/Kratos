/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Authors: Matthias Ebert
//           based on code of
//             Thomas Oberbichler
//             Armin Geiser
*/

#if !defined(KRATOS_GEN_EIGENSYSTEM_SOLVER_H_INCLUDED)
#define KRATOS_GEN_EIGENSYSTEM_SOLVER_H_INCLUDED

// External includes
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/iterative_solver.h"
#include "utilities/openmp_utils.h"
#include "custom_utilities/ublas_wrapper.h"

namespace Kratos
{

template<
    //class TSolverType,
    class TSparseSpaceType, // = typename TSolverType::TGlobalSpace,
    class TDenseSpaceType, // = typename TSolverType::TLocalSpace,
    class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class GenEigensystemSolver
    : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
    Parameters mParam;

  public:
    KRATOS_CLASS_POINTER_DEFINITION(GenEigensystemSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    GenEigensystemSolver(
        Parameters param
    ) : mParam(param)
    {   
        Parameters default_params(R"(
        {
            "number_of_eigenvalues": 1,
            "compute_eigenvectors": true,
            "normalize_eigenvectors": false,
            "max_iteration": 400,
            "tolerance": 1e-10,
            "echo_level": 1
        })");

        mParam.ValidateAndAssignDefaults(default_params);

        BaseType::SetTolerance(mParam["tolerance"].GetDouble());
        BaseType::SetMaxIterationsNumber(mParam["max_iteration"].GetInt());
    }

    ~GenEigensystemSolver() override {}

    /**
     * Solve the generalized eigenvalue problem.
     *
     * K is a (non-)symmetric matrix. M is a symmetric positive-definite matrix.
     */
    void Solve(
        SparseMatrixType& rK,
        SparseMatrixType& rM,
        VectorType& rEigenvalues,
        DenseMatrixType& rEigenvectors) override
    {
        using scalar_t = double;
        using vector_t = Eigen::VectorXd;
        using matrix_t = Eigen::MatrixXd;

        // --- get settings

        const int nroot = mParam["number_of_eigenvalues"].GetInt();
        const int max_iteration = BaseType::GetMaxIterationsNumber();
        const double tolerance = BaseType::GetTolerance();
        const int echo_level = mParam["echo_level"].GetInt();
        //const int normalization_needed = mParam["normalize_eigenvectors"].GetBool();
        const int eigenvectors_needed = mParam["compute_eigenvectors"].GetBool();

KRATOS_INFO("---this is the new one----"); //TODO: delete before final version
        // --- wrap ublas matrices

        UblasWrapper<scalar_t> a_wrapper(rK);
        UblasWrapper<scalar_t> b_wrapper(rM);

        const auto& a = a_wrapper.matrix();
        const auto& b = b_wrapper.matrix();


        // --- timer

        double start_time = OpenMPUtils::GetCurrentTime();

        KRATOS_INFO_IF("GenEigensystemSolver", echo_level > 0) << "Start"  << std::endl;

        // --- calculation

        // TODO: check which is best in terms of time/memory
        // constructor which does not allocate the sizes in eig; needs resizing when calling eig.compute()
        Eigen::GeneralizedEigenSolver<matrix_t> eig;
        eig.setMaxIterations(max_iteration);
        eig.compute(a, b, eigenvectors_needed);

        // constructor which immediately calls eig.compute and does not need resizing 
        // setting max iterations of Schur decomposition is not possible
        //Eigen::GeneralizedEigenSolver<matrix_t> eig(a,b,true);

        // constructor which allocates storage needed in eig
        // however based on the code it gets resized when calling eig.compute anyway
        // also if calculating eigenvectors is not necessary, resizing part is skipped
        // -> allocating before (at construction time) would just need additional time and memory
        //Eigen::GeneralizedEigenSolver<matrix_t> eig(rK.size1());


            if(eig.info() != Eigen::Success) {
                KRATOS_WARNING("GenEigensystemSolver") << "Eigen solution was not successful!" << std::endl;
            }

        // debug print, to delete
        // // KRATOS_WATCH(nroot);
        // // KRATOS_WATCH(eig.eigenvalues());
        // // KRATOS_WATCH(eig.eigenvectors());
        // // KRATOS_WATCH(eig.betas());

        // test if some beta is close to zero 
        // could yield a zero division error for right eigenvalues -> left eigenvalue should be build
        vector_t b_test = eig.betas();
        for(int i=0; i<nroot; i++){
            if(std::abs(b_test(i))<tolerance){
                KRATOS_WARNING("GenEigensystemSolver") << "Some beta value is close to zero! Decrease tolerance or consider left eigenvalues." << std::endl;
                // TODO: maybe even an error instead of warning
                // this code only supports right eigenvalues
                // either extend the code or use another solver
                break;
            }
        }


        // resize Eigenvalue vector and Eigenvector matrix
        // store complex result in double notation
        // Eigenvalues      r1 i1 r2 i2 r3 i3 ...      (real part of eigenvalue 1, imaginary part of eigenvalue 1, etc.)
        // Eigenvectors     same structure, each row corresponds to one Eigenvector //TODO: check, if this is correct
        if (static_cast<int>(rEigenvalues.size()) != nroot*2) {
            rEigenvalues.resize(nroot*2);
        }

        if (eigenvectors_needed){
            if (static_cast<int>(rEigenvectors.size1()) != nroot || static_cast<int>(rEigenvectors.size2()) != nroot*2) {
                rEigenvectors.resize(nroot, nroot*2);
            }
        }

 
        //TODO: consider OMP for size>11; otherwise too much overhead for parallelization
        // store the Eigenvalues
        Eigen::VectorXcd eigvals_tmp = eig.eigenvalues();
        for(int i=0; i<nroot; i++){
            rEigenvalues(2*i) = eigvals_tmp(i).real();
            rEigenvalues(2*i+1) = eigvals_tmp(i).imag();
        }

        // store the Eigenvectors
        if (eigenvectors_needed){
            Eigen::MatrixXcd eigvecs_tmp = eig.eigenvectors();
            for(int i=0; i<nroot; i++){
                for(int j=0; j<nroot; j++){
                    rEigenvectors(i,2*j) = eigvecs_tmp(i,j).real();
                    rEigenvectors(i,2*j+1) = eigvecs_tmp(i,j).imag();
                }
            }
        }
    



 











        
/* 
        eigvals = eig.eigenvalues().head(nroot);

        for (int i = 0; i != nroot; ++i) {
            tmp = r.col(i);
            eigvecs.row(i) = solver.solve(tmp).normalized();
        }
*/

// TODO: could be reused
/*

        // --- normalization
        // Given generalized eigenvalue problem (A - eigenvalue * B) * eigenvector = 0,
        // eigenvector is normalized such that eigenvector^T * B * eigenvector = 1
        if(mParam["normalize_eigenvectors"].GetBool())
        {
            for (int i = 0; i != nroot; ++i)
            {
                const double tmp = eigvecs.row(i) * b * eigvecs.row(i).transpose();
                const double factor = 1.0 / std::sqrt(tmp);
                eigvecs.row(i) *=  factor;
                KRATOS_INFO_IF("GenEigensystemSolver:", echo_level > 0) << "Eigenvector " << i+1 << " is normalized - used factor: " << factor << std::endl;
            }
        }
*/

        // --- output
        if (echo_level > 0) {
            double end_time = OpenMPUtils::GetCurrentTime();
            double duration = end_time - start_time;

            Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "[ ", " ]");

            KRATOS_INFO("GenEigensystemSolver") << "Completed in " << duration << " seconds" << std::endl
                      << "                      Eigenvalues = " << eig.eigenvalues().transpose().format(fmt) << std::endl;
        }
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "GenEigensystemSolver";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream &rOStream) const override
    {
    }

}; // class GenEigensystemSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >>(
    std::istream& rIStream,
    GenEigensystemSolver<TSparseSpaceType,
    TDenseSpaceType,
    TReordererType>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator <<(
    std::ostream& rOStream,
    const GenEigensystemSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

#endif // defined(KRATOS_GEN_EIGENSYSTEM_SOLVER_H_INCLUDED)
